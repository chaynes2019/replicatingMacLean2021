// Minimal MP4 muxer: packs H.264 frames that were already compressed (e.g. by the
// browser's WebCodecs VideoEncoder) into a playable .mp4 file. Single video track, no audio,
// no B-frames, constant frame rate. No dependencies and no DOM, so it can be unit-tested.
//
//   Mp4Mux.build({ width, height, fps, avcC, samples }) -> Uint8Array
//     avcC    Uint8Array  the AVCDecoderConfigurationRecord (SPS/PPS), i.e. the
//                         `decoderConfig.description` the encoder reports with its first frame
//     samples [{ data: Uint8Array (AVCC-format frame), key: boolean }]
//
// File layout:  ftyp | moov | mdat.  The moov box (the "index": where every frame is, how long
// it lasts, which are keyframes) goes before the frame data so players can start immediately.

(function (root) {
  'use strict';

  const u8 = (...bytes) => Uint8Array.from(bytes);
  const u16 = v => u8((v >>> 8) & 255, v & 255);
  const u32 = v => u8((v >>> 24) & 255, (v >>> 16) & 255, (v >>> 8) & 255, v & 255);
  const str = s => Uint8Array.from(s, ch => ch.charCodeAt(0));   // box names are ASCII
  const zeros = n => new Uint8Array(n);

  function concat(parts) {
    let n = 0; for (const p of parts) n += p.length;
    const out = new Uint8Array(n); let o = 0;
    for (const p of parts) { out.set(p, o); o += p.length; }
    return out;
  }
  // An MP4 "box": 4-byte size (including this header), 4-char type, payload.
  const box = (type, ...payload) => {
    const body = concat(payload);
    return concat([u32(8 + body.length), str(type), body]);
  };
  const fullBox = (type, version, flags, ...payload) =>
    box(type, u8(version, (flags >>> 16) & 255, (flags >>> 8) & 255, flags & 255), ...payload);

  const MATRIX = concat([u32(0x00010000), u32(0), u32(0), u32(0), u32(0x00010000), u32(0), u32(0), u32(0), u32(0x40000000)]);

  function moov(o, chunkOffset) {
    const { width, height, avcC, samples } = o;
    const n = samples.length;
    const timescale = o.fps * 1000, sampleDelta = 1000;       // exactly 1/fps per frame
    const durationMs = Math.round(n * 1000 / o.fps);
    const sizes = concat(samples.map(s => u32(s.data.length)));
    const keys = []; samples.forEach((s, i) => { if (s.key) keys.push(i + 1); });

    const mvhd = fullBox('mvhd', 0, 0, u32(0), u32(0), u32(1000), u32(durationMs), u32(0x00010000), u16(0x0100),
      zeros(10), MATRIX, zeros(24), u32(2));
    const tkhd = fullBox('tkhd', 0, 3, u32(0), u32(0), u32(1), u32(0), u32(durationMs), zeros(8),
      u16(0), u16(0), u16(0), u16(0), MATRIX, u32(width << 16), u32(height << 16));
    const mdhd = fullBox('mdhd', 0, 0, u32(0), u32(0), u32(timescale), u32(n * sampleDelta), u16(0x55c4), u16(0));
    const hdlr = fullBox('hdlr', 0, 0, u32(0), str('vide'), zeros(12), str('VideoHandler\0'));
    const vmhd = fullBox('vmhd', 0, 1, u16(0), zeros(6));
    const dinf = box('dinf', fullBox('dref', 0, 0, u32(1), fullBox('url ', 0, 1)));
    const avc1 = box('avc1', zeros(6), u16(1), u16(0), u16(0), zeros(12), u16(width), u16(height),
      u32(0x00480000), u32(0x00480000), u32(0), u16(1), zeros(32), u16(0x0018), u16(0xffff), box('avcC', avcC));
    const stsd = fullBox('stsd', 0, 0, u32(1), avc1);
    const stts = fullBox('stts', 0, 0, u32(1), u32(n), u32(sampleDelta));
    const stss = fullBox('stss', 0, 0, u32(keys.length), concat(keys.map(u32)));
    const stsc = fullBox('stsc', 0, 0, u32(1), u32(1), u32(n), u32(1));          // all samples in one chunk
    const stsz = fullBox('stsz', 0, 0, u32(0), u32(n), sizes);
    const stco = fullBox('stco', 0, 0, u32(1), u32(chunkOffset));
    const stbl = box('stbl', stsd, stts, stss, stsc, stsz, stco);
    const minf = box('minf', vmhd, dinf, stbl);
    const mdia = box('mdia', mdhd, hdlr, minf);
    return box('moov', mvhd, box('trak', tkhd, mdia));
  }

  function build(o) {
    if (!o.samples.length) throw new Error('Mp4Mux.build: no samples');
    const ftyp = box('ftyp', str('isom'), u32(0x200), str('isom'), str('iso2'), str('avc1'), str('mp41'));
    // stco holds an absolute file offset, which depends on the moov size, which does not
    // depend on the offset value: build once to measure, then again with the real offset.
    const moovLen = moov(o, 0).length;
    const mdatHeader = 8;
    const moovBox = moov(o, ftyp.length + moovLen + mdatHeader);
    let total = 0; for (const s of o.samples) total += s.data.length;
    const mdat = concat([u32(mdatHeader + total), str('mdat'), ...o.samples.map(s => s.data)]);
    return concat([ftyp, moovBox, mdat]);
  }

  const api = { build };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.Mp4Mux = api;
})(typeof window !== 'undefined' ? window : globalThis);
