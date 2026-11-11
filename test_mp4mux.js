// Run: node test_mp4mux.js   (or via the osascript harness if Node is not installed)
// Builds an MP4 from synthetic frames and checks the container's structure.
const Mp4Mux = require('./mp4mux.js');
let fail = 0;
const check = (name, ok, detail) => { console.log((ok ? 'PASS ' : 'FAIL ') + name + (detail ? '  ' + detail : '')); if (!ok) fail++; };

const FPS = 30, W = 640, H = 360, N = 70;
const samples = [];
for (let i = 0; i < N; i++) samples.push({ data: Uint8Array.from({ length: 10 + i }, (_, k) => (i * 7 + k) & 255), key: i % 30 === 0 });
const avcC = Uint8Array.from([1, 0x64, 0, 0x33, 0xff, 0xe1, 0, 2, 0x67, 0x64, 0x01, 0, 1, 0x68]);
const file = Mp4Mux.build({ width: W, height: H, fps: FPS, avcC, samples });

const dv = new DataView(file.buffer, file.byteOffset, file.byteLength);
const u32 = o => dv.getUint32(o), type = o => String.fromCharCode(file[o], file[o + 1], file[o + 2], file[o + 3]);
const CONTAINERS = ['moov', 'trak', 'mdia', 'minf', 'dinf', 'stbl'];
function parse(start, end) { // -> [{type, start, size, kids}]
  const out = []; let o = start;
  while (o < end) {
    const size = u32(o), t = type(o + 4);
    if (size < 8 || o + size > end) throw new Error('bad box at ' + o + ' (' + t + ', size ' + size + ')');
    out.push({ type: t, start: o, size, kids: CONTAINERS.includes(t) ? parse(o + 8, o + size) : [] });
    o += size;
  }
  return out;
}
const find = (list, ...path) => path.reduce((l, p) => l && l.kids ? l.kids.find(b => b.type === p) : (l.find(b => b.type === p)), list);

let top;
try { top = parse(0, file.length); check('boxes tile the file exactly', true); } catch (e) { check('boxes tile the file exactly', false, e.message); process.exit(1); }
check('top-level order ftyp, moov, mdat', top.map(b => b.type).join() === 'ftyp,moov,mdat', top.map(b => b.type).join());

const stbl = find(top, 'moov', 'trak', 'mdia', 'minf', 'stbl'), mdat = top[2];
const body = b => b.start + 8;                       // start of a box's payload
const stts = find(stbl.kids, 'stts'), stsz = find(stbl.kids, 'stsz'), stss = find(stbl.kids, 'stss'), stco = find(stbl.kids, 'stco');
check('stts: N samples of 1000 ticks', u32(body(stts) + 8) === N && u32(body(stts) + 12) === 1000);
check('stsz: N sizes match input', u32(body(stsz) + 8) === N && samples.every((s, i) => u32(body(stsz) + 12 + 4 * i) === s.data.length));
check('stss: keyframes at samples 1,31,61', u32(body(stss) + 4) === 3 && [1, 31, 61].every((v, i) => u32(body(stss) + 8 + 4 * i) === v));
const total = samples.reduce((a, s) => a + s.data.length, 0);
check('mdat payload == sum of sample sizes', mdat.size - 8 === total);
const off = u32(body(stco) + 8);
check('stco points at first frame data', off === body(mdat) && samples[0].data.every((v, k) => file[off + k] === v), 'offset ' + off);
const lastOff = off + total - samples[N - 1].data.length;
check('last frame data is where sizes say', samples[N - 1].data.every((v, k) => file[lastOff + k] === v));
const mdhd = find(top, 'moov', 'trak', 'mdia', 'mdhd');
check('mdhd: timescale 30000, duration N*1000 (=N/fps s)', u32(body(mdhd) + 12) === FPS * 1000 && u32(body(mdhd) + 16) === N * 1000);
const mvhd = find(top, 'moov', 'mvhd');
check('mvhd: duration in ms ~ 2333', u32(body(mvhd) + 16) === Math.round(N * 1000 / FPS), String(u32(body(mvhd) + 16)));
const avcCBox = find(stbl.kids, 'stsd');
check('avcC bytes embedded', (() => { for (let o = avcCBox.start; o < avcCBox.start + avcCBox.size - 4; o++) if (type(o) === 'avcC') return avcC.every((v, k) => file[o + 4 + k] === v); return false; })());
process.exit(fail ? 1 : 0);
