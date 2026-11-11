// Multi-cell model of Rommelfanger & MacLean (2021): many copies of the GATA1/PU.1 cell
// (Eqns 1-3, integrated with Model.rk4Step) coupled by Poisson-process signals (Eqns 4-5).
// Pure logic, no page code, so it can be unit-tested.
//
// Signaling rules implemented here (the paper's Methods, plus two interpretations marked [*]):
//   - Every directed edge sender -> receiver has its OWN Poisson clock: waits are Exp(mean mu),
//     rounded to whole numbers (min 1).  The wait for the next signal starts when the previous
//     signal's pulse has ended, so the cycle is: wait -> signal -> pulse of length tau -> wait ... [*]
//     (Checked against the paper: with this rule a 2-cell chain at lambda=1 has its cell-2 midpoint at
//     A0 = 0.696 and a 2-cell loop at 0.558, versus ~0.69 and ~0.57 in Figs S11 / S2A. Starting the
//     wait at the signal instead gave 0.667 and 0.506.)
//   - A signal at time ts starts a pulse lasting tau time units.  Its factor is frozen at ts:
//       consensus (Eqn 4):  G_sender / (lambda + P_sender)
//       dissensus (Eqn 5):  P_sender / (kappa  + G_sender)
//   - A receiver's A is A0 times the product of the factors of all its currently active pulses
//     (A = A0 with none).  One pulse gives Eqn 4/5 exactly; pulses from DIFFERENT senders that overlap
//     compound, which is how the paper's "A2(ts)" (the value of A2 at the signal time) behaves. [*]
//     (Pulses on the same edge cannot overlap under the rule above.)

(function (root) {
  'use strict';
  const Model = (typeof require === 'function') ? require('./model.js') : root.Model;

  // ---------- topologies ----------
  // { name, n, edges:[{from,to,type:'consensus'|'dissensus'}], layout:[{x,y}] }   (cells numbered
  // from 0 internally; the text format and the page use 1-based numbers.)  layout is in the unit
  // square with y pointing down.
  const C = 'consensus', D = 'dissensus';
  const edge = (from, to, type) => ({ from, to, type: type || C });
  const line = n => Array.from({ length: n }, (_, i) => ({ x: n === 1 ? 0.5 : 0.08 + 0.84 * i / (n - 1), y: 0.5 }));
  const ring = n => Array.from({ length: n }, (_, i) => ({ x: 0.5 + 0.38 * Math.cos(-Math.PI / 2 + 2 * Math.PI * i / n), y: 0.5 + 0.38 * Math.sin(-Math.PI / 2 + 2 * Math.PI * i / n) }));
  const TRI = [{ x: 0.2, y: 0.72 }, { x: 0.5, y: 0.28 }, { x: 0.8, y: 0.72 }];          // cell 1 bottom-left, 2 top, 3 bottom-right (as in Fig. 3)
  const check = n => { if (!(n >= 1 && Number.isInteger(n))) throw new Error('cell count must be a whole number >= 1'); return n; };

  const Topology = {
    single: () => ({ name: 'Single cell', n: 1, edges: [], layout: [{ x: 0.5, y: 0.5 }] }),
    chain: n => { check(n); return { name: `Chain of ${n}`, n, edges: range(n - 1).map(i => edge(i, i + 1)), layout: line(n) }; },
    loop: n => { check(n); return { name: `Loop of ${n}`, n, edges: range(n).filter(() => n > 1).map(i => edge(i, (i + 1) % n)), layout: ring(n) }; },
    bidirectionalChain: n => { check(n); return { name: `Bidirectional chain of ${n}`, n, edges: range(n - 1).flatMap(i => [edge(i, i + 1), edge(i + 1, i)]), layout: line(n) }; },
    fullyConnected: n => { check(n); const e = []; for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) if (i !== j) e.push(edge(i, j)); return { name: `Fully connected, ${n}`, n, edges: e, layout: ring(n) }; },
    // Fig. 3 motifs (panel letters) and Fig. S3
    dissensusChain2: () => ({ name: 'Two cells, dissensus (Fig. 3C)', n: 2, edges: [edge(0, 1, D)], layout: line(2) }),
    loop2Dissensus: () => ({ name: 'Two-cell loop, one dissensus (Fig. 3D)', n: 2, edges: [edge(0, 1, C), edge(1, 0, D)], layout: line(2) }),
    feedForward3: () => ({ name: 'Feed-forward (Fig. 3E)', n: 3, edges: [edge(0, 1), edge(1, 2), edge(0, 2)], layout: TRI }),
    incoherentFeedForward3: () => ({ name: 'Incoherent feed-forward (Fig. 3F)', n: 3, edges: [edge(0, 1), edge(1, 2), edge(0, 2, D)], layout: TRI }),
    loop3TwoDissensus: () => ({ name: 'Three-cell loop, two dissensus (Fig. 3G)', n: 3, edges: [edge(0, 1, D), edge(1, 2, D), edge(2, 0, C)], layout: TRI }),
    loop3OneInverse: () => ({ name: 'Three-cell loop, one inverse (Fig. S3)', n: 3, edges: [edge(0, 1), edge(1, 2), edge(2, 0, D)], layout: TRI }),
  };
  function range(n) { return Array.from({ length: Math.max(0, n) }, (_, i) => i); }

  // "1>2, 2>3, 3-|1"  ->  { n, edges }.  ">" is a consensus arrow, "-|" a dissensus arrow.
  // Separators: commas, semicolons or newlines.  Cells are numbered from 1.
  function parseEdges(text, minCells) {
    const edges = []; let maxIdx = Math.max(1, minCells || 1);
    const items = String(text).split(/[,;\n]+/).map(s => s.trim()).filter(Boolean);
    for (const it of items) {
      const m = /^(\d+)\s*(>|-\|)\s*(\d+)$/.exec(it);
      if (!m) throw new Error(`Cannot read "${it}". Use e.g. 1>2 (consensus) or 3-|1 (dissensus).`);
      const a = +m[1], b = +m[3];
      if (a < 1 || b < 1) throw new Error(`Cell numbers start at 1 (in "${it}").`);
      if (a === b) throw new Error(`A cell cannot signal itself ("${it}").`);
      edges.push({ from: a - 1, to: b - 1, type: m[2] === '>' ? C : D });
      maxIdx = Math.max(maxIdx, a, b);
    }
    if (maxIdx > 200) throw new Error('At most 200 cells.');
    return { n: maxIdx, edges };
  }
  function customTopology(text, minCells) {
    const { n, edges } = parseEdges(text, minCells);
    return { name: 'Custom', n, edges, layout: ring(n) };
  }

  // ---------- random numbers ----------
  // mulberry32: a small, well-known seeded generator; the same seed always gives the same run.
  function mulberry32(seed) {
    let a = seed >>> 0;
    return function () {
      a = (a + 0x6D2B79F5) >>> 0;
      let t = a; t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  // ---------- the network ----------
  const EPS = 1e-9;
  function createNetwork(topology, params, opts) {
    opts = opts || {};
    // Signaling constants. If the caller passes an object it is used directly (not copied), so the
    // page's sliders can edit lambda/kappa/mu/tau while the simulation runs; missing keys get defaults.
    const sig = opts.signal || {};
    for (const [k, v] of Object.entries({ lambda: 1, kappa: 1, mu: 50, tau: 5 })) if (sig[k] == null) sig[k] = v;
    const net = {
      topology, params, sig,
      onA: null,        // optional hook onA(cell, time, newA), called when a pulse start/end changes a cell's A
      t: 0, seed: opts.seed == null ? 1 : opts.seed,
      y: [], A: [], cellParams: [], pulses: [], next: [], events: [], rng: null,
    };
    const n = topology.n, E = topology.edges;

    const wait = () => Math.max(1, Math.round(-net.sig.mu * Math.log(1 - net.rng())));

    net.reset = function (init, seed) {
      if (seed != null) net.seed = seed;
      net.rng = mulberry32(net.seed); net.t = 0; net.events.length = 0;
      net.y = range(n).map(() => (init || [20, 1, 20]).slice());          // paper: (G,P,X) = (20,1,20)
      net.pulses = range(n).map(() => []);
      net.cellParams = range(n).map(() => Object.create(params));        // inherits live slider values; own A
      net.next = E.map(() => wait());
      applyA();
    };

    function applyA() {                                                   // A_k = A0 x product of active pulse factors
      for (let k = 0; k < n; k++) {
        let a = params.A; for (const p of net.pulses[k]) a *= p.factor;
        const changed = a !== net.A[k];
        net.A[k] = a; net.cellParams[k].A = a;
        if (changed && net.onA) net.onA(k, net.t, a);
      }
    }
    net.refreshA = applyA;                                                // re-read A0 (e.g. after the user moves the A slider)

    // Process everything due at (or before) the current time: pulses ending, then signals firing.
    function fireDue() {
      let changed = false;
      for (let k = 0; k < n; k++) {
        const ps = net.pulses[k];
        for (let i = ps.length - 1; i >= 0; i--) if (ps[i].end <= net.t + EPS) { ps.splice(i, 1); changed = true; }
      }
      for (let i = 0; i < E.length; i++) {
        while (net.next[i] <= net.t + EPS) {
          const e = E[i], s = net.y[e.from], G = s[0], P = s[1];
          const factor = e.type === C ? G / (net.sig.lambda + P) : P / (net.sig.kappa + G);
          net.pulses[e.to].push({ end: net.next[i] + net.sig.tau, factor, edge: i });
          net.events.push({ edge: i, t: net.next[i], factor });
          net.next[i] += net.sig.tau + wait(); changed = true;     // next wait starts when this pulse ends
        }
      }
      if (changed) applyA();
    }
    function nextEventTime() {
      let m = Infinity;
      for (const t of net.next) if (t < m) m = t;
      for (const ps of net.pulses) for (const p of ps) if (p.end < m) m = p.end;
      return m;
    }

    // Advance by `len` time units, stopping exactly at every signal start / pulse end.
    // With no event inside the span this is a single RK4 step of size len (so a network with no
    // edges reproduces the single-cell integrator exactly).
    net.advanceBy = function (len) {
      let rem = len;
      while (rem > EPS) {
        fireDue();
        const toEvent = nextEventTime() - net.t;
        const seg = toEvent >= rem - EPS ? rem : toEvent;
        const aNow = params.A;                                            // A0 may have been moved by the user
        for (let k = 0; k < n; k++) {
          let a = aNow; for (const p of net.pulses[k]) a *= p.factor;
          net.A[k] = a; net.cellParams[k].A = a;
          Model.rk4Step(net.y[k], seg, net.cellParams[k]);
        }
        net.t += seg; rem -= seg;
        const te = nextEventTime(); if (Math.abs(net.t - te) < EPS) net.t = te;   // land exactly on the event
      }
      fireDue();
    };
    net.advanceSteps = function (count, h) { for (let i = 0; i < count; i++) net.advanceBy(h); };
    net.advanceTo = function (tEnd, h) { while (net.t < tEnd - EPS) net.advanceBy(Math.min(h, tEnd - net.t)); };
    net.drainEvents = function () { const e = net.events.slice(); net.events.length = 0; return e; };

    net.reset(opts.init, net.seed);
    return net;
  }

  const api = { Topology, parseEdges, customTopology, createNetwork, mulberry32 };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.Network = api;
})(typeof window !== 'undefined' ? window : globalThis);
