// Tests for network.js.  Run via the harness (Node is optional): node test_network.js
const Model = require('./model.js');
const Network = require('./network.js');
const { Topology, parseEdges, customTopology, createNetwork } = Network;
let fail = 0;
const check = (name, ok, detail) => { console.log((ok ? 'PASS ' : 'FAIL ') + name + (detail ? '  ' + detail : '')); if (!ok) fail++; };
const throwsMsg = f => { try { f(); return null; } catch (e) { return e.message; } };
const close = (a, b, tol) => Math.abs(a - b) <= tol * Math.max(1, Math.abs(b));
const P0 = () => ({ ...Model.DEFAULTS });

// ---- text edge lists ----
{
  const r = parseEdges('1>2, 2>3, 3-|1');
  check('parse "1>2, 2>3, 3-|1"', r.n === 3 && r.edges.length === 3 && r.edges[0].type === 'consensus' && r.edges[2].type === 'dissensus' && r.edges[2].from === 2 && r.edges[2].to === 0);
  check('parse infers n from the largest number; minCells raises it', parseEdges('1>4').n === 4 && parseEdges('1>2', 6).n === 6);
  check('parse accepts newlines/semicolons/spaces', parseEdges('1 > 2;\n2 -| 1').edges.length === 2);
  check('parse rejects garbage', /Cannot read/.test(throwsMsg(() => parseEdges('1>>2')) || ''));
  check('parse rejects self-signals', /itself/.test(throwsMsg(() => parseEdges('2>2')) || ''));
  check('parse rejects cell 0', /start at 1/.test(throwsMsg(() => parseEdges('0>1')) || ''));
  check('customTopology has a layout entry per cell', customTopology('1>2,2>3').layout.length === 3);
}
// ---- builders ----
{
  const ec = t => t.edges.filter(e => e.type === 'consensus').length, ed = t => t.edges.filter(e => e.type === 'dissensus').length;
  check('single: 1 cell, no edges', Topology.single().n === 1 && Topology.single().edges.length === 0);
  check('chain(5): 4 edges 1->2..4->5', (t => t.edges.length === 4 && t.edges.every((e, i) => e.from === i && e.to === i + 1))(Topology.chain(5)));
  check('loop(5): 5 edges incl. 5->1', (t => t.edges.length === 5 && t.edges[4].from === 4 && t.edges[4].to === 0)(Topology.loop(5)));
  check('loop(2): 1->2 and 2->1', (t => t.edges.length === 2)(Topology.loop(2)));
  check('bidirectionalChain(4): 6 edges', Topology.bidirectionalChain(4).edges.length === 6);
  check('fullyConnected(5): 20 edges', Topology.fullyConnected(5).edges.length === 20);
  check('Fig 3C: 1 dissensus edge', ed(Topology.dissensusChain2()) === 1 && ec(Topology.dissensusChain2()) === 0);
  check('Fig 3D: one consensus 1->2, one dissensus 2->1', (t => ec(t) === 1 && ed(t) === 1 && t.edges[1].from === 1 && t.edges[1].type === 'dissensus')(Topology.loop2Dissensus()));
  check('Fig 3E: 3 consensus edges incl. 1->3', (t => ec(t) === 3 && t.edges.some(e => e.from === 0 && e.to === 2))(Topology.feedForward3()));
  check('Fig 3F: 1-|3 dissensus', (t => ec(t) === 2 && t.edges.some(e => e.from === 0 && e.to === 2 && e.type === 'dissensus'))(Topology.incoherentFeedForward3()));
  check('Fig 3G: two dissensus (1-|2, 2-|3), 3->1', (t => ed(t) === 2 && ec(t) === 1 && t.edges[2].from === 2 && t.edges[2].to === 0)(Topology.loop3TwoDissensus()));
  check('every layout has one position per cell, inside the unit square', [Topology.chain(7), Topology.loop(9), Topology.fullyConnected(6), Topology.feedForward3()].every(t => t.layout.length === t.n && t.layout.every(p => p.x >= 0 && p.x <= 1 && p.y >= 0 && p.y <= 1)));
  check('cell count must be a whole number >= 1', !!throwsMsg(() => Topology.chain(0)) && !!throwsMsg(() => Topology.loop(2.5)));
}
// ---- regression: a network with no edges is exactly the old single-cell integrator ----
{
  const p = P0(), net = createNetwork(Topology.single(), p, { seed: 7 });
  const y = [20, 1, 20], q = P0();
  net.advanceSteps(8000, 0.5); for (let i = 0; i < 8000; i++) Model.rk4Step(y, 0.5, q);
  check('no-edge network == plain RK4, bit for bit', net.y[0][0] === y[0] && net.y[0][1] === y[1] && net.y[0][2] === y[2] && net.t === 4000, `G=${net.y[0][0]}`);
  net.advanceSteps(500, 0.5); p.A = 1.5; net.advanceSteps(500, 0.5); q.A = 1.5;
  for (let i = 0; i < 500; i++) Model.rk4Step(y, 0.5, q);    // note: first 500 above ran at A=1.0
  const y2 = [20, 1, 20], q2 = P0(); for (let i = 0; i < 8500; i++) Model.rk4Step(y2, 0.5, q2); for (let i = 0; i < 500; i++) { q2.A = 1.5; Model.rk4Step(y2, 0.5, q2); }
  check('moving A0 live changes the cell on the very next step (same as before)', net.y[0][0] === y2[0] && net.A[0] === 1.5);
}
// ---- signal formulas (Eqns 4 and 5) ----
{
  const p = P0(), t = Topology.chain(2), net = createNetwork(t, p, { seed: 3, signal: { lambda: 38, mu: 1e9 } });
  net.next[0] = 0; net.advanceBy(0.5);          // fires at t=0 with the sender at (20, 1, 20)
  const ev = net.drainEvents();
  check('Eqn 4 at start state, lambda=38: factor = 20/39 = 0.5128 (cf. Fig. S9 ~0.55)', ev.length === 1 && close(ev[0].factor, 20 / 39, 1e-12), String(ev[0] && ev[0].factor));
  check('receiver A = A0 x factor during the pulse', close(net.A[1], p.A * 20 / 39, 1e-12), String(net.A[1]));
}
{
  const p = P0(), net = createNetwork(Topology.dissensusChain2(), p, { seed: 3, signal: { kappa: 28, mu: 1e9 } });
  net.next[0] = 0; net.advanceBy(0.5); const ev = net.drainEvents();
  check('Eqn 5 at start state, kappa=28: factor = P/(kappa+G) = 1/48', ev.length === 1 && close(ev[0].factor, 1 / 48, 1e-12), String(ev[0] && ev[0].factor));
}
// ---- pulse timing: [ts, ts+tau) ----
{
  const p = P0(), net = createNetwork(Topology.chain(2), p, { seed: 3, signal: { lambda: 1, tau: 5, mu: 1e9 } });
  net.next[0] = 10;
  net.advanceTo(9.5, 0.5); const before = net.A[1];
  net.advanceTo(10, 0.5); const atStart = net.A[1];
  net.advanceTo(14.5, 0.5); const mid = net.A[1];
  net.advanceTo(15, 0.5); const atEnd = net.A[1];
  check('A = A0 before the signal', before === p.A);
  check('A is changed from t = ts', atStart !== p.A && Math.abs(net.events.length) === 1);
  check('A still changed just before ts + tau', mid === atStart);
  check('A back to exactly A0 at ts + tau (half-open interval)', atEnd === p.A);
  const net2 = createNetwork(Topology.chain(2), P0(), { seed: 3, signal: { tau: 5, mu: 1e9 } });
  net2.next[0] = 10; net2.advanceTo(10, 0.37);
  check('events land exactly on the signal time even when h does not divide it', net2.t === 10 && net2.events.length === 1);
}
// ---- overlapping pulses multiply ----
{
  const p = P0(), t = { name: 't', n: 3, edges: [{ from: 0, to: 2, type: 'consensus' }, { from: 1, to: 2, type: 'consensus' }], layout: [] };
  const net = createNetwork(t, p, { seed: 3, signal: { lambda: 1, tau: 5, mu: 1e9 } });
  net.next[0] = 10; net.next[1] = 12; net.advanceTo(13, 0.5);
  const ev = net.drainEvents();
  check('two pulses overlap: A = A0 x f1 x f2', ev.length === 2 && close(net.A[2], p.A * ev[0].factor * ev[1].factor, 1e-12), net.A[2].toFixed(5));
  net.advanceTo(15, 0.5); check('first ends at 15, second still active', close(net.A[2], p.A * ev[1].factor, 1e-12));
  net.advanceTo(17, 0.5); check('both ended: A = A0', net.A[2] === p.A);
}
// ---- Poisson clock statistics ----
{
  const p = P0(), net = createNetwork(Topology.chain(2), p, { seed: 11, signal: { mu: 50, tau: 5 } });
  const ts = []; for (let i = 0; i < 400; i++) { net.advanceTo(net.t + 250, 1); for (const e of net.drainEvents()) ts.push(e.t); }
  // cycle = pulse (tau) + wait, so the wait before signal i is (ts[i] - ts[i-1] - tau); the first wait starts at t = 0
  const waits = ts.map((t, i) => i ? t - ts[i - 1] - 5 : t);
  const mean = waits.reduce((a, b) => a + b, 0) / waits.length;
  check('signal times are whole numbers', ts.every(t => Number.isInteger(t)));
  check('waits >= 1 and mean ~ mu = 50', waits.every(w => w >= 1) && Math.abs(mean - 50) < 3, `${waits.length} signals, mean wait ${mean.toFixed(2)}`);
  const sd = Math.sqrt(waits.reduce((a, w) => a + (w - mean) ** 2, 0) / waits.length);
  check('waits look exponential (sd ~ mean)', Math.abs(sd / mean - 1) < 0.1, `sd/mean = ${(sd / mean).toFixed(2)}`);
  check('pulses on one edge never overlap (next wait starts after the pulse ends)', ts.every((t, i) => !i || t - ts[i - 1] >= 5 + 1));
}
// ---- independent clock per edge ----
{
  const p = P0(), t = { name: 't', n: 3, edges: [{ from: 0, to: 1, type: 'consensus' }, { from: 0, to: 2, type: 'consensus' }], layout: [] };
  const net = createNetwork(t, p, { seed: 5, signal: { mu: 50 } }); net.advanceTo(20000, 1);
  const ev = net.drainEvents(), a = ev.filter(e => e.edge === 0).map(e => e.t), b = ev.filter(e => e.edge === 1).map(e => e.t);
  check('one sender, two targets: the two edges have different signal times', a.length > 100 && b.length > 100 && a.join() !== b.join());
}
// ---- reproducibility and step-size independence of event times ----
{
  const run = (seed, h) => { const n = createNetwork(Topology.loop(3), P0(), { seed, signal: { lambda: 1, mu: 50 } }); n.advanceTo(3000, h); return n; };
  const a = run(9, 0.5), b = run(9, 0.5), c = run(10, 0.5);
  check('same seed -> identical run (exactly)', a.y.every((s, k) => s.every((v, i) => v === b.y[k][i])));
  check('different seed -> different run', a.y.some((s, k) => s.some((v, i) => v !== c.y[k][i])));
  const fine = run(9, 0.1), coarse = run(9, 0.5);
  const ts = n => (n.drainEvents(), n.next.join());           // next signal times depend only on the random stream
  check('signal times do not depend on the step size', ts(fine) === ts(coarse));
  check('end state agrees across step sizes (RK4 error only)', fine.y.every((s, k) => s.every((v, i) => close(v, coarse.y[k][i], 1e-4))), `G1: ${fine.y[0][0].toFixed(6)} vs ${coarse.y[0][0].toFixed(6)}`);
}
// ---- dynamics: a cell receiving consensus signals is pulled toward the sender's fate ----
{
  const p = P0(); p.A = 0.8;        // single cell at A=0.8 settles G-low from (20,1,20) ... checked below
  const alone = createNetwork(Topology.single(), p, { seed: 1 }); alone.advanceTo(40000, 1);
  check('single cell at A0=0.8 from (20,1,20) ends G-low', alone.y[0][0] < 20, 'G=' + alone.y[0][0].toFixed(2));
}
// ---- hook and live signal constants (used by the page) ----
{
  const p = P0(), live = { lambda: 38, mu: 1e9, tau: 5 }, net = createNetwork(Topology.chain(2), p, { seed: 3, signal: live });
  const log = []; net.onA = (k, t, a) => log.push([k, t, a]);
  net.next[0] = 10; net.advanceTo(20, 0.5);
  check('onA fires at pulse start and pulse end, with exact times', log.length === 2 && log[0][0] === 1 && log[0][1] === 10 && log[1][1] === 15 && log[1][2] === p.A, JSON.stringify(log.map(l => [l[0], l[1]])));
  live.lambda = 1; net.next[0] = 30; net.advanceTo(30, 0.5); const ev = net.drainEvents().filter(e => e.t === 30)[0];
  check('editing the live signal object changes the next pulse', ev && close(ev.factor, net.y[0][0] / (1 + net.y[0][1]), 1e-9));
  const before = net.A[1]; p.A = 0.5; net.refreshA();
  check('refreshA picks up a moved A0 slider', net.A[1] !== before && close(net.A[1], 0.5 * (net.pulses[1].reduce((a, q) => a * q.factor, 1)), 1e-12));
}
process.exit(fail ? 1 : 0);
