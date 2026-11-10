// Run: node test_model.js
// Checks the RK4 integrator against an exact solution and the model against the paper.
const { DEFAULTS, derivs, rk4Step } = require('./model.js');
let fail = 0;
const check = (name, ok, detail) => { console.log((ok ? 'PASS ' : 'FAIL ') + name + (detail ? '  ' + detail : '')); if (!ok) fail++; };

// 1. Order of convergence on a linear case: all production off => dG/dt = -b1 G, exact G0*exp(-b1 t).
//    Halving h should cut the error by ~2^4 = 16 for RK4.
function err(h) {
  const p = { ...DEFAULTS, a1: 0, a2: 0, a3: 0, a4: 0, a5: 0, b1: 0.5, A: 0, B: 0 };
  const y = [10, 0, 0];
  const T = 4; for (let i = 0; i < Math.round(T / h); i++) rk4Step(y, h, p);
  return Math.abs(y[0] - 10 * Math.exp(-p.b1 * T));
}
const e1 = err(0.4), e2 = err(0.2);
check('RK4 is 4th order', Math.abs(Math.log2(e1 / e2) - 4) < 0.3, `error ratio = ${(e1 / e2).toFixed(1)} (expect ~16)`);

// 2. Bistability with Table S1 values.
function settle(y0) {
  const y = y0.slice(); for (let i = 0; i < 4e5; i++) rk4Step(y, 0.5, DEFAULTS); return y;
}
const hi = settle([60, 0.2, 0]), lo = settle([1, 5, 0]);
check('G-high state ~ (53, 0.22)  [paper Fig. S1]', Math.abs(hi[0] - 53) < 1 && Math.abs(hi[1] - 0.22) < 0.02,
  hi.map(v => v.toFixed(3)).join(', '));
check('G-low / P-high state exists and differs', lo[0] < 20 && lo[1] > 1, lo.map(v => v.toFixed(3)).join(', '));

// 3. Both are true fixed points.
for (const [n, s] of [['high', hi], ['low', lo]]) {
  const d = derivs(s, DEFAULTS, [0, 0, 0]);
  check(`${n} state is a fixed point`, d.every(v => Math.abs(v) < 1e-6), d.map(v => v.toExponential(1)).join(', '));
}
// 4. Bifurcation curve: every point is a fixed point, and at A=1 it contains the two
//    attractors the simulation found plus an unstable branch between them.
const { steadyStateCurve } = require('./model.js');
const curve = steadyStateCurve(DEFAULTS, 3000);
let maxRes = 0;
for (const s of curve) { const d = derivs([s.G, s.P, s.X], { ...DEFAULTS, A: s.A }, [0, 0, 0]); maxRes = Math.max(maxRes, ...d.map(Math.abs)); }
check('curve points are fixed points', maxRes < 1e-8, 'max |f| = ' + maxRes.toExponential(1));
function atA(A0) { // all (G, stable) where the curve crosses A = A0
  const r = [];
  for (let i = 1; i < curve.length; i++) {
    const u = curve[i - 1], v = curve[i];
    if ((u.A - A0) * (v.A - A0) < 0) { const w = (A0 - u.A) / (v.A - u.A); r.push({ G: u.G + w * (v.G - u.G), stable: u.stable && v.stable }); }
  }
  return r;
}
const at1 = atA(1.0);
console.log('   branches at A=1: ' + at1.map(b => `G=${b.G.toFixed(2)}${b.stable ? ' stable' : ' unstable'}`).join(', '));
check('A=1: stable G~13.4 and G~52.7 exist',
  at1.some(b => b.stable && Math.abs(b.G - 13.38) < 0.1) && at1.some(b => b.stable && Math.abs(b.G - 52.68) < 0.1));
check('A=1: an unstable branch lies between them', at1.some(b => !b.stable && b.G > 14 && b.G < 52));
process.exit(fail ? 1 : 0);
