// Core model of Rommelfanger & MacLean (2021), Eqns 1-3 (after Chickarmane et al. 2009),
// plus a hand-written classical 4th-order Runge-Kutta integrator. No external solver.
//
//   dG/dt = -b1 G + (a1 A + a2 G) / (1 + g1 A + g2 G + g3 G P)                      (1)
//   dP/dt = -b2 P + (a3 B + a4 P) / (1 + g4 B + g5 P + g6 G P + g7 G X)             (2)
//   dX/dt = -b3 X +  a5 G         / (1 + g8 G + g9 C)                               (3)
//
// G = GATA1, P = PU.1, X = cofactor recruited by GATA1 to inhibit PU.1.
// A, B, C summarise the extracellular environment.

(function (root) {
  'use strict';

  // Table S1 of the supplement. A is not in Table S1 (it is the signalled parameter,
  // resting value A0); 1.0 is the value used in Fig. S1.
  const DEFAULTS = {
    a1: 1.0, a2: 0.25, a3: 1.0, a4: 0.25, a5: 0.01,
    b1: 0.01, b2: 0.01, b3: 0.01,
    g1: 1.0, g2: 0.25, g3: 1.0, g4: 1.0, g5: 0.25, g6: 1.0, g7: 0.13, g8: 0.01, g9: 10.0,
    A: 1.0, B: 0.5, C: 0.0,
  };

  // Right-hand side. Reads p on every call, so parameter edits take effect immediately.
  // Writes into `out` to avoid allocating in the inner loop.
  function derivs(y, p, out) {
    const G = y[0], P = y[1], X = y[2];
    out[0] = -p.b1 * G + (p.a1 * p.A + p.a2 * G) /
             (1 + p.g1 * p.A + p.g2 * G + p.g3 * G * P);
    out[1] = -p.b2 * P + (p.a3 * p.B + p.a4 * P) /
             (1 + p.g4 * p.B + p.g5 * P + p.g6 * G * P + p.g7 * G * X);
    out[2] = -p.b3 * X + (p.a5 * G) /
             (1 + p.g8 * G + p.g9 * p.C);
    return out;
  }

  // Scratch buffers, so rk4Step does no allocation.
  const k1 = [0, 0, 0], k2 = [0, 0, 0], k3 = [0, 0, 0], k4 = [0, 0, 0], tmp = [0, 0, 0];

  // One classical RK4 step of size h, in place on y.
  function rk4Step(y, h, p) {
    derivs(y, p, k1);
    for (let i = 0; i < 3; i++) tmp[i] = y[i] + 0.5 * h * k1[i];
    derivs(tmp, p, k2);
    for (let i = 0; i < 3; i++) tmp[i] = y[i] + 0.5 * h * k2[i];
    derivs(tmp, p, k3);
    for (let i = 0; i < 3; i++) tmp[i] = y[i] + h * k3[i];
    derivs(tmp, p, k4);
    for (let i = 0; i < 3; i++) y[i] += (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    // Concentrations cannot go negative; clamp tiny numerical undershoot.
    for (let i = 0; i < 3; i++) if (y[i] < 0) y[i] = 0;
    return y;
  }

  // ---------- Steady states (bifurcation curve) ----------
  // Set dG/dt = dP/dt = dX/dt = 0 and parametrise by G:
  //   X(G) = a5 G / (b3 (1 + g8 G + g9 C))                       (Eqn 3)
  //   P(G): b2 P (1+g4 B+g5 P+g6 G P+g7 G X) = a3 B + a4 P        (Eqn 2, quadratic in P)
  //   A(G): Eqn 1 is LINEAR in A, so A is solved explicitly.      (Eqn 1)
  // This traces every branch (stable and unstable, including folds) exactly, with no
  // time stepping. A in p is ignored; it is the output.
  function pOfG(G, X, p) {
    const a = p.b2 * (p.g5 + p.g6 * G);
    const b = p.b2 * (1 + p.g4 * p.B + p.g7 * G * X) - p.a4;
    const c = -p.a3 * p.B;
    if (a > 1e-14) return (-b + Math.sqrt(Math.max(0, b * b - 4 * a * c))) / (2 * a);
    if (b > 1e-14) return -c / b;
    return NaN;
  }
  function stableAt(s, p) {
    // Numeric 3x3 Jacobian; stable iff all eigenvalues have Re < 0 (Routh-Hurwitz on
    // s^3 + a1 s^2 + a2 s + a3).
    const q = { ...p, A: s.A }, J = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    const f1 = [0, 0, 0], f0 = [0, 0, 0], y = [s.G, s.P, s.X];
    for (let j = 0; j < 3; j++) {
      const h = 1e-6 * Math.max(1, Math.abs(y[j])), yj = y[j];
      y[j] = yj + h; derivs(y, q, f1); y[j] = yj - h; derivs(y, q, f0); y[j] = yj;
      for (let i = 0; i < 3; i++) J[i][j] = (f1[i] - f0[i]) / (2 * h);
    }
    const tr = J[0][0] + J[1][1] + J[2][2];
    const m2 = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) + (J[0][0] * J[2][2] - J[0][2] * J[2][0]) + (J[1][1] * J[2][2] - J[1][2] * J[2][1]);
    const det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
    const a1 = -tr, a2 = m2, a3 = -det;
    return a1 > 0 && a3 > 0 && a1 * a2 > a3;
  }
  // Returns an array of {A, G, P, X, stable}, ordered by increasing G.
  function steadyStateCurve(p, n) {
    n = n || 1500;
    const out = [];
    if (p.b3 <= 0) return out;
    const k = p.b1 * p.g1;
    const Gmax = k > 0 ? Math.min(150, 0.999 * p.a1 / k) : 150;
    for (let i = 0; i <= n; i++) {
      const G = Gmax * i / n;
      const X = p.a5 * G / (p.b3 * (1 + p.g8 * G + p.g9 * p.C));
      const P = pOfG(G, X, p);
      if (!isFinite(P)) continue;
      const A = (p.b1 * G * (1 + p.g2 * G + p.g3 * G * P) - p.a2 * G) / (p.a1 - k * G);
      if (!isFinite(A)) continue;
      const s = { A, G, P, X, stable: false };
      s.stable = stableAt(s, p);
      out.push(s);
    }
    return out;
  }

  const api = { DEFAULTS, derivs, rk4Step, steadyStateCurve };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.Model = api;
})(typeof window !== 'undefined' ? window : globalThis);
