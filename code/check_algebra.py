"""Checks of the algebraic identities used in the paper
"Weak-Type Endpoint Bounds for Bochner-Riesz Means in Every Dimension".

Every identity below is proved by hand in the paper; this program confirms
them symbolically (sympy) or numerically on random data (numpy).  It prints
one line per check and a final count.

Run:  python3 -I check_algebra.py
"""

import itertools
import random
from fractions import Fraction

import numpy as np
import sympy as sp

PASSED = 0
FAILED = 0


def report(name, ok):
    global PASSED, FAILED
    if ok:
        PASSED += 1
        print(f"[pass] {name}")
    else:
        FAILED += 1
        print(f"[FAIL] {name}")


# ---------------------------------------------------------------------------
# Quaternions as 4-tuples (a, b, c, d) = a + b i + c j + d k.
# ---------------------------------------------------------------------------

def qmul(p, q):
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return (a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
            a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
            a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
            a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2)


def qconj(p):
    return (p[0], -p[1], -p[2], -p[3])


def qadd(p, q):
    return tuple(x + y for x, y in zip(p, q))


def qscale(t, p):
    return tuple(t * x for x in p)


ZERO = (0, 0, 0, 0)
UNITS = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]


def qmat_mul(A, B):
    n, m, l = len(A), len(B), len(B[0])
    C = [[ZERO for _ in range(l)] for _ in range(n)]
    for a in range(n):
        for b in range(l):
            s = ZERO
            for c in range(m):
                s = qadd(s, qmul(A[a][c], B[c][b]))
            C[a][b] = s
    return C


def qmat_star(A):
    return [[qconj(A[b][a]) for b in range(len(A))] for a in range(len(A[0]))]


def complex_rep(A):
    """2m x 2m complex matrix of a quaternionic m x m matrix.

    q = z + w j with z = a + b i, w = c + d i  ->  [[z, w], [-conj w, conj z]].
    """
    m = len(A)
    M = np.zeros((2 * m, 2 * m), dtype=complex)
    for a in range(m):
        for b in range(m):
            x = A[a][b]
            z = complex(x[0], x[1])
            w = complex(x[2], x[3])
            M[2 * a, 2 * b] = z
            M[2 * a, 2 * b + 1] = w
            M[2 * a + 1, 2 * b] = -np.conj(w)
            M[2 * a + 1, 2 * b + 1] = np.conj(z)
    return M


def moore_det_numeric(A):
    """Moore determinant of a Hermitian quaternionic matrix: product of the m
    real eigenvalues, each of which appears twice in the complex representation."""
    ev = np.sort(np.linalg.eigvalsh(complex_rep(A)))
    return float(np.prod(ev[::2]))


def random_hermitian(m, rng, pos=False):
    A = [[ZERO for _ in range(m)] for _ in range(m)]
    for a in range(m):
        A[a][a] = (rng.normal(), 0.0, 0.0, 0.0)
        for b in range(a + 1, m):
            x = tuple(rng.normal(size=4))
            A[a][b] = x
            A[b][a] = qconj(x)
    if pos:
        B = qmat_mul(A, qmat_star(A))
        for a in range(m):
            B[a][a] = qadd(B[a][a], (0.5, 0, 0, 0))
        return B
    return A


def moore_det3(A):
    a, b, c = A[0][0][0], A[1][1][0], A[2][2][0]
    x, y, z = A[0][1], A[0][2], A[1][2]
    n2 = lambda q: sum(t * t for t in q)
    re_xzy = qmul(qmul(x, z), qconj(y))[0]
    return a * b * c - a * n2(z) - b * n2(y) - c * n2(x) + 2 * re_xzy


# ---------------------------------------------------------------------------
# 1. Moore determinant identities (Lemma on matrix identities).
# ---------------------------------------------------------------------------

def check_moore():
    rng = np.random.default_rng(1)
    ok3 = True
    for _ in range(50):
        A = random_hermitian(3, rng)
        ok3 &= abs(moore_det3(A) - moore_det_numeric(A)) < 1e-8 * (1 + abs(moore_det3(A)))
    report("3x3 Moore determinant polynomial = product of real eigenvalues", ok3)

    # cofactor C(A) = det(A) A^{-1} as the derivative of det(A + tB)
    okc = True
    for m in (2, 3, 4):
        for _ in range(20):
            A = random_hermitian(m, rng, pos=True)
            B = random_hermitian(m, rng)
            t = 1e-6
            Ap = [[qadd(A[a][b], qscale(t, B[a][b])) for b in range(m)] for a in range(m)]
            Am = [[qadd(A[a][b], qscale(-t, B[a][b])) for b in range(m)] for a in range(m)]
            deriv = (moore_det_numeric(Ap) - moore_det_numeric(Am)) / (2 * t)
            M = complex_rep(A)
            inv = np.linalg.inv(M)
            # (det A) A^{-1} : B = det(A) Re tr(A^{-1} B); Re tr is half the complex trace.
            val = moore_det_numeric(A) * np.trace(inv @ complex_rep(B)).real / 2
            okc &= abs(deriv - val) < 1e-4 * (1 + abs(val))
    report("cofactor (det A) A^{-1} is the derivative of the Moore determinant", okc)

    # det(A + t v v*) = det A + t v* C(A) v and v* C(A) v = det(A | v-perp)
    okr = True
    for m in (2, 3, 4):
        for _ in range(20):
            A = random_hermitian(m, rng, pos=True)
            v = [tuple(rng.normal(size=4)) for _ in range(m)]
            nv = np.sqrt(sum(sum(x * x for x in q) for q in v))
            v = [qscale(1 / nv, q) for q in v]
            t = rng.normal()
            vv = [[qmul(v[a], qconj(v[b])) for b in range(m)] for a in range(m)]
            At = [[qadd(A[a][b], qscale(t, vv[a][b])) for b in range(m)] for a in range(m)]
            M = complex_rep(A)
            cof = moore_det_numeric(A) * np.linalg.inv(M)
            vc = np.zeros(2 * m, dtype=complex)
            for a in range(m):
                vc[2 * a] = complex(v[a][0], v[a][1])
                vc[2 * a + 1] = complex(v[a][2], v[a][3])
            # v as a column: q = z + w j -> rows (z, w); the quadratic form is
            # read off the first column of the 2x2 blocks.
            col = np.zeros((2 * m, 2), dtype=complex)
            for a in range(m):
                z = complex(v[a][0], v[a][1])
                w = complex(v[a][2], v[a][3])
                col[2 * a, 0] = z
                col[2 * a + 1, 0] = -np.conj(w)
                col[2 * a, 1] = w
                col[2 * a + 1, 1] = np.conj(z)
            vCv = (col.conj().T @ cof @ col)[0, 0].real
            lhs = moore_det_numeric(At)
            okr &= abs(lhs - (moore_det_numeric(A) + t * vCv)) < 1e-7 * (1 + abs(lhs))
            # compression to the orthogonal complement of v: complete v to an
            # orthonormal basis with the complex representation
            Q, _ = np.linalg.qr(np.concatenate([col, np.eye(2 * m)[:, : 2 * m - 2]], axis=1))
            rest = Q[:, 2:]
            comp = rest.conj().T @ M @ rest
            ev = np.sort(np.linalg.eigvalsh(comp))
            okr &= abs(np.prod(ev[::2]) - vCv) < 1e-7 * (1 + abs(vCv))
    report("rank-one update and compression formula for the cofactor", okr)

    # Schur complement: det M = det(D) sigma and M - sigma e0 e0* >= 0
    oks = True
    for m in (2, 3, 4):
        for _ in range(20):
            A = random_hermitian(m, rng, pos=True)
            A[0][0] = (rng.normal() * 5, 0, 0, 0)    # time entry of either sign
            M = complex_rep(A)
            D = M[2:, 2:]
            c = M[2:, :2]
            ell = M[:2, :2]
            S = ell - c.conj().T @ np.linalg.inv(D) @ c
            sigma = S[0, 0].real
            evD = np.sort(np.linalg.eigvalsh(D))
            detD = np.prod(evD[::2])
            oks &= abs(moore_det_numeric(A) - detD * sigma) < 1e-6 * (1 + abs(detD * sigma))
            Ms = M.copy()
            Ms[0, 0] -= sigma
            Ms[1, 1] -= sigma
            oks &= np.linalg.eigvalsh(Ms).min() > -1e-8
    report("Schur complement: det M = (det D) sigma, M - sigma e0 e0* >= 0", oks)


# ---------------------------------------------------------------------------
# 2. Quaternionic Hessian calculus and the cofactor divergence (symbolic).
# ---------------------------------------------------------------------------

def hessian_q(u, coords, slots):
    """coords[a][alpha] are sympy symbols; slots = number of real slots (3 or 4)."""
    m = len(coords)
    H = [[None] * m for _ in range(m)]
    for a in range(m):
        for b in range(m):
            s = (0, 0, 0, 0)
            for al in range(slots):
                for be in range(slots):
                    d2 = sp.diff(u, coords[a][al], coords[b][be])
                    e = qmul(UNITS[al], qconj(UNITS[be]))
                    s = qadd(s, qscale(d2, e))
            H[a][b] = tuple(sp.expand(x) for x in s)
    return H


def gradient_q(u, coords, slots):
    return [tuple(sp.expand(x) for x in
                  [sum((qscale(sp.diff(u, coords[a][al]), UNITS[al])[t] for al in range(slots)), 0)
                   for t in range(4)])
            for a in range(len(coords))]


def check_hessian_calculus():
    slots = 3
    X = [[sp.Symbol(f"z{a}{al}", real=True) for al in range(slots)] for a in range(2)]
    allv = [x for row in X for x in row]
    rng = random.Random(3)

    def rand_poly(deg):
        mons = [sp.Mul(*[v ** e for v, e in zip(allv, exps)])
                for exps in itertools.product(range(deg + 1), repeat=len(allv))
                if sum(exps) <= deg]
        return sum(rng.randint(-3, 3) * mm for mm in rng.sample(mons, 12))

    u, v = rand_poly(3), rand_poly(3)
    Hu, Hv = hessian_q(u, X, slots), hessian_q(v, X, slots)
    Huv = hessian_q(sp.expand(u * v), X, slots)
    Gu, Gv = gradient_q(u, X, slots), gradient_q(v, X, slots)
    ok = True
    for a in range(2):
        for b in range(2):
            rhs = qadd(qadd(qscale(u, Hv[a][b]), qscale(v, Hu[a][b])),
                       qadd(qmul(Gu[a], qconj(Gv[b])), qmul(Gv[a], qconj(Gu[b]))))
            ok &= all(sp.expand(Huv[a][b][t] - rhs[t]) == 0 for t in range(4))
            ok &= all(sp.expand(Hu[a][b][t] - qconj(Hu[b][a])[t]) == 0 for t in range(4))
    tr = sum(Hu[a][a][0] for a in range(2))
    lap = sum(sp.diff(u, x, x) for x in allv)
    ok &= sp.expand(tr - lap) == 0
    report("Hessian calculus: Hermitian, trace = Laplacian, product rule", ok)

    # Hessian of |F|^2 for holomorphic F in the variables P_a = p_a + k q_a
    # (slot i carries p_a, slot j carries q_a): H(|F|^2) = 4 J* J.
    p = [X[a][1] for a in range(2)]
    q = [X[a][2] for a in range(2)]
    # holomorphic polynomials in P with C_k coefficients, C_k identified with C
    P = [p[a] + sp.I * q[a] for a in range(2)]
    F = [P[0] ** 2 + 2 * P[0] * P[1] + 3, (1 + 2 * sp.I) * P[1] ** 3 - P[0]]
    absF2 = sp.expand(sum(sp.expand(f * sp.conjugate(f)) for f in F))
    absF2 = sp.expand(absF2.subs({sp.conjugate(x): x for x in allv}))
    H = hessian_q(absF2, X, slots)
    ok2 = True
    for a in range(2):
        for b in range(2):
            # 4 sum_nu conj(dF_nu/dP_a) dF_nu/dP_b, a C_k number written with k
            # compute complex derivative via d/dp (holomorphic: d/dP = d/dp)
            val = sum(4 * sp.conjugate(sp.diff(f, p[a])) * sp.diff(f, p[b]) for f in F)
            val = sp.expand(val.subs({sp.conjugate(x): x for x in allv}))
            re, im = sp.expand(sp.re(val)), sp.expand(sp.im(val))
            target = (re, 0, 0, im)          # C_k number re + im k
            ok2 &= all(sp.simplify(H[a][b][t] - target[t]) == 0 for t in range(4))
    report("Hessian of |F|^2 equals 4 J* J for holomorphic F", ok2)

    # cofactor divergence: sum_{a,alpha} conj(e_alpha) d_{a alpha} C(H(u))_{ab} = 0
    def cof2(H):
        return [[H[1][1], qscale(-1, H[0][1])], [qscale(-1, H[1][0]), H[0][0]]]

    ok3 = True
    for _ in range(3):
        u = rand_poly(4)
        H = hessian_q(u, X, slots)
        C = cof2(H)
        for b in range(2):
            s = (0, 0, 0, 0)
            for a in range(2):
                for al in range(slots):
                    dC = tuple(sp.diff(C[a][b][t], X[a][al]) for t in range(4))
                    s = qadd(s, qmul(qconj(UNITS[al]), dC))
            ok3 &= all(sp.expand(x) == 0 for x in s)
    report("cofactor of the Hessian has divergence-free rows (m = 2, three slots)", ok3)

    # m = 3: the cofactor of the 3x3 Moore determinant, three slots
    X3 = [[sp.Symbol(f"y{a}{al}", real=True) for al in range(slots)] for a in range(3)]
    allv3 = [x for row in X3 for x in row]
    rng3 = random.Random(5)
    mons = [sp.Mul(*[v ** e for v, e in zip(allv3, exps)])
            for exps in itertools.product(range(3), repeat=len(allv3)) if sum(exps) <= 3]
    u3 = sum(rng3.randint(-2, 2) * mm for mm in rng3.sample(mons, 10))
    H3 = hessian_q(u3, X3, slots)

    def det3(A):
        a, b, c = A[0][0][0], A[1][1][0], A[2][2][0]
        x, y, z = A[0][1], A[0][2], A[1][2]
        n2 = lambda q: sum(t * t for t in q)
        return a * b * c - a * n2(z) - b * n2(y) - c * n2(x) + 2 * qmul(qmul(x, z), qconj(y))[0]

    # cofactor by polarization: C(A):B = d/dt det(A + tB)
    tt = sp.Symbol('tt')
    C3 = [[None] * 3 for _ in range(3)]
    for a in range(3):
        for b in range(3):
            comps = []
            for t in range(4):
                if a == b and t > 0:
                    comps.append(0)
                    continue
                if a == b:
                    E = [[ZERO] * 3 for _ in range(3)]
                    E[a][a] = (1, 0, 0, 0)
                    Bm = E
                    val = sp.diff(det3([[qadd(H3[i][j], qscale(tt, Bm[i][j])) for j in range(3)]
                                        for i in range(3)]), tt).subs(tt, 0)
                    comps.append(sp.expand(val))
                else:
                    E = [[ZERO] * 3 for _ in range(3)]
                    E[a][b] = UNITS[t]
                    E[b][a] = qconj(UNITS[t])
                    val = sp.diff(det3([[qadd(H3[i][j], qscale(tt, E[i][j])) for j in range(3)]
                                        for i in range(3)]), tt).subs(tt, 0)
                    # d det(A + tE)/dt = Re tr(C E) = 2 <C_ab, e_t>
                    comps.append(sp.expand(val / 2))
            C3[a][b] = tuple(comps)
    ok4 = True
    for b in range(3):
        s = (0, 0, 0, 0)
        for a in range(3):
            for al in range(slots):
                dC = tuple(sp.diff(C3[a][b][t], X3[a][al]) for t in range(4))
                s = qadd(s, qmul(qconj(UNITS[al]), dC))
        ok4 &= all(sp.expand(x) == 0 for x in s)
    report("cofactor of the Hessian has divergence-free rows (m = 3, three slots)", ok4)


# ---------------------------------------------------------------------------
# 3. Identities of the wave current (null frame).
# ---------------------------------------------------------------------------

def check_wave_current():
    d = 3
    kappa = sp.Symbol('kappa', positive=True)
    t, s0 = sp.symbols('t s0', real=True)
    x = sp.symbols('x1:4', real=True)
    v = sp.symbols('v1:4', real=True)
    vv = sum(c * c for c in v)
    u = sp.exp(sp.I * kappa * (s0 + sum(a * b for a, b in zip(v, x)) - vv * t / 2))
    wave = 2 * sp.diff(u, t, s0) + sum(sp.diff(u, c, c) for c in x)
    report("plane mode solves 2 d_t d_s0 + Laplacian = 0", sp.simplify(wave / u) == 0)

    # physical frequency: T = (t - s0)/sqrt2, b = (t + s0)/sqrt2
    xib = kappa * (1 - vv / 2) / sp.sqrt(2)
    norm2 = kappa ** 2 * vv + xib ** 2
    ok = sp.simplify(norm2 - (kappa * (1 + vv / 2) / sp.sqrt(2)) ** 2) == 0
    ok &= sp.simplify(xib + kappa * (1 + vv / 2) / sp.sqrt(2) - sp.sqrt(2) * kappa) == 0
    report("frequency of a mode: |xi| = kappa(1+|v|^2/2)/sqrt2, xi_b + |xi| = sqrt2 kappa", ok)

    # symbol of the system and the Gaussian
    n = d + 1
    E = lambda a, b: sp.Matrix(n, n, lambda i, j: 1 if (i, j) == (a, b) else 0)
    Ct = 2 * E(0, 0)
    Cx = [E(0, j) + E(j, 0) for j in range(1, n)]
    Cs = -sum((E(j, j) for j in range(1, n)), sp.zeros(n, n))
    pcov = [kappa * (-vv / 2)] + [kappa * c for c in v] + [kappa]
    Cp = pcov[0] * Ct + sum((pcov[1 + j] * Cx[j] for j in range(d)), sp.zeros(n, n)) + pcov[-1] * Cs
    target = -kappa * sp.Matrix(n, n, lambda i, j:
                                 vv if (i, j) == (0, 0) else
                                 (-v[j - 1] if i == 0 and j > 0 else
                                  (-v[i - 1] if j == 0 and i > 0 else (1 if i == j else 0))))
    ok = sp.simplify(Cp - target) == sp.zeros(n, n)
    w = sp.Matrix(sp.symbols('w0:4', real=True))
    quad = sp.expand(-(w.T * Cp * w)[0])
    ok &= sp.expand(quad - kappa * sum((w[j] - v[j - 1] * w[0]) ** 2 for j in range(1, n))) == 0
    Ev = sp.Matrix([1] + list(v))
    ok &= sp.simplify(Cp * Ev) == sp.zeros(n, 1)
    report("symbol of the wave system, its kernel E(v) and the quadratic form", ok)

    # holomorphic flux map: Theta(R, zeta, eta) = (R^2, R zeta, eta - zeta.zeta/2)
    R = sp.Symbol('R')
    zeta = sp.symbols('zeta1:4')
    eta = sp.Symbol('eta')
    Theta = [R ** 2] + [R * z for z in zeta] + [eta - sum(z * z for z in zeta) / 2]
    Jm = sp.Matrix([[sp.diff(f, var) for var in [R, *zeta, eta]] for f in Theta])
    report("complex Jacobian of the flux map is 2 R^{d+1}",
           sp.simplify(Jm.det() - 2 * R ** (d + 1)) == 0)

    # M E(v) = 2 R V^v and the two differentials of Phi and W
    z = sp.Matrix([R] + list(zeta))
    Ms = [(Ct * z), *[(Cx[j] * z) for j in range(d)], (Cs * z)]
    MEv = sp.Matrix([(m.T * Ev)[0] for m in Ms])
    W = [zz / R for zz in zeta]
    Vv = sp.Matrix([1] + [(v[j] + W[j]) / 2 for j in range(d)] + [-sum(a * b for a, b in zip(W, v)) / 2])
    ok = sp.simplify(MEv - 2 * R * Vv) == sp.zeros(n + 1, 1)
    Y, S1 = sp.symbols('Y S1')
    Pv = sp.symbols('P1:4')
    Phi = S1 + sum(c * c for c in Pv) / (2 * Y)
    Wy = [c / Y for c in Pv]
    dirn = [Vv[0]] + [Vv[1 + j] for j in range(d)] + [Vv[-1]]
    sub = {Y: R ** 2, **{Pv[j]: R * zeta[j] for j in range(d)}}
    dPhi = sum(sp.diff(Phi, var) * c for var, c in zip([Y, *Pv, S1], dirn))
    ok &= sp.simplify(dPhi.subs(sub)) == 0
    for j in range(d):
        dW = sum(sp.diff(Wy[j], var) * c for var, c in zip([Y, *Pv, S1], dirn))
        ok &= sp.simplify(dW.subs(sub) - (v[j] - W[j]) / (2 * R ** 2)) == 0
    report("M E(v) = 2 R V^v; (dPhi) V^v = 0; (dW) V^v = (v - W)/(2Y)", ok)

    # l(w) V^v = -(w - W).(w - v)/2
    wv = sp.symbols('u1:4')
    lw = sp.Matrix([-sum(c * c for c in wv) / 2] + list(wv) + [1])
    lhs = (lw.T * Vv)[0]
    rhs = -sum((wv[j] - W[j]) * (wv[j] - v[j]) for j in range(d)) / 2
    report("l(w) V^v = -(w - W).(w - v)/2", sp.simplify(lhs - rhs) == 0)

    # Gaussian weight: |zeta - v R|^2 = |R|^2 |W - v|^2 for real data
    Rr = sp.Symbol('Rr', positive=True)
    zr = sp.symbols('zr1:4', real=True)
    lhs = sum((zr[j] - v[j] * Rr) ** 2 for j in range(d))
    rhs = Rr ** 2 * sum((zr[j] / Rr - v[j]) ** 2 for j in range(d))
    report("Gaussian weight in the flux coordinates", sp.simplify(lhs - rhs) == 0)


# ---------------------------------------------------------------------------
# 4. Energy comparison, power counting and the radial comparison.
# ---------------------------------------------------------------------------

def check_energy_and_exponents():
    th = sp.Symbol('theta', positive=True)
    n = sp.Symbol('n', positive=True, integer=True)
    ok = True
    for nn in range(2, 9):
        a1 = sp.integrate(th ** nn, (th, 0, 1))
        a2 = sp.integrate((1 - th) * th ** (nn - 1), (th, 0, 1))
        ok &= a1 == sp.Rational(1, nn + 1) and a2 == sp.Rational(1, nn * (nn + 1))
        # F = X0 + E, E = (1/lam) int u (theta^n detDu + (1-theta)theta^(n-1) C:B)
        # with  lam X0 = int (-u) det Du  and  I = int (-u) C:B  gives
        # F = n/(n+1) X0 - I/(lam n (n+1)); and n lam X0 - I = int G* C G.
        X0, I, lam = sp.symbols('X0 I lam', positive=True)
        F = X0 - a1 * X0 - a2 * I / lam
        G = n * lam * X0 - I
        ok &= sp.simplify(F.subs(n, nn) - (sp.Rational(nn - 1, nn + 1) * X0
                                           + G.subs(n, nn) / (lam * nn * (nn + 1)))) == 0
    report("energy comparison: F = (n-1)/(n+1) X0 + (1/(lam n(n+1))) int G* C G", ok)

    # power counting in the interior lower bound: h^2 s^{2d} (h s)^{-d q} = sigma^{n-2},
    # with q = n/d, sigma = s/h, d = n - 1
    h, s = sp.symbols('h s', positive=True)
    ok = True
    for nn in range(2, 9):
        dd = nn - 1
        q = sp.Rational(nn, dd)
        expr = h ** 2 * s ** (2 * dd) * (h * s) ** (-dd * q)
        ok &= sp.simplify(expr - (s / h) ** (nn - 2)) == 0
        ok &= 2 * q == sp.Rational(2 * nn, nn - 1)
    report("power counting: h^2 s^(2d) (hs)^(-dq) = (s/h)^(n-2), 2q = 2n/(n-1)", ok)

    # quaternionic determinant of a rank-one background plus alpha I
    ok = True
    for nn in range(2, 7):
        beta, alpha = sp.symbols('beta alpha', positive=True)
        D = sp.diag(beta + alpha, *([alpha] * (nn - 1)))
        ok &= sp.simplify(D.det() - (beta + alpha) * alpha ** (nn - 1)) == 0
    report("rank-one background: det(B + alpha I) = (beta + alpha) alpha^(n-1)", ok)

    # Remark on sharpness of the radial comparison: E = (1-cos T) p + L - g
    T, sv, b, D = sp.symbols('T s b D', positive=True)
    p = sp.log(3 / sv) / D
    g = sv ** 2 * (1 + sp.log(3 / sv)) / (4 * D)
    lap = lambda f: sp.diff(f, sv, 2) + sp.diff(f, sv) / sv
    ok = sp.simplify(lap(p)) == 0 and sp.simplify(lap(g) - p) == 0
    Lfun = p * g.subs(sv, b) + (1 - p) * g.subs(sv, 3)
    E = (1 - sp.cos(T)) * p + Lfun - g
    lapE = sp.diff(E, T, 2) + lap(E)
    ok &= sp.simplify(lapE - (sp.cos(T) - 1) * p) == 0
    ok &= sp.simplify(g.subs(sv, 3) - sp.Rational(9, 4) / D) == 0
    report("sharpness example for the radial comparison", ok)

    # Gronwall representation: psi(x)/x = 1 + int z(x-z)/x P psi/z, psi'' = P psi
    xs, zs = sp.symbols('x z', positive=True)
    Pc = sp.Rational(1, 3)
    psi = sp.sinh(sp.sqrt(Pc) * xs) / sp.sqrt(Pc)
    rhs = 1 + sp.integrate(zs * (xs - zs) / xs * Pc * psi.subs(xs, zs) / zs, (zs, 0, xs))
    ok = sp.simplify(psi / xs - rhs) == 0
    report("integral identity behind the radial comparison", ok)


def check_numerology():
    """Exponent identities for the epsilon-removal and the sparse bounds."""
    ok = True
    for d in range(2, 30):
        for pn in range(1, 60):
            # p ranges over (1, 2d/(d+1)) on a grid
            p = 1 + Fraction(pn, 60) * (Fraction(2 * d, d + 1) - 1)
            if p >= Fraction(2 * d, d + 1):
                continue
            A = Fraction(d + 3, 2 * (d + 1))
            x = Fraction(d + 1, 2 * d)            # 1/p0 at the limit p0 = 2d/(d+1)
            inv_rstar = (x * (A - 1 / p) + Fraction(1, 2) * (1 / p - x)) / (A - x)
            ok &= inv_rstar == Fraction(d + 1, d - 1) * (1 - 1 / p)
            q_opt = Fraction((d - 1)) * p / (d + 1 - 2 * p)
            ok &= 1 - 1 / q_opt == inv_rstar if p > Fraction(2 * (d + 1), d + 3) else True
            # weighted exponent: sigma = (q'/p)' at q = q_opt
            qp = q_opt / (q_opt - 1)
            s = qp / p
            sigma = s / (s - 1)
            ok &= 1 - 1 / sigma == Fraction(d + 1, d - 1) * (p - 1)
            # critical index lambda(p) equals Bochner-Riesz order max(0, d|1/p-1/2|-1/2)
            lam = d * (1 / p - Fraction(1, 2)) - Fraction(1, 2)
            ok &= lam > 0
    report("limit of r_*(p, p0, p0) as p0 -> 2d/(d+1); q_opt; weighted exponent", ok)

    # monotonicity of q_* in p0 at r0 = p0 (the limit is approached from above)
    ok = True
    for d in range(2, 15):
        for pn in range(1, 20):
            lo = Fraction(2 * (d + 1), d + 3)
            pc = Fraction(2 * d, d + 1)
            p = lo + Fraction(pn, 21) * (pc - lo)
            prev = None
            for k in range(1, 40):
                p0 = p + Fraction(k, 40) * (pc - p)
                A = Fraction(d + 3, 2 * (d + 1))
                x = 1 / p0
                inv_r = (x * (A - 1 / p) + Fraction(1, 2) * (1 / p - x)) / (A - x)
                qstar = 1 / (1 - inv_r)
                if prev is not None:
                    ok &= qstar < prev
                prev = qstar
            q_opt = Fraction(d - 1) * p / (d + 1 - 2 * p)
            ok &= prev > q_opt
    report("q_*(p, p0, p0) decreases to q_opt as p0 increases", ok)


def main():
    check_moore()
    check_hessian_calculus()
    check_wave_current()
    check_energy_and_exponents()
    check_numerology()
    print(f"{PASSED} checks passed, {FAILED} failed")


if __name__ == "__main__":
    main()
