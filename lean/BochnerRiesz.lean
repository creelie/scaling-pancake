/-
  Exact arithmetic behind the weak-type endpoint Bochner-Riesz theorem.

  Every statement is an identity or inequality between rational numbers
  (exponents, interpolation parameters and the coefficients of the energy
  comparison).  The file uses only the core library of Lean 4 and is checked
  with `lean BochnerRiesz.lean`.
-/

namespace BochnerRiesz

/-- The Bochner-Riesz index `lambda(p) = d (1/p - 1/2) - 1/2`, written with
`x = 1/p`. -/
def lam (d x : Rat) : Rat := d * (x - 1/2) - 1/2

/-- The index vanishes exactly at the critical exponent `1/p = (d+1)/(2d)`. -/
theorem lam_critical (d : Rat) (hd : d ≠ 0) : lam d ((d + 1) / (2 * d)) = 0 := by
  unfold lam
  grind

/-- At `p = 1` the index is `(d-1)/2`. -/
theorem lam_one (d : Rat) : lam d 1 = (d - 1) / 2 := by
  unfold lam
  grind

/-- The constant `A = (d+3)/(2(d+1))` of the interpolation formula and the
endpoint value `x = (d+1)/(2d)` of `1/p0` differ by `(d-1)/(2d(d+1))`, so the
denominator of the formula for `1/r_*` does not vanish at the limit when
`d >= 2`. -/
theorem A_minus_x (d : Rat) (hd : d ≠ 0) (hd1 : d + 1 ≠ 0) :
    (d + 3) / (2 * (d + 1)) - (d + 1) / (2 * d) = (d - 1) / (2 * d * (d + 1)) := by
  grind

/-- The interpolation formula returns `1/r0` at `p = p0`. -/
theorem rstar_at_p0 (A x y : Rat) (h : A - x ≠ 0) :
    (y * (A - x) + (1/2) * (x - x)) / (A - x) = y := by
  grind

/-- The limit of `1/r_*(p, p0, p0)` as `p0 -> 2d/(d+1)`. -/
theorem rstar_limit (d p : Rat) (hd : d - 1 ≠ 0) (hp : p ≠ 0) (h3 : d + 1 ≠ 0)
    (h4 : 2 * d ≠ 0) :
    let x := (d + 1) / (2 * d)
    let A := (d + 3) / (2 * (d + 1))
    (x * (A - 1/p) + (1/2) * (1/p - x)) / (A - x) = (d + 1) / (d - 1) * (1 - 1/p) := by
  intro x A
  simp only [x, A]
  grind

/-- The optimal sparse exponent `q_opt = (d-1)p/(d+1-2p)` is the dual exponent
of the limit `r_*`: `1 - 1/q_opt = (d+1)/(d-1) (1 - 1/p)`. -/
theorem qopt_dual (d p : Rat) (hd : d - 1 ≠ 0) (hp : p ≠ 0) (hq : d + 1 - 2 * p ≠ 0) :
    1 - 1 / ((d - 1) * p / (d + 1 - 2 * p)) = (d + 1) / (d - 1) * (1 - 1/p) := by
  grind

/-- The reverse Holder exponent produced by a sparse form `Sp(p, q)`:
with `s = q'/p`, `sigma = s'`, one has `1 - 1/sigma = p/q'`. -/
theorem sigma_formula (p qq : Rat) (hp : p ≠ 0) (hq : qq ≠ 0) (hs : qq - p ≠ 0) :
    let s := qq / p
    let sigma := s / (s - 1)
    1 - 1 / sigma = p / qq := by
  intro s sigma
  simp only [s, sigma]
  grind

/-- With `q = q_opt`, `p/q' = (d+1)/(d-1) (p - 1)`; this is the weighted range
`p < 1 + (d-1)/(d+1) (1 - 1/sigma)`. -/
theorem weighted_exponent (d p : Rat) (hd : d - 1 ≠ 0) (hp : p ≠ 0)
    (hq : d + 1 - 2 * p ≠ 0) (hq1 : (d - 1) * p - (d + 1 - 2 * p) ≠ 0) :
    let q := (d - 1) * p / (d + 1 - 2 * p)
    let qq := q / (q - 1)
    p / qq = (d + 1) / (d - 1) * (p - 1) := by
  intro q qq
  simp only [q, qq]
  grind

/-- Inverting the previous relation. -/
theorem weighted_inverse (d p t : Rat) (hd : d - 1 ≠ 0) (hd1 : d + 1 ≠ 0)
    (h : t = (d + 1) / (d - 1) * (p - 1)) :
    p = 1 + (d - 1) / (d + 1) * t := by
  grind

/-- The energy comparison: with `int theta^n = 1/(n+1)` and
`int (1-theta) theta^(n-1) = 1/(n(n+1))`,
`X0 - X0/(n+1) - I/(lam n (n+1)) = (n-1)/(n+1) X0 + (n lam X0 - I)/(lam n (n+1))`. -/
theorem energy_comparison (n l X I : Rat) (hn : n ≠ 0) (hn1 : n + 1 ≠ 0) (hl : l ≠ 0) :
    X - X / (n + 1) - I / (l * (n * (n + 1)))
      = (n - 1) / (n + 1) * X + (n * l * X - I) / (l * n * (n + 1)) := by
  grind

/-- Power counting in the interior lower bound: with `d = n - 1` and
`q = n/d`, the exponent of `h` in `h^2 s^(2d) (hs)^(-dq)` is `2 - n` and the
exponent of `s` is `n - 2`. -/
theorem power_counting (n : Rat) (hd : n - 1 ≠ 0) :
    2 - (n - 1) * (n / (n - 1)) = -(n - 2) ∧
    2 * (n - 1) - (n - 1) * (n / (n - 1)) = n - 2 := by
  constructor <;> grind

/-- The conversion between the multipliers `(1 - rho^2)` and `(1 - rho)`. -/
theorem multiplier_factor (rho : Rat) : 1 - rho ^ 2 = (1 - rho) * (1 + rho) := by
  grind

/-- Quaternions with rational coordinates. -/
structure Quat where
  re : Rat
  i : Rat
  j : Rat
  k : Rat

def Quat.mul (p q : Quat) : Quat where
  re := p.re * q.re - p.i * q.i - p.j * q.j - p.k * q.k
  i := p.re * q.i + p.i * q.re + p.j * q.k - p.k * q.j
  j := p.re * q.j - p.i * q.k + p.j * q.re + p.k * q.i
  k := p.re * q.k + p.i * q.j - p.j * q.i + p.k * q.re

def Quat.normSq (p : Quat) : Rat := p.re * p.re + p.i * p.i + p.j * p.j + p.k * p.k

def Quat.conj (p : Quat) : Quat := ⟨p.re, -p.i, -p.j, -p.k⟩

/-- The norm is multiplicative (Euler's four-square identity). -/
theorem normSq_mul (p q : Quat) : (p.mul q).normSq = p.normSq * q.normSq := by
  cases p; cases q
  simp only [Quat.mul, Quat.normSq]
  grind

/-- Conjugation reverses products. -/
theorem conj_mul (p q : Quat) : (p.mul q).conj = q.conj.mul p.conj := by
  cases p; cases q
  simp only [Quat.mul, Quat.conj, Quat.mk.injEq]
  refine ⟨?_, ?_, ?_, ?_⟩ <;> grind

/-- `p * conj p` is the real number `|p|^2`. -/
theorem mul_conj (p : Quat) : p.mul p.conj = ⟨p.normSq, 0, 0, 0⟩ := by
  cases p
  simp only [Quat.mul, Quat.conj, Quat.normSq, Quat.mk.injEq]
  refine ⟨?_, ?_, ?_, ?_⟩ <;> grind

/-- Multiplication is associative. -/
theorem mul_assoc' (p q r : Quat) : (p.mul q).mul r = p.mul (q.mul r) := by
  cases p; cases q; cases r
  simp only [Quat.mul, Quat.mk.injEq]
  refine ⟨?_, ?_, ?_, ?_⟩ <;> grind

/-- The Moore determinant of a Hermitian quaternionic `3 x 3` matrix with
diagonal `a, b, c` and upper entries `x, y, z`. -/
def moore3 (a b c : Rat) (x y z : Quat) : Rat :=
  a * b * c - a * z.normSq - b * y.normSq - c * x.normSq + 2 * ((x.mul z).mul y.conj).re

/-- For real off-diagonal entries the Moore determinant is the usual one. -/
theorem moore3_real (a b c x y z : Rat) :
    moore3 a b c ⟨x, 0, 0, 0⟩ ⟨y, 0, 0, 0⟩ ⟨z, 0, 0, 0⟩
      = a * (b * c - z * z) - x * (x * c - z * y) + y * (x * z - b * y) := by
  simp only [moore3, Quat.mul, Quat.conj, Quat.normSq]
  grind

/-- The rank-one background used for the Dirichlet problem:
`det(diag(beta + alpha, alpha, alpha)) = (beta + alpha) alpha^2`. -/
theorem moore3_rank_one (alpha beta : Rat) :
    moore3 (beta + alpha) alpha alpha ⟨0, 0, 0, 0⟩ ⟨0, 0, 0, 0⟩ ⟨0, 0, 0, 0⟩
      = (beta + alpha) * alpha ^ 2 := by
  simp only [moore3, Quat.mul, Quat.conj, Quat.normSq]
  grind

end BochnerRiesz
