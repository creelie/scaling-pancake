# Weak-Type Endpoint Bounds for Bochner–Riesz Means in Every Dimension

Author: Deep Bhattacharjee. LaTeX source and verification code for the paper
of that title (54 pages, amsart). Release v1.0.0.

## What the paper proves

For n ≥ 2, a > 0 and t > 0 let R^λ_{a,t} be the Fourier multiplier with symbol
(1 − |ξ|^a / t^a)^λ_+, and let λ(p) = n(1/p − 1/2) − 1/2 be the critical index.

- **Theorem A (weak-type endpoint).** For 1 < p < 2n/(n+1) the operators
  R^{λ(p)}_{a,t} map L^p(R^n) to L^{p,∞}(R^n) uniformly in t, and
  R^{λ(p)}_{a,t} f → f in L^{p,∞} as t → ∞.
- **Theorem B (sparse bounds).** For 1 ≤ p < 2n/(n+1) and
  q > (n−1)p/(n+1−2p) the operators R^{λ(p)}_{a,t} have (p,q) sparse bounds
  uniform in t. This covers the critical side of the trapezoid of Beltran,
  Roos and Seeger except its upper corner, which stays open.
- **Theorem C (weights).** If w ∈ A_1 ∩ RH_σ with σ > 1, then for
  1 < p < 1 + (n−1)/(n+1) (1 − 1/σ) the operators map L^p(w) to L^{p,∞}(w)
  uniformly in t, with convergence.

For n = 2 these were known (Seeger; Beltran, Roos and Seeger); the new content
is n ≥ 3. The theorems are consequences of OpenAI's preprint *Local smoothing
for the Euclidean wave equation* (25 October 2026), which proves Sogge's local
smoothing conjecture and with it the Bochner–Riesz conjecture in every
dimension, combined with the ε-removal of Tao in the vector-valued form of
Beltran, Roos and Seeger, their sparse domination theorem (Math. Ann. 2025)
and the weighted theorem of Frey and Nieraeth (J. Geom. Anal. 2019).

The two unrefereed inputs are reproduced with complete proofs and with credit:
Part I (Sections 2–7 and Appendix A) gives OpenAI's proof of the critical
annular estimate, of local smoothing and of the Bochner–Riesz conjecture,
including the Dirichlet problem it uses; Part II (Section 8) gives the
ε-removal of the Oberwolfach preprint of Beltran, Roos and Seeger. Part III
(Section 9) proves Theorems A–C, and Section 10 lists what remains open: the
corner of the trapezoid, its top side, the sharp range of the weights, other
convex domains and the maximal operator.

## Layout

- `paper/`: LaTeX source (`main.tex` and the section files); figures in
  `paper/figures/`. Build with `latexmk -pdf main.tex` (0 errors, no bad
  boxes, no undefined references).
- `figures/make_figures.py`: regenerates the three figures (matplotlib).
- `code/check_algebra.py`: symbolic and numerical checks of the algebraic
  identities and exponent computations (Python with sympy and numpy). Its
  recorded output is `code/check_output.txt`; the last line is
  `22 checks passed, 0 failed`.
- `lean/BochnerRiesz.lean`: eighteen theorems in Lean 4 (core library only,
  toolchain in `lean/lean-toolchain`): the exponent arithmetic of Sections 8
  and 9 and the quaternion identities of Part I. Each theorem depends only on
  the standard axioms `propext`, `Classical.choice` and `Quot.sound`
  (`lean/axioms.txt`); there is no `sorry`.
- `dist/`: the compiled paper, the TeX zip and the arXiv tarball attached to
  the release.
- `.github/workflows/check.yml` runs the checks, the Lean file and the LaTeX
  build; `.github/workflows/release.yml` publishes a release by hand.

No statement of the paper depends on the programs; every step is proved in the
text.

## Citation

See `CITATION.cff`. The Zenodo DOI of each release is listed below.

| Release | DOI |
| --- | --- |
| v1.0.0 | [10.5281/zenodo.23701803](https://doi.org/10.5281/zenodo.23701803) |

Concept DOI (all versions): [10.5281/zenodo.23701802](https://doi.org/10.5281/zenodo.23701802).

## Licence

MIT, see `LICENSE`.
