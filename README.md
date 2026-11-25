# Guiders Mathematics (`guiders-mathematics`)

**Projective Computer Algebra System over $\mathbb{Q}(i, \sqrt{d-1})$ and the Riemann Sphere $\mathbb{P}^1(\mathbb{C})$**

`guiders-mathematics` replaces the classical Liouville–Risch transcendental integration and floating-point series bottleneck with an exact algebraic-projective engine over the biquadratic field $\mathbb{Q}(i, \sqrt{d-1})$ (default $d = 3 \implies \mathbb{Q}(i, \sqrt{2})$ with eighth root of unity $\zeta_8 = \frac{1+i}{\sqrt{2}}$) and homogeneous projective line $\mathbb{P}^1(\mathbb{C})$:

1. **Exact Biquadratic Field $\mathbb{Q}(i, \sqrt{d-1})$ & Riemann Sphere $\mathbb{P}^1$ (`guiders_math.field`):**
   - Exact elements $a + b\,i + c\,\sqrt{d-1} + e\,i\sqrt{d-1}$ with $a,b,c,e \in \mathbb{Q}$.
   - Homogeneous projective points $[U : V] \in \mathbb{P}^1(\mathbb{Q}(i, \sqrt{d-1}))$ where infinity $\infty = [1 : 0]$ is a regular first-class coordinate (no `NaN`/`zoo`), stereographic circle $(\mathcal{C}([1:0]), \mathcal{S}([1:0])) = (-1, 0)$, Möbius group $\mathrm{PGL}_2(\mathbb{Q}(i, \sqrt{d-1}))$, Cayley transform, and projective cross-ratio $(z_1, z_2; z_3, z_4)$.
2. **Projective Functors $\mathcal{P}_{\mathrm{in}}, \mathcal{P}_{\mathrm{out}}$ (`guiders_math.functor`):**
   - Lifts trigonometric, hyperbolic, exponential, quadratic circle radicals, and Fresnel–Airy–Gaussian polynomial phase kernels onto $\mathbb{P}^1$.
   - Polynomial-time trigonometric and hyperbolic identity decision procedure (`p_simplify`, `p_prove_identity`).
3. **Complete Projective Quadrature (`guiders_math.quadrature`, `guiders_math.poly_engine`):**
   - **Ostrowski–Hermite reduction** extracts the rational pole part $R_{\mathrm{hermite}}(t)$ purely via polynomial GCD without root finding.
   - **Lagrange–Rothstein–Trager residues** decompose the square-free part over $\mathbb{Q}(i, \sqrt{d-1})$ (and Aberth–Ehrlich rational root solver in $\mathbb{C}$ for irreducible degree $\ge 5$ factors):
     $$
     F(t) = P_{\mathrm{poly}}(t) + R_{\mathrm{hermite}}(t) + \sum_{k=1}^{m} c_k\,\mathrm{Log}_{\mathbb{P}}(t - \alpha_k).
     $$
   - **Definite & Improper integrals on $\mathbb{R}$ and $\mathbb{C}$:** handles finite complex segments $[z_1, z_2] \subset \mathbb{C}$ without branch-cut jumps, improper integrals on $[0, \infty)$ and $(-\infty, \infty)$ via residue collapse at $\infty = [1 : 0] \in \mathbb{P}^1$, closed contour integrals $\oint_{|z-z_0|=r} f(z)\,dz$, and proper/improper/complex Fresnel integrals $\int \sin(x^2)\,dx$, $\int_0^\infty \sin(x^2)\,dx = \frac{\sqrt{2\pi}}{4}$ via $\zeta_8$ contour rotation.
4. **Pure-Rational Projective Evaluator on $\mathbb{R}$ and $\mathbb{C}$ (`guiders_math.eval_rc`):**
   - Zero `mpmath` imports and zero `math.sin`/`cos`/`atan`/`log`/`exp` calls: arbitrary-precision evaluation ($N$ decimal digits on $\mathbb{R}$ or $\mathbb{C}$) uses quadratic Heron–Newton rational square roots, complex square-root (half-angle) descent $z \mapsto \sqrt{z}$, Cayley mapping $w = \frac{z-1}{z+1}$, and the diagonal $(5, 5)$ Padé rational approximant:
     $$
     \mathrm{Log}_{\mathbb{P}}(z_m) \approx \frac{2 w_m \bigl(945 - 735 w_m^2 + 64 w_m^4\bigr)}{945 - 1050 w_m^2 + 225 w_m^4} \in \mathbb{Q}(i).
     $$
5. **Cavalieri–Lagrange Algebraic Differential & Padé Approximants (`guiders_math.calculus`):**
   - Limit-free derivative `p_diff` in the dual ring $\mathbb{Q}(i, \sqrt{d-1})(t)[\varepsilon]/(\varepsilon^2)$, exact Toeplitz Padé approximant `p_pade(expr, x, m, n)`, and projective ODE solver `p_dsolve`.

---

## Quick Start & Verification (auto-bootstraps `.venv` on a fresh clone)

```bash
python3 verify_principles.py
```

## CLI Examples

```bash
# 1. Improper rational integral on P^1: int_0^oo dx / (1 + x^4) = pi / (2*sqrt(2)) to 60 digits
python3 -m guiders_math integrate "1/(1+x^4)" --var x --from 0 --to oo --digits 60

# 2. Improper Fresnel integral via zeta_8 rotation in Q(i, sqrt(2)): int_0^oo sin(x^2) dx = sqrt(2*pi)/4
python3 -m guiders_math integrate "sin(x^2)" --var x --from 0 --to oo --digits 60

# 3. Proper Fresnel integral on [0, 1] and complex Fresnel integral on [0, 1+I]
python3 -m guiders_math integrate "sin(x^2)" --var x --from 0 --to 1 --digits 40
python3 -m guiders_math integrate "sin(z^2)" --var z --from 0 --to "1+I" --digits 40

# 4. Complex path integral in C without branch-cut jumps
python3 -m guiders_math integrate "1/(1+z^2)" --var z --from "1+I" --to "2+3*I" --digits 50

# 5. Closed contour integral on the Riemann sphere P^1(C)
python3 -m guiders_math contour "1/(1+z^4)" --var z --center "I" --radius 1 --digits 50
```
