"""Complete Projective Quadrature Engine over Q(i, sqrt(d-1)) and Riemann Sphere P^1(C).

Implements:
  1. Indefinite projective integration (`p_integrate(expr, x)`) via Ostrowski-Hermite
     reduction + Lagrange-Rothstein-Trager residues over Q(i, sqrt(d-1)).
  2. Definite and improper integration (`p_integrate(expr, (x, a, b))`) on R and C,
     including a, b in C and improper endpoints (-oo, oo) and [a, oo) on P^1(C).
  3. Complex contour integration (`p_contour_integrate(expr, z, center, radius)`) on P^1(C).
  4. Fresnel-Airy-Gaussian polynomial phase integration (`int sin(x^2) dx`, `int_0^oo sin(x^2) dx`, etc.)
     via biquadratic zeta_8 = (1+i)/sqrt(2) contour rotation in Q(i, sqrt(2)) and pure-rational
     projective evaluation on R and C.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
import math
from typing import Any

import sympy as sp

from guiders_math.eval_rc import (
    ComplexDecimal,
    eval_to_complex_decimal,
    projective_atan_R,
    projective_cos_sin_C,
    projective_exp_C,
    projective_log_C,
    projective_pi,
    rational_sqrt_decimal,
)
from guiders_math.field import P1Point, QQuadExt
from guiders_math.functor import ProjectiveLift, lift_integrand_to_projective
from guiders_math.poly_engine import (
    HermiteDecomposition,
    PoleResidue,
    eval_sympy_alg_to_complex_decimal,
    hermite_reduction,
    lagrange_residues,
)


def _eval_poly_rational_at(expr: sp.Expr, var: sp.Symbol, z_val: ComplexDecimal, digits: int) -> ComplexDecimal:
    """Evaluate a rational function expr(var) at z_val in C using pure-rational arithmetic."""
    work_prec = digits + 20
    expr = sp.cancel(expr)
    if expr == 0:
        return ComplexDecimal(Decimal(0), Decimal(0), digits=digits)
    num, den = sp.fraction(expr)

    def _eval_p(p_expr: sp.Expr) -> ComplexDecimal:
        if not p_expr.has(var):
            return eval_sympy_alg_to_complex_decimal(p_expr, digits=work_prec)
        poly = sp.Poly(p_expr, var)
        acc = ComplexDecimal(Decimal(0), Decimal(0), digits=work_prec)
        for coeff in poly.all_coeffs():
            c_val = eval_sympy_alg_to_complex_decimal(coeff, digits=work_prec)
            acc = acc * z_val + c_val
        return acc

    num_v = _eval_p(num)
    den_v = _eval_p(den)
    res = num_v / den_v
    with localcontext() as ctx:
        ctx.prec = digits
        return ComplexDecimal(+res.re, +res.im, digits=digits)


def _continuous_segment_log_ratio(
    z_a: ComplexDecimal,
    z_b: ComplexDecimal,
    pole: ComplexDecimal,
    digits: int = 50,
) -> ComplexDecimal:
    """Compute int_{z_a}^{z_b} dz / (z - pole) along the straight complex segment [z_a, z_b] in C.

    Eliminates branch-cut jumps by tracking the continuous phase rotation (winding angle)
    of (z(s) - pole) from s = 0 to s = 1 using projective cross-ratio steps.
    """
    work_prec = digits + 25
    with localcontext() as ctx:
        ctx.prec = work_prec
        da = z_a - pole
        db = z_b - pole
        if (da.re == 0 and da.im == 0) or (db.re == 0 and db.im == 0):
            raise ZeroDivisionError(f"Integration endpoint coincides with pole {pole}.")

        # Subdivide segment [z_a, z_b] if the pole is close to the segment so each step changes phase by < pi/2
        seg_vec = z_b - z_a
        seg_len2 = seg_vec.re * seg_vec.re + seg_vec.im * seg_vec.im
        min_d2 = min(da.re * da.re + da.im * da.im, db.re * db.re + db.im * db.im)
        if seg_len2 > 0:
            # Projection of pole onto line z_a + t*(z_b - z_a)
            t_proj = -((da.re * seg_vec.re + da.im * seg_vec.im) / seg_len2)
            if Decimal(0) < t_proj < Decimal(1):
                closest_re = da.re + t_proj * seg_vec.re
                closest_im = da.im + t_proj * seg_vec.im
                d2_line = closest_re * closest_re + closest_im * closest_im
                if d2_line < min_d2:
                    min_d2 = d2_line

        if min_d2 == 0:
            raise ZeroDivisionError(f"Integration path [z_a, z_b] passes directly through pole {pole}.")

        ratio_est = float(seg_len2 / min_d2) if min_d2 > 0 else 1.0
        n_sub = max(1, min(64, int(math.ceil(math.sqrt(ratio_est) * 2.0))))

        total_log = ComplexDecimal(Decimal(0), Decimal(0), digits=work_prec)
        prev_w = da
        for k in range(1, n_sub + 1):
            frac_k = Decimal(k) / Decimal(n_sub)
            cur_z = ComplexDecimal(
                z_a.re + frac_k * seg_vec.re,
                z_a.im + frac_k * seg_vec.im,
                digits=work_prec,
            )
            cur_w = cur_z - pole
            step_ratio = cur_w / prev_w
            total_log = total_log + projective_log_C(step_ratio, digits=work_prec)
            prev_w = cur_w

    with localcontext() as ctx:
        ctx.prec = digits
        return ComplexDecimal(+total_log.re, +total_log.im, digits=digits)


@dataclass(frozen=True)
class ProjectiveAntiderivative:
    """Exact projective antiderivative on P^1(Q(i, sqrt(d-1))):

        F(t) = poly_antideriv(t) + hermite_rat(t) + sum_{k=1}^m c_k * Log_P(t - alpha_k).
    """

    lift: ProjectiveLift
    hermite: HermiteDecomposition
    residues: list[PoleResidue]
    d: int = 3

    @property
    def t_var(self) -> sp.Symbol:
        return self.hermite.var

    @property
    def orig_var(self) -> sp.Symbol:
        return self.lift.orig_var

    def to_projective_sympy(self) -> sp.Expr:
        """Return closed-form symbolic antiderivative over Q(i, sqrt(d-1)) using Log_P terms."""
        t = self.t_var
        expr = self.hermite.poly_antideriv + self.hermite.hermite_rat
        for pr in self.residues:
            if pr.pole_sym is not None and pr.residue_sym is not None:
                expr = expr + pr.residue_sym * sp.log(t - pr.pole_sym)
        if self.lift.forward_sym is not None and self.lift.forward_sym != t:
            expr = expr.subs(t, self.lift.forward_sym)
        return expr

    def to_real_sympy(self) -> sp.Expr:
        """Pair conjugate poles alpha_k, conj(alpha_k) in Q(i, sqrt(d-1)) to produce a manifest
        real expression in terms of ln(...) and atan(...).
        """
        t = self.t_var
        expr = self.hermite.poly_antideriv + self.hermite.hermite_rat
        used = [False] * len(self.residues)

        for i, pr in enumerate(self.residues):
            if used[i]:
                continue
            if pr.pole_sym is None or pr.residue_sym is None:
                continue
            p_re, p_im = sp.expand(pr.pole_sym).as_real_imag()
            c_re, c_im = sp.expand(pr.residue_sym).as_real_imag()
            if p_im == 0:
                used[i] = True
                expr = expr + pr.residue_sym * sp.log(t - p_re)
                continue

            # Look for complex conjugate partner pole p_re - I*p_im
            partner_idx = None
            for j in range(i + 1, len(self.residues)):
                if used[j] or self.residues[j].pole_sym is None:
                    continue
                q_re, q_im = sp.expand(self.residues[j].pole_sym).as_real_imag()
                if sp.simplify(p_re - q_re) == 0 and sp.simplify(p_im + q_im) == 0:
                    partner_idx = j
                    break

            if partner_idx is not None and p_im > 0:
                used[i] = True
                used[partner_idx] = True
                # Pair: c * log(t - (u + i*v)) + conj(c) * log(t - (u - i*v))
                # = c_re * ln((t - u)^2 + v^2) - 2 * c_im * atan((t - u) / v)
                q_res = self.residues[partner_idx].residue_sym
                d_re, d_im = sp.expand(q_res).as_real_imag()
                if sp.simplify(c_re - d_re) == 0 and sp.simplify(c_im + d_im) == 0:
                    quad_arg = sp.expand((t - p_re) ** 2 + p_im**2)
                    atan_arg = sp.simplify((t - p_re) / p_im)
                    if c_re != 0:
                        expr = expr + sp.simplify(c_re) * sp.log(quad_arg)
                    if c_im != 0:
                        expr = expr - sp.simplify(2 * c_im) * sp.atan(atan_arg)
                    continue

            used[i] = True
            expr = expr + pr.residue_sym * sp.log(t - pr.pole_sym)

        if self.lift.forward_sym is not None and self.lift.forward_sym != t:
            expr = expr.subs(t, self.lift.forward_sym)
        return expr

    def to_latex(self, real_form: bool = True) -> str:
        target = self.to_real_sympy() if real_form else self.to_projective_sympy()
        return sp.latex(target)

    def eval_at_t(self, t_val: ComplexDecimal | int | Fraction | float | complex, digits: int = 50) -> ComplexDecimal:
        """Evaluate F(t_val) at a finite point t_val in R or C using pure-rational projective descent."""
        work_prec = digits + 20
        t_c = eval_to_complex_decimal(t_val, digits=work_prec)
        rat_expr = self.hermite.poly_antideriv + self.hermite.hermite_rat
        acc = _eval_poly_rational_at(rat_expr, self.t_var, t_c, digits=work_prec)

        for pr in self.residues:
            pole_c = pr.eval_pole(digits=work_prec)
            res_c = pr.eval_residue(digits=work_prec)
            diff_c = t_c - pole_c
            log_c = projective_log_C(diff_c, digits=work_prec)
            acc = acc + (res_c * log_c)

        with localcontext() as ctx:
            ctx.prec = digits
            return ComplexDecimal(+acc.re, +acc.im, digits=digits)

    def eval_interval_t(
        self,
        t_a: ComplexDecimal | str,
        t_b: ComplexDecimal | str,
        digits: int = 50,
    ) -> ComplexDecimal:
        """Evaluate definite or improper integral int_{t_a}^{t_b} R(t) dt on P^1(C).

        Handles:
          - Finite complex or real endpoints t_a, t_b in C,
          - Full projective line (-oo, +oo) on P^1:
                int_{-oo}^{+oo} R(t) dt = 2*pi*i * sum_{Im(alpha_k) > 0} c_k,
          - Semi-infinite projective ray [t_a, +oo) on P^1:
                int_{t_a}^{+oo} R(t) dt = -R_hermite(t_a) - sum_{k=1}^m c_k * Log_P(t_a - alpha_k),
          - Semi-infinite projective ray (-oo, t_b] on P^1.
        """
        work_prec = digits + 25

        # Case 1: Full projective circle (-oo, +oo) on P^1
        if t_a == "-oo" and t_b == "oo":
            if self.hermite.poly_part != 0:
                raise ValueError("Improper integral diverges on (-oo, +oo): non-zero polynomial part.")
            pi_c = ComplexDecimal(projective_pi(work_prec), Decimal(0), digits=work_prec)
            two_pi_i = ComplexDecimal(Decimal(0), Decimal(2) * pi_c.re, digits=work_prec)
            res_sum = ComplexDecimal(Decimal(0), Decimal(0), digits=work_prec)
            thr = Decimal(10) ** (-(work_prec - 6))
            for pr in self.residues:
                p_c = pr.eval_pole(digits=work_prec)
                c_c = pr.eval_residue(digits=work_prec)
                if abs(p_c.im) <= thr:
                    raise ValueError(f"Pole {p_c} lies on the real integration axis (-oo, +oo).")
                if p_c.im > 0:
                    res_sum = res_sum + c_c
            ans = two_pi_i * res_sum
            with localcontext() as ctx:
                ctx.prec = digits
                return ComplexDecimal(+ans.re, +ans.im, digits=digits)

        # Case 2: Semi-infinite ray [t_a, +oo) on P^1
        if t_b == "oo" and isinstance(t_a, ComplexDecimal):
            if self.hermite.poly_part != 0:
                raise ValueError("Improper integral diverges at +oo: non-zero polynomial part.")
            # At t = oo = [1 : 0] in P^1: R_hermite(oo) = 0 and sum(c_k) = 0, so
            # lim_{R -> +oo} [F(R) - F(t_a)] = -R_hermite(t_a) + sum_k c_k * (-Log_P(t_a - alpha_k)).
            acc = ComplexDecimal(Decimal(0), Decimal(0), digits=work_prec) - _eval_poly_rational_at(
                self.hermite.hermite_rat, self.t_var, t_a, digits=work_prec
            )
            for pr in self.residues:
                p_c = pr.eval_pole(digits=work_prec)
                c_c = pr.eval_residue(digits=work_prec)
                diff_c = t_a - p_c
                log_c = projective_log_C(diff_c, digits=work_prec)
                acc = acc - (c_c * log_c)
            with localcontext() as ctx:
                ctx.prec = digits
                return ComplexDecimal(+acc.re, +acc.im, digits=digits)

        # Case 3: Semi-infinite ray (-oo, t_b] on P^1
        if t_a == "-oo" and isinstance(t_b, ComplexDecimal):
            # int_{-oo}^{t_b} = int_{-oo}^{+oo} - int_{t_b}^{+oo}
            full_val = self.eval_interval_t("-oo", "oo", digits=work_prec)
            tail_val = self.eval_interval_t(t_b, "oo", digits=work_prec)
            ans = full_val - tail_val
            with localcontext() as ctx:
                ctx.prec = digits
                return ComplexDecimal(+ans.re, +ans.im, digits=digits)

        # Case 4: Both endpoints finite in R or C
        assert isinstance(t_a, ComplexDecimal) and isinstance(t_b, ComplexDecimal)
        rat_expr = self.hermite.poly_antideriv + self.hermite.hermite_rat
        val_b = _eval_poly_rational_at(rat_expr, self.t_var, t_b, digits=work_prec)
        val_a = _eval_poly_rational_at(rat_expr, self.t_var, t_a, digits=work_prec)
        acc = val_b - val_a

        for pr in self.residues:
            p_c = pr.eval_pole(digits=work_prec)
            c_c = pr.eval_residue(digits=work_prec)
            seg_log = _continuous_segment_log_ratio(t_a, t_b, p_c, digits=work_prec)
            acc = acc + (c_c * seg_log)

        with localcontext() as ctx:
            ctx.prec = digits
            return ComplexDecimal(+acc.re, +acc.im, digits=digits)


def _projective_gamma_rational(s: Fraction, digits: int = 50) -> Decimal:
    """Evaluate Gamma(p/q) on R>0 via pure-rational projective reflection/lanczos-free Euler product
    and exact closed forms for half-integers and standard fractions.
    """
    work_prec = digits + 25
    with localcontext() as ctx:
        ctx.prec = work_prec
        if s <= 0:
            raise ValueError(f"Gamma requires s > 0, got {s}")
        # Shift s into [1, 2) using Gamma(s + 1) = s * Gamma(s)
        mult = Decimal(1)
        cur = s
        while cur > 1:
            cur -= 1
            mult *= Decimal(cur.numerator) / Decimal(cur.denominator)
        while cur < 1:
            mult /= Decimal(cur.numerator) / Decimal(cur.denominator)
            cur += 1

        # Exact closed form for s = 1
        if cur == 1:
            res = mult
        # Exact closed form for s = 1/2 ( wait: shifted to 1 + 1/2 or 1/2 -> cur == 1 not reached, cur == 1 is integer)
        elif s.denominator == 2:
            # All half-integers reduce to Gamma(1/2) = sqrt(pi)!
            pi_val = projective_pi(work_prec)
            sqrt_pi = rational_sqrt_decimal(pi_val, work_prec)
            # Recompute shift to 1/2
            m_half = Decimal(1)
            c_h = s
            half = Fraction(1, 2)
            while c_h > half:
                c_h -= 1
                m_half *= Decimal(c_h.numerator) / Decimal(c_h.denominator)
            res = m_half * sqrt_pi
        else:
            # General rational s in (0, 1]: Euler-Stirling asymptotic with rising factorial shift N
            # Gamma(x) = Gamma(x + N) / (x * (x+1) * ... * (x+N-1)),
            # ln Gamma(z) = (z - 1/2)*ln(z) - z + (1/2)*ln(2*pi) + sum_{k=1}^K B_{2k} / (2k*(2k-1)*z^{2k-1})
            n_shift = max(digits + 15, 40)
            x_d = Decimal(s.numerator) / Decimal(s.denominator)
            rising = Decimal(1)
            for k in range(n_shift):
                rising *= x_d + Decimal(k)
            z = x_d + Decimal(n_shift)
            pi_val = projective_pi(work_prec)
            ln_z = projective_log_C(ComplexDecimal(z, Decimal(0), work_prec), work_prec).re
            ln_2pi = projective_log_C(ComplexDecimal(Decimal(2) * pi_val, Decimal(0), work_prec), work_prec).re
            stirling = (z - Decimal("0.5")) * ln_z - z + Decimal("0.5") * ln_2pi

            # Bernoulli numbers B_{2k} via exact rational recurrence, scaled automatically with `digits`
            bernoulli_even: list[Fraction] = []
            b_all: list[Fraction] = [Fraction(1, 1)]
            max_2k = max(32, 2 * ((digits // 3) + 8))
            for m_idx in range(1, max_2k + 1):
                if m_idx > 1 and m_idx % 2 == 1:
                    b_all.append(Fraction(0, 1))
                    continue
                s_b = Fraction(0, 1)
                from math import comb

                for j in range(m_idx):
                    s_b += comb(m_idx + 1, j) * b_all[j]
                b_m = -s_b / Fraction(m_idx + 1, 1)
                b_all.append(b_m)
                if m_idx % 2 == 0:
                    bernoulli_even.append(b_m)

            inv_z = Decimal(1) / z
            inv_z2 = inv_z * inv_z
            z_pow = inv_z
            for k_idx, b_2k in enumerate(bernoulli_even, start=1):
                two_k = 2 * k_idx
                term_c = b_2k / Fraction(two_k * (two_k - 1), 1)
                term_d = (Decimal(term_c.numerator) / Decimal(term_c.denominator)) * z_pow
                stirling += term_d
                if abs(term_d) < Decimal(10) ** (-(work_prec - 2)):
                    break
                z_pow *= inv_z2

            gamma_z = projective_exp_C(ComplexDecimal(stirling, Decimal(0), work_prec), work_prec).re
            res = gamma_z / rising

    with localcontext() as ctx:
        ctx.prec = digits
        return +res


def _gauss_legendre_nodes_weights(n_pts: int, digits: int) -> list[tuple[Decimal, Decimal]]:
    """Compute n_pts Gauss-Legendre quadrature nodes and weights on [-1, 1] using pure-rational
    Newton iteration on Legendre polynomials P_n(x).
    """
    work_prec = digits + 20
    with localcontext() as ctx:
        ctx.prec = work_prec
        pi_val = projective_pi(work_prec)
        nodes_weights: list[tuple[Decimal, Decimal]] = []

        def eval_legendre(x_val: Decimal) -> tuple[Decimal, Decimal]:
            p0 = Decimal(1)
            p1 = x_val
            for k in range(2, n_pts + 1):
                pk = (Decimal(2 * k - 1) * x_val * p1 - Decimal(k - 1) * p0) / Decimal(k)
                p0, p1 = p1, pk
            # Derivative P'_n(x) = n * (x * P_n(x) - P_{n-1}(x)) / (x^2 - 1)
            dp = Decimal(n_pts) * (x_val * p1 - p0) / (x_val * x_val - Decimal(1))
            return p1, dp

        half_n = (n_pts + 1) // 2
        tol = Decimal(10) ** (-(digits + 14))
        for i in range(1, half_n + 1):
            # Initial angle theta = pi * (i - 0.25) / (n + 0.5)
            theta = pi_val * (Decimal(4 * i - 1) / Decimal(4 * n_pts + 2))
            c_init, _ = projective_cos_sin_C(ComplexDecimal(theta, Decimal(0), work_prec), 25)
            x_cur = c_init.re
            for _ in range(60):
                pn, dpn = eval_legendre(x_cur)
                dx = pn / dpn
                x_cur -= dx
                if abs(dx) <= tol:
                    break
            _, dpn = eval_legendre(x_cur)
            w_cur = Decimal(2) / ((Decimal(1) - x_cur * x_cur) * dpn * dpn)
            nodes_weights.append((-x_cur, w_cur))
            if 2 * i - 1 != n_pts:
                nodes_weights.append((x_cur, w_cur))

        nodes_weights.sort(key=lambda pair: pair[0])
        return nodes_weights


def _eval_fresnel_phase_integral(
    meta: dict[str, Any],
    a_lim: ComplexDecimal | str,
    b_lim: ComplexDecimal | str,
    digits: int = 50,
) -> tuple[sp.Expr | None, ComplexDecimal]:
    """Evaluate proper or improper Fresnel-Airy-Gaussian phase integral:
        int_{a_lim}^{b_lim} coeff * x^m * fn(a * x^n) dx
    via biquadratic contour rotation zeta_{2n} in P^1(Q(i, sqrt(d-1))) for improper [0, oo),
    or adaptive projective rational quadrature for finite intervals [a_lim, b_lim] in R or C.
    """
    work_prec = digits + 20
    fn = meta["fn"]
    coeff = meta["coeff"]
    m_pow = int(meta["m"])
    a_param = meta["a"]
    n_deg = int(meta["n"])

    coeff_c = eval_sympy_alg_to_complex_decimal(coeff, work_prec)

    # Case 1: Improper integral on [0, oo) or (-oo, oo) via exact projective contour rotation!
    if (a_lim == "-oo" and b_lim == "oo") or (
        isinstance(a_lim, ComplexDecimal) and a_lim.re == 0 and a_lim.im == 0 and b_lim == "oo"
    ):
        is_full_line = a_lim == "-oo"
        # By parity on (-oo, oo): if integrand is odd, integral is 0; if even, 2 * int_0^oo.
        if is_full_line:
            if fn == "sin":
                # sin(a * x^n) has parity (-1)^n, so x^m * sin(a * x^n) has parity (-1)^(m + n)
                if (m_pow + n_deg) % 2 == 1:
                    return sp.Integer(0), ComplexDecimal(Decimal(0), Decimal(0), digits=digits)
            else:
                if m_pow % 2 == 1:
                    return sp.Integer(0), ComplexDecimal(Decimal(0), Decimal(0), digits=digits)

        s_frac = Fraction(m_pow + 1, n_deg)
        s_sym = sp.Rational(m_pow + 1, n_deg)
        gamma_val = _projective_gamma_rational(s_frac, work_prec)
        pi_val = projective_pi(work_prec)

        if fn in ("sin", "cos"):
            # Sector rotation x = zeta_{2n} * u with zeta_{2n} = exp(i*pi/(2n)):
            # int_0^oo x^m * exp(i*a*x^n) dx = (1 / (n * |a|^s)) * Gamma(s) * exp(i * sgn(a) * pi * s / 2)
            a_c = eval_sympy_alg_to_complex_decimal(a_param, work_prec)
            if not a_c.is_real() or a_c.re == 0:
                raise ValueError("Improper Fresnel integral requires real non-zero frequency parameter a.")
            sgn_a = Decimal(1) if a_c.re > 0 else Decimal(-1)
            abs_a_sym = sp.Abs(a_param)
            angle_sym = sp.pi * s_sym / 2

            if fn == "cos":
                trig_factor_sym = sp.cos(angle_sym)
            else:
                trig_factor_sym = sp.sign(a_param) * sp.sin(angle_sym)

            sym_ans = sp.simplify(
                coeff * (sp.Integer(2) if is_full_line else sp.Integer(1))
                * sp.gamma(s_sym)
                * trig_factor_sym
                / (sp.Integer(n_deg) * (abs_a_sym ** s_sym))
            )

            with localcontext() as ctx:
                ctx.prec = work_prec
                angle_d = pi_val * Decimal(s_frac.numerator) / Decimal(2 * s_frac.denominator)
                c_ang, s_ang = projective_cos_sin_C(ComplexDecimal(angle_d, Decimal(0), work_prec), work_prec)
                trig_d = c_ang.re if fn == "cos" else sgn_a * s_ang.re
                # |a|^s = exp_P(s * Log_P(|a|))
                abs_a_d = abs(a_c.re)
                if abs_a_d == Decimal(1):
                    a_pow_s = Decimal(1)
                else:
                    ln_a = projective_log_C(ComplexDecimal(abs_a_d, Decimal(0), work_prec), work_prec).re
                    a_pow_s = projective_exp_C(
                        ComplexDecimal(ln_a * Decimal(s_frac.numerator) / Decimal(s_frac.denominator), Decimal(0), work_prec),
                        work_prec,
                    ).re
                raw_re = gamma_val * trig_d / (Decimal(n_deg) * a_pow_s)
                if is_full_line:
                    raw_re *= Decimal(2)
                num_ans = coeff_c * ComplexDecimal(raw_re, Decimal(0), digits=work_prec)

            with localcontext() as ctx:
                ctx.prec = digits
                return sym_ans, ComplexDecimal(+num_ans.re, +num_ans.im, digits=digits)

        if fn == "exp":
            # int_0^oo x^m * exp(a * x^n) dx where Re(a) < 0 or Re(a) == 0 (purely imaginary i*b -> Fresnel!)
            a_c = eval_sympy_alg_to_complex_decimal(a_param, work_prec)
            if a_c.re == 0 and a_c.im != 0:
                # Complex Fresnel integral int_0^oo x^m * exp(i * b * x^n) dx!
                b_sym = sp.simplify(a_param / sp.I)
                abs_b_sym = sp.Abs(b_sym)
                sgn_b = sp.sign(b_sym)
                angle_sym = sp.pi * s_sym / 2
                sym_ans = sp.simplify(
                    coeff
                    * (sp.Integer(2) if is_full_line else sp.Integer(1))
                    * sp.gamma(s_sym)
                    * (sp.cos(angle_sym) + sp.I * sgn_b * sp.sin(angle_sym))
                    / (sp.Integer(n_deg) * (abs_b_sym ** s_sym))
                )
                num_ans = eval_sympy_alg_to_complex_decimal(
                    sym_ans.subs(sp.gamma(sp.Rational(1, 2)), sp.sqrt(sp.pi)),
                    digits=digits,
                )
                return sym_ans, num_ans
            if a_c.re < 0:
                neg_a = sp.simplify(-a_param)
                sym_ans = sp.simplify(
                    coeff
                    * (sp.Integer(2) if is_full_line else sp.Integer(1))
                    * sp.gamma(s_sym)
                    / (sp.Integer(n_deg) * (neg_a ** s_sym))
                )
                with localcontext() as ctx:
                    ctx.prec = work_prec
                    neg_a_c = ComplexDecimal(-a_c.re, -a_c.im, digits=work_prec)
                    log_na = projective_log_C(neg_a_c, work_prec)
                    s_dec = Decimal(s_frac.numerator) / Decimal(s_frac.denominator)
                    na_pow_s = projective_exp_C(
                        ComplexDecimal(log_na.re * s_dec, log_na.im * s_dec, digits=work_prec),
                        work_prec,
                    )
                    factor_c = ComplexDecimal(
                        gamma_val * (Decimal(2) if is_full_line else Decimal(1)) / Decimal(n_deg),
                        Decimal(0),
                        digits=work_prec,
                    )
                    num_ans = (coeff_c * factor_c) / na_pow_s
                with localcontext() as ctx:
                    ctx.prec = digits
                    return sym_ans, ComplexDecimal(+num_ans.re, +num_ans.im, digits=digits)

    # Case 2: Finite endpoints [a_lim, b_lim] in R or C
    assert isinstance(a_lim, ComplexDecimal) and isinstance(b_lim, ComplexDecimal)
    a_c = eval_sympy_alg_to_complex_decimal(a_param, work_prec)

    # Evaluate int_{a_lim}^{b_lim} coeff * z^m * fn(a * z^n) dz along straight segment in C
    # using composite high-order projective Gauss-Legendre quadrature
    n_pts = max(16, int(math.ceil(digits * 0.45)) + 6)
    if n_pts % 2 == 1:
        n_pts += 1
    gl_nodes = _gauss_legendre_nodes_weights(n_pts, work_prec)

    with localcontext() as ctx:
        ctx.prec = work_prec
        dz_tot = b_lim - a_lim
        dist_est = float(abs(dz_tot.re) + abs(dz_tot.im))
        n_panels = max(2, int(math.ceil(dist_est * 2.0)))
        panel_dz = ComplexDecimal(
            dz_tot.re / Decimal(n_panels),
            dz_tot.im / Decimal(n_panels),
            digits=work_prec,
        )
        half_panel = ComplexDecimal(
            panel_dz.re / Decimal(2),
            panel_dz.im / Decimal(2),
            digits=work_prec,
        )

        total_int = ComplexDecimal(Decimal(0), Decimal(0), digits=work_prec)
        for p_idx in range(n_panels):
            z_mid = ComplexDecimal(
                a_lim.re + (Decimal(p_idx) + Decimal("0.5")) * panel_dz.re,
                a_lim.im + (Decimal(p_idx) + Decimal("0.5")) * panel_dz.im,
                digits=work_prec,
            )
            panel_acc = ComplexDecimal(Decimal(0), Decimal(0), digits=work_prec)
            for xi, wi in gl_nodes:
                z_pt = ComplexDecimal(
                    z_mid.re + xi * half_panel.re,
                    z_mid.im + xi * half_panel.im,
                    digits=work_prec,
                )
                # Compute z_pt^m
                z_m = ComplexDecimal(Decimal(1), Decimal(0), digits=work_prec)
                for _ in range(m_pow):
                    z_m = z_m * z_pt
                # Compute a * z_pt^n
                z_n = ComplexDecimal(Decimal(1), Decimal(0), digits=work_prec)
                for _ in range(n_deg):
                    z_n = z_n * z_pt
                arg_c = a_c * z_n
                if fn == "cos":
                    f_v, _ = projective_cos_sin_C(arg_c, work_prec)
                elif fn == "sin":
                    _, f_v = projective_cos_sin_C(arg_c, work_prec)
                else:
                    f_v = projective_exp_C(arg_c, work_prec)
                val_pt = z_m * f_v
                w_c = ComplexDecimal(wi, Decimal(0), digits=work_prec)
                panel_acc = panel_acc + (w_c * val_pt)
            total_int = total_int + (half_panel * panel_acc)

        final_ans = coeff_c * total_int

    with localcontext() as ctx:
        ctx.prec = digits
        return None, ComplexDecimal(+final_ans.re, +final_ans.im, digits=digits)


@dataclass(frozen=True)
class DefiniteIntegralResult:
    """Result of a definite or improper projective integral on R or C."""

    integrand: sp.Expr
    var: sp.Symbol
    lower: Any
    upper: Any
    antiderivative: ProjectiveAntiderivative | None
    exact_sympy: sp.Expr | None
    numeric_value: ComplexDecimal

    def eval_R(self, digits: int | None = None) -> Decimal:
        if not self.numeric_value.is_real():
            raise ValueError(f"Integral value {self.numeric_value} has non-zero imaginary part.")
        if digits is not None and digits != self.numeric_value.digits:
            with localcontext() as ctx:
                ctx.prec = digits
                return +self.numeric_value.re
        return self.numeric_value.re

    def eval_C(self, digits: int | None = None) -> ComplexDecimal:
        if digits is not None and digits != self.numeric_value.digits:
            with localcontext() as ctx:
                ctx.prec = digits
                return ComplexDecimal(+self.numeric_value.re, +self.numeric_value.im, digits=digits)
        return self.numeric_value


def _parse_limit(lim: Any, digits: int) -> ComplexDecimal | str:
    if lim is sp.oo or lim == "oo" or lim == "+oo" or lim == float("inf"):
        return "oo"
    if lim is -sp.oo or lim == "-oo" or lim == float("-inf"):
        return "-oo"
    if isinstance(lim, ComplexDecimal):
        return eval_to_complex_decimal(lim, digits=digits)
    if isinstance(lim, QQuadExt):
        return eval_to_complex_decimal(lim, digits=digits)
    sym_lim = sp.sympify(lim)
    if sym_lim is sp.oo:
        return "oo"
    if sym_lim is -sp.oo:
        return "-oo"
    return eval_sympy_alg_to_complex_decimal(sym_lim, digits=digits)


def p_antiderivative(
    expr: sp.Expr | str,
    var: sp.Symbol | str | None = None,
    d: int = 3,
    digits: int = 50,
) -> ProjectiveAntiderivative:
    """Construct the exact ProjectiveAntiderivative of `expr` over Q(i, sqrt(d-1))."""
    sym_expr = sp.sympify(expr)
    if var is None:
        free = sorted(sym_expr.free_symbols, key=lambda s: s.name)
        x_sym = free[0] if free else sp.Symbol("x")
    elif isinstance(var, str):
        x_sym = sp.Symbol(var)
    else:
        x_sym = var

    lift = lift_integrand_to_projective(sym_expr, x_sym)
    herm = hermite_reduction(lift.rat_num, lift.rat_den, lift.t_var)
    residues = lagrange_residues(herm.rem_num, herm.sqfree_den, lift.t_var, d=d, digits=digits)
    return ProjectiveAntiderivative(lift=lift, hermite=herm, residues=residues, d=d)


def p_integrate(
    expr: sp.Expr | str,
    var_or_limits: sp.Symbol | str | tuple[Any, Any, Any] | None = None,
    *,
    d: int = 3,
    digits: int = 50,
) -> ProjectiveAntiderivative | DefiniteIntegralResult:
    """Unified entry point for indefinite, definite, and improper projective integrals on R and C.

    Examples:
      - Indefinite:
            F = p_integrate("1 / (1 + x^4)", "x")
      - Definite on R or C:
            res = p_integrate("1 / (1 + x^4)", ("x", 0, 1), digits=50)
            res = p_integrate("1 / (1 + z^2)", ("z", "1 + I", "2 + 3*I"), digits=50)
      - Improper on P^1:
            res = p_integrate("1 / (1 + x^4)", ("x", 0, "oo"), digits=50)
            res = p_integrate("sin(x^2)", ("x", 0, "oo"), digits=50)
    """
    sym_expr = sp.sympify(expr)

    # Indefinite integral
    if not isinstance(var_or_limits, (tuple, list)):
        return p_antiderivative(sym_expr, var=var_or_limits, d=d, digits=digits)

    if len(var_or_limits) != 3:
        raise ValueError("Definite integral limits must be a 3-tuple (var, lower, upper).")

    raw_var, raw_a, raw_b = var_or_limits
    x_sym = sp.Symbol(raw_var) if isinstance(raw_var, str) else raw_var

    work_prec = digits + 20
    lim_a = _parse_limit(raw_a, work_prec)
    lim_b = _parse_limit(raw_b, work_prec)

    lift = lift_integrand_to_projective(sym_expr, x_sym)

    # Special handling for Fresnel / Airy / Gaussian phase kernels
    if lift.kind == "fresnel_phase" and lift.phase_meta is not None:
        sym_ans, num_ans = _eval_fresnel_phase_integral(lift.phase_meta, lim_a, lim_b, digits=digits)
        return DefiniteIntegralResult(
            integrand=sym_expr,
            var=x_sym,
            lower=raw_a,
            upper=raw_b,
            antiderivative=None,
            exact_sympy=sym_ans,
            numeric_value=num_ans,
        )

    antideriv = p_antiderivative(sym_expr, var=x_sym, d=d, digits=digits)

    # Special case: full periodic trigonometric cycle [0, 2*pi] or [-pi, pi]
    if lift.kind == "trig" and isinstance(lim_a, ComplexDecimal) and isinstance(lim_b, ComplexDecimal):
        pi_d = projective_pi(work_prec)
        span = lim_b - lim_a
        if span.is_real() and abs(abs(span.re) - Decimal(2) * pi_d) <= Decimal(10) ** (-(digits - 4)):
            # Full circle on S^1 corresponds to full projective line (-oo, +oo) on P^1!
            sgn = Decimal(1) if span.re > 0 else Decimal(-1)
            val_full = antideriv.eval_interval_t("-oo", "oo", digits=digits)
            val_signed = ComplexDecimal(sgn * val_full.re, sgn * val_full.im, digits=digits)
            return DefiniteIntegralResult(
                integrand=sym_expr,
                var=x_sym,
                lower=raw_a,
                upper=raw_b,
                antiderivative=antideriv,
                exact_sympy=None,
                numeric_value=val_signed,
            )

    t_a = lift.forward_eval(lim_a, work_prec) if lift.forward_eval is not None else lim_a
    t_b = lift.forward_eval(lim_b, work_prec) if lift.forward_eval is not None else lim_b
    num_val = antideriv.eval_interval_t(t_a, t_b, digits=digits)

    # Construct exact symbolic value for standard improper rational integrals when poles are symbolic
    exact_sym: sp.Expr | None = None
    if lift.kind == "rational" and lim_a == "-oo" and lim_b == "oo":
        if all(pr.pole_sym is not None and pr.residue_sym is not None for pr in antideriv.residues):
            res_acc = sp.Integer(0)
            for pr in antideriv.residues:
                p_c = pr.eval_pole(digits=20)
                if p_c.im > 0:
                    res_acc += pr.residue_sym
            exact_sym = sp.simplify(2 * sp.pi * sp.I * res_acc)
    elif (
        lift.kind == "rational"
        and isinstance(lim_a, ComplexDecimal)
        and lim_a.re == 0
        and lim_a.im == 0
        and lim_b == "oo"
    ):
        # For even rational integrand R(-t) == R(t), int_0^oo = (1/2) * int_{-oo}^{oo}!
        t_v = antideriv.t_var
        rat_t = lift.rat_num / lift.rat_den
        if sp.simplify(rat_t - rat_t.subs(t_v, -t_v)) == 0:
            if all(pr.pole_sym is not None and pr.residue_sym is not None for pr in antideriv.residues):
                res_acc = sp.Integer(0)
                for pr in antideriv.residues:
                    p_c = pr.eval_pole(digits=20)
                    if p_c.im > 0:
                        res_acc += pr.residue_sym
                exact_sym = sp.simplify(sp.pi * sp.I * res_acc)

    return DefiniteIntegralResult(
        integrand=sym_expr,
        var=x_sym,
        lower=raw_a,
        upper=raw_b,
        antiderivative=antideriv,
        exact_sympy=exact_sym,
        numeric_value=num_val,
    )


def p_contour_integrate(
    expr: sp.Expr | str,
    var: sp.Symbol | str,
    center: Any = 0,
    radius: Any = 1,
    *,
    d: int = 3,
    digits: int = 50,
) -> tuple[sp.Expr | None, ComplexDecimal]:
    """Compute closed contour integral oint_{|z - center| = radius} f(z) dz on P^1(C).

    By Cauchy-Lagrange residue theorem on the Riemann sphere, equals:
        2 * pi * i * sum_{|alpha_k - center| < radius} Res_{z = alpha_k} f(z).
    """
    sym_expr = sp.sympify(expr)
    z_sym = sp.Symbol(var) if isinstance(var, str) else var
    work_prec = digits + 20

    c_center = _parse_limit(center, work_prec)
    c_radius = _parse_limit(radius, work_prec)
    if not isinstance(c_center, ComplexDecimal) or not isinstance(c_radius, ComplexDecimal):
        raise ValueError("Contour center and radius must be finite numbers.")
    r_dec = c_radius.re
    if r_dec <= 0:
        raise ValueError("Contour radius must be positive.")

    # Use partial fraction / Hermite reduction: note that d/dz [R_hermite(z)] has ZERO residue
    # everywhere (since it has a single-valued rational antiderivative R_hermite(z))!
    # Only the simple-pole part A(z) / Q_sqfree(z) contributes to contour residues!
    if not sym_expr.is_rational_function(z_sym):
        raise ValueError("p_contour_integrate currently requires a rational function on P^1(C).")

    cancelled = sp.cancel(sym_expr)
    num, den = sp.fraction(cancelled)
    herm = hermite_reduction(num, den, z_sym)
    residues = lagrange_residues(herm.rem_num, herm.sqfree_den, z_sym, d=d, digits=work_prec)

    r2 = r_dec * r_dec
    thr = Decimal(10) ** (-(digits - 4))
    sym_sum: sp.Expr | None = sp.Integer(0)
    num_sum = ComplexDecimal(Decimal(0), Decimal(0), digits=work_prec)

    for pr in residues:
        p_c = pr.eval_pole(digits=work_prec)
        c_c = pr.eval_residue(digits=work_prec)
        diff = p_c - c_center
        dist2 = diff.re * diff.re + diff.im * diff.im
        if abs(dist2 - r2) <= thr:
            raise ValueError(f"Pole {p_c} lies directly on the contour |z - {c_center}| = {r_dec}.")
        if dist2 < r2:
            num_sum = num_sum + c_c
            if sym_sum is not None and pr.residue_sym is not None:
                sym_sum = sym_sum + pr.residue_sym
            else:
                sym_sum = None

    pi_val = projective_pi(work_prec)
    two_pi_i = ComplexDecimal(Decimal(0), Decimal(2) * pi_val, digits=work_prec)
    final_num = two_pi_i * num_sum
    final_sym = sp.simplify(2 * sp.pi * sp.I * sym_sum) if sym_sum is not None else None

    with localcontext() as ctx:
        ctx.prec = digits
        return final_sym, ComplexDecimal(+final_num.re, +final_num.im, digits=digits)
