"""Exact polynomial engine: Ostrowski-Hermite reduction and Lagrange-Rothstein-Trager residues.

Given any rational function R(t) = P(t) / Q(t) over Q(i, sqrt(d-1)):
  1. Polynomial division extracts the polynomial part P_poly(t) so that deg(R_rem) < deg(Q).
  2. Square-free factorization + Ostrowski-Hermite reduction extracts the rational pole part
     R_hermite(t) purely via polynomial GCD / Bezout identities (zero root-finding!),
     reducing the remaining integral to a square-free denominator:
         P(t) / Q(t) = P_poly(t) + d/dt [ R_hermite(t) ] + A(t) / Q_sqfree(t).
  3. Lagrange-Rothstein-Trager residue extraction decomposes A(t) / Q_sqfree(t) into
     simple poles alpha_k with residues c_k = A(alpha_k) / Q'_sqfree(alpha_k):
         A(t) / Q_sqfree(t) = sum_{k=1}^m c_k / (t - alpha_k).
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
import math

import sympy as sp

from guiders_math.eval_rc import (
    ComplexDecimal,
    _c_add,
    _c_div,
    _c_mul,
    _c_sub,
    eval_qquad_to_complex_decimal,
    eval_to_complex_decimal,
    projective_cos_sin_C,
    projective_pi,
    rational_sqrt_decimal,
)
from guiders_math.field import QQuadExt


@dataclass(frozen=True)
class PoleResidue:
    """A simple pole alpha_k and its Lagrange residue c_k = A(alpha_k) / Q'(alpha_k) on P^1(C).

    Stores:
      - `pole_sym` / `residue_sym`: exact symbolic representation (in Q(i, sqrt(d-1)) or radicals)
        when available,
      - `pole_q` / `residue_q`: exact `QQuadExt` elements when the pole lies in Q(i, sqrt(d-1)),
      - `pole_num` / `residue_num`: arbitrary-precision `ComplexDecimal` fallback/cache.
    """

    pole_sym: sp.Expr | None
    residue_sym: sp.Expr | None
    pole_q: QQuadExt | None = None
    residue_q: QQuadExt | None = None
    pole_num: ComplexDecimal | None = None
    residue_num: ComplexDecimal | None = None

    def eval_pole(self, digits: int = 50) -> ComplexDecimal:
        if self.pole_q is not None:
            return eval_qquad_to_complex_decimal(self.pole_q, digits=digits)
        if self.pole_num is not None and self.pole_num.digits >= digits:
            return eval_to_complex_decimal(self.pole_num, digits=digits)
        if self.pole_sym is not None:
            return eval_sympy_alg_to_complex_decimal(self.pole_sym, digits=digits)
        if self.pole_num is not None:
            return self.pole_num
        raise ValueError("PoleResidue has neither symbolic nor numeric pole.")

    def eval_residue(self, digits: int = 50) -> ComplexDecimal:
        if self.residue_q is not None:
            return eval_qquad_to_complex_decimal(self.residue_q, digits=digits)
        if self.residue_num is not None and self.residue_num.digits >= digits:
            return eval_to_complex_decimal(self.residue_num, digits=digits)
        if self.residue_sym is not None:
            return eval_sympy_alg_to_complex_decimal(self.residue_sym, digits=digits)
        if self.residue_num is not None:
            return self.residue_num
        raise ValueError("PoleResidue has neither symbolic nor numeric residue.")


@dataclass(frozen=True)
class HermiteDecomposition:
    """Result of Ostrowski-Hermite reduction on R(t) = P(t) / Q(t).

    Satisfies:
        int (P(t) / Q(t)) dt = poly_antideriv(t) + hermite_rat(t) + int (rem_num(t) / sqfree_den(t)) dt.
    """

    var: sp.Symbol
    poly_part: sp.Expr
    poly_antideriv: sp.Expr
    hermite_rat: sp.Expr
    rem_num: sp.Expr
    sqfree_den: sp.Expr


def eval_sympy_alg_to_complex_decimal(expr: sp.Expr, digits: int = 50) -> ComplexDecimal:
    """Evaluate an algebraic SymPy expression (built from rationals, I, +, *, /, fractional powers)
    to `ComplexDecimal` using pure projective Heron-Newton & Pade evaluation (NO mpmath/evalf).
    """
    work_prec = digits + 20
    expr = sp.sympify(expr)

    # Try fast exact extraction into Q(i, sqrt(2)) or Q(i, sqrt(3)) first
    for d_try in (3, 4, 6):
        try:
            q_elem = QQuadExt.from_sympy(expr, d=d_try)
            return eval_qquad_to_complex_decimal(q_elem, digits=digits)
        except Exception:
            pass

    def _eval_node(node: sp.Expr) -> ComplexDecimal:
        if node.is_Integer:
            return ComplexDecimal(Decimal(int(node)), Decimal(0), digits=work_prec)
        if node.is_Rational:
            with localcontext() as ctx:
                ctx.prec = work_prec
                val = Decimal(int(node.p)) / Decimal(int(node.q))
            return ComplexDecimal(val, Decimal(0), digits=work_prec)
        if node is sp.I:
            return ComplexDecimal(Decimal(0), Decimal(1), digits=work_prec)
        if node is sp.pi:
            return ComplexDecimal(projective_pi(work_prec), Decimal(0), digits=work_prec)
        if node.is_Add:
            acc = ComplexDecimal(Decimal(0), Decimal(0), digits=work_prec)
            for arg in node.args:
                acc = acc + _eval_node(arg)
            return acc
        if node.is_Mul:
            acc = ComplexDecimal(Decimal(1), Decimal(0), digits=work_prec)
            for arg in node.args:
                acc = acc * _eval_node(arg)
            return acc
        if node.is_Pow:
            base_val = _eval_node(node.base)
            exp = node.exp
            if exp.is_Integer:
                n = int(exp)
                if n == 0:
                    return ComplexDecimal(Decimal(1), Decimal(0), digits=work_prec)
                if n < 0:
                    one = ComplexDecimal(Decimal(1), Decimal(0), digits=work_prec)
                    return one / _eval_node(node.base ** (-n))
                res = ComplexDecimal(Decimal(1), Decimal(0), digits=work_prec)
                b = base_val
                while n > 0:
                    if n & 1:
                        res = res * b
                    b = b * b
                    n >>= 1
                return res
            if isinstance(exp, sp.Rational) and int(exp.q) == 2:
                # Half-integer power base^(p/2) via pure complex Heron square root
                from guiders_math.eval_rc import complex_sqrt_heron

                re_s, im_s = complex_sqrt_heron(base_val.re, base_val.im, work_prec)
                sqrt_c = ComplexDecimal(re_s, im_s, digits=work_prec)
                p = int(exp.p)
                if p == 1:
                    return sqrt_c
                if p == -1:
                    return ComplexDecimal(Decimal(1), Decimal(0), digits=work_prec) / sqrt_c
                return _eval_node(sp.Dummy() ** p) if False else (
                    sqrt_c * _eval_node(node.base ** ((p - 1) // 2))
                    if p > 0
                    else ComplexDecimal(Decimal(1), Decimal(0), digits=work_prec)
                    / (sqrt_c * _eval_node(node.base ** ((-p - 1) // 2)))
                )
            if exp.is_Rational:
                # General rational power z^(p/q) = exp_P((p/q) * Log_P(z))
                from guiders_math.eval_rc import projective_exp_C, projective_log_C

                log_b = projective_log_C(base_val, work_prec)
                with localcontext() as ctx:
                    ctx.prec = work_prec
                    ratio = Decimal(int(exp.p)) / Decimal(int(exp.q))
                scaled = ComplexDecimal(log_b.re * ratio, log_b.im * ratio, digits=work_prec)
                return projective_exp_C(scaled, work_prec)
        if isinstance(node, sp.cos):
            arg_val = _eval_node(node.args[0])
            c_v, _ = projective_cos_sin_C(arg_val, work_prec)
            return c_v
        if isinstance(node, sp.sin):
            arg_val = _eval_node(node.args[0])
            _, s_v = projective_cos_sin_C(arg_val, work_prec)
            return s_v
        if isinstance(node, sp.exp):
            from guiders_math.eval_rc import projective_exp_C

            return projective_exp_C(_eval_node(node.args[0]), work_prec)
        if isinstance(node, sp.log):
            from guiders_math.eval_rc import projective_log_C

            return projective_log_C(_eval_node(node.args[0]), work_prec)
        if isinstance(node, sp.atan):
            from guiders_math.eval_rc import projective_atan_R

            arg_val = _eval_node(node.args[0])
            if arg_val.is_real():
                return ComplexDecimal(projective_atan_R(arg_val.re, work_prec), Decimal(0), digits=work_prec)
        raise ValueError(f"Unsupported symbolic node in pure-rational evaluator: {node!r}")

    res = _eval_node(expr)
    with localcontext() as ctx:
        ctx.prec = digits
        return ComplexDecimal(+res.re, +res.im, digits=digits)


def hermite_reduction(num: sp.Expr, den: sp.Expr, var: sp.Symbol) -> HermiteDecomposition:
    """Perform exact Ostrowski-Hermite reduction on P(t)/Q(t).

    Decomposes:
        P(t) / Q(t) = poly_part(t) + d/dt [ hermite_rat(t) ] + rem_num(t) / sqfree_den(t),
    where `sqfree_den(t)` is square-free (all roots simple) and deg(rem_num) < deg(sqfree_den).
    """
    num_p = sp.Poly(sp.expand(num), var)
    den_p = sp.Poly(sp.expand(den), var)
    if den_p.is_zero:
        raise ZeroDivisionError("Denominator polynomial Q(t) is zero.")

    # Normalize denominator to be monic
    lc = den_p.LC()
    num_expr = sp.simplify(num_p.as_expr() / lc)
    den_expr = sp.simplify(den_p.as_expr() / lc)

    num_p = sp.Poly(num_expr, var)
    den_p = sp.Poly(den_expr, var)

    # 1. Polynomial quotient and remainder: P = q_poly * Q + r_rem
    q_poly, r_poly = sp.div(num_p, den_p, var)
    poly_part = q_poly.as_expr()
    poly_antideriv = sp.integrate(poly_part, var)

    if r_poly.is_zero:
        return HermiteDecomposition(
            var=var,
            poly_part=poly_part,
            poly_antideriv=poly_antideriv,
            hermite_rat=sp.Integer(0),
            rem_num=sp.Integer(0),
            sqfree_den=sp.Integer(1),
        )

    # 2.Mack / Yun square-free reduction:
    # Let V = Q, A = r_poly. While gcd(V, V') has degree > 0:
    # Let V_minus = gcd(V, V'), V_star = V / V_minus.
    # Actually standard Hermite-Mack iterative step:
    # Write Q = U * V^m via sqf_list.
    _, sqf_factors = sp.sqf_list(den_p.as_expr(), var)
    # Build dictionary {multiplicity: poly}
    mult_map: dict[int, sp.Expr] = {}
    max_m = 1
    for fac, mult in sqf_factors:
        mult_map[mult] = sp.expand(mult_map.get(mult, sp.Integer(1)) * fac)
        if mult > max_m:
            max_m = mult

    a_cur = r_poly.as_expr()
    hermite_acc = sp.Integer(0)

    # Reduce highest multiplicity m down to 1 step-by-step
    while max_m > 1:
        v_m = mult_map.pop(max_m)
        # Product of all other factors U = Q_cur / V_m^max_m
        u_other = sp.Integer(1)
        for k_mult, k_fac in mult_map.items():
            u_other = sp.expand(u_other * (k_fac ** k_mult))

        v_m_diff = sp.diff(v_m, var)
        w_poly = sp.expand(u_other * v_m_diff)
        # Since v_m is square-free and coprime to u_other, gcd(w_poly, v_m) = 1.
        # Solve Bezout identity: B(t)*w_poly + C(t)*v_m = a_cur / (1 - max_m)
        scaled_a = sp.expand(a_cur / sp.Rational(1 - max_m))
        s_bez, t_bez, g_bez = sp.gcdex(sp.Poly(w_poly, var), sp.Poly(v_m, var))
        # s_bez * w_poly + t_bez * v_m = g_bez (a constant != 0)
        g_expr = g_bez.as_expr()
        b_full = sp.Poly(sp.expand(s_bez.as_expr() * scaled_a / g_expr), var)
        _, b_rem = sp.div(b_full, sp.Poly(v_m, var), var)
        b_expr = b_rem.as_expr()

        # Exact algebraic quotient: (scaled_a - b_expr * w_poly) / v_m
        c_quo, _ = sp.div(
            sp.Poly(sp.expand(scaled_a - b_expr * w_poly), var),
            sp.Poly(v_m, var),
            var,
        )
        c_expr = c_quo.as_expr()

        # Update Hermite rational part: add B(t) / V_m^(m - 1)
        hermite_acc = sp.together(hermite_acc + b_expr / (v_m ** (max_m - 1)))

        # Update remaining numerator for denominator U * V_m^(m - 1):
        a_cur = sp.expand((1 - max_m) * c_expr - u_other * sp.diff(b_expr, var))

        # Merge v_m into multiplicity (max_m - 1)
        mult_map[max_m - 1] = sp.expand(mult_map.get(max_m - 1, sp.Integer(1)) * v_m)
        max_m -= 1
        while max_m > 1 and max_m not in mult_map:
            max_m -= 1

    sqfree_den = sp.expand(mult_map.get(1, sp.Integer(1)))
    rem_num = sp.expand(a_cur)
    return HermiteDecomposition(
        var=var,
        poly_part=poly_part,
        poly_antideriv=poly_antideriv,
        hermite_rat=sp.cancel(hermite_acc),
        rem_num=rem_num,
        sqfree_den=sqfree_den,
    )


def _aberth_roots_pure_rational(coeffs: list[ComplexDecimal], digits: int = 50) -> list[ComplexDecimal]:
    """Pure-rational Aberth-Ehrlich simultaneous root solver in C for high-degree polynomials.

    Finds all n complex roots of P(z) = sum_{k=0}^n coeffs[n-k] * z^k = 0 (coeffs[0] = leading coeff)
    to `digits` precision without mpmath.
    """
    work_prec = digits + 25
    n = len(coeffs) - 1
    if n <= 0:
        return []
    with localcontext() as ctx:
        ctx.prec = work_prec
        lc = (coeffs[0].re, coeffs[0].im)
        norm_c = [_c_div((c.re, c.im), lc) for c in coeffs]

        # Cauchy upper bound on root radius: 1 + max |a_k|
        max_a = Decimal(1)
        for c in norm_c[1:]:
            mag = rational_sqrt_decimal(c[0] * c[0] + c[1] * c[1], 25)
            if mag + Decimal(1) > max_a:
                max_a = mag + Decimal(1)

        # Initialize roots on an asymmetrically shifted circle in Q(i) using rational stereographic points
        roots: list[tuple[Decimal, Decimal]] = []
        for k in range(n):
            # Rational parameter t_k = (2k + 1) / (n + 1) + 1/(3n+7)
            t_f = Fraction(2 * k + 1, n + 1) + Fraction(1, 3 * n + 7)
            t_d = Decimal(t_f.numerator) / Decimal(t_f.denominator)
            den = Decimal(1) + t_d * t_d
            c_t = (Decimal(1) - t_d * t_d) / den
            s_t = (Decimal(2) * t_d) / den
            r_scale = max_a * (Decimal(1) - Decimal(k) / Decimal(2 * n + 3))
            if k % 2 == 1:
                c_t = -c_t
            roots.append((r_scale * c_t, r_scale * s_t))

        tol = Decimal(10) ** (-(digits + 16))
        max_aberth_iter = max(120, digits * 2)
        for _ in range(max_aberth_iter):
            max_delta = Decimal(0)
            new_roots = list(roots)
            for i in range(n):
                zi = roots[i]
                # Horner evaluation of P(zi) and P'(zi)
                p_val = norm_c[0]
                dp_val = (Decimal(0), Decimal(0))
                for c_k in norm_c[1:]:
                    dp_val = _c_add(_c_mul(dp_val, zi), p_val)
                    p_val = _c_add(_c_mul(p_val, zi), c_k)

                p_mag2 = p_val[0] * p_val[0] + p_val[1] * p_val[1]
                if p_mag2 <= tol * tol:
                    continue

                newton_step = _c_div(p_val, dp_val)
                # Aberth sum: sum_{j != i} 1 / (z_i - z_j)
                ab_sum = (Decimal(0), Decimal(0))
                for j in range(n):
                    if i == j:
                        continue
                    diff = _c_sub(zi, roots[j])
                    ab_sum = _c_add(ab_sum, _c_div((Decimal(1), Decimal(0)), diff))

                one_minus = _c_sub((Decimal(1), Decimal(0)), _c_mul(newton_step, ab_sum))
                offset = _c_div(newton_step, one_minus)
                new_roots[i] = _c_sub(zi, offset)
                off_mag = abs(offset[0]) + abs(offset[1])
                if off_mag > max_delta:
                    max_delta = off_mag
            roots = new_roots
            if max_delta <= tol:
                break

    with localcontext() as ctx:
        ctx.prec = digits
        return [ComplexDecimal(+r[0], +r[1], digits=digits) for r in roots]


def _exact_roots_of_factor(fac: sp.Expr, var: sp.Symbol) -> list[sp.Expr] | None:
    """Return exact symbolic roots in Q(i, sqrt(d-1)) or radicals for degree <= 4 and cyclotomic/biquadratic polynomials."""
    poly = sp.Poly(fac, var)
    deg = poly.degree()
    if deg == 0:
        return []
    if deg == 1:
        a1, a0 = poly.all_coeffs()
        return [sp.simplify(-a0 / a1)]
    if deg == 2:
        a2, a1, a0 = poly.all_coeffs()
        disc = sp.simplify(a1 * a1 - 4 * a2 * a0)
        s_disc = sp.sqrt(disc)
        return [
            sp.simplify((-a1 + s_disc) / (2 * a2)),
            sp.simplify((-a1 - s_disc) / (2 * a2)),
        ]
    # Check binomial a*t^n + b = 0
    terms = poly.terms()
    if len(terms) == 2 and terms[0][0][0] == deg and terms[1][0][0] == 0:
        an = terms[0][1]
        a0 = terms[1][1]
        # t^n = -a0 / an
        rhs = sp.simplify(-a0 / an)
        # If deg == 4 and rhs == -1 (i.e. t^4 + 1 = 0), return exact Q(i, sqrt(2)) roots!
        if deg == 4 and rhs == -1:
            s2 = sp.sqrt(2) / 2
            return [
                s2 * (1 + sp.I),
                s2 * (-1 + sp.I),
                s2 * (-1 - sp.I),
                s2 * (1 - sp.I),
            ]
        if deg == 4 and rhs.is_positive:
            r4 = sp.simplify(rhs ** sp.Rational(1, 4))
            return [r4, -r4, sp.I * r4, -sp.I * r4]
        if deg == 4 and rhs.is_negative:
            r4 = sp.simplify((-rhs) ** sp.Rational(1, 4))
            s2 = sp.sqrt(2) / 2
            return [
                r4 * s2 * (1 + sp.I),
                r4 * s2 * (-1 + sp.I),
                r4 * s2 * (-1 - sp.I),
                r4 * s2 * (1 - sp.I),
            ]
    # Check biquadratic a4*t^4 + a2*t^2 + a0 = 0
    if deg == 4:
        c4, c3, c2, c1, c0 = poly.all_coeffs()
        if c3 == 0 and c1 == 0:
            u_var = sp.Dummy("u")
            u_roots = _exact_roots_of_factor(c4 * u_var**2 + c2 * u_var + c0, u_var)
            if u_roots is not None:
                out: list[sp.Expr] = []
                for ur in u_roots:
                    sr = sp.sqrt(ur)
                    out.append(sp.simplify(sr))
                    out.append(sp.simplify(-sr))
                return out

    # General irreducible degree >= 3 factors (Cardano/Ferrari/quintic):
    # delegate directly to the pure-rational Aberth-Ehrlich solver in Q(i)
    return None


def lagrange_residues(
    rem_num: sp.Expr,
    sqfree_den: sp.Expr,
    var: sp.Symbol,
    d: int = 3,
    digits: int = 50,
) -> list[PoleResidue]:
    """Compute all poles alpha_k and Lagrange residues c_k = A(alpha_k) / Q'_sqfree(alpha_k)
    of a square-free proper rational function A(t) / Q_sqfree(t).
    """
    rem_num = sp.expand(rem_num)
    sqfree_den = sp.expand(sqfree_den)
    if rem_num == 0 or sqfree_den == 1:
        return []

    den_diff = sp.diff(sqfree_den, var)
    # Factor over Q (or Q(i) if I is present in coefficients) — instant O(1) factorization
    if sqfree_den.has(sp.I):
        _, fac_list = sp.factor_list(sqfree_den, var, extension=[sp.I])
    else:
        _, fac_list = sp.factor_list(sqfree_den, var)

    results: list[PoleResidue] = []

    for fac, _mult in fac_list:
        exact_roots = _exact_roots_of_factor(fac, var)
        if exact_roots is not None:
            for r_sym in exact_roots:
                res_sym = sp.simplify(rem_num.subs(var, r_sym) / den_diff.subs(var, r_sym))
                pole_q: QQuadExt | None = None
                res_q: QQuadExt | None = None
                try:
                    pole_q = QQuadExt.from_sympy(r_sym, d=d)
                    res_q = QQuadExt.from_sympy(res_sym, d=d)
                except Exception:
                    pass
                results.append(
                    PoleResidue(
                        pole_sym=r_sym,
                        residue_sym=res_sym,
                        pole_q=pole_q,
                        residue_q=res_q,
                    )
                )
        else:
            # High-degree irreducible factor: solve roots in C via pure-rational Aberth-Newton
            work_prec = digits + 25
            f_poly = sp.Poly(fac, var)
            c_list = [
                eval_sympy_alg_to_complex_decimal(c, digits=work_prec)
                for c in f_poly.all_coeffs()
            ]
            num_poly = sp.Poly(rem_num, var)
            dd_poly = sp.Poly(den_diff, var)
            num_coeffs = [
                eval_sympy_alg_to_complex_decimal(c, digits=work_prec)
                for c in num_poly.all_coeffs()
            ]
            dd_coeffs = [
                eval_sympy_alg_to_complex_decimal(c, digits=work_prec)
                for c in dd_poly.all_coeffs()
            ]

            def eval_horner(cs: list[ComplexDecimal], z: ComplexDecimal) -> ComplexDecimal:
                acc = ComplexDecimal(Decimal(0), Decimal(0), digits=work_prec)
                for c_item in cs:
                    acc = acc * z + c_item
                return acc

            roots_num = _aberth_roots_pure_rational(c_list, digits=work_prec)
            for r_c in roots_num:
                a_val = eval_horner(num_coeffs, r_c)
                qp_val = eval_horner(dd_coeffs, r_c)
                res_c = a_val / qp_val
                results.append(
                    PoleResidue(
                        pole_sym=None,
                        residue_sym=None,
                        pole_q=None,
                        residue_q=None,
                        pole_num=r_c,
                        residue_num=res_c,
                    )
                )

    return results
