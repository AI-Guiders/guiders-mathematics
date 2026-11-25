"""Input and output projective functors P_in and P_out over Q(i, sqrt(d-1)).

Transforms transcendental expressions (trigonometric, hyperbolic, exponential,
quadratic radicals, and Fresnel-Airy polynomial phase kernels) into exact
rational forms on P^1(Q(i, sqrt(d-1))).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import sympy as sp

from guiders_math.eval_rc import (
    ComplexDecimal,
    eval_to_complex_decimal,
    projective_atan_R,
    projective_cos_sin_C,
    projective_exp_C,
    projective_log_C,
    rational_sqrt_decimal,
)


@dataclass(frozen=True)
class ProjectiveLift:
    """Result of lifting an integrand f(x) dx onto the projective line P^1.

    Attributes:
        kind: "rational" | "trig" | "exp" | "euler_circle" | "fresnel_phase"
        t_var: Projective coordinate symbol t on P^1
        rat_num: Numerator polynomial P(t) of the lifted rational integrand
        rat_den: Denominator polynomial Q(t) of the lifted rational integrand
        forward_sym: Symbolic forward map t(x)
        forward_eval: Pure-rational evaluator mapping x in R or C (or oo) to t in C (or oo)
        phase_meta: Metadata dictionary for Fresnel / Gaussian sector rotation integrals
    """

    kind: str
    orig_var: sp.Symbol
    t_var: sp.Symbol
    rat_num: sp.Expr
    rat_den: sp.Expr
    forward_sym: sp.Expr | None = None
    forward_eval: Callable[[ComplexDecimal | str, int], ComplexDecimal | str] | None = None
    phase_meta: dict | None = None


def _replace_trig_with_projective(expr: sp.Expr, x: sp.Symbol, t: sp.Symbol) -> sp.Expr:
    """Replace all trig functions of integer or half-integer multiples of x using the Cayley rational power:
        t = tan(x/2),
        cos(n*x) + i*sin(n*x) = ((1 + i*t) / (1 - i*t))^n.
    """
    u_unit = (1 + sp.I * t) / (1 - sp.I * t)
    u_inv = (1 - sp.I * t) / (1 + sp.I * t)

    def _walk(node: sp.Expr) -> sp.Expr:
        if not node.has(x):
            return node
        if isinstance(node, (sp.sin, sp.cos, sp.tan, sp.cot, sp.sec, sp.csc)):
            arg = sp.expand(node.args[0])
            coeff = sp.simplify(arg / x)
            if coeff == sp.Rational(1, 2):
                if isinstance(node, sp.tan):
                    return t
                if isinstance(node, sp.cot):
                    return 1 / t
            if coeff.is_Integer:
                n = int(coeff)
                un = u_unit ** n
                u_minus_n = u_inv ** n
                c_nx = (un + u_minus_n) / 2
                s_nx = (un - u_minus_n) / (2 * sp.I)
                if isinstance(node, sp.cos):
                    return sp.cancel(c_nx)
                if isinstance(node, sp.sin):
                    return sp.cancel(s_nx)
                if isinstance(node, sp.tan):
                    return sp.cancel(s_nx / c_nx)
                if isinstance(node, sp.cot):
                    return sp.cancel(c_nx / s_nx)
                if isinstance(node, sp.sec):
                    return sp.cancel(1 / c_nx)
                if isinstance(node, sp.csc):
                    return sp.cancel(1 / s_nx)
        if node.args:
            return node.func(*[_walk(a) for a in node.args])
        return node

    return sp.cancel(_walk(expr))


def _replace_hyp_with_projective(expr: sp.Expr, x: sp.Symbol, u: sp.Symbol) -> sp.Expr:
    """Replace hyperbolic and exponential functions of integer multiples of x with u = exp(x)."""
    def _walk(node: sp.Expr) -> sp.Expr:
        if not node.has(x):
            return node
        if isinstance(node, (sp.sinh, sp.cosh, sp.tanh, sp.coth, sp.sech, sp.csch, sp.exp)):
            arg = sp.expand(node.args[0])
            coeff = sp.simplify(arg / x)
            if coeff.is_Integer:
                n = int(coeff)
                un = u ** n
                u_minus_n = u ** (-n)
                if isinstance(node, sp.exp):
                    return un
                ch_nx = (un + u_minus_n) / 2
                sh_nx = (un - u_minus_n) / 2
                if isinstance(node, sp.cosh):
                    return sp.cancel(ch_nx)
                if isinstance(node, sp.sinh):
                    return sp.cancel(sh_nx)
                if isinstance(node, sp.tanh):
                    return sp.cancel(sh_nx / ch_nx)
                if isinstance(node, sp.coth):
                    return sp.cancel(ch_nx / sh_nx)
                if isinstance(node, sp.sech):
                    return sp.cancel(1 / ch_nx)
                if isinstance(node, sp.csch):
                    return sp.cancel(1 / sh_nx)
        if node.args:
            return node.func(*[_walk(a) for a in node.args])
        return node

    return sp.cancel(_walk(expr))


def p_simplify(expr: sp.Expr | str, var: sp.Symbol | None = None) -> sp.Expr:
    """Decision-procedure projective simplification over Q(i, sqrt(d-1))(t).

    Lifts trigonometric and hyperbolic functions to rational functions on P^1,
    performs exact polynomial GCD cancellation in Q(i, sqrt(d-1))[t],
    and descends back if transcendental variables remain.
    """
    sym_expr = sp.sympify(expr)
    free_syms = sorted(sym_expr.free_symbols, key=lambda s: s.name)
    if not free_syms:
        return sp.simplify(sym_expr)

    x = var if var is not None else free_syms[0]
    t = sp.Symbol("t_proj", real=True)

    has_trig = sym_expr.has(sp.sin, sp.cos, sp.tan, sp.cot, sp.sec, sp.csc)
    has_hyp = sym_expr.has(sp.sinh, sp.cosh, sp.tanh, sp.coth, sp.sech, sp.csch, sp.exp)

    work = sym_expr
    if has_trig and not has_hyp:
        proj_rat = _replace_trig_with_projective(work, x, t)
        if not proj_rat.has(x):
            # Pure rational function in t on P^1!
            cancelled = sp.cancel(proj_rat)
            if not cancelled.has(t):
                return sp.simplify(cancelled)
            return sp.simplify(cancelled.subs(t, sp.tan(x / 2)))
    if has_hyp and not has_trig:
        u = sp.Symbol("u_proj", positive=True)
        proj_rat = _replace_hyp_with_projective(work, x, u)
        if not proj_rat.has(x):
            cancelled = sp.cancel(proj_rat)
            if not cancelled.has(u):
                return sp.simplify(cancelled)
            return sp.simplify(cancelled.subs(u, sp.exp(x)))

    return sp.cancel(sp.simplify(sym_expr))


def p_prove_identity(lhs: sp.Expr | str, rhs: sp.Expr | str, var: sp.Symbol | None = None) -> bool:
    """Prove whether lhs == rhs identically by lifting lhs - rhs to P^1 and checking P(t) == 0."""
    diff_expr = sp.sympify(lhs) - sp.sympify(rhs)
    simp = p_simplify(diff_expr, var=var)
    return sp.simplify(simp) == 0


def _detect_fresnel_or_gaussian_phase(expr: sp.Expr, x: sp.Symbol) -> dict | None:
    """Detect oscillating Fresnel-Airy or Gaussian phase integrands of the form:
        c * x^m * sin(a * x^n),
        c * x^m * cos(a * x^n),
        c * x^m * exp(alpha * x^n)
    with integer n >= 2 and integer m >= 0.
    """
    expr = sp.expand(expr)
    # Handle single term c * x^m * func(a * x^n)
    coeff, dep = expr.as_independent(x, as_Add=False)
    # Separate power x^m from trig/exp factor
    m_pow = 0
    core_func = dep
    if dep.is_Mul:
        other_factors = []
        for f in dep.args:
            if f == x:
                m_pow += 1
            elif f.is_Pow and f.base == x and f.exp.is_Integer and int(f.exp) >= 1:
                m_pow += int(f.exp)
            else:
                other_factors.append(f)
        if len(other_factors) == 1:
            core_func = other_factors[0]
        else:
            return None

    if isinstance(core_func, (sp.sin, sp.cos, sp.exp)):
        arg = sp.expand(core_func.args[0])
        a_coeff, x_part = arg.as_independent(x, as_Add=False)
        n_deg = 0
        if x_part.is_Pow and x_part.base == x and x_part.exp.is_Integer:
            n_deg = int(x_part.exp)
        if n_deg >= 2 and a_coeff != 0:
            kind_fn = (
                "sin"
                if isinstance(core_func, sp.sin)
                else "cos"
                if isinstance(core_func, sp.cos)
                else "exp"
            )
            return {
                "fn": kind_fn,
                "coeff": coeff,
                "m": m_pow,
                "a": a_coeff,
                "n": n_deg,
            }
    return None


def lift_integrand_to_projective(expr: sp.Expr | str, var: sp.Symbol) -> ProjectiveLift:
    """Lift an integrand f(x) dx to a rational integrand (P(t)/Q(t)) dt on P^1(Q(i, sqrt(d-1)))."""
    sym_expr = sp.sympify(expr)
    t = sp.Symbol("t", real=True)

    # 1. Check Fresnel / Airy / Gaussian polynomial phase kernel x^m * sin(a*x^n), etc.
    phase_info = _detect_fresnel_or_gaussian_phase(sym_expr, var)
    if phase_info is not None:
        # Spinor half-angle lift on P^1: u = a*x^n, stereographic carrier on phase u
        return ProjectiveLift(
            kind="fresnel_phase",
            orig_var=var,
            t_var=t,
            rat_num=sp.Integer(2) * t,
            rat_den=1 + t**4,
            forward_sym=var,
            forward_eval=None,
            phase_meta=phase_info,
        )

    # 2. Check pure rational function in var over Q(i, sqrt(d-1))
    if sym_expr.is_rational_function(var):
        cancelled = sp.cancel(sym_expr)
        num, den = sp.fraction(cancelled)
        num_t = num.subs(var, t)
        den_t = den.subs(var, t)

        def _eval_identity(val: ComplexDecimal | str, digits: int) -> ComplexDecimal | str:
            return val

        return ProjectiveLift(
            kind="rational",
            orig_var=var,
            t_var=t,
            rat_num=num_t,
            rat_den=den_t,
            forward_sym=var,
            forward_eval=_eval_identity,
        )

    # 3. Check trigonometric integrand R(sin(x), cos(x), tan(x), ...)
    has_trig = sym_expr.has(sp.sin, sp.cos, sp.tan, sp.cot, sp.sec, sp.csc)
    if has_trig and not sym_expr.has(var**2):
        # Weierstrass-stereographic lift: t = tan(x/2), dx = 2 / (1 + t^2) dt
        lifted_f = _replace_trig_with_projective(sym_expr, var, t)
        if not lifted_f.has(var) and lifted_f.is_rational_function(t):
            full_rat = sp.cancel(lifted_f * sp.Integer(2) / (1 + t**2))
            num_t, den_t = sp.fraction(full_rat)

            def _eval_tan_half(val: ComplexDecimal | str, digits: int) -> ComplexDecimal | str:
                if isinstance(val, str):
                    return val
                work_prec = digits + 20
                from decimal import Decimal, localcontext

                with localcontext() as ctx:
                    ctx.prec = work_prec
                    half_z = ComplexDecimal(val.re / Decimal(2), val.im / Decimal(2), digits=work_prec)
                    c_h, s_h = projective_cos_sin_C(half_z, work_prec)
                    if abs(c_h.re) + abs(c_h.im) <= Decimal(10) ** (-(digits + 5)):
                        return "oo"
                    return s_h / c_h

            return ProjectiveLift(
                kind="trig",
                orig_var=var,
                t_var=t,
                rat_num=num_t,
                rat_den=den_t,
                forward_sym=sp.tan(var / 2),
                forward_eval=_eval_tan_half,
            )

    # 4. Check hyperbolic / exponential integrand R(sinh(x), cosh(x), exp(x))
    has_hyp = sym_expr.has(sp.sinh, sp.cosh, sp.tanh, sp.coth, sp.sech, sp.csch, sp.exp)
    if has_hyp:
        lifted_f = _replace_hyp_with_projective(sym_expr, var, t)
        if not lifted_f.has(var) and lifted_f.is_rational_function(t):
            # u = exp(x) => dx = du / u
            full_rat = sp.cancel(lifted_f / t)
            num_t, den_t = sp.fraction(full_rat)

            def _eval_exp(val: ComplexDecimal | str, digits: int) -> ComplexDecimal | str:
                if val == "oo":
                    return "oo"
                if val == "-oo":
                    from decimal import Decimal

                    return ComplexDecimal(Decimal(0), Decimal(0), digits=digits)
                assert isinstance(val, ComplexDecimal)
                return projective_exp_C(val, digits=digits)

            return ProjectiveLift(
                kind="exp",
                orig_var=var,
                t_var=t,
                rat_num=num_t,
                rat_den=den_t,
                forward_sym=sp.exp(var),
                forward_eval=_eval_exp,
            )

    # 5. Check quadratic circle radical sqrt(a^2 - x^2) or sqrt(x^2 + a^2)
    # Example: Substitution x = a * 2*t / (1 + t^2), sqrt(a^2 - x^2) = a*(1 - t^2)/(1 + t^2), dx = 2*a*(1 - t^2)/(1 + t^2)^2 dt
    pows = [
        p
        for p in sym_expr.atoms(sp.Pow)
        if p.exp == sp.Rational(1, 2) and p.base.is_polynomial(var) and sp.Poly(p.base, var).degree() == 2
    ]
    if len(pows) == 1:
        rad = pows[0]
        q_poly = sp.Poly(rad.base, var)
        c2, c1, c0 = q_poly.all_coeffs()
        if c1 == 0 and c2 == -1 and c0 > 0:
            a_rad = sp.sqrt(c0)
            x_sub = a_rad * (2 * t) / (1 + t**2)
            rad_sub = a_rad * (1 - t**2) / (1 + t**2)
            dx_sub = a_rad * (2 * (1 - t**2)) / ((1 + t**2) ** 2)
            sub_expr = sp.cancel(sym_expr.subs({rad: rad_sub, var: x_sub}) * dx_sub)
            if sub_expr.is_rational_function(t):
                num_t, den_t = sp.fraction(sub_expr)

                def _eval_euler_circle(val: ComplexDecimal | str, digits: int) -> ComplexDecimal | str:
                    if isinstance(val, str):
                        return val
                    # t = x / (a + sqrt(a^2 - x^2))
                    work_prec = digits + 20
                    from decimal import Decimal, localcontext
                    from guiders_math.poly_engine import eval_sympy_alg_to_complex_decimal

                    a_c = eval_sympy_alg_to_complex_decimal(a_rad, work_prec)
                    c0_c = eval_sympy_alg_to_complex_decimal(c0, work_prec)
                    diff_c = c0_c - (val * val)
                    from guiders_math.eval_rc import complex_sqrt_heron

                    sr, si = complex_sqrt_heron(diff_c.re, diff_c.im, work_prec)
                    s_c = ComplexDecimal(sr, si, digits=work_prec)
                    return val / (a_c + s_c)

                return ProjectiveLift(
                    kind="euler_circle",
                    orig_var=var,
                    t_var=t,
                    rat_num=num_t,
                    rat_den=den_t,
                    forward_sym=var / (a_rad + sp.sqrt(c0 - var**2)),
                    forward_eval=_eval_euler_circle,
                )

    raise ValueError(
        f"Integrand {sym_expr!r} could not be lifted to a rational form on P^1(Q(i, sqrt(d-1)))."
    )
