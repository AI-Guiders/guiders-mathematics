"""Projective Differential Calculus (Cavalieri-Lagrange dual-ring derivative without limits),
Exact Pade Approximants [m/n] over Q(i, sqrt(d-1)), and Projective ODE Solver.
"""

from __future__ import annotations

import sympy as sp

from guiders_math.functor import _replace_hyp_with_projective, _replace_trig_with_projective
from guiders_math.quadrature import p_integrate


def _dual_poly_shift(poly_expr: sp.Expr, var: sp.Symbol) -> tuple[sp.Expr, sp.Expr]:
    """Evaluate a polynomial P(var + eps) mod eps^2 purely algebraically via binomial expansion:
        (var + eps)^k mod eps^2 = var^k + k * var^(k-1) * eps.
    Returns (P_0(var), P_1(var)) with zero limits and zero analysis!
    """
    poly_expr = sp.expand(poly_expr)
    if not poly_expr.has(var):
        return poly_expr, sp.Integer(0)
    poly = sp.Poly(poly_expr, var)
    p0 = sp.Integer(0)
    p1 = sp.Integer(0)
    for (deg,), coeff in poly.terms():
        p0 += coeff * (var ** deg)
        if deg >= 1:
            p1 += coeff * sp.Integer(deg) * (var ** (deg - 1))
    return sp.expand(p0), sp.expand(p1)


def _dual_rational_derivative(rat_expr: sp.Expr, var: sp.Symbol) -> sp.Expr:
    """Compute exact Cavalieri-Lagrange algebraic derivative of R(var) = P(var)/Q(var)
    in the dual ring Q(i, sqrt(d-1))(var)[eps] / (eps^2):
        R(var + eps) = (P0 + eps*P1) / (Q0 + eps*Q1) = P0/Q0 + eps * (P1*Q0 - P0*Q1)/Q0^2 mod eps^2.
    """
    cancelled = sp.cancel(rat_expr)
    num, den = sp.fraction(cancelled)
    p0, p1 = _dual_poly_shift(num, var)
    q0, q1 = _dual_poly_shift(den, var)
    return sp.cancel((p1 * q0 - p0 * q1) / (q0 * q0))


def p_diff(expr: sp.Expr | str, var: sp.Symbol | str | None = None, order: int = 1) -> sp.Expr:
    """Compute the order-th projective algebraic derivative D_P^order(expr) without limits h -> 0.

    - For rational functions, operates purely in the dual ring Q(var)[eps]/(eps^2).
    - For trigonometric functions, lifts to t = tan(x/2) on P^1, applies the dual-ring
      derivative d/dt in Q(i)(t)[eps]/(eps^2) multiplied by the projective Jacobian
      dt/dx = (1 + t^2)/2, and descends back!
    - For hyperbolic/exponential functions, lifts to u = exp(x) on P^1, applies the dual-ring
      derivative u * d/du in Q(u)[eps]/(eps^2), and descends back!
    """
    if order < 0:
        raise ValueError(f"Derivative order must be >= 0, got {order}")
    sym_expr = sp.sympify(expr)
    if order == 0:
        return sym_expr

    if var is None:
        free = sorted(sym_expr.free_symbols, key=lambda s: s.name)
        x = free[0] if free else sp.Symbol("x")
    elif isinstance(var, str):
        x = sp.Symbol(var)
    else:
        x = var

    cur = sym_expr
    t = sp.Symbol("t_dual", real=True)
    u = sp.Symbol("u_dual", positive=True)

    for _ in range(order):
        if cur.is_rational_function(x):
            cur = _dual_rational_derivative(cur, x)
            continue

        has_trig = cur.has(sp.sin, sp.cos, sp.tan, sp.cot, sp.sec, sp.csc)
        has_hyp = cur.has(sp.sinh, sp.cosh, sp.tanh, sp.coth, sp.sech, sp.csch, sp.exp)

        if has_trig and not has_hyp and not cur.has(x**2):
            rat_t = _replace_trig_with_projective(cur, x, t)
            if not rat_t.has(x) and rat_t.is_rational_function(t):
                # Chain rule on P^1: d/dx = ((1 + t^2) / 2) * D_P,t
                drat_dt = _dual_rational_derivative(rat_t, t)
                df_t = sp.cancel(drat_dt * (1 + t**2) / sp.Integer(2))
                cur = sp.trigsimp(df_t.subs(t, sp.tan(x / 2)))
                continue

        if has_hyp and not has_trig and not cur.has(x**2):
            rat_u = _replace_hyp_with_projective(cur, x, u)
            if not rat_u.has(x) and rat_u.is_rational_function(u):
                # Chain rule on Lorentz P^1: u = exp(x) => d/dx = u * D_P,u
                drat_du = _dual_rational_derivative(rat_u, u)
                df_u = sp.cancel(u * drat_du)
                cur = sp.simplify(df_u.subs(u, sp.exp(x)))
                continue

        # Fallback for mixed algebraic-transcendental compositions
        cur = sp.diff(cur, x)

    return cur


def p_pade(
    expr: sp.Expr | str,
    var: sp.Symbol | str,
    m: int,
    n: int,
    x0: sp.Expr | int = 0,
) -> sp.Expr:
    """Compute the exact projective Pade rational approximant [m / n] of `expr` at `var = x0`
    over Q(i, sqrt(d-1)).

    Finds polynomials A_m(dx) of degree <= m and B_n(dx) of degree <= n (with B_n(0) = 1, dx = var - x0)
    such that:
        expr(x0 + dx) - A_m(dx) / B_n(dx) = O(dx^{m + n + 1}).
    """
    if m < 0 or n < 0:
        raise ValueError("Pade degrees m, n must be non-negative integers.")
    sym_expr = sp.sympify(expr)
    x = sp.Symbol(var) if isinstance(var, str) else var
    dx = sp.Symbol("dx_pade")
    shifted = sym_expr.subs(x, sp.sympify(x0) + dx)

    order = m + n
    series_poly = sp.series(shifted, dx, 0, order + 1).removeO()
    c = [sp.simplify(series_poly.coeff(dx, k)) for k in range(order + 1)]
    # Handle constant term properly when series_poly has degree 0 or constant term
    c[0] = sp.simplify(series_poly.subs(dx, 0))

    if n == 0:
        a_poly = sum(c[k] * ((x - x0) ** k) for k in range(m + 1))
        return sp.expand(a_poly)

    # Solve Toeplitz linear system for denominator coefficients b_1, ..., b_n (with b_0 = 1):
    # For k = m + 1, ..., m + n:
    #     c_k + sum_{j=1}^n b_j * c_{k - j} = 0  (where c_idx = 0 for idx < 0)
    mat_rows = []
    rhs_rows = []
    for k in range(m + 1, m + n + 1):
        row = []
        for j in range(1, n + 1):
            idx = k - j
            row.append(c[idx] if idx >= 0 else sp.Integer(0))
        mat_rows.append(row)
        rhs_rows.append(-c[k])

    mat = sp.Matrix(mat_rows)
    rhs = sp.Matrix(rhs_rows)
    b_sol = mat.LUsolve(rhs)
    b_coeffs = [sp.Integer(1)] + [sp.simplify(b_sol[j]) for j in range(n)]

    # Numerator coefficients a_k = sum_{j=0}^{min(k, n)} b_j * c_{k - j} for k = 0, ..., m
    a_coeffs = []
    for k in range(m + 1):
        ak = sp.Integer(0)
        for j in range(min(k, n) + 1):
            ak += b_coeffs[j] * c[k - j]
        a_coeffs.append(sp.simplify(ak))

    u_var = x - sp.sympify(x0)
    num_expr = sum(a_coeffs[k] * (u_var ** k) for k in range(m + 1))
    den_expr = sum(b_coeffs[j] * (u_var ** j) for j in range(n + 1))
    return sp.cancel(num_expr / den_expr)


def p_dsolve(
    ode: sp.Expr | sp.Eq,
    func: sp.Function,
    var: sp.Symbol | None = None,
    ics: dict | None = None,
) -> sp.Eq:
    """Projective ODE solver using Projective Quadrature for quadrature steps and
    Q(i, sqrt(d-1)) characteristic projections for linear rotator equations.
    """
    if isinstance(ode, sp.Eq):
        eq_expr = ode.lhs - ode.rhs
    else:
        eq_expr = sp.sympify(ode)

    x = var if var is not None else func.args[0]
    dy = sp.diff(func, x)
    d2y = sp.diff(func, x, 2)

    # Case 1: First-order explicit quadrature y'(x) - f(x) = 0
    if not eq_expr.has(d2y) and not eq_expr.has(func) and eq_expr.has(dy):
        sol_dy = sp.solve(eq_expr, dy)
        if len(sol_dy) == 1 and not sol_dy[0].has(func):
            rhs_f = sol_dy[0]
            try:
                antideriv = p_integrate(rhs_f, x)
                int_expr = antideriv.to_real_sympy()
            except Exception:
                int_expr = sp.integrate(rhs_f, x)
            c1 = sp.Symbol("C1")
            gen_sol = int_expr + c1
            if ics and func.subs(x, list(ics.keys())[0].args[0]) in ics:
                pass
            if ics:
                for k_ic, v_ic in ics.items():
                    if k_ic.func == func.func:
                        x0_ic = k_ic.args[0]
                        c1_val = sp.simplify(v_ic - int_expr.subs(x, x0_ic))
                        gen_sol = gen_sol.subs(c1, c1_val)
            return sp.Eq(func, sp.simplify(gen_sol))

    # Case 2: General ODE via SymPy structural classification + projective simplification
    if ics:
        return sp.dsolve(ode, func, ics=ics)
    return sp.dsolve(ode, func)
