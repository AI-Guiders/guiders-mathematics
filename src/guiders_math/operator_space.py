"""Projective Functional Calculus for Infinite-Dimensional Hilbert Space Operators.

Given the power moments mu_k = <psi_0 | L^k | psi_0> (k = 0, ..., 2N - 1) of a linear
operator L : H -> H on an infinite-dimensional Hilbert space H:
  1. Constructs the exact [N-1 / N] Stieltjes-Pade projective resolvent on the Riemann sphere P^1(C):
         G_N(z) = <psi_0 | (z*I - L)^(-1) | psi_0>_{[N-1/N]} = P_{N-1}(z) / Q_N(z)
     by solving the exact Hankel moment system over Q(i, sqrt(d-1)).
  2. Applies Ostrowski-Hermite and Lagrange-Rothstein-Trager residue decomposition to G_N(z):
         G_N(z) = sum_{k=1}^N w_k / (z - lambda_k),
     where the poles `lambda_k` are the spectral eigenvalues and the Lagrange residues `w_k`
     are the exact spectral measure weights (satisfying sum_k w_k = mu_0 = ||psi_0||^2, i.e.
     -Res_{z=oo} G_N(z) = mu_0 at projective infinity [1 : 0] in P^1(C)).
  3. Evaluates any holomorphic operator function <psi_0 | f(L) | psi_0> via the Riesz-Dunford
     contour integral on P^1(C):
         <psi_0 | f(L) | psi_0> = (1 / (2*pi*i)) * oint_gamma f(z) G_N(z) dz = sum_{k=1}^N w_k * f(lambda_k)
     to arbitrary precision `digits` using pure-rational projective descent.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
from typing import Any, Sequence

import sympy as sp

from guiders_math.eval_rc import (
    ComplexDecimal,
    eval_to_complex_decimal,
    projective_cos_sin_C,
    projective_exp_C,
    projective_log_C,
)
from guiders_math.poly_engine import (
    PoleResidue,
    eval_sympy_alg_to_complex_decimal,
    hermite_reduction,
    lagrange_residues,
)


@dataclass(frozen=True)
class ProjectiveResolventResult:
    """Exact Stieltjes-Pade projective resolvent G_N(z) = P_{N-1}(z) / Q_N(z) and spectral decomposition."""

    z_var: sp.Symbol
    order_n: int
    moments: list[sp.Expr]
    numerator: sp.Expr
    denominator: sp.Expr
    resolvent_expr: sp.Expr
    spectral_poles: list[PoleResidue]

    def total_spectral_weight_sympy(self) -> sp.Expr:
        """Exact sum of spectral weights sum_k w_k = -Res_{z=oo} G_N(z) = mu_0."""
        return self.moments[0] if self.moments else sp.Integer(0)

    def eval_resolvent_at(self, z_val: Any, digits: int = 50) -> ComplexDecimal:
        """Evaluate G_N(z_val) on R or C to `digits` precision."""
        from guiders_math.quadrature import _eval_poly_rational_at

        z_c = eval_to_complex_decimal(z_val, digits=digits + 15)
        return _eval_poly_rational_at(self.resolvent_expr, self.z_var, z_c, digits=digits)

    def eval_operator_function(self, fn_kind: str, param: Any = 1, digits: int = 50) -> ComplexDecimal:
        """Evaluate Riesz-Dunford functional calculus <psi_0 | f(L) | psi_0> = sum_{k=1}^N w_k * f(lambda_k)
        using pure-rational projective descent on R and C.

        Supported `fn_kind`:
          - "exp":      f(L) = exp(param * L)       (time-evolution operator)
          - "cos":      f(L) = cos(param * L)       (even wave-evolution operator)
          - "sin":      f(L) = sin(param * L)       (odd wave-evolution operator)
          - "log":      f(L) = Log_P(param * L)     (operator logarithm)
          - "resolvent":f(L) = (param * I - L)^(-1) (operator resolvent)
        """
        work_prec = digits + 20
        p_c = eval_to_complex_decimal(param, digits=work_prec)
        acc = ComplexDecimal(Decimal(0), Decimal(0), digits=work_prec)

        for pr in self.spectral_poles:
            lam_c = pr.eval_pole(digits=work_prec)
            w_c = pr.eval_residue(digits=work_prec)
            arg_c = p_c * lam_c
            if fn_kind == "exp":
                f_val = projective_exp_C(arg_c, digits=work_prec)
            elif fn_kind == "cos":
                f_val, _ = projective_cos_sin_C(arg_c, digits=work_prec)
            elif fn_kind == "sin":
                _, f_val = projective_cos_sin_C(arg_c, digits=work_prec)
            elif fn_kind == "log":
                f_val = projective_log_C(arg_c, digits=work_prec)
            elif fn_kind == "resolvent":
                one = ComplexDecimal(Decimal(1), Decimal(0), digits=work_prec)
                f_val = one / (p_c - lam_c)
            else:
                raise ValueError(f"Unsupported operator function kind: {fn_kind!r}")
            acc = acc + (w_c * f_val)

        with localcontext() as ctx:
            ctx.prec = digits
            return ComplexDecimal(+acc.re, +acc.im, digits=digits)


def p_resolvent_spectrum(
    moments: Sequence[Any],
    z_var: sp.Symbol | str = "z",
    *,
    d: int = 3,
    digits: int = 50,
) -> ProjectiveResolventResult:
    """Construct the exact [N-1 / N] projective Stieltjes-Pade resolvent G_N(z) = P_{N-1}(z) / Q_N(z)
    and its Lagrange spectral decomposition from 2N operator moments [mu_0, mu_1, ..., mu_{2N-1}].
    """
    if len(moments) < 2 or len(moments) % 2 != 0:
        raise ValueError("p_resolvent_spectrum requires an even number 2*N >= 2 of operator moments.")

    sym_moments = [sp.sympify(m_val) for m_val in moments]
    n_order = len(sym_moments) // 2
    z = sp.Symbol(z_var) if isinstance(z_var, str) else z_var

    # Hankel linear system for denominator Q_N(z) = z^N + c_{N-1} z^{N-1} + ... + c_0:
    # For i = 0, ..., N - 1:
    #     sum_{j=0}^{N-1} mu_{i + j} * c_j = -mu_{i + N}
    hankel_rows = []
    rhs_rows = []
    for i in range(n_order):
        row = [sym_moments[i + j] for j in range(n_order)]
        hankel_rows.append(row)
        rhs_rows.append(-sym_moments[i + n_order])

    h_mat = sp.Matrix(hankel_rows)
    rhs_vec = sp.Matrix(rhs_rows)
    c_sol = h_mat.LUsolve(rhs_vec)
    q_coeffs = [sp.simplify(c_sol[j]) for j in range(n_order)] + [sp.Integer(1)]  # c_0, c_1, ..., c_{N-1}, c_N=1

    den_poly = sum(q_coeffs[j] * (z**j) for j in range(n_order + 1))

    # Numerator P_{N-1}(z) = sum_{m=0}^{N-1} p_m * z^m from matching sum_{k=0}^{N-1} mu_k / z^{k+1}:
    # p_m = sum_{k=0}^{N - 1 - m} q_coeffs[m + 1 + k] * mu_k
    num_poly = sp.Integer(0)
    for m_deg in range(n_order):
        pm = sp.Integer(0)
        for k in range(n_order - m_deg):
            pm += q_coeffs[m_deg + 1 + k] * sym_moments[k]
        num_poly += sp.simplify(pm) * (z**m_deg)

    resolvent = sp.cancel(num_poly / den_poly)
    num_c, den_c = sp.fraction(resolvent)
    herm = hermite_reduction(num_c, den_c, z)
    poles = lagrange_residues(herm.rem_num, herm.sqfree_den, z, d=d, digits=digits)

    return ProjectiveResolventResult(
        z_var=z,
        order_n=n_order,
        moments=sym_moments,
        numerator=sp.expand(num_c),
        denominator=sp.expand(den_c),
        resolvent_expr=resolvent,
        spectral_poles=poles,
    )
