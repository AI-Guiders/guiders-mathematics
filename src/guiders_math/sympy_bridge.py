"""Drop-in SymPy kernel bridge for external CAS frontends (SymvaCAS, SymPad, Jupyter, etc.).

Provides:
  - `projective_sympy_primitive(expression, variable)`: exact projective antiderivative over
    Q(i, sqrt(d-1)) with Ostrowski-Hermite and Lagrange-Rothstein-Trager step breakdown.
  - `projective_sympy_definite(expression, variable, lower, upper)`: exact definite and
    improper integration on P^1(C) (including Fresnel-Airy sector rotations zeta_{2n}).
  - `install_sympy_kernel()`: monkey-patches `sympy.integrate` and `sympy.trigsimp` so any
    existing SymPy wrapper UI transparently uses the projective Q(i, sqrt(d-1)) kernel.
"""

from __future__ import annotations

from typing import Any

import sympy as sp

from guiders_math.functor import p_simplify
from guiders_math.quadrature import DefiniteIntegralResult, ProjectiveAntiderivative, p_antiderivative, p_integrate

_ORIGINAL_SYMPY_INTEGRATE = sp.integrate
_ORIGINAL_SYMPY_TRIGSIMP = sp.trigsimp
_KERNEL_INSTALLED = False


def projective_sympy_primitive(
    expression: sp.Expr,
    variable: sp.Symbol,
    d: int = 3,
) -> tuple[sp.Expr, list[sp.Basic], str, dict[str, Any]] | None:
    """Compute an exact closed-form antiderivative over Q(i, sqrt(d-1)) and return
    (primitive, conditions, note, step_metadata) compatible with SymvaCAS `calc_integrals.py`.
    """
    try:
        antideriv = p_antiderivative(expression, variable, d=d)
    except Exception:
        return None

    # Ensure all poles were resolved symbolically
    if any(pr.pole_sym is None or pr.residue_sym is None for pr in antideriv.residues):
        return None

    real_prim = antideriv.to_real_sympy()
    if real_prim.has(sp.Integral):
        return None

    # Verify derivative algebraically
    try:
        diff_check = p_simplify(sp.diff(real_prim, variable) - expression, var=variable)
        if sp.simplify(diff_check) != 0:
            return None
    except Exception:
        return None

    t_var = antideriv.t_var
    rat_integrand = sp.cancel(antideriv.lift.rat_num / antideriv.lift.rat_den)
    hermite_part = antideriv.hermite.poly_antideriv + antideriv.hermite.hermite_rat
    sqfree_part = sp.cancel(antideriv.hermite.rem_num / antideriv.hermite.sqfree_den)

    step_meta = {
        "lift_kind": antideriv.lift.kind,
        "t_var": t_var,
        "forward_sym": antideriv.lift.forward_sym,
        "rat_integrand": rat_integrand,
        "hermite_part": hermite_part,
        "sqfree_part": sqfree_part,
        "poles_residues": [
            (pr.pole_sym, pr.residue_sym) for pr in antideriv.residues
        ],
        "projective_form": antideriv.to_projective_sympy(),
        "real_form": real_prim,
    }
    note = (
        f"Integrated via projective calculus over Q(i, sqrt({d - 1})) on P^1: "
        "Ostrowski-Hermite rational reduction + Lagrange-Rothstein-Trager conjugate pole pairing."
    )
    return real_prim, [], note, step_meta


def projective_sympy_definite(
    expression: sp.Expr,
    variable: sp.Symbol,
    lower: Any,
    upper: Any,
    d: int = 3,
    digits: int = 50,
) -> tuple[sp.Expr, str, dict[str, Any]] | None:
    """Compute an exact or arbitrary-precision definite/improper integral on P^1(C)
    (including Fresnel-Airy sector rotations zeta_{2n} in Q(i, sqrt(d-1))).
    """
    try:
        res = p_integrate(expression, (variable, lower, upper), d=d, digits=digits)
    except Exception:
        return None

    assert isinstance(res, DefiniteIntegralResult)
    if res.exact_sympy is not None and not res.exact_sympy.has(sp.Integral):
        note = (
            f"Evaluated on projective line P^1(Q(i, sqrt({d - 1}))) via residue collapse / "
            "biquadratic sector rotation."
        )
        return res.exact_sympy, note, {"numeric_C": res.numeric_value}

    return None


def install_sympy_kernel(d: int = 3) -> None:
    """Monkey-patch `sympy.integrate` and `sympy.trigsimp` to use the Projective CAS kernel first."""
    global _KERNEL_INSTALLED
    if _KERNEL_INSTALLED:
        return

    def _patched_integrate(*args: Any, **kwargs: Any) -> Any:
        if len(args) == 2 and not kwargs:
            expr = sp.sympify(args[0])
            spec = args[1]
            if isinstance(spec, sp.Symbol):
                proj_res = projective_sympy_primitive(expr, spec, d=d)
                if proj_res is not None:
                    return proj_res[0]
            elif isinstance(spec, (tuple, list)) and len(spec) == 3 and isinstance(spec[0], sp.Symbol):
                proj_def = projective_sympy_definite(expr, spec[0], spec[1], spec[2], d=d)
                if proj_def is not None:
                    return proj_def[0]
        return _ORIGINAL_SYMPY_INTEGRATE(*args, **kwargs)

    def _patched_trigsimp(expr: Any, **kwargs: Any) -> Any:
        if not kwargs:
            try:
                return p_simplify(expr)
            except Exception:
                pass
        return _ORIGINAL_SYMPY_TRIGSIMP(expr, **kwargs)

    sp.integrate = _patched_integrate
    sp.trigsimp = _patched_trigsimp
    _KERNEL_INSTALLED = True
