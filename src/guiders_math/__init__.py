"""Guiders Mathematics — Projective Computer Algebra System over Q(i, sqrt(d-1)) and P^1(C)."""

from __future__ import annotations

from guiders_math.calculus import p_diff, p_dsolve, p_pade
from guiders_math.eval_rc import (
    ComplexDecimal,
    complex_sqrt_heron,
    eval_qquad_to_complex_decimal,
    eval_to_complex_decimal,
    projective_atan_R,
    projective_cos_sin_C,
    projective_exp_C,
    projective_log_C,
    projective_pi,
    rational_sqrt_decimal,
)
from guiders_math.field import MoebiusTransform, P1Point, QQuadExt, cross_ratio
from guiders_math.functor import ProjectiveLift, lift_integrand_to_projective, p_prove_identity, p_simplify
from guiders_math.poly_engine import (
    HermiteDecomposition,
    PoleResidue,
    eval_sympy_alg_to_complex_decimal,
    hermite_reduction,
    lagrange_residues,
)
from guiders_math.quadrature import (
    DefiniteIntegralResult,
    ProjectiveAntiderivative,
    p_antiderivative,
    p_contour_integrate,
    p_integrate,
)
from guiders_math.sympy_bridge import (
    install_sympy_kernel,
    projective_sympy_definite,
    projective_sympy_primitive,
)

__all__ = [
    "QQuadExt",
    "P1Point",
    "MoebiusTransform",
    "cross_ratio",
    "ComplexDecimal",
    "rational_sqrt_decimal",
    "complex_sqrt_heron",
    "projective_atan_R",
    "projective_pi",
    "projective_log_C",
    "projective_exp_C",
    "projective_cos_sin_C",
    "eval_qquad_to_complex_decimal",
    "eval_to_complex_decimal",
    "eval_sympy_alg_to_complex_decimal",
    "PoleResidue",
    "HermiteDecomposition",
    "hermite_reduction",
    "lagrange_residues",
    "ProjectiveLift",
    "lift_integrand_to_projective",
    "p_simplify",
    "p_prove_identity",
    "ProjectiveAntiderivative",
    "DefiniteIntegralResult",
    "p_antiderivative",
    "p_integrate",
    "p_contour_integrate",
    "p_diff",
    "p_pade",
    "p_dsolve",
    "projective_sympy_primitive",
    "projective_sympy_definite",
    "install_sympy_kernel",
]
