"""Comprehensive test suite for Guiders Mathematics (Projective CAS over Q(i, sqrt(d-1)))."""

from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction

import sympy as sp

from guiders_math import (
    ComplexDecimal,
    MoebiusTransform,
    P1Point,
    QQuadExt,
    cross_ratio,
    hermite_reduction,
    p_contour_integrate,
    p_diff,
    p_dsolve,
    p_integrate,
    p_pade,
    p_prove_identity,
    p_simplify,
    projective_atan_R,
    projective_exp_C,
    projective_log_C,
    projective_pi,
    rational_sqrt_decimal,
)


def test_biquadratic_field_and_zeta8() -> None:
    """Verify exact arithmetic in Q(i, sqrt(2)) and 8th root of unity zeta_8 = (1+i)/sqrt(2)."""
    z8 = QQuadExt.zeta8(d=3)
    i_u = QQuadExt.i_unit(d=3)
    one = QQuadExt.one(d=3)

    assert z8 * z8 == i_u
    assert z8**4 == -one
    assert z8**8 == one
    assert z8 * z8.conj() == one

    # Inverse test on generic element 3 + 2*i - 5*sqrt(2) + 7*i*sqrt(2)
    x = QQuadExt(Fraction(3), Fraction(2), Fraction(-5), Fraction(7), d=3)
    assert x * x.inv() == one


def test_riemann_sphere_p1_and_infinity() -> None:
    """Verify homogeneous P^1 coordinates, infinity [1:0], stereographic circle, and cross-ratio."""
    inf_pt = P1Point.infinity(d=3)
    zero_pt = P1Point.from_affine(0, d=3)
    one_pt = P1Point.from_affine(1, d=3)
    minus_one_pt = P1Point.from_affine(-1, d=3)

    # Stereographic circle at t = 0 -> (1, 0), t = 1 -> (0, 1), t = oo -> (-1, 0) without limits!
    assert zero_pt.stereographic_circle() == (QQuadExt.one(), QQuadExt.zero())
    assert one_pt.stereographic_circle() == (QQuadExt.zero(), QQuadExt.one())
    assert inf_pt.stereographic_circle() == (-QQuadExt.one(), QQuadExt.zero())

    # Angle doubling at t = 1 (45 deg) reaches infinity [1 : 0] (90 deg half-angle) smoothly
    assert one_pt.double_angle() == inf_pt

    # Cross-ratio (oo, 0; 1, -1) = -1 on Riemann sphere
    cr = cross_ratio(inf_pt, zero_pt, one_pt, minus_one_pt)
    assert not cr.is_infinity
    assert cr.affine() == -QQuadExt.one()


def test_pure_rational_pi_log_exp_on_R_and_C() -> None:
    """Verify pure-rational Pade (5, 5) + half-angle descent to 55 decimal places on R and C."""
    with localcontext() as ctx:
        ctx.prec = 55
        pi_55 = projective_pi(digits=55)
        ref_pi = Decimal("3.1415926535897932384626433832795028841971693993751058209")
        assert abs(pi_55 - ref_pi) < Decimal("1e-52")

        # Log_P(1 + i) = (1/2)*ln(2) + i*(pi/4)
        z = QQuadExt(Fraction(1), Fraction(1), Fraction(0), Fraction(0), d=3)
        log_z = projective_log_C(z, digits=50)
        ln2_half = projective_log_C(2, digits=52).re / Decimal(2)
        pi_quarter = pi_55 / Decimal(4)
        assert abs(log_z.re - ln2_half) < Decimal("1e-47")
        assert abs(log_z.im - pi_quarter) < Decimal("1e-47")

        # Round-trip exp_P(Log_P(3 + 4i)) == 3 + 4i
        w = ComplexDecimal(Decimal(3), Decimal(4), digits=50)
        rt = projective_exp_C(projective_log_C(w, digits=52), digits=50)
        assert abs(rt.re - Decimal(3)) < Decimal("1e-47")
        assert abs(rt.im - Decimal(4)) < Decimal("1e-47")


def test_improper_rational_integrals_on_P1() -> None:
    """Verify improper rational integrals on [0, oo) and (-oo, oo) via P^1 residue collapse."""
    with localcontext() as ctx:
        ctx.prec = 55
        # 1. int_0^oo dx / (1 + x^4) = pi / (2 * sqrt(2))
        res = p_integrate("1 / (1 + x**4)", ("x", 0, "oo"), digits=50)
        assert res.exact_sympy == sp.pi / (2 * sp.sqrt(2))
        expected = projective_pi(52) / (Decimal(2) * rational_sqrt_decimal(2, 52))
        assert abs(res.eval_R() - expected) < Decimal("1e-47")

        # 2. int_{-oo}^{oo} dx / (1 + x^2)^2 = pi / 2 (Ostrowski-Hermite reduction + P^1 residues)
        res2 = p_integrate("1 / (1 + x**2)**2", ("x", "-oo", "oo"), digits=50)
        assert res2.exact_sympy == sp.pi / 2
        assert abs(res2.eval_R() - projective_pi(52) / Decimal(2)) < Decimal("1e-47")


def test_fresnel_proper_improper_and_complex() -> None:
    """Verify proper, improper, and complex Fresnel integrals on R and C."""
    with localcontext() as ctx:
        ctx.prec = 55
        # 1. Improper Fresnel integrals: int_0^oo sin(x^2) dx = int_0^oo cos(x^2) dx = sqrt(2*pi) / 4
        s_inf = p_integrate("sin(x**2)", ("x", 0, "oo"), digits=50)
        c_inf = p_integrate("cos(x**2)", ("x", 0, "oo"), digits=50)
        assert sp.simplify(s_inf.exact_sympy - sp.sqrt(2 * sp.pi) / 4) == 0
        assert sp.simplify(c_inf.exact_sympy - sp.sqrt(2 * sp.pi) / 4) == 0

        ref_inf = rational_sqrt_decimal(Decimal(2) * projective_pi(52), 52) / Decimal(4)
        assert abs(s_inf.eval_R() - ref_inf) < Decimal("1e-47")
        assert abs(c_inf.eval_R() - ref_inf) < Decimal("1e-47")

        # 2. Proper Fresnel integrals on [0, 1]:
        # S(1) = int_0^1 sin(x^2) dx = sum_{k=0}^oo (-1)^k / ((2k+1)! * (4k+3)) = 0.310268301723381101808152423165...
        # C(1) = int_0^1 cos(x^2) dx = sum_{k=0}^oo (-1)^k / ((2k)! * (4k+1))   = 0.904524237900272081474788366833...
        s_1 = p_integrate("sin(x**2)", ("x", 0, 1), digits=35)
        c_1 = p_integrate("cos(x**2)", ("x", 0, 1), digits=35)
        ref_s1 = Decimal("0.31026830172338110180815242316541880")
        ref_c1 = Decimal("0.90452423790027208147478836683255715")
        assert abs(s_1.eval_R() - ref_s1) < Decimal("1e-28")
        assert abs(c_1.eval_R() - ref_c1) < Decimal("1e-28")

        # 3. Complex Fresnel integral in C: int_0^{1+I} sin(z^2) dz
        s_c = p_integrate("sin(z**2)", ("z", 0, "1 + I"), digits=30)
        val_c = s_c.eval_C()
        z0 = ComplexDecimal(Decimal(1), Decimal(1), digits=35)
        series_acc = ComplexDecimal(Decimal(0), Decimal(0), digits=35)
        fact = 1
        for k in range(25):
            if k > 0:
                fact *= (2 * k) * (2 * k + 1)
            p_exp = 4 * k + 3
            z_pow = ComplexDecimal(Decimal(1), Decimal(0), digits=35)
            for _ in range(p_exp):
                z_pow = z_pow * z0
            denom = Decimal(fact * p_exp) * (Decimal(-1) if k % 2 == 1 else Decimal(1))
            series_acc = series_acc + (z_pow / ComplexDecimal(denom, Decimal(0), digits=35))
        assert abs(val_c.re - series_acc.re) < Decimal("1e-25")
        assert abs(val_c.im - series_acc.im) < Decimal("1e-25")


def test_complex_path_and_contour_integrals_in_C() -> None:
    """Verify complex path integrals and closed Riemann sphere contour integrals without branch-cut bugs."""
    with localcontext() as ctx:
        ctx.prec = 52
        # 1. Complex path integral: int_{1+i}^{2+3i} dz / (1 + z^2)
        res_c = p_integrate("1 / (1 + z**2)", ("z", "1 + I", "2 + 3*I"), digits=45)
        val = res_c.eval_C()
        # By projective identity: int dz/(1+z^2) = (1/(2i)) * Log_P((z-i)/(z+i)) = (-i/2) * Log_P((z-i)/(z+i))
        z1 = ComplexDecimal(Decimal(1), Decimal(1), digits=48)
        z2 = ComplexDecimal(Decimal(2), Decimal(3), digits=48)
        i_c = ComplexDecimal(Decimal(0), Decimal(1), digits=48)
        ratio = ((z2 - i_c) / (z2 + i_c)) / ((z1 - i_c) / (z1 + i_c))
        exact_c = (ComplexDecimal(Decimal(0), Decimal("-0.5"), digits=48)) * projective_log_C(ratio, digits=48)
        assert abs(val.re - exact_c.re) < Decimal("1e-42")
        assert abs(val.im - exact_c.im) < Decimal("1e-42")

        # 2. Closed contour integral enclosing two upper-half-plane poles of 1/(1 + z^4):
        # Circle |z - i| = 1 encloses alpha_1 = (1+i)/sqrt(2) and alpha_2 = (-1+i)/sqrt(2)
        sym_cnt, num_cnt = p_contour_integrate("1 / (1 + z**4)", "z", center="I", radius=1, digits=45)
        assert sp.simplify(sym_cnt - sp.pi / sp.sqrt(2)) == 0
        ref_cnt = projective_pi(48) / rational_sqrt_decimal(2, 48)
        assert abs(num_cnt.re - ref_cnt) < Decimal("1e-42")
        assert abs(num_cnt.im) < Decimal("1e-42")


def test_cavalieri_dual_diff_pade_and_identities() -> None:
    """Verify Cavalieri-Lagrange dual-ring derivative, Pade [m/n], and P^1 identity prover."""
    x = sp.Symbol("x")
    # 1. Dual-ring derivative of rational and trigonometric functions
    f_rat = (x**3 - 2 * x + 1) / (x**2 + 1)
    assert sp.simplify(p_diff(f_rat, x) - sp.diff(f_rat, x)) == 0

    f_trig = sp.sin(2 * x) / (1 + sp.cos(x) ** 2)
    assert p_prove_identity(p_diff(f_trig, x), sp.diff(f_trig, x), var=x)

    # 2. Pade [2/2] of exp(x) = (12 + 6x + x^2) / (12 - 6x + x^2)
    pade_22 = p_pade(sp.exp(x), x, 2, 2)
    expected_pade = (12 + 6 * x + x**2) / (12 - 6 * x + x**2)
    assert sp.simplify(pade_22 - expected_pade) == 0

    # 3. Projective trigonometric identity proof on P^1
    assert p_prove_identity("sin(3*x)", "3*sin(x) - 4*sin(x)**3")
    assert p_prove_identity("cos(4*x)", "8*cos(x)**4 - 8*cos(x)**2 + 1")


def test_quintic_aberth_and_cubic_airy_fresnel() -> None:
    """Verify irreducible degree-5 rational quadrature and cubic Airy-Fresnel improper integral."""
    with localcontext() as ctx:
        ctx.prec = 45
        # 1. Irreducible degree-5 factor int_0^1 dx / (x^5 + x + 2)
        res_q = p_integrate("1 / (x**5 + x + 2)", ("x", 0, 1), digits=35)
        val_q = res_q.eval_R()
        # Compare with 2-panel 24-point Gauss-Legendre reference on [0, 1/2] and [1/2, 1]
        from guiders_math.quadrature import _gauss_legendre_nodes_weights

        gl = _gauss_legendre_nodes_weights(24, 42)
        ref_q = Decimal(0)
        for panel_start in (Decimal(0), Decimal("0.5")):
            for xi, wi in gl:
                x_pt = panel_start + (xi + Decimal(1)) / Decimal(4)
                ref_q += (wi / Decimal(4)) / (x_pt**5 + x_pt + Decimal(2))
        assert abs(val_q - ref_q) < Decimal("1e-32")

        # 2. Cubic Airy-Fresnel improper integral: int_0^oo cos(x^3) dx = Gamma(1/3) / (2 * sqrt(3))
        res_airy = p_integrate("cos(x**3)", ("x", 0, "oo"), digits=35)
        assert sp.simplify(res_airy.exact_sympy - sp.gamma(sp.Rational(1, 3)) / (2 * sp.sqrt(3))) == 0
        # Gamma(1/3) / (2 * sqrt(3)) = 0.77334294207798985019610161129521734...
        ref_airy = Decimal(str(sp.N(sp.gamma(sp.Rational(1, 3)) / (2 * sp.sqrt(3)), 40)))
        assert abs(res_airy.eval_R() - ref_airy) < Decimal("1e-32")

