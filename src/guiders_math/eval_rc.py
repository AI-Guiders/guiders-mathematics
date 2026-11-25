"""Pure-rational projective evaluator on R and C to arbitrary precision (N digits).

ZERO usage of math.sin, math.cos, math.atan, math.log, math.exp, or mpmath!
Every transcendental value on R and C is computed strictly from:
  1. Quadratic Heron-Newton rational iteration for square roots in Q and Q(i),
  2. Projective half-angle / square-root descent z -> sqrt(z) on P^1(Q(i)),
  3. Cayley transform w = (z - 1)/(z + 1) and diagonal (5, 5) Pade rational approximant:
         Log_P(z_m) = 2*w * (945 - 735*w^2 + 64*w^4) / (945 - 1050*w^2 + 225*w^4),
     which achieves O(w^11) = O(2^{-11*m}) accuracy in m half-angle steps.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
import math
from typing import Union

from guiders_math.field import QQuadExt


def _decimal_from_scaled_int(val: int, scale_exp: int, digits: int) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = max(digits + 25, 50)
        d = Decimal(val) / (Decimal(10) ** scale_exp)
        # Round cleanly to `digits` significant/decimal places
    with localcontext() as ctx:
        ctx.prec = digits
        return +d


@dataclass(frozen=True)
class ComplexDecimal:
    """Arbitrary-precision complex number re + i*im on C produced by projective descent."""

    re: Decimal
    im: Decimal
    digits: int = 50

    def is_real(self, tol_exp: int | None = None) -> bool:
        thr = Decimal(10) ** (-(tol_exp if tol_exp is not None else max(self.digits - 4, 10)))
        return abs(self.im) <= thr

    def to_real(self) -> Decimal:
        return self.re

    def __add__(self, other: ComplexDecimal) -> ComplexDecimal:
        prec = max(self.digits, other.digits)
        with localcontext() as ctx:
            ctx.prec = prec + 10
            re = self.re + other.re
            im = self.im + other.im
        with localcontext() as ctx:
            ctx.prec = prec
            return ComplexDecimal(+re, +im, digits=prec)

    def __sub__(self, other: ComplexDecimal) -> ComplexDecimal:
        prec = max(self.digits, other.digits)
        with localcontext() as ctx:
            ctx.prec = prec + 10
            re = self.re - other.re
            im = self.im - other.im
        with localcontext() as ctx:
            ctx.prec = prec
            return ComplexDecimal(+re, +im, digits=prec)

    def __mul__(self, other: ComplexDecimal) -> ComplexDecimal:
        prec = max(self.digits, other.digits)
        with localcontext() as ctx:
            ctx.prec = prec + 10
            re = self.re * other.re - self.im * other.im
            im = self.re * other.im + self.im * other.re
        with localcontext() as ctx:
            ctx.prec = prec
            return ComplexDecimal(+re, +im, digits=prec)

    def __truediv__(self, other: ComplexDecimal) -> ComplexDecimal:
        prec = max(self.digits, other.digits)
        with localcontext() as ctx:
            ctx.prec = prec + 15
            denom = other.re * other.re + other.im * other.im
            if denom == 0:
                raise ZeroDivisionError("ComplexDecimal division by zero.")
            re = (self.re * other.re + self.im * other.im) / denom
            im = (self.im * other.re - self.re * other.im) / denom
        with localcontext() as ctx:
            ctx.prec = prec
            return ComplexDecimal(+re, +im, digits=prec)

    def conj(self) -> ComplexDecimal:
        return ComplexDecimal(self.re, -self.im, digits=self.digits)

    def format(self, digits: int | None = None) -> str:
        d = digits if digits is not None else self.digits
        with localcontext() as ctx:
            ctx.prec = d
            r = +self.re
            i = +self.im
            thr = Decimal(10) ** (-(d - 2))
            if abs(i) <= thr:
                return f"{r:f}" if abs(r) >= Decimal("1e-4") or r == 0 else str(r)
            sign = "+" if i >= 0 else "-"
            abs_i = abs(i)
            return f"{r} {sign} {abs_i}*i"

    def __str__(self) -> str:
        return self.format()


def rational_sqrt_decimal(x: Decimal | int | Fraction, digits: int = 50) -> Decimal:
    """Compute sqrt(x) on R>=0 via pure Heron-Newton rational iteration r_{k+1} = (r_k + x/r_k)/2."""
    work_prec = digits + 25
    with localcontext() as ctx:
        ctx.prec = work_prec
        if isinstance(x, Fraction):
            val = Decimal(x.numerator) / Decimal(x.denominator)
        else:
            val = Decimal(x)
        if val < 0:
            raise ValueError(f"rational_sqrt_decimal requires x >= 0, got {val}")
        if val == 0:
            return Decimal(0)

        # Initial rational seed using exponent halving
        exp10 = val.adjusted()
        half_exp = exp10 // 2
        mant = val.scaleb(-exp10)
        # 3-step integer Heron seed on mantissa in [1, 10)
        r = Decimal(2) * (Decimal(10) ** half_exp)
        r = (r + val / r) / Decimal(2)

        tol = Decimal(10) ** (-(digits + 18))
        max_steps = int(math.ceil(math.log2(work_prec))) + 10
        two = Decimal(2)
        for _ in range(max_steps):
            r_next = (r + val / r) / two
            if abs(r_next - r) <= abs(r_next) * tol:
                r = r_next
                break
            r = r_next

    with localcontext() as ctx:
        ctx.prec = digits
        return +r


def complex_sqrt_heron(re: Decimal, im: Decimal, digits: int = 50) -> tuple[Decimal, Decimal]:
    """Principal square root sqrt(re + i*im) on C using pure algebraic Heron-Newton radicals.

    For z = a + b*i:
        |z| = sqrt(a^2 + b^2),
        u = sqrt((|z| + a) / 2),
        v = sgn(b) * sqrt((|z| - a) / 2).
    To avoid catastrophic cancellation when a < 0 and |b| << |a|, uses the exact relation
    u * v = b / 2.
    """
    work_prec = digits + 25
    with localcontext() as ctx:
        ctx.prec = work_prec
        if im == 0:
            if re >= 0:
                return rational_sqrt_decimal(re, work_prec), Decimal(0)
            return Decimal(0), rational_sqrt_decimal(-re, work_prec)

        r_abs = rational_sqrt_decimal(re * re + im * im, work_prec)
        two = Decimal(2)
        if re >= 0:
            u = rational_sqrt_decimal((r_abs + re) / two, work_prec)
            v = im / (two * u)
        else:
            v_abs = rational_sqrt_decimal((r_abs - re) / two, work_prec)
            v = v_abs if im >= 0 else -v_abs
            u = im / (two * v)
        return u, v


def _c_mul(a: tuple[Decimal, Decimal], b: tuple[Decimal, Decimal]) -> tuple[Decimal, Decimal]:
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def _c_add(a: tuple[Decimal, Decimal], b: tuple[Decimal, Decimal]) -> tuple[Decimal, Decimal]:
    return (a[0] + b[0], a[1] + b[1])


def _c_sub(a: tuple[Decimal, Decimal], b: tuple[Decimal, Decimal]) -> tuple[Decimal, Decimal]:
    return (a[0] - b[0], a[1] - b[1])


def _c_div(a: tuple[Decimal, Decimal], b: tuple[Decimal, Decimal]) -> tuple[Decimal, Decimal]:
    den = b[0] * b[0] + b[1] * b[1]
    if den == 0:
        raise ZeroDivisionError("Complex division by zero in projective Pade evaluator.")
    return ((a[0] * b[0] + a[1] * b[1]) / den, (a[1] * b[0] - a[0] * b[1]) / den)


def _pade55_arctanh_complex(w: tuple[Decimal, Decimal]) -> tuple[Decimal, Decimal]:
    """Diagonal (5, 5) Pade rational approximant for 2 * arctanh(w) = Log_P((1+w)/(1-w)).

    Formula:
        2 * arctanh(w) = 2 * w * (945 - 735*w^2 + 64*w^4) / (945 - 1050*w^2 + 225*w^4) + O(w^11).
    """
    w2 = _c_mul(w, w)
    w4 = _c_mul(w2, w2)

    c945 = (Decimal(945), Decimal(0))
    num_poly = _c_add(
        _c_sub(c945, (Decimal(735) * w2[0], Decimal(735) * w2[1])),
        (Decimal(64) * w4[0], Decimal(64) * w4[1]),
    )
    den_poly = _c_add(
        _c_sub(c945, (Decimal(1050) * w2[0], Decimal(1050) * w2[1])),
        (Decimal(225) * w4[0], Decimal(225) * w4[1]),
    )
    two_w = (Decimal(2) * w[0], Decimal(2) * w[1])
    return _c_div(_c_mul(two_w, num_poly), den_poly)


def projective_atan_R(t: Decimal | int | Fraction, digits: int = 50) -> Decimal:
    """Compute arctan(t) on R via projective half-angle descent + diagonal (5, 5) Pade rational.

    Uses the exact projective half-angle relation:
        t_{k+1} = t_k / (1 + sqrt(1 + t_k^2))
    followed by the Cayley-Pade (5, 5) rational approximant for arctan(t_m):
        arctan(u) = u * (945 + 735*u^2 + 64*u^4) / (945 + 1050*u^2 + 225*u^4) + O(u^11).
    """
    work_prec = digits + 25
    with localcontext() as ctx:
        ctx.prec = work_prec
        if isinstance(t, Fraction):
            u = Decimal(t.numerator) / Decimal(t.denominator)
        else:
            u = Decimal(t)
        if u == 0:
            return Decimal(0)

        sign = Decimal(-1) if u < 0 else Decimal(1)
        u = abs(u)

        # If u > 1, use projective inversion on P^1: arctan(u) = 2*arctan(1) - arctan(1/u)
        # Or simply perform 2 extra half-angle steps since t/(1+sqrt(1+t^2)) < 1 after 1 step for ANY u > 0!
        # Specifically, for ANY u in (0, oo), after 1 step u_1 = u/(1+sqrt(1+u^2)) < 1!
        # To guarantee O(u_m^11) < 10^{-(digits+15)}, we need u_m < 10^{-(digits+15)/10}.
        target_log2 = int(math.ceil((digits + 20) * 3.321928 / 9.5)) + 3
        m_steps = max(target_log2, 10)

        one = Decimal(1)
        for _ in range(m_steps):
            u = u / (one + rational_sqrt_decimal(one + u * u, work_prec))

        u2 = u * u
        u4 = u2 * u2
        num = u * (Decimal(945) + Decimal(735) * u2 + Decimal(64) * u4)
        den = Decimal(945) + Decimal(1050) * u2 + Decimal(225) * u4
        res = sign * (Decimal(2) ** m_steps) * (num / den)

    with localcontext() as ctx:
        ctx.prec = digits
        return +res


def projective_pi(digits: int = 50) -> Decimal:
    """Compute pi to `digits` decimal places via Machin's projective identity in Q(i):
        (5 + i)^4 * (239 - i) = 114244 * (1 + i)
        => pi = 16 * arctan_P(1/5) - 4 * arctan_P(1/239).
    """
    work_prec = digits + 20
    with localcontext() as ctx:
        ctx.prec = work_prec
        a1 = projective_atan_R(Fraction(1, 5), work_prec)
        a2 = projective_atan_R(Fraction(1, 239), work_prec)
        pi_val = Decimal(16) * a1 - Decimal(4) * a2
    with localcontext() as ctx:
        ctx.prec = digits
        return +pi_val


def projective_log_C(z: ComplexDecimal | complex | int | Fraction | QQuadExt, digits: int = 50) -> ComplexDecimal:
    """Compute principal projective logarithm Log_P(z) = ln|z| + i*Arg(z) on C\\{0}.

    Uses m steps of complex square-root descent z_{k+1} = sqrt(z_k) in Q(i),
    followed by Cayley map w_m = (z_m - 1)/(z_m + 1) and diagonal (5, 5) Pade rational:
        Log_P(z) = 2^m * 2*w_m*(945 - 735*w_m^2 + 64*w_m^4)/(945 - 1050*w_m^2 + 225*w_m^4).
    """
    work_prec = digits + 25
    c_val = eval_to_complex_decimal(z, work_prec)
    with localcontext() as ctx:
        ctx.prec = work_prec
        re, im = c_val.re, c_val.im
        if re == 0 and im == 0:
            raise ZeroDivisionError("Log_P(0) is the projective south pole infinity [1 : 0].")

        # Handle negative real axis (branch cut Arg(z) = +pi) exactly:
        if im == 0 and re < 0:
            pi_val = projective_pi(work_prec)
            ln_abs = projective_log_C(ComplexDecimal(-re, Decimal(0), work_prec), work_prec).re
            with localcontext() as out_ctx:
                out_ctx.prec = digits
                return ComplexDecimal(+ln_abs, +pi_val, digits=digits)

        # Normalize modulus by power of 2 to keep z near unit circle before square-root descent
        r2 = re * re + im * im
        # Choose number of complex square-root halving steps m_steps
        target_steps = int(math.ceil((digits + 20) * 3.321928 / 9.5)) + 5
        # Add extra steps if |z| is very large or very small
        mag_exp = abs(int(r2.adjusted())) // 2
        m_steps = max(target_steps, 12) + int(math.ceil(math.log2(mag_exp + 2))) + 2

        cur_re, cur_im = re, im
        for _ in range(m_steps):
            cur_re, cur_im = complex_sqrt_heron(cur_re, cur_im, work_prec)

        # Cayley transform w = (z_m - 1) / (z_m + 1)
        one = (Decimal(1), Decimal(0))
        zm = (cur_re, cur_im)
        w = _c_div(_c_sub(zm, one), _c_add(zm, one))

        log_zm = _pade55_arctanh_complex(w)
        scale = Decimal(2) ** m_steps
        out_re = scale * log_zm[0]
        out_im = scale * log_zm[1]

    with localcontext() as ctx:
        ctx.prec = digits
        return ComplexDecimal(+out_re, +out_im, digits=digits)


def projective_exp_C(z: ComplexDecimal | complex | int | Fraction | QQuadExt, digits: int = 50) -> ComplexDecimal:
    """Compute exp_P(z) on C via argument halving z_m = z / 2^m, diagonal (5, 5) Pade Cayley rational,
    and m squaring steps z_{k+1} = z_k^2.

    Diagonal (5, 5) Pade approximant to exp(u):
        P_5(u) = 30240 + 15120*u + 3360*u^2 + 420*u^3 + 30*u^4 + u^5,
        Q_5(u) = 30240 - 15120*u + 3360*u^2 - 420*u^3 + 30*u^4 - u^5,
        exp(u) = P_5(u) / Q_5(u) + O(u^11).
    """
    work_prec = digits + 25
    c_val = eval_to_complex_decimal(z, work_prec)
    with localcontext() as ctx:
        ctx.prec = work_prec
        re, im = c_val.re, c_val.im
        if re == 0 and im == 0:
            return ComplexDecimal(Decimal(1), Decimal(0), digits=digits)

        mag = max(abs(re), abs(im))
        extra = int(math.ceil(math.log2(float(mag) + 1.0))) if mag > 1 else 0
        target_steps = int(math.ceil((digits + 20) * 3.321928 / 9.5)) + 4 + extra
        m_steps = max(target_steps, 12)

        inv_scale = Decimal(1) / (Decimal(2) ** m_steps)
        u = (re * inv_scale, im * inv_scale)
        u2 = _c_mul(u, u)
        u3 = _c_mul(u2, u)
        u4 = _c_mul(u2, u2)
        u5 = _c_mul(u4, u)

        even_part = (
            Decimal(30240) + Decimal(3360) * u2[0] + Decimal(30) * u4[0],
            Decimal(3360) * u2[1] + Decimal(30) * u4[1],
        )
        odd_part = (
            Decimal(15120) * u[0] + Decimal(420) * u3[0] + u5[0],
            Decimal(15120) * u[1] + Decimal(420) * u3[1] + u5[1],
        )
        cur = _c_div(_c_add(even_part, odd_part), _c_sub(even_part, odd_part))

        for _ in range(m_steps):
            cur = _c_mul(cur, cur)

    with localcontext() as ctx:
        ctx.prec = digits
        return ComplexDecimal(+cur[0], +cur[1], digits=digits)


def projective_cos_sin_C(z: ComplexDecimal | int | Fraction | QQuadExt, digits: int = 50) -> tuple[ComplexDecimal, ComplexDecimal]:
    """Compute (cos_P(z), sin_P(z)) on R or C via projective exp_P(i*z) and exp_P(-i*z)."""
    work_prec = digits + 20
    c_val = eval_to_complex_decimal(z, work_prec)
    with localcontext() as ctx:
        ctx.prec = work_prec
        iz = ComplexDecimal(-c_val.im, c_val.re, digits=work_prec)
        minus_iz = ComplexDecimal(c_val.im, -c_val.re, digits=work_prec)
        e_pos = projective_exp_C(iz, work_prec)
        e_neg = projective_exp_C(minus_iz, work_prec)
        two = ComplexDecimal(Decimal(2), Decimal(0), digits=work_prec)
        two_i = ComplexDecimal(Decimal(0), Decimal(2), digits=work_prec)
        cos_val = (e_pos + e_neg) / two
        sin_val = (e_pos - e_neg) / two_i
    with localcontext() as ctx:
        ctx.prec = digits
        return (
            ComplexDecimal(+cos_val.re, +cos_val.im, digits=digits),
            ComplexDecimal(+sin_val.re, +sin_val.im, digits=digits),
        )


def eval_qquad_to_complex_decimal(elem: QQuadExt, digits: int = 50) -> ComplexDecimal:
    """Evaluate an exact element a + b*i + c*sqrt(m) + e*i*sqrt(m) of Q(i, sqrt(d-1)) on C to `digits` places."""
    work_prec = digits + 20
    with localcontext() as ctx:
        ctx.prec = work_prec
        sqrt_m = rational_sqrt_decimal(elem.m, work_prec) if (elem.c != 0 or elem.e != 0) else Decimal(0)
        a_d = Decimal(elem.a.numerator) / Decimal(elem.a.denominator)
        b_d = Decimal(elem.b.numerator) / Decimal(elem.b.denominator)
        c_d = Decimal(elem.c.numerator) / Decimal(elem.c.denominator)
        e_d = Decimal(elem.e.numerator) / Decimal(elem.e.denominator)
        re = a_d + c_d * sqrt_m
        im = b_d + e_d * sqrt_m
    with localcontext() as ctx:
        ctx.prec = digits
        return ComplexDecimal(+re, +im, digits=digits)


def eval_to_complex_decimal(
    val: ComplexDecimal | QQuadExt | int | Fraction | Decimal | float | complex,
    digits: int = 50,
) -> ComplexDecimal:
    """Coerce any exact or numeric input to a ComplexDecimal with `digits` precision."""
    if isinstance(val, ComplexDecimal):
        if val.digits >= digits:
            with localcontext() as ctx:
                ctx.prec = digits
                return ComplexDecimal(+val.re, +val.im, digits=digits)
        return ComplexDecimal(val.re, val.im, digits=digits)
    if isinstance(val, QQuadExt):
        return eval_qquad_to_complex_decimal(val, digits=digits)
    with localcontext() as ctx:
        ctx.prec = digits + 15
        if isinstance(val, Fraction):
            re = Decimal(val.numerator) / Decimal(val.denominator)
            im = Decimal(0)
        elif isinstance(val, (int, Decimal)):
            re = Decimal(val)
            im = Decimal(0)
        elif isinstance(val, float):
            re = Decimal(str(val))
            im = Decimal(0)
        elif isinstance(val, complex):
            re = Decimal(str(val.real))
            im = Decimal(str(val.imag))
        else:
            raise TypeError(f"Unsupported type for projective evaluation: {type(val)!r}")
    with localcontext() as ctx:
        ctx.prec = digits
        return ComplexDecimal(+re, +im, digits=digits)


def projective_e(digits: int = 50) -> Decimal:
    """Compute Euler's number e = exp_P(1) to `digits` decimal places via dyadic halving + (5, 5) Cayley-Pade."""
    return projective_exp_C(1, digits=digits).re


def projective_ln2(digits: int = 50) -> Decimal:
    """Compute ln(2) = Log_P(2) to `digits` decimal places via Heron square-root descent + Cayley-Pade (5, 5)."""
    return projective_log_C(2, digits=digits).re


def projective_golden_ratio(digits: int = 50) -> Decimal:
    """Compute the golden ratio phi = (1 + sqrt(5)) / 2 via pure Heron-Newton rational iteration."""
    work_prec = digits + 15
    with localcontext() as ctx:
        ctx.prec = work_prec
        val = (Decimal(1) + rational_sqrt_decimal(5, work_prec)) / Decimal(2)
    with localcontext() as ctx:
        ctx.prec = digits
        return +val


def projective_lemniscate_varpi(digits: int = 50) -> Decimal:
    """Compute Gauss's lemniscate constant varpi = 2 * int_0^1 dx / sqrt(1 - x^4) = pi / AGM(1, sqrt(2))
    using quadratic Arithmetic-Geometric Mean descent over Q(sqrt(2)) (d = 3).
    """
    work_prec = digits + 20
    with localcontext() as ctx:
        ctx.prec = work_prec
        pi_val = projective_pi(work_prec)
        a = Decimal(1)
        b = rational_sqrt_decimal(2, work_prec)
        n_steps = int(math.ceil(math.log2(work_prec))) + 4
        two = Decimal(2)
        for _ in range(n_steps):
            a_next = (a + b) / two
            b_next = rational_sqrt_decimal(a * b, work_prec)
            if abs(a_next - b_next) <= Decimal(10) ** (-(digits + 15)):
                a = a_next
                break
            a, b = a_next, b_next
        varpi = pi_val / a
    with localcontext() as ctx:
        ctx.prec = digits
        return +varpi


def projective_gamma_quarter(digits: int = 50) -> Decimal:
    """Compute Gamma(1/4) = sqrt(2 * varpi * sqrt(2 * pi)) purely via Gauss AGM(1, sqrt(2)) and Heron radicals."""
    work_prec = digits + 20
    with localcontext() as ctx:
        ctx.prec = work_prec
        pi_val = projective_pi(work_prec)
        varpi = projective_lemniscate_varpi(work_prec)
        sqrt_2pi = rational_sqrt_decimal(Decimal(2) * pi_val, work_prec)
        g14 = rational_sqrt_decimal(Decimal(2) * varpi * sqrt_2pi, work_prec)
    with localcontext() as ctx:
        ctx.prec = digits
        return +g14


def projective_catalan_G(digits: int = 50) -> Decimal:
    """Compute Catalan's constant G = sum_{k=0}^oo (-1)^k / (2k+1)^2 = int_0^1 arctan(x)/x dx
    via Ramanujan's projective identity over Q(sqrt(3)) (d = 4):
        G = (pi / 8) * Log_P(2 + sqrt(3)) + (3 / 8) * sum_{k=0}^oo (k!)^2 / ((2k+1)! * (2k+1)).
    """
    work_prec = digits + 25
    with localcontext() as ctx:
        ctx.prec = work_prec
        pi_val = projective_pi(work_prec)
        sqrt3 = rational_sqrt_decimal(3, work_prec)
        ln_unit = projective_log_C(Decimal(2) + sqrt3, work_prec).re
        s = Decimal(0)
        term = Decimal(1)
        tol = Decimal(10) ** (-(digits + 20))
        max_k = int(work_prec * 1.8) + 20
        for k in range(max_k):
            contrib = term / Decimal(2 * k + 1)
            s += contrib
            if abs(contrib) <= tol:
                break
            term = term * Decimal(k + 1) / Decimal(2 * (2 * k + 3))
        g_val = (pi_val / Decimal(8)) * ln_unit + (Decimal(3) / Decimal(8)) * s
    with localcontext() as ctx:
        ctx.prec = digits
        return +g_val


def projective_apery_zeta3(digits: int = 50) -> Decimal:
    """Compute Apery's constant zeta(3) = sum_{n=1}^oo 1/n^3 via Amdeberhan-Zeilberger pure-rational
    hypergeometric Pade acceleration (~3.01 decimal digits per step):
        zeta(3) = (1/64) * sum_{k=0}^oo (-1)^k * (205*k^2 + 250*k + 77) * (k!)^10 / ((2k+1)!)^5.
    """
    work_prec = digits + 25
    with localcontext() as ctx:
        ctx.prec = work_prec
        zeta3 = Decimal(0)
        term_num = Decimal(1)
        tol = Decimal(10) ** (-(digits + 20))
        max_k = (work_prec // 3) + 12
        for k in range(max_k):
            poly = Decimal(205 * k * k + 250 * k + 77)
            sign = Decimal(-1) if (k & 1) else Decimal(1)
            contrib = sign * poly * term_num
            zeta3 += contrib
            if abs(contrib) <= tol:
                break
            num_step = Decimal((k + 1) ** 5)
            den_step = Decimal(32) * Decimal((2 * k + 3) ** 5)
            term_num = term_num * num_step / den_step
        zeta3 = zeta3 / Decimal(64)
    with localcontext() as ctx:
        ctx.prec = digits
        return +zeta3


def projective_euler_gamma(digits: int = 50) -> Decimal:
    """Compute the Euler-Mascheroni constant gamma via Brent-McMillan rational Bessel acceleration
    paired with pure-rational projective logarithm Log_P(n).
    """
    work_prec = digits + 25
    with localcontext() as ctx:
        ctx.prec = work_prec
        n = int(digits * 1.15) + 10
        ln_n = projective_log_C(n, work_prec).re
        a_term = -ln_n
        b_term = Decimal(1)
        u_sum = a_term
        v_sum = b_term
        n2 = Decimal(n * n)
        tol = Decimal(10) ** (-(digits + 20))
        for k in range(1, 5 * n + 30):
            k_dec = Decimal(k)
            b_term = b_term * n2 / (k_dec * k_dec)
            a_term = (a_term * n2 / k_dec + b_term) / k_dec
            u_sum += a_term
            v_sum += b_term
            if abs(a_term) < tol and abs(b_term) < tol:
                break
        gamma_val = u_sum / v_sum
    with localcontext() as ctx:
        ctx.prec = digits
        return +gamma_val
