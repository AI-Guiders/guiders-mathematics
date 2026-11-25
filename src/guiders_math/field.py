"""Exact arithmetic over the biquadratic field Q(i, sqrt(d-1)) and Riemann sphere P^1.

In projective calculus over Q(i, sqrt(d-1)):
  - Every element x in Q(i, sqrt(m)) (with m = d - 1, default d = 3 -> m = 2)
    is uniquely represented as:
        x = a + b*i + c*sqrt(m) + e*i*sqrt(m),    a, b, c, e in Q.
  - Projective line / Riemann sphere P^1(Q(i, sqrt(m))) consists of homogeneous pairs
    [U : V] with (U, V) != (0, 0). The point at infinity is the regular element
    [1 : 0], eliminating division-by-zero exceptions at infinity.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Union

import sympy as sp

RationalLike = Union[int, Fraction, "QQuadExt"]


def _to_frac(val: int | Fraction | sp.Basic) -> Fraction:
    if isinstance(val, Fraction):
        return val
    if isinstance(val, int):
        return Fraction(val, 1)
    if isinstance(val, sp.Rational):
        return Fraction(int(val.p), int(val.q))
    if isinstance(val, sp.Integer):
        return Fraction(int(val), 1)
    raise TypeError(f"Expected exact rational (int or Fraction), got {type(val)!r}")


@dataclass(frozen=True)
class QQuadExt:
    """Exact element of Q(i, sqrt(m)) where m = d - 1 (default d = 3, m = 2).

    Represents:
        x = a + b*i + c*sqrt(m) + e*i*sqrt(m)
    with exact rational coordinates a, b, c, e in Q.
    """

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)
    c: Fraction = Fraction(0)
    e: Fraction = Fraction(0)
    d: int = 3

    def __post_init__(self) -> None:
        if self.d < 2:
            raise ValueError(f"Dimension parameter d must be >= 2, got {self.d}")
        object.__setattr__(self, "a", _to_frac(self.a))
        object.__setattr__(self, "b", _to_frac(self.b))
        object.__setattr__(self, "c", _to_frac(self.c))
        object.__setattr__(self, "e", _to_frac(self.e))

    @property
    def m(self) -> int:
        """Radicand m = d - 1 of the quadratic extension Q(sqrt(d-1))."""
        return self.d - 1

    @classmethod
    def zero(cls, d: int = 3) -> QQuadExt:
        return cls(Fraction(0), Fraction(0), Fraction(0), Fraction(0), d=d)

    @classmethod
    def one(cls, d: int = 3) -> QQuadExt:
        return cls(Fraction(1), Fraction(0), Fraction(0), Fraction(0), d=d)

    @classmethod
    def i_unit(cls, d: int = 3) -> QQuadExt:
        return cls(Fraction(0), Fraction(1), Fraction(0), Fraction(0), d=d)

    @classmethod
    def sqrt_m(cls, d: int = 3) -> QQuadExt:
        return cls(Fraction(0), Fraction(0), Fraction(1), Fraction(0), d=d)

    @classmethod
    def zeta8(cls, d: int = 3) -> QQuadExt:
        """Eighth root of unity zeta_8 = (1 + i)/sqrt(2) = sqrt(2)/2 + i*sqrt(2)/2 in Q(i, sqrt(2))."""
        if d != 3:
            raise ValueError("zeta8 = (1+i)/sqrt(2) lives in d=3 (Q(i, sqrt(2))).")
        half = Fraction(1, 2)
        return cls(Fraction(0), Fraction(0), half, half, d=3)

    @classmethod
    def coerce(cls, val: RationalLike, d: int = 3) -> QQuadExt:
        if isinstance(val, QQuadExt):
            if val.d != d:
                if val.c == 0 and val.e == 0:
                    return cls(val.a, val.b, Fraction(0), Fraction(0), d=d)
                raise ValueError(f"Cannot mix QQuadExt with d={val.d} and d={d}")
            return val
        return cls(_to_frac(val), Fraction(0), Fraction(0), Fraction(0), d=d)

    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0 and self.c == 0 and self.e == 0

    def is_rational(self) -> bool:
        return self.b == 0 and self.c == 0 and self.e == 0

    def is_real(self) -> bool:
        return self.b == 0 and self.e == 0

    def is_gaussian(self) -> bool:
        return self.c == 0 and self.e == 0

    def real_part(self) -> QQuadExt:
        """Return Re(x) = a + c*sqrt(m) in Q(sqrt(m))."""
        return QQuadExt(self.a, Fraction(0), self.c, Fraction(0), d=self.d)

    def imag_part(self) -> QQuadExt:
        """Return Im(x) = b + e*sqrt(m) in Q(sqrt(m))."""
        return QQuadExt(self.b, Fraction(0), self.e, Fraction(0), d=self.d)

    def conj(self) -> QQuadExt:
        """Complex conjugate: i -> -i."""
        return QQuadExt(self.a, -self.b, self.c, -self.e, d=self.d)

    def radical_conj(self) -> QQuadExt:
        """Galois conjugate over Q(i): sqrt(m) -> -sqrt(m)."""
        return QQuadExt(self.a, self.b, -self.c, -self.e, d=self.d)

    def norm_over_q_sqrt_m(self) -> QQuadExt:
        """Return x * conj(x) = Re(x)^2 + Im(x)^2 in Q(sqrt(m))."""
        return self * self.conj()

    def norm_Q(self) -> Fraction:
        """Absolute field norm N_{Q(i,sqrt(m))/Q}(x) in Q."""
        n1 = self.norm_over_q_sqrt_m()
        # n1 is in Q(sqrt(m)): n1 = A + C*sqrt(m), its norm in Q is A^2 - m*C^2
        return n1.a * n1.a - Fraction(self.m) * n1.c * n1.c

    def __neg__(self) -> QQuadExt:
        return QQuadExt(-self.a, -self.b, -self.c, -self.e, d=self.d)

    def __add__(self, other: RationalLike) -> QQuadExt:
        o = QQuadExt.coerce(other, d=self.d)
        return QQuadExt(self.a + o.a, self.b + o.b, self.c + o.c, self.e + o.e, d=self.d)

    def __radd__(self, other: RationalLike) -> QQuadExt:
        return self.__add__(other)

    def __sub__(self, other: RationalLike) -> QQuadExt:
        o = QQuadExt.coerce(other, d=self.d)
        return QQuadExt(self.a - o.a, self.b - o.b, self.c - o.c, self.e - o.e, d=self.d)

    def __rsub__(self, other: RationalLike) -> QQuadExt:
        return QQuadExt.coerce(other, d=self.d).__sub__(self)

    def __mul__(self, other: RationalLike) -> QQuadExt:
        o = QQuadExt.coerce(other, d=self.d)
        m = Fraction(self.m)
        a1, b1, c1, e1 = self.a, self.b, self.c, self.e
        a2, b2, c2, e2 = o.a, o.b, o.c, o.e
        # Basis: 1, i, s=sqrt(m), i*s with i^2 = -1, s^2 = m, (i*s)^2 = -m
        a = a1 * a2 - b1 * b2 + m * (c1 * c2 - e1 * e2)
        b = a1 * b2 + b1 * a2 + m * (c1 * e2 + e1 * c2)
        c = a1 * c2 - b1 * e2 + c1 * a2 - e1 * b2
        e = a1 * e2 + b1 * c2 + c1 * b2 + e1 * a2
        return QQuadExt(a, b, c, e, d=self.d)

    def __rmul__(self, other: RationalLike) -> QQuadExt:
        return self.__mul__(other)

    def inv(self) -> QQuadExt:
        """Exact multiplicative inverse in Q(i, sqrt(m))."""
        if self.is_zero():
            raise ZeroDivisionError("Cannot invert 0 in Q(i, sqrt(d-1)).")
        # Step 1: multiply by complex conjugate conj(self) -> n1 in Q(sqrt(m))
        c_conj = self.conj()
        n1 = self * c_conj  # n1 = A + C*sqrt(m) (with b=0, e=0)
        # Step 2: multiply by radical conjugate n1.radical_conj() -> norm in Q
        r_conj = n1.radical_conj()
        denom_q = n1.a * n1.a - Fraction(self.m) * n1.c * n1.c
        if denom_q == 0:
            raise ZeroDivisionError(f"Singular norm in Q(i, sqrt({self.m})) for {self}")
        num = c_conj * r_conj
        inv_d = Fraction(1, 1) / denom_q
        return QQuadExt(num.a * inv_d, num.b * inv_d, num.c * inv_d, num.e * inv_d, d=self.d)

    def __truediv__(self, other: RationalLike) -> QQuadExt:
        o = QQuadExt.coerce(other, d=self.d)
        return self * o.inv()

    def __rtruediv__(self, other: RationalLike) -> QQuadExt:
        return QQuadExt.coerce(other, d=self.d) / self

    def __pow__(self, exp: int) -> QQuadExt:
        if exp == 0:
            return QQuadExt.one(d=self.d)
        if exp < 0:
            return self.inv() ** (-exp)
        result = QQuadExt.one(d=self.d)
        base = self
        n = exp
        while n > 0:
            if n & 1:
                result = result * base
            base = base * base
            n >>= 1
        return result

    def to_sympy(self) -> sp.Expr:
        """Convert to an exact SymPy expression in Q(I, sqrt(d-1))."""
        rad = sp.sqrt(sp.Integer(self.m))
        return (
            sp.Rational(self.a.numerator, self.a.denominator)
            + sp.Rational(self.b.numerator, self.b.denominator) * sp.I
            + sp.Rational(self.c.numerator, self.c.denominator) * rad
            + sp.Rational(self.e.numerator, self.e.denominator) * sp.I * rad
        )

    @classmethod
    def from_sympy(cls, expr: sp.Expr, d: int = 3) -> QQuadExt:
        """Extract exact QQuadExt coordinates from a SymPy expression in Q(i, sqrt(d-1))."""
        m = d - 1
        rad = sp.sqrt(sp.Integer(m))
        expanded = sp.expand(sp.sympify(expr))
        re_expr, im_expr = expanded.as_real_imag()

        def split_rad(part: sp.Expr) -> tuple[Fraction, Fraction]:
            part = sp.expand(part)
            c_coeff = part.coeff(rad)
            a_coeff = sp.simplify(part - c_coeff * rad)
            if a_coeff.has(rad) or not a_coeff.is_rational or not c_coeff.is_rational:
                raise ValueError(f"Expression {expr!r} is not in Q(i, sqrt({m}))")
            return _to_frac(sp.Rational(a_coeff)), _to_frac(sp.Rational(c_coeff))

        a_val, c_val = split_rad(re_expr)
        b_val, e_val = split_rad(im_expr)
        return cls(a_val, b_val, c_val, e_val, d=d)

    def to_latex(self) -> str:
        return sp.latex(self.to_sympy())

    def __str__(self) -> str:
        return str(self.to_sympy())


@dataclass(frozen=True)
class P1Point:
    """Point [U : V] on the projective line / Riemann sphere P^1(Q(i, sqrt(d-1))).

    - Affine point z corresponds to [z : 1].
    - Point at infinity corresponds to [1 : 0].
    - Canonical normalization: if V != 0, normalized to [U/V : 1];
      if V == 0, normalized to [1 : 0].
    """

    u: QQuadExt
    v: QQuadExt

    def __post_init__(self) -> None:
        if self.u.is_zero() and self.v.is_zero():
            raise ValueError("Homogeneous coordinates [0 : 0] are invalid in P^1.")
        if self.v.is_zero():
            object.__setattr__(self, "u", QQuadExt.one(d=self.u.d))
            object.__setattr__(self, "v", QQuadExt.zero(d=self.u.d))
        else:
            object.__setattr__(self, "u", self.u / self.v)
            object.__setattr__(self, "v", QQuadExt.one(d=self.u.d))

    @property
    def d(self) -> int:
        return self.u.d

    @classmethod
    def from_affine(cls, z: RationalLike, d: int = 3) -> P1Point:
        return cls(QQuadExt.coerce(z, d=d), QQuadExt.one(d=d))

    @classmethod
    def infinity(cls, d: int = 3) -> P1Point:
        return cls(QQuadExt.one(d=d), QQuadExt.zero(d=d))

    @property
    def is_infinity(self) -> bool:
        return self.v.is_zero()

    def affine(self) -> QQuadExt:
        if self.is_infinity:
            raise ZeroDivisionError("Point [1 : 0] is at infinity on P^1.")
        return self.u

    def reciprocal(self) -> P1Point:
        """Projective inversion [U : V] -> [V : U] (swaps 0 and infinity)."""
        return P1Point(self.v, self.u)

    def stereographic_circle(self) -> tuple[QQuadExt, QQuadExt]:
        """Exact projective circle parameterization (C(t), S(t)) on P^1 -> S^1.

        In homogeneous coordinates t = [U : V]:
            C([U : V]) = (V^2 - U^2) / (V^2 + U^2),
            S([U : V]) = (2 * U * V) / (V^2 + U^2).
        Evaluates at infinity [1 : 0] directly to (-1, 0) without any limit!
        """
        u2 = self.u * self.u
        v2 = self.v * self.v
        denom = v2 + u2
        if denom.is_zero():
            raise ZeroDivisionError("Isotropic pole U^2 + V^2 = 0 on complex circle.")
        c_val = (v2 - u2) / denom
        s_val = (QQuadExt.coerce(2, d=self.d) * self.u * self.v) / denom
        return c_val, s_val

    def stereographic_hyperbola(self) -> tuple[QQuadExt, QQuadExt]:
        """Exact projective Lorentz parameterization (cosh_P(t), sinh_P(t)) on P^1.

        In homogeneous coordinates t = [U : V]:
            cosh_P([U : V]) = (V^2 + U^2) / (V^2 - U^2),
            sinh_P([U : V]) = (2 * U * V) / (V^2 - U^2).
        """
        u2 = self.u * self.u
        v2 = self.v * self.v
        denom = v2 - u2
        if denom.is_zero():
            raise ZeroDivisionError("Light-cone pole V^2 - U^2 = 0 on Lorentz hyperbola.")
        ch_val = (v2 + u2) / denom
        sh_val = (QQuadExt.coerce(2, d=self.d) * self.u * self.v) / denom
        return ch_val, sh_val

    def double_angle(self) -> P1Point:
        """Projective angle doubling [U : V] -> [2*U*V : V^2 - U^2].

        Maps t = [1 : 1] (45 deg) to infinity [1 : 0] (90 deg half-angle / 180 deg) smoothly!
        """
        two = QQuadExt.coerce(2, d=self.d)
        return P1Point(two * self.u * self.v, self.v * self.v - self.u * self.u)

    def cayley_to_disk(self) -> P1Point:
        """Cayley map z -> w = (z - 1)/(z + 1), or in homogeneous coords [U - V : U + V].

        Maps z = 1 -> 0, z = infinity -> 1, z = 0 -> -1, z = -1 -> infinity.
        """
        return P1Point(self.u - self.v, self.u + self.v)

    def cayley_from_disk(self) -> P1Point:
        """Inverse Cayley map w -> z = (1 + w)/(1 - w), or [V + U : V - U]."""
        return P1Point(self.v + self.u, self.v - self.u)

    def __str__(self) -> str:
        if self.is_infinity:
            return "[1 : 0] (oo)"
        return f"[{self.u} : 1]"


@dataclass(frozen=True)
class MoebiusTransform:
    """Fractional linear transformation M = [[a, b], [c, d]] in PGL_2(Q(i, sqrt(d-1))).

    Acts on P^1 via M([U : V]) = [a*U + b*V : c*U + d*V].
    """

    a: QQuadExt
    b: QQuadExt
    c: QQuadExt
    d_elem: QQuadExt

    @property
    def d(self) -> int:
        return self.a.d

    def det(self) -> QQuadExt:
        return self.a * self.d_elem - self.b * self.c

    def __post_init__(self) -> None:
        if self.det().is_zero():
            raise ValueError("Moebius transformation must have non-zero determinant.")

    def apply(self, pt: P1Point) -> P1Point:
        new_u = self.a * pt.u + self.b * pt.v
        new_v = self.c * pt.u + self.d_elem * pt.v
        return P1Point(new_u, new_v)

    def compose(self, other: MoebiusTransform) -> MoebiusTransform:
        return MoebiusTransform(
            a=self.a * other.a + self.b * other.c,
            b=self.a * other.b + self.b * other.d_elem,
            c=self.c * other.a + self.d_elem * other.c,
            d_elem=self.c * other.b + self.d_elem * other.d_elem,
        )

    def inverse(self) -> MoebiusTransform:
        return MoebiusTransform(
            a=self.d_elem,
            b=-self.b,
            c=-self.c,
            d_elem=self.a,
        )


def cross_ratio(z1: P1Point, z2: P1Point, z3: P1Point, z4: P1Point) -> P1Point:
    """Projective cross-ratio (z1, z2; z3, z4) on the Riemann sphere P^1.

    In homogeneous coordinates z_k = [U_k : V_k], using the determinant
    det(z_j, z_k) = U_j*V_k - V_j*U_k:
        (z1, z2; z3, z4) = [ det(z1, z3)*det(z2, z4) : det(z1, z4)*det(z2, z3) ].
    Works seamlessly when any of z1, z2, z3, z4 is infinity [1 : 0]!
    """
    def det2(p: P1Point, q: P1Point) -> QQuadExt:
        return p.u * q.v - p.v * q.u

    num = det2(z1, z3) * det2(z2, z4)
    den = det2(z1, z4) * det2(z2, z3)
    return P1Point(num, den)
