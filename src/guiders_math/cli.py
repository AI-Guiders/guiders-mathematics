"""Command-line interface for Guiders Mathematics (Projective CAS over Q(i, sqrt(d-1)))."""

from __future__ import annotations

import argparse
import json
import sys

import sympy as sp

from guiders_math.calculus import p_diff, p_pade
from guiders_math.functor import p_prove_identity, p_simplify
from guiders_math.poly_engine import eval_sympy_alg_to_complex_decimal
from guiders_math.quadrature import DefiniteIntegralResult, ProjectiveAntiderivative, p_contour_integrate, p_integrate


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="guiders-math",
        description="Projective Computer Algebra System over Q(i, sqrt(d-1)) and Riemann Sphere P^1(C)",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    # 1. integrate
    p_int = sub.add_parser("integrate", help="Indefinite, definite, or improper projective integral on R or C")
    p_int.add_argument("expr", help="Integrand expression, e.g. '1/(1+x^4)' or 'sin(x^2)'")
    p_int.add_argument("--var", default="x", help="Integration variable (default: x)")
    p_int.add_argument("--from", dest="lower", default=None, help="Lower limit (e.g. 0, -oo, '1+I')")
    p_int.add_argument("--to", dest="upper", default=None, help="Upper limit (e.g. 1, oo, '2+3*I')")
    p_int.add_argument("--digits", type=int, default=50, help="Decimal precision on R or C (default: 50)")
    p_int.add_argument("--json", action="store_true", help="Output JSON")

    # 2. contour
    p_cnt = sub.add_parser("contour", help="Closed contour integral oint_{|z-center|=radius} f(z) dz on P^1(C)")
    p_cnt.add_argument("expr", help="Rational integrand f(z)")
    p_cnt.add_argument("--var", default="z", help="Complex variable (default: z)")
    p_cnt.add_argument("--center", default="0", help="Contour center in C (default: 0)")
    p_cnt.add_argument("--radius", default="1", help="Contour radius > 0 (default: 1)")
    p_cnt.add_argument("--digits", type=int, default=50, help="Decimal precision (default: 50)")

    # 3. eval
    p_ev = sub.add_parser("eval", help="Pure-rational projective evaluation of an expression on R or C")
    p_ev.add_argument("expr", help="Symbolic expression, e.g. 'pi/(2*sqrt(2))' or 'log(1+I)'")
    p_ev.add_argument("--digits", type=int, default=50, help="Decimal precision (default: 50)")

    # 4. simplify
    p_sim = sub.add_parser("simplify", help="Projective P^1 simplification and identity prover")
    p_sim.add_argument("expr", help="Expression to simplify")
    p_sim.add_argument("--equals", default=None, help="Optional RHS expression to prove identity lhs == rhs")

    # 5. diff
    p_df = sub.add_parser("diff", help="Cavalieri-Lagrange dual-ring projective derivative (no limits)")
    p_df.add_argument("expr", help="Expression to differentiate")
    p_df.add_argument("--var", default="x", help="Variable (default: x)")
    p_df.add_argument("--order", type=int, default=1, help="Derivative order (default: 1)")

    # 6. pade
    p_pd = sub.add_parser("pade", help="Exact projective Pade rational approximant [m/n]")
    p_pd.add_argument("expr", help="Expression to approximate")
    p_pd.add_argument("-m", type=int, default=3, help="Numerator degree m")
    p_pd.add_argument("-n", type=int, default=3, help="Denominator degree n")
    p_pd.add_argument("--var", default="x", help="Variable (default: x)")

    args = parser.parse_args(argv)

    if args.cmd == "integrate":
        if args.lower is not None and args.upper is not None:
            res = p_integrate(args.expr, (args.var, args.lower, args.upper), digits=args.digits)
            assert isinstance(res, DefiniteIntegralResult)
            if args.json:
                print(
                    json.dumps(
                        {
                            "integrand": str(res.integrand),
                            "var": str(res.var),
                            "lower": str(res.lower),
                            "upper": str(res.upper),
                            "exact": str(res.exact_sympy) if res.exact_sympy is not None else None,
                            "numeric_C": res.numeric_value.format(args.digits),
                            "is_real": res.numeric_value.is_real(),
                        },
                        indent=2,
                    )
                )
            else:
                if res.exact_sympy is not None:
                    print(f"Exact (Q(i, sqrt(d-1))): {res.exact_sympy}")
                print(f"Value ({args.digits} digits):   {res.numeric_value.format(args.digits)}")
        else:
            antideriv = p_integrate(args.expr, args.var, digits=args.digits)
            assert isinstance(antideriv, ProjectiveAntiderivative)
            if args.json:
                print(
                    json.dumps(
                        {
                            "projective_form": str(antideriv.to_projective_sympy()),
                            "real_form": str(antideriv.to_real_sympy()),
                            "latex": antideriv.to_latex(real_form=True),
                        },
                        indent=2,
                    )
                )
            else:
                print(f"Projective Q(i, sqrt(d-1)) form: {antideriv.to_projective_sympy()}")
                print(f"Real R form:                     {antideriv.to_real_sympy()}")
                print(f"LaTeX:                           {antideriv.to_latex(real_form=True)}")
        return 0

    if args.cmd == "contour":
        sym_ans, num_ans = p_contour_integrate(
            args.expr, args.var, center=args.center, radius=args.radius, digits=args.digits
        )
        if sym_ans is not None:
            print(f"Exact residue sum:     {sym_ans}")
        print(f"Value ({args.digits} digits): {num_ans.format(args.digits)}")
        return 0

    if args.cmd == "eval":
        val = eval_sympy_alg_to_complex_decimal(sp.sympify(args.expr), digits=args.digits)
        print(val.format(args.digits))
        return 0

    if args.cmd == "simplify":
        if args.equals is not None:
            ok = p_prove_identity(args.expr, args.equals)
            print(f"Identity holds: {ok}")
            return 0 if ok else 1
        print(p_simplify(args.expr))
        return 0

    if args.cmd == "diff":
        print(p_diff(args.expr, args.var, order=args.order))
        return 0

    if args.cmd == "pade":
        print(p_pade(args.expr, args.var, args.m, args.n))
        return 0

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
