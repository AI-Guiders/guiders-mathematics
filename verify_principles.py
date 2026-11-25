#!/usr/bin/env python3
"""Verification gate for Guiders Mathematics (Projective CAS over Q(i, sqrt(d-1))).

Enforces:
  1. Pure-rational evaluation gate (`eval_rc_purity`): `src/guiders_math/eval_rc.py`
     must NOT import `mpmath` and must NOT call `math.sin`, `math.cos`, `math.atan`,
     `math.log`, or `math.exp` — all evaluations on R and C must go through pure-rational
     Heron-Newton radicals and half-angle Pade (5, 5) descent.
  2. Full pytest verification suite across field Q(i, sqrt(d-1)), Riemann sphere P^1(C),
     Ostrowski-Hermite reduction, Lagrange-Rothstein-Trager quadrature, improper integrals,
     Fresnel-Airy phase integrals, contour integration, Cavalieri dual derivatives, and Pade [m/n].
"""

import os
import subprocess
import sys

if sys.version_info[0] < 3:
    os.execvp("python3", ["python3", os.path.abspath(__file__)] + sys.argv[1:])

import ast
from pathlib import Path

_REPO = Path(__file__).resolve().parent
_SRC = _REPO / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))


def _ensure_venv_python(repo_root):
    """Return a Python executable with `sympy` and `pytest` installed, auto-creating `.venv` if needed."""
    check = subprocess.call(
        [sys.executable, "-c", "import sympy, pytest"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if check == 0:
        return sys.executable
    venv_dir = repo_root / ".venv"
    venv_py = (
        venv_dir / "Scripts" / "python.exe"
        if os.name == "nt"
        else venv_dir / "bin" / "python3"
    )
    req_file = repo_root / "requirements.txt"
    if not venv_py.is_file():
        print("[guiders-mathematics] Creating virtual environment in {} ...".format(venv_dir), file=sys.stderr)
        subprocess.check_call([sys.executable, "-m", "venv", str(venv_dir)])
    if subprocess.call([str(venv_py), "-c", "import sympy, pytest"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
        print("[guiders-mathematics] Installing dependencies from {} into .venv ...".format(req_file.name), file=sys.stderr)
        subprocess.check_call([str(venv_py), "-m", "pip", "install", "--upgrade", "pip"])
        subprocess.check_call([str(venv_py), "-m", "pip", "install", "-r", str(req_file)])
    return str(venv_py)


FORBIDDEN_MATH_FUNCS = {"sin", "cos", "tan", "atan", "atan2", "asin", "acos", "log", "exp", "sqrt"}


def check_eval_rc_purity(eval_rc_path):
    """Verify via AST that eval_rc.py never imports mpmath or calls floating transcendental math.* functions."""
    text = eval_rc_path.read_text(encoding="utf-8")
    tree = ast.parse(text, filename=str(eval_rc_path))
    violations = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            mod = getattr(node, "module", "") or ""
            names = [alias.name for alias in node.names]
            if "mpmath" in mod or "mpmath" in names:
                violations.append("line {}: forbidden import of mpmath in eval_rc.py".format(node.lineno))
        if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name):
            if node.value.id == "math" and node.attr in FORBIDDEN_MATH_FUNCS:
                violations.append("line {}: forbidden call to math.{} in eval_rc.py".format(node.lineno, node.attr))
    return violations


def main():
    eval_rc_file = _REPO / "src" / "guiders_math" / "eval_rc.py"
    violations = check_eval_rc_purity(eval_rc_file)
    if violations:
        print("Purity check FAILED:")
        for v in violations:
            print("  -", v)
        return 1

    print("[PASS] eval_rc_purity: 0 mpmath imports, 0 floating math.sin/cos/atan/log/exp/sqrt calls in eval_rc.py.")

    py_exe = _ensure_venv_python(_REPO)
    env = os.environ.copy()
    env["PYTHONPATH"] = str(_SRC) + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    res = subprocess.run([py_exe, "-m", "pytest", "-q"], cwd=str(_REPO), env=env, check=False)
    return int(res.returncode)


if __name__ == "__main__":
    raise SystemExit(main())
