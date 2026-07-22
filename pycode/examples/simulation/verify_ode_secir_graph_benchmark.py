#!/usr/bin/env python3
"""Compare final states from both benchmark runners."""
import argparse

import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument("cpp")
parser.add_argument("python")
parser.add_argument("nodes", type=int)
parser.add_argument("age_groups", type=int)
args = parser.parse_args()

cpp = np.fromfile(args.cpp, dtype=np.float64)
python = np.fromfile(args.python, dtype=np.float64)
expected = args.nodes * args.age_groups * 10
if cpp.size != expected or python.size != expected:
    raise SystemExit(
        f"wrong state size: expected {expected}, got {cpp.size}/{python.size}")
finite = np.isfinite(cpp).all() and np.isfinite(python).all()
difference = np.abs(cpp - python)
max_abs = difference.max(initial=0.0)
scale = np.maximum(1.0, np.maximum(np.abs(cpp), np.abs(python)))
max_rel = np.max(difference / scale, initial=0.0)
bitwise = np.array_equal(cpp.view(np.uint64), python.view(np.uint64))
close = finite and np.allclose(cpp, python, rtol=1e-12, atol=1e-9)
print(f"{args.nodes},{args.age_groups},{expected},{int(bitwise)},"
      f"{max_abs:.17g},{max_rel:.17g},{int(close)}")
if not close:
    raise SystemExit("C++ and Python results differ")
