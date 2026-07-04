"""Exact-rational certification of two PUBLISHED signed step vectors for the
third autocorrelation inequality (C3, signed max semantics):

  (a) Matolcsi--Vinuesa closing-remark example, n=150, printed coefficients
      (arXiv:0907.1379 source, appendix; 8-significant-digit decimals).
  (b) AlphaEvolve public n=400 vector (google-deepmind/alphaevolve_results,
      mathematical_results.ipynb height_sequence_3; candidate JSON copy).

Semantics (arena board 4 / third-autocorrelation-inequality):
  dx = 0.5/n on [-1/4, 1/4];  C3 = |max( conv(f,f)*dx )| / (sum(f)*dx)^2.
For a signed step function f on the uniform grid, g = f*f is piecewise linear
with breakpoints on the half-grid and knot values dx*conv[k]; hence
  sup_t g(t) = dx * max_k conv[k]  (exactly), and
  C3(f) = 2n * max_k conv[k] / (sum f)^2 (dx and any common scale cancel).
Coefficients are parsed as exact decimal rationals (the printed values),
convolution is exact big-integer; NO floating point on the certified path.
Also reports max_k |conv| to expose the max|g| vs max g distinction.
"""
import re, json, sys
from fractions import Fraction
from decimal import Decimal, getcontext
getcontext().prec = 50


def exact_c3(dec_strings, label):
    n = len(dec_strings)
    fr = [Fraction(s) for s in dec_strings]  # exact decimal parse
    # common denominator -> integers
    from math import lcm
    D = 1
    for x in fr:
        D = lcm(D, x.denominator)
    F = [int(x * D) for x in fr]
    assert all(Fraction(a, D) == x for a, x in zip(F, fr))
    S = sum(F)
    # exact integer autoconvolution, direct O(n^2)
    L = 2 * n - 1
    conv = [0] * L
    for i in range(n):
        Fi = F[i]
        for j in range(n):
            conv[i + j] += Fi * F[j]
    maxc = max(conv)
    minc = min(conv)
    amax = conv.index(maxc)
    maxabs = max(maxc, -minc)
    assert maxc > 0 and S != 0
    C3 = Fraction(2 * n * maxc, S * S)          # signed-max semantics (board)
    C3abs = Fraction(2 * n * maxabs, S * S)     # abs-inside-max cousin
    dC3 = Decimal(C3.numerator) / Decimal(C3.denominator)
    dAbs = Decimal(C3abs.numerator) / Decimal(C3abs.denominator)
    print(f"=== {label} ===")
    print(f"n = {n}   sum-of-coeffs (float) = {float(Fraction(S, D)):.10f}")
    print(f"exact C3 (signed max)  num = {C3.numerator}")
    print(f"                       den = {C3.denominator}")
    print(f"exact C3 decimal (40dp) = {dC3}")
    print(f"float(C3_exact)         = {float(C3)!r}")
    print(f"argmax k = {amax} of 0..{L-1}   min conv (float scale) = {float(Fraction(minc, D*D)) * (0.5/n):.6f}")
    print(f"max_t |f*f| exact decimal (40dp) = {dAbs}")
    print(f"float(max|f*f| ratio)   = {float(C3abs)!r}")
    print(f"abs()-inert (max|g| attained at positive max)? {maxabs == maxc}")
    print()
    return C3, C3abs


# (a) Matolcsi--Vinuesa n=150 printed vector
import os as _os
_VEC = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "vectors")
_mv_path = sys.argv[1] if len(sys.argv) > 1 else _os.path.join(_VEC, "mv150_printed.txt")
_ae_path = sys.argv[2] if len(sys.argv) > 2 else _os.path.join(_VEC, "AE-2025-third-400.json")
mv = open(_mv_path).read().split()
assert len(mv) == 150
c3_mv, abs_mv = exact_c3(mv, "MATOLCSI-VINUESA printed n=150 signed vector (arXiv:0907.1379 appendix)")

# (b) AlphaEvolve public n=400 vector: parse raw JSON text for decimal literals
raw = open(_ae_path).read()
d = json.loads(raw)
# Print only whitelisted, publication-safe metadata keys (never echo an
# arbitrary 'source' path that a locally-wrapped input JSON might carry).
_SAFE_META = ("note", "upstream_url", "upstream_sha256", "retrieved")
meta = {k: d[k] for k in _SAFE_META if k in d}
print("AE candidate metadata:", meta)
m = re.search(r'"values"\s*:\s*\[(.*?)\]', raw, re.S)
ae = re.findall(r'-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?', m.group(1))
assert len(ae) == len(d["values"]) == 400
# cross-check decimal parse matches the JSON floats to float precision
for s, v in zip(ae, d["values"]):
    assert abs(float(s) - v) == 0.0
c3_ae, abs_ae = exact_c3(ae, "ALPHAEVOLVE public n=400 vector (height_sequence_3)")

# margins vs our exact leader certificate Q3
Q3 = Fraction(11753128449293701953238517385067272445617294540800000,
              8092744874989952471246071559466128309374865340943729)
for name, c in (("MV", c3_mv), ("AE", c3_ae)):
    diff = c - Q3
    print(f"{name} cert minus Q3 = {float(diff)!r}  (positive => Q3 is lower/better)")
    print(f"  40dp: {Decimal(diff.numerator)/Decimal(diff.denominator)}")
