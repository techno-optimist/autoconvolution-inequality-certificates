#!/usr/bin/env python3
"""Enclosure script (cited in the proofs of Propositions 5 and 7).

Computes the directed decimal enclosures of every re-certified value and every
margin quoted in the note, by EXACT Fraction floor/ceil arithmetic -- no
floating point on any printed bound. Output byte-matches certs/recert_enclosures.txt.

The Jaech-Joseph and Boyer-Li exact ratios are recomputed here from the frozen
input vectors (vectors/); the theorem/proposition fractions Q2, Q3 and the
Matolcsi-Vinuesa / AlphaEvolve certified ratios are the exact rationals proved
in the note (reproduced by certify.py / cert_published_c3.py).
"""
import os
from fractions import Fraction
from decimal import Decimal, getcontext, ROUND_FLOOR, ROUND_CEILING

getcontext().prec = 90
_VEC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "vectors")


def exact_c2_ratio(int_heights):
    """Exact ||g||_2^2 / (||g||_1 ||g||_inf) for a nonneg step vector, via the
    piecewise-linear closed forms (Lemma 1(iii)-(iv)). Integer O(n^2) conv."""
    v = list(int_heights)
    n = len(v)
    L = [0] * (2 * n - 1)
    for i, a in enumerate(v):
        if a:
            for j, b in enumerate(v):
                L[i + j] += a * b
    Lpad = [0] + L + [0]
    S2x3 = sum(Lpad[j] ** 2 + Lpad[j] * Lpad[j + 1] + Lpad[j + 1] ** 2
               for j in range(len(Lpad) - 1))
    LInf = max(L)
    L1x2 = sum(Lpad[j] + Lpad[j + 1] for j in range(len(Lpad) - 1))
    return Fraction(2 * S2x3, 3 * LInf * L1x2)


def sig_enclosure(x, sig=25):
    """Directed [floor, ceil] enclosure of Fraction x at `sig` significant digits."""
    val = Decimal(x.numerator) / Decimal(x.denominator)
    q = Decimal(1).scaleb(val.adjusted() - (sig - 1))
    return val.quantize(q, rounding=ROUND_FLOOR), val.quantize(q, rounding=ROUND_CEILING)


def floor_dp(x, dp):
    """Lower bound of Fraction x to `dp` decimal places (for '> ...' margins)."""
    val = Decimal(x.numerator) / Decimal(x.denominator)
    return val.quantize(Decimal(1).scaleb(-dp), rounding=ROUND_FLOOR)


# ---- recompute JJ, BL exact ratios from frozen vectors ----
import ast
jj = ast.literal_eval(open(os.path.join(_VEC, "jj_function.py")).read())
jj_fr = [Fraction(*x.as_integer_ratio()) for x in jj]
D = max(f.denominator for f in jj_fr)
jj_int = [int(f * D) for f in jj_fr]
JJ = exact_c2_ratio(jj_int)

import re
bl = [int(x) for x in re.findall(r"-?\d+", open(os.path.join(_VEC, "coeffBL.txt")).read())]
BL = exact_c2_ratio(bl)

# ---- certified fractions proved in the note ----
Q2 = Fraction(140651861665566489683881393353250795846281833,
              146070932420211259869783468438333325818535926)
Q3 = Fraction(11753128449293701953238517385067272445617294540800000,
              8092744874989952471246071559466128309374865340943729)
MV = Fraction(174972681773148828, 119999999644045921)          # max_t (v*v)/(int v)^2
AE = Fraction(116451423625336345095200, 79999999996821408075961)

lines = []
lo, hi = sig_enclosure(JJ);  lines.append(f"JJ exact C2 ratio in [ {lo} , {hi} ]   (25 sig digits, floor/ceil)")
lines.append(f"BL exact C2 ratio    = {sig_enclosure(BL)[0]}   (== published 96086089410408840289339769/106577013451431545242354944)")
lo, hi = sig_enclosure(MV);  lines.append(f"MV cert (max_t) in [ {lo} , {hi} ]")
lo, hi = sig_enclosure(AE);  lines.append(f"AE cert (max_t) in [ {lo} , {hi} ]")
lines.append("")
lines.append(f"Q2 - JJ  > {floor_dp(Q2 - JJ, 16)}   (exact begins {floor_dp(Q2 - JJ, 20)})")
lines.append(f"Q2 - BL  > {floor_dp(Q2 - BL, 16)}")
lines.append(f"MV - Q3  > {floor_dp(MV - Q3, 16)}")
lines.append(f"AE - Q3  > {floor_dp(AE - Q3, 16)}")

out = "\n".join(lines) + "\n"
outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certs", "recert_enclosures.txt")
open(outpath, "w").write(out)
print(out, end="")
