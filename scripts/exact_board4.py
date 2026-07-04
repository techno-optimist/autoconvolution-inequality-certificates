#!/usr/bin/env python3
"""
EXACT certifier for the third autocorrelation inequality (signed variant),
EinsteinArena problem `third-autocorrelation-inequality`.
MINIMIZE. f MAY BE SIGNED, supported [-1/4,1/4], n=100000, leader OrganonAgent 1.4523043331831582.

A certified evaluation of the leader gives a PROVEN UPPER BOUND on the inf-constant C3*.
Exact big-integer / Fraction arithmetic ONLY on the certified path (each float64 is an
exact dyadic rational). Signed Kronecker: F = P - Q, conv(F,F)=conv(P,P)-2conv(P,Q)+conv(Q,Q).

Functional (pinned verifier): C3(f) = |max_p (conv(f,f)*dx)[p]| / (sum(f)*dx)^2.
Lemma: scaled_conv[p] = g(t_{p+1}) exactly (g = continuum autoconvolution, piecewise
linear, max attained at an interior vertex, positive), so score = exact continuum ratio
up to fp rounding, ZERO discretization loss. dx cancels: C3 = 2n*max_p c_p / (sum F)^2.
"""
import json, os, random, sys, time
from fractions import Fraction
import numpy as np

t0 = time.time()
LEADER = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "vectors", "ac_board4_leader.json")
d = json.load(open(LEADER))
vals = d["solution"]["values"]
n = len(vals)
reported = d.get("score")
print(f"[load] n={n} reported={reported!r}", flush=True)

# ---------------------------------------------------------------------------
# VERIFIER-FIRST: reproduce the published verifier's EXACT functional.
# Pinned verifier: score = |max_p (conv(f,f)*dx)[p]| / (sum(f)*dx)^2 with
# conv = np.convolve(f,f,'full'), dx = 0.5/n. We reproduce with scipy fftconvolve
# (O(n log n), fp-equivalent to np.convolve; the box is at load ~90 with production
# kissat/LP jobs, so the O(n^2) np.convolve/ArenaClient path is starved). This
# reproduction has been shown to match the board score to all 17 digits.
# ---------------------------------------------------------------------------
from scipy.signal import fftconvolve
f64 = np.array(vals, dtype=np.float64)
dx = 0.5 / n
scaled = fftconvolve(f64, f64) * dx          # length 2n-1, == np.convolve(f,f,'full')*dx
arena_score = float(abs(np.max(scaled)) / ((np.sum(f64) * dx) ** 2))
print(f"[verify] fftconvolve reproduction of pinned verifier = {arena_score!r}", flush=True)
print(f"[verify] matches board leader {reported!r}: {arena_score == reported}", flush=True)
print(f"[verify] fp argmax p = {int(np.argmax(scaled))}  signed max = {float(np.max(scaled))}"
      f"  signed min = {float(np.min(scaled))}", flush=True)

# ---------------------------------------------------------------------------
# EXACT RATIONAL PATH (no floats on the certified path)
# ---------------------------------------------------------------------------
fr = [Fraction(v) for v in vals]           # exact: float64 IS a dyadic rational
maxexp = 0
for x in fr:
    den = x.denominator
    assert (den & (den - 1)) == 0, "denominator not a power of two"
    maxexp = max(maxexp, den.bit_length() - 1)
D = 1 << maxexp
F = [int(x * D) for x in fr]               # exact signed integers, F_k = f_k * D
Ssum = sum(F)
Mabs = max(abs(v) for v in F)
print(f"[rat] D=2^{maxexp}  |F|max bitlen={Mabs.bit_length()}  Ssum bitlen={Ssum.bit_length()}", flush=True)
assert Ssum != 0, "integral is zero -- inadmissible"

# signed split for nonneg Kronecker packing
P = [v if v > 0 else 0 for v in F]
Q = [-v if v < 0 else 0 for v in F]

# slot width: 8W > bitlen(n * Mabs^2) + 1
bound = n * Mabs * Mabs
W = (bound.bit_length() + 2 + 7) // 8
L = 2 * n - 1
print(f"[kron] slot bytes W={W}  L={L}", flush=True)

def pack(A):
    buf = bytearray(W * len(A))
    for k, a in enumerate(A):
        if a:
            nb = (a.bit_length() + 7) // 8
            buf[k * W:k * W + nb] = a.to_bytes(nb, "little")
    return int.from_bytes(buf, "little")

def unpack(x):
    nb = (x.bit_length() + 7) // 8
    b = x.to_bytes(nb, "little").ljust(W * L, b"\x00")
    return [int.from_bytes(b[i * W:(i + 1) * W], "little") for i in range(L)]

pP, pQ = pack(P), pack(Q)
cPP = unpack(pP * pP); print("[kron] cPP done", flush=True)
cQQ = unpack(pQ * pQ); print("[kron] cQQ done", flush=True)
cPQ = unpack(pP * pQ); print("[kron] cPQ done", flush=True)
cint = [cPP[i] - 2 * cPQ[i] + cQQ[i] for i in range(L)]   # exact conv(F,F)[p]
assert sum(cint) == Ssum * Ssum, "checksum sum_m c_m = S^2 failed (slot overflow?)"

# spot-check Kronecker against direct O(n) sums at random positions
for _ in range(8):
    pp = random.randint(0, L - 1)
    s = 0
    for i in range(max(0, pp - (n - 1)), min(pp, n - 1) + 1):
        s += F[i] * F[pp - i]
    assert s == cint[pp], (pp, s, cint[pp])
print("[kron] spot-checks vs direct O(n): OK", flush=True)

maxc = max(cint)
amax = cint.index(maxc)
assert maxc > 0, "max autoconvolution not positive -- lemma abs() inertness fails"

# dx cancels: C3 = 2n * maxc / Ssum^2  (exact Fraction)
C3 = Fraction(2 * n * maxc, Ssum * Ssum)

print(f"[exact] integer argmax p = {amax}", flush=True)
print(f"[exact] C3 numerator digits = {len(str(C3.numerator))}", flush=True)
print(f"[exact] C3 denominator digits = {len(str(C3.denominator))}", flush=True)
print(f"[exact] float(C3) = {float(C3)!r}", flush=True)
print(f"[exact] reported  = {reported!r}", flush=True)
print(f"[exact] rel diff exact-vs-reported = {abs(float(C3) - reported) / reported:.3e}", flush=True)

# 25-significant-digit enclosure, ROUNDED UP (safe for an UPPER bound):
from decimal import Decimal, getcontext, ROUND_CEILING, ROUND_FLOOR
getcontext().prec = 60
dnum = Decimal(C3.numerator); dden = Decimal(C3.denominator)
val = dnum / dden
up25 = val.quantize(Decimal("1." + "0" * 24), rounding=ROUND_CEILING)
dn25 = val.quantize(Decimal("1." + "0" * 24), rounding=ROUND_FLOOR)
print(f"[exact] 25-digit UP  (safe upper bound) = {up25}", flush=True)
print(f"[exact] 25-digit DOWN                    = {dn25}", flush=True)

out = f"""third-autocorrelation-inequality (signed variant)
MINIMIZE. f MAY BE SIGNED. leader OrganonAgent, n={n}, reported={reported!r}
Direction: construction -> PROVEN UPPER BOUND on inf-constant C3*.

board_score_reproduced = {arena_score!r}   (matches board leader: {arena_score == reported})

EXACT C3(leader) = num/den (exact rational, dx cancels, D^2 cancels):
num={C3.numerator}
den={C3.denominator}
float(C3_exact) = {float(C3)!r}
integer argmax p = {amax}
maxc>0 (abs() inert) = {maxc > 0}
Ssum bitlen = {Ssum.bit_length()}  D=2^{maxexp}  W={W} bytes

25-sig-digit enclosure:
  UP  (report as upper bound): {up25}
  DOWN                       : {dn25}

PROVEN BOUND:  C3* <= {C3.numerator}/{C3.denominator}  =  {up25}  (rounded UP)

rel diff exact-vs-reported = {abs(float(C3) - reported) / reported:.3e}
NOTE: exact 1.45230433318315... is a HAIR BELOW reported float (fp rounding of the
SAME leader vector). This is an EXACT EVALUATION of the existing leader, NOT an
improvement; UPPER bound only.
runtime_seconds = {time.time() - t0:.1f}
"""
open("exact_board4_out.txt", "w").write(out)
print("[done] wrote exact_board4_out.txt (cwd)", flush=True)
print(f"[done] runtime {time.time() - t0:.1f}s", flush=True)
