import json, os, sys, time
from fractions import Fraction
from decimal import Decimal, getcontext

getcontext().prec = 80
t_start = time.time()

_DEF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "vectors", "ac_board3_leader.json")
d = json.load(open(sys.argv[1] if len(sys.argv) > 1 else _DEF))
vals = d["solution"]["values"]
n = len(vals)
reported = d["score"]
print(f"n = {n}  reported = {reported!r}", flush=True)

# STEP 1: exact dyadic scaling to integers
# each float is a dyadic rational num/den, den a power of 2. D = max den.
ratios = [v.as_integer_ratio() for v in vals]   # (num, den)
D = max(den for _, den in ratios)
assert (D & (D-1)) == 0, "D not power of 2"
A = [num * (D // den) for num, den in ratios]
assert all(isinstance(a, int) for a in A)
maxA = max(A)
minA = min(A)
print(f"D = {D}  maxA bitlen = {maxA.bit_length()}  minA = {minA}", flush=True)
assert minA >= 0, "nonneg board but negative scaled int"

# STEP 2: Kronecker squaring
maxbits = maxA.bit_length()
b_needed = 2*maxbits + n.bit_length() + 2
b = ((b_needed + 7)//8)*8
if b < 96: b = 96
nb = b // 8
print(f"b_needed = {b_needed}  chosen b = {b}  nb = {nb} bytes", flush=True)

buf = bytearray(nb * n)
for i, a in enumerate(A):
    buf[i*nb:(i+1)*nb] = a.to_bytes(nb, 'little')
P = int.from_bytes(buf, 'little')
print(f"P bits = {P.bit_length()}  squaring...", flush=True)
t0 = time.time()
P2 = P * P
print(f"squaring done in {time.time()-t0:.1f}s", flush=True)
P2bytes = P2.to_bytes(nb*(2*n), 'little')  # enough room
c = [int.from_bytes(P2bytes[k*nb:(k+1)*nb], 'little') for k in range(2*n-1)]
assert len(c) == 2*n - 1, f"len(c)={len(c)}"

# checksum: S1 == (sum A)^2
sumA = sum(A)
S1 = sum(c)
assert S1 == sumA*sumA, f"S1 checksum FAIL {S1} != {sumA*sumA}"
print("S1 == (sum A)^2 checksum: PASS", flush=True)

# STEP 3: exact norms (h cancels)
Linf = max(c)
# S2 = sum_{m}(y_m^2 + y_m y_{m+1} + y_{m+1}^2) with y=[0]+c+[0]
# = 2*sum(c_k^2) + sum_k c_k c_{k+1}   (endpoints 0)
sum_sq = sum(ck*ck for ck in c)
sum_adj = sum(c[k]*c[k+1] for k in range(len(c)-1))
S2 = 2*sum_sq + sum_adj

C = Fraction(S2, 3*S1*Linf)
print(f"S1 = {S1}", flush=True)
print(f"Linf = {Linf}", flush=True)
print(f"S2 = {S2}", flush=True)
print(f"exact_num = {C.numerator}", flush=True)
print(f"exact_den = {C.denominator}", flush=True)
dec = Decimal(C.numerator) / Decimal(C.denominator)
print(f"exact_decimal = {dec}", flush=True)
print(f"float(exact) = {float(C)!r}", flush=True)
print(f"reported     = {reported!r}", flush=True)
print(f"abs diff float(exact) vs reported = {abs(float(C)-reported):.3e}", flush=True)

out = {
  "n": n,
  "exact_num": str(C.numerator), "exact_den": str(C.denominator),
  "exact_decimal": str(dec)[:62],
  "float_exact": float(C), "reported": reported,
  "board_verifier_rerun": reported,
  "S1": str(S1), "Linf": str(Linf), "S2": str(S2),
  "kron_bits": b, "elapsed_s": round(time.time()-t_start,1),
}
json.dump(out, open("board3_exact_result_v2.json","w"), indent=1)
print("saved board3_exact_result_v2.json (cwd)", flush=True)
print(f"TOTAL {time.time()-t_start:.1f}s", flush=True)
