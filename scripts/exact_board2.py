# Exact certifier for EinsteinArena board 2 (first-autocorrelation-inequality), MINIMIZE.
# Proven UPPER bound on C1* = inf_f (sup_t (f*f)(t))/(int f)^2.
# Certified path uses ONLY exact integer/rational arithmetic. No float on the certified value.
import json, os, sys, time
from fractions import Fraction
t0=time.time()

_DEF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "vectors", "ac_board2_leader.json")
d=json.load(open(sys.argv[1] if len(sys.argv) > 1 else _DEF))
vals=d["solution"]["values"]; n=len(vals)
claimed=d.get("score")

# 1. exact conversion (Fraction(float) = exact dyadic value of submitted vector)
fr=[Fraction(v) for v in vals]
assert all(x>=0 for x in fr), "verifier guard: all f_i >= 0"
Ssum_frac=sum(fr)
assert Ssum_frac>0, "verifier guard: integral > 0"

# 2. common denominator D=2^E (floats are dyadic), integer numerators a_i
maxden=1
for x in fr:
    if x.denominator>maxden: maxden=x.denominator
E=maxden.bit_length()-1
assert maxden==(1<<E), "denominator not a power of two"
D=maxden
a=[x.numerator*(D//x.denominator) for x in fr]
S=sum(a)                      # S = D * sum(f_i)
maxbits=max(x.bit_length() for x in a)

# 3. exact integer autoconvolution A[m]=sum_{i+j=m} a_i a_j via Kronecker substitution.
#    slot width b bits, byte-aligned; b > 2*maxbits + log2(n) guarantees no slot overflow.
bbits=2*maxbits + (n-1).bit_length() + 4
bbytes=(bbits+7)//8
b=bbytes*8
# pack a_i little-endian into b-bit slots
buf=bytearray(n*bbytes)
for i,ai in enumerate(a):
    if ai:
        buf[i*bbytes:(i+1)*bbytes]=ai.to_bytes(bbytes,'little')
X=int.from_bytes(buf,'little')
X2=X*X                        # exact big-int square  (Karatsuba)
L=2*n-1                       # number of conv coefficients m=0..2n-2
sq_bytes=X2.to_bytes(L*bbytes+bbytes,'little')  # generous length
mask=(1<<b)-1
M=0; argm=-1; checksum=0
for m in range(L):
    chunk=int.from_bytes(sq_bytes[m*bbytes:(m+1)*bbytes],'little')
    # chunk must fit in slot (no carry into next slot) if b chosen correctly
    checksum+=chunk
    if chunk>M:
        M=chunk; argm=m
assert checksum==S*S, "Kronecker slot overflow / packing error (sum A != S^2)"
assert M < (1<<b), "slot overflow at max coefficient"

# 4. exact functional. C = dx*max(conv)/(dx*sum f)^2 ; dx=1/(2n).
#    conv[m]=A[m]/D^2 , (sum f)=S/D  =>  C = (1/(2n)) * (A/D^2) / (S/D)^2 = 2n*A / (D^2 * ... )
#    careful: C = (1/(2n))^{-1}?  dx*max / (dx*sum)^2 = max/(dx*sum^2). dx=1/(2n) => = 2n*max/sum^2 (in conv/D units)
#    = 2n*(M/D^2)/((S/D)^2) = 2n*(M/D^2)/(S^2/D^2) = 2n*M/S^2.
C_exact=Fraction(2*n*M, S*S)

# 5. cross-check float(exact) vs verifier score
Cf=float(C_exact)

# 25-digit enclosure ROUNDED UP (safe direction for an UPPER bound)
num,den=C_exact.numerator, C_exact.denominator
scale=10**25
q,r=divmod(num*scale,den)
if r!=0: q+=1                 # ceil -> guaranteed >= true value
s=str(q).rjust(26,'0')
dec_up=s[:-25]+'.'+s[-25:]

print("n",n)
print("claimed_score",repr(claimed))
print("E_denom_exp",E,"maxbits_a",maxbits,"slot_b_bits",b)
print("S_bits",S.bit_length(),"M_bits",M.bit_length(),"argmax_m",argm)
print("C_exact_num_digits",len(str(C_exact.numerator)),"den_digits",len(str(C_exact.denominator)))
print("C_exact_num", C_exact.numerator)
print("C_exact_den", C_exact.denominator)
print("float_C_exact",repr(Cf))
print("absdiff_vs_claimed",abs(Cf-claimed))
print("C_exact_25dp_roundUP",dec_up)
print("below_leader?", C_exact < Fraction(claimed))
print("secs",round(time.time()-t0,1))
