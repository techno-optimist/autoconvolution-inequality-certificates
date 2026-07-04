import re, ast, os
from fractions import Fraction
from decimal import Decimal, getcontext
getcontext().prec = 80
_VEC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "vectors")

def exact_ratio_from_int_heights(vint):
    """vint: list of nonneg integers = heights * S (S a common scale).
    Returns exact Fraction of the SCALE-INVARIANT continuum ratio
        ||f*f||_2^2 / (||f*f||_1 * ||f*f||_inf)
    computed with the EXACT closed forms:
      f*f is piecewise-linear with knot values L[j] = discrete autocorrelation.
      SquareL2 = (1/3) sum_j (L[j]^2 + L[j]L[j+1] + L[j+1]^2)   [exact integral of (pw-linear)^2]
      LInf     = max_j L[j]                                     [max of pw-linear is at a knot]
      L1       = (1/2) sum_j (L[j] + L[j+1])                    [exact trapezoid integral]
    Scale S cancels in the ratio (degree-4 numerator / degree-4 denominator)."""
    n = len(vint)
    # autocorrelation L[k] = sum_i v[i]*v[k-i], k=0..2n-2  (all nonneg)
    assert all(x >= 0 for x in vint), "packing trick needs nonneg heights"
    SLOT = 256  # bits per slot; max L < n * max(v)^2, keep well under 2^256
    packed = 0
    for i, x in enumerate(vint):
        packed |= x << (SLOT * i)
    sq = packed * packed
    mask = (1 << SLOT) - 1
    L = [ (sq >> (SLOT*k)) & mask for k in range(2*n - 1) ]
    # sanity: no slot overflow (next slot boundary clean) -> check reconstruct sum equals total^2 later
    # SquareL2 numerator*3
    S2x3 = 0
    for j in range(len(L)-1):
        a, b = L[j], L[j+1]
        S2x3 += a*a + a*b + b*b
    # include the last knot pair with implicit 0 beyond support: L extends to 0, contributes a^2 for last
    aL = L[-1]
    S2x3 += aL*aL  # pair (L[last], 0): a^2 + a*0 + 0
    # also leading pair (0, L[0]) contributes L[0]^2
    a0 = L[0]
    S2x3 += a0*a0
    LInf = max(L)
    # L1*2 = sum over all unit knot-gaps of (L[j]+L[j+1]); with zero endpoints
    # total integral = sum_j (L[j]) since trapezoid with 0 ends telescopes to sum of interior... 
    # Do it explicitly with zero padding:
    Lpad = [0] + L + [0]
    L1x2 = 0
    for j in range(len(Lpad)-1):
        L1x2 += Lpad[j] + Lpad[j+1]
    # ratio = (S2x3/3) / (LInf * (L1x2/2)) = (2*S2x3) / (3*LInf*L1x2)
    num = 2*S2x3
    den = 3*LInf*L1x2
    return Fraction(num, den), L, LInf

def dec(fr, p=40):
    return (Decimal(fr.numerator)/Decimal(fr.denominator)).quantize(Decimal(10)**-p)

# ---------- Boyer-Li: 575 integers ----------
bl = [int(x) for x in re.findall(r'-?\d+', open(os.path.join(_VEC, 'coeffBL.txt')).read())]
assert len(bl)==575
r_bl, L_bl, inf_bl = exact_ratio_from_int_heights(bl)
print("=== BOYER-LI (575 integers, exact) ===")
print("our exact ratio =", r_bl.numerator, "/", r_bl.denominator)
print("our exact decimal =", dec(r_bl, 40))
published = Fraction(96086089410408840289339769, 106577013451431545242354944)
print("published fraction=", published.numerator, "/", published.denominator)
print("published decimal =", dec(published, 40))
print("EXACT MATCH to published fraction?", r_bl == published)
print(">= 0.901564 ?", r_bl >= Fraction(901564, 10**6))
print(">= 0.901562 (JJ's citation) ?", r_bl >= Fraction(901562, 10**6))

# ---------- Jaech-Joseph high-res: 2540 floats ----------
vf = ast.literal_eval(open(os.path.join(_VEC, 'jj_function.py')).read())
print("\n=== JAECH-JOSEPH github function.py (%d floats) ===" % len(vf))
print("min=%.6g max=%.6g any_negative=%s" % (min(vf), max(vf), any(x<0 for x in vf)))
# exact: each float -> exact dyadic rational; common scale = 2^maxexp
ratios = [Fraction(x).limit_denominator(1) if False else Fraction(*x.as_integer_ratio()) for x in vf]
maxden = max(fr.denominator for fr in ratios)
# all denominators are powers of two; scale to integers
vint = [ (fr.numerator * (maxden // fr.denominator)) for fr in ratios ]
assert all(fr.numerator*(maxden//fr.denominator)==fr*maxden for fr in ratios)
r_jj, L_jj, inf_jj = exact_ratio_from_int_heights(vint)
print("our EXACT continuum ratio =", dec(r_jj, 40))
print("as fraction num/den bits:", r_jj.numerator.bit_length(), "/", r_jj.denominator.bit_length())
print(">= 0.94136 (their claim) ?", r_jj >= Fraction(94136,10**5))
print(">= 0.926529 ?", r_jj >= Fraction(926529,10**6))
print(">= our C2* 0.9629011010961756912885013 ?", r_jj >= Fraction(Decimal('0.9629011010961756912885013')))

# Also reproduce their FLOAT method-ish sanity: naive float ratio using same exact-knot closed form but in float
import math
n=len(vf)
Lf=[0.0]*(2*n-1)
for i,a in enumerate(vf):
    for k,b in enumerate(vf):
        Lf[i+k]+=a*b
S2=0.0
Lp=[0.0]+Lf+[0.0]
for j in range(len(Lp)-1):
    S2 += (Lp[j]**2+Lp[j]*Lp[j+1]+Lp[j+1]**2)/3.0
LInf=max(Lf)
L1=sum((Lp[j]+Lp[j+1])/2.0 for j in range(len(Lp)-1))
print("float closed-form ratio =", S2/(LInf*L1))
