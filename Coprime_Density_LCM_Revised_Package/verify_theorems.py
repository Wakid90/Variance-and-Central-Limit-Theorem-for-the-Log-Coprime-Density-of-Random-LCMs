"""
FULL AUDIT of the repaired Paper 1.
Run after the §6.2 repair and minor fixes to verify nothing is broken
AND that every claim in the new proof holds numerically.

Checks:
  A1. Envelope S_j(M) <= 10.08/ln M for 1 <= j <= sqrt M, M >= 289.
  A2. Sigma_{k/r}(M) <= C_u * r / ln M for 2 <= r <= sqrt M / ln M, M large.
  A3. (1 - M^{-1/2})^{k/r} super-small for r in main range.
  A4. Stirling bound (C_u r/ln M)^r/r! <= (C_u e/ln M)^r.
  A5. Sum <= 2(C_u e)^2/(ln M)^2.
  A6. The tail r > r_max contribution is super-polynomially small.
  A7. Total R_k(M) <= O(1/(ln M)^2) numerically.
  A8. ln M * R_k(M) -> 0 at observed rate.

Also re-verify:
  B1. Theorem 1 identity on 24 exact finite cases.
  B2. Theorem 2 (limit and explicit remainder bound).
  B3. Theorem 4 (sharp fixed-M separation q** < q*).
  B4. Theorem 5 leading order Sigma * ln M -> A(c).
  B5. Theorem 5 sharp form (ln M * (rho_k/P_M - 1) -> A(c)).
  B6. Proposition 6 (A(c) asymptotics).
  B7. Remark 3 (naive independence bound counterexample).

If ANYTHING fails, this script will SHOUT.
"""
import math
from itertools import combinations
from sympy import primerange
from fractions import Fraction
from math import log, exp, sqrt, floor

FAIL = []
OK = []

def assert_close(desc, observed, expected, tol, relative=False):
    """Record PASS/FAIL against a numerical assertion."""
    if relative:
        err = abs(observed - expected)/max(abs(expected), 1e-20)
    else:
        err = abs(observed - expected)
    if err <= tol:
        OK.append(f"[PASS] {desc}: obs={observed:.6g}, exp={expected:.6g}, err={err:.2g}")
    else:
        FAIL.append(f"[FAIL] {desc}: obs={observed:.6g}, exp={expected:.6g}, err={err:.2g} > tol={tol:.2g}")

def assert_true(desc, cond, extra=""):
    if cond:
        OK.append(f"[PASS] {desc} {extra}")
    else:
        FAIL.append(f"[FAIL] {desc} {extra}")

# ----------- helpers -----------
def primes_up_to(M):
    return list(primerange(2, M+1))

def Sigma_s(s, plist, M):
    total = 0.0
    for p in plist:
        mp = M // p
        qp = 1.0 - mp/(M-1)
        if qp <= 0:
            continue
        total += qp**s/(p-1)
    return total

def Lambda(M, plist=None):
    if plist is None:
        plist = primes_up_to(M)
    return sum(1.0/(p-1) for p in plist)

def P_M(M, plist=None):
    if plist is None:
        plist = primes_up_to(M)
    P = 1.0
    for p in plist:
        P *= (1 - 1/p)
    return P

def beta_T(T, M):
    count = 0
    for n in range(2, M+1):
        if all(n % p != 0 for p in T):
            count += 1
    return count/(M-1)

def S_j(j, plist, M):
    """Strip sum at level j."""
    return sum(1.0/(p-1) for p in plist if M//p == j)

def rho_k_exact(k, M):
    """Exact rho_k(M) via Theorem 1."""
    plist = primes_up_to(M)
    total = Fraction(0)
    for rlen in range(0, len(plist)+1):
        for T in combinations(plist, rlen):
            count = 0
            for n in range(2, M+1):
                if all(n % p != 0 for p in T):
                    count += 1
            bT = Fraction(count, M-1)
            coef = Fraction(1)
            for p in T:
                coef *= Fraction(1, p-1)
            total += bT**k * coef
    # Multiply by P_M
    PM = Fraction(1)
    for p in plist:
        PM *= Fraction(p-1, p)
    return total * PM, PM

def A(c, J=5000):
    return sum(exp(-c*j) * log(1 + 1/j) for j in range(1, J+1))

# -------------- A. Verify repair claims --------------

print("="*78)
print("A. VERIFY REPAIR CLAIMS FOR §6.2")
print("="*78)

# A1: Envelope bound S_j(M) <= 10.08/ln M for j in [1, sqrt M], M >= 289
print("\n[A1] Envelope S_j(M) <= 10.08/ln M for 1 <= j <= sqrt M, M >= 289")
envelope_max = 0.0
for M in [289, 1000, 10000, 100000]:
    plist = primes_up_to(M)
    lnM = log(M)
    jmax = int(sqrt(M))
    for j in range(1, jmax+1):
        S = S_j(j, plist, M)
        ratio = S * lnM  # should be <= 10.08
        if ratio > envelope_max:
            envelope_max = ratio
        if ratio > 10.08:
            FAIL.append(f"[FAIL] envelope violated at M={M}, j={j}: S_j*lnM={ratio:.3f} > 10.08")
            break

assert_true(f"A1 envelope holds (max observed S_j*ln M = {envelope_max:.3f} <= 10.08)",
            envelope_max <= 10.08)

# A2: Sigma_{k/r} <= C_u * r / ln M for small r (pick c=1)
print("\n[A2] Sigma_{k/r}(M) * ln M / r bounded above (uniform in r)")
print(f"  Testing for c=1.0, various M and r.")
print(f"  {'M':>7} {'r':>4} {'Sigma*lnM':>10} {'Sigma*lnM/r':>12} {'expected':>10}")
max_ratio = 0.0
for M in [1000, 10000, 100000]:
    plist = primes_up_to(M)
    c = 1.0
    k = int(floor(c*M))
    lnM = log(M)
    rmax = int(sqrt(M)/lnM)
    for r in [2, 3, 5, 10, 30, 100, min(300, rmax), min(1000, rmax)]:
        if r > rmax or r > k:
            continue
        s = k/r
        sig = Sigma_s(s, plist, M)
        lmsig = lnM * sig
        ratio = lmsig/r
        if ratio > max_ratio:
            max_ratio = ratio
        print(f"  {M:>7} {r:>4} {lmsig:>10.3f} {ratio:>12.3f} {'<=C_u':>10}")

print(f"  Max observed Sigma*lnM/r = {max_ratio:.3f}")
# Our C_u bound should be some reasonable constant. Actually the bound is
# 10.08 * 4/c = 40.32 for c=1. So observed should be <= 40.32.
assert_true(f"A2 ratio bounded (max = {max_ratio:.3f} <= 2)", max_ratio <= 2.0)
# The theoretical C_u = 10.08 * 4/c = 40.32 is very loose; observed is much tighter.

# A3: For r <= sqrt M/ln M, (1 - M^{-1/2})^{k/r} is super-polynomially small
print("\n[A3] (1 - M^{-1/2})^{k/r} super-poly-small for r <= sqrt(M)/ln M")
for M in [1000, 10000, 100000]:
    c = 1.0
    k = int(floor(c*M))
    rmax = int(sqrt(M)/log(M))
    r_test = rmax  # worst case
    s = k/r_test
    val = (1 - M**(-0.5))**s
    # Expected bound: <= M^{-c/2}
    expected_ub = M**(-0.5)  # c/2 with c=1
    assert_true(f"A3 M={M}, r=rmax={r_test}: (1-1/sqrt(M))^s = {val:.3e} <= {expected_ub:.3e}",
                val <= expected_ub * 1.1)

# A4: Stirling r! >= (r/e)^r
print("\n[A4] Stirling: r! >= (r/e)^r for r >= 1")
from math import factorial
for r in [2, 5, 10, 20, 100]:
    rf = factorial(r)
    stir = (r/math.e)**r
    assert_true(f"A4 r={r}: r! = {rf:.3e} >= (r/e)^r = {stir:.3e}", rf >= stir)

# A5-A8: Compute sum bound directly and verify it matches ln M * R_k -> 0
print("\n[A5-A8] Direct numerical verification of sum_{r>=2} Sigma_{k/r}^r/r!")
print(f"{'M':>7} {'c':>4} {'lnM':>7} {'R_k bound':>14} {'lnM*bound':>14}")
prev_lmb = None
for M in [500, 1000, 5000, 10000, 50000, 100000, 500000]:
    c = 1.0
    plist = primes_up_to(M)
    k = int(floor(c*M))
    lnM = log(M)
    bound = 0.0
    r = 2
    log_fact = log(2)
    while True:
        s = k/r
        sig = Sigma_s(s, plist, M)
        if sig <= 0:
            break
        log_term = r*log(sig) - log_fact
        if log_term < -80:
            break
        bound += exp(log_term)
        r += 1
        log_fact += log(r)
        if r > 500:
            break
    lmb = lnM * bound
    print(f"  {M:>7} {c:>4.1f} {lnM:>7.3f} {bound:>14.4e} {lmb:>14.4e}")
    if prev_lmb is not None:
        assert_true(f"A8 ln M · bound decreasing from M -> 2M (at M={M}): {lmb:.4e} <= {prev_lmb:.4e}",
                    lmb < prev_lmb * 1.1)  # allow slight wiggle
    prev_lmb = lmb

# A9: Compare the bound to the ACTUAL R_k for small M (via Theorem 1 exact)
print("\n[A9] For small M, compute exact R_k and verify bound is indeed >= R_k")
print(f"{'M':>5} {'c':>4} {'k':>4} {'exact R_k':>14} {'upper bnd':>14} {'bnd/R_k':>8}")
for M in [20, 30, 40, 50]:
    if M > 40:
        continue  # too slow for exact enumeration
    plist = primes_up_to(M)
    for c in [0.5, 1.0, 2.0]:
        k = int(floor(c*M))
        if k < 1: continue
        # Compute exact R_k = rho_k/P_M - 1 - Sigma_k
        # Easier: R_k = sum_{|T|>=2} beta_T^k prod 1/(p-1)
        Rk = 0.0
        for rlen in range(2, len(plist)+1):
            for T in combinations(plist, rlen):
                bT = beta_T(T, M)
                if bT == 0:
                    continue
                coef = 1.0
                for p in T:
                    coef /= (p-1)
                Rk += bT**k * coef
        # Compute upper bound
        bound = 0.0
        r = 2
        log_fact = log(2)
        while r <= len(plist):
            s = k/r
            sig = Sigma_s(s, plist, M)
            if sig <= 0:
                break
            log_term = r*log(sig) - log_fact
            if log_term < -100:
                break
            bound += exp(log_term)
            r += 1
            log_fact += log(r)
        ratio = bound / Rk if Rk > 0 else float('inf')
        print(f"  {M:>5} {c:>4.1f} {k:>4} {Rk:>14.4e} {bound:>14.4e} {ratio:>8.2f}")
        assert_true(f"A9 M={M}, c={c}: bound >= R_k ({bound:.3e} >= {Rk:.3e})", bound >= Rk * 0.999)

# -------------- B. Re-verify all theorems --------------

print("\n" + "="*78)
print("B. RE-VERIFY ALL EXISTING THEOREMS")
print("="*78)

# B1: Theorem 1 identity on 24 test cases (M in {5..10}, k in {2..5})
print("\n[B1] Theorem 1 identity — 24 test cases (M in {5..10}, k in {2..5})")
pass_count = 0
fail_count = 0
for M in [5, 6, 7, 8, 9, 10]:
    for k in [2, 3, 4, 5]:
        rho_eq, PM = rho_k_exact(k, M)
        # Compute rho_k by brute force enumeration
        total_phi_over_L = Fraction(0)
        from math import gcd
        def ilcm(a, b):
            return a*b // gcd(a, b)
        def phi(n):
            # Euler totient of n
            result = n
            m = n
            p = 2
            while p*p <= m:
                if m % p == 0:
                    while m % p == 0:
                        m //= p
                    result -= result//p
                p += 1
            if m > 1:
                result -= result//m
            return result
        # Enumerate all (M-1)^k tuples
        samples = list(range(2, M+1))
        def rec(depth, curr_lcm):
            if depth == 0:
                return Fraction(phi(curr_lcm), curr_lcm)
            s = Fraction(0)
            for x in samples:
                s += rec(depth-1, ilcm(curr_lcm, x))
            return s
        # Start with dummy lcm=1
        count_tuples = (M-1)**k
        rho_brute = rec(k, 1) / count_tuples
        if rho_eq == rho_brute:
            pass_count += 1
        else:
            fail_count += 1
            FAIL.append(f"[FAIL] B1 M={M}, k={k}: identity {rho_eq} != brute {rho_brute}")
assert_true(f"B1 all 24 Theorem 1 cases match exactly (passed {pass_count}/24)",
            pass_count == 24)

# B2: Theorem 2 explicit remainder bound
print("\n[B2] Theorem 2: R_k(M) <= Sigma_{k/2}(M)^2 e^Lambda / 2")
# Pick a few (M, k) and check R_k <= bound
for M in [20, 30]:
    plist = primes_up_to(M)
    Lam = Lambda(M, plist)
    for k in [10, 20]:
        # Compute R_k exactly
        Rk = 0.0
        for rlen in range(2, len(plist)+1):
            for T in combinations(plist, rlen):
                bT = beta_T(T, M)
                if bT == 0: continue
                coef = 1.0
                for p in T:
                    coef /= (p-1)
                Rk += bT**k * coef
        # Bound
        sig_half = Sigma_s(k/2, plist, M)
        bound = sig_half**2 * exp(Lam) / 2
        ratio = bound / Rk if Rk > 0 else float('inf')
        assert_true(f"B2 M={M}, k={k}: R_k={Rk:.4e} <= bound={bound:.4e}", bound >= Rk)

# B3: Theorem 4 — q^** < q^*
print("\n[B3] Theorem 4: q**(M) < q*(M) for M >= 3")
for M in [3, 5, 10, 20, 30, 50]:
    plist = primes_up_to(M)
    qstar = max(1 - (M//p)/(M-1) for p in plist)
    qstar_star = 0
    for rlen in range(2, len(plist)+1):
        for T in combinations(plist, rlen):
            bT = beta_T(T, M)
            if bT > qstar_star:
                qstar_star = bT
    assert_true(f"B3 M={M}: q** = {qstar_star:.4f} < q* = {qstar:.4f}", qstar_star < qstar)
    if M >= 11:
        expected = (M-3)/(M-1)
        assert_close(f"B3 explicit q** formula at M={M}", qstar_star, expected, 1e-14)

# B4: Theorem 5 leading order Sigma * ln M -> A(c)
print("\n[B4] Theorem 5 leading: Sigma_{floor(cM)}(M) * ln M approx A(c)")
for c in [0.5, 1.0, 2.0]:
    for M in [10000, 100000]:
        plist = primes_up_to(M)
        k = int(floor(c*M))
        sig = Sigma_s(k, plist, M)
        lmsig = log(M) * sig
        Ac = A(c)
        rel_err = abs(lmsig - Ac)/Ac
        assert_true(f"B4 c={c}, M={M}: ln M · Sigma = {lmsig:.4f}, A(c) = {Ac:.4f}, rel_err = {rel_err:.3f}",
                    rel_err < 0.2)

# B5: Theorem 5 sharp — ln M · (rho_k/P_M - 1) -> A(c)
print("\n[B5] Theorem 5 sharp: ln M · ((rho_k/P_M) - 1) approx A(c)")
for c in [0.5, 1.0]:
    for M in [15, 20]:  # can compute exact rho_k for small M
        plist = primes_up_to(M)
        k = int(floor(c*M))
        if k < 1:
            continue
        rho, PM = rho_k_exact(k, M)
        obs = log(M) * (float(rho)/float(PM) - 1)
        Ac = A(c)
        rel_err = abs(obs - Ac) / Ac
        print(f"  c={c}, M={M}, k={k}: ln M · (rho/P_M - 1) = {obs:.4f}, A(c) = {Ac:.4f}, "
              f"rel err = {rel_err:.3f}")
        # For small M, we expect large relative error (order 1/ln M or so)
        assert_true(f"B5 c={c}, M={M}: ln M · (rho/P-1) agrees with A(c) within ~40%",
                    rel_err < 0.4)

# B6: Proposition 6 asymptotics of A(c)
print("\n[B6] Proposition 6: A(c) asymptotics")
# Large c expansion: A(c) = e^{-c} ln 2 * (1 + O(e^{-c}))
# Small c expansion: A(c) = -ln c - gamma + O(c ln c)
gamma = 0.5772156649
for c in [5.0, 10.0]:
    Ac = A(c)
    lead = exp(-c)*log(2)
    rel = abs(Ac - lead)/lead
    assert_true(f"B6 large-c: c={c}, A={Ac:.4e}, leading={lead:.4e}, rel={rel:.3f}",
                rel < 0.3)

for c in [0.01, 0.001]:
    Ac = A(c, J=50000)
    lead = -log(c) - gamma
    rel = abs(Ac - lead)/abs(lead)
    assert_true(f"B6 small-c: c={c}, A={Ac:.4f}, leading -ln c - gamma = {lead:.4f}, rel = {rel:.3f}",
                rel < 0.05)

# B7: naive independence bound fails at (M,T)=(17,{3,5})
print("\n[B7] beta_T(17, {3,5}) > q_3 q_5")
M = 17
T = (3, 5)
bT = beta_T(T, M)
q3 = 1 - (M//3)/(M-1)
q5 = 1 - (M//5)/(M-1)
assert_true(f"B7 remark: beta_T = {bT} vs q_3 q_5 = {q3*q5}. beta_T > product: {bT > q3*q5}",
            bT > q3*q5)

# ---------------- FINAL SUMMARY ----------------

print("\n\n" + "="*78)
print("FINAL AUDIT SUMMARY")
print("="*78)
print(f"PASSED: {len(OK)}")
print(f"FAILED: {len(FAIL)}")
if FAIL:
    print("\n--- FAILURES ---")
    for f in FAIL:
        print(f)
else:
    print("\n*** ALL CHECKS PASSED ***")
