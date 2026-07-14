"""Exact finite-M evaluation of the variance of the log coprime density.

For X_1, ..., X_k i.i.d. uniform on {2, ..., M} and L_k = lcm(X_1, ..., X_k),

    W_k(M) = ln(phi(L_k)/L_k) - E ln(phi(L_k)/L_k),
    s_k^2  = Var W_k(M) = V1 + V2,

where, writing m_p = floor(M/p), pi_p = m_p/(M-1), q_p = 1 - pi_p, and
beta_{p,q} = P(p does not divide X, q does not divide X) for a single draw X,

    V1 = sum_p        ln(1 - 1/p)^2 * q_p^k * (1 - q_p^k),
    V2 = sum_{p != q} ln(1 - 1/p) ln(1 - 1/q) * (beta_{p,q}^k - q_p^k q_q^k).

V2 is the ORDERED off-diagonal sum, so each unordered prime pair contributes
twice. The full prime-pair formula is used, including pairs with p*q <= M.

Theorem 1 of the paper states M ln M * s_k^2 -> B(c) as M -> infinity along
k = floor(cM), with B(c) = 1/(e^c - 1) - 1/(e^{2c} - 1).

Two independent evaluations are provided:

  * rescaled_variance(...)  -- extended precision (numpy.longdouble, 80-bit on
    x86-64). This is what populates Table 1 and Figure 1.

  * rescaled_variance_mp(...) -- arbitrary precision (mpmath, configurable
    number of decimal digits). V2 is a difference of nearly equal quantities,
    so this guards against catastrophic cancellation in the fast path.

Both are fully deterministic; no Monte Carlo is involved.

Command line
------------
    python code/exact_variance.py --M 100 --c 1
    python code/exact_variance.py --M 100 --c 1 --mpmath --dps 80
"""

from __future__ import annotations

import argparse
import math

import numpy as np


def primes_upto(n: int) -> np.ndarray:
    """Primes <= n, by sieve of Eratosthenes."""
    if n < 2:
        return np.array([], dtype=np.int64)
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    return np.nonzero(sieve)[0].astype(np.int64)


def B(c: float) -> float:
    """Limiting constant B(c) = 1/(e^c - 1) - 1/(e^{2c} - 1) of Theorem 1."""
    return 1.0 / math.expm1(c) - 1.0 / math.expm1(2.0 * c)


def variance_parts(M: int, k: int, dtype=np.longdouble):
    """Return (V1, V2) exactly, in extended-precision floating point."""
    if M < 3:
        raise ValueError("M must be at least 3")
    P = primes_upto(M)
    m = M // P                                        # levels m_p
    pi = m.astype(dtype) / dtype(M - 1)
    q = dtype(1) - pi
    a = np.log1p(-dtype(1) / P.astype(dtype))         # ln(1 - 1/p) < 0
    qk = np.exp(k * np.log(q))

    V1 = np.sum(a * a * qk * (dtype(1) - qk))

    # beta_{p,q} = (M - 1 - m_p - m_q + floor(M/(pq))) / (M - 1)
    PQ = np.outer(P, P)
    mpq = (M // PQ).astype(dtype)
    beta = (dtype(M - 1) - m[:, None] - m[None, :] + mpq) / dtype(M - 1)
    np.fill_diagonal(beta, dtype(1))                  # diagonal excluded below
    betak = np.exp(k * np.log(beta))
    cross = betak - np.outer(qk, qk)
    np.fill_diagonal(cross, dtype(0))

    V2 = np.sum(np.outer(a, a) * cross)
    return V1, V2


def rescaled_variance(M: int, c: float, dtype=np.longdouble):
    """Return (M ln M * s_k^2, M ln M * V1, M ln M * V2) at k = floor(cM)."""
    k = int(math.floor(c * M))
    V1, V2 = variance_parts(M, k, dtype=dtype)
    scale = dtype(M) * np.log(dtype(M))
    return float(scale * (V1 + V2)), float(scale * V1), float(scale * V2)


def rescaled_variance_mp(M: int, c: float, dps: int = 80):
    """Independent arbitrary-precision recomputation of M ln M * s_k^2.

    Slow (quadratic in pi(M) with full-precision powers); intended for spot
    checks at moderate M, not for populating the whole table.
    """
    from mpmath import mp, mpf, log, power

    mp.dps = dps
    k = int(math.floor(c * M))
    P = [int(p) for p in primes_upto(M)]
    one = mpf(1)
    Mm1 = mpf(M - 1)

    m = {p: M // p for p in P}
    q = {p: one - mpf(m[p]) / Mm1 for p in P}
    a = {p: log(one - one / mpf(p)) for p in P}
    qk = {p: power(q[p], k) for p in P}

    V1 = sum(a[p] ** 2 * qk[p] * (one - qk[p]) for p in P)

    V2 = mpf(0)
    for i, p in enumerate(P):
        for r in P[i + 1 :]:
            beta = (Mm1 - m[p] - m[r] + (M // (p * r))) / Mm1
            # factor 2: the ordered sum counts each unordered pair twice
            V2 += 2 * a[p] * a[r] * (power(beta, k) - qk[p] * qk[r])

    return mpf(M) * log(mpf(M)) * (V1 + V2)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--M", type=int, required=True, help="support bound M")
    ap.add_argument("--c", type=float, required=True, help="ratio c, with k = floor(cM)")
    ap.add_argument("--mpmath", action="store_true", help="use arbitrary precision")
    ap.add_argument("--dps", type=int, default=80, help="decimal digits for --mpmath")
    args = ap.parse_args()

    if args.mpmath:
        val = rescaled_variance_mp(args.M, args.c, dps=args.dps)
        from mpmath import nstr

        print(nstr(val, 30))
    else:
        total, v1, v2 = rescaled_variance(args.M, args.c)
        k = int(math.floor(args.c * args.M))
        print(f"M = {args.M}, c = {args.c}, k = {k}")
        print(f"  M ln M * V1    = {v1:.12f}")
        print(f"  M ln M * V2    = {v2:.12e}")
        print(f"  M ln M * s_k^2 = {total:.12f}")
        print(f"  B(c)           = {B(args.c):.12f}")


if __name__ == "__main__":
    main()
