#!/usr/bin/env python3
"""Exact finite-M variance formula for random-LCM coprime density.

The combinatorial formula is exact.  The default numerical evaluation uses
NumPy long-double arithmetic for speed.  The optional ``--mpmath`` mode
re-evaluates one parameter pair with arbitrary precision for validation.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import mpmath as mp
import numpy as np


@dataclass(frozen=True)
class VarianceResult:
    M: int
    c: float
    k: int
    V1: float
    V2: float
    variance: float
    rescaled_variance: float
    B_c: float


def primes_up_to(n: int) -> list[int]:
    if n < 2:
        return []
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, math.isqrt(n) + 1):
        if sieve[p]:
            start = p * p
            sieve[start : n + 1 : p] = b"\x00" * (((n - start) // p) + 1)
    return [i for i, flag in enumerate(sieve) if flag]


def variance_components(M: int, c: float) -> VarianceResult:
    """Vectorized long-double evaluation of the exact finite formula."""
    if M < 2:
        raise ValueError("M must be at least 2")
    if c <= 0:
        raise ValueError("c must be positive")
    k = math.floor(c * M)
    if k < 1:
        raise ValueError("floor(c*M) must be at least 1")

    dtype = np.longdouble
    primes = np.asarray(primes_up_to(M), dtype=np.int64)
    levels = M // primes
    denom = dtype(M - 1)
    q = (denom - levels.astype(dtype)) / denom
    weights = np.log1p(-dtype(1) / primes.astype(dtype))
    qk = np.power(q, k)
    V1 = np.sum(weights * weights * qk * (dtype(1) - qk), dtype=dtype)

    p_col = primes[:, None]
    level_col = levels[:, None]
    common = M // (p_col * p_col.T)
    beta_num = (M - 1) - level_col - level_col.T + common
    beta = beta_num.astype(dtype) / denom
    covariance = np.power(beta, k) - np.power(q[:, None] * q[None, :], k)
    weight_matrix = weights[:, None] * weights[None, :]
    np.fill_diagonal(covariance, dtype(0))
    V2 = np.sum(weight_matrix * covariance, dtype=dtype)

    variance = V1 + V2
    rescaled = dtype(M) * np.log(dtype(M)) * variance
    c_ld = dtype(str(c))
    B_c = dtype(1) / np.expm1(c_ld) - dtype(1) / np.expm1(dtype(2) * c_ld)
    return VarianceResult(
        M=M,
        c=c,
        k=k,
        V1=float(V1),
        V2=float(V2),
        variance=float(variance),
        rescaled_variance=float(rescaled),
        B_c=float(B_c),
    )


def variance_components_mpmath(M: int, c: float, dps: int = 80) -> dict[str, mp.mpf]:
    """Slow arbitrary-precision validation for a single (M,c) pair."""
    mp.mp.dps = dps
    k = math.floor(c * M)
    denom = M - 1
    primes = primes_up_to(M)
    levels = [M // p for p in primes]
    q_values = [mp.mpf(denom - j) / denom for j in levels]
    weights = [mp.log1p(-mp.mpf(1) / p) for p in primes]
    V1 = mp.fsum(
        a * a * mp.power(q, k) * (1 - mp.power(q, k))
        for a, q in zip(weights, q_values)
    )
    terms = []
    for i, p in enumerate(primes):
        for j in range(i + 1, len(primes)):
            qprime = primes[j]
            beta = mp.mpf(
                denom - levels[i] - levels[j] + M // (p * qprime)
            ) / denom
            cov = mp.power(beta, k) - mp.power(q_values[i] * q_values[j], k)
            terms.append(2 * weights[i] * weights[j] * cov)
    V2 = mp.fsum(terms)
    variance = V1 + V2
    c_mp = mp.mpf(str(c))
    return {
        "V1": V1,
        "V2": V2,
        "variance": variance,
        "rescaled_variance": M * mp.log(M) * variance,
        "B(c)": 1 / mp.expm1(c_mp) - 1 / mp.expm1(2 * c_mp),
    }


def write_table(
    output_csv: Path,
    M_values: Iterable[int],
    c_values: Iterable[float],
) -> list[VarianceResult]:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    results: list[VarianceResult] = []
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["M", "c", "k", "V1", "V2", "variance", "M_log_M_variance", "B(c)"]
        )
        for M in M_values:
            for c in c_values:
                result = variance_components(M, c)
                results.append(result)
                writer.writerow(
                    [
                        result.M,
                        f"{result.c:g}",
                        result.k,
                        f"{result.V1:.17g}",
                        f"{result.V2:.17g}",
                        f"{result.variance:.17g}",
                        f"{result.rescaled_variance:.17g}",
                        f"{result.B_c:.17g}",
                    ]
                )
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--M", type=int)
    parser.add_argument("--c", type=float)
    parser.add_argument("--mpmath", action="store_true", help="arbitrary-precision single case")
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if args.M is not None or args.c is not None:
        if args.M is None or args.c is None:
            parser.error("--M and --c must be supplied together")
        if args.mpmath:
            result = variance_components_mpmath(args.M, args.c, args.dps)
            for key, value in result.items():
                print(f"{key}={mp.nstr(value, 30)}")
        else:
            result = variance_components(args.M, args.c)
            for key in ("V1", "V2", "variance", "rescaled_variance", "B_c"):
                print(f"{key}={getattr(result, key):.17g}")
        return

    output = args.output or Path("data/variance_table.csv")
    results = write_table(
        output,
        M_values=(100, 200, 400, 800, 1600, 3200),
        c_values=(0.5, 1.0, 2.0, 3.0),
    )
    print(f"Wrote {len(results)} rows to {output}")


if __name__ == "__main__":
    main()
