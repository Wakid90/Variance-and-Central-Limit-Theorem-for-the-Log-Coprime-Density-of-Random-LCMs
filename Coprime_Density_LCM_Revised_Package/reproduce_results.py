#!/usr/bin/env python3
"""Reproduce the numerical figures and tables in the preprint.

Outputs
-------
fig1_convergence.pdf
fig2_joint_scaling.pdf
figure1_points.csv
table1_monte_carlo.csv
table2_fixed_M.csv
table3_joint_scaling.csv
numerical_report.txt

The Monte Carlo calculations are deterministic. Each (M,k) cell uses seed
2026 + 1000*M + k and N=2000 replicates unless stated otherwise.
"""
from __future__ import annotations

import csv
import math
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

BASE_SEED = 2026
N_MC = 2000
A_SERIES_TERM_TOL = 1e-17
ROOT = Path(__file__).resolve().parent

TABLE1_CELLS = [
    (100, 50), (100, 200), (300, 100), (300, 500),
    (1000, 500), (1000, 2000), (1000, 5000), (3000, 5000),
]
TABLE3_CS = [0.05, 0.10, 0.15, 0.25, 0.50, 0.75, 1.00, 1.25,
             1.50, 2.00, 3.00, 5.00, 7.50, 10.00, 15.00]
TABLE3_MS = [10_000, 100_000, 1_000_000, 10_000_000]


def sieve_primes(n: int) -> np.ndarray:
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    limit = int(n**0.5)
    for p in range(2, limit + 1):
        if sieve[p]:
            sieve[p*p:n+1:p] = False
    return np.flatnonzero(sieve)


def prime_masks(M: int):
    """Return masks[n], factors, primes, and P_M for 2 <= n <= M."""
    spf = list(range(M + 1))
    for p in range(2, int(M**0.5) + 1):
        if spf[p] == p:
            for n in range(p*p, M+1, p):
                if spf[n] == n:
                    spf[n] = p
    primes = [p for p in range(2, M+1) if spf[p] == p]
    index = {p: i for i, p in enumerate(primes)}
    masks = [0] * (M + 1)
    for n in range(2, M + 1):
        x = n
        mask = 0
        while x > 1:
            p = spf[x]
            mask |= 1 << index[p]
            while x % p == 0:
                x //= p
        masks[n] = mask
    factors = [1.0 - 1.0/p for p in primes]
    return masks, factors, primes, math.prod(factors)


def phi_ratio_from_mask(mask: int, factors: list[float]) -> float:
    value = 1.0
    while mask:
        low = mask & -mask
        idx = low.bit_length() - 1
        value *= factors[idx]
        mask ^= low
    return value


def sigma_k(M: int, k: float, primes) -> float:
    p = np.asarray(primes, dtype=np.int64)
    q = 1.0 - (M // p) / (M - 1)
    return float(np.sum(np.power(q, k) / (p - 1)))


def A(c: float) -> float:
    total = 0.0
    j = 1
    while True:
        term = math.exp(-c*j) * math.log1p(1.0/j)
        total += term
        if term < A_SERIES_TERM_TOL:
            return total
        j += 1


def monte_carlo_cell(M: int, k: int, N: int = N_MC):
    masks, factors, primes, PM = prime_masks(M)
    seed = BASE_SEED + 1000*M + k
    rng = np.random.default_rng(seed)
    vals = np.empty(N, dtype=float)
    for i in range(N):
        xs = rng.integers(2, M+1, size=k)
        union = 0
        for x in np.unique(xs):
            union |= masks[int(x)]
        vals[i] = phi_ratio_from_mask(union, factors)
    rho = float(vals.mean())
    se = float(vals.std(ddof=1) / math.sqrt(N))
    pred = PM * sigma_k(M, k, primes)
    gap = rho - PM
    return {
        "M": M, "k": k, "N": N, "seed": seed,
        "rho_mc": rho, "standard_error": se, "P_M": PM,
        "gap_observed": gap, "gap_predicted": pred,
        "ratio": gap / pred,
    }


def write_csv(path: Path, rows: list[dict]):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def make_table1(report):
    rows = []
    for M, k in TABLE1_CELLS:
        row = monte_carlo_cell(M, k)
        rows.append(row)
        report.append(
            f"Table 1 M={M}, k={k}: rho={row['rho_mc']:.8f}, "
            f"SE={row['standard_error']:.8f}, P_M={row['P_M']:.8f}, "
            f"gap={row['gap_observed']:.8f}, pred={row['gap_predicted']:.8f}, "
            f"ratio={row['ratio']:.5f}"
        )
    write_csv(ROOT / "table1_monte_carlo.csv", rows)
    return rows


def make_figure1(report):
    Ms = [100, 300, 1000]
    k_values = [5, 8, 12, 20, 35, 50, 100, 200, 500, 1000, 2000]
    rows = []
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    markers = {100: "o", 300: "s", 1000: "^"}
    for M in Ms:
        means, ses = [], []
        PM = None
        for k in k_values:
            row = monte_carlo_cell(M, k)
            rows.append(row)
            means.append(row["rho_mc"])
            ses.append(row["standard_error"])
            PM = row["P_M"]
        ax.errorbar(k_values, means, yerr=ses, marker=markers[M],
                    markersize=5, capsize=2, linewidth=1.1, label=f"M={M}")
        ax.axhline(PM, linestyle="--", linewidth=0.9, alpha=0.65)
    ax.set_xscale("log")
    ax.set_xlabel("k (number of samples)")
    ax.set_ylabel(r"$\rho_k(M)=\mathbb{E}[\varphi(L_k)/L_k]$")
    ax.set_title("Convergence of coprime density to Mertens product")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(ROOT / "fig1_convergence.pdf", bbox_inches="tight")
    plt.close(fig)
    write_csv(ROOT / "figure1_points.csv", rows)
    report.append("Figure 1 regenerated from deterministic Monte Carlo points.")


def beta_T(M: int, T: tuple[int, ...]) -> float:
    return sum(all(n % p for p in T) for n in range(2, M+1)) / (M-1)


def make_table2(report):
    rows = []
    for M in [20, 30]:
        primes = sieve_primes(M).tolist()
        betas = []
        for r in range(2, len(primes)+1):
            for T in combinations(primes, r):
                coef = math.prod(1.0/(p-1) for p in T)
                betas.append((beta_T(M, T), coef))
        for k in [20, 40, 80, 160]:
            R = sum((b**k)*coef for b, coef in betas)
            Sigma = sigma_k(M, k, primes)
            rate = ((M-3)/(M-2))**k
            rows.append({"M": M, "k": k, "R_over_Sigma": R/Sigma,
                         "qstarstar_over_qstar_power_k": rate})
    write_csv(ROOT / "table2_fixed_M.csv", rows)
    report.append("Table 2 regenerated by exact finite subset enumeration.")


def make_table3_and_figure2(report):
    all_primes = sieve_primes(max(TABLE3_MS))
    rows = []
    for c in TABLE3_CS:
        row = {"c": c, "A_c": A(c)}
        for M in TABLE3_MS:
            p = all_primes[all_primes <= M]
            row[f"M_{M}"] = math.log(M) * sigma_k(M, math.floor(c*M), p)
        rows.append(row)
    write_csv(ROOT / "table3_joint_scaling.csv", rows)

    plot_Ms = [1_000, 10_000, 100_000, 1_000_000]
    cs = np.concatenate([np.array([0.10,0.15,0.20,0.25,0.30,0.40,0.50,0.60,0.75]),
                         np.arange(1.0, 5.01, 0.25)])
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    markers = {1_000:"o", 10_000:"s", 100_000:"^", 1_000_000:"d"}
    for M in plot_Ms:
        p = all_primes[all_primes <= M]
        ys = [math.log(M)*sigma_k(M, math.floor(float(c)*M), p) for c in cs]
        ax.plot(cs, ys, linestyle="", marker=markers[M], markersize=5,
                alpha=0.85, label=f"M={M}")
    c_fine = np.linspace(0.05, 5.0, 500)
    ax.plot(c_fine, [A(float(c)) for c in c_fine], linewidth=1.4,
            label=r"$A(c)=\sum_{j\geq1}e^{-cj}\ln(1+1/j)$")
    ax.set_xlabel(r"$c=k/M$")
    ax.set_ylabel(r"$\ln M\,\Sigma_{\lfloor cM\rfloor}(M)$")
    ax.set_title(r"Joint scaling: convergence to $A(c)$")
    ax.set_xlim(0, 5)
    ax.set_ylim(bottom=-0.05)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(ROOT / "fig2_joint_scaling.pdf", bbox_inches="tight")
    plt.close(fig)
    report.append("Table 3 and Figure 2 regenerated from deterministic direct prime sums in double precision.")
    report.append(f"A(c) series truncated adaptively when the next term is below {A_SERIES_TERM_TOL:.0e}.")


def main():
    report = [
        "Numerical reproduction report",
        "Base Monte Carlo seed: 2026",
        "Per-cell seed: 2026 + 1000*M + k",
        f"Monte Carlo replicates per cell: {N_MC}",
        "",
    ]
    make_table1(report)
    make_figure1(report)
    make_table2(report)
    make_table3_and_figure2(report)
    (ROOT / "numerical_report.txt").write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report))


if __name__ == "__main__":
    main()
