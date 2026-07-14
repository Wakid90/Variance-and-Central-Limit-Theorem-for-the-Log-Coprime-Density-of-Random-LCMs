"""Regenerate Figure 2: Monte Carlo diagnostic for the central limit theorem.

Writes:
    paper/figures/fig_clt_distribution.png
    data/simulation_diagnostics.csv

For each panel we draw N independent replicates of (X_1, ..., X_k), i.i.d.
uniform on {2, ..., M}, form L_k = lcm(X_1, ..., X_k) implicitly by collecting
the set of primes dividing at least one X_i, and evaluate the coprime density

    phi(L_k)/L_k = prod_{p | L_k} (1 - 1/p).

The standardised statistic (phi(L_k)/L_k - rho_hat) / sigma_hat is histogrammed
against the standard normal density, where rho_hat and sigma_hat are the sample
mean and sample standard deviation (ddof = 1).

Reproducibility: each panel uses numpy.random.default_rng(SEED), freshly
initialised, so the figure is bit-for-bit reproducible.

Reported skewness and excess kurtosis (Fisher convention) are descriptive moment
statistics, not hypothesis tests.
"""

from __future__ import annotations

import csv
import pathlib

import numpy as np
import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 9})
import matplotlib.pyplot as plt  # noqa: E402

from exact_variance import primes_upto  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
FIG = ROOT / "paper" / "figures" / "fig_clt_distribution.png"
CSV = ROOT / "data" / "simulation_diagnostics.csv"

SEED = 2026
PANELS = [(100, 3000), (300, 3000), (1000, 1000)]  # (M, N); k = M throughout


def smallest_prime_factor_table(M: int) -> np.ndarray:
    """spf[n] = smallest prime factor of n, for n <= M."""
    spf = np.zeros(M + 1, dtype=np.int64)
    for i in range(2, M + 1):
        if spf[i] == 0:
            spf[i::i] = np.where(spf[i::i] == 0, i, spf[i::i])
    return spf


def prime_index_masks(M: int):
    """For each n in {2..M}, a boolean row over primes <= M marking p | n."""
    P = primes_upto(M)
    idx = {int(p): j for j, p in enumerate(P)}
    masks = np.zeros((M + 1, len(P)), dtype=bool)
    spf = smallest_prime_factor_table(M)
    for n in range(2, M + 1):
        x = n
        while x > 1:
            p = int(spf[x])
            masks[n, idx[p]] = True
            while x % p == 0:
                x //= p
    return P, masks


def coprime_densities(M: int, k: int, N: int, rng) -> np.ndarray:
    """N independent samples of phi(L_k)/L_k."""
    P, masks = prime_index_masks(M)
    log_factor = np.log1p(-1.0 / P.astype(np.float64))  # ln(1 - 1/p)

    out = np.empty(N, dtype=np.float64)
    for r in range(N):
        draws = rng.integers(2, M + 1, size=k)
        hit = masks[draws].any(axis=0)      # primes dividing L_k
        out[r] = np.exp(log_factor[hit].sum())
    return out


def moments(x: np.ndarray):
    """Sample mean, sd (ddof=1), skewness, excess kurtosis (Fisher)."""
    n = x.size
    mean = x.mean()
    sd = x.std(ddof=1)
    z = (x - mean) / x.std(ddof=0)
    skew = float((z**3).mean())
    exkurt = float((z**4).mean() - 3.0)
    return float(mean), float(sd), skew, exkurt, n


def main() -> None:
    FIG.parent.mkdir(parents=True, exist_ok=True)
    CSV.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, len(PANELS), figsize=(11.0, 3.2))
    grid = np.linspace(-4, 4, 400)
    normal = np.exp(-0.5 * grid**2) / np.sqrt(2 * np.pi)

    rows = []
    for ax, (M, N) in zip(axes, PANELS):
        k = M
        rng = np.random.default_rng(SEED)   # freshly initialised per panel
        samples = coprime_densities(M, k, N, rng)
        mean, sd, skew, exkurt, n = moments(samples)
        z = (samples - mean) / sd

        ax.hist(z, bins=35, range=(-4, 4), density=True, alpha=0.75, color="tab:blue")
        ax.plot(grid, normal, color="tab:orange", lw=1.8, label=r"$\mathcal{N}(0,1)$")
        ax.set_title(
            f"M = {M}, k = {k}\nskew = {skew:.2f}, excess kurt = {exkurt:.2f}",
            fontsize=9,
        )
        ax.set_xlabel(r"$(\varphi(L_k)/L_k - \hat\rho_k)/\hat\sigma_k$")
        ax.legend(fontsize=8, loc="upper left")
        if ax is axes[0]:
            ax.set_ylabel("density")

        rows.append(
            {
                "M": M,
                "k": k,
                "N": n,
                "seed": SEED,
                "sample_mean": f"{mean:.10f}",
                "sample_sd": f"{sd:.10e}",
                "skewness": f"{skew:.4f}",
                "excess_kurtosis": f"{exkurt:.4f}",
            }
        )
        print(f"M = {M:4d}  N = {n:4d}  skew = {skew:+.4f}  excess kurt = {exkurt:+.4f}")

    fig.tight_layout()
    fig.savefig(FIG, dpi=200)
    print(f"wrote {FIG.relative_to(ROOT)}")

    with CSV.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {CSV.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
