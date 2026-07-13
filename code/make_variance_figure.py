#!/usr/bin/env python3
"""Generate the rescaled-variance figure and its machine-readable CSV."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from exact_variance import variance_components, write_table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--figure", type=Path, default=Path("paper/figures/fig_variance.pdf")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("data/variance_table.csv")
    )
    args = parser.parse_args()

    M_values = (100, 200, 400, 800, 1600, 3200)
    c_values = (0.5, 1.0, 2.0, 3.0)
    results = write_table(args.output, M_values, c_values)

    by_c = {c: [] for c in c_values}
    limits = {}
    for result in results:
        by_c[result.c].append((result.M, float(result.rescaled_variance)))
        limits[result.c] = float(result.B_c)

    args.figure.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.4, 5.1), constrained_layout=True)
    markers = ("o", "s", "^", "D")
    for marker, c in zip(markers, c_values):
        pairs = sorted(by_c[c])
        x = [item[0] for item in pairs]
        y = [item[1] for item in pairs]
        ax.plot(x, y, marker=marker, label=f"c = {c:g}")
        ax.axhline(limits[c], linestyle="--", linewidth=1.0, alpha=0.65)
        ax.text(x[-1] * 1.03, limits[c], f"B({c:g}) = {limits[c]:.3f}", va="center", fontsize=8)

    ax.set_xscale("log")
    ax.set_xlabel("M")
    ax.set_ylabel(r"$M\ln M\,\mathrm{Var}[\ln(\varphi(L_k)/L_k)]$")
    ax.set_title(r"Convergence of the rescaled variance to $B(c)$, at $k=\lfloor cM\rfloor$")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.25)
    fig.savefig(args.figure)
    plt.close(fig)
    print(f"Wrote {args.figure}")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
