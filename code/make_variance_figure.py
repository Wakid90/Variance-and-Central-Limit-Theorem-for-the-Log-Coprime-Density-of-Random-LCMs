"""Regenerate Figure 1 and the Table 1 data from the exact finite formula.

Writes:
    paper/figures/fig_variance.pdf
    data/variance_table.csv

Fully deterministic: the rescaled variance M ln M * s_k^2 is evaluated from the
exact finite formula s_k^2 = V1 + V2 by direct summation over the primes up to M
and over all ordered prime pairs. No Monte Carlo is involved.

Fonts are embedded as TrueType (pdf.fonttype = 42) rather than as Type 3 bitmap
subsets, which several journals reject.
"""

from __future__ import annotations

import csv
import pathlib

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "pdf.fonttype": 42,   # TrueType, not Type 3
        "ps.fonttype": 42,
        "font.size": 10,
    }
)
import matplotlib.pyplot as plt  # noqa: E402

from exact_variance import B, rescaled_variance  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
FIG = ROOT / "paper" / "figures" / "fig_variance.pdf"
CSV = ROOT / "data" / "variance_table.csv"

CS = [0.5, 1.0, 2.0, 3.0]
MS = [100, 200, 400, 800, 1600, 3200]

STYLE = {
    0.5: ("tab:blue", "o"),
    1.0: ("tab:orange", "s"),
    2.0: ("tab:green", "^"),
    3.0: ("tab:red", "D"),
}


def main() -> None:
    FIG.parent.mkdir(parents=True, exist_ok=True)
    CSV.parent.mkdir(parents=True, exist_ok=True)

    values = {}
    rows = []
    for M in MS:
        for c in CS:
            total, v1, v2 = rescaled_variance(M, c)
            values[(M, c)] = total
            rows.append(
                {
                    "M": M,
                    "c": c,
                    "k": int(c * M),
                    "MlnM_s2": f"{total:.10f}",
                    "MlnM_V1": f"{v1:.10f}",
                    "MlnM_V2": f"{v2:.10e}",
                    "B_c": f"{B(c):.10f}",
                }
            )
            assert v2 < 0, f"V2 should be negative at (M, c) = ({M}, {c})"

    with CSV.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {CSV.relative_to(ROOT)}  ({len(rows)} rows)")

    fig, ax = plt.subplots(figsize=(8.4, 5.1))
    for c in CS:
        col, mk = STYLE[c]
        ax.plot(
            MS,
            [values[(M, c)] for M in MS],
            marker=mk,
            ms=5,
            lw=1.4,
            color=col,
            label=f"c = {c:g}",
        )
        ax.axhline(B(c), color=col, ls="--", lw=0.9, alpha=0.65)
        ax.annotate(
            f"B({c:g}) = {B(c):.3f}",
            xy=(1.005, B(c)),
            xycoords=("axes fraction", "data"),
            va="center",
            ha="left",
            fontsize=8.5,
            color=col,
        )

    ax.set_xscale("log")
    ax.set_xlabel("M")
    ax.set_ylabel(r"$M\,\ln M\;\mathrm{Var}[\ln(\varphi(L_k)/L_k)]$")
    ax.set_title(
        r"Convergence of the rescaled variance to $B(c)$, at $k=\lfloor cM\rfloor$",
        fontsize=11,
    )
    ax.set_ylim(0.0, 1.05)
    ax.legend(loc="center left", fontsize=9, framealpha=0.9)
    fig.subplots_adjust(right=0.86)
    fig.savefig(FIG)
    print(f"wrote {FIG.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
