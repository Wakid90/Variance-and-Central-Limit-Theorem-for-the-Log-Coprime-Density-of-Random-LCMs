"""Reproduce every numerical output in the paper, from scratch.

    python code/reproduce_all.py

Rewrites:
    paper/figures/fig_variance.pdf          (Figure 1)
    paper/figures/fig_clt_distribution.png  (Figure 2)
    data/variance_table.csv                 (Table 1)
    data/simulation_diagnostics.csv         (Figure 2 diagnostics)

and runs an independent 80-decimal arbitrary-precision validation of the
extended-precision variance values.

Everything is deterministic: the variance computation is exact, and the Monte
Carlo panels are seeded (numpy.random.default_rng(2026), freshly initialised
per panel).
"""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import make_variance_figure  # noqa: E402
import simulate_clt  # noqa: E402
from exact_variance import rescaled_variance, rescaled_variance_mp  # noqa: E402

# (M, c) pairs checked at 80 decimal digits. Kept small: the arbitrary-precision
# path is quadratic in pi(M) with full-precision powers.
SPOT_CHECKS = [(100, c) for c in (0.5, 1.0, 2.0, 3.0)] + [
    (400, c) for c in (0.5, 1.0, 2.0, 3.0)
]
TOL = 1e-12


def validate() -> None:
    print("== arbitrary-precision validation (80 decimal digits) ==")
    worst = 0.0
    for M, c in SPOT_CHECKS:
        fast, _, _ = rescaled_variance(M, c)
        slow = float(rescaled_variance_mp(M, c, dps=80))
        delta = abs(fast - slow)
        worst = max(worst, delta)
        status = "OK" if delta < TOL else "FAIL"
        print(f"  M = {M:4d}, c = {c:>3}:  |delta| = {delta:.2e}  {status}")
        if delta >= TOL:
            raise AssertionError(
                f"extended precision disagrees with 80-digit value at (M, c) = ({M}, {c}): "
                f"{fast!r} vs {slow!r}"
            )
    print(f"  worst discrepancy: {worst:.2e}  (tolerance {TOL:.0e})")


def main() -> None:
    print("== Table 1 and Figure 1 (exact, deterministic) ==")
    make_variance_figure.main()

    print()
    print("== Figure 2 (Monte Carlo, seed 2026) ==")
    simulate_clt.main()

    print()
    validate()

    print()
    print("All outputs regenerated.")


if __name__ == "__main__":
    main()
