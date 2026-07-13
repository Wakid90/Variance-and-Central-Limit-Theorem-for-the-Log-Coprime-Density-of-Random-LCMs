#!/usr/bin/env python3
"""Reproduce all numerical outputs used by the manuscript."""

from pathlib import Path

from make_variance_figure import main as make_variance_main
from simulate_clt import run_experiment


if __name__ == "__main__":
    # The variance script uses its own CLI defaults.  Calling it directly here
    # would parse this script's arguments, so reproduce the deterministic output
    # with an explicit lightweight wrapper instead.
    import subprocess
    import sys

    root = Path(__file__).resolve().parents[1]
    subprocess.run(
        [
            sys.executable,
            str(root / "code" / "make_variance_figure.py"),
            "--figure",
            str(root / "paper" / "figures" / "fig_variance.pdf"),
            "--output",
            str(root / "data" / "variance_table.csv"),
        ],
        check=True,
        cwd=root,
    )
    run_experiment(
        root / "paper" / "figures" / "fig_clt_distribution.png",
        root / "data" / "simulation_diagnostics.csv",
        seed=2026,
    )
