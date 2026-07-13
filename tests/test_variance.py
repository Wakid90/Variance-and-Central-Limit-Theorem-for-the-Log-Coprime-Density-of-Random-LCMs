from __future__ import annotations

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from exact_variance import variance_components  # noqa: E402


def test_table_value_M100_c1() -> None:
    result = variance_components(100, 1.0)
    assert math.isclose(result.rescaled_variance, 0.3876825327460782, rel_tol=0, abs_tol=2e-15)
    assert result.V2 < 0


def test_table_value_M3200_c05() -> None:
    result = variance_components(3200, 0.5)
    assert math.isclose(result.rescaled_variance, 0.9826, rel_tol=0, abs_tol=5e-5)


def test_limit_constant() -> None:
    result = variance_components(100, 2.0)
    expected = 1 / math.expm1(2.0) - 1 / math.expm1(4.0)
    assert math.isclose(result.B_c, expected, rel_tol=0, abs_tol=2e-16)
