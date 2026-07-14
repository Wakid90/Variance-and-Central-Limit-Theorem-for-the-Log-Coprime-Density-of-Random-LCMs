"""Numerical regression tests for the exact variance computation.

Run with:  pytest -q
"""

from __future__ import annotations

import math
import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "code"))

from exact_variance import (  # noqa: E402
    B,
    primes_upto,
    rescaled_variance,
    rescaled_variance_mp,
    variance_parts,
)

# Table 1 of the paper: M ln M * s_k^2 at k = floor(cM).
TABLE_1 = {
    (100, 0.5): 0.8977, (100, 1.0): 0.3877, (100, 2.0): 0.1261, (100, 3.0): 0.0454,
    (200, 0.5): 0.9702, (200, 1.0): 0.4478, (200, 2.0): 0.1506, (200, 3.0): 0.0550,
    (400, 0.5): 0.9785, (400, 1.0): 0.4205, (400, 2.0): 0.1303, (400, 3.0): 0.0464,
    (800, 0.5): 0.9585, (800, 1.0): 0.4259, (800, 2.0): 0.1407, (800, 3.0): 0.0514,
    (1600, 0.5): 0.9767, (1600, 1.0): 0.4342, (1600, 2.0): 0.1422, (1600, 3.0): 0.0518,
    (3200, 0.5): 0.9826, (3200, 1.0): 0.4323, (3200, 2.0): 0.1403, (3200, 3.0): 0.0510,
}


def test_primes_upto():
    assert list(primes_upto(30)) == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    assert len(primes_upto(1000)) == 168  # pi(1000)


def test_B_closed_form():
    # B(c) = sum_{j>=1} e^{-cj}(1 - e^{-cj}), summed directly
    for c in (0.5, 1.0, 2.0, 3.0):
        series = sum(math.exp(-c * j) * (1 - math.exp(-c * j)) for j in range(1, 4000))
        assert abs(series - B(c)) < 1e-12


@pytest.mark.parametrize("key,expected", sorted(TABLE_1.items()))
def test_table_1_reproduces(key, expected):
    M, c = key
    total, _, _ = rescaled_variance(M, c)
    assert round(total, 4) == expected


@pytest.mark.parametrize("M,c", [(100, 0.5), (100, 1.0), (100, 2.0), (100, 3.0),
                                 (400, 0.5), (400, 1.0), (400, 2.0), (400, 3.0)])
def test_extended_precision_matches_80_digits(M, c):
    """Guards against catastrophic cancellation in V2, a difference of near-equal terms."""
    fast, _, _ = rescaled_variance(M, c)
    slow = float(rescaled_variance_mp(M, c, dps=80))
    assert abs(fast - slow) < 1e-12


@pytest.mark.parametrize("M", [100, 200, 400, 800])
@pytest.mark.parametrize("c", [0.5, 1.0, 2.0, 3.0])
def test_V2_is_negative(M, c):
    """Remark 9: the off-diagonal contribution is strictly negative."""
    k = int(math.floor(c * M))
    _, V2 = variance_parts(M, k)
    assert V2 < 0


@pytest.mark.parametrize("M", [100, 400])
@pytest.mark.parametrize("c", [0.5, 1.0, 2.0, 3.0])
def test_V2_is_an_order_of_magnitude_below_V1(M, c):
    k = int(math.floor(c * M))
    V1, V2 = variance_parts(M, k)
    assert abs(V2) < 0.2 * abs(V1)


def test_documented_spot_value():
    """The value quoted in the README."""
    val = float(rescaled_variance_mp(100, 1.0, dps=80))
    assert abs(val - 0.387682532746078228541128165354) < 1e-25


def test_variance_is_positive():
    for M in (100, 400, 1600):
        for c in (0.5, 1.0, 2.0, 3.0):
            total, _, _ = rescaled_variance(M, c)
            assert total > 0
