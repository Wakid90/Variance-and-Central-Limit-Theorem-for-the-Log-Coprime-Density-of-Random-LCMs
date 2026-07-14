# Variance and Central Limit Theorem for the Log Coprime Density of Random LCMs

Manuscript, LaTeX source, figures, exact finite-variance calculations, Monte Carlo diagnostics, and machine-readable numerical outputs for:

**Muhamad Wakid, *Variance and Central Limit Theorem for the Log Coprime Density of Random LCMs* (2026).**
DOI: [10.5281/zenodo.21339769](https://doi.org/10.5281/zenodo.21339769)

The companion first-order paper is:

**Muhamad Wakid, *On the Expected Coprime Density of the Least Common Multiple of Random Integers*, Zenodo (2026).**
DOI: [10.5281/zenodo.21325880](https://doi.org/10.5281/zenodo.21325880) · [source repository](https://github.com/Wakid90/coprime-density-lcm)

## Results

Let `X_1, X_2, ...` be i.i.d. uniform on `{2, ..., M}`, let `L_k = lcm(X_1, ..., X_k)`, and let

```
W_k(M) = ln(phi(L_k)/L_k) - E ln(phi(L_k)/L_k)
```

be the centred log coprime density, with `s_k^2 = Var W_k(M)`. In the joint regime `k = floor(cM)`, `M -> infinity`, with `c > 0` fixed, the paper proves three results.

**1. Variance asymptotic.**

```
M ln M * s_k^2  ->  B(c) = 1/(e^c - 1) - 1/(e^{2c} - 1)
```

proved by a diagonal / off-diagonal decomposition, the off-diagonal term being controlled by an exact covariance identity valid for prime pairs with `pq > M`.

**2. Central limit theorem.** `W_k(M)/s_k -> N(0, 1)` in distribution, via the Hall–Heyde martingale CLT applied to the Doob filtration. The two nontrivial ingredients are a uniform increment bound combining a coupon-collector estimate with a sieve-geometric observation, and a conditional-variance concentration resting on explicit hitting-time moment bounds together with an `L^2` Minkowski summation.

**3. Berry–Esseen rate.** The Kolmogorov distance to the standard normal is `O_c((ln M)^{-2/5})`, on both the logarithmic and the linear scale, via the Heyde–Brown martingale inequality together with a quantitative refinement of result 1,

```
M ln M * s_k^2  =  B(c) + O_c(1 / ln M),
```

obtained from the effective prime-counting bounds of Rosser and Schoenfeld.

## Layout

```
paper/main.tex                       LaTeX manuscript
paper/main.pdf                       compiled manuscript
paper/figures/fig_variance.pdf       Figure 1 (deterministic)
paper/figures/fig_clt_distribution.png   Figure 2 (Monte Carlo)
code/exact_variance.py               exact finite formula for V1, V2, s_k^2
code/make_variance_figure.py         regenerates Figure 1 and Table 1
code/simulate_clt.py                 regenerates Figure 2 and its diagnostics
code/reproduce_all.py                regenerates everything, then validates
data/variance_table.csv              machine-readable Table 1
data/simulation_diagnostics.csv      sample moments behind Figure 2
tests/test_variance.py               numerical regression tests
```

## Reproduce the numerical results

Python 3.11 or newer.

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\Activate.ps1
pip install -r requirements.txt
python code/reproduce_all.py       # or: make repro
```

This rewrites both figures and both CSVs, then runs the arbitrary-precision validation described below. Runtime is a few minutes, dominated by the 80-digit checks.

## Numerical methodology

**Table 1 and Figure 1 are exact, not sampled.** They evaluate the finite formula `s_k^2 = V1 + V2` by direct summation over the primes up to `M` and over all ordered prime pairs, in extended-precision floating point (`numpy.longdouble`). No Monte Carlo is involved.

**Cancellation is checked, not assumed.** `V2` is a difference of nearly equal quantities, so `exact_variance.py` carries a second, independent implementation at 80 decimal digits (`mpmath`). `reproduce_all.py` compares the two across the grid and *fails loudly* on any disagreement above `1e-12`. Observed agreement is to within `1e-16`. A single pair can be checked directly:

```bash
python code/exact_variance.py --M 100 --c 1 --mpmath --dps 80
# 0.387682532746078228541128165354
```

**Figure 2 is Monte Carlo and seeded.** Each panel uses `numpy.random.default_rng(2026)`, freshly initialised, so the figure is bit-for-bit reproducible. Sample sizes are `N = 3000` at `M = 100` and `M = 300`, and `N = 1000` at `M = 1000`, with `k = M` throughout. Standardisation uses the sample mean and sample standard deviation with `ddof = 1`. Reported skewness and excess kurtosis (Fisher convention) are descriptive moment statistics, not hypothesis tests.

**Conventions.** `V2` is the ordered off-diagonal sum, so each unordered prime pair contributes twice. The full prime-pair formula is used, including pairs with `pq <= M`. `V2` is negative at every tested `(k, M)`, consistent with Remark 9 of the paper; the test suite asserts this.

## Tests

```bash
pytest -q          # or: make test
```

60 tests, covering: the prime sieve, the closed form for `B(c)`, exact reproduction of all 24 Table 1 entries, agreement between the extended-precision and 80-digit evaluations, the sign of `V2`, and its magnitude relative to `V1`.

## Compile the paper

Requires a TeX installation with `amsmath`, `amsthm`, `hyperref`, `mathpazo`, `helvet`, `courier`, `eso-pic`, `orcidlink`, and `placeins`.

```bash
cd paper && latexmk -pdf main.tex     # or: make paper
```

## Citation

```bibtex
@misc{wakid2026clt,
  author = {Wakid, Muhamad},
  title  = {Variance and Central Limit Theorem for the Log Coprime Density of Random {LCM}s},
  year   = {2026},
  doi    = {10.5281/zenodo.21339769},
  note   = {Preprint}
}
```

## License

Code: MIT. Manuscript and figures: CC BY 4.0. See `LICENSE`.
