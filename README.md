# Variance and Central Limit Theorem for the Log Coprime Density of Random LCMs

This repository contains the revised manuscript, LaTeX source, figures, exact finite-variance calculations, Monte Carlo diagnostics, and machine-readable numerical outputs for:

**Muhamad Wakid, _Variance and Central Limit Theorem for the Log Coprime Density of Random LCMs_ (2026).**

The companion first-order paper is:

**Muhamad Wakid, _On the Expected Coprime Density of the Least Common Multiple of Random Integers_, Zenodo (2026), DOI: 10.5281/zenodo.21325880.**

## Repository layout

- `paper/main.tex` - revised LaTeX manuscript
- `paper/main.pdf` - compiled manuscript
- `paper/figures/fig_variance.pdf` - deterministic variance figure
- `paper/figures/fig_clt_distribution.png` - Monte Carlo CLT diagnostic
- `code/exact_variance.py` - exact finite formula for `V1`, `V2`, and the total variance
- `code/make_variance_figure.py` - regenerates Figure 1 and Table 1 data
- `code/simulate_clt.py` - regenerates Figure 2 and its diagnostics
- `code/reproduce_all.py` - reproduces all numerical outputs
- `data/variance_table.csv` - machine-readable values behind Table 1
- `data/simulation_diagnostics.csv` - sample means, standard deviations, skewness, and excess kurtosis
- `tests/test_variance.py` - numerical regression tests

## Reproduce the numerical results

Python 3.11 or newer is recommended.

```bash
python -m venv .venv
```

Linux or macOS:

```bash
source .venv/bin/activate
```

Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Install dependencies and regenerate everything:

```bash
python -m pip install --upgrade pip
pip install -r requirements.txt
python code/reproduce_all.py
```

The script rewrites:

```text
paper/figures/fig_variance.pdf
paper/figures/fig_clt_distribution.png
data/variance_table.csv
data/simulation_diagnostics.csv
```

## Arbitrary-precision validation

The table grid is evaluated with NumPy long-double arithmetic. A parameter pair can be independently checked with `mpmath`:

```bash
python code/exact_variance.py --M 100 --c 1 --mpmath --dps 80
```

Expected rescaled variance:

```text
0.387682532746078228541128165354
```

## Compile the paper

A TeX installation containing `orcidlink`, `mathpazo`, `hyperref`, `placeins`, and the AMS packages is required.

```bash
cd paper
latexmk -pdf main.tex
```

## Run tests

```bash
pytest -q
```

## Numerical conventions

- The exact finite variance is evaluated from the full prime-pair formula, including pairs with `pq <= M`.
- `V2` is the ordered off-diagonal sum, so each unordered prime pair contributes twice.
- Monte Carlo panels use `numpy.random.default_rng(2026)`, freshly initialised for each panel.
- Standardisation uses the sample mean and sample standard deviation with `ddof=1`.
- Reported skewness and excess kurtosis are descriptive moment statistics, not hypothesis tests.

## Public release checklist

Before creating the public GitHub repository:

1. Choose a repository name.
2. Add the preferred code and manuscript licenses.
3. Push the complete directory.
4. Create a GitHub release.
5. Connect the repository to Zenodo and archive the release.
6. Insert the resulting repository DOI into the manuscript's data-and-code statement in a later version.
