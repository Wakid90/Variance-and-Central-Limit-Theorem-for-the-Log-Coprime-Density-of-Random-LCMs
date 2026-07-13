# On the Expected Coprime Density of the Least Common Multiple of Random Integers

**Author:** Muhamad Wakid  
**ORCID:** [0009-0008-6274-778X](https://orcid.org/0009-0008-6274-778X)  
**Preprint date:** 12 July 2026  
**DOI:** [10.5281/zenodo.21325880](https://doi.org/10.5281/zenodo.21325880)

## Overview

This repository contains the LaTeX source, figures, numerical reproduction code, raw numerical output, and theorem-audit code for the preprint.

The paper studies

\[
\rho_k(M)=\mathbb E\!\left[\frac{\varphi(L_k)}{L_k}\right],
\qquad
L_k=\operatorname{lcm}(X_1,\ldots,X_k),
\]

where the variables are independent and uniform on \(\{2,\ldots,M\}\). It proves an exact subset identity, sharp fixed-support convergence, and the joint-scaling limit

\[
\log M\,\frac{\rho_{\lfloor cM\rfloor}(M)-P_M}{P_M}
\longrightarrow
A(c)=\sum_{j\ge1}e^{-cj}\log\!\left(1+\frac1j\right).
\]

## Repository contents

- `Coprime_Density_LCM_Preprint_Revised.tex` - revised LaTeX source.
- `Coprime_Density_LCM_Preprint_Revised.pdf` - compiled revised preprint.
- `fig1_convergence.pdf` - deterministic Monte Carlo convergence figure.
- `fig2_joint_scaling.pdf` - deterministic direct-prime-sum joint-scaling figure.
- `reproduce_results.py` - regenerates both figures and all three numerical tables.
- `figure1_points.csv` - raw points and standard errors used in Figure 1.
- `table1_monte_carlo.csv` - raw Table 1 estimates, seeds, and standard errors.
- `table2_fixed_M.csv` - exact fixed-\(M\) correction ratios used in Table 2.
- `table3_joint_scaling.csv` - exact joint-scaling values used in Table 3.
- `numerical_report.txt` - concise numerical reproduction log.
- `verify_theorems.py` - independent numerical audit of the main identities and bounds.
- `AUDIT_REPORT.txt` - output of the theorem audit.
- `requirements.txt` - Python dependencies.

## Compile the paper

A standard TeX Live installation with `pdflatex` is sufficient.

```bash
pdflatex -interaction=nonstopmode -halt-on-error Coprime_Density_LCM_Preprint_Revised.tex
pdflatex -interaction=nonstopmode -halt-on-error Coprime_Density_LCM_Preprint_Revised.tex
pdflatex -interaction=nonstopmode -halt-on-error Coprime_Density_LCM_Preprint_Revised.tex
```

The figures must remain in the same directory as the `.tex` file.

## Reproduce the numerical results

Create an environment and install the dependencies:

```bash
python -m venv .venv
source .venv/bin/activate            # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Regenerate the figures and tables:

```bash
python reproduce_results.py
```

Monte Carlo cells use NumPy's default generator with

```text
base seed = 2026
cell seed = 2026 + 1000*M + k
replicates = 2000
```

The Table 3 calculation sieves primes up to \(10^7\), so memory use can reach several hundred megabytes. Deterministic prime sums use double precision, and the positive series for `A(c)` is truncated when the next term is below `1e-17`.

## Run the theorem audit

```bash
python verify_theorems.py > AUDIT_REPORT.txt
```

Regenerate the audit report after any source or code change. The included revised report records **54 passed checks and 0 failures**. The audit checks exact finite cases of the subset identity, the fixed-\(M\) bounds, the formula for \(q^{**}(M)\), the joint-scaling numerics, and the asymptotics of \(A(c)\). These computational checks supplement, but do not replace, the proofs in the paper.

## Citation

```bibtex
@misc{wakid2026coprime,
  author       = {Muhamad Wakid},
  title        = {On the Expected Coprime Density of the Least Common Multiple of Random Integers},
  year         = {2026},
  month        = jul,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.21325880},
  url          = {https://doi.org/10.5281/zenodo.21325880},
  note         = {Preprint, 12 July 2026}
}
```
