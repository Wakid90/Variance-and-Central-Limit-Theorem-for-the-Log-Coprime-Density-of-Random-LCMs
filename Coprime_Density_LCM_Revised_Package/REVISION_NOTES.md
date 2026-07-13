# Revision notes

This revision addresses the paper-only findings from the audit.

## Mathematical and literature corrections

- Corrected the classical two-variable expected-LCM asymptotic from \(\zeta(3)n^2/\zeta(2)\) to \(\zeta(3)n^2/(4\zeta(2))\).
- Corrected the Bostan--Marynych--Raschel normalization from \(n^{r_k}\) to \(n^{rk}\).
- Added the directly relevant Fernández--Fernández 2013 random gcd/lcm reference.
- Qualified the extension to general sampling distributions: the limiting product contains only primes that have positive probability of appearing.
- Replaced the inaccurate description of \(q^{**}/q^*\) as a ratio of two single-prime densities with its correct multi-prime/single-prime interpretation.
- Clarified the real-exponent definition of \(\Sigma_s\), including the \(s=0\) convention, and corrected the endpoint range in the geometric-mean proof.

## Numerical and reproducibility corrections

- Replaced claims of exact floating-point computation with precise descriptions of deterministic double-precision evaluation.
- Updated the Figure 2 caption to match the adaptive truncation used by the released code.
- Removed the unarchived small-case Monte Carlo agreement claim while retaining the exact 24-case rational verification.
- Corrected Table 2 wording: the displayed quantity is the geometric rate factor, while the theorem's full bound includes an \(M\)-dependent prefactor.
- Updated the theorem-audit numbering and changed its Theorem 2 test to check the bound stated in the paper.
- Updated the build instructions to three LaTeX passes.

The revised theorem audit reports 54 passed checks and 0 failures.
