# Unimodal Optimization Playground (Svenn + classic 1D search)

This mini project demonstrates **bracketing and 1D minimization** for a unimodal function:
- **Svenn’s method** to bracket an uncertainty interval `[a, b]`;
- **Dichotomy** (bisection for minimization);
- **Interval halving**;
- **Golden-section search**;
- **Adaptive-step search**.

The target function used in examples is:

\[
f(x) = 2 - \frac{1}{\log_2(x^4 + 2x^2 + 2)}, \quad
\text{with } x^4 + 2x^2 + 2 > 0 \text{ and } \neq 1.
\]

A small wrapper `FCounter` counts the number of function evaluations.

---

## Requirements

- Python 3.10+
- No external packages required (only the standard library).

---

## Run

```bash
python unimodal_search.py
