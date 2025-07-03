# Lepage Analysis: A Renormalization Approach to the Schrödinger Equation

## Overview and Project Objective

This project solves the time-independent radial Schrödinger equation using the Numerov method on a logarithmic grid. The physical system under study is a hydrogen-like atom perturbed by a short-range interaction. The main goal is to investigate how different approximations — including effective theories derived from renormalization techniques — can reproduce the phase shifts and energy data of an unknown short-range potential.

The inspiration comes from renormalization techniques explained by G. P. Lepage in [these lectures](https://arxiv.org/abs/nucl-th/9706029).

---

## Getting Started

### Requirements

This code is written in Python and depends on the following libraries, imported at the beginning of the program:

- `numpy`
- `pandas`
- `matplotlib`
- `scipy`
- `os`, `sys`

The code auto-generates several directories (`potentials/`, `eigenfunctions/`, `energy/`, `phase_shift/`, `relative_errors/`) to store results and plots.

---

### Libraries Used and Their Role

| Library                   | Purpose                                                                 |
|--------------------------|-------------------------------------------------------------------------|
| `numpy`                  | Efficient array operations and mathematical functions                   |
| `pandas`                 | Organizing numerical data into readable tables to print and/or save to `.csv` files |
| `matplotlib`             | Plotting wavefunctions, potentials, and relative errors                 |
| `scipy.special.erf`      | Used in smoothing the Coulomb potential (cutoff regularization)         |
| `scipy.optimize.minimize`| Minimization of cost functions to tune parameters in effective theory   |
| `os`, `sys`              | Directory creation and clean error messages                             |
| `matplotlib.cm`          | Color map for the plots                                                 |

---

## Project Structure

The program is divided into various cells which can be run separately.
The main executive cells aim to:
1. Numerically solve the radial Schrödinger equation for:
   - Pure Coulomb potential (estimate the intrinsic error of the numerical method)
   - Coulomb + short-range exponential potential (real theory whose data we want to reproduce with the effective theory)
   - Regularized Coulomb + smeared delta function potential (effective theory)

2. Compare eigenvalues and phase shifts from different potentials through the study of relative errors with respect to the real theory.

They exploit functions and parameters defined in the preceding cells of the program, which the user must therefore run preliminarily.
The global parameters are defined and commented at the beginning of the program, in the cell named Choice of the Global Parameters, and can be freely modified by the user.

---

### Main Functions

1. **Grid Initialization**
   - Uses a logarithmic grid x = log(r) to better handle the singularity at `r = 0` while keeping constant spacing for Numerov.

2. **Potential Setup**
   - Three main potentials are defined, and saved as `.csv` and `.pdf` files:
     - `v_coulomb = -1/r`
     - `v_real` (used to generate "true" data)
     - `v_eff`: regularized potential with parameters tuned via phase shifts

   The user can modify the form of the short-range potential and of the effective potential as they wish.
   This must be done inside the functions which initialize the two potentials, in the first part of the code, namely     `short_range_initialization` and `effective_initialization`. One example for a certain choice of these potentials is carried out in the presentation TNANPpresentation.pdf.

3. **Solving Schrödinger Equation Using Numerov Integration**
   - Solves the Schrödinger equation using an outward + inward integration scheme
   - The eigenvalue is adjusted iteratively to match the number of nodes (`n - l - 1`) and wavefunction smoothness

4. **Analysis of the Given Potential**
   - Solves the Schrödinger equation for the energy levels under consideration
   - Prints the results for the energies and phase shifts

5. **Cost Functions**
   - Define the quantity to minimize in order to tune the effective theories

A more detailed analysis of the code is carried out in the file TNANPpresentation.pdf present in this repository.

### Parameters the User May Wish to Modify

- `r_max`, `x_min`, `dx`: control grid resolution and range
- `n_max`: maximum quantum number `n` to compute
- `a`: cutoff parameter for effective theory
- `a_values`: set of values for the cutoff if the user wishes to compare the results

---

## Changes Implemented in Different Versions of the Code

| Code             | `lepage_analysis.py`                                          | `lepage_analysis_new.py`                                      |
|---------------------|---------------------------------------------------------------|----------------------------------------------------------------|
| **Real Potential** (can be modified by the user)  | $V_{\text{real}}(r) = -\dfrac{1}{r} - \dfrac{e^{-r}}{r}$        | $V_{\text{real}}(r) = -\dfrac{1}{r} - \dfrac{1}{1 + e^r}$        |
| **Effective Potential** (can be modified by the user)  | $V_{\text{eff}}(r) = -\dfrac{1}{r}  \text{erf}\left( \dfrac{r}{\sqrt{2} a} \right) + c  a^2  \dfrac{ e^{-r/a}}{8 \pi a^3} + d_1  a^4  \left( \dfrac{r}{a} - 2 \right) \dfrac{  e^{-r/a}}{8 \pi a^4 r}$ | $V_{\text{eff}}(r) = -\dfrac{1}{r}  \text{erf}\left( \dfrac{r}{\sqrt{2} a} \right) + c  a^2 \dfrac{ e^{-r^2 / (2 a^2)}}{(2\pi)^{3/2} a^3} + d_1  a^4 \left( \dfrac {r^2} {a^2} - 3 \right) \dfrac{ e^{-r^2 / (2 a^2)}}{(2\pi)^{3/2} a^5}$ |
| **Improvements** | | - Improved plots of relative errors for the effective theory across cutoff values.<br>- Enhanced eigenvalue-finding algorithm in the Schrödinger solver for stability and accuracy. |

---

## Inspiration and References

This project is realised for the course "Theoretical and Numerical Aspects of Nuclear Physics", by Prof. Paolo Finelli, Alma Mater Studiorum, Bologna.

Core references:
- For the theory: G. P. Lepage, *How to Renormalize the Schrödinger Equation* ([arXiv:nucl-th/9706029](https://arxiv.org/abs/nucl-th/9706029))
- For the code: inspired from C++ coding pieces from Paolo Giannozzi, *Numerical Methods in Quantum Mechanics, Lecture Notes* (2021)

---

## Results & Known Issues

- The `a^2` and `a^4` effective theories both approximate the full potential well
- However, the difference between `a^2` and `a^4` results is extremely small, even though higher-order corrections were expected to improve the results
- For high-energy states (lower `n`), the short-range potential has a more significant impact and the relative errors are bigger. This result was, however, expected and is a limitation of the renormalization method
- Relative errors are lower at small binding energy (high `n`), as expected from the renormalization framework
- Relative errors vary also at the varying of the cutoff parameter, as the user may verify and as shown in the TNANPpresentation.pdf, where all important results are plotted as a reference example

---

## Possible Future Improvements

- Higher-order contact terms: Include `a^6` and beyond in the effective theory for better accuracy
- Generic angular momentum `l ≠ 0`: Currently only `l = 0` (S-wave) is implemented

