# Lepage Analysis: A Renormalization Approach to the Schrödinger Equation

## 📌 Overview

This project solves the time-independent radial Schrödinger equation using the Numerov method on a logarithmic grid. The physical system under study is a hydrogen-like atom perturbed by a short-range interaction. The main goal is to investigate how different approximations — including effective theories derived from renormalization techniques — reproduce the phase shifts and energy eigenvalues of a more "realistic" potential.

The inspiration comes from renormalization techniques explained by G. P. Lepage in [this lecture](https://arxiv.org/abs/nucl-th/9706029).

---

## 🚀 Getting Started

### Requirements

This code is written in Python 3 and depends on the following libraries:

- numpy
- pandas
- matplotlib
- scipy
- os, sys

Install the required dependencies via pip:

```bash
pip install numpy pandas matplotlib scipy
```

The code auto-generates several directories (`potentials/`, `eigenfunctions/`, `energy/`, `phase_shift/`, `relative_errors/`) to store intermediate results and plots.

---

## 📚 Libraries Used and Their Role

| Library        | Purpose                                                                 |
|----------------|-------------------------------------------------------------------------|
| numpy          | Efficient array operations and mathematical functions                   |
| pandas         | Saving and organizing numerical data into .csv files                    |
| matplotlib     | Plotting wavefunctions, potentials, and relative errors                 |
| scipy.special.erf | Used in smoothing the Coulomb potential (cutoff regularization)         |
| scipy.optimize.minimize | Minimization of cost functions to tune parameters in effective theory |
| os, sys        | Directory creation and clean error messages                             |

---

## 🎯 Project Objective

The aim is to:

1. Numerically solve the radial Schrödinger equation for:
   - Pure Coulomb potential.
   - Coulomb + short-range exponential potential.
   - Regularized Coulomb + smeared delta function (effective theory).
2. Compare eigenvalues and phase shifts from different potentials.
3. Extract physical insight from how effective potentials (a², a⁴ truncations) reproduce low-energy observables.

---

## ⚙️ How the Code Works

The code executes the following main steps:

1. Grid Initialization
   - Uses a logarithmic grid x = log(r) to better handle the singularity at r = 0 while keeping constant spacing for Numerov.

2. Potential Setup
   - Three main potentials:
     - V_coulomb = -1/r
     - V_real = -1/r - e^{-r}/r (used to generate "true" data)
     - V_eff: regularized potential with parameters tuned via phase shifts.

3. Numerov Integration
   - Solves the Schrödinger equation using an outward + inward integration scheme.
   - The eigenvalue is adjusted iteratively to match the number of nodes (n-l-1) and phase shift smoothness.

4. Effective Theory
   - Adds cutoff and smeared delta terms (up to order a⁴) to mimic the short-range physics.
   - Parameters are tuned via minimization of phase shift discrepancy (cost_function_1, cost_function_2).

5. Plots and Data Output
   - Saves energy levels, phase shifts, eigenfunctions, and various comparisons in .csv and .pdf.

### 🔧 Modifiable Parameters

- r_max, x_min, dx: control grid resolution and range.
- n_max: maximum quantum number n to compute.
- a: cutoff parameter for effective theory.
- Initial guesses in c_a2, init_parameters.

---

## 🧠 Inspiration and References

This project is inspired by:

> 📌 **[Insert exact course, professor, or assignment description here — e.g., "Lecture notes by Prof. XYZ for the TNANP course"]**

Core references:
- G. P. Lepage, "How to Renormalize the Schrödinger Equation" (https://arxiv.org/abs/nucl-th/9706029)
- Paolo Giannozzi, Numerical Methods in Quantum Mechanics, Lecture Notes (2021)

---

## 📊 Results & Known Issues

- The a² and a⁴ effective theories both approximate the full potential very well.
- However, the difference between a² and a⁴ results is extremely small, even though higher-order corrections were expected to improve the results.
- For high-energy states (lower n), the short-range potential has a more significant impact.
- Relative errors are lowest at small binding energy (high n), as expected from the renormalization framework.

---

## 🛠️ Future Improvements

- Higher-order contact terms: Include a⁶ and beyond in the effective theory for better accuracy at smaller cutoff.
- Generic angular momentum l ≠ 0: Currently only l = 0 (S-wave) is implemented.
- More robust root-finding: Current implementation uses bisection with perturbative refinement; may be improved with eigenvalue bracketing strategies.
- Modular design: Functions could be refactored into modules/classes for cleaner structure.

---

## 🧾 License

This project is for academic and educational purposes.

---

*Author: Caterina Vagnoni*
