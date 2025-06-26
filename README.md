# Lepage Analysis: A Renormalization Approach to the Schrödinger Equation

## Overview and Project Objective

This project solves the time-independent radial Schrödinger equation using the Numerov method on a logarithmic grid. The physical system under study is a hydrogen-like atom perturbed by a short-range interaction. The main goal is to investigate how different approximations — including effective theories derived from renormalization techniques — can reproduce the phase shifts and energy data of an unknown short-range potential.

The inspiration comes from renormalization techniques explained by G. P. Lepage in [this lectures](https://arxiv.org/abs/nucl-th/9706029).

---

## Getting Started

### Requirements

This code is written in Python and depends on the following libraries, imported at the beginning of the program:

- numpy
- pandas
- matplotlib
- scipy
- os, sys

The code auto-generates several directories (`potentials/`, `eigenfunctions/`, `energy/`, `phase_shift/`, `relative_errors/`) to store results and plots.

---

### Libraries Used and Their Role

| Library        | Purpose                                                                 |
|----------------|-------------------------------------------------------------------------|
| numpy          | Efficient array operations and mathematical functions                   |
| pandas         | Organizing numerical data into readable tables to print and/or save to .csv files,                    |
| matplotlib     | Plotting wavefunctions, potentials, and relative errors                 |
| scipy.special.erf | Used in smoothing the Coulomb potential (cutoff regularization)         |
| scipy.optimize.minimize | Minimization of cost functions to tune parameters in effective theory |
| os, sys        | Directory creation and clean error messages                             |

---

## Project Structure

The program is divided into various cells which can be run separately.
The main executive cells aim to:
1. Numerically solve the radial Schrödinger equation for:
   - Pure Coulomb potential (estimate the intrinsic error of the numerical method)
   - Coulomb + short-range exponential potential.
   - Regularized Coulomb + smeared delta function (effective theory).
2. Compare eigenvalues and phase shifts from different potentials.
They exploit functions and parameters defined in the preceding cells of the program.
The global parameters are defined and commented at the beginning of the program, in the cell named "Choice of the Global Parametes", and can be modified as the user wishes.


---

### Main Functions

1. Grid Initialization
   - Uses a logarithmic grid x = log(r) to better handle the singularity at r = 0 while keeping constant spacing for Numerov.

2. Potential Setup
   - Three main potentials are defined, and saved as .csv and .pdf files:
     - V_coulomb = -1/r
     - V_real = -1/r - e^{-r}/r (used to generate "true" data)
     - V_eff: regularized potential with parameters tuned via phase shifts.

3. Solving Schrödinger equation using Numerov integration
   - Solves the Schrödinger equation using an outward + inward integration scheme.
   - The eigenvalue is adjusted iteratively to match the number of nodes (n-l-1) and wavefunction smoothness.
    
4. Analysis of the given potential
   - Solves the Schrödinger equation for the energy levels under consideration.
   - Prints the results for the energies and phase shifts.

5. Cost functions
   - Define the quantity to minimize in order to tune the effective theories.
   

A more detailed analysis of the code is carried out in the file TNANPpresentation.pdf present in this repository.

### Parameters the User May Wish to Modify
- r_max, x_min, dx: control grid resolution and range.
- n_max: maximum quantum number n to compute.
- a: cutoff parameter for effective theory.


---

## Inspiration and References

This project is realised for the course "Theoretical and Numerical Aspects of Nuclear Physics", by Prof. Paolo Finelli, Alma Mater Studiorum, Bologna.

Core references:
- For the theory: G. P. Lepage, "How to Renormalize the Schrödinger Equation" (https://arxiv.org/abs/nucl-th/9706029)
- For the code: inspired from c++ coding pieces from Paolo Giannozzi, Numerical Methods in Quantum Mechanics, Lecture Notes (2021)

---

## Results & Known Issues

- The a² and a⁴ effective theories both approximate the full potential well.
- However, the difference between a² and a⁴ results is extremely small, even though higher-order corrections were expected to improve the results.
- For high-energy states (lower n), the short-range potential has a more significant impact and the relative errors are bigger. This result was, however, expected and is a limitation of the renormalization method.
- Relative errors are lowest at small binding energy (high n), as expected from the renormalization framework.
- Relative errors vary also at the varying of the cutoff parameter, as the user may verify and as shown in the TNANPpresentation.pdf, where all important results are plotted.

---

## Possible Future Improvements

- Higher-order contact terms: Include a⁶ and beyond in the effective theory for better accuracy at smaller cutoff.
- Generic angular momentum l ≠ 0: Currently only l = 0 (S-wave) is implemented.

---


