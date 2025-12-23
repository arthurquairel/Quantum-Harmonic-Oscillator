# Quantum Harmonic Oscillator & Potential Wells Solver ⚛️

![Python](https://img.shields.io/badge/Python-3.x-blue?style=flat&logo=python)
![NumPy](https://img.shields.io/badge/NumPy-Scientific%20Computing-blue?style=flat&logo=numpy)
![Physics](https://img.shields.io/badge/Physics-Quantum%20Mechanics-purple)
![Status](https://img.shields.io/badge/Status-Completed-success)

## 📖 Overview

This repository contains a numerical implementation for solving the **1D Schrödinger equation** in various potential configurations. The project focuses on the **Quantum Harmonic Oscillator** confined within an infinite potential well, as well as the study of **Finite Potential Wells** in semiconductor heterostructures.

The core approach involves using the eigenfunctions of the infinite well as a basis to diagonalize the Hamiltonian matrix of more complex systems (F. Marsiglio's method).

---

## 📂 Project Structure & Theoretical Background

The project is divided into two main physical studies:

### 1. Harmonic Oscillator in an Infinite Well
We analyze a particle subject to a harmonic potential $V_{OH}(x) = \frac{1}{2}m\omega^2x^2$ confined within an infinite well of width $a$.

* **Basis Expansion:** The wavefunctions $|\psi\rangle$ are expanded on the basis of the infinite well eigenstates $|\phi_n\rangle$.
* **Matrix Diagonalization:** We compute the Hamiltonian matrix elements $H_{nm} = \langle \phi_n | \hat{H} | \phi_m \rangle$.
* **The $R$ Parameter:** We study the system behavior based on the ratio $R = \frac{\hbar\omega}{E_1^0}$, which dictates the dominance between the harmonic potential and the well boundaries.
* **Numerical Cut-off:** The infinite basis is truncated to the first $N$ modes (e.g., $N=100$) to allow numerical resolution.

### 2. Finite Potential Well (Semiconductor Heterostructures)
We extend the method to a finite potential well structure, modeling electron confinement in semiconductors.

* **Potential Profile:** A well of width $2b$ and depth $V_0$.
* **Bound States:** The algorithm specifically searches for states with energy $E < 0$ (bound states).

---

## 🛠 Technologies & Libraries

* **Language:** Python 3.x
* **Linear Algebra:** `numpy.linalg` (Eigenvalue decomposition)
* **Integration:** `scipy.integrate` (Computing matrix elements $H_{nm}$)
* **Visualization:** `matplotlib.pyplot` (Wavefunctions, Potentials, and Energy levels)

---


## 📊 Key Findings & Results

The numerical simulations yielded the following physical insights, confirming theoretical predictions:

### 1. Matrix Sparsity & Parity
The Hamiltonian matrix $H$ exhibits a checkerboard pattern where every other off-diagonal term is zero. [cite_start]This numerically confirms the **decoupling of even and odd parity states**.

### 2. Energy Regimes & The Constant $C$
[cite_start]For the simulation with parameter $R=24$:
* [cite_start]**Harmonic Regime:** For $n \le 15$, the eigenvalues match the standard harmonic oscillator levels.
* [cite_start]**Infinite Well Regime:** For $n > 15$, the energies diverge from the harmonic ladder and follow $n^2 + C$.
* [cite_start]**Physical Meaning:** We determined numerically that $C \approx 118.34$, which matches the **average value of the reduced harmonic potential**.

### 3. Semiconductor Finite Well
[cite_start]For the finite well configuration ($a=1nm$, $b=3nm$, $V_0=0.5eV$), the simulation successfully isolated exactly **3 bound states** ($E < 0$).
[cite_start]The computed eigenenergies are:
1.  **Ground State:** $E_1 \approx -0.4437 \text{ eV}$
2.  **1st Excited State:** $E_2 \approx -0.2817 \text{ eV}$
3.  **2nd Excited State:** $E_3 \approx -0.0507 \text{ eV}$


   ## 🚀 Installation & Usage

### Prerequisites
Ensure you have the required scientific Python libraries installed:

```bash
pip install numpy scipy matplotlib*
