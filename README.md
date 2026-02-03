# Abaqus Hypercomplex UEL Benchmarking

This repository contains a comprehensive framework for benchmarking different **Forward-Mode Automatic Differentiation** and **Numerical Differentiation** methods applied to the calculation of the **Tangent Stiffness Matrix** in Abaqus User Elements (UELs).

The primary test case is a **3D Large-Deformation Beam** (formulated with Hermite Cubic shape functions for bending and Linear shape functions for axial/torsion components).

## Overview

In nonlinear Finite Element Analysis, the convergence rate of the Newton-Raphson solver depends critically on the accuracy of the Tangent Stiffness Matrix (Jacobian). This project implements and compares five strategies for computing this matrix:

1.  **Exact Analytical (`Exact`)**:
    *   Manually derived Jacobian.
    *   Serves as the "Ground Truth" for accuracy and performance speed.
2.  **Hyperdual Numbers (`HD`)**:
    *   Uses a custom `hyperdual_mod` Fortran module.
    *   Second-derivative tracking AD.
    *   **Optimization**: Computes 2 columns of the stiffness matrix simultaneously per residual pass using dual imaginary channels (`e1`, `e2`).
3.  **Complex Step (`CS`)**:
    *   Uses the complex-step derivative approximation ($f'(x) \approx \text{Im}(f(x+ih))/h$).
    *   Avoids subtractive cancellation errors common in Finite Difference.
    *   Highly robust and often as fast as analytical methods for complex formulations.
4.  **Finite Difference (`FD`)**:
    *   Traditional numerical perturbation.
    *   Susceptible to rounding errors (step size dilemma).
5.  **Otis (`Otis`)**:
    *   A generalized multicomplex AD implementation using the `OTIM6N1` module.
    *   Supports higher-order derivatives and multiple directions (up to 6).

## Directory Structure

```plaintext
.
├── scripts/                # Python & Bash automation tools
│   ├── run_comparison.sh   # MAIN ENTRY POINT: Orchestrates the entire benchmark
│   ├── generate_inp.py     # Generates Abaqus input files (.inp)
│   ├── extract_tip_history.py # Parses .odb/.dat files for displacement results
│   └── visualize_comparison.py # Generates interactive HTML reports
├── src/                    # Fortran Source Code
│   ├── beam_uel.f90        # Hyperdual implementation
│   ├── beam_uel_fd.f90     # Finite Difference implementation
│   ├── beam_uel_complex.f90# Complex Step implementation
│   ├── beam_uel_exact.f90  # Analytic implementation
│   ├── beam_uel_otis.f90   # Generalized AD implementation
│   └── ...                 # Helper modules (real_utils, master_parameters)
├── sim/                    # Simulation Artifacts (generated at runtime)
│   ├── results/            # CSV logs and summary stats
│   ├── plots/              # HTML animations and graphs
│   └── runs/               # Abaqus working directories
└── README.md
```

## Prerequisites

*   **Abaqus Standard** (Solver and Fortran compiler setup).
*   **Python 3.x**:
    *   `numpy`
    *   `plotly` (for visualization)
*   **Bash** (WSL/Linux/Git Bash for Windows).

## Usage

The `run_comparison.sh` script handles compilation, execution, and post-processing in one go.

### Standard Run
Run a simulation with 10 elements and a tip force of -10,000.

```bash
./scripts/run_comparison.sh 10 -10000.0
```

### Advanced Run
Run with customized element count, force, number of increments, and a custom tag for the output folder.

```bash
# Syntax: ./scripts/run_comparison.sh <n_elems> <force> <n_increments> <tag>
./scripts/run_comparison.sh 20 -50000.0 50 "high_load_test"
```

## Outputs

After a successful run, check the `sim/` directory:

1.  ** Interactive Report**: `sim/beam_comparison.html`
    *   **3D Animation**: Visualizes the beam deformation for all methods overlaid.
    *   **Error Plots**: L2 Norm error vs. Step time relative to the Exact solution.
    *   **Runtime Analysis**: CPU time comparisons.

2.  **Data Logs**:
    *   `sim/timings.csv`: Aggregated execution times for performance analysis.
    *   `sim/tip_history.csv`: Displacement history of the beam tip.

## Formulation Details

*   **Element Type**: 2-Node 3D Beam (12 DOFs).
*   **Interpolation**:
    *   Axial/Torsion: Linear.
    *   Bending (v, w): Hermite Cubic.
*   **Material**: Linear Elastic (E=210 GPa, $\nu$=0.3).
*   **Geometric Nonlinearity**: `NLGEOM=YES` (Large rotation/displacement).

## License

This project is open-source. Please attribute the original author (santiagarcia) if used for derivative works.
