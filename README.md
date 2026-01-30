# Abaqus Hypercomplex UEL Comparison Project

This project implements and compares four different differentiation methods for calculating the Tangent Stiffness Matrix in an Abaqus User Element (UEL) subroutine (Fortran).

## Implemented Methods
1.  **Exact Analytical**: Manually derived Jacobian (Fastest & Reference).
2.  **Complex Step**: Uses a complex imaginary perturbation step (Robust & Fast).
3.  **Finite Difference (FD)**: Traditional numerical perturbation. Includes toggle for Forward vs Central difference.
4.  **Hyperdual Numbers**: Second-derivative tracking automatic differentiation (Most Robust, Slower).

## Running Simulations

Use the `run_comparison.sh` script to compile the UELs, run the Abaqus simulations, and generate comparison plots.

### Usage
```bash
./scripts/run_comparison.sh <number_of_elements> <force_magnitude> [number_of_increments]
```

### Examples
Run with 10 elements:
```bash
./scripts/run_comparison.sh 10 -10000.0
```

Run a larger model with 100 elements and 100 increments:
```bash
./scripts/run_comparison.sh 100 -100000.0 100
```

## Results
The script generates an interactive HTML comparison report:
*   `sim/beam_comparison.html`: Contains 3D beam deflection animations and error analysis plots.

## Project Structure
*   `src/`: Fortran UEL source code (`beam_uel.f90`, `beam_uel_exact.f90`, etc.).
*   `scripts/`: Python generation, visualization, and shell runner scripts.
*   `sim/`: Abaqus working directory (logs, input files, output databases).
