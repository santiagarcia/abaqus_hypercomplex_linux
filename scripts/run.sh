#!/bin/bash
# scripts/run.sh

# Ensure we are in the scripts directory
cd "$(dirname "$0")"

# Environment check
echo "Checking compilers..."
which ifort
which ifx
gfortran --version

# Abaqus Environment Info
echo "Checking Abaqus environment..."
abaqus information=environment
echo "Verifying Abaqus User Subroutines..."
abaqus verify -user_std

# Compilation of User Subroutine
echo "Compiling User Subroutine..."
# Create a joined source file because Abaqus likes single file inputs usually or we link objects.
# We will cat them together. Order matters: module first.
# cat ../src/hyperdual_mod.f90 ../src/beam_uel.f90 > beam_uel_complete.f90

# Run Abaqus Job
echo "Running Abaqus Datacheck..."
cd ../sim
abaqus job=cantilever_uel user=../src/beam_uel.f90 datacheck interactive ask_delete=OFF

echo "Running Abaqus Job relative to sim/..."
abaqus job=cantilever_uel user=../src/beam_uel.f90 interactive ask_delete=OFF

echo "Generating Visualization..."
python3 ../scripts/visualize_results.py

echo "Done."
