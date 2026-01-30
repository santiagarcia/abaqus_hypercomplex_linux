#!/bin/bash
# scripts/run_comparison.sh

cd "$(dirname "$0")"

if [ "$#" -ge 2 ]; then
    NUM_ELEMS=$1
    FORCE_VAL=$2
    if [ "$#" -eq 3 ]; then
        NUM_INCS=$3
    else
        NUM_INCS=10
    fi
else
    NUM_ELEMS=1
    FORCE_VAL=-1000.0
    NUM_INCS=10
fi

echo "============================================="
echo "Running Comparison: Hyperdual vs. Finite Diff"
echo "Elements: $NUM_ELEMS | Force: $FORCE_VAL | Increments: $NUM_INCS"
echo "============================================="

echo "Generating Input Files..."
python3 ../scripts/generate_inp.py ../sim/cantilever_uel.inp $NUM_ELEMS $FORCE_VAL $NUM_INCS
# FD uses the same geometry, just copy or generate again with same params
cp ../sim/cantilever_uel.inp ../sim/cantilever_fd.inp
# Fix HEADING in FD file
sed -i 's/Cantilever Beam/Cantilever Beam FD/' ../sim/cantilever_fd.inp

# 1. Run Hyperdual UEL
echo "[1/2] Running Hyperdual UEL..."
cd ../sim
abaqus job=cantilever_uel user=../src/beam_uel.f90 interactive ask_delete=OFF > uel_run.log 2>&1

# 2. Run Finite Difference UEL
echo "[2/3] Running FD UEL..."
abaqus job=cantilever_fd user=../src/beam_uel_fd.f90 interactive ask_delete=OFF > fd_run.log 2>&1

cp ../sim/cantilever_uel.inp ../sim/cantilever_complex.inp
sed -i 's/Cantilever Beam/Cantilever Beam Complex/' ../sim/cantilever_complex.inp

# 3. Run Hypercomplex UEL
echo "[3/4] Running Hypercomplex UEL..."
abaqus job=cantilever_complex user=../src/beam_uel_complex.f90 interactive ask_delete=OFF > complex_run.log 2>&1

cp ../sim/cantilever_uel.inp ../sim/cantilever_exact.inp
sed -i 's/Cantilever Beam/Cantilever Beam Exact/' ../sim/cantilever_exact.inp

# 4. Run Exact UEL
echo "[4/4] Running Exact Analytical UEL..."
abaqus job=cantilever_exact user=../src/beam_uel_exact.f90 interactive ask_delete=OFF > exact_run.log 2>&1

echo "============================================="
echo "Done. Checking Logs..."

# Minimal Timing Check
echo "Hyperdual Analysis Time:"
grep "TOTAL CPU TIME" cantilever_uel.dat | tail -n 1
echo "FD Analysis Time:"
grep "TOTAL CPU TIME" cantilever_fd.dat | tail -n 1
echo "Hypercomplex Analysis Time:"
grep "TOTAL CPU TIME" cantilever_complex.dat | tail -n 1
echo "Exact Analysis Time:"
grep "TOTAL CPU TIME" cantilever_exact.dat | tail -n 1

echo "Generating Comparison Visualization..."

python3 ../scripts/visualize_comparison.py
