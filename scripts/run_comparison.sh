#!/bin/bash
# scripts/run_comparison.sh

# Resolve absolute path to the script and repo root
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
SIM_DIR="$REPO_ROOT/sim"
SRC_DIR="$REPO_ROOT/src"

# Argument Parsing
if [ "$#" -ge 2 ]; then
    NUM_ELEMS=$1
    FORCE_VAL=$2
    if [ "$#" -ge 3 ]; then
        NUM_INCS=$3
    else
        NUM_INCS=10
    fi
    if [ "$#" -ge 4 ]; then
        TAG=$4
    else
        TAG="default"
    fi
else
    # Defaults
    NUM_ELEMS=1
    FORCE_VAL=-1000.0
    NUM_INCS=10
    TAG="default"
fi

# Define Run Directory
if [ "$TAG" == "default" ]; then
    RUN_DIR="$SIM_DIR"
else
    RUN_DIR="$SIM_DIR/runs/$TAG"
fi

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

echo "============================================="
echo "Running Comparison: Hyperdual vs. FD vs. Complex vs. Exact"
echo "Elements: $NUM_ELEMS | Force: $FORCE_VAL | Increments: $NUM_INCS | Tag: $TAG"
echo "Dir: $RUN_DIR"
echo "============================================="

# Copy env file if needed (for UEL compilation flags)
if [ -f "$SIM_DIR/abaqus_v6.env" ]; then
    cp "$SIM_DIR/abaqus_v6.env" .
fi

# Define Base Names
JOB_UEL="cantilever_${TAG}_uel"
JOB_FD="cantilever_${TAG}_fd"
JOB_CS="cantilever_${TAG}_complex"
JOB_EXACT="cantilever_${TAG}_exact"
JOB_OTIS="cantilever_${TAG}_otis"

# Copy OTIM module source to run dir for compilation dependencies
if [ -f "$SRC_DIR/otim6n1.f90" ]; then
    cp "$SRC_DIR/otim6n1.f90" .
fi
if [ -f "$SRC_DIR/master_parameters.f90" ]; then
    cp "$SRC_DIR/master_parameters.f90" .
fi
if [ -f "$SRC_DIR/real_utils.f90" ]; then
    cp "$SRC_DIR/real_utils.f90" .
fi

# 1. Generate Input Files
echo "Generating Input Files..."
# We use the python script but pass the specific output filename
python3 "$SCRIPT_DIR/generate_inp.py" "${JOB_UEL}.inp" $NUM_ELEMS $FORCE_VAL $NUM_INCS

# FD Duplicate
cp "${JOB_UEL}.inp" "${JOB_FD}.inp"
sed -i "s/Cantilever Beam/Cantilever Beam FD/" "${JOB_FD}.inp"

# Complex Duplicate
cp "${JOB_UEL}.inp" "${JOB_CS}.inp"
sed -i "s/Cantilever Beam/Cantilever Beam Complex/" "${JOB_CS}.inp"

# Exact Duplicate
cp "${JOB_UEL}.inp" "${JOB_EXACT}.inp"
sed -i "s/Cantilever Beam/Cantilever Beam Exact/" "${JOB_EXACT}.inp"

# Otis Duplicate
cp "${JOB_UEL}.inp" "${JOB_OTIS}.inp"
sed -i "s/Cantilever Beam/Cantilever Beam Otis/" "${JOB_OTIS}.inp"

# 2. Run Jobs
# Note: abaqus uses 'user=' for the subroutine. We need absolute path or relative to RUN_DIR.
# SRC_DIR is absolute.

run_abaqus() {
    JOBNAME=$1
    USERFILE=$2
    echo "Running $JOBNAME..."
    # Clean old files
    rm -f ${JOBNAME}.lck ${JOBNAME}.odb ${JOBNAME}.dat ${JOBNAME}.msg ${JOBNAME}.sta ${JOBNAME}.prt ${JOBNAME}.com ${JOBNAME}.sim
    
    abaqus job=$JOBNAME user=$USERFILE interactive ask_delete=OFF > "${JOBNAME}.log" 2>&1
}

run_abaqus $JOB_UEL "$SRC_DIR/beam_uel.f90"
run_abaqus $JOB_FD "$SRC_DIR/beam_uel_fd.f90"
run_abaqus $JOB_CS "$SRC_DIR/beam_uel_complex.f90"
run_abaqus $JOB_EXACT "$SRC_DIR/beam_uel_exact.f90"
run_abaqus $JOB_OTIS "$SRC_DIR/beam_uel_otis.f90"


# 3. Extract Timings
extract_time() {
    JOB=$1
    METHOD=$2
    DATFILE="${JOB}.dat"
    TAG_ARG=$3
    NE=$4
    FV=$5
    NI=$6
    
    # Extract time using grep and awk. Matches: " TOTAL CPU TIME (SEC) =  1.20"
    if [ -f "$DATFILE" ]; then
        TIME=$(grep "TOTAL CPU TIME" "$DATFILE" | tail -n 1 | awk -F'=' '{print $2}' | xargs)
    else
        TIME=""
    fi
    
    if [ -z "$TIME" ]; then TIME="NaN"; fi
    
    # Append to global timings csv
    # Use flock if possible or just append (not concurrent here)
    echo "$TAG_ARG,$NE,$FV,$NI,$METHOD,$TIME" >> "$SIM_DIR/timings.csv"
}

# Init header if not exists
TIMINGS_FILE="$SIM_DIR/timings.csv"
if [ ! -f "$TIMINGS_FILE" ]; then
    echo "tag,nelems,force,nincs,method,cpu_time_seconds" > "$TIMINGS_FILE"
fi

extract_time $JOB_UEL "HD" "$TAG" $NUM_ELEMS $FORCE_VAL $NUM_INCS
extract_time $JOB_FD "FD" "$TAG" $NUM_ELEMS $FORCE_VAL $NUM_INCS
extract_time $JOB_CS "CS" "$TAG" $NUM_ELEMS $FORCE_VAL $NUM_INCS
extract_time $JOB_EXACT "Exact" "$TAG" $NUM_ELEMS $FORCE_VAL $NUM_INCS
extract_time $JOB_OTIS "Otis" "$TAG" $NUM_ELEMS $FORCE_VAL $NUM_INCS

# 4. Extract History Data (Tip Displacement)
echo "Running Extraction..."
abaqus python "$SCRIPT_DIR/extract_tip_history.py" "$RUN_DIR" "$TAG"

# 5. Generate Visualization (Plotly HTML)
echo "Running Visualization..."
# Ensure plotly is available or handle error
python3 "$SCRIPT_DIR/visualize_comparison.py" "$RUN_DIR" "$TAG"

echo "Done with run $TAG."
