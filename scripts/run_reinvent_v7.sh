#!/bin/bash
#SBATCH --job-name=REINVENT_V7
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=reinvent_v7_%j.out
#SBATCH --error=reinvent_v7_%j.err

set -euo pipefail

# Repository location
if [[ -n "${REINVENT_ROOT:-}" ]]; then
    ROOT="$REINVENT_ROOT"
elif [[ -n "${SLURM_JOB_ID:-}" ]]; then
    ROOT="${SLURM_SUBMIT_DIR:?Missing Slurm submission directory}"
else
    ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

ROOT="$(realpath "$ROOT")"
cd "$ROOT"

export REINVENT_ROOT="$ROOT"

# Reuse the V6 runtime configuration.
# Keep this separate from the commercial runtime JSON.
RUNTIME_CONFIG="${RUNTIME_CONFIG:-config/runtime.conf}"

if [[ ! -f "$RUNTIME_CONFIG" ]]; then
    RUNTIME_CONFIG="config/runtime.example.conf"
fi

RUNTIME_CONFIG="$(realpath "$RUNTIME_CONFIG")"

source "$RUNTIME_CONFIG"

export RUNTIME_CONFIG

resolve_python() {
    local env_name="$1"

    conda run --no-capture-output -n "$env_name" \
        python -c 'import sys; print(sys.executable)'
}

REINVENT_PY="${REINVENT_PY:-$(resolve_python "$REINVENT_ENV")}"
PYTHON="$REINVENT_PY"

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-${CPU_THREADS:-4}}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-${CPU_THREADS:-4}}

trap 'echo "V7 PIPELINE FAILED at line $LINENO" >&2' ERR

stage() {
    echo
    echo "=================================================="
    echo "$1"
    echo "Time: $(date)"
    echo "=================================================="
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

usage() {
    cat <<'EOF'

REINVENT V7 FLEXIBLE PIPELINE

FULL PIPELINE:

  bash scripts/run_reinvent_v7.sh full \
      INPUT N OUTPUT CMOS_PROFILE SELECT [OPTIONS]

COMMERCIAL-ONLY:

  bash scripts/run_reinvent_v7.sh commercial-only \
      SELECTED_CSV OUTPUT [OPTIONS]

COMMERCIAL OPTIONS:

  --profile PROFILE
      default, vendors, validated, all,
      enamine, chembridge, chemdiv, lifechemicals

  --libraries LIBRARY [LIBRARY ...]
      Select individual commercial sources.

  --skip-commercial
      Skip commercial screening.

  --mode both
      Integrated screening mode.

  --minimum-similarity FLOAT
      Default: 0.5

  --maximum-analogs INTEGER
      Default: 20

  --allow-remote
      Explicitly permit remote screening.

  --maximum-new-jobs INTEGER
      Required for remote screening.

  --help
      Show this help message.

DEFAULT BEHAVIOR:

  Libraries: Enamine + ChemBridge
  Screening mode: both
  Remote screening: disabled
  Existing cache: reused

EOF
}

# --------------------------------------------------
# Defaults
# --------------------------------------------------

COMMERCIAL_PROFILE=default
SCREENING_MODE=both
MINIMUM_SIMILARITY=0.5
MAXIMUM_ANALOGS=20

ALLOW_REMOTE=false
MAXIMUM_NEW_JOBS=0
SKIP_COMMERCIAL=false

PROFILE_EXPLICIT=false

LIBRARIES=()

MODE=${1:-help}

case "$MODE" in

    help|--help|-h)
        usage
        exit 0
        ;;

    full)
        shift

        INPUT=${1:?Provide input SMILES file}
        shift

        N=1000
        OUT=results/v7/REINVENT_V7_Results
        PROFILE=cns
        SELECT=30

        if [[ $# -gt 0 && "$1" != --* ]]; then
            N=$1
            shift
        fi

        if [[ $# -gt 0 && "$1" != --* ]]; then
            OUT=$1
            shift
        fi

        if [[ $# -gt 0 && "$1" != --* ]]; then
            PROFILE=$1
            shift
        fi

        if [[ $# -gt 0 && "$1" != --* ]]; then
            SELECT=$1
            shift
        fi

        [[ -s "$INPUT" ]] || \
            die "Input file missing or empty: $INPUT"

        INPUT=$(realpath "$INPUT")
        ;;

    commercial-only)
        shift

        SELECTED=${1:?Provide CMOS selected CSV}
        shift

        OUT=results/v7/Commercial_Results

        if [[ $# -gt 0 && "$1" != --* ]]; then
            OUT=$1
            shift
        fi

        [[ -s "$SELECTED" ]] || \
            die "Selected CSV missing or empty: $SELECTED"

        SELECTED=$(realpath "$SELECTED")
        ;;

    *)
        die "Unknown execution mode: $MODE"
        ;;

esac

# --------------------------------------------------
# Backward-compatible positional commercial profile
# --------------------------------------------------

if [[ $# -gt 0 && "$1" != --* ]]; then

    COMMERCIAL_PROFILE=$1
    PROFILE_EXPLICIT=true

    shift

fi

# --------------------------------------------------
# Flexible commercial options
# --------------------------------------------------

while [[ $# -gt 0 ]]; do

    case "$1" in

        --profile)
            [[ $# -ge 2 ]] || die "--profile requires a value"

            COMMERCIAL_PROFILE=$2
            PROFILE_EXPLICIT=true

            shift 2
            ;;

        --libraries)
            shift

            [[ $# -gt 0 && "$1" != --* ]] || \
                die "--libraries requires at least one vendor"

            while [[ $# -gt 0 && "$1" != --* ]]; do

                LIBRARIES+=("$1")
                shift

            done
            ;;

        --skip-commercial)
            SKIP_COMMERCIAL=true
            shift
            ;;

        --mode)
            [[ $# -ge 2 ]] || die "--mode requires a value"

            SCREENING_MODE=$2
            shift 2
            ;;

        --minimum-similarity)
            [[ $# -ge 2 ]] || \
                die "--minimum-similarity requires a value"

            MINIMUM_SIMILARITY=$2
            shift 2
            ;;

        --maximum-analogs)
            [[ $# -ge 2 ]] || \
                die "--maximum-analogs requires a value"

            MAXIMUM_ANALOGS=$2
            shift 2
            ;;

        --allow-remote)
            ALLOW_REMOTE=true
            shift
            ;;

        --maximum-new-jobs)
            [[ $# -ge 2 ]] || \
                die "--maximum-new-jobs requires a value"

            MAXIMUM_NEW_JOBS=$2
            shift 2
            ;;

        --help|-h)
            usage
            exit 0
            ;;

        *)
            die "Unknown option: $1"
            ;;

    esac

done

# --------------------------------------------------
# Validate execution options
# --------------------------------------------------

[[ -x "$PYTHON" ]] || die "Python executable not found"

if [[ "$MODE" == "full" ]]; then

    [[ "$N" =~ ^[1-9][0-9]*$ ]] || \
        die "N must be a positive integer"

    [[ "$SELECT" =~ ^[1-9][0-9]*$ ]] || \
        die "SELECT must be a positive integer"

fi

[[ "$MAXIMUM_ANALOGS" =~ ^(0|[1-9][0-9]*)$ ]] || \
    die "Maximum analogs must be a nonnegative integer"

[[ "$MAXIMUM_NEW_JOBS" =~ ^(0|[1-9][0-9]*)$ ]] || \
    die "Maximum new jobs must be a nonnegative integer"

if [[ "$SCREENING_MODE" != "both" ]]; then

    die "CMOS integration currently requires --mode both"

fi

if [[ "$PROFILE_EXPLICIT" == true && ${#LIBRARIES[@]} -gt 0 ]]; then

    die "Choose either --profile or --libraries, not both"

fi

if [[ "$ALLOW_REMOTE" == true && "$MAXIMUM_NEW_JOBS" -eq 0 ]]; then

    die "--allow-remote requires a positive --maximum-new-jobs"

fi

if [[ "$ALLOW_REMOTE" == false && "$MAXIMUM_NEW_JOBS" -gt 0 ]]; then

    die "--maximum-new-jobs requires --allow-remote"

fi

if [[ "$SKIP_COMMERCIAL" == true && "$ALLOW_REMOTE" == true ]]; then

    die "Cannot combine --skip-commercial with --allow-remote"

fi

if [[ "$MODE" == "commercial-only" && "$SKIP_COMMERCIAL" == true ]]; then

    die "--skip-commercial is only applicable to full mode"

fi

# --------------------------------------------------
# Output directories
# --------------------------------------------------

mkdir -p "$OUT"

OUT=$(realpath "$OUT")

COMMERCIAL="$OUT/commercial"
FINAL="$OUT/final"

BASE_CONFIG="$ROOT/config/commercial_sources.json"

COMMERCIAL_RUNTIME_CONFIG="$OUT/commercial_runtime_config.json"

# Prevent accidental reuse of previous commercial outputs.

if [[ "$SKIP_COMMERCIAL" == false ]]; then

    [[ ! -e "$COMMERCIAL" && ! -e "$FINAL" ]] || \
        die "Commercial output already exists. Choose a new output directory."

fi

# --------------------------------------------------
# Create and validate runtime configuration
# --------------------------------------------------

stage "REINVENT V7 PIPELINE INITIALIZATION"

echo "Execution mode: $MODE"
echo "Output directory: $OUT"
echo "Skip commercial: $SKIP_COMMERCIAL"
echo "Remote screening: $ALLOW_REMOTE"
echo "Maximum new jobs: $MAXIMUM_NEW_JOBS"

if [[ "$SKIP_COMMERCIAL" == false ]]; then

    "$PYTHON" - \
        "$ROOT" \
        "$BASE_CONFIG" \
        "$COMMERCIAL_RUNTIME_CONFIG" \
        "$COMMERCIAL_PROFILE" \
        "$ALLOW_REMOTE" \
        "$MAXIMUM_NEW_JOBS" \
        "$MINIMUM_SIMILARITY" \
        "$MAXIMUM_ANALOGS" \
        "${LIBRARIES[@]}" <<'PY'

import sys
import json
import math
from pathlib import Path

root = Path(sys.argv[1])
base = Path(sys.argv[2])
runtime = Path(sys.argv[3])

profile = sys.argv[4]
allow_remote = sys.argv[5] == "true"
maximum_new_jobs = int(sys.argv[6])

minimum_similarity = float(sys.argv[7])
maximum_analogs = int(sys.argv[8])

requested = list(dict.fromkeys(sys.argv[9:]))

sys.path.insert(
    0,
    str(root / "scripts/commercial"),
)

from commercial_config import CommercialConfig

original_bytes = base.read_bytes()

source = CommercialConfig(base)

config = json.loads(original_bytes)

if config["api"]["allow_remote_submission"] is not False:
    raise RuntimeError(
        "Production configuration must remain offline"
    )

if not math.isfinite(minimum_similarity):
    raise ValueError("Similarity must be finite")

if not 0 <= minimum_similarity <= 1:
    raise ValueError("Similarity must be between 0 and 1")

if requested:
    libraries = requested
else:
    libraries = source.get_libraries(profile)

if not libraries:
    raise ValueError("No enabled libraries selected")

for library in libraries:

    if library not in config["libraries"]:
        raise ValueError(
            f"Unknown library: {library}"
        )

    settings = config["libraries"][library]

    if not settings["enabled"]:
        raise ValueError(
            f"Library disabled: {library}"
        )

    if not settings["validated"]:
        raise ValueError(
            f"Library not validated: {library}"
        )

config["profiles"]["runtime_selected"] = libraries

config["default_profile"] = "runtime_selected"

config["screening"]["mode"] = "both"

config["screening"]["minimum_analog_similarity"] = (
    minimum_similarity
)

config["screening"]["maximum_analogs_per_library"] = (
    maximum_analogs
)

config["api"]["allow_remote_submission"] = allow_remote

config["api"]["maximum_new_jobs_per_run"] = (
    maximum_new_jobs
)

config["api"]["parallel_jobs"] = 1

config["cache"]["database"] = str(
    source.database_path()
)

config["cache"]["enabled"] = True
config["cache"]["reuse_completed"] = True
config["cache"]["refresh"] = False

# Do not resume remote sessions during offline execution.

config["cache"]["resume_sessions"] = allow_remote

runtime.write_text(
    json.dumps(config, indent=4) + "\n"
)

# Validate the newly created configuration.

CommercialConfig(runtime)

assert base.read_bytes() == original_bytes

print("\nSelected commercial libraries:")

for library in libraries:
    print(" ", library)

print("\nNumber of libraries:", len(libraries))
print("Minimum similarity:", minimum_similarity)
print("Maximum analogs:", maximum_analogs)
print("Remote submission:", allow_remote)
print("Maximum new jobs:", maximum_new_jobs)
print("Runtime configuration:", runtime)
print("Production configuration preserved")

PY

fi

# --------------------------------------------------
# Stages 1-7: Existing V6 pipeline
# --------------------------------------------------

if [[ "$MODE" == "full" ]]; then

    stage "STAGES 1-7: MOLECULAR GENERATION AND OPTIMIZATION"

    bash scripts/run_reinvent_v6.sh \
        "$INPUT" \
        "$N" \
        "$OUT" \
        "$PROFILE" \
        "$SELECT"

    SELECTED="$OUT/cmos/CMOS_V2_selected.csv"

else

    stage "STAGES 1-7: USING EXISTING CMOS V2 LIBRARY"

    echo "Skipping molecular generation"
    echo "Skipping ADMET prediction"
    echo "Skipping CMOS optimization"

fi

[[ -s "$SELECTED" ]] || \
    die "Selected molecular library missing: $SELECTED"

echo "Selected molecular library: $SELECTED"

# --------------------------------------------------
# Optional commercial screening
# --------------------------------------------------

if [[ "$SKIP_COMMERCIAL" == true ]]; then

    stage "COMMERCIAL SCREENING SKIPPED"

    echo "Final molecular library: $SELECTED"

    stage "REINVENT V7 PIPELINE COMPLETE"

    exit 0

fi

mkdir -p "$COMMERCIAL" "$FINAL"

# --------------------------------------------------
# Stage 8: Commercial screening
# --------------------------------------------------

stage "STAGE 8: COMMERCIAL SCREENING"

SCREENING_ARGS=(

    --input "$SELECTED"

    --smiles-column SMILES

    --config "$COMMERCIAL_RUNTIME_CONFIG"

    --profile runtime_selected

    --mode "$SCREENING_MODE"

    --minimum-similarity "$MINIMUM_SIMILARITY"

    --maximum-analogs "$MAXIMUM_ANALOGS"

    --output "$COMMERCIAL"

)

if [[ "$ALLOW_REMOTE" == true ]]; then

    SCREENING_ARGS+=(

        --allow-remote

        --maximum-new-jobs "$MAXIMUM_NEW_JOBS"

    )

fi

"$PYTHON" -u \
    scripts/commercial/commercial_screening.py \
    "${SCREENING_ARGS[@]}"

test -s "$COMMERCIAL/Commercial_Availability.csv"
test -s "$COMMERCIAL/screening_summary.json"

# --------------------------------------------------
# Stage 9: CMOS integration
# --------------------------------------------------

stage "STAGE 9: CMOS COMMERCIAL INTEGRATION"

"$PYTHON" -u \
    scripts/commercial/commercial_integration.py \
    --input "$SELECTED" \
    --screening-dir "$COMMERCIAL" \
    --output "$FINAL"

test -s "$FINAL/CMOS_V2_Commercial_Annotated.csv"
test -s "$FINAL/Commercial_Screening_Summary.json"

# --------------------------------------------------
# Final validation
# --------------------------------------------------

stage "FINAL V7 DATA INTEGRITY VALIDATION"

"$PYTHON" - \
    "$SELECTED" \
    "$FINAL" \
    "$COMMERCIAL_RUNTIME_CONFIG" \
    "$BASE_CONFIG" \
    "$ALLOW_REMOTE" \
    "$MAXIMUM_NEW_JOBS" <<'PY'

import sys
import json
from pathlib import Path

import pandas as pd

source = Path(sys.argv[1])
final = Path(sys.argv[2])

runtime = json.loads(
    Path(sys.argv[3]).read_text()
)

base = json.loads(
    Path(sys.argv[4]).read_text()
)

allow_remote = sys.argv[5] == "true"
maximum_new_jobs = int(sys.argv[6])

original = pd.read_csv(
    source,
    dtype=str,
    keep_default_na=False,
)

annotated = pd.read_csv(
    final / "CMOS_V2_Commercial_Annotated.csv",
    dtype=str,
    keep_default_na=False,
)

pd.testing.assert_frame_equal(
    annotated[original.columns],
    original,
)

assert len(original) == len(annotated)

assert len(annotated.columns) == len(original.columns) + 21

assert annotated[
    "Commercial_Stock_Status"
].eq("UNVERIFIED").all()

summary = json.loads(
    (
        final / "Commercial_Screening_Summary.json"
    ).read_text()
)

screening = summary["Screening_Summary"]

expected_libraries = runtime["profiles"]["runtime_selected"]

assert set(screening["Libraries"]) == set(expected_libraries)

assert screening["Mode"] == "both"

assert screening["New_Submission_Attempts"] <= maximum_new_jobs

assert base["api"]["allow_remote_submission"] is False

if not allow_remote:

    assert screening["New_Submission_Attempts"] == 0
    assert screening["Resumed_Jobs"] == 0

print("Original molecules:", len(original))
print("Annotated molecules:", len(annotated))
print("Original columns:", len(original.columns))
print("Final columns:", len(annotated.columns))

print("\nCommercial libraries:", expected_libraries)

print(
    "New submissions:",
    screening["New_Submission_Attempts"],
)

print("Resumed sessions:", screening["Resumed_Jobs"])
print("Cache hits:", screening["Cache_Hits"])

print("\nScreening status:")

print(
    annotated["Commercial_Screening_Status"]
    .value_counts()
    .to_string()
)

print("\nScreening coverage:")

print(
    annotated["Commercial_Coverage_Status"]
    .value_counts()
    .to_string()
)

problem_statuses = {
    "SUBMISSION_UNCERTAIN",
    "SCREENING_ERROR",
    "REVIEW_REQUIRED",
    "CLASSIFICATION_ERROR",
    "RESUME_DISABLED",
}

problems = {
    status: count
    for status, count in screening["Job_Status_Counts"].items()
    if status in problem_statuses
}

if problems:
    print("\nWARNING: Screening jobs require review:")
    print(json.dumps(problems, indent=4))

print("\nALL V7 DATA INTEGRITY TESTS PASSED")

PY

stage "REINVENT V7 PIPELINE COMPLETE"

echo "Results: $OUT"
echo "Runtime configuration: $COMMERCIAL_RUNTIME_CONFIG"
echo "Commercial screening: $COMMERCIAL"
echo "Final molecular library: $FINAL/CMOS_V2_Commercial_Annotated.csv"
echo "Exact candidates: $FINAL/CMOS_V2_Commercial_Exact.csv"
echo "Strict exact candidates: $FINAL/CMOS_V2_Commercial_Strict_Exact.csv"
echo "Analog candidates: $COMMERCIAL/commercial_analogs.csv"
echo "Per-parent summary: $FINAL/Commercial_Per_Parent_Summary.csv"
echo "Final summary: $FINAL/Commercial_Screening_Summary.json"

