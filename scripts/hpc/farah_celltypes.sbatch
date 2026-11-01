#!/bin/bash
#SBATCH --job-name=farah-celltypes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=0-04:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --mail-type=END,FAIL
# --mail-user is deliberately NOT a #SBATCH line: Slurm reads directives
# literally and does not expand env vars, so a `--mail-user=$VAR` directive
# silently mails nobody. submit_farah_celltypes.sh passes it on the CLI.
#
# Pseudobulks the Farah et al. 2024 heart atlas with hvantk, twice (by Age, Rep,
# cell class and Region; then without Region), and flattens both into
# $OUT_DIR/pseudobulk.tsv for scripts/convert_farah2024.py.
# The memory and time are NOT yet measured on this input: the first real run
# should be read with `seff` and these numbers corrected from it.

set -euo pipefail
PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_DIR"

# --- Site (overridable; no absolute path is committed, $WORK and $HOME expand per user) ---
HPC_ENV_MODULE="${HPC_ENV_MODULE:-hpc-env/13.1}"
HVANTK_PYTHON="${HVANTK_PYTHON:-$HOME/projects/pyvatk/.venv/bin/python}"
HVANTK_REPO="${HVANTK_REPO:-$HOME/projects/pyvatk}"
# tables.MIN_CELLS_PER_HEART: the site says a heart counts only with 10 or more
# cells, and convert_farah2024.py refuses a group under it, so a lower value here
# fails there.
MIN_CELLS="${MIN_CELLS:-10}"
INPUT_H5AD="${INPUT_H5AD:-${WORK:?WORK is unset; pass INPUT_H5AD explicitly}/hvantk-datasets/expression/cardiac/farah_2024.h5ad}"

# --- Required ---
REGISTRY="${REGISTRY:?Set REGISTRY to the genes.tsv shipped with this job}"
OUT_DIR="${OUT_DIR:?Set OUT_DIR to a writable output directory}"

BY_REGION="$OUT_DIR/farah_by_region.h5ad"
POOLED="$OUT_DIR/farah_pooled.h5ad"
PSEUDOBULK="$OUT_DIR/pseudobulk.tsv"

# --- Prove the inputs and the output path before the expensive passes ---
for f in "$INPUT_H5AD" "$REGISTRY" "$PROJECT_DIR/farah_extract.py" ; do
  [[ -r "$f" ]] || { echo "ERROR: not readable: $f" >&2; exit 1; }
done
[[ -x "$HVANTK_PYTHON" ]] || { echo "ERROR: HVANTK_PYTHON not executable: $HVANTK_PYTHON" >&2; exit 1; }
mkdir -p "$OUT_DIR"
: > "$OUT_DIR/.writetest" 2>/dev/null || { echo "ERROR: $OUT_DIR not writable" >&2; exit 1; }
rm -f "$OUT_DIR/.writetest"

# --- Modules: only the site environment; python comes from the hvantk venv ---
if command -v ml >/dev/null 2>&1; then
  ml purge || true
  ml "$HPC_ENV_MODULE" || echo "WARNING: could not load $HPC_ENV_MODULE" >&2
fi
# A cluster module must not shadow the venv's own packages.
unset PYTHONPATH

# Captured rather than echoed inline: a command substitution inside `echo`
# fails silently under `set -e` and would print an empty commit. The output's
# provenance needs this line, so an empty one stops the job.
HVANTK_HEAD="$(git -C "$HVANTK_REPO" rev-parse HEAD)"
: "${HVANTK_HEAD:?could not read the hvantk commit from $HVANTK_REPO}"

echo "== farah celltypes =="
echo "date:         $(date -Iseconds)"
echo "host:         $(hostname)"
echo "job:          ${SLURM_JOB_ID:-<none>}"
echo "python:       $("$HVANTK_PYTHON" -V 2>&1)"
echo "hvantk HEAD:  $HVANTK_HEAD"
echo "min cells:    $MIN_CELLS"
echo "input:        $INPUT_H5AD"
echo "out:          $OUT_DIR"
echo "input sha256: $(sha256sum "$INPUT_H5AD" | cut -d' ' -f1)"

"$HVANTK_PYTHON" -m hvantk expression summarize -m "$INPUT_H5AD" \
  --group-by Age --group-by Rep --group-by major_cell_class --group-by Region \
  --min-cells "$MIN_CELLS" -o "$BY_REGION"

"$HVANTK_PYTHON" -m hvantk expression summarize -m "$INPUT_H5AD" \
  --group-by Age --group-by Rep --group-by major_cell_class \
  --min-cells "$MIN_CELLS" -o "$POOLED"

"$HVANTK_PYTHON" "$PROJECT_DIR/farah_extract.py" \
  --by-region "$BY_REGION" --pooled "$POOLED" \
  --registry "$REGISTRY" --out "$PSEUDOBULK"

echo "-- output --"
wc -l "$PSEUDOBULK"
sha256sum "$PSEUDOBULK"
echo "done: $(date -Iseconds)"
