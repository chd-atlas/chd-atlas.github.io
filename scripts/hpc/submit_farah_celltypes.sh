#!/bin/bash
# Submit wrapper for farah_celltypes.sbatch.
#
# Sets --mail-user from $SLURM_EMAIL_NOTIFICATION (never a literal address) and
# forwards job vars through the environment. Override anything inline:
#   OUT_DIR=/scratch/me/farah bash submit_farah_celltypes.sh
# Run it from a directory that holds this script, farah_celltypes.sbatch,
# farah_extract.py and the registry as genes.tsv (or set REGISTRY).
set -euo pipefail
cd "$(dirname "$0")"

CPU_PARTITION="${CPU_PARTITION:-all_cpu.p}"
REGISTRY="${REGISTRY:-$PWD/genes.tsv}"
OUT_DIR="${OUT_DIR:?Set OUT_DIR to a writable output directory}"

if [[ -z "${SLURM_EMAIL_NOTIFICATION:-}" && -r "$HOME/.bashrc" ]]; then
  SLURM_EMAIL_NOTIFICATION="$(
    sed -n "s/^[[:space:]]*\\(export[[:space:]]\\+\\)\\?SLURM_EMAIL_NOTIFICATION=[\"']\\?\\([^\"'#[:space:]]\\+\\).*/\\2/p" \
      "$HOME/.bashrc" | tail -n1
  )"
fi
export SLURM_EMAIL_NOTIFICATION="${SLURM_EMAIL_NOTIFICATION:-}"
if [[ -n "$SLURM_EMAIL_NOTIFICATION" ]]; then
  echo 'mail: ON  (address from $SLURM_EMAIL_NOTIFICATION; not printed)'
else
  echo 'mail: OFF (export SLURM_EMAIL_NOTIFICATION=you@example.org to enable)'
fi

export REGISTRY OUT_DIR
for v in INPUT_H5AD HVANTK_PYTHON HVANTK_REPO MIN_CELLS; do
  [[ -n "${!v:-}" ]] && export "$v"
done

SBATCH_ARGS=(--partition="$CPU_PARTITION" --export=ALL)
[[ -n "$SLURM_EMAIL_NOTIFICATION" ]] && SBATCH_ARGS+=(--mail-user="$SLURM_EMAIL_NOTIFICATION")

echo "partition: $CPU_PARTITION"
echo "registry:  $REGISTRY"
echo "out dir:   $OUT_DIR"
# Build a real array: expanding an unset one as "${ARR[@]:-}" passes an empty
# filename to sbatch.
EXTRA=()
if [[ -n "${EXTRA_SBATCH_ARGS:-}" ]]; then
  # shellcheck disable=SC2206  # deliberate word-splitting: caller passes a flag string
  EXTRA=($EXTRA_SBATCH_ARGS)
fi

sbatch "${SBATCH_ARGS[@]}" ${EXTRA[@]+"${EXTRA[@]}"} farah_celltypes.sbatch
