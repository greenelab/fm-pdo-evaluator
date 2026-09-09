#!/usr/bin/env bash
# Submit rung 0's three stages as a dependency chain: assign -> slice array -> combine.
#
# Run from the repository root on the LOCAL machine; every submission goes through ralpine so
# the boundary it enforces holds (PROCESS section 2). Each job id is read from sbatch's own
# output and passed as --dependency=afterok to the next, then the dependency is read back with
# `ralpine jobinfo` -- a chain that silently drops its dependency runs immediately and out of
# order (PROCESS section 2, "Chained jobs").
#
#   scripts/alpine/submit_rung0_chain.sh            # all three stages
#   scripts/alpine/submit_rung0_chain.sh --from slice   # the array and the combine only
#   scripts/alpine/submit_rung0_chain.sh --from combine # the combine only
set -euo pipefail

RALPINE="$(dirname "${BASH_SOURCE[0]}")/ralpine"
FROM="assign"
if [[ "${1:-}" == "--from" ]]; then
  FROM="${2:?--from needs a stage: assign, slice or combine}"
fi

job_id() { grep -oE '[0-9]+$' <<<"$1" | tail -1; }
check_dependency() {
  local id="$1" want="$2"
  local dep
  dep="$("$RALPINE" jobinfo "$id" | grep -oE 'Dependency=[^ ]+' || true)"
  if [[ "$dep" != *"$want"* ]]; then
    echo "job $id was submitted without its dependency ($dep, wanted $want); cancel it" >&2
    exit 1
  fi
  echo "  job $id: $dep"
}

DEP=""
if [[ "$FROM" == "assign" ]]; then
  out="$("$RALPINE" submit scripts/alpine/rung0_assign.sbatch)"
  ASSIGN="$(job_id "$out")"
  echo "assign:  job $ASSIGN"
  DEP="--dependency=afterok:$ASSIGN"
fi
if [[ "$FROM" == "assign" || "$FROM" == "slice" ]]; then
  out="$("$RALPINE" submit scripts/alpine/rung0_slice.sbatch ${DEP:+"$DEP"})"
  SLICE="$(job_id "$out")"
  echo "slices:  array job $SLICE"
  [[ -n "$DEP" ]] && check_dependency "$SLICE" "afterok:$ASSIGN"
  DEP="--dependency=afterok:$SLICE"
fi
out="$("$RALPINE" submit scripts/alpine/rung0_combine.sbatch ${DEP:+"$DEP"})"
COMBINE="$(job_id "$out")"
echo "combine: job $COMBINE"
[[ -n "$DEP" ]] && check_dependency "$COMBINE" "afterok:"
echo "watch with: scripts/alpine/ralpine sq"
