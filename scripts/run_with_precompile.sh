#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./scripts/run_with_precompile.sh <job> [job args...]

Examples:
  ./scripts/run_with_precompile.sh src/code/jobs/run_reversal_transition_job.jl --config config/default.toml
  ./scripts/run_with_precompile.sh run_reversal_transition_job --config config/default.toml
  ./scripts/run_with_precompile.sh reversal_transition --config config/default.toml

Notes:
  - Always runs precompilation first (scripts/precompile.jl).
  - Forwards all extra arguments to the job (e.g. --config, --set key=value).
EOF
}

resolve_job_path() {
  local input="$1"
  local candidate
  local -a candidates=(
    "$input"
    "src/code/jobs/$input"
    "src/code/jobs/${input}.jl"
    "src/code/jobs/run_${input}_job.jl"
  )

  for candidate in "${candidates[@]}"; do
    if [[ -f "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  return 1
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || "$#" -lt 1 ]]; then
  usage
  exit 0
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
cd "$ROOT_DIR"

job_input="$1"
shift

if ! job_path="$(resolve_job_path "$job_input")"; then
  echo "Error: job not found for input '$job_input'." >&2
  echo "Try one of: src/code/jobs/<file>.jl, <file>, <task>, or run_<task>_job." >&2
  exit 1
fi

julia_exe="${JULIA_EXE:-julia}"
common_flags=(--startup-file=no --project=.)

echo "[1/2] Precompile environment"
"$julia_exe" "${common_flags[@]}" scripts/precompile.jl

echo "[2/2] Run job: $job_path"
"$julia_exe" "${common_flags[@]}" "$job_path" "$@"

