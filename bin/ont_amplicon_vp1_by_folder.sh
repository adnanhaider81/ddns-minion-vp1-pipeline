#!/usr/bin/env bash
set -euo pipefail

die(){ echo "ERROR: $*" >&2; exit 1; }

usage(){
  cat <<'USAGE'
Standalone ONT VP1 pipeline for DDNS stool-culture sequencing

This wrapper always runs the bundled vp1 pipeline with:
  - references/references.all.unique.fasta
  - resources/BarcodedPrimers.xlsx
  - piranha-vp1 behavior only

Required:
  --barcode-map /path/to/barcodes.csv
  --fastq-pass /path/to/fastq_pass
  --out /path/to/output_dir

Common optional overrides:
  --threads N
  --run-name NAME
  --run-report-html PATH
  --min-len N
  --max-len N
  --min-cov N

Any additional options supported by the internal vp1 engine may also be passed,
except --refs, --primers-xlsx, --preset, and --piranha-vp1-mode, which are
fixed by this wrapper.
USAGE
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ENGINE="$SCRIPT_DIR/vp1_pipeline_internal.sh"
REFS="$REPO_DIR/references/references.all.unique.fasta"
PRIMERS="$REPO_DIR/resources/BarcodedPrimers.xlsx"

[[ -f "$ENGINE" ]] || die "Missing internal engine: $ENGINE"
[[ -f "$REFS" ]] || die "Missing bundled reference FASTA: $REFS"
[[ -f "$PRIMERS" ]] || die "Missing bundled primers workbook: $PRIMERS"

if [[ $# -eq 0 ]]; then
  usage
  exit 1
fi

for arg in "$@"; do
  case "$arg" in
    -h|--help)
      usage
      exit 0
      ;;
    --refs|--primers-xlsx|--preset|--piranha-vp1-mode)
      die "Do not pass $arg to this standalone vp1 wrapper; it is fixed internally."
      ;;
  esac
done

bash "$ENGINE" \
  --piranha-vp1-mode \
  --refs "$REFS" \
  --primers-xlsx "$PRIMERS" \
  "$@"
