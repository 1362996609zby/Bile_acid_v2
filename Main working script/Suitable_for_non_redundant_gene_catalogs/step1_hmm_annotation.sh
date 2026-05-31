#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Step 1: Annotate a non-redundant protein catalog with custom HMM profiles
#
# This script runs HMMER `hmmscan` against each HMM profile in a directory and
# writes one domain table per HMM family. It is designed for reusable pipelines:
#   - no project-specific hard-coded paths
#   - configurable E-value thresholds
#   - automatic resume: finished HMMs are skipped
#   - temporary output files are used to avoid incomplete results after crashes
#   - optional automatic hmmpress for hmmscan-compatible HMM databases
#
# Typical usage:
#   bash step1_hmm_annotation.sh \
#     --protein-fasta /path/to/Unigenes.protein.fa \
#     --hmm-dir /path/to/hmm_mod \
#     --out-dir /path/to/hmm_results \
#     --cpu 40 \
#     --evalue 1e-3 \
#     --dom-evalue 1e-3
#
# Outputs:
#   <out-dir>/<hmm_name>_hmmscan.txt          # HMMER --domtblout
#   <out-dir>/<hmm_name>_hmmscan.stdout.txt   # full hmmscan stdout
#   <out-dir>/<hmm_name>.done                 # completion flag for resuming
###############################################################################

usage() {
  cat <<'EOF'
Usage:
  bash step1_hmm_annotation.sh \
    --protein-fasta FILE \
    --hmm-dir DIR \
    --out-dir DIR \
    [--cpu INT] \
    [--evalue FLOAT] \
    [--dom-evalue FLOAT] \
    [--hmm-pattern "*.hmm"] \
    [--force] \
    [--no-auto-press]

Required arguments:
  --protein-fasta FILE   Protein FASTA file to annotate.
  --hmm-dir DIR          Directory containing HMM profile files.
  --out-dir DIR          Output directory for hmmscan results.

Optional arguments:
  --cpu INT              Number of CPU threads for hmmscan. Default: 8.
  --evalue FLOAT         Sequence-level E-value cutoff (-E). Default: 1e-3.
  --dom-evalue FLOAT     Domain-level E-value cutoff (--domE). Default: 1e-3.
  --hmm-pattern PATTERN  File pattern for HMM profiles. Default: "*.hmm".
  --force                Re-run all HMMs even if output/done files already exist.
  --no-auto-press        Do not run hmmpress automatically when .h3* files are missing.
  -h, --help             Show this help message.

Notes:
  - For hmmscan, HMM files usually need to be pressed with hmmpress.
  - The script skips a family if <hmm_name>.done exists and the domtblout file is non-empty.
  - Output is first written to *.tmp and only renamed after hmmscan succeeds.
EOF
}

PROTEIN_FASTA=""
HMM_DIR=""
OUTPUT_DIR=""
CPU=8
EVALUE="1e-3"
DOMEVALUE="1e-3"
HMM_PATTERN="*.hmm"
FORCE=0
AUTO_PRESS=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --protein-fasta) PROTEIN_FASTA="$2"; shift 2 ;;
    --hmm-dir) HMM_DIR="$2"; shift 2 ;;
    --out-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --cpu) CPU="$2"; shift 2 ;;
    --evalue) EVALUE="$2"; shift 2 ;;
    --dom-evalue) DOMEVALUE="$2"; shift 2 ;;
    --hmm-pattern) HMM_PATTERN="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    --no-auto-press) AUTO_PRESS=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "${PROTEIN_FASTA}" || -z "${HMM_DIR}" || -z "${OUTPUT_DIR}" ]]; then
  echo "ERROR: --protein-fasta, --hmm-dir and --out-dir are required." >&2
  usage
  exit 1
fi

if [[ ! -s "${PROTEIN_FASTA}" ]]; then
  echo "ERROR: Protein FASTA not found or empty: ${PROTEIN_FASTA}" >&2
  exit 1
fi

if [[ ! -d "${HMM_DIR}" ]]; then
  echo "ERROR: HMM directory not found: ${HMM_DIR}" >&2
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"

command -v hmmscan >/dev/null 2>&1 || {
  echo "ERROR: hmmscan not found in PATH. Please activate/install HMMER." >&2
  exit 1
}

if [[ "${AUTO_PRESS}" -eq 1 ]]; then
  command -v hmmpress >/dev/null 2>&1 || {
    echo "ERROR: hmmpress not found in PATH, but auto-press is enabled." >&2
    echo "       Install HMMER or rerun with --no-auto-press if HMM files are already pressed." >&2
    exit 1
  }
fi

shopt -s nullglob
HMM_FILES=( "${HMM_DIR}"/${HMM_PATTERN} )
shopt -u nullglob

if [[ ${#HMM_FILES[@]} -eq 0 ]]; then
  echo "ERROR: No HMM files found in ${HMM_DIR} with pattern ${HMM_PATTERN}" >&2
  exit 1
fi

echo "============================================================"
echo "Step 1: HMM annotation with hmmscan"
echo "Protein FASTA : ${PROTEIN_FASTA}"
echo "HMM directory : ${HMM_DIR}"
echo "Output dir    : ${OUTPUT_DIR}"
echo "HMM files     : ${#HMM_FILES[@]}"
echo "CPU           : ${CPU}"
echo "E-value       : ${EVALUE}"
echo "Domain E-value: ${DOMEVALUE}"
echo "Force rerun   : ${FORCE}"
echo "Auto hmmpress : ${AUTO_PRESS}"
echo "============================================================"

for hmmfile in "${HMM_FILES[@]}"; do
  hmmname=$(basename "${hmmfile}" .hmm)

  domtbl="${OUTPUT_DIR}/${hmmname}_hmmscan.txt"
  stdout_log="${OUTPUT_DIR}/${hmmname}_hmmscan.stdout.txt"
  done_flag="${OUTPUT_DIR}/${hmmname}.done"

  domtbl_tmp="${domtbl}.tmp"
  stdout_tmp="${stdout_log}.tmp"

  if [[ "${FORCE}" -eq 0 && -f "${done_flag}" && -s "${domtbl}" ]]; then
    echo "SKIP: ${hmmname} (done flag and output found)"
    continue
  fi

  if [[ "${FORCE}" -eq 0 && -s "${domtbl}" ]]; then
    echo "SKIP: ${hmmname} (existing output found; creating done flag)"
    touch "${done_flag}"
    continue
  fi

  if [[ "${FORCE}" -eq 1 ]]; then
    rm -f "${domtbl}" "${stdout_log}" "${done_flag}"
  fi

  if [[ "${AUTO_PRESS}" -eq 1 ]]; then
    if [[ ! -s "${hmmfile}.h3f" || ! -s "${hmmfile}.h3i" || ! -s "${hmmfile}.h3m" || ! -s "${hmmfile}.h3p" ]]; then
      echo "PRESS: ${hmmname} (missing .h3* files)"
      hmmpress -f "${hmmfile}" >/dev/null
    fi
  fi

  echo "RUN : ${hmmname} (E=${EVALUE}, domE=${DOMEVALUE})"

  rm -f "${domtbl_tmp}" "${stdout_tmp}"

  hmmscan --cpu "${CPU}" \
    -E "${EVALUE}" --domE "${DOMEVALUE}" \
    --domtblout "${domtbl_tmp}" \
    "${hmmfile}" "${PROTEIN_FASTA}" > "${stdout_tmp}"

  if [[ ! -s "${domtbl_tmp}" ]]; then
    echo "ERROR: Empty domtblout for ${hmmname}. Not marking as done." >&2
    exit 1
  fi

  mv -f "${domtbl_tmp}" "${domtbl}"
  mv -f "${stdout_tmp}" "${stdout_log}"
  touch "${done_flag}"

  echo "DONE: ${hmmname}"
done

echo "============================================================"
echo "Step 1 completed successfully."
echo "Results saved in: ${OUTPUT_DIR}"
echo "============================================================"
