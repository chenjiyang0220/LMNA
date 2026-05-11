#!/usr/bin/env bash

set -euo pipefail

BW_DIR="/media/ggj/ggj/CJY/cuttag/scATAC/bw"
OUT_DIR="/media/ggj/ggj/CJY/cuttag/scATAC/TE/TE_pantissue_up"
PROJ_RDS="/media/ggj/ggj/CJY/cuttag/scATAC/proj_callpeak.rds"
TE_IN_PEAKSET_BED="${OUT_DIR}/mm10.te.peakolp.bed"

WT_BW="${BW_DIR}/WT.10kcell.bw"
KO_BW="${BW_DIR}/KO.10kcell.bw"

mkdir -p "${OUT_DIR}"

zcat "${TE_TXT_GZ}" \
    | awk 'BEGIN{OFS="\t"} NR>1 {print $6,$7,$8,$10,$11,$12}' \
    > "${TE_BED}"

RAW_COUNTS="${OUT_DIR}/panTissue_TE_ATAC.rawCounts.txt"
NPZ_FILE="${OUT_DIR}/panTissue_TE_ATAC.npz"

multiBigwigSummary BED-file \
  --bwfiles "${WT_BW}" "${KO_BW}" \
  --BED "${TE_IN_PEAKSET_BED}" \
  --outFileName "${NPZ_FILE}" \
  --outRawCounts "${RAW_COUNTS}" \
  -p 8

echo "Done. Raw counts (per TE in peakset, WT/KO) in: ${RAW_COUNTS}"
