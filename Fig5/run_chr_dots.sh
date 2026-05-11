#!/usr/bin/env bash
set -uo pipefail

COOL_PATH="/media/ggj/ggj/CJY/nature_WXY/Hi-c/file/KO/KO_merged_noXY_1kb_10bin.cool"
BED_DIR="/media/ggj/ggj/CJY/nature_WXY/Hi-c/file/RE_sites"
OUT_DIR="/media/ggj/ggj/CJY/nature_WXY/Hi-c/analysis/dots/KO_10bin"
WEIGHT_COL="weight"   

MASTER_LOG="${OUT_DIR}/KO_noXY_10kb_all_chr_dots_report.txt"
echo "=== $(date) === Starting all chromosomes dots calling" | tee "${MASTER_LOG}"
echo "COOL: ${COOL_PATH}" | tee -a "${MASTER_LOG}"
echo "BED_DIR: ${BED_DIR}" | tee -a "${MASTER_LOG}"
echo "OUT_DIR: ${OUT_DIR}" | tee -a "${MASTER_LOG}"
echo "" | tee -a "${MASTER_LOG}"

SUCCESS_CHR=""
FAILED_CHR=""

# loop
for i in $(seq 1 19); do
    CHR="chr${i}"
    VIEW_BED="${BED_DIR}/mm10_${CHR}.bed"

    if [[ ! -f "${VIEW_BED}" ]]; then
        echo "[WARNING] ${VIEW_BED} not found, skipping ${CHR}..." | tee -a "${MASTER_LOG}"
        FAILED_CHR="${FAILED_CHR} ${CHR}(no_bed)"
        continue
    fi

    EXPECTED_TSV="${OUT_DIR}/KO_noXY_10kb_${CHR}_expected.tsv"
    DOTS_BEDPE="${OUT_DIR}/KO_noXY_10kb_${CHR}_dots.bedpe"
    REPORT_TXT="${OUT_DIR}/KO_noXY_10kb_${CHR}_dots_report.txt"
    
    echo "========================================" | tee -a "${MASTER_LOG}"
    echo "[${CHR}] Starting at $(date)" | tee "${REPORT_TXT}" | tee -a "${MASTER_LOG}"
    echo "VIEW: ${VIEW_BED}" | tee -a "${REPORT_TXT}"
    echo "" | tee -a "${REPORT_TXT}"
    
    # Step 1: expected-cis
    echo "[1] Running expected-cis on ${CHR} ..." | tee -a "${REPORT_TXT}"
    if ! cooltools expected-cis \
      --view "${VIEW_BED}" \
      --clr-weight-name "${WEIGHT_COL}" \
      -p 4 \
      "${COOL_PATH}" \
      -o "${EXPECTED_TSV}" 2>&1 | tee -a "${REPORT_TXT}"; then
        echo "[ERROR] expected-cis failed for ${CHR}, skipping..." | tee -a "${REPORT_TXT}" | tee -a "${MASTER_LOG}"
        FAILED_CHR="${FAILED_CHR} ${CHR}(expected-cis)"
        continue
    fi
    
    echo "" | tee -a "${REPORT_TXT}"
    
    # Step 2: dots
    echo "[2] Running dots on ${CHR} ..." | tee -a "${REPORT_TXT}"
    if ! cooltools dots \
      --view "${VIEW_BED}" \
      --clr-weight-name "${WEIGHT_COL}" \
      --max-loci-separation 2000000 \
      -p 4 \
      "${COOL_PATH}" \
      "${EXPECTED_TSV}::balanced.avg" \
      -o "${DOTS_BEDPE}" 2>&1 | tee -a "${REPORT_TXT}"; then
        echo "[ERROR] dots failed for ${CHR}, continuing to next chromosome..." | tee -a "${REPORT_TXT}" | tee -a "${MASTER_LOG}"
        FAILED_CHR="${FAILED_CHR} ${CHR}(dots)"
        continue
    fi
    
    echo "" | tee -a "${REPORT_TXT}"
    echo "[DONE] ${CHR} dots calling finished successfully." | tee -a "${REPORT_TXT}" | tee -a "${MASTER_LOG}"
    echo "Expected file: ${EXPECTED_TSV}" | tee -a "${REPORT_TXT}"
    echo "Dots file    : ${DOTS_BEDPE}" | tee -a "${REPORT_TXT}"
    SUCCESS_CHR="${SUCCESS_CHR} ${CHR}"
done

echo "" | tee -a "${MASTER_LOG}"
echo "========================================" | tee -a "${MASTER_LOG}"
echo "=== $(date) === All chromosomes processing completed!" | tee -a "${MASTER_LOG}"
echo "" | tee -a "${MASTER_LOG}"
echo "SUCCESS:${SUCCESS_CHR:-" (none)"}" | tee -a "${MASTER_LOG}"
echo "FAILED:${FAILED_CHR:-" (none)"}" | tee -a "${MASTER_LOG}"