# ----- plot Heatmap & Profile example -----

# ---- H3K4me3up promoter heatmaps ----
BASE_DIR="./scATAC/copair/only_epi"
OUT_DIR="${BASE_DIR}/histoned_EP_tables"
WT_K4me3="./bulkCUT/bw/WT.H3K4me3.merge.bw"
KO_K4me3="./bulkCUT/bw/KO.H3K4me3.merge.bw"

WT_BED="${OUT_DIR}/WT_unique_promoter_H3K4me3up.bed"
KO_BED="${OUT_DIR}/KO_unique_promoter_H3K4me3up.bed"

# compare WT and KO signal on WT-up regions
computeMatrix reference-point \
  -S "${WT_K4me3}" "${KO_K4me3}" \
  -R "${WT_BED}" \
  --referencePoint center \
  --beforeRegionStartLength 5000 \
  --afterRegionStartLength 5000 \
  --binSize 10 \
  --sortRegions descend \
  --sortUsing mean \
  --missingDataAsZero \
  -o "${OUT_DIR}/H3K4me3_wtup_matrix.mat.gz"

plotProfile \
  -m "${OUT_DIR}/H3K4me3_wtup_matrix.mat.gz" \
  -out "${OUT_DIR}/H3K4me3_wtup_profile.pdf" \
  --perGroup \
  --colors blue red \
  --plotTitle "H3K4me3 signal around WT promoter peaks" \
  --refPointLabel "Peak Center" \
  --regionsLabel "ATAC Peaks" \
  --samplesLabel "WT" "KO" \
  --yAxisLabel "Normalized read density"

plotHeatmap \
  -m "${OUT_DIR}/H3K4me3_wtup_matrix.mat.gz" \
  -out "${OUT_DIR}/H3K4me3_wtup_heatmap.pdf" \
  --colorMap GnBu \
  --whatToShow 'heatmap and colorbar' \
  --refPointLabel "Peak Center" \
  --regionsLabel "ATAC Peaks" \
  --samplesLabel "WT" "KO" \
  --boxAroundHeatmaps no \
  --legendLocation upper-left
