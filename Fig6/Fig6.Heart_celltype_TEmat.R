library(ArchR)
library(GenomicRanges)
library(rtracklayer)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)
library(ComplexHeatmap)
library(circlize)

# Import RepeatMasker annotation (UCSC or custom)
# Option A: from UCSC table download
rmsk <- read.table("/path/scATAC/TE/Mouse_TE_annotation/TE_LTR_LINE_SINE.bed", header = FALSE, sep = "\t")
colnames(rmsk) <- c("chr", "start", "end","repClass")

# Convert to GRanges
te_gr <- GRanges(
  seqnames = rmsk$chr,
  ranges   = IRanges(start = rmsk$start, end = rmsk$end),
  repClass = rmsk$repClass
)

# Optional: further split by family
# e.g., SINE/Alu, SINE/B1, ERV1, ERVK, L1, etc.
setwd("/path/scATAC")
proj <- readRDS("/path/scATAC/tissue/Heart/proj_anno.rds")

# Heart single-cell TE accessibility analysis from fragments
# This workflow does not add a matrix back into ArchRProj.
# Instead, it extracts heart-project fragments, counts overlaps on merged
# non-overlapping TE regions, normalizes the TE-region x cell matrix, and
# summarizes TE accessibility bias by cell type and genotype.

out_dir <- "/path/scATAC/tissue/Heart/TE_matrix"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Create a non-overlapping TE feature set to avoid double-counting nearby TEs
te_by_class <- split(te_gr, te_gr$repClass)
te_merged <- lapply(te_by_class, reduce)

# Combine merged regions and annotate class + unique region id
te_features <- do.call(c, lapply(names(te_merged), function(cls) {
  gr <- te_merged[[cls]]
  gr$repClass <- cls
  gr$region_id <- paste0(cls, "_", seq_along(gr))
  gr
}))
names(te_features) <- te_features$region_id
te_classes <- sort(unique(as.character(te_features$repClass)))

# Cell metadata from the heart single-cell ArchR project
cell_md <- as.data.frame(getCellColData(proj))
cell_md$Cell <- rownames(cell_md)

if (!all(c("celltype", "type", "nFrags") %in% colnames(cell_md))) {
  stop("proj cellColData must contain Lineage, Genotype, and nFrags columns.")
}

cells_use <- cell_md$Cell[!is.na(cell_md$celltype) & !is.na(cell_md$type)]
cell_md <- cell_md[match(cells_use, cell_md$Cell), , drop = FALSE]

# Extract fragment-level data without adding any matrix to ArchRProj
frags_list <- getFragmentsFromProject(
  ArchRProj = proj,
  cellNames = cells_use
)

frag_gr <- unlist(frags_list, use.names = FALSE)
frag_gr$Cell <- as.character(mcols(frag_gr)$RG)
frag_gr <- frag_gr[frag_gr$Cell %in% cells_use]

# Build TE-region x cell sparse count matrix from fragment overlaps
ov <- findOverlaps(frag_gr, te_features, ignore.strand = TRUE)
cell_index <- match(frag_gr$Cell[queryHits(ov)], cells_use)
te_index <- subjectHits(ov)
valid <- !is.na(cell_index)

te_region_counts <- sparseMatrix(
  i = te_index[valid],
  j = cell_index[valid],
  x = 1L,
  dims = c(length(te_features), length(cells_use)),
  dimnames = list(te_features$region_id, cells_use)
)

# Standardize by sequencing depth: TE counts per 10k total fragments
depth_scale <- cell_md$nFrags / 1e4
depth_scale[depth_scale == 0 | is.na(depth_scale)] <- NA_real_
te_region_norm <- te_region_counts %*% Diagonal(
  x = ifelse(is.na(depth_scale), 0, 1 / depth_scale)
)
colnames(te_region_norm) <- cells_use
rownames(te_region_norm) <- te_features$region_id

# Save region-level matrices for downstream reuse
saveRDS(te_region_counts, file.path(out_dir, "Heart_TE_region_counts_raw.rds"))
saveRDS(te_region_norm, file.path(out_dir, "Heart_TE_region_counts_norm_per10k.rds"))

# Per-cell TE scores
te_total_score <- Matrix::colSums(te_region_norm)

te_class_scores <- do.call(cbind, lapply(te_classes, function(cls) {
  idx <- which(te_features$repClass == cls)
  Matrix::colSums(te_region_norm[idx, , drop = FALSE])
}))
colnames(te_class_scores) <- paste0("TE_score_", te_classes)
rownames(te_class_scores) <- cells_use

df <- cell_md %>%
  select(Cell, celltype, type, nFrags) %>%
  mutate(TE_score_total = te_total_score[Cell]) %>%
  bind_cols(as.data.frame(te_class_scores[.$Cell, , drop = FALSE]))

write.csv(df, file.path(out_dir, "Heart_TE_scores_per_cell.csv"), row.names = FALSE)

# ---- Cell type x genotype TE accessibility bias ----
stat_res <- df %>%
  group_by(celltype) %>%
  wilcox_test(TE_score_total ~ type) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

write.csv(stat_res, file.path(out_dir, "Heart_TE_total_WT_vs_KO_wilcox.csv"), row.names = FALSE)

fc_summary <- df %>%
  group_by(celltype, type) %>%
  summarise(
    median_TE = median(TE_score_total, na.rm = TRUE),
    mean_TE = mean(TE_score_total, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = type,
    values_from = c(median_TE, mean_TE, n_cells)
  ) %>%
  mutate(
    log2FC = log2((median_TE_KO + 0.01) / (median_TE_WT + 0.01))
  ) %>%
  arrange(desc(log2FC))

write.csv(fc_summary, file.path(out_dir, "Heart_TE_total_log2FC_by_Lineage.csv"), row.names = FALSE)
print(fc_summary)

# Order lineages by KO/WT TE change
lineage_order <- fc_summary$celltype

p1 <- df %>%
  mutate(celltype = factor(celltype, levels = lineage_order)) %>%
  ggplot(aes(x = celltype, y = TE_score_total, fill = type)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3) +
  scale_fill_manual(values = c("WT" = "#4A90D9", "KO" = "#D94A4A")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(y = "TE accessibility score (per 10k fragments)", x = NULL)
p1
ggsave(
  filename = file.path(out_dir, "Heart_TE_total_boxplot.pdf"),
  plot = p1, width = 9, height = 4.5
)

class_fc_mat <- class_fc_long %>%
  select(celltype, TE_class, log2FC) %>%
  pivot_wider(names_from = TE_class, values_from = log2FC)

mat <- as.matrix(class_fc_mat[, -1, drop = FALSE])
rownames(mat) <- class_fc_mat$celltype

pdf(file.path(out_dir, "Heart_TE_class_log2FC_heatmap.pdf"), width = 6, height = 5)
Heatmap(
  mat,
  name = "log2FC",
  col = circlize::colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")),
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_title = "TE class accessibility change (KO vs WT)"
)
dev.off()

# ---- TE region-level summary by lineage and genotype ----
group_ids <- paste(cell_md$celltype, cell_md$type, sep = "__")
group_levels <- unique(group_ids)

group_mean_mat <- do.call(cbind, lapply(group_levels, function(grp) {
  idx <- which(group_ids == grp)
  Matrix::rowMeans(te_region_norm[, idx, drop = FALSE])
}))
colnames(group_mean_mat) <- group_levels
rownames(group_mean_mat) <- te_features$region_id

group_mean_df <- as.data.frame(group_mean_mat)
group_mean_df$region_id <- rownames(group_mean_df)
group_mean_df$repClass <- te_features$repClass[match(group_mean_df$region_id, te_features$region_id)]
write.csv(group_mean_df, file.path(out_dir, "Heart_TE_region_mean_norm_per10k_celltype.csv"), row.names = FALSE)

celltype <- sort(unique(cell_md$celltype))
region_fc_list <- lapply(celltype, function(ct) {
  ko_col <- paste0(ct, "__KO")
  wt_col <- paste0(ct, "__WT")
  if (!(ko_col %in% colnames(group_mean_mat) && wt_col %in% colnames(group_mean_mat))) {
    return(NULL)
  }
  data.frame(
    region_id = rownames(group_mean_mat),
    repClass = te_features$repClass,
    celltype = ct,
    mean_KO = group_mean_mat[, ko_col],
    mean_WT = group_mean_mat[, wt_col],
    log2FC = log2((group_mean_mat[, ko_col] + 0.01) / (group_mean_mat[, wt_col] + 0.01)),
    stringsAsFactors = FALSE
  )
})

region_fc_df <- bind_rows(region_fc_list) %>%
  arrange(celltype, desc(log2FC))

write.csv(region_fc_df, file.path(out_dir, "Heart_TE_region_log2FC_KO_vs_WT.csv"), row.names = FALSE)

# select log2FC > 0.5 and do heeatmap
region_fc_sub <- region_fc_df %>% filter(log2FC > 0.5)
top_n_per_celltype <- 80L
region_fc_top <- region_fc_sub %>%
  group_by(celltype) %>%
  slice_head(n = top_n_per_celltype) %>%
  ungroup()

hm_wide <- region_fc_top %>%
  select(celltype, region_id, log2FC) %>%
  pivot_wider(names_from = region_id, values_from = log2FC, values_fill = NA)
hm_mat <- as.matrix(hm_wide[, -1])
rownames(hm_mat) <- hm_wide$celltype

fc_range <- range(region_fc_top$log2FC, na.rm = TRUE)
fc_breaks <- seq(fc_range[1], fc_range[2], length.out = 3)

order <- c("CM" , "Endothelial cell"  ,  "Fibroblast"  ,       "Atrial cardiomyocyte"   ,    "Macrophage" ,     "Immune cell"    )
hm_mat <- hm_mat[order,]

library(grid)
pdf(file.path(out_dir, "Heart_TE_region_log2FC_heatmap.pdf"), width = 6, height = 6)
ht <- Heatmap(
  hm_mat,
  name = "log2FC",
  col = circlize::colorRamp2(c(fc_range[1], mean(fc_range), fc_range[2]), c("#F5F8FA", "#408FBF", "#073162")),
  cluster_columns = TRUE,
  cluster_rows = F,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 9),
  row_names_max_width = max_text_width(rownames(hm_mat), gp = gpar(fontsize = 9)),
  show_column_names = FALSE,
  show_row_names = TRUE,
  column_title = "TE region accessibility change (KO vs WT), log2FC > 0.5",
  na_col = "grey90"
)
draw(ht, padding = unit(c(4, 4, 4, 8), "mm"))
dev.off()

