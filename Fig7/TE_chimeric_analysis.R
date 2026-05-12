library(dplyr)
library(readr)
library(Seurat)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(GenomicRanges)
library(data.table)
library(ggplot2)

####################### celltype zscore heatmap########################

log2fc_dir <- "/path/scRNA/analysis/TE_chimeric/tissue/results/log2FC"
log2fc_files <- list.files(log2fc_dir, pattern = "^Heart_.*_TE-derived_log2FC\\.csv$", full.names = TRUE)

gene_sets <- list()
celltypes <- character(0)
for (f in log2fc_files) {
  ct <- sub("^Heart_(.*)_TE-derived_log2FC\\.csv$", "\\1", basename(f))
  d <- read.csv(f, stringsAsFactors = FALSE)
  if (!"Gene" %in% colnames(d)) next
  genes_ct <- unique(na.omit(d$Gene))
  if (length(genes_ct) == 0) next
  gene_sets[[ct]] <- genes_ct
  celltypes <- c(celltypes, ct)
}
celltypes <- unique(celltypes)

cat("Heart celltypes with TE-derived genes:", paste(celltypes, collapse = ", "), "\n")

load("/path/scRNA/merge/Tissue/Heart/seob_anno_v4.rda")
meta <- seob@meta.data

normalize_ct_name <- function(x) gsub("[ /]", "_", x)

meta$cell_type_norm <- normalize_ct_name(meta$cell_type)
celltypes_norm <- normalize_ct_name(celltypes)

cells_use <- rownames(meta)[meta$cell_type_norm %in% celltypes_norm]
seob_sub <- subset(seob, cells = cells_use)
meta_sub <- seob_sub@meta.data

meta_sub$cell_type_norm <- normalize_ct_name(meta_sub$cell_type)
seob_sub$cell_type_norm <- meta_sub$cell_type_norm

exp_mat <- as.matrix(GetAssayData(seob_sub, slot = "data"))  # log-normalized expr
all_genes <- rownames(exp_mat)

## 3. calculate (source_celltype, gene) expression in Heart celltype
avg_list <- list()
for (ct in celltypes) {
  ct_norm <- normalize_ct_name(ct)
  genes_ct <- intersect(gene_sets[[ct]], all_genes)
  if (length(genes_ct) == 0) next

  cells_ct <- colnames(seob_sub)[seob_sub$cell_type_norm == ct_norm]
  if (length(cells_ct) == 0) next

  mat_ct <- exp_mat[genes_ct, cells_ct, drop = FALSE]
  avg_expr <- rowMeans(mat_ct, na.rm = TRUE)

  avg_list[[ct]] <- data.frame(
    source = ct,                
    gene   = names(avg_expr),
    celltype = ct,             
    value  = as.numeric(avg_expr),
    stringsAsFactors = FALSE
  )
}

if (length(avg_list) == 0) {
  stop("No valid gene expression data found for Heart celltypes.")
}

df_avg <- bind_rows(avg_list)

## 4.(source,gene) × celltype matrix
mat_wide <- dcast(df_avg, source + gene ~ celltype, value.var = "value", fill = 0)
mat_wide$gene_id <- paste(mat_wide$source, mat_wide$gene, sep = "|")
rownames(mat_wide) <- mat_wide$gene_id
mat_wide$source <- NULL
mat_wide$gene   <- NULL
mat_wide$gene_id <- NULL
mat_mat <- as.matrix(mat_wide)

row_zscore <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  mu <- mean(x, na.rm = TRUE)
  sdv <- sd(x, na.rm = TRUE)
  if (sdv == 0 || is.na(sdv)) {
    return(rep(0, length(x)))
  } else {
    (x - mu) / sdv
  }
}
mat_z <- t(apply(mat_mat, 1, row_zscore))
rownames(mat_z) <- rownames(mat_mat)
colnames(mat_z) <- colnames(mat_mat)

keep_rows <- apply(mat_z, 1, function(x) any(!is.na(x) & x != 0))
mat_z <- mat_z[keep_rows, , drop = FALSE]

cat("Heatmap matrix dimension (genes × celltypes):", nrow(mat_z), "×", ncol(mat_z), "\n")

row_source <- sub("\\|.*$", "", rownames(mat_z))
row_source_factor <- factor(row_source, levels = unique(row_source))

n_src <- length(levels(row_source_factor))
src_cols <- circlize::colorRamp2(
  seq(0, 1, length.out = n_src),
  grDevices::hcl(
    h = seq(0, 300, length.out = n_src),
    c = 40,   
    l = 80   
  )
)(seq(0, 1, length.out = n_src))
names(src_cols) <- levels(row_source_factor)

row_ha <- rowAnnotation(
  source = row_source_factor,
  col = list(source = src_cols),
  show_annotation_name = FALSE  
)

## 6.do heatmap
col_fun <- colorRamp2(c(-2, 0, 2), c("#4575B4", "white", "#D73027"))

ht <- Heatmap(
  mat_z,
  name = "z-score",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_title = "TE-derived genes mean expression (z-score) across Heart celltypes",
  row_split = row_source_factor,
  left_annotation = row_ha
)

pdf("/path/scRNA/analysis/TE_chimeric/tissue/analysis/Heart_TEderived_genes_celltype_zscore_heatmap.pdf",
    width = 6, height = 6)
draw(ht)
dev.off()


######################ORF analysis########################

library(data.table)
library(dplyr)
library(ggplot2)

gene_stat_file <- "/path/scRNA/analysis/TE_chimeric/tissue/analysis/gene_count_df_stat.csv"
kozac_cpc2_file <- "/path/scRNA/analysis/TE_chimeric/tissue/Heart/CM/Kozak_CPC2_final_peptide.csv"
out_dir_orf <- "/path/scRNA/analysis/TE_chimeric/tissue/analysis/CM_ORF"

dir.create(out_dir_orf, recursive = TRUE, showWarnings = FALSE)

gene_stat_dt <- fread(gene_stat_file)
ko_ids <- unique(gene_stat_dt$Name[gene_stat_dt$type == "KO"])

kozac_dt <- fread(kozac_cpc2_file)
kozac_ko <- kozac_dt[uniqid %in% ko_ids]
kozac_ko=kozac_ko[kozac_ko$type!="None",]

message("Total TE-chimeric KO transcripts with ORF/CPC2 annotation: ", nrow(kozac_ko))

## TE-chimeric transcript coding vs noncoding
if ("cpc2_label" %in% colnames(kozac_ko)) {
  kozac_ko <- kozac_ko %>%
    mutate(coding_status = cpc2_label)
} else if ("cpc2_is_coding" %in% colnames(kozac_ko)) {
  kozac_ko <- kozac_ko %>%
    mutate(coding_status = ifelse(cpc2_is_coding, "coding", "noncoding"))
} else {
  stop("Cannot find CPC2 coding label columns (cpc2_label or cpc2_is_coding) in Kozak_CPC2_final_peptide.csv")
}

df_coding <- kozac_ko %>%
  filter(!is.na(coding_status)) %>%
  count(coding_status, name = "n") %>%
  mutate(
    fraction = n / sum(n),
    label = paste0(coding_status, "\n", n)
  )


##  coding transcript: in-frame / out-of-frame fraction
df_coding_frame <- kozac_ko %>%
  filter(!is.na(type), type != "None", !is.na(cpc2_label), cpc2_label == "coding") %>%
  mutate(frame_status = ifelse(type == "out-of-frame", "out-of-frame", "in-frame")) %>%
  count(frame_status, name = "n") %>%
  mutate(
    frac = n / sum(n),
    label = paste0(frame_status, "\n", n)
  )

if (nrow(df_coding_frame) > 0) {
  df_coding_frame$frame_status <- factor(
    df_coding_frame$frame_status,
    levels = c("in-frame", "out-of-frame")
  )

  pdf(file.path(out_dir_orf, "Heart_CM_coding_frame_pie.pdf"), width = 4, height = 4)
  ggplot(df_coding_frame, aes(x = "", y = frac, fill = frame_status)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 3
    ) +
    scale_fill_manual(
      values = c("in-frame" = "#D94A4A", "out-of-frame" = "grey70"),
      name = NULL,
      drop = FALSE
    ) +
    theme_void() +
    ggtitle("Frame status in CPC2 coding TE-chimeric transcripts\n(CM, KO)")
  dev.off()
}

## NMD sensitivity in in-frame coding transcripts
if ("nmdcall" %in% colnames(kozac_ko)) {
  df_nmd <- kozac_ko %>%
    mutate(
      NMD_status = case_when(
        nmdcall %in% c("Yes", "yes", "NMD_sensitive", "NMD-sensitive") ~ "NMD-sensitive",
        nmdcall %in% c("No", "no", "NMD_insensitive", "NMD-insensitive") ~ "NMD-insensitive",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(
      frame_group == "in-frame",
      coding_status == "coding",
      !is.na(NMD_status)
    ) %>%
    count(NMD_status, name = "n") %>%
    mutate(
      frac = n / sum(n),
      label = paste0(NMD_status, "\n", n)
    )

  nmd_order <- c("NMD-sensitive", "NMD-insensitive")
  df_nmd$NMD_status <- factor(df_nmd$NMD_status, levels = nmd_order)

  pdf(file.path(out_dir_orf, "Heart_CM_TEchimeric_NMD_inframe_coding_pie.pdf"), width = 4.5, height = 4.5)
  ggplot(df_nmd, aes(x = "", y = frac, fill = NMD_status)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5),
              size = 3) +
    scale_fill_manual(
      values = c("NMD-sensitive" = "#4E79A7", "NMD-insensitive" = "#A0CBE8"),
      name = NULL
    ) +
    theme_void() +
    ggtitle("NMD sensitivity in in-frame coding TE-chimeric transcripts\n(CM, KO)")
  dev.off()
}


## coding transcripts type, with NMD sensitive
pdf(file.path(out_dir_orf, "Heart_CM_TEchimeric_type_in_coding_inframe2.pdf"), width = 5.5, height = 4)
ggplot(df_type_coding, aes(x = type, y = n, fill = NMD_status)) +
  geom_col() +
  geom_text(
    data = df_type_coding_total,
    inherit.aes = FALSE,
    aes(x = type, y = total, label = total),
    vjust = -0.2,
    size = 3
  ) +
  scale_fill_manual(
    values = c("NMD-sensitive" = "#4E79A7", "NMD-insensitive" = "#A0CBE8"),
    name = NULL
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Chimeric transcript type",
    y = "Number of coding TE-chimeric transcripts",
    title = "Chimeric transcript types among coding TE-chimeric transcripts by NMD status\n(CM, KO)"
  )
dev.off()

#--------------- junction reads statistic ------------------


gene_stat_df <- read_csv(
  "/path/scRNA/analysis/TE_chimeric/tissue/analysis/gene_count_df_stat.csv",
  show_col_types = FALSE
)

gene_stat_df <- gene_stat_df %>%
  transmute(
    celltype = ct.x,
    transcript_name = Name,
    genotype = type
  ) %>%
  filter(genotype %in% c("KO", "WT"))

isoform_files <- list.files(
  "/path/scRNA/analysis/TE_chimeric/tissue/results",
  pattern = "^Heart_.*_TE-derived_Alternative_Isoforms\\.csv$",
  full.names = TRUE
)

isoform_plot_list <- list()

for (f in isoform_files) {
  df_iso <- read_csv(f, show_col_types = FALSE)

  req_cols <- c(
    "Transcript Name",
    "Mean Intron Read Count (Treatment)",
    "Mean Intron Read Count (Normal)"
  )
  if (!all(req_cols %in% colnames(df_iso))) next

  celltype_label <- sub(
    "^Heart_(.*)_TE-derived_Alternative_Isoforms\\.csv$",
    "\\1",
    basename(f)
  )

  df_long <- df_iso %>%
    transmute(
      celltype = celltype_label,
      transcript_name = `Transcript Name`,
      tumor_intronjuncount_mean = `Mean Intron Read Count (Treatment)`,
      normal_intronjuncount_mean = `Mean Intron Read Count (Normal)`
    ) %>%
    left_join(gene_stat_df, by = c("celltype", "transcript_name")) %>%
    filter(genotype %in% c("KO", "WT")) %>%
    mutate(
      intron_read_count = ifelse(
        genotype == "KO",
        tumor_intronjuncount_mean,
        normal_intronjuncount_mean
      )
    ) %>%
    select(celltype, genotype, intron_read_count)

  isoform_plot_list[[length(isoform_plot_list) + 1L]] <- df_long
}

if (length(isoform_plot_list) > 0) {
  df_isoform_plot <- bind_rows(isoform_plot_list) %>%
    filter(!is.na(intron_read_count))

  cell_order <- df_isoform_plot %>%
    group_by(celltype) %>%
    summarise(mean_count = mean(intron_read_count, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_count)) %>%
    pull(celltype)

  df_isoform_plot$celltype <- factor(df_isoform_plot$celltype, levels = cell_order)
  df_isoform_plot$genotype <- factor(df_isoform_plot$genotype, levels = c("KO", "WT"))

  p_intron_box <- ggplot(
    df_isoform_plot,
    aes(x = celltype, y = intron_read_count, fill = genotype)
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.72),
      width = 0.62,
      outlier.size = 0.4,
      alpha = 0.9
    ) +
    scale_fill_manual(values = c("KO" = "#c9796f", "WT" = "#6c8ebf")) +
    coord_cartesian(ylim = c(0, 40)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank()
    ) +
    labs(
      x = "Cell type",
      y = "Mean intron junction count",
      title = "Mean intron junction count by source genotype in Heart TE-derived isoforms"
    )

  pdf(
    "/path/scRNA/analysis/TE_chimeric/tissue/analysis/TE_isoform_intronjuncount_by_source_boxplot.pdf",
    width = 8,
    height = 4.5
  )
  print(p_intron_box)
  dev.off()
}


