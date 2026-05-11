#-----------diffHic for chr13, 30kb binsize, no replicate--------

options(scipen=999)
library(diffHic)
library(edgeR)
library(csaw)
library(statmod)
library(limma)
library(ggplot2)
library(dplyr)

bin.size <- 30000  # 30kb binsize
target.chr <- "chr13" 

wt_file <- "/media/ggj/ggj/CJY/nature_WXY/Hi-c/file/WT/matrix_file/chr13_only/wt_chr13_30bin.tsv.gz"
ko_file <- "/media/ggj/ggj/CJY/nature_WXY/Hi-c/file/WT/matrix_file/chr13_only/ko_chr13_30bin.tsv.gz"

wt_data <- read.delim(wt_file, header = FALSE)
ko_data <- read.delim(ko_file, header = FALSE)


convertToGI <- function(df){
  row.regions <- GenomicRanges::GRanges(df$V1, IRanges::IRanges(df$V2, df$V3))
  col.regions <- GenomicRanges::GRanges(df$V4, IRanges::IRanges(df$V5, df$V6))
  gi <- InteractionSet::GInteractions(row.regions, col.regions)
  gi$counts <- df$V7  # Interaction frequencies
  return(gi)
}

# turn to GInteractions 
wt.gi <- convertToGI(wt_data)
ko.gi <- convertToGI(ko_data)

cat("WT interactions:", length(wt.gi), "\n")
cat("KO interactions:", length(ko.gi), "\n")

# merge regions
hic.gi <- list(wt.gi, ko.gi)
names(hic.gi) <- c("wt.gi", "ko.gi")
combined <- unique(c(hic.gi$wt.gi, hic.gi$ko.gi))
cat(length(combined), "\n")

# replace original regions to combined regions
InteractionSet::replaceRegions(hic.gi$wt.gi) <- InteractionSet::regions(combined)
InteractionSet::replaceRegions(hic.gi$ko.gi) <- InteractionSet::regions(combined)

matched <- lapply(hic.gi, function(x) {
  match(x, combined)
})

# counts matrix
counts <- matrix(0, ncol = 2, nrow = length(combined))
counts[matched$wt.gi, 1] <- hic.gi$wt.gi$counts
counts[matched$ko.gi, 2] <- hic.gi$ko.gi$counts

# InteractionSet object
iset <- InteractionSet::InteractionSet(counts, combined)
names(SummarizedExperiment::assays(iset)) <- "counts"
InteractionSet::interactions(iset) <- methods::as(InteractionSet::interactions(iset), "ReverseStrictGInteractions")
iset$totals <- colSums(SummarizedExperiment::assay(iset))

# filter low expression interactions
keep <- edgeR::aveLogCPM(diffHic::asDGEList(iset)) > 0
iset <- iset[keep,]
cat("kept interactions:", nrow(iset), "\n")

output_dir <- "/media/ggj/ggj/CJY/nature_WXY/Hi-c/analysis/diffhic/chr13_30kb_noRep"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ========== un-norm MA plot ==========
pdf(file.path(output_dir, "chr13_30kb_prenorm_MAplot.pdf"))
ab <- edgeR::aveLogCPM(diffHic::asDGEList(iset))
o <- order(ab)
adj.counts <- edgeR::cpm(diffHic::asDGEList(iset), log=TRUE)

# WT vs KO
mval <- adj.counts[,2] - adj.counts[,1]
smoothScatter(ab, mval, xlab="A (Average logCPM)", ylab="M (logFC)", 
              main=paste("WT vs KO (pre-normalization)\n", target.chr, ", 30kb bins"))
fit <- limma::loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red", lwd=2)
abline(h=0, col="blue", lty=2)
dev.off()

# ========== normalize ==========

neardiag <- diffHic::filterDiag(iset, by.dist=2*bin.size)
nb.off <- matrix(0, nrow=nrow(iset), ncol=ncol(iset))
nb.off[neardiag] <- diffHic::normOffsets(iset[neardiag,], method="loess", se.out=F)
nb.off[!neardiag] <- diffHic::normOffsets(iset[!neardiag,], method="loess", se.out=F)
SummarizedExperiment::assay(iset, "offset") <- nb.off

saveRDS(iset, file.path(output_dir, "normalized_iset_chr13_30kb.rds"))
saveRDS(nb.off, file.path(output_dir, "normalization_offsets_chr13_30kb.rds"))

# ========== normlized MA plot ==========
pdf(file.path(output_dir, "chr13_30kb_postnorm_MAplot.pdf"))
adj.counts <- log2(SummarizedExperiment::assay(iset) + 0.5) - nb.off/log(2)
mval <- adj.counts[,2] - adj.counts[,1]
smoothScatter(ab, mval, xlab="A (Average logCPM)", ylab="M (logFC)", 
              main=paste("WT vs KO (post-normalization)\n", target.chr, ", 30kb bins"))
fit <- limma::loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red", lwd=2)
abline(h=0, col="blue", lty=2)
dev.off()

# ========== diffrential analysis ==========
iset_dgeList <- diffHic::asDGEList(iset)

group <- factor(c("WT", "KO"))
design <- model.matrix(~group)
colnames(design) <- c("Intercept", "KO_vs_WT")

fixed_disp <- 0.1

# using glmFit + glmLRT
fit.iset_dgeList <- edgeR::glmFit(iset_dgeList, design, dispersion = fixed_disp)
result.iset_dgeList <- edgeR::glmLRT(fit.iset_dgeList, coef = 2)

results_df <- edgeR::topTags(result.iset_dgeList, n = Inf, adjust.method = "BH", sort.by = "logFC")$table

cat("total interactions:", nrow(results_df), "\n")


FDR_threshold <- 0.2
logFC_threshold <- 0.4
logCPM_threshold <- 0.3

# filter: FDR -> logFC -> logCPM
results_df_FDR <- results_df[results_df$FDR < FDR_threshold & !is.na(results_df$FDR),]
results_df_logFC <- results_df_FDR[abs(results_df_FDR$logFC) > logFC_threshold,]
results_df_logFC_logCPM <- results_df_logFC[results_df_logFC$logCPM > logCPM_threshold,]
sigInts <- rownames(results_df_logFC_logCPM)

dimnames(iset_dgeList@.Data[[1]])[2][[1]] <- c("WT", "KO")

# ========== output ==========
if(length(sigInts) > 0) {
  int_coords <- as.data.frame(InteractionSet::interactions(iset))[sigInts,]
  counts_data <- iset_dgeList@.Data[[1]][sigInts,]
  bedpe <- data.frame(int_coords[,1:3], int_coords[,6:8], ".", ".", 
                     int_coords$strand1, int_coords$strand2, ".", 
                     counts_data, results_df_logFC_logCPM, stringsAsFactors=F)
  colnames(bedpe) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", 
                      "name", "score", "strand1", "strand2", "color",
                      "WT", "KO",
                      "logFC", "logCPM", "LR", "PValue", "FDR")
  bedpe[bedpe$logFC<0,"color"] <- "0,0,139" 
  bedpe[bedpe$logFC>0,"color"] <- "255,140,0"  
  
  write.table(bedpe, file.path(output_dir, "diffHic_chr13_30kb_noRep.bedpe"), 
              sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
  cat("BEDPE saved:", file.path(output_dir, "diffHic_chr13_30kb_noRep.bedpe"), "\n")

  cat("WT > KO (logFC < 0):", sum(bedpe$logFC < 0), "\n")
  cat("KO > WT (logFC > 0):", sum(bedpe$logFC > 0), "\n")
  cat("mean |logFC|:", mean(abs(bedpe$logFC)), "\n")
  cat("median |logFC|:", median(abs(bedpe$logFC)), "\n")
} else {
  cat("\nwarning: no differential interactions\n")
}

# ========== visualize ==========
# 1. Volcano plot
pdf(file.path(output_dir, "chr13_30kb_volcano_plot.pdf"), width=8, height=6)
plot(results_df$logFC, -log10(results_df$FDR), 
     xlab="log2 Fold Change (KO vs WT)", 
     ylab="-log10(FDR)",
     main=paste("Volcano Plot\n", target.chr, ", 30kb bins"),
     pch=20, cex=0.5, col=ifelse(results_df$FDR < FDR_threshold & abs(results_df$logFC) > logFC_threshold, "red", "gray"))
abline(h=-log10(FDR_threshold), col="blue", lty=2)
abline(v=c(-logFC_threshold, logFC_threshold), col="blue", lty=2)
abline(v=0, col="black", lty=1)
dev.off()

# 2. MA plot (post-normalization with significance)
pdf(file.path(output_dir, "chr13_30kb_MAplot_significant.pdf"), width=8, height=6)
plot(ab, mval, 
     xlab="A (Average logCPM)", 
     ylab="M (logFC)",
     main=paste("MA Plot with Significant Interactions\n", target.chr, ", 30kb bins"),
     pch=20, cex=0.3, col="gray")
if(length(sigInts) > 0) {
  sig_ab <- ab[sigInts]
  sig_mval <- mval[sigInts]
  points(sig_ab, sig_mval, pch=20, cex=0.5, 
         col=ifelse(sig_mval > 0, "red", "blue"))
}
abline(h=0, col="black", lty=1)
abline(h=c(-logFC_threshold, logFC_threshold), col="blue", lty=2)
legend("topright", legend=c("KO > WT", "WT > KO", "Not significant"), 
       col=c("red", "blue", "gray"), pch=20, cex=0.8)
dev.off()

# 3. distance diff
if(length(sigInts) > 0) {
  sig.df <- as.data.frame(InteractionSet::interactions(iset))[sigInts, ]
  sig.df$logFC <- results_df_logFC_logCPM$logFC
  sig.df$FDR <- results_df_logFC_logCPM$FDR

  sig.df$dist_kb <- abs(
    (sig.df$start1 + sig.df$end1)/2 -
      (sig.df$start2 + sig.df$end2)/2
  ) / 1000

  dist_bins <- seq(0, max(sig.df$dist_kb, na.rm=TRUE), by=100)
  sig.df$dist_bin <- cut(sig.df$dist_kb, breaks=dist_bins, include.lowest=TRUE)
  
  dist.stats <- sig.df %>%
    filter(!is.na(dist_bin)) %>%
    group_by(dist_bin) %>%
    summarise(
      n = n(),
      mean_logFC = mean(logFC, na.rm=TRUE),
      mean_abs_logFC = mean(abs(logFC), na.rm=TRUE),
      .groups = "drop"
    ) %>%
    mutate(dist_center = (as.numeric(sub("\\(", "", sub(",.*", "", as.character(dist_bin)))) + 
                          as.numeric(sub(".*,", "", sub("\\]", "", as.character(dist_bin))))) / 2)

  pdf(file.path(output_dir, "chr13_30kb_distance_dependent_logFC.pdf"), width=10, height=6)
  p <- ggplot(dist.stats, aes(x=dist_center, y=mean_logFC)) +
    geom_line(color="#316879", linewidth=1) +
    geom_point(color="#316879", size=2) +
    geom_hline(yintercept=0, linetype="dashed", color="gray50", linewidth=0.5) +
    labs(x="Interaction Distance (kb)",
         y="Mean log2FC (KO vs WT)",
         title=paste("Distance-dependent Interaction Difference\n", target.chr, ", 30kb bins")) +
    theme_classic(base_size=14) +
    theme(
      axis.line = element_line(linewidth=0.5),
      axis.text = element_text(size=12, color="black"),
      axis.title = element_text(size=14, face="bold"),
      plot.title = element_text(size=14, face="bold", hjust=0.5)
    )
  print(p)
  dev.off()
  
  write.table(dist.stats, file.path(output_dir, "chr13_30kb_distance_stats.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
}

all_results <- data.frame(
  as.data.frame(InteractionSet::interactions(iset)),
  results_df,
  stringsAsFactors=FALSE
)
write.table(all_results, file.path(output_dir, "chr13_30kb_all_results.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

## stremely down contact
down_threshold <- -2  
strong_down_idx <- which(mval <= down_threshold)

if (length(strong_down_idx) > 0) {
  strong_down_coords <- as.data.frame(InteractionSet::interactions(iset))[strong_down_idx, ]

  strong_down_res <- data.frame(
    strong_down_coords,
    A = ab[strong_down_idx],
    M = mval[strong_down_idx],
    stringsAsFactors = FALSE
  )
  write.table(
    strong_down_res,
    file.path(output_dir, "chr13_30kb_strong_down_from_MA.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}
