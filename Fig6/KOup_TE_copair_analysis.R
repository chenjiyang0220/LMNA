#############Fig6.5
proj_callpeak <- readRDS("./scATAC/tissue/Heart/proj_integration_Co_addToArrow.rds")

#### ATAC marker gene heatmap

markersGS <- getMarkerFeatures(
    ArchRProj = proj_callpeak, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "celltype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerGenes_to_label <- c(
  # Atrial CM
  "Nppa", "Myl4",
  # CM
  "Myh7", "Myh6", "Ttn",
  # Endothelial
  "Pecam1", "Egfl7",
  # Fibroblast
  "Col1a1", "Col3a1", "Dcn",
  # Macrophage
  "Csf1r", "C1qc",
  # Immune
  "Cd79a", "Cd3e"
)
heatmapMarkers <- plotMarkerHeatmap(
  seMarker   = markersGS, 
  cutOff     = "FDR <= 0.05 & Log2FC >= 0.5",
  transpose  = TRUE,
  plotLog2FC = TRUE,labelMarkers = markerGenes_to_label , clusterCols = TRUE
)
heatmapMarkers <-plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  transpose = TRUE,labelMarkers = markerGenes_to_label ,  clusterCols = F
)
pdf("./scATAC/tissue/Heart/marker_gene_heatmap2.pdf",width = 8,height = 6)
print(heatmapMarkers)
dev.off()

## ------- marker peak -------
for (i in unique(proj_callpeak$celltype)) {
  proj <- proj_callpeak[proj_callpeak$celltype == i]
  markersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    maxCells = 2000,
    groupBy = "type",
    bias = c("TSSEnrichment", "log10(nFrags)"), 
    testMethod = "wilcoxon"
  )
  markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 0.5", returnGR = TRUE) # return=FALSE -> generate a data frame
  markerList <- as.data.frame(markerList)
  write.csv(markerList,paste0('./scATAC/tissue/Heart/marker_celltype/',i,'_markerPeaks.csv'),quote = F)
}

setwd("./scATAC/tissue/Heart/marker_celltype")

# ========== export TE overlap's KO up marker ==========
library(GenomicRanges)
library(data.table)

te_bed   <- "./scATAC/TE/TE_pantissue_up/mm10.te.bed"
marker_dir <- "./scATAC/tissue/Heart/marker_celltype"
out_te   <- "./scATAC/tissue/Heart/marker_celltype/TE_KOup"

dir.create(out_te, recursive = TRUE, showWarnings = FALSE)

te <- fread(te_bed, header = FALSE)
te_gr <- GRanges(seqnames = te[[1]], IRanges(as.numeric(te[[2]]), as.numeric(te[[3]])))

marker_files <- list.files(marker_dir, pattern = "_markerPeaks\\.csv$", full.names = TRUE)
for (f in marker_files) {
  celltype <- sub("_markerPeaks\\.csv$", "", basename(f))
  d <- read.csv(f)
  if (!"group_name" %in% names(d)) next
  ko <- d[d$group_name == "KO", ]
  if (nrow(ko) == 0) next
  ko_gr <- GRanges(seqnames = ko$seqnames, IRanges(ko$start, ko$end), mcols = ko)
  # ol <- findOverlaps(te_gr,ko_gr,  ignore.strand = TRUE)
  # idx <- unique(queryHits(ol))
  # ko_te <- ko[idx, ]
  ko_te <- subsetByOverlaps(te_gr,ko_gr)
  ko_te <- as.data.frame(ko_te)
  if (nrow(ko_te) == 0) next
  write.csv(ko_te, file.path(out_te, paste0(celltype, "_KOup_TEoverlap_markerPeaks.csv")), quote = FALSE)
  write.table(ko_te[, c("seqnames", "start", "end")], file.path(out_te, paste0(celltype, "_KOup_TEoverlap_markerPeaks.bed")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  message(celltype, ": KO-up ", nrow(ko), " -> TE overlap ", nrow(ko_te))
}

########anno and export bed

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

getTxDbGenes <- function(txdb, orgdb, gr = NULL, ignore.strand = TRUE) {
  if (is.null(gr)) {
    genes <- GenomicFeatures::genes(txdb)
  } else {
    genes <- suppressWarnings(
      subsetByOverlaps(
        GenomicFeatures::genes(txdb),
        gr,
        ignore.strand = ignore.strand
      )
    )
  }

  if (length(genes) > 1) {
    mcols(genes)$symbol <- suppressMessages(
      mapIds(
        orgdb,
        keys = mcols(genes)$gene_id,
        column = "SYMBOL",
        keytype = "ENTREZID",
        multiVals = "first"
      )
    )
    genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)
    names(genes) <- NULL
    out <- genes
  } else {
    out <- GRanges(seqnames(gr), ranges = IRanges(0, 0), gene_id = 0, symbol = "none")[-1]
  }
  out
}

genes_txdb <- getTxDbGenes(txdb = txdb, orgdb = org.Mm.eg.db)
id_symbol  <- data.table(geneId = genes_txdb$gene_id, symbol = genes_txdb$symbol)

te_koup_peak_dir <- "./scATAC/tissue/Heart/marker_celltype/TE_KOup"
te_peak_files <- list.files(te_koup_peak_dir, pattern = "_KOup_TEoverlap_markerPeaks\\.bed$", full.names = TRUE)

for (bf in te_peak_files) {
  celltype <- sub("_KOup_TEoverlap_markerPeaks\\.bed$", "", basename(bf))
  message("Annotating TE_KOup peaks for cell type: ", celltype)

  peak_gr <- readPeakFile(bf)
  peak_anno <- annotatePeak(
    peak_gr,
    TxDb = txdb,
    tssRegion = c(-2000, 200),
    verbose = FALSE
  )
  anno_df <- as.data.frame(peak_anno)
  if ("geneId" %in% names(anno_df)) {
    anno_df$symbol <- id_symbol$symbol[match(anno_df$geneId, id_symbol$geneId)]
  }

  write.table(
      anno_df,
      file = file.path(te_koup_peak_dir, paste0(celltype, "_KOupTE_anno.txt")),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  is_promoter <- grepl("Promoter", anno_df$annotation, ignore.case = TRUE)

  prom_df <- anno_df[is_promoter, ]
  enh_df  <- anno_df[!is_promoter, ]

  if (nrow(prom_df) > 0) {
    prom_bed <- prom_df[, c("seqnames", "start", "end")]
    write.table(
      prom_bed,
      file = file.path(te_koup_peak_dir, paste0(celltype, "_KOupTE_promoter_anno.bed")),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
  }

  if (nrow(enh_df) > 0) {
    enh_bed <- enh_df[, c("seqnames", "start", "end")]
    write.table(
      enh_bed,
      file = file.path(te_koup_peak_dir, paste0(celltype, "_KOupTE_Enhancer_anno.bed")),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
  }
}
########################################################
##############run HOMER(only Promoter)###############
########################################################

#################### Enhancer with Hic valid overlap with KOup TE ####################
## Hi-C validated Cicero copair 

hic_base_dir <- "./Hi-c/analysis/dots/result"
cicero_base_dir <- "./scATAC/tissue/Heart/cicero_downsample/result/process"
# cicero_base_dir <- "./scATAC/tissue/Heart/cicero_celltype/result/process"
out_dir <- "./scATAC/tissue/Heart/marker_celltype/TE_KOup_enhancer"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# samples_hic_cicero <- c("Atrial cardiomyocyte_KO", "Endothelial cell_KO", "Immune cell_KO", "Macrophage_KO")
samples_hic_cicero <- c("CM_KO", "Fibroblast_KO")

cicero_colnames <- c(
  "cicero_id", "peak1", "peak2", "coaccess",
  "is_promoter1", "is_promoter2", "interaction_count",
  "gene1", "gene2", "unique_id",
  "center1", "center2", "distance", "dist_class",
  "annot1_detail", "annot1_type", "annot1_gene",
  "annot2_detail", "annot2_type", "annot2_gene",
 "pair_type"
)

for (sample in samples_hic_cicero) {
  message("Merge Hi-C valid copairs with Cicero annotations for sample: ", sample)
  hic_file <- file.path(
    hic_base_dir,
    sample,
    "cicero_hic_validation_allChrom_interaction_regions_merged.tsv"
  )
  cicero_file <- file.path(
    cicero_base_dir,
    paste0(sample, "_2kcell.txt")
  )
  cicero_dt$cicero_id=as.numeric(cicero_dt$cicero_id)
  hic_dt <- fread(hic_file)
  cicero_dt <- fread(cicero_file, header = FALSE, col.names = cicero_colnames)
  cicero_dt$cicero_id=as.numeric(cicero_dt$cicero_id)
  if (!("cicero_id" %in% colnames(hic_dt))) {
    stop("Hi-C validation file does not contain 'cicero_id' column: ", hic_file)
  }

  merged_dt <- merge(
    hic_dt,
    cicero_dt,
    by = "cicero_id",
    all.x = TRUE,
    sort = FALSE
  )
  merged_dt <- merged_dt[merged_dt$is_valid == 'TRUE',]
  message("  Valid copairs: ", nrow(merged_dt))

  ## only EP（Promoter-Enhancer） remain
  if ("pair_type.x" %in% colnames(merged_dt)) {
    merged_dt <- merged_dt[pair_type.x == "EP", ]
  } else if ("pair_type" %in% colnames(merged_dt)) {
    merged_dt <- merged_dt[pair_type == "EP", ]
  }
  message("  EP-only copairs: ", nrow(merged_dt))
  is_prom1 <- !is.na(merged_dt$annot1_type) &
    grepl("Promoter", merged_dt$annot1_type, ignore.case = TRUE)
  is_prom2 <- !is.na(merged_dt$annot2_type) &
    grepl("Promoter", merged_dt$annot2_type, ignore.case = TRUE)

  merged_dt[, promoter_peak := fifelse(
    is_prom1, peak1.x,
    fifelse(is_prom2, peak2.x, NA_character_)
  )]
  merged_dt[, enhancer_peak := fifelse(
    is_prom1, peak2.x,
    fifelse(is_prom2, peak1.x, NA_character_)
  )]
  merged_dt[, promoter_gene := fifelse(
    is_prom1, annot1_gene,
    fifelse(is_prom2, annot2_gene, NA_character_)
  )]
  merged_dt[, enhancer_gene := fifelse(
    is_prom1, annot2_gene,
    fifelse(is_prom2, annot1_gene, NA_character_)
  )]

  dup_y_cols <- c(
    "peak1.y", "peak2.y", "coaccess.y", "gene1.y", "gene2.y",
    "center1.y", "center2.y", "distance.y", "pair_type.y"
  )
  dup_y_cols <- intersect(dup_y_cols, colnames(merged_dt))
  if (length(dup_y_cols) > 0) {
    merged_dt[, (dup_y_cols) := NULL]
  }
  rename_map <- c(
    "peak1.x"    = "peak1",
    "peak2.x"    = "peak2",
    "center1.x"  = "center1",
    "center2.x"  = "center2",
    "coaccess.x" = "coaccess",
    "distance.x" = "distance",
    "pair_type.x"= "pair_type",
    "gene1.x"    = "gene1",
    "gene2.x"    = "gene2"
  )
  rename_old <- intersect(names(rename_map), colnames(merged_dt))
  if (length(rename_old) > 0) {
    setnames(merged_dt, old = rename_old, new = unname(rename_map[rename_old]))
  }

  out_file <- paste0(
    out_dir,'/',sample,
    "_cicero_hic_validation_prom_enh.tsv"
  )
  fwrite(merged_dt, out_file, sep = "\t")
  message("  Merged table written to: ", out_file)
}


te_koup_dir <- "./scATAC/tissue/Heart/marker_celltype/TE_KOup"
enh_out_dir <- "./scATAC/tissue/Heart/marker_celltype/TE_KOup_enhancer"
h3k27ac_peak_file <- "./bulkCUT/peak_shuf/KO/H3K27ac/Heart/Heart_peaks.narrowPeak"
dir.create(enh_out_dir, recursive = TRUE, showWarnings = FALSE)
samples_hic_cicero <- c("CM_KO", "Fibroblast_KO","Atrial cardiomyocyte_KO", "Endothelial cell_KO", "Immune cell_KO", "Macrophage_KO")

for (sample in samples_hic_cicero) {
  message(sample)

  celltype_name <- sub("_KO$", "", sample)
  te_bed_file <- file.path(
    te_koup_dir,
    paste0(celltype_name, "_KOup_TEoverlap_markerPeaks.bed")
  )
  prom_enh_file <- file.path(
    enh_out_dir,
    paste0(sample, "_cicero_hic_validation_prom_enh.tsv")
  )

  te_dt <- fread(te_bed_file, header = FALSE)
  setnames(te_dt, c("V1", "V2", "V3"), c("chr", "start", "end"))

  prom_enh_dt <- fread(prom_enh_file)
  enh_peaks <- unique(na.omit(prom_enh_dt$enhancer_peak))


  enh_split <- tstrsplit(enh_peaks, "_", fixed = TRUE)
  enh_dt <- data.table(
    chr   = enh_split[[1]],
    start = suppressWarnings(as.integer(enh_split[[2]])),
    end   = suppressWarnings(as.integer(enh_split[[3]]))
  )
  enh_dt <- enh_dt[!is.na(start) & !is.na(end)]
  if (nrow(enh_dt) == 0) {
    message("  All enhancer_peak entries could not be parsed into coordinates for sample: ", sample)
    next
  }

  te_gr <- GRanges(
    seqnames = te_dt$chr,
    ranges   = IRanges(start = te_dt$start, end = te_dt$end)
  )
  enh_gr <- GRanges(
    seqnames = enh_dt$chr,
    ranges   = IRanges(start = enh_dt$start, end = enh_dt$end)
  )

  ov <- findOverlaps(te_gr, enh_gr, ignore.strand = TRUE)
  if (length(ov) == 0) {
    message("  No overlaps between TE_KOup peaks and Hi-C enhancers for sample: ", sample)
    next
  }
  te_idx <- sort(unique(queryHits(ov)))
  te_overlaps <- te_dt[te_idx, ]
  message("  TE_KOup enhancer-overlapping peaks: ", nrow(te_overlaps))

  enh_bed_out <- file.path(
    enh_out_dir,
    paste0(sample, "_TE_KOup_enhancer_overlap.bed")
  )
  fwrite(te_overlaps[, .(chr, start, end)], enh_bed_out,
         sep = "\t", col.names = FALSE)
  message("  TE_KOup enhancer-overlapping peaks written to: ", enh_bed_out)

  ## bulk KO H3K27ac peaks overlaps
  if (!file.exists(h3k27ac_peak_file)) {
    message("  H3K27ac peak file not found, skip H3K27ac overlap: ", h3k27ac_peak_file)
  } else {
    h3_dt <- fread(h3k27ac_peak_file, header = FALSE)
    if (nrow(h3_dt) > 0) {
      setnames(h3_dt, c("V1", "V2", "V3"), c("chr", "start", "end"))
      h3_gr <- GRanges(
        seqnames = h3_dt$chr,
        ranges   = IRanges(start = h3_dt$start, end = h3_dt$end)
      )
      te_gr2 <- GRanges(
        seqnames = te_overlaps$chr,
        ranges   = IRanges(start = te_overlaps$start, end = te_overlaps$end)
      )
      ov_h3 <- findOverlaps(te_gr2, h3_gr, ignore.strand = TRUE)
      if (length(ov_h3) > 0) {
        te_h3_idx <- sort(unique(queryHits(ov_h3)))
        te_h3_overlaps <- te_overlaps[te_h3_idx, ]
        h3_bed_out <- file.path(
          enh_out_dir,
          paste0(sample, "_TE_KOup_enhancer_H3K27ac_overlap.bed")
        )
        fwrite(te_h3_overlaps[, .(chr, start, end)], h3_bed_out,
               sep = "\t", col.names = FALSE)
        message("  TE_KOup enhancer peaks overlapping H3K27ac written to: ", h3_bed_out)
      } else {
        message("  No overlaps between TE_KOup enhancer peaks and H3K27ac for sample: ", sample)
      }
    }
  }
}

WTandKO=fread('./Hi-c/analysis/dots/result/compare/WT_KO_EP_promoter_pairs_with_overlap.tsv')
head(WTandKO)

wtko_sub <- WTandKO[overlap_status %in% c("wt_unique", "shared")]

co_split <- tstrsplit(wtko_sub$co_peak_region, "_", fixed = TRUE)
bed_dt <- data.table(
  chr   = co_split[[1]],
  start = as.integer(co_split[[2]]),
  end   = as.integer(co_split[[3]])
)
out_dir <- "./scATAC/tissue/Heart/marker_celltype/TE_KOup_anno_homer/CM"
out_bed <- file.path(out_dir, "CM_WT_enhancer_region.bed")
fwrite(bed_dt, out_bed, sep = "\t", col.names = FALSE)

message("BED written to: ", out_bed)

#################run HOMER######################
