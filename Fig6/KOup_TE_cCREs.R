library(ArchR)
library(readr)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)


## ----------- TE KOup cCREs ------------
addArchRThreads(threads = 1)
proj_callpeak <- readRDS("./scATAC/proj_callpeak_all.rds")
proj_callpeak <- addPeakMatrix(proj_callpeak)
setwd("./scATAC")
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj_callpeak,
  useMatrix = "PeakMatrix",
  maxCells = 10000,
  groupBy = "type",
  bias = c("TSSEnrichment", "log10(nFrags)"), #让ArchR考虑细胞群之间数据质量的差异
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 0.5", returnGR = TRUE) # return=FALSE -> generate a data frame
markerList <- as.data.frame(markerList)
write.csv(markerList,'./diffpeak/all_KOvsWT/pan_tissue_marker_202603.csv',quote = F)
marker.KO <- markerList[markerList$group_name=='KO',]
marker_KO.gr <- GRanges(seqnames = marker.KO$seqnames,IRanges(marker.KO$start,marker.KO$end))

TE_anno <- fread("/media/ggj/ggj/CJY/tools/mm10/mm10.te.txt.gz")
TE_gr <- GRanges(seqnames = TE_anno$V6, ranges = IRanges(start = TE_anno$V7, end = TE_anno$V8))
TE_gr$type <- TE_anno$V12  
TE_gr$names <- TE_anno$V11 

markerKO_TE.gr <- subsetByOverlaps(marker_KO.gr,TE_gr)
ol <- findOverlaps(markerKO_TE.gr, TE_gr)

markerKO_TE.gr$te_family <- NA
qh <- queryHits(ol)
sh <- subjectHits(ol)
te_family_list <- tapply(TE_gr$type[sh], qh, function(x) paste(unique(x), collapse = ";"))
markerKO_TE.gr$te_family[as.integer(names(te_family_list))] <- unname(te_family_list)
markerKO_TE <- as_tibble(markerKO_TE.gr)

markerKO_TE$id <- paste(markerKO_TE$seqnames,markerKO_TE$start,markerKO_TE$end,sep = '_')
marker.KO$id <- paste(marker.KO$seqnames,marker.KO$start,marker.KO$end,sep = '_')
markerKO_TE_merge <- merge(markerKO_TE,marker.KO,by='id')
write.csv(markerKO_TE_merge,'./TE/TE_pantissue_up/KOup_TE_cCREs.csv',quote = F)


###### KO up TE cCREs for each histone modify

peak_dir <- './bulkCUT/peak_noshuf/KO/merge'

peak_file <- list.files(path = peak_dir)
his <- c('CTCF','H3K27ac','H3K27me3','H3K4me1','H3K4me3')
df <- data.frame()
for (i in peak_file) {
  message(i)
  histone <- gsub('.fragments.sorted_peaks.narrowPeak','',i)
  histone <- gsub('KO.','',histone)
  peak <- fread(paste(peak_dir,i,sep='/'))
  peak.gr <- GRanges(seqnames = peak$V1,IRanges(peak$V2,peak$V3))
  hist_KOupTE.gr <- subsetByOverlaps(markerKO_TE.gr,peak.gr)
  hist_KOupTE <- as_tibble(hist_KOupTE.gr)
  hist_KOupTE$histone_modify <- histone
  write.csv(hist_KOupTE,file = paste0('./scATAC/TE/TE_diffpeak_mat/histone_KOup_TE_family/',histone,'_KOupTE_info.csv'),quote = F)
  for (j in c("LTR",'LINE','SINE')) {
    sep <- hist_KOupTE[grepl(j,hist_KOupTE$te_family),]
    write.table(sep,file = paste0('./scATAC/TE/TE_diffpeak_mat/histone_KOup_TE_family/',j,'/',histone,'_',j,'_KOupTE.bed'),quote = F,sep = '\t',col.names = F,row.names = F)
  }

  tmp <- data.frame(histone=histone,
                    number_of_KOupTE=length(hist_KOupTE.gr))
  df <- rbind(df,tmp)
  family <- hist_KOupTE.gr$te_family
  fam_count <- as.data.frame(sort(table(family), decreasing = TRUE))
  fam_count <- fam_count[fam_count$family %in% fam_count$family[grepl('SINE|LINE|LTR',fam_count$family)],]
  fam_count <- fam_count[fam_count$family %ni% fam_count$family[grepl('Simple_repeat|DNA|Unknown|Low_complexity|snRNA|Satellite|Other|scRNA|tRNA',fam_count$family)],]

  ggplot(fam_count, aes(x = family, y = Freq)) +
    geom_col(fill = "steelblue") +
    xlab("TE family") +
    ylab("Freq") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0('./scATAC/TE/TE_diffpeak_mat/histone_KOup_TE_family/',histone,'_KOupTE_famliy_barplot.pdf'),width = 5,height = 5)
}
write.csv(df,'./scATAC/TE/TE_diffpeak_mat/histone_KOup_TE_family/each_histone_TE_number.csv',quote = F)

########### go to HOMER analysis#############################
