setwd("J:/WXY/scATAC/")
library(ArchR)
library(Matrix)
library(reshape2)
addArchRThreads(threads = 6)
proj_callpeak<-readRDS("proj_callpeak_all.rds")

## ------------- ACD coaccess heatmap ---------------
Anno_type<-unique(proj_callpeak$Anno_type)
Anno_type<-Anno_type[-grep("APC",Anno_type)]
Anno<-unique(colsplit(Anno_type,"_",c(1,2))$`1`)
i=28
for (i in 1:30) {
tissue_counts <- table(proj_callpeak$tissue[proj_callpeak$Anno_type == Anno_type[2*i-1]])
selected_tissues <- names(tissue_counts[tissue_counts >= 10])
proj_callpeak_c_wt <- proj_callpeak[proj_callpeak$Anno_type == Anno_type[2*i-1], ]
proj_callpeak_c_wt <- proj_callpeak_c_wt[proj_callpeak_c_wt$tissue %in% selected_tissues, ]
if (length(rownames(proj_callpeak_c_wt@cellColData)) >= 1000) {
  proj_callpeak_c_wt <- proj_callpeak_c_wt[sample(1:length(rownames(proj_callpeak_c_wt@cellColData)), 1000), ]
}
mm<-getMatrixFromProject(ArchRProj=proj_callpeak_c_wt,useMatrix = "PeakMatrix")
test<-mm@assays@data
test<-test[["PeakMatrix"]]
acd_wt<-read.csv(paste0("ACD/",Anno_type[2*i-1],"_ACD_peak.csv"))

library(dplyr)
library(GenomicRanges)

# WT
peak_gr<-as.data.frame(mm@rowRanges)
peak_gr <- GRanges(seqnames = peak_gr$seqnames,
                   ranges = IRanges(start = peak_gr$start,
                                    end = peak_gr$end))

segment_df <- data.frame(
  chr = acd_wt$chr,
  start = acd_wt$start,
  end =acd_wt$end
)

peak_matrix <-test
segment_gr <- GRanges(
  seqnames = segment_df$chr,
  ranges = IRanges(start = segment_df$start, end = segment_df$end)
)

overlaps <- findOverlaps(segment_gr, peak_gr)
library(data.table)

dt <- data.table(
  seg_id = queryHits(overlaps),
  peak_id = subjectHits(overlaps)
)

dt <- cbind(dt, as.data.table(as.matrix(peak_matrix[dt$peak_id, ])))
result_matrix <- dt[, lapply(.SD, sum), by = seg_id, .SDcols = colnames(peak_matrix)]
final_matrix <- matrix(0, nrow = nrow(segment_df), ncol = ncol(peak_matrix))
final_matrix[result_matrix$seg_id, ] <- as.matrix(result_matrix[, -1])

rownames(final_matrix) <- paste(segment_df$chr, segment_df$start, segment_df$end, sep = "_")
colnames(final_matrix) <- colnames(peak_matrix)
WT_result_matrix<-wt_result_matrix

library(Seurat)
WT_result_matrix_n<-NormalizeData(WT_result_matrix)
WT_ACD_cor<-cor(t(WT_result_matrix_n),method = "pearson")
WT_ACD_cor[is.na(WT_ACD_cor)]<-0

#save(KO_result_matrix,WT_result_matrix,KO_ACD_cor,WT_ACD_cor,Var_ACD_cor,file = "acd_result_x.rdata")
#BiocManager::install("pheatmap")
library(pheatmap)
library(circlize)
library(RColorBrewer)

dir.create(paste0("coaccess_celltype/",Anno[i]))

WT_ACD_cor_adjusted<-WT_ACD_cor
data_trimmed <- WT_ACD_cor[upper.tri(WT_ACD_cor)]
q1 <- as.numeric(quantile(data_trimmed, 0.001))
q99 <- as.numeric(quantile(data_trimmed, 0.999))
WT_ACD_cor_adjusted[WT_ACD_cor_adjusted>q99]<-q99
WT_ACD_cor_adjusted[WT_ACD_cor_adjusted<q1]<-q1
min_val <- q1
max_val <- q99
total_colors <- 100  
n_blue_white <- round((0 - min_val) / (max_val - min_val) * total_colors)
n_white_red <- total_colors - n_blue_white
color1 <- colorRampPalette(c("#4682B4", "white"))(n_blue_white)
color2 <- colorRampPalette(c("white", "#B22222"))(n_white_red + 1)[-1]  
colors <- c(color1, color2)
breaks <- seq(min_val, max_val, length.out = total_colors + 1)

# plot Heatmap---------------
p<-pheatmap(WT_ACD_cor_adjusted,
         color = colors,
         breaks = breaks,
         cluster_rows = FALSE,   
         cluster_cols = FALSE,
         show_rownames = FALSE,  
         show_colnames = FALSE,
         legend = TRUE)  
pdf(paste0("coaccess_celltype/",Anno[i],"/wt_acd_cor.pdf"),w=10,h=10)
print(p)       
dev.off()


## -------- ACD co-access distance ---------
for (j in 1:19) {
  tmp<-KO_ACD_cor[grep(paste0("chr",j,"_"),rownames(KO_ACD_cor)),grep(paste0("chr",j,"_"),rownames(KO_ACD_cor))]
  rownames(tmp)<-(colsplit(rownames(tmp),"_",c(1,2,3))$`2`-1)/1000000
  colnames(tmp)<-rownames(tmp)
  result <- outer(as.numeric(colnames(tmp)), as.numeric(colnames(tmp)), function(x, y) abs(x - y))
  unique_d<-unique(as.vector(result))
  data<-data.table(pos=unique_d,ave=NA)
  for (k in 1:length(unique_d)) {
    data$ave[k] <- mean(KO_ACD_cor[which(result == unique_d[k], arr.ind = TRUE)])
  }
  dataall<-rbind(dataall,data)
}
dataall$group<-"KO"
dataall<-rbind(dataall,dataallwt)
dataall<-dataall[!dataall$ave==1,]
p<-ggplot(dataall, aes(x = pos, y = ave, color = group, fill = group)) +
  geom_point(size = 0.1, shape = 21, alpha = 0.9) +  
  geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "Test", x = "Distance", y = "Average Correlation") +
  theme_minimal() +  
  theme(panel.background = element_blank()) +
  scale_color_manual(values = c("WT" = "#1f77b4", "KO" = "#ff7f0e")) + 
  scale_fill_manual(values = c("WT" = "#1f77b4", "KO" = "#ff7f0e"))
pdf(paste0("coaccess_celltype/",Anno[i],"/cor_distance.pdf"),w=10,h=10)
print(p)
dev.off()
}


