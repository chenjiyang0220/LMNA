
## ---------- KO loss copair GO function -------------
data_KO <- fread("./scATAC/copair/All_KO_results_processed.csv")
data_WT <- fread("./scATAC/copair/All_WT_results_processed.csv")
data_promoter <- fread("/media/ggj/ggj/CJY/cuttag/scATAC/coAccessibility/allcell_copair_promoter_count.csv")
data_promoter[, WT_KO_diff := as.numeric(WT_count) - as.numeric(KO_count)]
setorder(data_promoter, -WT_KO_diff)    
Top10gene <- unique(data_promoter[1:10, gene])
Top300gene <- unique(data_promoter[1:300, gene])

# Top10 promoter gene selection
WT_top10_promoter <- data_WT[
  (peak1_typeuse == "Promoter" & gene1 %in% Top10gene) |
  (peak2_typeuse == "Promoter" & gene2 %in% Top10gene)
]
WT_top10_promoter[, Toptype := NA_character_]
WT_top10_promoter[
  peak1_typeuse == "Promoter" & gene1 %in% Top10gene,
  Toptype := gene1
]
WT_top10_promoter[
  is.na(Toptype) & peak2_typeuse == "Promoter" & gene2 %in% Top10gene,
  Toptype := gene2
]

KO_top10_promoter <- data_KO[
  (peak1_typeuse == "Promoter" & gene1 %in% Top10gene) |
  (peak2_typeuse == "Promoter" & gene2 %in% Top10gene)
]
KO_top10_promoter[, Toptype := NA_character_]
KO_top10_promoter[
  peak1_typeuse == "Promoter" & gene1 %in% Top10gene,
  Toptype := gene1
]
KO_top10_promoter[
  is.na(Toptype) & peak2_typeuse == "Promoter" & gene2 %in% Top10gene,
  Toptype := gene2
]

# Top300 gene GO enrichment
library(clusterProfiler)
library(org.Mm.eg.db)
library(forcats) 
library(scales)

go_database <- "org.Mm.eg.db"
top300_gene_id <- bitr(
  Top300gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = go_database
)

Top300_GO <- enrichGO(
  gene = unique(top300_gene_id$ENTREZID),
  OrgDb = go_database,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

#### tree plot
library(enrichplot)
Top300_GO_simplified <- simplify(
  Top300_GO,
  cutoff = 0.67,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang"
)

Top300_GO_simplified <- pairwise_termsim(Top300_GO_simplified)
p_top300_tree <- treeplot(
  Top300_GO_simplified,
  nCluster = 6,
  showCategory = 12
) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8)
  ) +
  ggtitle("Top300 gene GO semantic clustering")
ggsave(
  plot = p_top300_tree,
  filename = "./scATAC/copair/WT-KO_diff_Top300gene_GO_treeplot.pdf",
  width = 13,
  height = 6
)


# ------ candidate Enhancer validation --------
co_WT<- fread("./scATAC/copair/All_WT_results_TE.csv")
co_KO<- fread("./scATAC/copair/All_KO_results_TE.csv")
co_WT <- co_WT[,-1]
co_KO <- co_KO[,-1]
Peak_deg<- fread("./scATAC/copair/Promoter_wt_ko_DEG.csv")
WT_id <- Peak_deg %>%          
  filter(type == "WT") %>%    
  pull(peakid) %>%            
  unique()               
KO_id<- Peak_deg %>%        
  filter(type == "KO") %>%    
  pull(peakid) %>%            
  unique()   
co_WT_deg <- co_WT%>%filter(Peak1 %in%WT_id | Peak2 %in% WT_id)%>% filter(EP=="EP")
co_KO_deg <- co_KO%>%filter(Peak1 %in%KO_id | Peak2 %in% KO_id)%>% filter(EP=="EP")

WT_df <- bind_rows(
  co_WT_deg %>% filter(peak1_typeuse == "Enhancer") %>%
    transmute(type = "WT",
              peak_id   = Peak1,
              TE = Peak1_TE),
  co_WT_deg %>% filter(peak2_typeuse == "Enhancer") %>%
    transmute(type = "WT",
              peak_id   = Peak2,
              TE = Peak2_TE))%>% distinct()  
KO_df <- bind_rows(
  co_KO_deg %>% filter(peak1_typeuse == "Enhancer") %>%
    transmute(type = "KO",
              peak_id   = Peak1,
              TE = Peak1_TE),
  co_KO_deg %>% filter(peak2_typeuse == "Enhancer") %>%
    transmute(type = "KO",
              peak_id   = Peak2,
              TE = Peak2_TE))%>% distinct()  
df<-data.frame()
df<-rbind(WT_df,KO_df)
write.csv(df,"./scATAC/copair/co_wt_ko_nonP_deg.csv")

k27_wt <- fread("./bulk_cut/H3K27ac/WT_H3K27ac_peaks_20.narrowPeak")
k27_ko <- fread("./bulk_cut/H3K27ac/KO_H3K27ac_peaks_20.narrowPeak")
k27.wt.gr<- GRanges(seqnames = k27_wt$V1, IRanges(start = k27_wt$V2, end = k27_wt$V3))
k27.ko.gr<- GRanges(seqnames = k27_ko$V1, IRanges(start = k27_ko$V2, end = k27_ko$V3))
df <- fread("./scATAC/copair/co_wt_ko_nonP_deg.csv")
df <- df %>%
  separate(peak_id,
           into      = c("seqnames", "start", "end"),
           sep       = "_",          
           remove    = FALSE,      
           convert   = TRUE)       
wt<- df%>%filter(type=="WT")
ko<- df%>%filter(type=="KO")
wt.gr<- GRanges(seqnames = wt$seqnames,IRanges(start = wt$start,end = wt$end),peakid=wt$peak_id, TE=wt$TE)
ko.gr<- GRanges(seqnames = ko$seqnames,IRanges(start = ko$start,end = ko$end),peakid=ko$peak_id, TE=ko$TE)
ol.wt <- findOverlaps(k27.wt.gr,wt.gr)
ol.ko <- findOverlaps(k27.ko.gr,ko.gr)
wt.gr$k27 <- 0
ko.gr$k27 <- 0
wt.gr$k27[unique(subjectHits(ol.wt))] <- 1
ko.gr$k27[unique(subjectHits(ol.ko))] <- 1
wt_df<-data.frame(wt.gr)
ko_df<-data.frame(ko.gr)
wt_k27<- wt_df%>%filter(k27==1)
ko_k27 <- ko_df%>%filter(k27==1)
write_tsv(wt_k27,"./scATAC/copair/H3K27ac/wt_co_k27ac_deg.bed", col_names = FALSE, quote = "none")
write_tsv(ko_k27,"./scATAC/copair/H3K27ac/ko_co_k27ac_deg.bed", col_names = FALSE, quote = "none")

## --- turn to HOMER motif analysis ---
