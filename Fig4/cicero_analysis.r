setwd("/media/ggj/ggj/CJY/cuttag/scATAC")
library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)
addArchRThreads(threads = 1)

proj_callpeak$Anno_type<-paste0(proj_callpeak$celltype,"_",proj_callpeak$type)
table(proj_callpeak$Anno_type)

tissue_type<-unique(proj_callpeak$Anno_type)
table(proj_callpeak$Anno_type)


## ------------ get peak matrix -------------
i=1
for (i in 2) {
  use<-tissue_type[i]
  proj<-proj_callpeak[proj_callpeak$Anno_type==use]
  
  proj<-proj[sample(nrow(proj@cellColData),2000)]
  peak_matrix <- getMatrixFromProject(proj, useMatrix = 'PeakMatrix', binarize = T)
  saveRDS(peak_matrix,file = file.path(use,'pmat_2k.rds'))

  ## get peak matrix
  pmat <- peak_matrix@assays@data@listData[["PeakMatrix"]]
  dim(pmat)
  peakset <- as.data.frame(proj@peakSet, row.names = 1:length(proj@peakSet))
  gr <- peak_matrix@rowRanges
  gr.df <- as.data.frame(gr)
  dim(gr.df)
  rownames(pmat) <- paste(gr.df$seqnames,
                          gr.df$start,gr.df$end,sep = '_')

  pmat_sub <- pmat[rowSums(pmat)>0,]
  dim(pmat_sub) # 5226903    1000
  gr.df$id <- paste(gr.df$seqnames,
                    gr.df$start,gr.df$end,sep = '_')
  gr.df.sub <- gr.df[gr.df$id %in% rownames(pmat_sub),]

  ## save data
  class(pmat_sub) #"dgCMatrix"
  # writeMM(pmat_sub,file = file.path(outDir,'packed','matrix.mtx'))
  saveRDS(pmat_sub,file = file.path(use,'matrix.rds'))

  bc.df <- data.frame(cells = colnames(pmat_sub))
  write.csv(bc.df, file = file.path(use,'barcodes.csv'), quote = F,
            row.names = F)

  peak.df <- data.frame(chr = gr.df.sub$seqnames, bp1 = gr.df.sub$start, bp2 = gr.df.sub$end)
  write.csv(peak.df, file = file.path(use,'peaks.csv'),quote = F, row.names = F)

}


## ------------ cicero analysis -------------
library(cicero)
library(monocle3)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(reshape2)

output_folder <- "cicero_output"
dir.create(output_folder, showWarnings = FALSE)

indata <- readRDS('matrix_2k.rds')
indata@x[indata@x > 0] <- 1

cellinfo <- read.csv('barcodes_2k.csv', header = TRUE)
row.names(cellinfo) <- cellinfo$cells

peakinfo <- read.csv('peaks_2k.csv', header = TRUE)
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
row.names(peakinfo) <- peakinfo$site_name

row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)

# Make CDS
input_cds <-  suppressWarnings(new_cell_data_set(indata,
                                                 cell_metadata = cellinfo,
                                                 gene_metadata = peakinfo))
input_cds <- monocle3::detect_genes(input_cds)
# Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]

# 4. Process cicero-CDS object
# Data preprocessing
set.seed(666)
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI", num_dim = 10)

# Dimensional reduction with umap
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP',
                              preprocess_method = "LSI")
umap_coords <- reducedDims(input_cds)$UMAP
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

saveRDS(cicero_cds, file = file.path(output_folder,'cicero_cds_2k.rds'))
## reload
# cicero_cds <- readRDS(file.path(output_folder,'cicero_cds.rds'))

chromosome_length <- read.table("/media/ggj/ggj/CJY/tools/mm10/mm10.main.chrom.sort.sizes")

conns <- run_cicero(cicero_cds, chromosome_length) # Takes a few minutes to run

# Check results
head(conns)
all_peaks <- row.names(exprs(input_cds))
# write.csv(x = all_peaks, file = paste0(output_folder, "/all_peaks.csv"))
# write.csv(x = conns, file = paste0(output_folder, "/cicero_connections.csv"))
saveRDS(all_peaks,file = paste0(output_folder, "/all_peaks_2k.rds"))
saveRDS(conns,file = paste0(output_folder, "/cicero_connections_2k.rds"))
conns_filter <- conns[conns$coaccess>0.25,]
saveRDS(conns_filter,file = paste0(output_folder, "/cicero_connections_coacess0.25_2k.rds"))

conns_filter<-readRDS(paste0(output_folder, "/cicero_connections_coacess0.25_2k.rds"))


conns_mouse <- na.omit(conns_filter)
write.csv(conns_mouse,"results_raw_2k.csv")

## ------------- EP loop --------------
library(reshape2)
library(GenomicRanges)

gene_anno_mouse <- rtracklayer::readGFF('/media/ggj/ggj/CJY/tools/mm10/Ensembl_release98/Mus_musculus.GRCm38.98.gtf.gz')
gene_anno_mouse$chromosome <- paste0("chr", gene_anno_mouse$seqid)
gene_anno_mouse$gene <- gene_anno_mouse$gene_id
gene_anno_mouse$transcript <- gene_anno_mouse$transcript_id
gene_anno_mouse$symbol <- gene_anno_mouse$gene_name

## get promoter
gene_mouse <- gene_anno_mouse[gene_anno_mouse$type == 'gene',]
gene_mouse.gr <- GRanges(gene_mouse$chromosome, IRanges(gene_mouse$start, gene_mouse$end), strand = gene_mouse$strand, genes = gene_mouse$symbol)
promoter_mouse.gr <- promoters(gene_mouse.gr,upstream = 2000,downstream = 200) # default: [-2000,+200]

extract_regulatory_loop <- function(promoter_gr = NULL, cicero_conns = NULL){
  peak1.df <- colsplit(cicero_conns$Peak1,pattern = '_',names = c('chr','start','end'))
  peak1.gr <- GRanges(peak1.df$chr,IRanges(peak1.df$start,peak1.df$end))
  peak2.df <- colsplit(cicero_conns$Peak2,pattern = '_',names = c('chr','start','end'))
  peak2.gr <- GRanges(peak2.df$chr,IRanges(peak2.df$start,peak2.df$end))
  
  peak1.gr$olp <- 0
  idx1 <- queryHits(findOverlaps(peak1.gr,promoter_gr))
  peak1.gr$olp[idx1] <- 1
  peak2.gr$olp <- 0
  idx2 <- queryHits(findOverlaps(peak2.gr,promoter_gr))
  peak2.gr$olp[idx2] <- 1
  
  cicero_conns$Peak1_olp <- peak1.gr$olp
  cicero_conns$Peak2_olp <- peak2.gr$olp
  cicero_conns$Peak_olp <- cicero_conns$Peak1_olp + cicero_conns$Peak2_olp
  
  conns_regulatory <- cicero_conns[cicero_conns$Peak_olp>0,]
  return(conns_regulatory)
}
conns_mouse<-read.csv("results_raw_2k.csv",row.names = 1)

conns_EP_mouse <- extract_regulatory_loop(promoter_gr = promoter_mouse.gr, cicero_conns = conns_mouse)
library(rtracklayer)
library(GenomicRanges)

gtf <- import('/media/ggj/ggj/CJY/tools/mm10/Ensembl_release98/Mus_musculus.GRCm38.98.gtf.gz')
genes <- gtf[gtf$type == "gene"]
start(genes) <- pmax(1, start(genes) - 2000)
end(genes) <- end(genes) + 200
conns_EP_mouse_P1<-colsplit(conns_EP_mouse$Peak1,"_",c(1,2,3))
conns_EP_mouse_P1$`1`<-gsub("chr","",conns_EP_mouse_P1$`1`)
gr_pos1 <- GRanges(seqnames = conns_EP_mouse_P1$`1`, ranges = IRanges(start = conns_EP_mouse_P1$`2`, end = conns_EP_mouse_P1$`3`))
conns_EP_mouse_P2<-colsplit(conns_EP_mouse$Peak2,"_",c(1,2,3))
conns_EP_mouse_P2$`1`<-gsub("chr","",conns_EP_mouse_P2$`1`)
gr_pos2 <- GRanges(seqnames = conns_EP_mouse_P2$`1`, ranges = IRanges(start = conns_EP_mouse_P2$`2`, end = conns_EP_mouse_P2$`3`))

overlaps_pos1 <- findOverlaps(gr_pos1, genes)
conns_EP_mouse$pos1_gene <- NA
conns_EP_mouse$pos1_gene[queryHits(overlaps_pos1)] <- genes$gene_name[subjectHits(overlaps_pos1)]


overlaps_pos2 <- findOverlaps(gr_pos2, genes)
conns_EP_mouse$pos2_gene <- NA
conns_EP_mouse$pos2_gene[queryHits(overlaps_pos2)] <- genes$gene_name[subjectHits(overlaps_pos2)]
write.csv(conns_EP_mouse,"results.csv")


## ------------- plot cicero pair ------------
gene_anno <- rtracklayer::readGFF('/media/ggj/ggj/CJY/tools/mm10/Mus_musculus.GRCm38.88.gtf')
# rename some columns to match requirements
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

df_filtered <- allcell_WT_results_processed_tri[grepl("Nfia", allcell_WT_results_processed_tri$gene1) | grepl("Nfia", allcell_WT_results_processed_tri$gene2), ]
df_filtered <- allcell_KO_results_processed_tri[grepl("Klf9", allcell_KO_results_processed_tri$gene1) | grepl("Klf9", allcell_KO_results_processed_tri$gene2), ]

df_use <- df_filtered[,c("Peak1","Peak2","coaccess")]

pdf('./KO.pdf',width = 6,height = 4)
plot_connections(df_use, "chr19", 22500000, 23200000,
                 gene_model = gene_anno,
                 coaccess_cutoff = .25,
                 connection_width = .5,
                 collapseTranscripts = "longest")
dev.off()
