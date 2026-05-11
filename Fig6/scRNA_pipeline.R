library(Seurat)
library(data.table)
library(stringr)
library(reshape2)
library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)
library(RColorBrewer)
library(harmony)
library(DoubletFinder)
library(dplyr)
library(clusterProfiler)

tissue <- "Heart"
setwd(paste0("./scRNA/merge/Tissue/",tissue))
#------------create seob--------------

count = list()
file.ls <- c('./LMNA_result/')

sample=basename(getwd())
for (i in 1:length(file.ls)){
  message(i)
  count[[i]] = data.frame(fread(paste0(file.ls[i],'LMNA_WT.',sample,'.dge.txt.gz')),row.names = 1)
}
seob_list = list()
for (i in 1:2){
  message(i)
  seob_list[[i]] = CreateSeuratObject(counts = count[[i]],
                                      min.cells = 3, 
                                      min.features=200) 
  seob_list[[i]]$batch = paste0(i)
}

seob = merge(seob_list[[1]],seob_list[-1])
save(seob,file = 'seob_raw.rda')

#--------------merge----------------
id <- c('KO1','KO2','KO3','WT1','WT2','WT3')
seob_list = list()
for (i in id) {
  message(i)
  seob <- readRDS(paste0('./scRNA/merge/',i,'/tissue/',tissue,'/seob_passQC.rds'))
  seob_list[[i]] <- seob
  seob_list[[i]]$MouseID = paste0(i)
}

seob <- merge(seob_list[[1]],seob_list[-1])
seob@meta.data$type[seob@meta.data$MouseID %in% c("KO1", "KO2", "KO3")] <- "KO"
seob@meta.data$type[seob@meta.data$MouseID %in% c("WT1", "WT2", "WT3")] <- "WT"
table(seob$type)
saveRDS(seob,paste0(tissue,'_seob_passQC_merge.rds'))
seob <- readRDS(paste0(tissue,'_seob_passQC_merge.rds'))

#-------DoubletFinder------
Idents(seob) <- seob$MouseID
Group=as.character(unique(Idents(seob)))
Group

seob_list <- list()
for(i in 1:length(Group)){
  temp<-subset(seob,idents=Group[i])
  temp<-NormalizeData(temp)
  temp<-FindVariableFeatures(temp)
  temp<-ScaleData(temp)
  temp<-RunPCA(temp)
  temp<-FindNeighbors(temp,dims = 1:15)
  temp<-FindClusters(temp,resolution = 0.5)
  
  doublettest <- paramSweep(temp, PCs = 1:30, sct = FALSE)
  doublettest2 <- summarizeSweep(doublettest, GT = FALSE)
  doublettest3 <- find.pK(doublettest2)
  mpK<-as.numeric(as.vector(doublettest3$pK[which.max(doublettest3$BCmetric)]))
  temp$use<-Idents(temp)
  annotations <- temp$use
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.05*length(temp$orig.ident))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  final <- doubletFinder(temp, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  cellname.use<-rownames(final@meta.data)[final@meta.data[11]=="Singlet"]
  temp<-subset(final,cells=cellname.use)
  seob_list[[i]]=temp
}
seob <- merge(seob_list[[1]],seob_list[-1])
table(seob$MouseID)

#-----------scDblFinder---------
Idents(seob) <- seob$MouseID
Group=as.character(unique(Idents(seob)))
Group

# i=1
seob_list <- list()
for(i in 1:length(Group)){
  temp=subset(seob,idents=Group[i])
  temp=as.SingleCellExperiment(temp)
  temp=scDblFinder(temp, BPPARAM = MulticoreParam(1))
  temp=as.Seurat(temp)
  seob_list[[i]]=temp
}
seob <- merge(seob_list[[1]],seob_list[-1])
seob.dedoub <- subset(seob,scDblFinder.class == "singlet") 
dim(seob.dedoub)
saveRDS(seob.dedoub,'./Heart_seob_dedoub2.rds')

#-------------normalize-------------
seob <- NormalizeData(seob,normalization.method = "LogNormalize")
seob <- FindVariableFeatures(seob, 
                             selection.method = "vst", 
                             nfeatures = 3000)
seob <- ScaleData(seob)#split.by = 'batch'

#------------cluster------------
# seob <- combined
# seob <- ScaleData(seob)
seob <- RunPCA(seob)
ElbowPlot(seob, ndims = 50)

seob<-RunUMAP(seob,dims = 1:14,n.neighbors = 100,min.dist = 0.005)
seob<-FindNeighbors(seob,dims = 1:14)
seob<-FindClusters(seob,resolution = 0.4)

colors <- colorRampPalette(brewer.pal(8,"Paired"))(length(levels(as.factor(seob$seurat_clusters))))
pdf('umap_use_14_100_0.4.pdf',width = 10,height = 5)
p1 <- DimPlot(seob, reduction = "umap",  group.by = "type", 
        label = F)
p2 <- DimPlot(seob, reduction = "umap", 
        group.by  = "seurat_clusters",
        cols = colors,
        label = T)
p2+p1
dev.off()
save(seob, file = paste0(tissue,'_use_14_100_0.4_v4.rda'))
# check if dedup is reasonable?
# DimPlot(seob, reduction = "umap",  group.by = "scDblFinder.class", 
#         label = F)
# FeaturePlot(seob,features = 'scDblFinder.cxds_score')

# check batch effect
pdf(file = paste0(tissue,'_umap_26_batch.pdf'),width = 6.6,height = 6)
DimPlot(seob, reduction = "umap",  group.by = "MouseID", 
        label = F)
dev.off()

DimPlot(seob, reduction = "umap",  group.by = "type", 
        label = F)

seob <- FindNeighbors(seob, dims = 1:26)
seob <- FindClusters(seob, resolution = 0.3, # 值越大，cluster 越多
                     random.seed = 1) 

#pdf(file = 'plot/umap_raw_18_0.6.pdf',width = 6.6,height = 6)

colors <- colorRampPalette(brewer.pal(8,"Paired"))(length(levels(as.factor(seob$seurat_clusters))))
DimPlot(seob, reduction = "umap", 
        group.by  = "seurat_clusters",
        cols = colors,
        label = T)

pdf(file = paste0(tissue,'_umap_26_0.3.pdf'),width = 13,height = 6)
p1 <- DimPlot(seob, reduction = "umap", 
              group.by  = "seurat_clusters",
              cols = colors,
              label = T)

p2 <- DimPlot(seob, reduction = "umap", 
              group.by  = "type",
              label = T)
p1+p2
dev.off()
save(seob,file = 'seob_harmony_30_0.3.rda')

table(seob$seurat_clusters, seob$type)

markers <- FindAllMarkers(seob, only.pos = TRUE,logfc.threshold = 0.5)
markers <- markers[order(markers$cluster, -markers$avg_log2FC, markers$p_val),]
write.csv(markers,file=paste0(tissue,'_use_14_100_0.4_markers.csv'),quote=F)

#----------anno---------
library(Seurat)
load("./scRNA/merge/Tissue/Heart/umi1000.rdata")

cluster2type <- c("0"='CM_WT',
                  "1"="CM_WT",
                  "2"="CM_KO",
                  "3"="CM_WT",
                  "4"="CM_KO",
                  "5"="Fibroblast_WT",
                  "6"="Fibroblast_KO",
                  "7"="CM_WT",
                  "8"='CM_KO',
                  "9"="Endothelial cell",
                  "10"="Atrial CM",
                  "11"="CM_WT",
                  "12"="CM_WT",
                  "13"="Macrophage",
                  "14"="CM_KO",
                  "15"="SMC",
                  "16"='Proliferating cell',
                  "17"="Unknown")

seob[['cell_type']] = unname(cluster2type[seob@meta.data$seurat_clusters])
source('/media/ggj/ggj/CJY/Code/source/cell_prop_LDY.R')
prop = cell_prop_LDY(seob = seob ,cell_type = 'cell_type',
                     group ='type',position = 'fill')

colors <- colorRampPalette(brewer.pal(8,"Paired"))(length(levels(as.factor(seob$cell_type))))
colors <- c('#C3DDEB','#769BC9','#CBE3AF','#86B96B','#F6BAB8','#EA7250','#FAD19C','#F5A661')
colors <- c('#9999B8','#957dad','#D4B3D2','#8CB4B1','#019092','#B8DBDB','#B091B4','#D3D3E9','#E9D5BA','#EAC09A')

pdf(file = paste0(tissue,'_umap_anno_separate_14_100_0.4.pdf'),width = 16,height = 6)
p1=DimPlot(seob, 
           reduction = "umap", repel = T,
           cols = colors,
           group.by = "cell_type",
           label = TRUE, pt.size = 0.3) + NoLegend()

p2=DimPlot(seob, reduction = "umap", 
           group.by  = "type",
           label = T)

p1+p2+prop[['p']]
dev.off()
save(seob, file = "seob_anno_separate.rda")

#-------marker average expression heatmap-----
Idents(seob) <- seob$cell_type
pbmc.markers <- FindAllMarkers(object =seob, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold  = 0.5)
pbmc.markers<-pbmc.markers[order(pbmc.markers$cluster,-pbmc.markers$avg_log2FC,pbmc.markers$p_val_adj  ),]
top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
DoHeatmap(object = seob, features  = top10$gene,  label = T, size=3)
ggsave('./marker_heatmap.pdf', width = 10,height = 14)

marker <- markerlist_heart_2_$marker
marker <- strsplit(marker,',')
marker <- unlist(marker)
marker <- gsub(' ', '',marker)
marker <- marker[!marker %in% c( "Sox10" , "Ntm" ,   "Gria2")]


gene_cell_exp <- AverageExpression(seob,
                                   group.by = 'cell_type') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)

library(ComplexHeatmap)
df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c( "Atrial CM" = '#9999B8',"CM_I"='#957dad', 'CM_II'='#D4B3D2',
                                                   "Endothelial cell"= '#8CB4B1', "Fibroblast_I"= '#019092','Fibroblast_II'='#B8DBDB',
                                                 "Macrophage" = '#B091B4',"Proliferating cell"='#D3D3E9',
                                                   "SMC" = '#E9D5BA',"Unknown" = '#EAC09A')))

top5$cluster <- as.character(top5$cluster)
top5 <- top5[order(top5$cluster),]
marker <- top5$gene
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
marker_exp_top <- marker_exp[rownames(marker_exp) %in% marker,]
marker_order <- unique(marker)
pdf('./average_expression_marker_anno2.pdf',width =6,height = 6)  
Heatmap(marker_exp_top,
        cluster_rows = F,
        cluster_columns = F,
        row_order = marker_order,
        column_order = c("Atrial CM" ,"CM_I", 'CM_II',"Endothelial cell","Fibroblast_I" , "Fibroblast_II" ,  "Macrophage", "Proliferating cell", "SMC",   "Unknown" ),
        show_column_names = F,
        show_row_names = T,
        column_title = NULL,
        heatmap_legend_param = list(
          title=' '),
        col = colorRampPalette(c("white","#BFDAE9","#6CA6CD"))(100),
        border = 'black',
        rect_gp = gpar(col = "black", lwd = 0.2),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        top_annotation = top_anno)
dev.off()       



