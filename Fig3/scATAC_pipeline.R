setwd("./scATAC/")
library(ArchR)
addArchRThreads(threads = 6)
addArchRGenome("mm10")
bed <- list.files(path = "./bed", pattern = "bed",full.names = T)
bed<-bed[-grep("tbi",bed)]
tissue<-gsub("./bed/","",bed)
tissue<-gsub(".sort.bed.gz","",tissue)

## ------------create ArrowFiles-------------
ArrowFiles <- createArrowFiles(
  inputFiles = bed,
  sampleNames = tissue,
  minTSS = 7,
  minFrags = 1000,
  threads =6, 
  addTileMat = TRUE,
  TileMatParams = list(tileSize = 500),
  addGeneScoreMat = TRUE)
ArrowFiles<-list.files(pattern = "arrow")
doubScore <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1,force = T
)

proj1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "Output",
  copyArrows = FALSE, #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

df <- getCellColData(ArchRProj = proj1, select = c('nFrags','TSSEnrichment'))
df$`log10(nFrags)` <- log10(df$nFrags)
dim(df) #41947
p <- ggPoint(
  x = df[,3], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight") +
  geom_hline(yintercept = 16, lty = "dashed") + 
  geom_vline(xintercept = 3.301, lty = "dashed")
p

UQ_cutoff = 2000
TSS_cutoff = 10

df_sub <- df[df$nFrags>UQ_cutoff & df$TSSEnrichment>TSS_cutoff,]
dim(df_sub)


cellsPass <- rownames(df_sub)
proj1 <- proj1[cellsPass,]

## ------------ Filter doublets ------------ 
proj1 <- filterDoublets(proj1)
proj1

# save
saveRDS(proj1,file = 'proj_TSS10_UQ2000_doubletsFiltered.rds')

## -------------- clustering -------------
proj2 <- addIterativeLSI(
  ArchRProj = proj1,
  useMatrix = 'TileMatrix',
  name = 'IterativeLSI',
  iterations = 2,
  clusterParams = list(
    resolution = c(1), # set res of LSI: 1/2
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 20000, # set varFeatures
  dimsToUse = 1:20, # set dims to use
  #filterQuantile = 0.995, # filter features
  # outlierQuantiles = c(0.05,0.95), # filter cells
  force = TRUE
)

proj2 <- addClusters(
  input = proj2,
  reducedDims = 'IterativeLSI',
  method = 'Seurat',
  name = 'Clusters',
  resolution = 0.5, # set res
  maxClusters = NULL, 
  #nOutlier = 50, # set the minimal number of cells
  force = TRUE
)

proj2 <- addUMAP(
  ArchRProj = proj2,
  reducedDims = 'IterativeLSI',
  name = 'UMAP',
  nNeighbors = 60,
  minDist = 0.2,
  metric = 'cosine',
  force = TRUE
)

plotEmbedding(ArchRProj = proj2, colorBy = 'cellColData', name = 'Clusters', embedding = 'UMAP')

markersGS <- getMarkerFeatures(
  ArchRProj = projHeme3, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList<-as.data.frame(markerList)
write.csv(markerList,file = "marker2_res0.5.csv")

save(projHeme3,file = "umap_v1.rdata")

## ------------- peak calling -------------
proj_callpeak@cellColData$Anno_type<-paste0(proj_callpeak@cellColData$Anno,"_",proj_callpeak@cellColData$type)
proj_callpeak <- addGroupCoverages(ArchRProj = proj_callpeak, groupBy = "Anno_type")
pathToMacs2 <- findMacs2()
proj_callpeak <- addReproduciblePeakSet(
  ArchRProj = proj_callpeak, 
  groupBy = "Anno_type", 
  pathToMacs2 = pathToMacs2
)
proj_callpeak <- addPeakMatrix(proj_callpeak)

tissue<-unique(proj_callpeak$Anno)
write.csv(proj_callpeak@peakSet,file = "peakset_all.csv")
saveRDS(proj_callpeak,file = "proj_callpeak_all.rds")


## --------------- motif annotation -----------------
proj_callpeak<-readRDS("proj_callpeak.rds")
metadata <- as.data.frame(proj_callpeak@cellColData)
metadata$cell_name <- rownames(metadata)
table(metadata$Anno)
metadata$Anno[metadata$Anno=="APC"]<-"Myeloid cell"
metadata$Anno[metadata$Anno=="Macrophage"]<-"Myeloid cell"
metadata$Annotype<-paste0(metadata$type,"_",metadata$Anno)
library(dplyr)
sampled_metadata <- metadata %>%
  group_by(Annotype) %>%
  slice_sample(n = 100) %>%
  ungroup()
sampled_cells <- sampled_metadata$cell_name

new_proj <- proj_callpeak[sampled_cells, ]

new_proj@projectMetadata$outputDirectory
proj_callpeak <- addMotifAnnotations(ArchRProj = proj_callpeak, motifSet = "homer", name = "Motif",force = TRUE)
proj_callpeak <- addBgdPeaks(proj_callpeak)
proj_callpeak <- addDeviationsMatrix(
  ArchRProj = proj_callpeak, 
  peakAnnotation = "Motif",
  force = TRUE
)

proj_callpeak$Anno[proj_callpeak$Anno=="APC"]<-"Myeloid cell"
proj_callpeak$Anno[proj_callpeak$Anno=="Macrophage"]<-"Myeloid cell"
proj_callpeak$Annotype<-paste0(proj_callpeak$type,"_",proj_callpeak$Anno)
seGroupMotif <- getGroupSE(ArchRProj = proj_callpeak, useMatrix = "MotifMatrix", groupBy = "Annotype")
motif <- as.matrix(assay(seGroupMotif))
anno <- as.data.frame(rowData(seGroupMotif))
anno <- anno[anno$seqnames=="z",]
rownames(motif) <- rowData(seGroupMotif)$name
library(viridis)

lineage_mapping <- list(
  Secretory = c("Acinar cell", "Pituitary cell", "Secretory cell"),
  Muscle = c("Cardiomyocyte", "Skeletal muscle cell", "Smooth muscle cell"),
  Immune = c("B cell", "Myeloid cell", "T cell"),
  Neural = c("Astrocyte", "Neural precursor cell", "Neural stem cell", 
             "Neuroendocrine cell", "Neuron", "Oligodendrocyte", 
             "Oligodendrocyte precursor cell", "Rod cell", "Schwann cell"),
  Endothelium = c("Endocardial endothelial cell", "Endothelial cell"),
  Epithelium = c("Epithelial cell", "Gastrointestinal epithelial cell", 
                 "Kenal epithelial cell", "Keratinocyte", "Urothelial cells", "Hepatocyte"),
  Stromal = c("Fibroblast", "Adiopcyte"),
  Erythroid = c("Erythroid cell")
)

label <- c('Nkx2.5','Tbx20','Gata4',
           'Mef2a', 'Mef2b', 'Mef2c', 'Mef2d','MyoD', 'MyoG', 'Pbx1', 'Pbx3',
           'Olig2', 'NeuroD1','NEUROG2','ASCL1','Sox10','PAX6',
           'CEBP', 'RUNX2', 'PPARγ',
           'Nr5a2', 'SF1', 'Six1', 'Six2', 'LHX3', 'MIST1','SF1', 'NEUROD1', 'ISL1',
           'RUNX1', 'RUNX3', 'IRF4', 'IRF8', 'IRF1',
           'GATA2', 'KLF3',
           'ELF3', 'FOXA1','FOXA2','FOXO1', 'GRHL2', 'KLF4',
           'ETS', 'Sox18', 'NRP1'
           )

label_motif <- rownames(motif)[grep(paste(label, collapse  = '|'),rownames(motif), ignore.case = T)]

data <- t(scale(t(motif)))

ht <- Heatmap(
  matrix = data,
  name = "value",  
  col = viridis(100),
  cluster_rows = T,
  cluster_columns = F,  
  show_row_names = F,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),  
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 10)
)

highlight_pos <- which(rownames(motif) %in% label_motif)
row_anno <- rowAnnotation(
  mark = anno_mark(at = highlight_pos,
                   labels = rownames(motif)[highlight_pos],
                   side = "right",
                   labels_gp = gpar(fontsize = 9))
)

pdf('heatmap.pdf',width = 7,height = 9)
ht+row_anno
dev.off()
