library(Seurat)
library(readxl)
library(reshape2)
setwd("J:\\WXY\\LMNA_species")
load("E:/MCA2/mca_rawdata_process.rdata")
mca<-data
load("E:/MCA2/mca_cellinfo.rdata")
mca<-mca[,row.names(cell)]
cell<-cell[colnames(mca),]
mca$anno<-cell$cell_type
mca$stage<-mca$Stage
mca$lineage<-cell$cell_lineage


load("E:/DMEL/dge/rawdmel.rdata")
dca<-data
dca$lineage<-dca$anno2
unique(dca$lineage)


load("E:/Xenopus/Xenopus _DGE/Xenopus_sct.rdata")
xca<-mergeuse
xca$stage<-xca$tissue
xca$stage[!grepl("^St", xca$stage)] <- "Adult"
xca$stage[xca$stage=="Stomach"]<-"Adult"
xca$lineage<-xca$anno2

load("D:/ZCA-whole/fish_5stage_20210926.rdata")
load("D:/ZCA-whole/zca_cellinfo.rdata")
rownames(fishcellinfo)<-fishcellinfo$barcodes
pbmc<-pbmc[,rownames(fishcellinfo)]
fishcellinfo<-fishcellinfo[colnames(pbmc),]
pbmc$anno<-fishcellinfo$cell_type
zca<-pbmc
zca$lineage<-fishcellinfo$cell_lineage

load("D:/HCLpseudocell/HCL_rawdata.Rdata")
data_pbmc<-UpdateSeuratObject(data_pbmc)
data_pbmc$stage<-"Adult"
data_pbmc$stage[grep("Fetal",data_pbmc$tissue)]<-"Fetal"
dge<-data_pbmc@assays$RNA@counts
colnames(dge)<-colsplit(colnames(data_pbmc),"-",c(1,2))$`1`
hcl<-CreateSeuratObject(dge)
hcl$stage<-data_pbmc@meta.data$stage
hcl$tissue<-data_pbmc@meta.data$tissue
data <- read_excel("D:/HCLpseudocell/HCL_Fig1_cell_Info.xlsx")
rownames(data)<-data$cellnames
setdiff(rownames(data),colnames(hcl))
data<-data[colnames(hcl),]
hcl$anno<-data$celltype
unique(data$celltype)

classify_lineage <- function(celltype) {
  if (grepl("stromal|fibroblast|mesenchymal", celltype, ignore.case = TRUE)) {
    return("Stromal")
  } else if (grepl("macrophage|monocyte|dendritic|B cell|T cell|mast|neutrophil|antigen presenting|proliferating T cell", celltype, ignore.case = TRUE)) {
    return("Immune")
  } else if (grepl("muscle|cardiomyocyte|skeletal muscle", celltype, ignore.case = TRUE)) {
    return("Muscle")
  } else if (grepl("endocrine|goblet|chief|pancreas|thyroid", celltype, ignore.case = TRUE)) {
    return("Secretory")
  } else if (grepl("endothelial", celltype, ignore.case = TRUE)) {
    return("Endothelial")
  } else if (grepl("erythroid", celltype, ignore.case = TRUE)) {
    return("Erythroid")
  } else if (grepl("neuron|astrocyte|oligodendrocyte", celltype, ignore.case = TRUE)) {
    return("Neuron")
  } else if (grepl("epithelial|hepatocyte|AT2|basal|stratified", celltype, ignore.case = TRUE)) {
    return("Epithelial")
  } else if (grepl("germ|sertoli", celltype, ignore.case = TRUE)) {
    return("Germ")
  } else if (grepl("proliferating", celltype, ignore.case = TRUE)) {
    return("Proliferating")
  } else {
    return("Other")
  }
}

data$lineage <- sapply(data$celltype, classify_lineage)
data$lineage<-as.character(data$lineage)
rownames(data)<-data$cellnames
hcl$lineage<-data$lineage

standard_lineage <- c(
  "Immune", "Stromal", "Muscle", "Erythroid", "Endothelial",
  "Epithelial", "Neuron", "Secretory", "Germ", "Proliferating", "Other"
)

xca$lineage <- as.character(xca$lineage) 
xca$lineage[xca$lineage == "Germline"] <- "Germ"     
xca$lineage[xca$lineage == "Hepatocyte"] <- "Epithelial"  
xca$lineage <- factor(xca$lineage, levels = standard_lineage)  

dca$lineage <- as.character(dca$lineage)
dca$lineage[dca$lineage %in% c("Rp_high", "Fatbody")] <- "Other"  
dca$lineage[dca$lineage == "Follicle"] <- "Stromal"  
dca$lineage[dca$lineage == "MAG"] <- "Secretory"   
dca$lineage <- factor(dca$lineage, levels = standard_lineage)

hcl$lineage <- factor(hcl$lineage, levels = standard_lineage)
mca$lineage <- factor(mca$lineage, levels = standard_lineage)
zca$lineage <- factor(zca$lineage, levels = standard_lineage)

cat("\nCheck unique lineage categories after standardization:")
cat("\nHCL:", toString(levels(hcl$lineage)), 
    "\nXCA:", toString(levels(xca$lineage)),
    "\nZCA:", toString(levels(zca$lineage)), 
    "\nDCA:", toString(levels(dca$lineage)),
    "\nMCA:", toString(levels(mca$lineage)))


mca<-NormalizeData(mca)
zca<-NormalizeData(zca)
dca<-NormalizeData(dca)
hcl<-NormalizeData(hcl)
xca<-NormalizeData(xca)



hcl$stage2<-"Adult"
hcl$stage2[hcl$stage=="Fetal"]<-"Fetal"

#mca<-mca[,!mca$stage=="StemCell"]
table(mca$stage2)
mca$stage2<-"Adult"
mca$stage2[mca$stage%in%c("E10.5","E12.5","Fetal")]<-"Fetal"
mca$stage2[mca$stage%in%c("Neonatal","TenDays","ThreeWeeks")]<-"Development"

xca$stage2<-"Tadpole"
xca$stage2[xca$stage=="Adult"]<-"Adult"

zca$stage2<-"Adult"
zca$stage2[zca$stage%in%c("24hpf")]<-"Pharyngula"
zca$stage2[zca$stage%in%c("21Day","72hpf")]<-"Larval"

table(mca$stage)
mca$stage[mca$stage=="Fetal"]<-"E14.5"
mca$stage[mca$stage=="6-8Weeks"]<-"SixToEightWeeks"

###################################time point
time_orders <- list(
  mca = c("E10.5", "E12.5", "E14.5", "Neonatal", "TenDays", "ThreeWeeks", "SixToEightWeeks", "OneYear", "EighteenMonths", "TwoYears"),
  zca = c("24hpf", "72hpf", "21Day", "3Month", "22Month"),
  xca = c("St48", "St54", "St59", "St66", "Adult"),
  dca = unique(dca$stage)
)


plot_avg_expression <- function(
    obj, 
    time_order, 
    title = "Average Expression by Stage",
    palette = "RdBu",angle = 45,     
    hjust = 1
) {
  Idents(obj) <- obj$stage
  
  expr <- as.data.frame(AverageExpression(
    object = obj,
    features = c("lmnb2", "Lmna", "lmnb2.L", "LMNB1", "Lam"),
    assays = "RNA"  
  )$RNA)
  colnames(expr) <- gsub("^g", "", colnames(expr))
  expr_data <- data.frame(
    data = unlist(t(expr)),  
    group = factor(colnames(expr), levels = time_order)
  )
  
  ggplot(expr_data, aes(x = group, y = data, fill = group)) +
    geom_col(width = 0.7) +
    labs(
      title = title,
      x = "Stage",
      y = "Average Expression Level",
      fill = "Stage"
    ) +
    scale_fill_brewer(palette = palette) +
    theme_classic() +
    theme(
      axis.text.x = element_text(
        angle = angle,    
        hjust = hjust,     
        vjust = 1       
      )  
    )
}



pdf("mca_lmna_stage.pdf",w=6,h=4)
plot_avg_expression(
  obj = mca,
  time_order = time_orders$mca,
  title = "Mouse Lamin Expression"
)
dev.off()



###############heatmap

features = c("lmnb1","Lmnb1","lmnb1.L","LMNB1","Lam")
features =c("lmnb2","Lmnb2","lmnb2.L","LMNB2","Lam")
features = c("lmna","Lmna","lmna.S","LMNA","LamC")

obj<-hcl[,hcl$stage2=="Adult"]
Idents(obj)<-obj$lineage
expr<-as.data.frame(AverageExpression(obj,features = features)$RNA)
expr_data<-data.frame(data=t(expr),group=colnames(expr))
colnames(expr_data)<-c("data","group")
expr_data_all<-expr_data
expr_data_all$sp<-"human"

obj<-mca[,mca$stage2=="Adult"]
Idents(obj)<-obj$lineage
expr<-as.data.frame(AverageExpression(obj,features = features)$RNA)
expr_data<-data.frame(data=t(expr),group=colnames(expr))
colnames(expr_data)<-c("data","group")
expr_data$sp<-"mouse"
expr_data_all<-rbind(expr_data_all,expr_data)

obj<-zca[,zca$stage2=="Adult"]
Idents(obj)<-obj$lineage
expr<-as.data.frame(AverageExpression(obj,features = features)$RNA)
expr_data<-data.frame(data=t(expr),group=colnames(expr))
colnames(expr_data)<-c("data","group")
expr_data$sp<-"zebrafish"
expr_data_all<-rbind(expr_data_all,expr_data)

obj<-xca[,xca$stage2=="Adult"]
Idents(obj)<-obj$lineage
expr<-as.data.frame(AverageExpression(obj,features = features)$RNA)
expr_data<-data.frame(data=t(expr),group=colnames(expr))
colnames(expr_data)<-c("data","group")
expr_data$sp<-"xenopus"
expr_data_all<-rbind(expr_data_all,expr_data)

obj<-dca
Idents(obj)<-obj$lineage
expr<-as.data.frame(AverageExpression(obj,features = features)$RNA)
expr_data<-data.frame(data=t(expr),group=colnames(expr))
colnames(expr_data)<-c("data","group")
expr_data$sp<-"drosophila"
expr_data_all<-rbind(expr_data_all,expr_data)

short_data <- dcast(expr_data_all, group ~ sp, value.var = "data")
rownames(short_data)<-short_data$group
short_data<-short_data[,-1]
short_data[is.na(short_data)]<-0
unique(colnames(short_data))
short_data<-short_data[,c("drosophila","zebrafish","xenopus","mouse","human")]
library(pheatmap)
short_data[short_data>3]<-3
pdf("heatmap_lmnb2(L).pdf",w=6,h=6)
pheatmap(short_data, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         color = colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(100)[50:1])
dev.off()

