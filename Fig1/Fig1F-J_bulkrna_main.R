library("FactoMineR")
library("factoextra")
library(clusterProfiler)
library("enrichplot")
library(ggplot2)
library(ggpubr)
library(limma)
library(org.Mm.eg.db)
library(pheatmap)
library(devtools)
library(stringr)
library(reshape2)
library("ggfortify")
library(dplyr)

setwd("./bulkRNA")
options(stringsAsFactors = F)

############# load dataset#################
files<-list.files(pattern = "count",recursive = T)

filename<-colsplit(files,"/",c(1,2))
names<-colsplit(filename$`2`,"\\.",c(1,2,3))
names$file<-files
names$type<-names$`1`
names$type[grep("HE",names$type)]<-"HE"
names$type[grep("WT",names$type)]<-"WT"
names$type[grep("KO",names$type)]<-"KO"
names$`2`[grep("Adrenal",names$`2`)]<-"AdrenalGland"
names$`2`[grep("Brain",names$`2`)]<-"BrainStem"
names$`2`[grep("Mammary Gland",names$`2`)]<-"MammaryGland"
names$`2`[grep("Fat",names$`2`)]<-"WhiteFat"

names<-names[-grep("Cranium",names$`2`),]
sample<-c("WT30.AdranalGland", "KO201.AdranalGland", 
            "WT202.Bladder", 
            "WT6.BoneMarrow", "WT202.BoneMarrow", "HE4.BoneMarrow", 
            "WT30.BrainStem", "WT202.BrainStem", "WT25.BrainStem", 
            "KO29.Cerebellum", "WT6.Cerebellum", 
            "WT202.Cerebrum", "KO2.Cerebrum", "WT22.Cerebrum", 
            "KO2.Eye", "WT202.Eye", 
            "KO8.Gonad", "WT101.Gonad", 
            "HE15.Heart", "WT9.Heart", 
            "KO29.Kidney", "WT6.Kidney", 
            "KO201.LargeIntestine", "KO207.LargeIntestine", "KO1.LargeIntestine", 
            "KO201.Liver", 
            "WT1.Lung", "WT6.Lung", 
            "WT25.LymphNode", "KO29.LymphNode", 
            "WT25.MammaryGland", 
            "WT25.Muscle", 
            "WT9.Pancreas", "HE15.Pancreas", 
            "KO201.PituitaryGland", "WT21.PituitaryGland", 
            "WT25.Skin", "KO20.Skin", "KO201.Skin", 
            "WT25.SmallIntestine", "WT30.SmallIntestine", "WT21.SmallIntestine", "KO7.SmallIntestine", 
            "WT202.SpinalCord", "WT22.SpinalCord", "WT21.SpinalCord", 
            "WT23.Spleen", 
            "KO201.Stomach", "KO27.Stomach", 
            "KO27.Testicles", "WT205.Testicles", 
            "WT25.Thymus", 
            "KO1.Vascular", "KO27.Vascular", "WT25.Vascular", "HE1.Vascular", 
            "WT30.WhiteFat")
names$sample<-paste0(names$`1`,".",names$`2`)

names<-names[!names$sample%in%sample,]
write.csv(names,file = "names.csv")
table(names$`2`)

names_use<-paste0(names$`1`,"_",names$`2`)
files<-names$file
gset1 <- read.table(file = files[1],header = T)

for (i in 2:length(files)) {
  gset <- read.table(file = files[i],header = T)
  gset1[,i+6] <- gset$starAligned.sortedByCoord.out.bam
}


bb<-c()
gsetuse<-gset1[,7:(length(files)+6)]
i=1
for (i in 1:length(files)) {
    aa2<-as.data.frame(table(gsetuse[,i]))
    bb[i]<-26214-aa2$Freq[1]
}

gsetuse$sum<-rowSums(gsetuse)
###gene
summary(bb)
names$genenum<-bb
write.csv(names,file = "sample_data.csv")

dge<-gset1[,7:(length(files)+6)]
rownames(dge)<-gset1$Geneid
colnames(dge)<-names_use
meta<-names
library(Seurat)
pbmc<-CreateSeuratObject(dge)
pbmc$type<-names$type
pbmc$tissue<-names$`2`
pbmc$mouse<-names$`1`
save(pbmc,file = "dge.rdata")


gene.datExpr <- log2(gset1[,7:(length(files)+6)]+1)
colnames(gene.datExpr)<-names_use
gene.datExpr$ENSEMBL <- gset1$Geneid

gene.datExpr$median <- apply(gene.datExpr[,1:length(files)],1,median)
gene.datExpr <- gene.datExpr[gene.datExpr$median>0,]
gene.datExpr <- gene.datExpr[order(gene.datExpr$ENSEMBL,gene.datExpr$median,decreasing = T),]
gene.datExpr=gene.datExpr[!duplicated(gene.datExpr$ENSEMBL),

row.names(gene.datExpr) <- gene.datExpr$ENSEMBL
gene.datExpr <- gene.datExpr[,1:length(files)]
group_list <- names$type
save(gene.datExpr,group_list,file = "ready.data.Rdata")

tissue<-names$`2`[!duplicated(names$`2`)]

load(file = 'ready.data.Rdata')
getwd();
if(!dir.exists('exports')){dir.create('export')}
if(!dir.exists('data')){dir.create('data')}
rownames(names)<-names_use

############### PCA and hierarchical clustering#################
for (i in tissue) {
  data = gene.datExpr
  namesuse<-names[grep(i,colnames(data)),]
  data = data[,grep(i,colnames(data))]
  nodePar <- list( lab.cex = 0.3, pch = c( NA, 19 ), cex = 0.5, col = "red" )
  hc = hclust(dist(t(data)))

  data=as.data.frame(t(data))
  dat.pca <- PCA(data, graph = F)
  p1<-fviz_pca_ind(dat.pca,
               geom.ind = "point", 
               col.ind = namesuse$type, # color by groups
               addEllipses = T, 
               legend.title = "Groups"
  )

  pdf(paste0('export/',i,'_simple_PCA.pdf'),w=6,h=4)
  print(p1)
  dev.off()
}

##############correlation between samples#################
load("dge.rdata")

pbmc<-FindVariableFeatures(pbmc)
pbmc$tissue_type<-paste0(pbmc$tissue,"_",pbmc$type)
Idents(pbmc)<-pbmc$tissue_type
pbmc<-NormalizeData(pbmc)

Ave<-AggregateExpression(pbmc,return.seurat = T)
Ave<-NormalizeData(Ave)
Ave<-FindVariableFeatures(Ave,nfeatures = 6000)

data<-as.data.frame(Ave@assays$RNA$data)
data<-data[VariableFeatures(Ave),]
cor<-cor(data)
library(RColorBrewer)
library(pheatmap)
metause<-colsplit(rownames(cor),"-",c(1,2))
colnames(metause)<-c("Tissue","Type")
rownames(metause)<-rownames(cor)
metause$Tissue<-as.factor(metause$Tissue)
metause$Type<-as.factor(metause$Type)

ann_colors=list(Tissue=c("Mesenchymal subtype"="#E41A1C","Differentiation subtype"="#377EB8","Immune infiltration subtype"="#4DAF4A","Metabolic subsype"="#984EA3","Muscle subtype"="#FF7F00","Undifferentiated-like subtype"="#FFFF33"),
                Type=c(WT="#619CFF",KO='#00BA38',HE="#F8766D"))
Tissue<-c("#FF8E9E","#FF597B","#C5484E","#C52A28","#EB5431","#EB6D37","#EE833E","#F5E670","#F1FD81","#F6FEAE","#E4E9BE","#CDDBB9","#A0BB72","#779289","#4F7068","#22533A","#314D91","#395D9E","#446FA7","#4BA7BB","#8CB9BD","#CFE1E3","#777BC8","#BE7FFF","#D3B4FA","#DFCCFB","#FDE2FD")

names(Tissue)<-unique(metause$Tissue)
ann_colors$Tissue<-Tissue
pheatmap(cor,clustering_method = "ward.D2",border_color = NA,annotation_col =metause,annotation_row = metause,annotation_colors = ann_colors)



############### DEGs analysis#################
group_listuse<-group_list
gene.datExpruse<-gene.datExpr
library("limma")
design <- model.matrix(~0 + factor(group_listuse))
colnames(design) <- levels(factor(group_listuse))
rownames(design) <- colnames(gene.datExpruse)
contrast.matrix <- makeContrasts("WT-KO",levels=design)
contrast.matrix
# start the DE analysis
fit <- lmFit(gene.datExpruse,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
nrDEG <- topTable(fit2,coef=1,n=Inf)

aa<-nrDEG
aa$gene<-rownames(aa)

logFC <-aa$logFC
adj <- aa$P.Value
gene<- aa$gene

data <- data.frame(logFC=logFC,padj=adj,gene=gene)
data$sig <- "no"
data$sig[data$padj <= 0.05 & data$logFC >= 0.5] <- "up"
data$sig[data$padj <= 0.05 & data$logFC <= -0.5] <- "down"

x_lim <- max(logFC,-logFC)
theme_set(theme_bw())
data$gene<-as.character(data$gene)
data$label=ifelse(data$padj < 0.001 &(data$logFC >= 1 | data$logFC <= -1)| data$gene=="GADD45G"|data$gene=="TGFBR2"|
                    data$gene=="MT2A"|data$gene=="ITM2A"|data$gene=="FNTA"|data$gene=="SPSB3"|data$gene=="LENG8"|data$gene=="TBC1D10C",
                  data$gene,"")

p <- ggplot(data,aes(logFC,-1*log10(padj),
                     color = sig))+geom_point()+
  labs(x="log2(FoldChange)",y="-log10(padj)")
p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))

p <- p +theme(panel.grid =element_blank())+
  theme(axis.line = element_line(size=0))
p <- p +guides(colour = FALSE)+theme_bw()
p<-p+geom_text_repel(data = data, aes(x = data$logFC, 
                                      y = -log10(data$padj), 
                                      label = label),
                     size = 3,box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.8, "lines"), 
                     segment.color = "black", 
                     show.legend = FALSE)
pdf(paste0("./volcano/volcano.pdf"),w=6,h=6)
print(p)
dev.off()
write.csv(data,file = paste0("./volcano/diffgene.csv"))

