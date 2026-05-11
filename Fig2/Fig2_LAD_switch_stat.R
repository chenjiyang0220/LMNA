library(GenomicRanges)
library(ChIPpeakAnno)
library(BiocGenerics)
library(reshape2)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(datarium)
library(tidyverse)
library(rstatix)

setwd("./bulkCUT/analysis/BLAD_trans")

tissue.ls <- c('AdrenalGland' ,  'LargeIntestine' , 'Pituitary' ,    'Thymus',
               'Bladder' ,    'Eye' ,   'Liver' ,         'Skin' ,       'Vascular',
               'BoneMarrow'  , 'Fat'  ,  'Lung',           'SmallIntestine',
               'BrainStem',   'Gonad' ,  'MammaryGland' ,   'SpinalCord',
               'Cerebellum' , 'Heart' , 'Muscle'   ,     'Spleen',
               'Cerebrum' ,     'Kidney', 'Pancreas'   ,     'Stomach','LymphNode','Testicle')

Gsize <- 2.72e9

#---- LAD-reorganized regions analysis----
dir <- './bulkCUT/LAD_calling/analysis/WT_region/WT_ABLADs/'
KO_dir <- './bulkCUT/LAD_calling/results/LMNB1/KO/'
for (tissue in tissue.ls) {
  WT_AB <- fread(paste0(dir,tissue,'_WT_ABLADs.bed'))
  onlyA <- WT_AB[which(WT_AB$V5 == 'ALAD'),]
  onlyB <- WT_AB[which(WT_AB$V5 == 'BLAD'),]
  A_B <- WT_AB[which(WT_AB$V5 == 'ALAD,BLAD'),]
  write.table(onlyA,paste0(dir,tissue,'_onlyA.bed'),quote = F,sep = '\t',col.names = F, row.names = F)
  write.table(onlyB,paste0(dir,tissue,'_onlyB.bed'),quote = F,sep = '\t',col.names = F, row.names = F)
  write.table(A_B,paste0(dir,tissue,'_A_B.bed'),quote = F,sep = '\t',col.names = F, row.names = F)
}

dir <- './bulkCUT/LAD_calling/analysis/KO_region/BLAD_switch/'
df <- data.frame()
for (tissue in tissue.ls) {
  AtoB <- fread(paste0(dir,tissue,'_AtoB.bed'))
  NontoB <- fread(paste0(dir,tissue,'_NontoB.bed'))
  ABtoNon <- fread(paste0(dir,tissue,'_ABtoNon.bed'))
  BtoNon <- fread(paste0(dir,tissue,'_BtoNon.bed'))
  
  AtoB.gr <- GRanges(seqnames = AtoB$V1, IRanges(start = AtoB$V2, end = AtoB$V3))
  NontoB.gr <- GRanges(seqnames = NontoB$V1, IRanges(start = NontoB$V2, end = NontoB$V3))
  ABtoNon.gr <- GRanges(seqnames = ABtoNon$V1, IRanges(start = ABtoNon$V2, end = ABtoNon$V3))
  BtoNon.gr <- GRanges(seqnames = BtoNon$V1, IRanges(start = BtoNon$V2, end = BtoNon$V3))
  
  tmp <- data.frame(tissue = tissue,
                    AtoB = -width(AtoB.gr) %>% sum() / Gsize,
                    NontoB = -width(NontoB.gr) %>% sum() / Gsize,
                    ABtoNon = -(width(ABtoNon.gr) %>% sum() / Gsize),
                    BtoNon = -(width(BtoNon.gr) %>% sum() / Gsize)
  )
  df <- rbind(df,tmp)
}


#####-------- LAD-reorganized regions GC content --------

library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome)
library(biovizBase)

dir <- './bulkCUT/LAD_calling/analysis/KO_region/BLAD_switch/'
tissue.ls <- c('AdrenalGland' ,  'LargeIntestine' , 'Pituitary' ,    'Thymus',
               'Bladder' ,    'Eye' ,   'Liver' ,         'Skin' ,       'Vascular',
               'BoneMarrow'  , 'Fat'  ,  'Lung',           'SmallIntestine',
               'BrainStem',   'Gonad' ,  'MammaryGland' ,   'SpinalCord',
               'Cerebellum' , 'Heart' , 'Muscle'   ,     'Spleen',
               'Cerebrum' ,     'Kidney', 'Pancreas'   ,     'Stomach','LymphNode','Testicle')
df <- data.frame()
for (tissue in tissue.ls) {
  AtoB <- fread(paste0(dir,tissue,'_AtoB.bed'))
  NontoB <- fread(paste0(dir,tissue,'_NontoB.bed'))
  ABtoNon <- fread(paste0(dir,tissue,'_ABtoNon.bed'))
  BtoNon <- fread(paste0(dir,tissue,'_BtoNon.bed'))
  
  AtoB.gr <- GRanges(seqnames = AtoB$V1, IRanges(start = AtoB$V2, end = AtoB$V3))
  NontoB.gr <- GRanges(seqnames = NontoB$V1, IRanges(start = NontoB$V2, end = NontoB$V3))
  ABtoNon.gr <- GRanges(seqnames = ABtoNon$V1, IRanges(start = ABtoNon$V2, end = ABtoNon$V3))
  BtoNon.gr <- GRanges(seqnames = BtoNon$V1, IRanges(start = BtoNon$V2, end = BtoNon$V3))
  
  seqlengths <- c(chr1=	195471971,chr10=	130694993,chr11=	122082543, chr12=	120129022, chr13=	120421639, chr14=124902244,
                              chr15=	104043685,chr16=98207768, chr17=94987271, chr18=90702639, chr19=	61431566,
                              chr2=	182113224, chr3=	160039680, chr4=	156508116,chr5=	151834684,
                              chr6=	149736546, chr7=	145441459, chr8=129401213, chr9=	124595110, chrX=	171031299, chrY=91744698 )
  chrom <- seqlengths(AtoB.gr)
  seqlength <- seqlengths[names(seqlengths) %in% names(chrom)]
  seqlengths(AtoB.gr) <- seqlength
  AtoB.gr <- trim(AtoB.gr)
  
  chrom <- seqlengths(NontoB.gr)
  seqlength <- seqlengths[names(seqlengths) %in% names(chrom)]
  seqlengths(NontoB.gr) <- seqlength
  NontoB.gr <- trim(NontoB.gr)
  
  chrom <- seqlengths(ABtoNon.gr)
  seqlength <- seqlengths[names(seqlengths) %in% names(chrom)]
  seqlengths(ABtoNon.gr) <- seqlength
  ABtoNon.gr <- trim(ABtoNon.gr)
  
  chrom <- seqlengths(BtoNon.gr)
  seqlength <- seqlengths[names(seqlengths) %in% names(chrom)]
  seqlengths(BtoNon.gr) <- seqlength
  BtoNon.gr <- trim(BtoNon.gr)
  
  AtoB_GC <- GCcontent(Mmusculus, AtoB.gr)
  NontoB_GC <- GCcontent(Mmusculus, NontoB.gr)
  ABtoNon_GC <- GCcontent(Mmusculus, ABtoNon.gr)
  BtoNon_GC <- GCcontent(Mmusculus, BtoNon.gr)
  
  tmp <-  data.frame(tissue = tissue,
                   AtoB_GCcontent = round(mean(AtoB_GC),2),
                   NontoB_GCcontent = round(mean(NontoB_GC),2),
                   ABtoNon_GCcontent = round(mean(ABtoNon_GC),2),
                   BtoNon_GCcontent = round(mean(BtoNon_GC),2)
                   )  
  df <- rbind(df,tmp)
  
}

dat <- melt(df)
colflg <- c('#C3D78E','#F6F2B8','#93DBDC','#9C84A3')
ggplot(dat, aes(x=variable, y=value, fill = variable)) + geom_violin(trim=FALSE,width=0.8,color="white")+
  stat_summary(fun.data = "mean_sdl",
                 fun.args = list(mult = 1), geom = "pointrange", color = "white") + 
  theme(panel.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF", size = 2))+
  theme(axis.line = element_line(linetype = 1, color = "black", size=0.3))+
  theme(axis.text = element_text(hjust = 1,angle=45)) +
  scale_fill_manual(values = colflg)
ggsave('./BLAD_trans_GCcontent.pdf',width = 6, height = 4)

######-------- motif enrichment analysis in LAD-reorganized regions (FigS2K) -----------
library(GenomicRanges)
library(rtracklayer)
library(chromVARmotifs)
library(BSgenome.Mmusculus.UCSC.mm10) 
library(data.table) 
library(TFBSTools)
library(JASPAR2020)
library(motifmatchr)
library(TFBSTools)
library(reshape2)

tissue.ls <- c('AdrenalGland', 'LargeIntestine', 'Pituitary', 'Thymus',
               'Bladder', 'Eye', 'Liver', 'Skin', 'Vascular',
               'BoneMarrow', 'Fat', 'Lung', 'SmallIntestine',
               'BrainStem', 'Gonad', 'MammaryGland', 'SpinalCord',
               'Cerebellum', 'Heart', 'Muscle', 'Spleen',
               'Cerebrum', 'Kidney', 'Pancreas', 'Stomach', 'LymphNode', 'Testicle')

df <- data.frame()

species <- "mouse"
source('./Code/source/Load_txdb_annodb.R')
Load_txdb_annodb(species = species)

# load Motif 
data("encode_pwms") # chromVARmotifs package
motifs_mouse <- encode_pwms
dir <- "./bulkCUT/LAD_calling/analysis/KO_region/BLAD_switch_new/" 

for (tissue in tissue.ls) {
  message(tissue)
  LADtoNon <- fread(paste0(dir, tissue, '_LADtoNon.bed'))
  NontoB <- fread(paste0(dir, tissue, '_NontoB.bed'))
  NontoNon <- fread(paste0(dir, tissue, '_NontoNon.bed'))
  
  LADtoNon.gr <- GRanges(seqnames = LADtoNon$V1, IRanges(start = LADtoNon$V2, end = LADtoNon$V3))
  NontoB.gr <- GRanges(seqnames = NontoB$V1, IRanges(start = NontoB$V2, end = NontoB$V3))
  NontoNon.gr <- GRanges(seqnames = NontoNon$V1, IRanges(start = NontoNon$V2, end = NontoNon$V3))
  seqlengths <- c(chr1=	195471971,chr10=	130694993,chr11=	122082543, chr12=	120129022, chr13=	120421639, chr14=124902244,
                  chr15=	104043685,chr16=98207768, chr17=94987271, chr18=90702639, chr19=	61431566,
                  chr2=	182113224, chr3=	160039680, chr4=	156508116,chr5=	151834684,
                  chr6=	149736546, chr7=	145441459, chr8=129401213, chr9=	124595110, chrX=	171031299, chrY=91744698 )
  chrom <- seqlengths(LADtoNon.gr)
  seqlength <- seqlengths[names(seqlengths) %in% names(chrom)]
  seqlengths(LADtoNon.gr) <- seqlength
  LADtoNon.gr <- trim(LADtoNon.gr)
  
  chrom <- seqlengths(NontoB.gr)
  seqlength <- seqlengths[names(seqlengths) %in% names(chrom)]
  seqlengths(NontoB.gr) <- seqlength
  NontoB.gr <- trim(NontoB.gr)
  
  chrom <- seqlengths(NontoNon.gr)
  seqlength <- seqlengths[names(seqlengths) %in% names(chrom)]
  seqlengths(NontoNon.gr) <- seqlength
  NontoNon.gr <- trim(NontoNon.gr)

  for (region_type in c("LADtoNon", "NontoB", "NontoNon")) {
    if (region_type == "LADtoNon") {
      subject.gr <- LADtoNon.gr
    } else if (region_type == "NontoB") {
      subject.gr <- NontoB.gr
    } else if (region_type == "NontoNon") {
      subject.gr <- NontoNon.gr
    }
    
    TF_pos <- matchMotifs(pwms = motifs_mouse, subject = subject.gr, genome = BSgenome.Mmusculus.UCSC.mm10,
                          out = "positions", p.cutoff = 5e-05, w = 7)
    
    scores_TF <- sapply(TF_pos, function(x) x$score)
    median_score_TF <- sapply(scores_TF, median)
    
    #TF_names <- colsplit(names(median_score_TF), pattern = '_', names = c('V1', 'V2', 'V3'))[,'V3']
    TF_names<-names(median_score_TF)

    df1 <- data.frame(TF = TF_names, scores = median_score_TF, tissue = tissue, region_type = region_type)
    df1 <- data.frame(
      TF = TF_names,
      scores = median_score_TF,
      tissue = rep(tissue, length(TF_names)),  
      region_type = rep(region_type, length(TF_names))  
    )

    df <- rbind(df, df1)
  }
}
df$TF <- colsplit(df$TF,pattern = '_',names = c('V1','V2'))$V1

df <- df[order(df$scores, decreasing = TRUE),]
df$rank <- 1:length(rownames(df))

