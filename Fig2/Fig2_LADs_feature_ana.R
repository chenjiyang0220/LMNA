library(ChIPpeakAnno)
library(BiocGenerics)
library(reshape2)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(ggplot2)
library(DescTools)
library(tidyr)
library(scales)
library(ggpubr)

setwd("./bulkCUT/analysis/LADs_feature/Fig")
tissue.ls <- c('AdrenalGland' , 'BrainStem'  , 'Cerebrum' ,'Fat', 'Kidney' , 'Lung' , 'Muscle' , 'Skin', 'Spleen' , 'Thymus',
               'Bladder' ,'Testicle' ,'Gonad','LargeIntestine' ,'LymphNode' , 'Pancreas', 'SmallIntestine','Stomach' ,  'Vascular',
               'BoneMarrow', 'Cerebellum' , 'Eye' ,  'Heart', 'Liver' , 'MammaryGland', 'Pituitary' , 'SpinalCord' )

###-----feature of A/B LADs in each tissue(Fig2SB)----

####width fraction of all genome region
chromesize <- fread('/media/ggj/ggj/CJY/tools/mm10/mm10.main.modi.chrom.sizes')
all_len <- sum(chromesize$V2)

## WT AB LADs length analysis
df <- data.frame() 
Gsize <- 2.72e9
for (tissue in tissue.ls) {
  multi <- fread(paste0('./bulkCUT/LAD_calling/analysis/WT_ABLADs/',tissue,'_WT_ABLADs.bed'))
  ALAD <- multi[which(multi$V5 == 'ALAD'),]
  BLAD <- multi[which(multi$V5 == 'BLAD'),]
  ola <- multi[which(multi$V5 == 'ALAD,BLAD'),]
  
  ALAD.gr <- GRanges(seqnames = ALAD$V1, IRanges(start = ALAD$V2, end = ALAD$V3))
  BLAD.gr <- GRanges(seqnames = BLAD$V1, IRanges(start = BLAD$V2, end = BLAD$V3))
  ola.gr <- GRanges(seqnames = ola$V1, IRanges(start = ola$V2, end = ola$V3))
  
  tmp <- data.frame(tissue=tissue,
                    A_only=width(ALAD.gr) %>% sum()/Gsize,
                    AB=width(ola.gr) %>% sum()/Gsize,
                    B_only=(width(BLAD.gr) %>% sum()/Gsize))
  df <- rbind(tmp,df)
}
WriteXLS::WriteXLS(df,'./Fig1/tissue_LADs_width.xlsx',col.names = T,row.names = F)

#round barplot
df$sum <- df$A_only + df$AB + abs(df$B_only)
df <- df[rev(order(df$sum)),]
df[28:29,] <- NA
df[29,1] <- 1
df$id <- c(1:29)
df$id <- as.character(df$id)
df$sum <- as.character(df$sum)
tissue <- factor(df$tissue,levels = df$tissue)
df$angle <-  as.character(cumsum(c(90,-rep(12.41379,28))))

colflg <- c('#A32A31','#DFB6BC','#407BD0')

data <- melt(df)
data$id <- as.numeric(data$id)
data$sum <- as.numeric(data$sum)
data$tissue <- factor(data$tissue, levels = levels(tissue))
p <- ggplot(data,aes(x = tissue, y = value, fill = variable))+
    geom_bar(stat = 'identity',position = 'stack',width = 0.7)+coord_polar()+
    ylim(-0.5, NA)+scale_fill_manual(values = colflg)+ #geom_segment(data = data, mapping =  aes(x = id-0.3 , y = 0.85, xend = id+0.3, yend = 0.85), color = colflg2 , size = 3)+
    theme(panel.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF", size = 2))+
    theme(axis.text.y = element_blank(), axis.ticks = element_blank(),axis.title = element_blank())
p

p <- ggplot(data,aes(x = tissue, y = value, fill = variable))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.8)+coord_polar()+
  ylim(-0.5, NA)+scale_fill_manual(values = colflg)+
  theme(panel.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF", size = 2))+
  geom_text(data=data, aes(x=tissue, y= sum +0.3, label=tissue, hjust=1), 
            color="black", alpha=0.4, size=2, angle=data$angle, inherit.aes = FALSE )+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank(),axis.title = element_blank())
p



axis <- data.frame(
     x=0.5,
     y=c(0, 0.2, 0.4, 0.6,0.8),label=c(0, 20, 40, 60,80))    
axis2 <- data.frame(
    x=0.5, xend=0.6,
    y=c(0, 0.2, 0.4, 0.6,0.8),yend=c(0, 0.2, 0.4, 0.6,0.8)
    )
 p2 <- p+geom_text(data = axis,aes(x=x,y=y,label = label),
                                 inherit.aes = FALSE, hjust=1,size=3)+
     geom_segment(data = axis2,aes(x=x,xend = xend,y=y,yend = yend),inherit.aes = FALSE)+
     annotate(geom = "segment",x=0.6,xend = 0.6,y=0,yend = 0.8)+
     annotate(geom = "text", x=28 ,y=0.6, label = "Width Fraction\nof Genome(%)",angle=90,size =3)
 p2


ggsave('./bulkCUT/analysis/Fig1/tissue_LADs_new.pdf',width = 10,height = 6)



###----------------LADs border cumulative density analysis(Fig2E)----------------
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(ChIPpeakAnno)

LAD_list <- list()
for (sample in tissue.ls) {
  tmp <- toGRanges(paste0('./bulkCUT/LAD_calling/results/LMNA/WT/',sample,'/BED_files_2states/',sample,'_LADs.bed'))
  LAD_list[[sample]] <- tmp
}

BLAD_list <- list()
for (sample in tissue.ls) {
  tmp <- toGRanges(paste0('./bulkCUT/LAD_calling/results/LMNB1/WT/',sample,'/BED_files_2states/',sample,'_LADs.bed'))
  BLAD_list[[sample]] <- tmp
}

KO_BLAD_list <- list()
for (sample in tissue.ls) {
  tmp <- toGRanges(paste0('./bulkCUT/LAD_calling/results/LMNB1/KO/',sample,'/BED_files_2states/',sample,'_LADs.bed'))
  KO_BLAD_list[[sample]] <- tmp
}

LADBorders <- function(LADs, bins, min.distance = 0) {
  # Given a (named) GRangesList of LADs, return a GRangesList with borders
  # Strand information defines left (+) or right (-) border
  # Also, require complete bin object to determine chromosome start and end

  cells <- names(LADs)
  LAD.borders <- GRangesList()

  for (cell in tissue.ls ) {
    # Get LADs and bins for cell
    cell.LADs <- LADs[[cell]]
    # cell.bins <- bins
    # Remove small iLADs and remove small LADs
    cell.LADs <- GenomicRanges::reduce(cell.LADs, min.gapwidth = min.distance)
    cell.LADs <- cell.LADs[width(cell.LADs) > min.distance]
    # Get borders
    cell.borders <- c(GRanges(seqnames = seqnames(cell.LADs),
                              ranges = IRanges(start = start(cell.LADs),
                                               end = start(cell.LADs)),
                              strand = "+"),
                      GRanges(seqnames = seqnames(cell.LADs),
                              ranges = IRanges(start = end(cell.LADs),
                                               end = end(cell.LADs)),
                              strand = "-"))
    cell.borders <- sort(cell.borders, ignore.strand = T)
    cell.borders$cell <- cell
    LAD.borders <- c(LAD.borders, GRangesList(cell.borders))
  }
  names(LAD.borders) <- tissue.ls
  LAD.borders
}

# Prepare borders

LAD_border_list <- LADBorders(LAD_list, bins)
BLAD_border_list <- LADBorders(BLAD_list, bins)
KO_BLAD_border_list <- LADBorders(KO_BLAD_list, bins)

# For all samples, determine distance to nearest border
border_distance <- tibble()
for (sample in tissue.ls) {
  # Get borders for sample
  LAD_sample <- BLAD_list[[sample]]
  LAD_border_sample <- BLAD_border_list[[sample]]
  LAD_borders <- LAD_border_list[[sample]]
  LADs <- LAD_list[[sample]]
  # Distance to nearest border
  dis <- as_tibble(distanceToNearest(LAD_borders, LAD_border_sample, ignore.strand = T)) %>%
    rename_all(~ c("border_idx", "sample_idx", "distance"))
  
  # Filter and add metadata
  dis <- dis %>%
    arrange(sample_idx, distance) %>%
    filter(! duplicated(sample_idx)) %>%
    mutate(sample = sample) %>%
    mutate(border_strand = as.character(strand(LAD_borders))[border_idx],
           sample_strand = as.character(strand(LAD_border_sample))[sample_idx],
           border_within_lad = overlapsAny(LAD_borders[border_idx], LAD_sample,
                                           ignore.strand = T),
           sample_within_lad = overlapsAny(LAD_border_sample[sample_idx], LADs,
                                           ignore.strand = T)) %>%
    mutate(distance = case_when(sample_within_lad == T ~ distance,
                                T ~ - distance))
  # Add to all border distances
  border_distance <- bind_rows(border_distance, 
                                dis)
}

border_distance2 <- tibble()
for (sample in tissue.ls) {
  # Get borders for sample
  LAD_sample <- KO_BLAD_list[[sample]]
  LAD_border_sample <- KO_BLAD_border_list[[sample]]
  LAD_borders <- LAD_border_list[[sample]]
  LADs <- LAD_list[[sample]]
  # Distance to nearest border
  dis <- as_tibble(distanceToNearest(LAD_borders, LAD_border_sample, ignore.strand = T)) %>%
    rename_all(~ c("border_idx", "sample_idx", "distance"))
  # Filter and add metadata
  dis <- dis %>%
    arrange(sample_idx, distance) %>%
    filter(! duplicated(sample_idx)) %>%
    mutate(sample = sample) %>%
    mutate(border_strand = as.character(strand(LAD_borders))[border_idx],
           sample_strand = as.character(strand(LAD_border_sample))[sample_idx],
           border_within_lad = overlapsAny(LAD_borders[border_idx], LAD_sample,
                                           ignore.strand = T),
           sample_within_lad = overlapsAny(LAD_border_sample[sample_idx], LADs,
                                           ignore.strand = T)) %>%
    mutate(distance = case_when(sample_within_lad == T ~ distance,
                                T ~ - distance))
  # Add to all border distances
  border_distance2 <- bind_rows(border_distance2, 
                            dis)
}

border_distance$LAD_type <- 'WT_B'
border_distance2$LAD_type <- 'KO_B'
border_distance_all <- rbind(border_distance,border_distance2)
# Prepare for plotting
border_distance_all <- border_distance_all %>%
  mutate(distance_kb = distance / 1e3,
         distance_kb = (distance_kb %/% 10) * 10) 
# Plot 
border_distance_all %>%
  ggplot(aes(x = distance_kb, col = LAD_type)) +
  stat_ecdf(geom = "line") +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, 
           fill = "grey", alpha = 0.3) +  
  geom_hline(yintercept = 0, col = "black", linetype = "dashed") +
  geom_vline(xintercept = 0, col = "black", linetype = "dashed") +
  facet_wrap(~ sample, ncol = 3, scales = "free_y") +
  coord_cartesian(xlim = c(-150, 150)) +
  xlab("Distance to LAD border (kb)") +
  ylab("Cum density") +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(aspect.ratio = 1, 
        axis.text.x = element_text(angle = 90, hjust = 1))
ggsave('./bulkCUT/analysis/BLAD_trans/BLAD_border_cum_density_tissue_150.pdf',width = 5,height = 12)

####----------CUT&Tag signal profile of LADs borders (Fig2D) ------------

library(rtracklayer)
library(ggplot2)
library(caTools)
library(purrr)
library(data.table)

setwd("./bulkCUT/analysis/profile/LAD_domain/WT/Plot")
####---- all tissue level ----
extend <- 100000
bin_size <- 1000
bins <- 1:(extend*3/bin_size)
antibody <- 'CTCF'
dir <- './bulkCUT/analysis/profile/LAD_domain/WT'
tib.all <- data.frame()
for (tissue in tissue.ls) {
  track_names <- tissue
  LAD <- fread(paste0('./bulkCUT/LAD_calling/analysis/WT_region/ALAD/',tissue,'_LADs.bed'))
  K9 <- fread(paste0('./bulkCUT/LAD_calling/analysis/WT_region/K9me2_3/',tissue,'.K9.LA.bed'))
  nonK9 <- fread(paste0('./bulkCUT/LAD_calling/analysis/WT_region/K9me2_3/',tissue,'.nonK9.LA.bed'))
  nonLAD <- fread(paste0('./bulkCUT/LAD_calling/analysis/WT_region/ALAD/',tissue,'_nonLADs.bed'))
  
  tib <- read_tsv(file.path(dir, paste0(antibody,'/',tissue,'.',antibody,".1kb.matrix.mat.gz")), 
                  skip = 1, 
                  col_names = c("seqnames", "start", "end", 
                                "name", "unknown", "strand",
                                paste(rep(track_names, 
                                          each = length(bins)), 
                                      bins, sep = "_")))
  tib_gather <- tib %>%
    dplyr::select(-seqnames, -start, -end, -name, -unknown, -strand) %>%
    add_column(border = 1:nrow(.),
               region = rep(c("LADs", "K9_LADs","nonK9_LADs","nonLADs"), 
                            times = c(length(rownames(LAD)),
                                      length(rownames(K9)),
                                      length(rownames(nonK9)),
                                      length(rownames(nonLAD))))) %>%
    gather(key, value, -border, -region) %>%
    separate(key, c("track", "bin"), remove = T) %>%
    mutate(track = factor(track, levels = track_names),
           bin = as.numeric(bin),
           bin = (bin * bin_size) - (extend + 0.5 * bin_size),
           bin = bin / 1e3) %>%
    drop_na()
  
  tib.all <- rbind(tib.all, tib_gather)
}  
tib.all$antibody <- antibody 
###all tissue
pdf(paste0('./',antibody,'_profile_1kb.pdf'),width =4,height =3)
tib.all %>%
  drop_na() %>%
  ggplot(aes(x = bin, y = value, col = region, fill = region,
             group = interaction(region,track))) +
  geom_vline(xintercept = 0, col = "black", linetype = "dashed") +
  geom_vline(xintercept = 100, col = "black", linetype = "dashed") + 
  stat_summary(aes(group = region, col = region),fun = mean, geom = "line", size = 0.7)+
  stat_summary(aes(group = region, col = region),fun.data = mean_se, geom = "ribbon", alpha = 0.25, col = NA,
               fun.args = list(mult = 1.96)) +
  xlab("Distance from region border (kb)") +
  ylab("Score") +
  coord_cartesian(xlim = c(-extend/1e3, 2*extend/1e3)) +
  scale_fill_brewer(palette = "Set2", name = "Region class") +
  scale_color_brewer(palette = "Set2", name = "Region class") +
  theme_bw() +
  theme(aspect.ratio = 0.6, 
        axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


