library(GenomicRanges)
library(data.table)
library(pracma)
library(reshape2)
setwd("./scATAC/ACD")
library(ArchR)
addArchRThreads(threads = 6)
proj_callpeak<-readRDS("../proj_callpeak_all.rds")

unique_annos<-unique(proj_callpeak$Anno_type)
anno<-unique_annos[1]

library(dplyr)
cell_meta <- as.data.frame(getCellColData(proj_callpeak))
cell_meta$cell<-rownames(cell_meta)

sampled_cells <- cell_meta %>%
  group_by(Anno_type) %>%
  slice_sample(n = 1000, replace = FALSE) 
sampled_cells_WT<-sampled_cells[sampled_cells$type=="WT",]
sampled_cells_KO<-sampled_cells[sampled_cells$type=="KO",]

setwd("../")

proj_callpeak<-proj_callpeak[sampled_cells_KO$cell,]
aa<-getFragmentsFromProject(
  ArchRProj = proj_callpeak,
  subsetBy = NULL,
  verbose = T,
  logFile = createLogFile("getFragmentsFromProject"))
a1<-aa[[1]]

total_frags <- length(a1)
target_frags <- 1e8  

if (total_frags > target_frags) {
  set.seed(123)  
  sampled_idx <- sample(seq_len(total_frags), size = target_frags, replace = FALSE)
  a1_downsampled <- a1[sampled_idx]
} else {
  message("GRanges already has <= 100M fragments, no downsampling needed.")
  a1_downsampled <- a1
}
a1_downsampled

chrom_sizes <- fread("./mm10_chromSize.txt", header = FALSE)
colnames(chrom_sizes) <- c("chromosome", "size")

split_chromosomes <- function(chrom_sizes, bin_size = 200000) {
  granges_list <- lapply(1:nrow(chrom_sizes), function(i) {
    chr <- chrom_sizes$chromosome[i]
    chr_len <- chrom_sizes$size[i]

    starts <- seq(1, chr_len, by = bin_size)
    ends <- pmin(starts + bin_size - 1, chr_len)

    GRanges(seqnames = chr, ranges = IRanges(start = starts, end = ends))
  })

  do.call(c, granges_list)
}

bins_gr <- split_chromosomes(chrom_sizes)
fragments_gr <- a1_downsampled
midpoints <- start(fragments_gr) + (width(fragments_gr) / 2)
midpoints_gr <- GRanges(seqnames = seqnames(fragments_gr), ranges = IRanges(start = midpoints, width = 1))

counts <- countOverlaps(bins_gr, midpoints_gr)

result <- data.frame(
  chromosome = as.character(seqnames(bins_gr)),
  range = paste0(start(bins_gr), "-", end(bins_gr)),
  count = counts
)
head(result)

result<-result[result$chromosome%in%paste0("chr",c(1:19,"X")),]

i=1
result_use<-result[result$chromosome==paste0('chr',i),]
result_use$x<-1:nrow(result_use)
#ggplot(data = result_use, aes(x = x, y = count,  color = chromosome, shape = chromosome)) + 
#  geom_point(size = 3) + 
#  geom_line(size = 1) + 
#  labs(x = "IGT Block", y = "Mean Net Score") + 
#  geom_hline(aes(yintercept=0), alpha = 0.65)

result_use$count<-as.integer(smooth(result_use$count))
peaks<-findpeaks(result_use$count,ndowns=0)
peaks<-as.data.frame(peaks)
peaks$chr<-paste0('chr',i)
peaks$start<-result_use$range[match(peaks$V3,result_use$x)]
peaks$start<-colsplit(peaks$start,"-",c("a","b"))$a
peaks$end<-result_use$range[match(peaks$V4,result_use$x)]
peaks$end<-colsplit(peaks$end,"-",c("a","b"))$b
peaks$height<-peaks$V1
peaksall<-peaks

for (i in c(2:19,"X")) {
  result_use<-result[result$chromosome==paste0('chr',i),]
  result_use$x<-1:nrow(result_use)
  result_use$count<-as.integer(smooth(result_use$count))
  peaks<-findpeaks(result_use$count,ndowns=0)
  peaks<-as.data.frame(peaks)
  peaks$chr<-paste0('chr',i)
  peaks$start<-result_use$range[match(peaks$V3,result_use$x)]
  peaks$start<-colsplit(peaks$start,"-",c("a","b"))$a
  peaks$end<-result_use$range[match(peaks$V4,result_use$x)]
  peaks$end<-colsplit(peaks$end,"-",c("a","b"))$b
  peaks$height<-peaks$V1
  peaksall<-rbind(peaksall,peaks)
  }

write.csv(peaksall,file = paste0("ACD/all_KO_ACD_peak.csv"))
wtresult<-result 
wtresult$type<-"WT"
koresult<-result 
koresult$type<-"KO"

result<-rbind(wtresult,koresult)
i="X"
for (i in c(1:19,"X")) {
  result_use<-result[result$chromosome==paste0('chr',i),]
result_use$x[result_use$type=="WT"]<-1:nrow(result_use[result_use$type=="WT",])
result_use$x[result_use$type=="KO"]<-1:nrow(result_use[result_use$type=="KO",])
result_use$count[result_use$type=="WT"]<-as.integer(smooth(result_use$count[result_use$type=="WT"]))
result_use$count[result_use$type=="KO"]<-as.integer(smooth(result_use$count[result_use$type=="KO"]))

p<-ggplot(data = result_use, aes(x = x, y = count,  color = type, shape = chromosome)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  labs(x = "IGT Block", y = "Mean Net Score") + 
  geom_hline(aes(yintercept=0), alpha = 0.6)
pdf(paste0("ACD/binplot_all_chr",i,".pdf"),w=10,h=4)
print(p)
dev.off()

}


##---- ACD covplot-----

covplot2_track <- function(
    peak,
    weightCol = NULL,
    xlab = "Chromosome Size (bp)",
    ylab = "",
    title = "ChIP Peaks over Chromosomes",
    chrs = NULL,
    xlim = NULL,
    lower = 1,
    fill_color = c("steelblue", "tomato"),
    track_height = 0.8,
    track_gap = 0.25
) {
  if (!requireNamespace("ChIPseeker", quietly = TRUE)) {
    stop("Package 'ChIPseeker' is required.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.")
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required.")
  }
  
  if (!is.list(peak)) {
    stop("`peak` must be a named list, e.g. list(WT = gr1, KO = gr2)")
  }
  
  if (length(peak) != 2) {
    stop("This function currently expects exactly 2 GRanges objects, e.g. WT and KO.")
  }
  
  if (is.null(names(peak)) || any(names(peak) == "")) {
    names(peak) <- paste0("peak", seq_along(peak))
    warning("Input list had no names, assigned automatically: ",
            paste(names(peak), collapse = ", "))
  }
  
  ltm <- lapply(
    peak,
    ChIPseeker:::getChrCov,
    weightCol = weightCol,
    chrs = chrs,
    xlim = xlim,
    lower = lower
  )
  
  tm <- dplyr::bind_rows(ltm, .id = ".id")
  
  if (nrow(tm) == 0 || length(unique(tm$chr)) == 0) {
    return(
      ggplot2::ggplot(data.frame(x = 1), ggplot2::aes(x)) +
        ggplot2::geom_blank() +
        ggplot2::theme_classic() +
        ggplot2::labs(x = xlab, y = ylab, title = title)
    )
  }
  
  chr.sorted <- ChIPseeker:::sortChrName(as.character(unique(tm$chr)))
  tm$chr <- factor(tm$chr, levels = chr.sorted)
  tm$.id <- factor(tm$.id, levels = names(peak))
  
  if (length(fill_color) == length(peak) &&
      all(ChIPseeker:::is_valid_color(fill_color))) {
    cols <- fill_color
    names(cols) <- names(peak)
  } else {
    cols <- ChIPseeker:::generate_colors(fill_color, n = length(peak))
    names(cols) <- names(peak)
  }

  max_df <- tm %>%
    dplyr::group_by(chr, .id) %>%
    dplyr::summarise(max_value = max(value, na.rm = TRUE), .groups = "drop")
  
  tm <- tm %>%
    dplyr::left_join(max_df, by = c("chr", ".id")) %>%
    dplyr::mutate(
      max_value = ifelse(is.na(max_value) | max_value == 0, 1, max_value),
      value_scaled = value / max_value
    )

  track_ids <- names(peak)
  top_track <- track_ids[1]
  bottom_track <- track_ids[2]
  
  tm <- tm %>%
    dplyr::mutate(
      track_base = dplyr::case_when(
        .id == bottom_track ~ 0,
        .id == top_track ~ track_height + track_gap
      ),
      ymin = track_base,
      ymax = track_base + value_scaled * track_height
    )

  label_df <- data.frame(
    .id = factor(c(top_track, bottom_track), levels = track_ids),
    y = c(track_height + track_gap + track_height / 2, track_height / 2),
    label = c(top_track, bottom_track)
  )
  
  p <- ggplot2::ggplot(tm) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = start,
        xmax = end,
        ymin = ymin,
        ymax = ymax,
        fill = .id,
        color = .id
      ),
      linewidth = 0.15
    ) +
    ggplot2::facet_grid(chr ~ ., scales = "free_x", space = "free_x", switch = "y") +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_x_continuous(
      labels = scales::label_number(scale_cut = scales::cut_si(""))
    ) +
    ggplot2::scale_y_continuous(
      breaks = label_df$y,
      labels = label_df$label,
      limits = c(0, 2 * track_height + track_gap),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      x = xlab,
      y = ylab,
      title = title,
      fill = NULL,
      color = NULL
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      strip.text.y = ggplot2::element_text(angle = 0),
      strip.background = ggplot2::element_blank(),
      legend.position = "none"
    )
  
  if (!is.null(xlim) && !all(is.na(xlim)) && is.numeric(xlim) && length(xlim) == 2) {
    p <- p + ggplot2::xlim(xlim)
  }
  
  return(p)
}

all_KO_ACD_peak <- read_csv("./ACD/all_KO_ACD_peak.csv")
all_WT_ACD_peak <- read_csv("./ACD/all_WT_ACD_peak.csv")
WT.gr <- GRanges(all_WT_ACD_peak$chr,IRanges(all_WT_ACD_peak$start,all_WT_ACD_peak$end))
KO.gr <- GRanges(all_KO_ACD_peak$chr,IRanges(all_KO_ACD_peak$start,all_KO_ACD_peak$end))

p <- covplot2_track(
  peak = list(WT = WT.gr, KO = KO.gr),
  fill_color = c(WT = "#3e6574", KO = "#db7e65"),
  title = "WT vs KO Peaks over Chromosomes"
)

print(p)
ggsave('./plot/WT_KO_all_covplot.pdf',width = 8,height = 6)

