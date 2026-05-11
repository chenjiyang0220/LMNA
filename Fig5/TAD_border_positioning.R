
## TAD_border_positioning.R
##
## Compare positional shifts of CTCF-enriched strong TAD boundaries
## between WT and KO on chr13, and visualize cumulative distributions
## of border-to-border distances (similar to LAD_boder_positioning_cum_density.R).
##
## Usage (from terminal, in an R-enabled environment):
##   Rscript TAD_border_positioning.R


library(GenomicRanges)
library(rtracklayer)
library(tidyverse)


## -------------------------------------------------------------------
## Input / output paths
## -------------------------------------------------------------------

base_dir <- "/media/ggj/ggj/CJY/nature_WXY/Hi-c/analysis/TAD"

wt_bed <- file.path(base_dir, "WT_chr13_ctcf_enriched_strong_boundaries_30kb.bed")
ko_bed <- file.path(base_dir, "KO_chr13_ctcf_enriched_strong_boundaries_30kb.bed")

out_dir <- base_dir
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

pdf_pairwise   <- file.path(out_dir, "WT_KO_chr13_TAD_border_cum_density_pairwise_150kb.pdf")

rda_out        <- file.path(out_dir, "WT_KO_chr13_TAD_border_distance.rda")
tsv_out        <- file.path(out_dir, "WT_KO_chr13_TAD_border_distance.tsv")

## -------------------------------------------------------------------
## Load WT / KO TAD borders (CTCF-enriched strong boundaries)
## -------------------------------------------------------------------

wt_borders_df <- read_tsv(
  wt_bed,
  col_names = c("chrom", "start", "end"),
  col_types = "cii"
)
ko_borders_df <- read_tsv(
  ko_bed,
  col_names = c("chrom", "start", "end"),
  col_types = "cii"
)

wt_borders <- GRanges(
  seqnames = wt_borders_df$chrom,
  ranges   = IRanges(start = wt_borders_df$start, end = wt_borders_df$end)
)
ko_borders <- GRanges(
  seqnames = ko_borders_df$chrom,
  ranges   = IRanges(start = ko_borders_df$start, end = ko_borders_df$end)
)

## -------------------------------------------------------------------
## Compute nearest-border distances (WT -> KO, and KO -> WT)
## -------------------------------------------------------------------

border_distance_all <- tibble()

if (length(wt_borders) > 0 && length(ko_borders) > 0) {
  hits_wt_to_ko <- distanceToNearest(wt_borders, ko_borders, ignore.strand = TRUE)
  
  df_wt_to_ko <- as_tibble(hits_wt_to_ko) %>%
    rename(
      query_idx   = queryHits,
      subject_idx = subjectHits,
      distance_abs = distance
    ) %>%
    mutate(
      wt_pos = start(wt_borders)[query_idx],
      ko_pos = start(ko_borders)[subject_idx],
      distance_signed = ko_pos - wt_pos,
      group = "WT_to_KO"
    )
  
  border_distance_all <- bind_rows(border_distance_all, df_wt_to_ko)
}

if (length(ko_borders) > 0 && length(wt_borders) > 0) {
  hits_ko_to_wt <- distanceToNearest(ko_borders, wt_borders, ignore.strand = TRUE)
  
  df_ko_to_wt <- as_tibble(hits_ko_to_wt) %>%
    rename(
      query_idx   = queryHits,
      subject_idx = subjectHits,
      distance_abs = distance
    ) %>%
    mutate(
      ko_pos = start(ko_borders)[query_idx],
      wt_pos = start(wt_borders)[subject_idx],
      distance_signed = wt_pos - ko_pos,
      group = "KO_to_WT"
    )
  
  border_distance_all <- bind_rows(border_distance_all, df_ko_to_wt)
}

if (nrow(border_distance_all) == 0) {
  stop("No border distances could be computed; check WT/KO border BEDs.")
}

border_distance_all <- border_distance_all %>%
  mutate(
    distance_kb = distance_signed / 1e3,
    distance_kb_binned = (distance_kb %/% 30) * 30
  )

write_tsv(border_distance_all, tsv_out)
save(border_distance_all, file = rda_out)

## -------------------------------------------------------------------
## Plotting: cumulative density of border distance
## -------------------------------------------------------------------
library(tidyverse)

base_dir <- "/media/ggj/ggj/CJY/nature_WXY/Hi-c/analysis/TAD"
out_dir  <- base_dir

tsv_out <- file.path(out_dir, "WT_KO_chr13_TAD_border_distance.tsv")
border_distance_all <- read_tsv(tsv_out, show_col_types = FALSE)

message("Plotting cumulative distance distributions ...")

border_distance_all_abs <- border_distance_all %>%
  dplyr::mutate(
    distance_kb_abs = abs(distance_kb),
    distance_kb_abs_binned = (distance_kb_abs %/% 10) * 10
  )

p_all <- border_distance_all_abs %>%
  ggplot(aes(x = distance_kb_abs_binned)) +
  stat_ecdf(geom = "line", col = "steelblue") +
  annotate(
    "rect",
    xmin = 0, xmax = Inf,
    ymin = -Inf, ymax = Inf,
    fill = "grey", alpha = 0.3
  ) +
  geom_hline(yintercept = 0, col = "black", linetype = "dashed") +
  coord_cartesian(xlim = c(0, 600)) +
  xlab("Absolute distance between TAD borders (kb)") +
  ylab("Cumulative density") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
pdf_all        <- file.path(out_dir, "WT_KO_chr13_TAD_border_cum_density_all_wide.pdf")
ggsave(pdf_all, p_all, width = 8, height = 4)

message("Finished. Outputs written to: ", out_dir)

