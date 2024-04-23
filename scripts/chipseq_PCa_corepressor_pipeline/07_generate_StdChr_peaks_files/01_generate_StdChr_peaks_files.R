setwd("/Users/chris/Desktop/20240422_PCa_ChIP_corepressor")

library(tidyverse)
library(GenomicRanges)

#
fastq_list_filename <- "chipseq_PCA_corepressor_fastq_list.txt"
df <- read_tsv(file.path("input", "chipseq_PCa_corepressor", fastq_list_filename)) %>% 
  dplyr::filter(type == "PAIRED_END",
                cell_line == "LNCaP",
                condition == "DMSO") 
output_pipeline_dir <- "chip-pipeline_PCA_corepressor-GRCh38"
# script_pipeline_dir <- "chipseq_H3K27ac_0_25m_pipeline"
  
#
ENCODE_elr <- rtracklayer::import("input/ENCODE_exclusion_list_regions_ENCFF356LFX.bed")
stdChr <- paste0("chr", c(seq(1:22), "X", "Y"))
peaks_dir <- file.path("output", output_pipeline_dir, "peak_call")

for (i in 1:nrow(df)) {
  sample_name <- df$sample_name[i]
  antibody_factor <- df$antibody[i]
  condition <- df$condition[i]
  # rep_i <- df$isogenic_replicate[i]
  
  # line <- df[i, ]
  # timepoint_i <- line %>% pull(time_point)
  # rep_i <- line %>% pull(isogenic_replicate)
  message("# ", sample_name)
  
  basename <- paste("LNCaP", condition, antibody_factor, sep = "_")
  message("   > ", basename)
  
  peaks_path <- file.path(peaks_dir, basename, antibody_factor, paste0(basename, ".", antibody_factor, "_peaks.narrowPeak"))
  
  #
  peaks_raw <- rtracklayer::import(peaks_path)
  message("\tRaw number of peaks : ", length(peaks_raw))
  
  #
  message("\t> Remove ", length(peaks_raw[!seqnames(peaks_raw) %in% stdChr]), " regions not on standard chromosomes")
  peaks_stdchr <- keepSeqlevels(peaks_raw, stdChr[stdChr %in% seqlevels(peaks_raw)], pruning.mode = "coarse")
  # message("\tNumber of peaks on standard chromosomes : ", length(peaks_stdchr))
  
  #
  message("\t> Remove ", length(subsetByOverlaps(peaks_stdchr, ENCODE_elr)), " excluded ENCODE regions")
  peaks_notbl <- subsetByOverlaps(peaks_stdchr, ENCODE_elr, invert = TRUE)
  
  # message("\tNumber of peaks not on the ENCODE exclusion list regions : ", length(peaks_notbl))
   
  #
  message("\tFinal number of peaks : ", length(peaks_notbl))
  output_filename <- paste0(basename, ".", antibody_factor, "_peaks.narrowPeak.stdchr.bed")
  output_path <- file.path(peaks_dir, basename, antibody_factor, output_filename)
  message("\t", output_path)
  rtracklayer::export(peaks_notbl, con = output_path, format = "bed")
}
