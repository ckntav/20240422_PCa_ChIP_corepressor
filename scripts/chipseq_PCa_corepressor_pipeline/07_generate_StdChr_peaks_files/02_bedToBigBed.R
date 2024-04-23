setwd("/Users/chris/Desktop/20240422_PCa_ChIP_corepressor")

library(tidyverse)

#
wdir_path <- "/Users/chris/Desktop/20240422_PCa_ChIP_corepressor"
output_pipeline_dir <- "chip-pipeline_PCA_corepressor-GRCh38"
tracks_dir <- file.path(wdir_path, "output", output_pipeline_dir, "bigbed_byReplicate")
chromsizeshg38 <- file.path(wdir_path, "input/genome_annot/hg38.chrom.sizes")

#
fastq_list_filename <- "chipseq_PCA_corepressor_fastq_list.txt"
df <- read_tsv(file.path("input", "chipseq_PCa_corepressor", fastq_list_filename)) %>% 
  dplyr::filter(type == "PAIRED_END",
                cell_line == "LNCaP",
                condition == "DMSO") 
peaks_antibody_dir <- file.path(wdir_path, "output", output_pipeline_dir, "peak_call")

# for (i in (which(df$antibody != "WCE"))) {
for (i in 1:nrow(df)) {
  sample_name <- df$sample_name[i]
  antibody_factor <- df$antibody[i]
  condition <- df$condition[i]
  # rep_i <- df$isogenic_replicate[i]
  message("# ", sample_name)
  
  #
  basename <- paste("LNCaP", condition, antibody_factor, sep = "_")
  
  input_bed <- paste0(basename, ".",  antibody_factor, "_peaks.narrowPeak.stdchr.bed")
  input_path <- file.path(peaks_antibody_dir, basename,  antibody_factor, input_bed)
  output_bb <- str_replace(input_bed, pattern = "bed", replacement = "bb")
  output_path <- file.path(tracks_dir, output_bb)
  
  # 
  call_bedToBigBed <- paste("/Users/chris/Documents/software/bedToBigBed",
                            "-type=bed3+3",
                            input_path, chromsizeshg38, output_path)
  # message(call_bedToBigBed)
  system(call_bedToBigBed)
}
