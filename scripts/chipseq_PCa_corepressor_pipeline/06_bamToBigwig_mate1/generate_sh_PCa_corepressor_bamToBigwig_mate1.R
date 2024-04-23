setwd("/Users/chris/Desktop/20240422_PCa_ChIP_corepressor")

library(tidyverse)

#
fastq_list_filename <- "chipseq_PCA_corepressor_fastq_list.txt"
df <- read_tsv(file.path("input", "chipseq_PCa_corepressor", fastq_list_filename)) %>% 
  dplyr::filter(type == "PAIRED_END",
                cell_line == "LNCaP",
                condition == "DMSO") 
fastq_folder <- "chipseq_PCa_corepressor"
output_pipeline_dir <- "chip-pipeline_PCA_corepressor-GRCh38"
script_pipeline_dir <- "chipseq_PCa_corepressor_pipeline"
workdir <- "/home/chris11/projects/def-stbil30/chris11/20240422_PCa_ChIP_corepressor"

alignment_dir <- file.path("output", output_pipeline_dir, "alignment")
tracks_dir <- file.path("output", output_pipeline_dir, "tracks_byReplicate")
message("mkdir -p ", file.path(workdir, tracks_dir))


header_sh <- c("#!/bin/sh",
               "#SBATCH --time=3:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --ntasks-per-node=1",
               "#SBATCH --cpus-per-task=16",
               "#SBATCH --mem-per-cpu=8G",
               "#SBATCH --account=def-stbil30",
               "#SBATCH --mail-user=christophe.tav@gmail.com",
               "#SBATCH --mail-type=ALL")

ENCODE_blacklist <- "input/ENCODE_exclusion_list_regions_ENCFF356LFX.bed"

norm <- "RPKM"

for (i in 1:nrow(df)) {
  sample_name <- df$sample_name[i]
  antibody_factor <- df$antibody[i]
  # message(" > ", sample_name)
  
  bam <- file.path(workdir, alignment_dir, sample_name, antibody_factor, paste0(sample_name, ".",  antibody_factor, ".sorted.dup.filtered.bam"))
  bw <- file.path(workdir, tracks_dir, paste0(sample_name, "_", norm, "_mate1.bw"))
  
  # print(paste0(" > bam : ", bam))
  # print(paste0(" > bigwig : ", bw))
  # message(basename(bw))
  call_bamCoverage <- paste("bamCoverage", "--extendReads", 225, "--binSize", 10,
                            "--smoothLength", 30,
                            "-p", 16,
                            "--samFlagInclude", 64,
                            "--normalizeUsing", norm,
                            "--blackListFileName", file.path(workdir, ENCODE_blacklist),
                            "-b", bam,
                            "-o", bw)
  # message(call_bamCoverage)
  # system(call_bamCoverage)
  
  file_sh <- file.path("scripts", script_pipeline_dir, "06_bamToBigwig_mate1", "batch_sh",
                       paste0("bamToBigwig_mate1_", sample_name, ".sh"))
  message("sbatch ", file_sh)
  fileConn <- file(file_sh)
  writeLines(c(header_sh, "\n", call_bamCoverage), fileConn)
  close(fileConn)
}
