setwd("/Users/chris/Desktop/20240422_PCa_ChIP_corepressor")

library(tidyverse)

##### module load fastqc

#
fastq_list_filename <- "chipseq_PCA_corepressor_fastq_list.txt"
df <- read_tsv(file.path("input", "chipseq_PCa_corepressor", fastq_list_filename))
fastq_folder <- "chipseq_PCa_corepressor/raw_fastq"
output_pipeline_dir <- "chip-pipeline_PCA_corepressor-GRCh38"
script_pipeline_dir <- "chipseq_PCa_corepressor_pipeline"
workdir <- "/home/chris11/projects/def-stbil30/chris11/20240422_PCa_ChIP_corepressor"

#
header_sh <- c("#!/bin/sh",
               "#SBATCH --time=3:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --ntasks-per-node=1",
               "#SBATCH --cpus-per-task=4",
               "#SBATCH --mem-per-cpu=8G",
               "#SBATCH --account=def-stbil30",
               "#SBATCH --mail-user=christophe.tav@gmail.com",
               "#SBATCH --mail-type=ALL")

fastqc_path <- "/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/fastqc/0.11.9/fastqc"

for (i in 1:nrow(df)) {
  sample_name <- df$sample_name[i]
  # message("# ", i, " | ", sample_name)
  
  ##### R1
  #
  # message(" > R1")
  fastq_R1_filename <- df$fastq_R1_filename[i]
  fastq_R1_filepath <- file.path(workdir, "raw", fastq_folder, fastq_R1_filename)
  output_fastqc_R1_path <- file.path(workdir, "output", output_pipeline_dir, "fastqc_beforetrim_output", paste(sep = "_", sample_name, "R1"))
  call_mkdir_R1 <- paste("mkdir", "-p", output_fastqc_R1_path)
  
  #
  call_fastqc_R1 <- paste(fastqc_path,
                          "--outdir", output_fastqc_R1_path,
                          "--format", "fastq",
                          fastq_R1_filepath)
  
  #
  file_sh1 <- file.path("scripts", script_pipeline_dir , "02_fastqc_beforetrim", "batch_sh",
                       paste0("fastqc_beforetrim_", sample_name, "_R1.sh"))
  message("sbatch ", file_sh1)
  fileConn1 <- file(file_sh1)
  writeLines(c(header_sh, "\n", call_mkdir_R1, "\n", call_fastqc_R1), fileConn1)
  close(fileConn1)

  
  ##### R2
  #
  # message(" > R2")
  fastq_R2_filename <- df$fastq_R2_filename[i]
  if (!is.na(fastq_R2_filename)) {
    fastq_R2_filepath <- file.path(workdir, "raw", fastq_folder, fastq_R2_filename)
    output_fastqc_R2_path <- file.path(workdir, "output", output_pipeline_dir, "fastqc_beforetrim_output", paste(sep = "_", sample_name, "R2"))
    call_mkdir_R2 <- paste("mkdir", "-p", output_fastqc_R2_path)
    
    #
    call_fastqc_R2 <- paste(fastqc_path,
                            "--outdir", output_fastqc_R2_path,
                            "--format", "fastq",
                            fastq_R2_filepath)
    
    #
    file_sh2 <- file.path("scripts", script_pipeline_dir , "02_fastqc_beforetrim", "batch_sh",
                          paste0("fastqc_beforetrim_", sample_name, "_R2.sh"))
    message("sbatch ", file_sh2)
    fileConn2 <- file(file_sh2)
    writeLines(c(header_sh, "\n", call_mkdir_R2, "\n", call_fastqc_R2), fileConn2)
    close(fileConn2)
  }
}


