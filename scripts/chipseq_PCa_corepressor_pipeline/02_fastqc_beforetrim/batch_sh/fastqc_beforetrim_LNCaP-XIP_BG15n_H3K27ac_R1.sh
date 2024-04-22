#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --account=def-stbil30
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL


mkdir -p /home/chris11/projects/def-stbil30/chris11/20240422_PCa_ChIP_corepressor/output/chip-pipeline_PCA_corepressor-GRCh38/fastqc_beforetrim_output/LNCaP-XIP_BG15n_H3K27ac_R1


/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/fastqc/0.11.9/fastqc --outdir /home/chris11/projects/def-stbil30/chris11/20240422_PCa_ChIP_corepressor/output/chip-pipeline_PCA_corepressor-GRCh38/fastqc_beforetrim_output/LNCaP-XIP_BG15n_H3K27ac_R1 --format fastq /home/chris11/projects/def-stbil30/chris11/20240422_PCa_ChIP_corepressor/raw/chipseq_PCa_corepressor/raw_fastq/SRR19639164_1.fastq.gz
