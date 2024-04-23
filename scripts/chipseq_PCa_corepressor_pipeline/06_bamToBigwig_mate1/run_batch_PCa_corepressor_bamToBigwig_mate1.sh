#!/bin/sh

mkdir -p /home/chris11/projects/def-stbil30/chris11/20240422_PCa_ChIP_corepressor/output/chip-pipeline_PCA_corepressor-GRCh38/tracks_byReplicate

sbatch scripts/chipseq_PCa_corepressor_pipeline/06_bamToBigwig_mate1/batch_sh/bamToBigwig_mate1_LNCaP_DMSO_SMRT.sh
sbatch scripts/chipseq_PCa_corepressor_pipeline/06_bamToBigwig_mate1/batch_sh/bamToBigwig_mate1_LNCaP_DMSO_NCOR1.sh
