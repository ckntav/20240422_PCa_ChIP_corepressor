#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --account=def-stbil30
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL


bamCoverage --extendReads 225 --binSize 10 --smoothLength 30 -p 16 --samFlagInclude 64 --normalizeUsing RPKM --blackListFileName /home/chris11/projects/def-stbil30/chris11/20240422_PCa_ChIP_corepressor/input/ENCODE_exclusion_list_regions_ENCFF356LFX.bed -b /home/chris11/projects/def-stbil30/chris11/20240422_PCa_ChIP_corepressor/output/chip-pipeline_PCA_corepressor-GRCh38/alignment/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.dup.filtered.bam -o /home/chris11/projects/def-stbil30/chris11/20240422_PCa_ChIP_corepressor/output/chip-pipeline_PCA_corepressor-GRCh38/tracks_byReplicate/LNCaP_DMSO_SMRT_RPKM_mate1.bw
