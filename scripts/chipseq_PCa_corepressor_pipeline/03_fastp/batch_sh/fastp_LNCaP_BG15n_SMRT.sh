#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --account=def-stbil30
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL


mkdir -p /home/chris11/projects/def-stbil30/chris11/20240422_PCa_ChIP_corepressor/output/chip-pipeline_PCA_corepressor-GRCh38/fastp_output/LNCaP_BG15n_SMRT


mkdir -p /home/chris11/projects/def-stbil30/chris11/20240422_PCa_ChIP_corepressor/raw/chipseq_PCa_corepressor/fastp_output


/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/fastp/0.23.4/bin/fastp --in1 /home/chris11/projects/def-stbil30/chris11/20240422_PCa_ChIP_corepressor/raw/chipseq_PCa_corepressor/raw_fastq/SRR19639160_1.fastq.gz --in2 /home/chris11/projects/def-stbil30/chris11/20240422_PCa_ChIP_corepressor/raw/chipseq_PCa_corepressor/raw_fastq/SRR19639160_2.fastq.gz --detect_adapter_for_pe --overrepresentation_analysis --overrepresentation_sampling 10 --thread 8 --length_required 74 --out1 /home/chris11/projects/def-stbil30/chris11/20240422_PCa_ChIP_corepressor/raw/chipseq_PCa_corepressor/fastp_output/LNCaP_BG15n_SMRT_1.fastq.gz --out2 /home/chris11/projects/def-stbil30/chris11/20240422_PCa_ChIP_corepressor/raw/chipseq_PCa_corepressor/fastp_output/LNCaP_BG15n_SMRT_2.fastq.gz --html /home/chris11/projects/def-stbil30/chris11/20240422_PCa_ChIP_corepressor/output/chip-pipeline_PCA_corepressor-GRCh38/fastp_output/LNCaP_BG15n_SMRT/LNCaP_BG15n_SMRT_fastp_report.html --json /home/chris11/projects/def-stbil30/chris11/20240422_PCa_ChIP_corepressor/output/chip-pipeline_PCA_corepressor-GRCh38/fastp_output/LNCaP_BG15n_SMRT/LNCaP_BG15n_SMRT_fastp_report.json
