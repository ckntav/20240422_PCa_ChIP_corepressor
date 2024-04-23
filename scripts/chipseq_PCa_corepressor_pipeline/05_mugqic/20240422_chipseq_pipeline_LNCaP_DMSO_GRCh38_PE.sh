#!/bin/bash
# Exit immediately on error

set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq SLURM Job Submission Bash script
# Version: 4.5.0
# Created on: 2024-04-22T21.26.12
# Steps:
#   picard_sam_to_fastq: 0 job... skipping
#   trimmomatic: 2 jobs
#   merge_trimmomatic_stats: 1 job
#   mapping_bwa_mem_sambamba: 2 jobs
#   sambamba_merge_bam_files: 2 jobs
#   sambamba_mark_duplicates: 3 jobs
#   sambamba_view_filter: 3 jobs
#   bedtools_blacklist_filter: 2 jobs
#   metrics: 5 jobs
#   homer_make_tag_directory: 2 jobs
#   qc_metrics: 1 job
#   homer_make_ucsc_file: 5 jobs
#   macs2_callpeak: 5 jobs
#   TOTAL: 33 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/project/6001942/chris11/20240422_PCa_ChIP_corepressor/output/chip-pipeline_PCA_corepressor-GRCh38
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/ChipSeq.chipseq.job_list.$TIMESTAMP
export CONFIG_FILES="/cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-4.5.0/pipelines/chipseq/chipseq.base.ini,/cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-4.5.0/pipelines/common_ini/cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini,input/cedar_customParameters.txt"
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

#-------------------------------------------------------------------------------
# STEP: trimmomatic
#-------------------------------------------------------------------------------
STEP=trimmomatic
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: trimmomatic_1_JOB_ID: trimmomatic.ChIP_LNCaP_DMSO_SMRT
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ChIP_LNCaP_DMSO_SMRT
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.ChIP_LNCaP_DMSO_SMRT.a64e33968b70d5b4c7b6ecded77f8df5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'trimmomatic.ChIP_LNCaP_DMSO_SMRT.a64e33968b70d5b4c7b6ecded77f8df5.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.39 && \
mkdir -p trim/LNCaP_DMSO_SMRT/SMRT && \
touch trim/LNCaP_DMSO_SMRT/SMRT && \
`cat > trim/LNCaP_DMSO_SMRT/SMRT/ChIP_LNCaP_DMSO_SMRT.trim.adapters.fa << END
>Prefix/1
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>Prefix/2
TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
END
` && \
java -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx20G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /project/6001942/chris11/20240422_PCa_ChIP_corepressor/raw/chipseq_PCa_corepressor/fastp_output/LNCaP_DMSO_SMRT_1.fastq.gz \
  /project/6001942/chris11/20240422_PCa_ChIP_corepressor/raw/chipseq_PCa_corepressor/fastp_output/LNCaP_DMSO_SMRT_2.fastq.gz \
  trim/LNCaP_DMSO_SMRT/SMRT/ChIP_LNCaP_DMSO_SMRT.trim.pair1.fastq.gz \
  trim/LNCaP_DMSO_SMRT/SMRT/ChIP_LNCaP_DMSO_SMRT.trim.single1.fastq.gz \
  trim/LNCaP_DMSO_SMRT/SMRT/ChIP_LNCaP_DMSO_SMRT.trim.pair2.fastq.gz \
  trim/LNCaP_DMSO_SMRT/SMRT/ChIP_LNCaP_DMSO_SMRT.trim.single2.fastq.gz \
  ILLUMINACLIP:trim/LNCaP_DMSO_SMRT/SMRT/ChIP_LNCaP_DMSO_SMRT.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/LNCaP_DMSO_SMRT/SMRT/ChIP_LNCaP_DMSO_SMRT.trim.log
trimmomatic.ChIP_LNCaP_DMSO_SMRT.a64e33968b70d5b4c7b6ecded77f8df5.mugqic.done
chmod 755 $COMMAND
trimmomatic_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 20G -c 6 -N 1    | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$trimmomatic_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: trimmomatic_2_JOB_ID: trimmomatic.ChIP_LNCaP_DMSO_NCOR1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ChIP_LNCaP_DMSO_NCOR1
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.ChIP_LNCaP_DMSO_NCOR1.55a5277313f9b42060f18f9d2b105923.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'trimmomatic.ChIP_LNCaP_DMSO_NCOR1.55a5277313f9b42060f18f9d2b105923.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.39 && \
mkdir -p trim/LNCaP_DMSO_NCOR1/NCOR1 && \
touch trim/LNCaP_DMSO_NCOR1/NCOR1 && \
`cat > trim/LNCaP_DMSO_NCOR1/NCOR1/ChIP_LNCaP_DMSO_NCOR1.trim.adapters.fa << END
>Prefix/1
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>Prefix/2
TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
END
` && \
java -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx20G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /project/6001942/chris11/20240422_PCa_ChIP_corepressor/raw/chipseq_PCa_corepressor/fastp_output/LNCaP_DMSO_NCOR1_1.fastq.gz \
  /project/6001942/chris11/20240422_PCa_ChIP_corepressor/raw/chipseq_PCa_corepressor/fastp_output/LNCaP_DMSO_NCOR1_2.fastq.gz \
  trim/LNCaP_DMSO_NCOR1/NCOR1/ChIP_LNCaP_DMSO_NCOR1.trim.pair1.fastq.gz \
  trim/LNCaP_DMSO_NCOR1/NCOR1/ChIP_LNCaP_DMSO_NCOR1.trim.single1.fastq.gz \
  trim/LNCaP_DMSO_NCOR1/NCOR1/ChIP_LNCaP_DMSO_NCOR1.trim.pair2.fastq.gz \
  trim/LNCaP_DMSO_NCOR1/NCOR1/ChIP_LNCaP_DMSO_NCOR1.trim.single2.fastq.gz \
  ILLUMINACLIP:trim/LNCaP_DMSO_NCOR1/NCOR1/ChIP_LNCaP_DMSO_NCOR1.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/LNCaP_DMSO_NCOR1/NCOR1/ChIP_LNCaP_DMSO_NCOR1.trim.log
trimmomatic.ChIP_LNCaP_DMSO_NCOR1.55a5277313f9b42060f18f9d2b105923.mugqic.done
chmod 755 $COMMAND
trimmomatic_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 20G -c 6 -N 1    | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$trimmomatic_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
STEP=merge_trimmomatic_stats
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: merge_trimmomatic_stats_1_JOB_ID: merge_trimmomatic_stats.
#-------------------------------------------------------------------------------
JOB_NAME=merge_trimmomatic_stats.
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$trimmomatic_2_JOB_ID
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats..9c133b2fed5937f045cb3d3d2c2c4588.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'merge_trimmomatic_stats..9c133b2fed5937f045cb3d3d2c2c4588.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/pandoc/2.16.2 && \
mkdir -p metrics && \
touch metrics && \

echo -e "Sample\tReadset\tMark Name\tRaw Paired Reads #\tSurviving Paired Reads #\tSurviving Paired Reads %" > metrics/trimReadsetTable.tsv && \
grep ^Input trim/LNCaP_DMSO_SMRT/SMRT/ChIP_LNCaP_DMSO_SMRT.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/LNCaP_DMSO_SMRT\tChIP_LNCaP_DMSO_SMRT\tSMRT\t\1\t\2/' | \
awk '{OFS="\t"; print $0, $5 / $4 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/LNCaP_DMSO_NCOR1/NCOR1/ChIP_LNCaP_DMSO_NCOR1.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/LNCaP_DMSO_NCOR1\tChIP_LNCaP_DMSO_NCOR1\tNCOR1\t\1\t\2/' | \
awk '{OFS="\t"; print $0, $5 / $4 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
cut -f1,3- metrics/trimReadsetTable.tsv | awk -F"\t" '{OFS="\t"; if (NR==1) {if ($3=="Raw Paired Reads #") {paired=1};print "Sample", "Mark Name", "Raw Reads #", "Surviving Reads #", "Surviving %"} else {if (paired) {$3=$3*2; $4=$4*2}; sample[$1$2]=$1; markname[$1$2]=$2; raw[$1$2]+=$3; surviving[$1$2]+=$4}}END{for (samplemark in raw){print sample[samplemark], markname[samplemark], raw[samplemark], surviving[samplemark], surviving[samplemark] / raw[samplemark] * 100}}' \
  > metrics/trimSampleTable.tsv && \
mkdir -p report && \
cp metrics/trimReadsetTable.tsv metrics/trimSampleTable.tsv report/ && \
trim_readset_table_md=`LC_NUMERIC=en_CA awk -F "\t" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----|-----|-----:|-----:|-----:"} else {print $1, $2, $3, sprintf("%\47d", $4), sprintf("%\47d", $5), sprintf("%.1f", $6)}}' metrics/trimReadsetTable.tsv` && \
pandoc \
  /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-4.5.0/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --template /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-4.5.0/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --variable trailing_min_quality=30 \
  --variable min_length=50 \
  --variable read_type=Paired \
  --variable trim_readset_table="$trim_readset_table_md" \
  --to markdown \
  > report/Illumina.merge_trimmomatic_stats.md
merge_trimmomatic_stats..9c133b2fed5937f045cb3d3d2c2c4588.mugqic.done
chmod 755 $COMMAND
merge_trimmomatic_stats_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=00:20:00 --mem-per-cpu 4000M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$merge_trimmomatic_stats_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$merge_trimmomatic_stats_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: mapping_bwa_mem_sambamba
#-------------------------------------------------------------------------------
STEP=mapping_bwa_mem_sambamba
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: mapping_bwa_mem_sambamba_1_JOB_ID: mapping_bwa_mem_sambamba.ChIP_LNCaP_DMSO_SMRT
#-------------------------------------------------------------------------------
JOB_NAME=mapping_bwa_mem_sambamba.ChIP_LNCaP_DMSO_SMRT
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID
JOB_DONE=job_output/mapping_bwa_mem_sambamba/mapping_bwa_mem_sambamba.ChIP_LNCaP_DMSO_SMRT.49ffdde5b33942989818b8c65e3af5fd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'mapping_bwa_mem_sambamba.ChIP_LNCaP_DMSO_SMRT.49ffdde5b33942989818b8c65e3af5fd.mugqic.done' > $COMMAND
module purge && \
module load mugqic/bwa/0.7.17 mugqic/sambamba/0.8.1 && \
mkdir -p alignment/LNCaP_DMSO_SMRT/SMRT/ChIP_LNCaP_DMSO_SMRT && \
touch alignment/LNCaP_DMSO_SMRT/SMRT/ChIP_LNCaP_DMSO_SMRT && \
bwa mem -K 100000000 -v 3 -t 12 -Y \
   \
  -R '@RG\tID:ChIP_LNCaP_DMSO_SMRT\tSM:LNCaP_DMSO_SMRT\tLB:LNCaP_DMSO_SMRT\tCN:McGill_Genome_Centre\tPL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/LNCaP_DMSO_SMRT/SMRT/ChIP_LNCaP_DMSO_SMRT.trim.pair1.fastq.gz \
  trim/LNCaP_DMSO_SMRT/SMRT/ChIP_LNCaP_DMSO_SMRT.trim.pair2.fastq.gz | \
sambamba view -S -f bam \
  /dev/stdin \
    | \
sambamba sort  \
  /dev/stdin \
  --tmpdir ${SLURM_TMPDIR} \
  --out alignment/LNCaP_DMSO_SMRT/SMRT/ChIP_LNCaP_DMSO_SMRT/ChIP_LNCaP_DMSO_SMRT.sorted.bam && \
sambamba index  \
  alignment/LNCaP_DMSO_SMRT/SMRT/ChIP_LNCaP_DMSO_SMRT/ChIP_LNCaP_DMSO_SMRT.sorted.bam \
  alignment/LNCaP_DMSO_SMRT/SMRT/ChIP_LNCaP_DMSO_SMRT/ChIP_LNCaP_DMSO_SMRT.sorted.bam.bai
mapping_bwa_mem_sambamba.ChIP_LNCaP_DMSO_SMRT.49ffdde5b33942989818b8c65e3af5fd.mugqic.done
chmod 755 $COMMAND
mapping_bwa_mem_sambamba_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu 4000M -c 12 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$mapping_bwa_mem_sambamba_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$mapping_bwa_mem_sambamba_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: mapping_bwa_mem_sambamba_2_JOB_ID: mapping_bwa_mem_sambamba.ChIP_LNCaP_DMSO_NCOR1
#-------------------------------------------------------------------------------
JOB_NAME=mapping_bwa_mem_sambamba.ChIP_LNCaP_DMSO_NCOR1
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID
JOB_DONE=job_output/mapping_bwa_mem_sambamba/mapping_bwa_mem_sambamba.ChIP_LNCaP_DMSO_NCOR1.428f636c10ff30465cf698c19c53b126.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'mapping_bwa_mem_sambamba.ChIP_LNCaP_DMSO_NCOR1.428f636c10ff30465cf698c19c53b126.mugqic.done' > $COMMAND
module purge && \
module load mugqic/bwa/0.7.17 mugqic/sambamba/0.8.1 && \
mkdir -p alignment/LNCaP_DMSO_NCOR1/NCOR1/ChIP_LNCaP_DMSO_NCOR1 && \
touch alignment/LNCaP_DMSO_NCOR1/NCOR1/ChIP_LNCaP_DMSO_NCOR1 && \
bwa mem -K 100000000 -v 3 -t 12 -Y \
   \
  -R '@RG\tID:ChIP_LNCaP_DMSO_NCOR1\tSM:LNCaP_DMSO_NCOR1\tLB:LNCaP_DMSO_NCOR1\tCN:McGill_Genome_Centre\tPL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/LNCaP_DMSO_NCOR1/NCOR1/ChIP_LNCaP_DMSO_NCOR1.trim.pair1.fastq.gz \
  trim/LNCaP_DMSO_NCOR1/NCOR1/ChIP_LNCaP_DMSO_NCOR1.trim.pair2.fastq.gz | \
sambamba view -S -f bam \
  /dev/stdin \
    | \
sambamba sort  \
  /dev/stdin \
  --tmpdir ${SLURM_TMPDIR} \
  --out alignment/LNCaP_DMSO_NCOR1/NCOR1/ChIP_LNCaP_DMSO_NCOR1/ChIP_LNCaP_DMSO_NCOR1.sorted.bam && \
sambamba index  \
  alignment/LNCaP_DMSO_NCOR1/NCOR1/ChIP_LNCaP_DMSO_NCOR1/ChIP_LNCaP_DMSO_NCOR1.sorted.bam \
  alignment/LNCaP_DMSO_NCOR1/NCOR1/ChIP_LNCaP_DMSO_NCOR1/ChIP_LNCaP_DMSO_NCOR1.sorted.bam.bai
mapping_bwa_mem_sambamba.ChIP_LNCaP_DMSO_NCOR1.428f636c10ff30465cf698c19c53b126.mugqic.done
chmod 755 $COMMAND
mapping_bwa_mem_sambamba_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu 4000M -c 12 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$mapping_bwa_mem_sambamba_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$mapping_bwa_mem_sambamba_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: sambamba_merge_bam_files
#-------------------------------------------------------------------------------
STEP=sambamba_merge_bam_files
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: sambamba_merge_bam_files_1_JOB_ID: symlink_readset_sample_bam.LNCaP_DMSO_SMRT.SMRT
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.LNCaP_DMSO_SMRT.SMRT
JOB_DEPENDENCIES=$mapping_bwa_mem_sambamba_1_JOB_ID
JOB_DONE=job_output/sambamba_merge_bam_files/symlink_readset_sample_bam.LNCaP_DMSO_SMRT.SMRT.26aded342af1320c06f8dd5031ecd5f1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'symlink_readset_sample_bam.LNCaP_DMSO_SMRT.SMRT.26aded342af1320c06f8dd5031ecd5f1.mugqic.done' > $COMMAND
mkdir -p alignment/LNCaP_DMSO_SMRT/SMRT && \
touch alignment/LNCaP_DMSO_SMRT/SMRT && \
ln -s -f \
  ChIP_LNCaP_DMSO_SMRT/ChIP_LNCaP_DMSO_SMRT.sorted.bam \
  alignment/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.bam && \
ln -s -f \
  ChIP_LNCaP_DMSO_SMRT/ChIP_LNCaP_DMSO_SMRT.sorted.bam.bai \
  alignment/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.bam.bai
symlink_readset_sample_bam.LNCaP_DMSO_SMRT.SMRT.26aded342af1320c06f8dd5031ecd5f1.mugqic.done
chmod 755 $COMMAND
sambamba_merge_bam_files_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=00:10:00 --mem-per-cpu 4000M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_merge_bam_files_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$sambamba_merge_bam_files_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: sambamba_merge_bam_files_2_JOB_ID: symlink_readset_sample_bam.LNCaP_DMSO_NCOR1.NCOR1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.LNCaP_DMSO_NCOR1.NCOR1
JOB_DEPENDENCIES=$mapping_bwa_mem_sambamba_2_JOB_ID
JOB_DONE=job_output/sambamba_merge_bam_files/symlink_readset_sample_bam.LNCaP_DMSO_NCOR1.NCOR1.01f2976ea1de22185421a87e6f5ab60b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'symlink_readset_sample_bam.LNCaP_DMSO_NCOR1.NCOR1.01f2976ea1de22185421a87e6f5ab60b.mugqic.done' > $COMMAND
mkdir -p alignment/LNCaP_DMSO_NCOR1/NCOR1 && \
touch alignment/LNCaP_DMSO_NCOR1/NCOR1 && \
ln -s -f \
  ChIP_LNCaP_DMSO_NCOR1/ChIP_LNCaP_DMSO_NCOR1.sorted.bam \
  alignment/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.sorted.bam && \
ln -s -f \
  ChIP_LNCaP_DMSO_NCOR1/ChIP_LNCaP_DMSO_NCOR1.sorted.bam.bai \
  alignment/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.sorted.bam.bai
symlink_readset_sample_bam.LNCaP_DMSO_NCOR1.NCOR1.01f2976ea1de22185421a87e6f5ab60b.mugqic.done
chmod 755 $COMMAND
sambamba_merge_bam_files_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=00:10:00 --mem-per-cpu 4000M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_merge_bam_files_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$sambamba_merge_bam_files_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: sambamba_mark_duplicates
#-------------------------------------------------------------------------------
STEP=sambamba_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: sambamba_mark_duplicates_1_JOB_ID: sambamba_mark_duplicates.LNCaP_DMSO_SMRT.SMRT
#-------------------------------------------------------------------------------
JOB_NAME=sambamba_mark_duplicates.LNCaP_DMSO_SMRT.SMRT
JOB_DEPENDENCIES=$sambamba_merge_bam_files_1_JOB_ID
JOB_DONE=job_output/sambamba_mark_duplicates/sambamba_mark_duplicates.LNCaP_DMSO_SMRT.SMRT.2cc14e261d1037722626b97a99976a36.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'sambamba_mark_duplicates.LNCaP_DMSO_SMRT.SMRT.2cc14e261d1037722626b97a99976a36.mugqic.done' > $COMMAND
module purge && \
module load mugqic/sambamba/0.8.1 && \
mkdir -p alignment/LNCaP_DMSO_SMRT/SMRT && \
touch alignment/LNCaP_DMSO_SMRT/SMRT && \
sambamba markdup -t 4 --sort-buffer-size=8192 --io-buffer-size=1024 \
  alignment/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.bam \
  --tmpdir ${SLURM_TMPDIR} \
  alignment/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.dup.bam
sambamba_mark_duplicates.LNCaP_DMSO_SMRT.SMRT.2cc14e261d1037722626b97a99976a36.mugqic.done
chmod 755 $COMMAND
sambamba_mark_duplicates_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu 4000M -c 4 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_mark_duplicates_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$sambamba_mark_duplicates_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: sambamba_mark_duplicates_2_JOB_ID: sambamba_mark_duplicates.LNCaP_DMSO_NCOR1.NCOR1
#-------------------------------------------------------------------------------
JOB_NAME=sambamba_mark_duplicates.LNCaP_DMSO_NCOR1.NCOR1
JOB_DEPENDENCIES=$sambamba_merge_bam_files_2_JOB_ID
JOB_DONE=job_output/sambamba_mark_duplicates/sambamba_mark_duplicates.LNCaP_DMSO_NCOR1.NCOR1.7edebd423934b4bcf898800c53b818dc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'sambamba_mark_duplicates.LNCaP_DMSO_NCOR1.NCOR1.7edebd423934b4bcf898800c53b818dc.mugqic.done' > $COMMAND
module purge && \
module load mugqic/sambamba/0.8.1 && \
mkdir -p alignment/LNCaP_DMSO_NCOR1/NCOR1 && \
touch alignment/LNCaP_DMSO_NCOR1/NCOR1 && \
sambamba markdup -t 4 --sort-buffer-size=8192 --io-buffer-size=1024 \
  alignment/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.sorted.bam \
  --tmpdir ${SLURM_TMPDIR} \
  alignment/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.sorted.dup.bam
sambamba_mark_duplicates.LNCaP_DMSO_NCOR1.NCOR1.7edebd423934b4bcf898800c53b818dc.mugqic.done
chmod 755 $COMMAND
sambamba_mark_duplicates_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu 4000M -c 4 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_mark_duplicates_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$sambamba_mark_duplicates_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: sambamba_mark_duplicates_3_JOB_ID: sambamba_mark_duplicates_report
#-------------------------------------------------------------------------------
JOB_NAME=sambamba_mark_duplicates_report
JOB_DEPENDENCIES=$sambamba_mark_duplicates_1_JOB_ID:$sambamba_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/sambamba_mark_duplicates/sambamba_mark_duplicates_report.d28efa532225145057f11c794b6760ec.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'sambamba_mark_duplicates_report.d28efa532225145057f11c794b6760ec.mugqic.done' > $COMMAND
mkdir -p report && \
cp \
  /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-4.5.0/bfx/report/ChipSeq.sambamba_mark_duplicates.md \
  report/ChipSeq.sambamba_mark_duplicates.md
sambamba_mark_duplicates_report.d28efa532225145057f11c794b6760ec.mugqic.done
chmod 755 $COMMAND
sambamba_mark_duplicates_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=04:00:00 --mem-per-cpu 4000M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_mark_duplicates_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$sambamba_mark_duplicates_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: sambamba_view_filter
#-------------------------------------------------------------------------------
STEP=sambamba_view_filter
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: sambamba_view_filter_1_JOB_ID: sambamba_view_filter.LNCaP_DMSO_SMRT.SMRT
#-------------------------------------------------------------------------------
JOB_NAME=sambamba_view_filter.LNCaP_DMSO_SMRT.SMRT
JOB_DEPENDENCIES=$sambamba_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/sambamba_view_filter/sambamba_view_filter.LNCaP_DMSO_SMRT.SMRT.1a104e1b6dfbdc38ac72ebd469758109.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'sambamba_view_filter.LNCaP_DMSO_SMRT.SMRT.1a104e1b6dfbdc38ac72ebd469758109.mugqic.done' > $COMMAND
module purge && \
module load mugqic/sambamba/0.8.1 && \
mkdir -p alignment/LNCaP_DMSO_SMRT/SMRT && \
touch alignment/LNCaP_DMSO_SMRT/SMRT && \
sambamba view -t 4 -f bam -F "not unmapped and not failed_quality_control and mapping_quality >= 20"  \
  alignment/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.dup.bam \
  -o alignment/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.dup.filtered.bam  && \
sambamba index  \
  alignment/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.dup.filtered.bam \
  alignment/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.dup.filtered.bam.bai
sambamba_view_filter.LNCaP_DMSO_SMRT.SMRT.1a104e1b6dfbdc38ac72ebd469758109.mugqic.done
chmod 755 $COMMAND
sambamba_view_filter_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu 4000M -c 4 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_view_filter_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$sambamba_view_filter_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: sambamba_view_filter_2_JOB_ID: sambamba_view_filter.LNCaP_DMSO_NCOR1.NCOR1
#-------------------------------------------------------------------------------
JOB_NAME=sambamba_view_filter.LNCaP_DMSO_NCOR1.NCOR1
JOB_DEPENDENCIES=$sambamba_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/sambamba_view_filter/sambamba_view_filter.LNCaP_DMSO_NCOR1.NCOR1.7491e19efd69910f4a87f086507bf1d9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'sambamba_view_filter.LNCaP_DMSO_NCOR1.NCOR1.7491e19efd69910f4a87f086507bf1d9.mugqic.done' > $COMMAND
module purge && \
module load mugqic/sambamba/0.8.1 && \
mkdir -p alignment/LNCaP_DMSO_NCOR1/NCOR1 && \
touch alignment/LNCaP_DMSO_NCOR1/NCOR1 && \
sambamba view -t 4 -f bam -F "not unmapped and not failed_quality_control and mapping_quality >= 20"  \
  alignment/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.sorted.dup.bam \
  -o alignment/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.sorted.dup.filtered.bam  && \
sambamba index  \
  alignment/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.sorted.dup.filtered.bam \
  alignment/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.sorted.dup.filtered.bam.bai
sambamba_view_filter.LNCaP_DMSO_NCOR1.NCOR1.7491e19efd69910f4a87f086507bf1d9.mugqic.done
chmod 755 $COMMAND
sambamba_view_filter_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu 4000M -c 4 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_view_filter_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$sambamba_view_filter_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: sambamba_view_filter_3_JOB_ID: sambamba_view_filter_report
#-------------------------------------------------------------------------------
JOB_NAME=sambamba_view_filter_report
JOB_DEPENDENCIES=$sambamba_view_filter_1_JOB_ID:$sambamba_view_filter_2_JOB_ID
JOB_DONE=job_output/sambamba_view_filter/sambamba_view_filter_report.bbb008b1f5430fc19d69da019005d621.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'sambamba_view_filter_report.bbb008b1f5430fc19d69da019005d621.mugqic.done' > $COMMAND
module purge && \
module load mugqic/pandoc/2.16.2 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-4.5.0/bfx/report/ChipSeq.sambamba_view_filter.md \
  --variable min_mapq="20" \
  /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-4.5.0/bfx/report/ChipSeq.sambamba_view_filter.md \
  > report/ChipSeq.sambamba_view_filter.md
sambamba_view_filter_report.bbb008b1f5430fc19d69da019005d621.mugqic.done
chmod 755 $COMMAND
sambamba_view_filter_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=04:00:00 --mem-per-cpu 4000M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_view_filter_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$sambamba_view_filter_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: bedtools_blacklist_filter
#-------------------------------------------------------------------------------
STEP=bedtools_blacklist_filter
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: bedtools_blacklist_filter_1_JOB_ID: bedtools_intersect.LNCaP_DMSO_SMRT.SMRT
#-------------------------------------------------------------------------------
JOB_NAME=bedtools_intersect.LNCaP_DMSO_SMRT.SMRT
JOB_DEPENDENCIES=$sambamba_view_filter_1_JOB_ID
JOB_DONE=job_output/bedtools_blacklist_filter/bedtools_intersect.LNCaP_DMSO_SMRT.SMRT.4624702e7b71a04752a2234d6bae04ef.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'bedtools_intersect.LNCaP_DMSO_SMRT.SMRT.4624702e7b71a04752a2234d6bae04ef.mugqic.done' > $COMMAND
module purge && \
module load mugqic/bedtools/2.30.0 mugqic/sambamba/0.8.1 && \
bedtools intersect -v -header \
  -a alignment/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.dup.filtered.bam \
  -b /cvmfs/soft.mugqic/CentOS6/genomes/blacklists/hg38.Kundaje.GRCh38_unified_Excludable.bed \
  > alignment/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.dup.filtered.cleaned.bam && \
sambamba index  \
  alignment/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.dup.filtered.cleaned.bam \
  alignment/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.dup.filtered.cleaned.bam.bai
bedtools_intersect.LNCaP_DMSO_SMRT.SMRT.4624702e7b71a04752a2234d6bae04ef.mugqic.done
chmod 755 $COMMAND
bedtools_blacklist_filter_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4000M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bedtools_blacklist_filter_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$bedtools_blacklist_filter_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: bedtools_blacklist_filter_2_JOB_ID: bedtools_intersect.LNCaP_DMSO_NCOR1.NCOR1
#-------------------------------------------------------------------------------
JOB_NAME=bedtools_intersect.LNCaP_DMSO_NCOR1.NCOR1
JOB_DEPENDENCIES=$sambamba_view_filter_2_JOB_ID
JOB_DONE=job_output/bedtools_blacklist_filter/bedtools_intersect.LNCaP_DMSO_NCOR1.NCOR1.bd85334922a86fa6b514ff9828cef8d1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'bedtools_intersect.LNCaP_DMSO_NCOR1.NCOR1.bd85334922a86fa6b514ff9828cef8d1.mugqic.done' > $COMMAND
module purge && \
module load mugqic/bedtools/2.30.0 mugqic/sambamba/0.8.1 && \
bedtools intersect -v -header \
  -a alignment/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.sorted.dup.filtered.bam \
  -b /cvmfs/soft.mugqic/CentOS6/genomes/blacklists/hg38.Kundaje.GRCh38_unified_Excludable.bed \
  > alignment/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.sorted.dup.filtered.cleaned.bam && \
sambamba index  \
  alignment/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.sorted.dup.filtered.cleaned.bam \
  alignment/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.sorted.dup.filtered.cleaned.bam.bai
bedtools_intersect.LNCaP_DMSO_NCOR1.NCOR1.bd85334922a86fa6b514ff9828cef8d1.mugqic.done
chmod 755 $COMMAND
bedtools_blacklist_filter_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4000M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bedtools_blacklist_filter_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$bedtools_blacklist_filter_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: metrics
#-------------------------------------------------------------------------------
STEP=metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: metrics_1_JOB_ID: picard_collect_multiple_metrics.LNCaP_DMSO_SMRT.SMRT
#-------------------------------------------------------------------------------
JOB_NAME=picard_collect_multiple_metrics.LNCaP_DMSO_SMRT.SMRT
JOB_DEPENDENCIES=$bedtools_blacklist_filter_1_JOB_ID
JOB_DONE=job_output/metrics/picard_collect_multiple_metrics.LNCaP_DMSO_SMRT.SMRT.3a654c518a5630a4076630a497e437a4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_collect_multiple_metrics.LNCaP_DMSO_SMRT.SMRT.3a654c518a5630a4076630a497e437a4.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.26.6 mugqic/R_Bioconductor/4.0.3_3.12 && \
mkdir -p metrics/LNCaP_DMSO_SMRT/SMRT && \
touch metrics/LNCaP_DMSO_SMRT/SMRT && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx8G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=${SLURM_TMPDIR} \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.dup.filtered.cleaned.bam \
 OUTPUT=metrics/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.dup.filtered.cleaned.all.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_collect_multiple_metrics.LNCaP_DMSO_SMRT.SMRT.3a654c518a5630a4076630a497e437a4.mugqic.done
chmod 755 $COMMAND
metrics_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=06:00:00 --mem-per-cpu 4000M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$metrics_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: metrics_2_JOB_ID: metrics_flagstat.LNCaP_DMSO_SMRT.SMRT
#-------------------------------------------------------------------------------
JOB_NAME=metrics_flagstat.LNCaP_DMSO_SMRT.SMRT
JOB_DEPENDENCIES=$sambamba_mark_duplicates_1_JOB_ID:$bedtools_blacklist_filter_1_JOB_ID
JOB_DONE=job_output/metrics/metrics_flagstat.LNCaP_DMSO_SMRT.SMRT.dd6f74bc666e980e132a38b50a19f85b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'metrics_flagstat.LNCaP_DMSO_SMRT.SMRT.dd6f74bc666e980e132a38b50a19f85b.mugqic.done' > $COMMAND
module purge && \
module load mugqic/sambamba/0.8.1 && \
mkdir -p metrics/LNCaP_DMSO_SMRT/SMRT && \
touch metrics/LNCaP_DMSO_SMRT/SMRT && \
sambamba flagstat  \
  alignment/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.dup.bam \
  > metrics/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.dup.flagstat && \
sambamba flagstat  \
  alignment/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.dup.filtered.cleaned.bam \
  > metrics/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.dup.filtered.cleaned.flagstat
metrics_flagstat.LNCaP_DMSO_SMRT.SMRT.dd6f74bc666e980e132a38b50a19f85b.mugqic.done
chmod 755 $COMMAND
metrics_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=06:00:00 --mem 8G -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$metrics_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: metrics_3_JOB_ID: picard_collect_multiple_metrics.LNCaP_DMSO_NCOR1.NCOR1
#-------------------------------------------------------------------------------
JOB_NAME=picard_collect_multiple_metrics.LNCaP_DMSO_NCOR1.NCOR1
JOB_DEPENDENCIES=$bedtools_blacklist_filter_2_JOB_ID
JOB_DONE=job_output/metrics/picard_collect_multiple_metrics.LNCaP_DMSO_NCOR1.NCOR1.584da01e21ad46185795a91fa2abb9eb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_collect_multiple_metrics.LNCaP_DMSO_NCOR1.NCOR1.584da01e21ad46185795a91fa2abb9eb.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.26.6 mugqic/R_Bioconductor/4.0.3_3.12 && \
mkdir -p metrics/LNCaP_DMSO_NCOR1/NCOR1 && \
touch metrics/LNCaP_DMSO_NCOR1/NCOR1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx8G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=${SLURM_TMPDIR} \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.sorted.dup.filtered.cleaned.bam \
 OUTPUT=metrics/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.sorted.dup.filtered.cleaned.all.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_collect_multiple_metrics.LNCaP_DMSO_NCOR1.NCOR1.584da01e21ad46185795a91fa2abb9eb.mugqic.done
chmod 755 $COMMAND
metrics_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=06:00:00 --mem-per-cpu 4000M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$metrics_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: metrics_4_JOB_ID: metrics_flagstat.LNCaP_DMSO_NCOR1.NCOR1
#-------------------------------------------------------------------------------
JOB_NAME=metrics_flagstat.LNCaP_DMSO_NCOR1.NCOR1
JOB_DEPENDENCIES=$sambamba_mark_duplicates_2_JOB_ID:$bedtools_blacklist_filter_2_JOB_ID
JOB_DONE=job_output/metrics/metrics_flagstat.LNCaP_DMSO_NCOR1.NCOR1.2f659b5bc57ad891ec5e4c840e7154a2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'metrics_flagstat.LNCaP_DMSO_NCOR1.NCOR1.2f659b5bc57ad891ec5e4c840e7154a2.mugqic.done' > $COMMAND
module purge && \
module load mugqic/sambamba/0.8.1 && \
mkdir -p metrics/LNCaP_DMSO_NCOR1/NCOR1 && \
touch metrics/LNCaP_DMSO_NCOR1/NCOR1 && \
sambamba flagstat  \
  alignment/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.sorted.dup.bam \
  > metrics/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.sorted.dup.flagstat && \
sambamba flagstat  \
  alignment/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.sorted.dup.filtered.cleaned.bam \
  > metrics/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.sorted.dup.filtered.cleaned.flagstat
metrics_flagstat.LNCaP_DMSO_NCOR1.NCOR1.2f659b5bc57ad891ec5e4c840e7154a2.mugqic.done
chmod 755 $COMMAND
metrics_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=06:00:00 --mem 8G -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$metrics_4_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: metrics_5_JOB_ID: metrics_report
#-------------------------------------------------------------------------------
JOB_NAME=metrics_report
JOB_DEPENDENCIES=$bedtools_blacklist_filter_1_JOB_ID:$bedtools_blacklist_filter_2_JOB_ID:$metrics_2_JOB_ID:$metrics_4_JOB_ID
JOB_DONE=job_output/metrics/metrics_report.8a670c7be95721050ff885c2e52addca.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'metrics_report.8a670c7be95721050ff885c2e52addca.mugqic.done' > $COMMAND
module purge && \
module load mugqic/pandoc/2.16.2 && \
module load mugqic/sambamba/0.8.1 && \
mkdir -p metrics
cp /dev/null metrics/SampleMetrics.tsv && \
declare -A samples_associative_array=(["LNCaP_DMSO_SMRT"]="SMRT" ["LNCaP_DMSO_NCOR1"]="NCOR1") && \
for sample in ${!samples_associative_array[@]}
do
  for mark_name in ${samples_associative_array[$sample]}
  do
    raw_flagstat_file=metrics/$sample/$mark_name/$sample.$mark_name.sorted.dup.flagstat
    filtered_flagstat_file=metrics/$sample/$mark_name/$sample.$mark_name.sorted.dup.filtered.flagstat
    bam_file=alignment/$sample/$mark_name/$sample.$mark_name.sorted.dup.filtered.cleaned.bam
    raw_supplementarysecondary_reads=`bc <<< $(grep "secondary" $raw_flagstat_file | sed -e 's/ + [[:digit:]]* secondary.*//')+$(grep "supplementary" $raw_flagstat_file | sed -e 's/ + [[:digit:]]* supplementary.*//')`
    mapped_reads=`bc <<< $(grep "mapped (" $raw_flagstat_file | sed -e 's/ + [[:digit:]]* mapped (.*)//')-$raw_supplementarysecondary_reads`
    filtered_supplementarysecondary_reads=`bc <<< $(grep "secondary" $filtered_flagstat_file | sed -e 's/ + [[:digit:]]* secondary.*//')+$(grep "supplementary" $filtered_flagstat_file | sed -e 's/ + [[:digit:]]* supplementary.*//')`
    filtered_reads=`bc <<< $(grep "in total" $filtered_flagstat_file | sed -e 's/ + [[:digit:]]* in total .*//')-$filtered_supplementarysecondary_reads`
    filtered_mapped_reads=`bc <<< $(grep "mapped (" $filtered_flagstat_file | sed -e 's/ + [[:digit:]]* mapped (.*)//')-$filtered_supplementarysecondary_reads`
    filtered_mapped_rate=`echo "scale=4; 100*$filtered_mapped_reads/$filtered_reads" | bc -l`
    filtered_dup_reads=`grep "duplicates" $filtered_flagstat_file | sed -e 's/ + [[:digit:]]* duplicates$//'`
    filtered_dup_rate=`echo "scale=4; 100*$filtered_dup_reads/$filtered_mapped_reads" | bc -l`
    filtered_dedup_reads=`echo "$filtered_mapped_reads-$filtered_dup_reads" | bc -l`
    if [[ -s metrics/trimSampleTable.tsv ]]
      then
        raw_reads=$(grep -P "${sample}\t${mark_name}" metrics/trimSampleTable.tsv | cut -f 3)
        raw_trimmed_reads=`bc <<< $(grep "in total" $raw_flagstat_file | sed -e 's/ + [[:digit:]]* in total .*//')-$raw_supplementarysecondary_reads`
        mapped_reads_rate=`echo "scale=4; 100*$mapped_reads/$raw_trimmed_reads" | bc -l`
        raw_trimmed_rate=`echo "scale=4; 100*$raw_trimmed_reads/$raw_reads" | bc -l`
        filtered_rate=`echo "scale=4; 100*$filtered_reads/$raw_trimmed_reads" | bc -l`
      else
        raw_reads=`bc <<< $(grep "in total" $raw_flagstat_file | sed -e 's/ + [[:digit:]]* in total .*//')-$raw_supplementarysecondary_reads`
        raw_trimmed_reads="NULL"
        mapped_reads_rate=`echo "scale=4; 100*$mapped_reads/$raw_reads" | bc -l`
        raw_trimmed_rate="NULL"
        filtered_rate=`echo "scale=4; 100*$filtered_reads/$raw_reads" | bc -l`
    fi
    filtered_mito_reads=$(sambamba view -F "not duplicate" -c $bam_file chrM)
    filtered_mito_rate=$(echo "scale=4; 100*$filtered_mito_reads/$filtered_mapped_reads" | bc -l)
    echo -e "$sample\t$mark_name\t$raw_reads\t$raw_trimmed_reads\t$raw_trimmed_rate\t$mapped_reads\t$mapped_reads_rate\t$filtered_reads\t$filtered_rate\t$filtered_mapped_reads\t$filtered_mapped_rate\t$filtered_dup_reads\t$filtered_dup_rate\t$filtered_dedup_reads\t$filtered_mito_reads\t$filtered_mito_rate" >> metrics/SampleMetrics.tsv
  done
done && \
sed -i -e "1 i\Sample\tMark Name\tRaw Reads #\tRemaining Reads after Trimming #\tRemaining Reads after Trimming %\tAligned Trimmed Reads #\tAligned Trimmed Reads %\tRemaining Reads after Filtering #\tRemaining Reads after Filtering %\tAligned Filtered Reads #\tAligned Filtered Reads %\tDuplicate Reads #\tDuplicate Reads %\tFinal Aligned Reads # without Duplicates\tMitochondrial Reads #\tMitochondrial Reads %" metrics/SampleMetrics.tsv && \
mkdir -p report && \
cp metrics/SampleMetrics.tsv report/SampleMetrics.tsv && \
sample_table=`LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:"} else {$1 = $1; print $0}}' report/SampleMetrics.tsv` && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-4.5.0/bfx/report/ChipSeq.metrics.md \
  --variable sample_table="$sample_table" \
  /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-4.5.0/bfx/report/ChipSeq.metrics.md \
  > report/ChipSeq.metrics.md

metrics_report.8a670c7be95721050ff885c2e52addca.mugqic.done
chmod 755 $COMMAND
metrics_5_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=03:00:00 --mem 16G -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$metrics_5_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: homer_make_tag_directory
#-------------------------------------------------------------------------------
STEP=homer_make_tag_directory
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_1_JOB_ID: homer_make_tag_directory.LNCaP_DMSO_SMRT.SMRT
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LNCaP_DMSO_SMRT.SMRT
JOB_DEPENDENCIES=$bedtools_blacklist_filter_1_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.LNCaP_DMSO_SMRT.SMRT.ab75ba3769bcf8ad5cf0428d836555c0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_tag_directory.LNCaP_DMSO_SMRT.SMRT.ab75ba3769bcf8ad5cf0428d836555c0.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.11 mugqic/samtools/1.14 && \
makeTagDirectory \
  tags/LNCaP_DMSO_SMRT/LNCaP_DMSO_SMRT.SMRT \
  alignment/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.dup.filtered.cleaned.bam \
  -genome hg38 \
  -checkGC \
 
homer_make_tag_directory.LNCaP_DMSO_SMRT.SMRT.ab75ba3769bcf8ad5cf0428d836555c0.mugqic.done
chmod 755 $COMMAND
homer_make_tag_directory_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 8G -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$homer_make_tag_directory_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$homer_make_tag_directory_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_2_JOB_ID: homer_make_tag_directory.LNCaP_DMSO_NCOR1.NCOR1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LNCaP_DMSO_NCOR1.NCOR1
JOB_DEPENDENCIES=$bedtools_blacklist_filter_2_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.LNCaP_DMSO_NCOR1.NCOR1.e06aeac269691be4dd266a74388c0155.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_tag_directory.LNCaP_DMSO_NCOR1.NCOR1.e06aeac269691be4dd266a74388c0155.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.11 mugqic/samtools/1.14 && \
makeTagDirectory \
  tags/LNCaP_DMSO_NCOR1/LNCaP_DMSO_NCOR1.NCOR1 \
  alignment/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.sorted.dup.filtered.cleaned.bam \
  -genome hg38 \
  -checkGC \
 
homer_make_tag_directory.LNCaP_DMSO_NCOR1.NCOR1.e06aeac269691be4dd266a74388c0155.mugqic.done
chmod 755 $COMMAND
homer_make_tag_directory_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 8G -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$homer_make_tag_directory_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$homer_make_tag_directory_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: qc_metrics
#-------------------------------------------------------------------------------
STEP=qc_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: qc_metrics_1_JOB_ID: qc_plots_R
#-------------------------------------------------------------------------------
JOB_NAME=qc_plots_R
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID:$homer_make_tag_directory_2_JOB_ID
JOB_DONE=job_output/qc_metrics/qc_plots_R.0e596d69f62d8ab42388b7790a2b8ee1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'qc_plots_R.0e596d69f62d8ab42388b7790a2b8ee1.mugqic.done' > $COMMAND
module purge && \
module load mugqic/mugqic_tools/2.12.4 mugqic/R_Bioconductor/4.0.3_3.12 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqGenerateQCMetrics.R \
  ../../raw/chipseq_PCa_corepressor/readset_chipseq_LNCaP_DMSO_PE_20240422.txt \
  /project/6001942/chris11/20240422_PCa_ChIP_corepressor/output/chip-pipeline_PCA_corepressor-GRCh38 && \
cp /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-4.5.0/bfx/report/ChipSeq.qc_metrics.md report/ChipSeq.qc_metrics.md && \
declare -A samples_associative_array=(["LNCaP_DMSO_SMRT"]="SMRT" ["LNCaP_DMSO_NCOR1"]="NCOR1") && \
for sample in ${!samples_associative_array[@]}
do
  for mark_name in ${samples_associative_array[$sample]}
  do
    cp --parents graphs/${sample}.${mark_name}_QC_Metrics.ps report/
    convert -rotate 90 graphs/${sample}.${mark_name}_QC_Metrics.ps report/graphs/${sample}.${mark_name}_QC_Metrics.png
    echo -e "----\n\n![QC Metrics for Sample $sample and Mark $mark_name ([download high-res image](graphs/${sample}.${mark_name}_QC_Metrics.ps))](graphs/${sample}.${mark_name}_QC_Metrics.png)\n" >> report/ChipSeq.qc_metrics.md
  done
done
qc_plots_R.0e596d69f62d8ab42388b7790a2b8ee1.mugqic.done
chmod 755 $COMMAND
qc_metrics_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu 4000M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$qc_metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$qc_metrics_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: homer_make_ucsc_file
#-------------------------------------------------------------------------------
STEP=homer_make_ucsc_file
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_1_JOB_ID: homer_make_ucsc_file.LNCaP_DMSO_SMRT.SMRT
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.LNCaP_DMSO_SMRT.SMRT
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.LNCaP_DMSO_SMRT.SMRT.b8fa5b52a274820384e98d90d18770d6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file.LNCaP_DMSO_SMRT.SMRT.b8fa5b52a274820384e98d90d18770d6.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.11 && \
mkdir -p tracks/LNCaP_DMSO_SMRT/SMRT && \
touch tracks/LNCaP_DMSO_SMRT/SMRT && \
makeUCSCfile \
  tags/LNCaP_DMSO_SMRT/LNCaP_DMSO_SMRT.SMRT \
  > tracks/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.ucsc.bedGraph && \
gzip -c tracks/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.ucsc.bedGraph \
  > tracks/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.ucsc.bedGraph.gz
homer_make_ucsc_file.LNCaP_DMSO_SMRT.SMRT.b8fa5b52a274820384e98d90d18770d6.mugqic.done
chmod 755 $COMMAND
homer_make_ucsc_file_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 32G -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$homer_make_ucsc_file_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$homer_make_ucsc_file_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_2_JOB_ID: homer_make_ucsc_file_bigWig.LNCaP_DMSO_SMRT.SMRT
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.LNCaP_DMSO_SMRT.SMRT
JOB_DEPENDENCIES=$homer_make_ucsc_file_1_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.LNCaP_DMSO_SMRT.SMRT.301fc7eecc9895bed772469c3983aad3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file_bigWig.LNCaP_DMSO_SMRT.SMRT.301fc7eecc9895bed772469c3983aad3.mugqic.done' > $COMMAND
module purge && \
module load mugqic/ucsc/v346 && \
mkdir -p tracks/LNCaP_DMSO_SMRT/SMRT/bigWig && \
touch tracks/LNCaP_DMSO_SMRT/SMRT/bigWig && \
(cat tracks/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.ucsc.bedGraph | head -n 1 > tracks/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.ucsc.bedGraph.head.tmp ; ec=$?; if [ "$ec" -eq 141 ]; then exit 0; else exit "$ec"; fi) && \
cat tracks/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -v "GL\|lambda\|pUC19\|KI\|\KN\|random"  | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.ucsc.bedGraph.body.tmp && \
cat tracks/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.ucsc.bedGraph.head.tmp tracks/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.ucsc.bedGraph.body.tmp > tracks/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.ucsc.bedGraph.sorted && \
rm tracks/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.ucsc.bedGraph.head.tmp tracks/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.ucsc.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/LNCaP_DMSO_SMRT/SMRT/bigWig/LNCaP_DMSO_SMRT.SMRT.bw
homer_make_ucsc_file_bigWig.LNCaP_DMSO_SMRT.SMRT.301fc7eecc9895bed772469c3983aad3.mugqic.done
chmod 755 $COMMAND
homer_make_ucsc_file_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4000M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$homer_make_ucsc_file_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$homer_make_ucsc_file_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_3_JOB_ID: homer_make_ucsc_file.LNCaP_DMSO_NCOR1.NCOR1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.LNCaP_DMSO_NCOR1.NCOR1
JOB_DEPENDENCIES=$homer_make_tag_directory_2_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.LNCaP_DMSO_NCOR1.NCOR1.1c292ddc6da3cfba1e2adadc24c1ab8c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file.LNCaP_DMSO_NCOR1.NCOR1.1c292ddc6da3cfba1e2adadc24c1ab8c.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.11 && \
mkdir -p tracks/LNCaP_DMSO_NCOR1/NCOR1 && \
touch tracks/LNCaP_DMSO_NCOR1/NCOR1 && \
makeUCSCfile \
  tags/LNCaP_DMSO_NCOR1/LNCaP_DMSO_NCOR1.NCOR1 \
  > tracks/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.ucsc.bedGraph && \
gzip -c tracks/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.ucsc.bedGraph \
  > tracks/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.ucsc.bedGraph.gz
homer_make_ucsc_file.LNCaP_DMSO_NCOR1.NCOR1.1c292ddc6da3cfba1e2adadc24c1ab8c.mugqic.done
chmod 755 $COMMAND
homer_make_ucsc_file_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 32G -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$homer_make_ucsc_file_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$homer_make_ucsc_file_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_4_JOB_ID: homer_make_ucsc_file_bigWig.LNCaP_DMSO_NCOR1.NCOR1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.LNCaP_DMSO_NCOR1.NCOR1
JOB_DEPENDENCIES=$homer_make_ucsc_file_3_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.LNCaP_DMSO_NCOR1.NCOR1.9b7bd6be84b5a3eb21fe9de51899d806.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file_bigWig.LNCaP_DMSO_NCOR1.NCOR1.9b7bd6be84b5a3eb21fe9de51899d806.mugqic.done' > $COMMAND
module purge && \
module load mugqic/ucsc/v346 && \
mkdir -p tracks/LNCaP_DMSO_NCOR1/NCOR1/bigWig && \
touch tracks/LNCaP_DMSO_NCOR1/NCOR1/bigWig && \
(cat tracks/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.ucsc.bedGraph | head -n 1 > tracks/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.ucsc.bedGraph.head.tmp ; ec=$?; if [ "$ec" -eq 141 ]; then exit 0; else exit "$ec"; fi) && \
cat tracks/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -v "GL\|lambda\|pUC19\|KI\|\KN\|random"  | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.ucsc.bedGraph.body.tmp && \
cat tracks/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.ucsc.bedGraph.head.tmp tracks/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.ucsc.bedGraph.body.tmp > tracks/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.ucsc.bedGraph.sorted && \
rm tracks/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.ucsc.bedGraph.head.tmp tracks/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.ucsc.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/LNCaP_DMSO_NCOR1/NCOR1/bigWig/LNCaP_DMSO_NCOR1.NCOR1.bw
homer_make_ucsc_file_bigWig.LNCaP_DMSO_NCOR1.NCOR1.9b7bd6be84b5a3eb21fe9de51899d806.mugqic.done
chmod 755 $COMMAND
homer_make_ucsc_file_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4000M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$homer_make_ucsc_file_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$homer_make_ucsc_file_4_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_5_JOB_ID: homer_make_ucsc_file_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_report
JOB_DEPENDENCIES=$homer_make_ucsc_file_1_JOB_ID:$homer_make_ucsc_file_3_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_report.6c18c57998ae36461b46104262e4b09e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file_report.6c18c57998ae36461b46104262e4b09e.mugqic.done' > $COMMAND
mkdir -p report && \
zip -r report/tracks.zip tracks/*/*/*.ucsc.bedGraph.gz && \
cp /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-4.5.0/bfx/report/ChipSeq.homer_make_ucsc_file.md report/
homer_make_ucsc_file_report.6c18c57998ae36461b46104262e4b09e.mugqic.done
chmod 755 $COMMAND
homer_make_ucsc_file_5_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4000M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$homer_make_ucsc_file_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$homer_make_ucsc_file_5_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# STEP: macs2_callpeak
#-------------------------------------------------------------------------------
STEP=macs2_callpeak
mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_1_JOB_ID: macs2_callpeak.LNCaP_DMSO_SMRT.SMRT
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LNCaP_DMSO_SMRT.SMRT
JOB_DEPENDENCIES=$bedtools_blacklist_filter_1_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LNCaP_DMSO_SMRT.SMRT.b6ec9de89fd9b27b68f13e68f2051ce5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'macs2_callpeak.LNCaP_DMSO_SMRT.SMRT.b6ec9de89fd9b27b68f13e68f2051ce5.mugqic.done' > $COMMAND
module purge && \
module load mugqic/python/3.7.3 mugqic/MACS2/2.2.7.1 && \
mkdir -p peak_call/LNCaP_DMSO_SMRT/SMRT && \
touch peak_call/LNCaP_DMSO_SMRT/SMRT && \
macs2 callpeak --format BAMPE --fix-bimodal  \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032 \
  --treatment \
  alignment/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.sorted.dup.filtered.cleaned.bam \
  --nolambda \
  --name peak_call/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT \
  >& peak_call/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT.diag.macs.out
macs2_callpeak.LNCaP_DMSO_SMRT.SMRT.b6ec9de89fd9b27b68f13e68f2051ce5.mugqic.done
chmod 755 $COMMAND
macs2_callpeak_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 32G -c 2 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$macs2_callpeak_1_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_2_JOB_ID: macs2_callpeak_bigBed.LNCaP_DMSO_SMRT.SMRT
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.LNCaP_DMSO_SMRT.SMRT
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.LNCaP_DMSO_SMRT.SMRT.433b8621985ae754d6b245029b631784.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'macs2_callpeak_bigBed.LNCaP_DMSO_SMRT.SMRT.433b8621985ae754d6b245029b631784.mugqic.done' > $COMMAND
module purge && \
module load mugqic/ucsc/v346 && \
awk '{if ($9 > 1000) {$9 = 1000}; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)}' peak_call/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT_peaks.narrowPeak > peak_call/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT_peaks.narrowPeak.bb
macs2_callpeak_bigBed.LNCaP_DMSO_SMRT.SMRT.433b8621985ae754d6b245029b631784.mugqic.done
chmod 755 $COMMAND
macs2_callpeak_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4000M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$macs2_callpeak_2_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_3_JOB_ID: macs2_callpeak.LNCaP_DMSO_NCOR1.NCOR1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LNCaP_DMSO_NCOR1.NCOR1
JOB_DEPENDENCIES=$bedtools_blacklist_filter_2_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LNCaP_DMSO_NCOR1.NCOR1.ea13f0cb2451351ec8f8714bf2ad52a6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'macs2_callpeak.LNCaP_DMSO_NCOR1.NCOR1.ea13f0cb2451351ec8f8714bf2ad52a6.mugqic.done' > $COMMAND
module purge && \
module load mugqic/python/3.7.3 mugqic/MACS2/2.2.7.1 && \
mkdir -p peak_call/LNCaP_DMSO_NCOR1/NCOR1 && \
touch peak_call/LNCaP_DMSO_NCOR1/NCOR1 && \
macs2 callpeak --format BAMPE --fix-bimodal  \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032 \
  --treatment \
  alignment/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.sorted.dup.filtered.cleaned.bam \
  --nolambda \
  --name peak_call/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1 \
  >& peak_call/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1.diag.macs.out
macs2_callpeak.LNCaP_DMSO_NCOR1.NCOR1.ea13f0cb2451351ec8f8714bf2ad52a6.mugqic.done
chmod 755 $COMMAND
macs2_callpeak_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem 32G -c 2 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$macs2_callpeak_3_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_4_JOB_ID: macs2_callpeak_bigBed.LNCaP_DMSO_NCOR1.NCOR1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.LNCaP_DMSO_NCOR1.NCOR1
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.LNCaP_DMSO_NCOR1.NCOR1.9cd5a1094de506dbb38b6ce10828f8ea.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'macs2_callpeak_bigBed.LNCaP_DMSO_NCOR1.NCOR1.9cd5a1094de506dbb38b6ce10828f8ea.mugqic.done' > $COMMAND
module purge && \
module load mugqic/ucsc/v346 && \
awk '{if ($9 > 1000) {$9 = 1000}; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)}' peak_call/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1_peaks.narrowPeak > peak_call/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.LNCaP_DMSO_NCOR1.NCOR1.9cd5a1094de506dbb38b6ce10828f8ea.mugqic.done
chmod 755 $COMMAND
macs2_callpeak_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4000M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$macs2_callpeak_4_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_5_JOB_ID: macs2_callpeak_report
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_report
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID:$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_report.f04a928a3a5fe70d28e3e1a16964166a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'macs2_callpeak_report.f04a928a3a5fe70d28e3e1a16964166a.mugqic.done' > $COMMAND
mkdir -p report && \
cp /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-4.5.0/bfx/report/ChipSeq.macs2_callpeak.md report/ && \
declare -A samples_associative_array=(["LNCaP_DMSO_SMRT"]="SMRT" ["LNCaP_DMSO_NCOR1"]="NCOR1") && \
for sample in ${!samples_associative_array[@]}
do
  for mark_name in ${samples_associative_array[$sample]}
  do
    cp -a --parents peak_call/$sample/$mark_name/ report/ && \
    echo -e "* [Peak Calls File for Sample $sample and Mark $mark_name](peak_call/$sample/$mark_name/${sample}.${mark_name}_peaks.xls)" >> report/ChipSeq.macs2_callpeak.md
  done
done
macs2_callpeak_report.f04a928a3a5fe70d28e3e1a16964166a.mugqic.done
chmod 755 $COMMAND
macs2_callpeak_5_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&      $COMMAND 
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE



if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu 4000M -c 1 -N 1    --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

echo "$macs2_callpeak_5_JOB_ID	$JOB_NAME submitted"
sleep 0.1

#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
LOG_MD5=$(echo $USER-'206.12.124.6-ChipSeq-LNCaP_DMSO_SMRT.ChIP_LNCaP_DMSO_SMRT,LNCaP_DMSO_NCOR1.ChIP_LNCaP_DMSO_NCOR1' | md5sum | awk '{ print $1 }')
if test -t 1; then ncolors=$(tput colors); if test -n "$ncolors" && test $ncolors -ge 8; then bold="$(tput bold)"; normal="$(tput sgr0)"; yellow="$(tput setaf 3)"; fi; fi
wget --quiet 'http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=cedar5.cedar.computecanada.ca&ip=206.12.124.6&pipeline=ChipSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,mapping_bwa_mem_sambamba,sambamba_merge_bam_files,sambamba_mark_duplicates,sambamba_view_filter,bedtools_blacklist_filter,metrics,homer_make_tag_directory,qc_metrics,homer_make_ucsc_file,macs2_callpeak&samples=2&md5=$LOG_MD5' -O /dev/null || echo "${bold}${yellow}Warning:${normal}${yellow} Genpipes ran successfully but was not send telemetry to mugqic.hpc.mcgill.ca. This error will not affect genpipes jobs you have submitted.${normal}"
