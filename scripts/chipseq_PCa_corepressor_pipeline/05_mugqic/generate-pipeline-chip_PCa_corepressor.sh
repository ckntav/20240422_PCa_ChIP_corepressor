mkdir -p $SCRATCH/20240422_PCa_ChIP_corepressor/output/chip-pipeline_PCA_corepressor-GRCh38

$MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.py -j slurm -s 1-13 \
    -l debug \
    -r raw/chipseq_PCa_corepressor/readset_chipseq_LNCaP_DMSO_PE_20240422.txt \
    -o output/chip-pipeline_PCA_corepressor-GRCh38 \
    --config $MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.base.ini \
        $MUGQIC_PIPELINES_HOME/pipelines/common_ini/cedar.ini \
        $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini \
        input/cedar_customParameters.txt