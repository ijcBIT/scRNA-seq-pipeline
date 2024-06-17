#!/bin/bash

# Record start time
start=$(date +%s)

# Load environment variables 
source /scRNAseq_data/samples.contig

# Create the output directory file 
mkdir -p ${OUTPUT_DIR}

cd ${OUTPUT_DIR}

#######################################
####### cellranger parameters #########
#######################################

# Set the paths
transcriptome="/opt/refdata-gex-GRCh38-2024-A"

# Cell Ranger count command with increased parallelism
cellranger count --id=$(basename ${OUTPUT_DIR}) \
                 --transcriptome=${transcriptome} \
                 --create-bam=true \
                 --fastqs=${FASTQS} \
                 --sample=${SAMPLE_NAME} \
                 --localcores=8 \
                 --localmem=128


# Record end time, total run time, and print it
end=$(date +%s)
runtime=$((end-start))
echo $runtime



