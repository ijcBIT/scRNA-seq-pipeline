#!/bin/bash


#SBATCH --job-name=build_out
#SBATCH --time=20:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8


# Load Singularity module 
module load singularity


# build the image and name it scRNAseq.sif
singularity build --remote --force scRNAseq_refData.sif scRNAseq_refData.def 2>&1 | tee singularity_build.log
