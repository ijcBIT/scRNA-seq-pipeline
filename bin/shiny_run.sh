#!/bin/bash

#SBATCH --job-name=shinyApp
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --time=12:00:00
#SBATCH --mem=110G


singularity exec \
  --bind /ijc/LABS/MERKEL/DATA/PROJECTS/mvilardell/:/scRNAseq_data \
  scRNAseq_refData.sif \
  Rscript -e "shiny::runApp('/scRNAseq_data/singularity/app.R')"
