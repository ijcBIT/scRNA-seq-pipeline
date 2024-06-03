#!/bin/bash

#SBATCH --job-name=scRNAseq_analysis
#SBATCH --output=scRNAseq_analysis.out
#SBATCH --error=scRNAseq_analysis.err
#SBATCH --time=12:00:00
#SBATCH --mem=110G

module load singularity

# Step 1: Run cellranger
echo "Running cellranger"

# Run the cellranger count script inside the Singularity container
singularity exec -C \
  --bind /ijc/LABS/MERKEL/RAW/NGS/scRNAseq_data:/scRNAseq_data \
  --bind /ijc/LABS/MERKEL/DATA/PROJECTS/mmerono/singularity:/scRNAseq_src \
  /ijc/LABS/MERKEL/DATA/PROJECTS/mmerono/singularity/scRNAseq.sif \
  /scRNAseq_src/cellranger_count.slm

# Check if cellranger was successful
if [ $? -ne 0 ]; then
  echo "cellranger failed."
  exit 1
fi


# Step 2: Generate HTML/pdf
echo "Generating HTML and Pdf"

singularity exec -C \
  --bind /ijc/LABS/MERKEL/RAW/NGS/scRNAseq_data:/scRNAseq_data \
  --bind /ijc/LABS/MERKEL/DATA/PROJECTS/mmerono/singularity:/scRNAseq_src \
  /ijc/LABS/MERKEL/DATA/PROJECTS/mmerono/singularity/scRNAseq.sif \
    Rscript -e "rmarkdown::render('/scRNAseq_src/scRNAseq.Rmd', 
                                output_file = '/scRNAseq_src/scRNAseq_singularity_output.html', 
                                params = list(file_10X_h5 = '/scRNAseq_data/filtered_feature_bc_matrix.h5'))"

# Check if HTML/Pdf generation was successful
if [ $? -ne 0 ]; then
  echo "HTML generation failed."
  exit 1
fi


# Step 3: Run the Shiny app
echo "Running Shiny app"

APP_DIR=/ijc/LABS/MERKEL/DATA/PROJECTS/mmerono/singularity
SIF_FILE=$APP_DIR/scRNAseq.sif
APP_FILE=$APP_DIR/app.R

singularity exec $SIF_FILE Rscript -e "shiny::runApp('$APP_FILE', host='0.0.0.0', port=3838)"


# Check if Shiny app was successful
if [ $? -ne 0 ]; then
  echo "Shiny app failed."
  exit 1
fi



echo "All the steps completed successfully!"
