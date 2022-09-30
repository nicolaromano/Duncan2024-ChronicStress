#!/bin/bash

STAR_IND_DIR=/home/nico/Storage/Genomes/Mus_musculus.GRCm39.107_star_indices_over150/ 
# Change these to correct paths
FASTQ_DIR=xxxx
OUT_DIR=xxxx

SAMPLE_NAMES=(
    "A1" "A2" "A3" 
    "B1" "B2" "B3" "B4" 
    "C1" "C2" "C3")

# Loop through the sample names
for SAMPLE_NAME in "${SAMPLE_NAMES[@]}"
do
    echo "Aligning sample $SAMPLE_NAME"
    echo ""

    STAR --genomeDir $STAR_IND_DIR --runThreadN 16 \
    --readFilesIn ${FASTQ_DIR}/${SAMPLE_NAME}_R1_001.fastq.gz ${FASTQ_DIR}/${SAMPLE_NAME}_R2_001.fastq.gz \
    --readFilesCommand zcat \
    --sjdbOverhang 150 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --outFileNamePrefix ${OUT_DIR}/$SAMPLE_NAME

done
