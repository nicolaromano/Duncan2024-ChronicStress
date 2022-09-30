#!/bin/bash

STAR_IND_DIR=/home/nico/Storage/Genomes/Mus_musculus.GRCm39.107_star_indices_over150/ 
# Change these to correct paths
FASTQ_DIR=xxxx
OUT_DIR=xxxx

SAMPLE_NAMES=(
    "CM-1" "CM-2" "CM-3" 
    "CS-1" "CS-2" "CS-3"
    "CSR-1" "CSR-2" "CSR-3")

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
