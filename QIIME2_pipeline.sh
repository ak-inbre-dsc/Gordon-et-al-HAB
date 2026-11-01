#!/bin/bash
# ==============================================================================
# Supplementary Script: QIIME 2 Bioinformatic Pipeline for 18S V4-V5 rRNA
# Description: Automated processing of paired-end 18S V4-V5 amplicon sequencing
#              data including demultiplexing import, primer trimming, DADA2 
#              denoising, taxonomic classification (PR2), and phylogeny.
# ==============================================================================

# ------------------------------------------------------------------------------
# SLURM Resource Allocation Directives
# ------------------------------------------------------------------------------
#SBATCH --job-name=18SV4-5_qiime
#SBATCH --partition=analysis
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=18SV4-5_qiime_%j.out
#SBATCH --error=18SV4-5_qiime_%j.err

# Exit immediately if any command fails
set -e

# ------------------------------------------------------------------------------
# Environment Setup
# ------------------------------------------------------------------------------
# Initialize Conda and activate the target QIIME 2 environment
eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2025.7

# Set thread utilization for parallelized steps
N_THREADS=${SLURM_CPUS_PER_TASK:-12}

# ------------------------------------------------------------------------------
# Path Definitions & Parameter Configuration (Relative Paths)
# ------------------------------------------------------------------------------
# Set base directory (defaults to current working directory if not supplied as $1)
BASE_DIR="${1:-.}"

# Input files and directories
INPUT_DIR="${BASE_DIR}/raw_data/18S_V4-5"
CLASSIFIER="${BASE_DIR}/classifiers/pr2-v4-classifier.qza"
METADATA="${BASE_DIR}/metadata_v4-5.tsv"

# Output directories
OUTPUT_DIR="${BASE_DIR}/qiime2/18SV4-5"
OUTPUT_DIR_VIS="${OUTPUT_DIR}/visualization"

# Create output directories if they do not exist
mkdir -p "$OUTPUT_DIR" "$OUTPUT_DIR_VIS"

# ------------------------------------------------------------------------------
# 1. Data Import
# ------------------------------------------------------------------------------
# Import raw demultiplexed, paired-end FASTQ files (Casava 1.8 format)
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$INPUT_DIR" \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path "$OUTPUT_DIR/paired-end-demux.qza"

# ------------------------------------------------------------------------------
# 2. Primer Trimming (Cutadapt)
# ------------------------------------------------------------------------------
# Trim 18S V4-V5 specific forward and reverse primer sequences
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences "$OUTPUT_DIR/paired-end-demux.qza" \
  --p-front-f "AAAAAGCTTCCATCCAACATCTCAGCATGATGAAA" \
  --p-front-r "AAACTGCAGCCCCTCAGAATGATATTTGTCCTCA" \
  --p-cores "$N_THREADS" \
  --o-trimmed-sequences "$OUTPUT_DIR/trimmed-demux.qza"

# ------------------------------------------------------------------------------
# 3. Quality Control Visualization
# ------------------------------------------------------------------------------
# Generate interactive sequence quality summaries before and after primer trimming
qiime demux summarize \
  --i-data "$OUTPUT_DIR/paired-end-demux.qza" \
  --o-visualization "$OUTPUT_DIR_VIS/paired-end-demux.qzv"

qiime demux summarize \
  --i-data "$OUTPUT_DIR/trimmed-demux.qza" \
  --o-visualization "$OUTPUT_DIR_VIS/trimmed-demux.qzv"

# ------------------------------------------------------------------------------
# 4. Denoising & Amplicon Sequence Variant (ASV) Inference (DADA2)
# ------------------------------------------------------------------------------
# Quality filtering, error modeling, dereplication, paired-end merging, and chimera removal
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$OUTPUT_DIR/trimmed-demux.qza" \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 220 \
  --p-trunc-len-r 220 \
  --p-n-threads "$N_THREADS" \
  --o-table "$OUTPUT_DIR/table.qza" \
  --o-representative-sequences "$OUTPUT_DIR/rep-seqs.qza" \
  --o-denoising-stats "$OUTPUT_DIR/denoising-stats.qza"

# Summarize feature table and sequence metadata
qiime feature-table summarize \
  --i-table "$OUTPUT_DIR/table.qza" \
  --o-visualization "$OUTPUT_DIR_VIS/table.qzv" \
  --m-sample-metadata-file "$METADATA"

# ------------------------------------------------------------------------------
# 5. Taxonomic Classification (PR2 Database)
# ------------------------------------------------------------------------------
# Assign taxonomy to ASVs using Naive Bayes classifier trained on PR2 database
qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER" \
  --i-reads "$OUTPUT_DIR/rep-seqs.qza" \
  --p-n-jobs "$N_THREADS" \
  --o-classification "$OUTPUT_DIR/taxonomy.qza"

# Tabulate taxonomy assignment metadata
qiime metadata tabulate \
  --m-input-file "$OUTPUT_DIR/taxonomy.qza" \
  --o-visualization "$OUTPUT_DIR_VIS/taxonomy.qzv"

# Generate interactive taxa abundance bar plots
qiime taxa barplot \
  --i-table "$OUTPUT_DIR/table.qza" \
  --i-taxonomy "$OUTPUT_DIR/taxonomy.qza" \
  --m-metadata-file "$METADATA" \
  --o-visualization "$OUTPUT_DIR_VIS/taxa-barplot.qzv"

# ------------------------------------------------------------------------------
# 6. Phylogenetic Tree Construction
# ------------------------------------------------------------------------------
# Align representative sequences with MAFFT, mask hypervariable sites, and build FastTree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "$OUTPUT_DIR/rep-seqs.qza" \
  --p-n-threads "$N_THREADS" \
  --o-alignment "$OUTPUT_DIR/aligned-rep-seqs.qza" \
  --o-masked-alignment "$OUTPUT_DIR/masked-aligned-rep-seqs.qza" \
  --o-tree "$OUTPUT_DIR/unrooted-tree.qza" \
  --o-rooted-tree "$OUTPUT_DIR/rooted-tree.qza"