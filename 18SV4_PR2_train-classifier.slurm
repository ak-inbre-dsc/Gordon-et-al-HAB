#!/bin/bash
#SBATCH --partition=analysis
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=18S_PR2_Train_%j.out

# Exit on error
set -e

# 1. Environment Setup
eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2025.7 

# Set directories
BASE_DIR="${1:-.}"
OUTPUT_DIR="${BASE_DIR%/}/classifiers"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Configuration
N_JOBS=24

# 18SV4 Primers
V4_FWD="CCAGCASCYGCGGTAATTCC"
V4_REV="ACTTTCGTTCTTGATYRA"

echo "Step 1: Downloading PR2 Database (v5.0.0)..."
qiime rescript get-pr2-data \
  --p-version '5.0.0' \
  --o-pr2-sequences pr2-raw-seqs.qza \
  --o-pr2-taxonomy pr2-raw-tax.qza

echo "Step 2: Culling sequences (Removing low quality)..."
qiime rescript cull-seqs \
  --i-sequences pr2-raw-seqs.qza \
  --o-clean-sequences pr2-culled-seqs.qza

echo "Step 3: Cleaning Taxonomy and Dereplicating..."
qiime rescript dereplicate \
  --i-sequences pr2-culled-seqs.qza \
  --i-taxa pr2-raw-tax.qza \
  --p-mode 'uniq' \
  --o-dereplicated-sequences pr2-final-seqs.qza \
  --o-dereplicated-taxa pr2-final-tax.qza

echo "Step 4a: Extracting V4 region..."
qiime feature-classifier extract-reads \
  --i-sequences pr2-final-seqs.qza \
  --p-f-primer "$V4_FWD" \
  --p-r-primer "$V4_REV" \
  --p-min-length 200 \
  --p-max-length 500 \
  --p-read-orientation 'forward' \
  --o-reads pr2-v4-extracts.qza

echo "Step 5a: Training V4 Naive Bayes..."
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads pr2-v4-extracts.qza \
  --i-reference-taxonomy pr2-final-tax.qza \
  --o-classifier pr2-v4-classifier.qza

echo "---------------------------------------------------------"
echo "18SV4  PR2 Training Complete at $(date)"
