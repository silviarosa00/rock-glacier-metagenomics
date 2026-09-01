#!/bin/bash
#SBATCH --job-name=fastqc_post_sickle
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --partition=bru-medium
#SBATCH --account=brusetti-lab
#SBATCH --time=04:00:00
#SBATCH --output=fastqc-%j.out
#SBATCH --error=fastqc-%j.err

# === ATTIVA CONDA se serve ===
module load fastqc/0.12.1-gcc-12.1.0

INPUT_DIR="/data/nanopore01/subsurfice/T2/DATA/Raw_data/sickle_out/"
OUTPUT_DIR="/data/nanopore01/subsurfice/T2/DATA/Raw_data/sickle_out/fastqc_results"
mkdir -p "$OUTPUT_DIR"

for fq in "$INPUT_DIR"/*.fastq; do
    [ -e "$fq" ] || continue
    echo "[INFO] Analizzo: $(basename "$fq")"
    fastqc "$fq" -o "$OUTPUT_DIR"
done

echo "[✅] FastQC post-Sickle completato."
