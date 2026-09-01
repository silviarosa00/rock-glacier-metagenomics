#!/bin/bash
#SBATCH --job-name=humann_T2-9
#SBATCH --output=humann_T2-9.log
#SBATCH --error=humann_T2-9.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=420000M
#SBATCH --time=24:00:00
#SBATCH --partition=bru-medium
#SBATCH --account=brusetti-lab

# Carica il modulo 
module load humann/3.9


# Percorsi
INPUT_FASTQ="/data/nanopore01/subsurfice/T2/DATA/Raw_data/sickle_out/T2-9_combined.clean.fastq"
OUTPUT_DIR="/data/nanopore01/subsurfice/T2/Humann_output/"

# Crea output dir se non esiste
mkdir -p "$OUTPUT_DIR"

# Comando HUMAnN
humann \
  --input "$INPUT_FASTQ" \
  --output "$OUTPUT_DIR" \
  --threads 8 \
  --pathways unipathway

