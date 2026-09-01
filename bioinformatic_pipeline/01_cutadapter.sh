#!/bin/bash
#SBATCH --job-name=cutadapt_batch
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --partition= ##add
#SBATCH --account= ##add
#SBATCH --time=00-10:00:00
#SBATCH --output=cutadapt-%j.out
#SBATCH --error=cutadapt-%j.err



# === CONFIGURA QUI ===
INPUT_DIR="/your/path" ##change
OUTPUT_DIR="${INPUT_DIR}/cutadapt_out"
ADAPTER_SEQ="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
THREADS=4  # if you want,  Cutadapt works on 1 CPU

# === CREA CARTELLA OUTPUT ===
mkdir -p "$OUTPUT_DIR"

# === ATTIVA CONDA ENV ===
source $(conda info --base)/etc/profile.d/conda.sh
conda activate cutadapt_env

# Adattatori Illumina TruSeq (KAPA HyperPrep)        ###report multiqc is written "illumina universal adapters"
ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# === LOOP SU TUTTE LE COPPIE *_R1.fastq.gz ===
cd "$INPUT_DIR"
for r1 in *_1*.fq.gz; do			##CHECK HOW IT'S WRITTEN: _R1 OR _1!!!"##
    r2="${r1/_1/_2}"
    base="${r1%%_1*}"

    echo "[INFO] Processando: $r1 e $r2"

    cutadapt \
      -a "$ADAPTER_R1" \
      -A "$ADAPTER_R2" \
      -q 20,20 \
      -O 5 \
      -e 0.1 \
      --trim-n \
      -m 30 \
      --pair-filter=any \
      -o "${OUTPUT_DIR}/T1-6_N2531_R1.trim.fq.gz" \
      -p "${OUTPUT_DIR}/T1-6_N2531_R2.trim.fq.gz" \
      "/data/nanopore01/subsurfice/T1/raw/Raw_data/T1-6_N2531_R1.fastq" "/data/nanopore01/subsurfice/T1/raw/Raw_data/T1-6_N2531_R2.fastq.gz"
done

echo "[✅] all samples have been processed."
