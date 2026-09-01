#!/bin/bash
#SBATCH --job-name=sickle
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --partition= ##add
#SBATCH --account= #add
#SBATCH --time=01-00:00:00
#SBATCH --output=sickle-%j.out
#SBATCH --error=sickle-%j.err

# === metti sickle nel PATH  ===
echo 'export PATH="$HOME/tools/sickle:$PATH"' >> ~/.bashrc
source ~/.bashrc


# === LOOP PAIRED-END ===
INPUT_DIR="/your/path/cutadapt_out/"		#seq after cutadapater
OUTPUT_DIR="your/path/sickle_out"
mkdir -p "$OUTPUT_DIR"

cd "$INPUT_DIR"
for r1 in *_R1*.fq.gz; do
    r2="${r1/_R1/_R2}"
    base="${r1%%_R1*}"
    sickle pe -f "$r1" -r "$r2" -t sanger -o "$OUTPUT_DIR/${base}_R1.clean.fastq." -p "$OUTPUT_DIR/${base}_R2.clean.fastq" -s "$OUTPUT_DIR/${base}_singles.fastq" -q 15 -l 35 
done
