#!/bin/bash
#SBATCH --job-name=kraken2_classify
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem=420000M
#SBATCH --partition= #add
#SBATCH --account= #add
#SBATCH --time=2-00:00:00
#SBATCH --output=kraken2-%j.out
#SBATCH --error=kraken2-%j.err

# --- CONFIGURAZIONE (modifica queste variabili) ---
INPUT_DIR="/your/path//MEGAHIT_assemblies"	#path to assemblies obtained with megahit
OUTPUT_DIR="/your/path/Kraken_contings/" 
DB_PATH="/your/path/databases/kraken2/2025-02-04/k2_pluspfp_2025040/"    #path to database Kraken2
THREADS=8


# add kraken to the path 
export PATH=/data/kraken2:$PATH

# Vai nella cartella input
cd "$INPUT_DIR"

echo "[??] Inizio analisi contigs con Kraken2: $(date)"


# Loop su ogni file final.contigs.fa
for contig in "$INPUT_DIR"/*/final.contigs.fa; do
  sample=$(basename "$(dirname "$contig")")

 echo "[??] $(date '+%Y-%m-%d %H:%M:%S') - Campione: $sample"
 
 start=$(date +%s)

  kraken2 \
    --db "$DB_PATH" \
    --threads $THREADS \
    --paired \
    --report "${OUTPUT_DIR}/${base}.report" \
    --output "${OUTPUT_DIR}/${base}.kraken" \
    "$contig"

end=$(date +%s)
  runtime=$((end - start))
done

echo "[🎉] Tutte le analisi Kraken2 completate: $(date '+%Y-%m-%d %H:%M:%S')"
