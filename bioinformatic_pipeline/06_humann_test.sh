#!/bin/bash
#SBATCH --job-name=humann_T2-9
#SBATCH --output=humann_T2-9.log
#SBATCH --error=humann_T2-9.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=420000M
#SBATCH --time=24:00:00
#SBATCH --partition= ##add
#SBATCH --account= ##add

#load module from server
module load humann/3.9

set -Eeuo pipefail

IN_DIR="/your/path/sickle_out/"
OUT_DIR="/your/path/sickle_out/humann_reads"
mkdir -p "$OUT_DIR"
shopt -s nullglob

# atteso: <sample>_R1.trimmed.fastq e <sample>_R2.trimmed.fastq

for R1 in "$IN_DIR"/*_R1.clean.fastq; do
  base=$(basename "$R1"); sample=${base%%_R1*}
  R2="$IN_DIR/${sample}_R2.clean.fastq"
  OUT="$OUT_DIR/${sample}.clean.concat.fastq.gz"

 if [[ ! -f "$R2" ]]; then
    echo "[WARN] manca R2 per $sample, salto."
    continue
  fi

  echo "[*] Unisco e comprimo $sample → $OUT"
  cat "$R1" "$R2" | gzip -1c > "$OUT"
done


# Percorsi
INPUT_FASTQ="/your/path/sickle_out/human_reads"  
OUTPUT_DIR="/your/path//Humann_output/"

# Crea output dir se non esiste
mkdir -p "$OUTPUT_DIR"

n_run=0
for fq in "$IN_DIR"/*.clean.concat.fastq.gz; do
  # questo è il prefisso che HUMAnN usa nei file di output
  sample_done=$(basename "$fq" .fastq.gz)  # es: T1-10_N2535.clean.concat
  # questo è un nome corto per i log/echo
  sample_short=$(basename "$fq" .clean.concat.fastq.gz)  # es: T1-10_N2535

  if [[ -f "$OUT_DIR/${sample_done}_genefamilies.tsv" && \
        -f "$OUT_DIR/${sample_done}_pathabundance.tsv" ]]; then
    echo "[SKIP] $sample_short già completato"
    continue
  fi

  echo "[RUN ] $(date '+%F %T') - HUMAnN: $sample"

# Comando HUMAnN
humann \
  --input "$fq" \
  --output "$OUT_DIR" \
  --threads 8 \
  --pathways unipathway \
  --remove-temp-output \

(( n_run++ ))
done

# Unisci tabelle solo se tutti i campioni hanno i file finali
all_ok=true
for fq in "$IN_DIR"/*.clean.concat.fastq.gz; do
  s=$(basename "$fq" .clean.concat.fastq.gz)
  if [[ ! ( -f "$OUT_DIR/${s}_genefamilies.tsv" && -f "$OUT_DIR/${s}_pathabundance.tsv" ) ]]; then
    all_ok=false; break
  fi
done

if $all_ok; then
  echo "[INFO] Tutti i campioni completi: unisco tabelle"
  humann_join_tables --input "$OUT_DIR" --file_name genefamilies  --output "$OUT_DIR/genefamilies.tsv"
  humann_join_tables --input "$OUT_DIR" --file_name pathabundance --output "$OUT_DIR/pathabundance.tsv"
  humann_renorm_table --input "$OUT_DIR/genefamilies.tsv"  --units cpm  --update-snames --output "$OUT_DIR/genefamilies_cpm.tsv"
  humann_renorm_table --input "$OUT_DIR/pathabundance.tsv" --units relab --update-snames --output "$OUT_DIR/pathabundance_rel.tsv"
  echo "[OK] Tabelle unite in $OUT_DIR"
else
  echo "[WARN] Non tutti i campioni sono completi. Salto la fase di join."
fi

