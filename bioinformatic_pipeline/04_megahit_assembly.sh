#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --partition= #add
#SBATCH --account=#add
#SBATCH --time=3-00:00:00
#SBATCH --output=megahit-%j.out
#SBATCH --error=megahit-%j.err

module load megahit/1.2.9

INPUT_DIR="/you/path/sickle_out"
OUTPUT_DIR="/you/path/MEGAHIT_assemblies"

mkdir -p "$OUTPUT_DIR"

echo "===== Inizio assemblaggi paired-end ====="

for R1 in ${INPUT_DIR}/*_R1.clean.fastq; do
    SAMPLE=$(basename "$R1" _R1.clean.fastq)
    R2=${INPUT_DIR}/${SAMPLE}_R2.clean.fastq
    OUTDIR="${OUTPUT_DIR}/${SAMPLE}"

    # ✅ SALTA se sample già fatto
    if [[ -f "${OUTDIR}/final.contigs.fa" ]]; then
        echo "✅ Assemblaggio già presente per $SAMPLE. Skipping..."
        continue
    fi

    echo "🔹 Assembling sample: $SAMPLE"

        megahit \
            -1 "$R1" \
            -2 "$R2" \
            -o "$OUTDIR" \
            --min-contig-len 1000 \
            --num-cpu-threads 8
done

echo "===== Assemblaggi individuali completati ====="

##Statistiche N50, contig , etc 

echo "===== Calcolo statistiche di assembly ====="

STATS_FILE=${OUTPUT_DIR}/assembly_stats.tsv
echo -e "Sample\tNum_contigs\tTotal_length(bp)\tN50(bp)" > "$STATS_FILE"

for fa in ${OUTPUT_DIR}/*/final.contigs.fa; do
    SAMPLE=$(basename $(dirname "$fa"))

    STATS=$(awk '
        /^>/ {if (seqlen){print seqlen}; seqlen=0; next}
        {seqlen += length($0)}
        END {print seqlen}' "$fa" | sort -nr |
        awk 'BEGIN{sum=0; total=0}
             {a[NR]=$1; total+=a[NR]}
             END{
                half=total/2;
                for(i=1;i<=NR;i++){sum+=a[i]; if(sum>=half){print NR"\t"total"\t"a[i]; exit}}
             }')

    NUM=$(echo "$STATS" | cut -f1)
    TOTAL=$(echo "$STATS" | cut -f2)
    N50=$(echo "$STATS" | cut -f3)

    echo -e "${SAMPLE}\t${NUM}\t${TOTAL}\t${N50}" >> "$STATS_FILE"
done

echo "===== Statistiche salvate in: $STATS_FILE ====="
echo "===== Pipeline completata con successo ====="
