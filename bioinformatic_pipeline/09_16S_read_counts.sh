#!/bin/bash
#SBATCH --job-name=cov16S
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=00-02:00:00
#SBATCH --partition=bru-medium
#SBATCH --account=brusetti-lab
#SBATCH --output=cov16S-%j.out
#SBATCH --error=cov16S-%j.err

############################################
# CONFIGURATION
############################################

GFF_DIR="/your/path/output_16S"
BAM_DIR="/your/path/read_mapping/bam"
OUT_DIR="/your/path/16S_coverage"

############################################
# MODULES
############################################

module purge
module load samtools/1.19.2-gcc-13.3.0-a2yhwkt
module load bedtools2/2.31.1-gcc-12.1.0

mkdir -p "$OUT_DIR"

echo "[INFO] GFF_DIR = $GFF_DIR"
echo "[INFO] BAM_DIR = $BAM_DIR"
echo "[INFO] OUT_DIR = $OUT_DIR"

############################################
# LOOP OVER 16S GFF FILES
############################################

for gff in "$GFF_DIR"/*_16S.gff; do
    [[ -e "$gff" ]] || {
        echo "[WARN] No *_16S.gff files found"
        break
    }

    base=$(basename "$gff")
    sample=${base%_16S.gff}
    bam="$BAM_DIR/${sample}.sorted.bam"

    if [[ ! -f "$bam" ]]; then
        echo "[WARN] BAM missing for $sample: $bam" >&2
        continue
    fi

    echo "[INFO] Sample $sample"

    bed="$OUT_DIR/${sample}_16S.bed"
    counts="$OUT_DIR/${sample}_16S_counts.tsv"

    # Extract 16S rRNA features from GFF and convert to BED
    awk 'BEGIN{OFS="\t"}
    $0 !~ /^#/ && $3=="rRNA" && $9 ~ /16S_rRNA/ {
        print $1, $4-1, $5, $9
    }' "$gff" > "$bed"

    if [[ ! -s "$bed" ]]; then
        echo "[WARN] No 16S_rRNA feature found in $gff, skipping sample." >&2
        continue
    fi

    {
        echo -e "Sample\tContig\tStart\tEnd\tLength_bp\tRead_count"
        bedtools coverage -a "$bed" -b "$bam" -counts \
          | awk -v S="$sample" 'BEGIN{OFS="\t"}
            {
                len=$3-$2;
                print S,$1,$2+1,$3,len,$NF
            }'
    } > "$counts"

    echo "[INFO] Written $counts"
done

echo "[INFO] 16S coverage completed."
