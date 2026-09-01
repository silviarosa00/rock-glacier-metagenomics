#!/bin/bash
#SBATCH --job-name=ARG_counts
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --partition=bru-medium
#SBATCH --account=brusetti-lab
#SBATCH --time=02-00:00:00
#SBATCH --output=ARG_counts-%j.out
#SBATCH --error=ARG_counts-%j.err

set -euo pipefail

############################################
# MODULES / ENVIRONMENT
############################################

module purge
module load samtools/1.19.2-gcc-13.3.0-a2yhwkt
module load bedtools2/2.31.1-gcc-12.1.0

echo "[INFO] samtools: $(samtools --version | head -n 1)"
echo "[INFO] bedtools: $(bedtools --version)"

############################################
# CONFIGURATION
############################################

THREADS=8
CHUNK_LINES=50000
CHECK_N=10

RGI_DIR="/your/path/output_rgi"
BAM_DIR="/your/path/ARG_abundance/bam"
BED_DIR="/your/path/ARG_abundance/bed_fixed"
COV_DIR="/your/path/ARG_abundance/coverage"

mkdir -p "$BED_DIR" "$COV_DIR"

############################################
# FUNCTIONS
############################################

make_bed_fixed() {
    local rgi_txt="$1"
    local bed_out="$2"

    awk -F'\t' 'BEGIN{OFS="\t"}
    NR==1{next}
    {
        orf=$1; aro=$11; hit=$9;

        split(orf,p,/ /);
        orf_name=p[1];
        contig=orf_name;
        sub(/_[0-9]+$/,"",contig);

        n=split(orf,a,/[[:space:]]*#[[:space:]]*/);
        s=a[2]; e=a[3];

        if(s=="" || e=="" || aro=="" || hit=="") next;

        split(hit,h,/[[:space:]]+/);
        hit_clean=h[1];
        arg=aro "|" hit_clean;

        if(s>e){tmp=s; s=e; e=tmp}

        print contig, s-1, e, arg
    }' "$rgi_txt" > "$bed_out"
}

bed_bam_ok_ratio() {
    local bed="$1"
    local bam="$2"
    local n="$3"

    local ok=0
    local tot=0
    local c

    while read -r c; do
        [[ -z "$c" ]] && continue
        tot=$((tot+1))

        if samtools idxstats "$bam" | awk -v C="$c" '$1==C{found=1} END{exit(found?0:1)}'; then
            ok=$((ok+1))
        fi

        [[ $tot -ge $n ]] && break
    done < <(cut -f1 "$bed")

    echo "${ok}/${tot}"
}

############################################
# LOOP OVER RGI FILES
############################################

shopt -s nullglob
RGI_FILES=("${RGI_DIR}"/*_rgi.txt)

if [[ ${#RGI_FILES[@]} -eq 0 ]]; then
    echo "[ERROR] No *_rgi.txt files found in $RGI_DIR"
    exit 1
fi

echo "[INFO] Found ${#RGI_FILES[@]} RGI files"

for rgi in "${RGI_FILES[@]}"; do
    sample=$(basename "$rgi" _rgi.txt)

    bam="${BAM_DIR}/${sample}.sorted.bam"
    bed="${BED_DIR}/${sample}_arg.fixed.bed"
    cov="${COV_DIR}/${sample}_arg_counts.tsv"

    if [[ ! -f "$bam" ]]; then
        echo "[WARN] ${sample}: BAM not found ($bam) -> skip"
        continue
    fi

    if [[ ! -f "${bam}.bai" ]]; then
        echo "[INFO] ${sample}: BAM index missing, creating it"
        samtools index "$bam"
    fi

    echo
    echo "[INFO] ${sample}: creating fixed BED from RGI"
    make_bed_fixed "$rgi" "$bed"

    nlines=$(wc -l < "$bed" || echo 0)
    if [[ "$nlines" -eq 0 ]]; then
        echo "[WARN] ${sample}: empty BED -> skip coverage"
        continue
    fi

    ratio=$(bed_bam_ok_ratio "$bed" "$bam" "$CHECK_N")
    echo "[INFO] ${sample}: BED contigs present in BAM (first ${CHECK_N}) = ${ratio}"

    echo "[INFO] ${sample}: coverage/counts -> $cov"
    echo -e "Sample\tContig\tStart\tEnd\tARG_ID\tLength_bp\tRead_count" > "$cov"

    split -l "$CHUNK_LINES" -d "$bed" "${bed}.chunk."

    for ch in "${bed}.chunk."*; do
        bedtools coverage -a "$ch" -b "$bam" -counts \
          | awk -v S="$sample" 'BEGIN{OFS="\t"}{len=$3-$2; print S,$1,$2+1,$3,$4,len,$(NF)}' >> "$cov"
    done

    rm -f "${bed}.chunk."*

    if awk -F'\t' 'NR>1 && $7>0{found=1; exit} END{exit(found?0:1)}' "$cov"; then
        echo "[INFO] ${sample}: OK, at least one ARG with Read_count > 0"
    else
        echo "[WARN] ${sample}: all Read_count values are 0"
    fi
done

echo
echo "[OK] Finished. Fixed BED files in $BED_DIR; counts in $COV_DIR"
