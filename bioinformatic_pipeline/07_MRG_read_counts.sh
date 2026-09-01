#!/bin/bash
#SBATCH --job-name=bacmet_counts
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --partition=bru-medium
#SBATCH --account=brusetti-lab
#SBATCH --time=02-00:00:00
#SBATCH --output=bacmet_counts-%j.out
#SBATCH --error=bacmet_counts-%j.err

############################
# 0) Conda + modules
############################

source ~/anaconda3/etc/profile.d/conda.sh
conda activate bowtie2-env

module purge
module load samtools/1.19.2-gcc-13.3.0-a2yhwkt
module load bedtools2/2.31.1-gcc-12.1.0

echo "PATH: $PATH"
which bowtie2
which samtools
which bedtools
bowtie2 --version
samtools --version
bedtools --version

############################################
# 1) CONFIGURATION
############################################

SCAFF_DIR="/your/path/scaffolds"
READS_DIR="/your/path/sickle_out"
BACMET_DIR="/your/path/output_bacmet"
OUT_DIR="/your/path/MRG_abundance"

THREADS=8

mkdir -p "${OUT_DIR}/bam" "${OUT_DIR}/bed" "${OUT_DIR}/coverage"

echo "[INFO] Starting BacMet read-count pipeline: $(date)"

############################################
# 2) LOOP OVER SCAFFOLDS
############################################

for contig in "${SCAFF_DIR}"/*.fa; do
    sample=$(basename "${contig}" .fa)

    echo
    echo "[INFO] Sample: ${sample} ($(date '+%Y-%m-%d %H:%M:%S'))"

    ########################################
    # 2.1) Reads
    ########################################

    R1_FASTQ="${READS_DIR}/${sample}_R1.clean.fastq"
    R2_FASTQ="${READS_DIR}/${sample}_R2.clean.fastq"

    R1_GZ="${READS_DIR}/${sample}_R1.clean.fastq.gz"
    R2_GZ="${READS_DIR}/${sample}_R2.clean.fastq.gz"

    if [[ -f "$R1_FASTQ" && -f "$R2_FASTQ" ]]; then
        R1="$R1_FASTQ"
        R2="$R2_FASTQ"
        echo "[INFO] Using uncompressed reads for ${sample}"
    elif [[ -f "$R1_GZ" && -f "$R2_GZ" ]]; then
        R1="$R1_GZ"
        R2="$R2_GZ"
        echo "[INFO] Using compressed reads for ${sample}"
    else
        echo "[WARN] Reads not found for ${sample}. Skipping."
        continue
    fi

    ########################################
    # 2.2) BacMet annotation file
    ########################################

    BACMET_TXT="${BACMET_DIR}/${sample}_bacmet2_output.txt"

    if [[ ! -f "$BACMET_TXT" ]]; then
        echo "[WARN] BacMet file not found for ${sample}: ${BACMET_TXT}. Skipping."
        continue
    fi

    ########################################
    # 2.3) Bowtie2 mapping
    ########################################

    idx_prefix="${OUT_DIR}/bam/${sample}_idx"
    bam="${OUT_DIR}/bam/${sample}.sorted.bam"

    if [[ ! -f "${bam}" ]]; then
        echo "[INFO] Building Bowtie2 index for ${sample}"
        bowtie2-build "${contig}" "${idx_prefix}"

        echo "[INFO] Mapping reads to scaffold for ${sample}"
        bowtie2 -x "${idx_prefix}" \
            -1 "${R1}" -2 "${R2}" \
            -p "${THREADS}" \
            -S "${OUT_DIR}/bam/${sample}.sam"

        echo "[INFO] Converting SAM to sorted BAM for ${sample}"
        samtools view -@ "${THREADS}" -bS "${OUT_DIR}/bam/${sample}.sam" \
            | samtools sort -@ "${THREADS}" -o "${bam}"

        samtools index "${bam}"
        rm -f "${OUT_DIR}/bam/${sample}.sam"
    else
        echo "[INFO] BAM already present for ${sample}, skipping mapping."
    fi

    ########################################
    # 2.4) DIAMOND -> BED
    ########################################

    bed="${OUT_DIR}/bed/${sample}_bacmet.bed"

    if [[ ! -f "${bed}" ]]; then
        echo "[INFO] Converting DIAMOND output to BED for ${sample}"

        awk 'BEGIN{OFS="\t"}
             {
               contig=$1; gene=$2; start=$7; end=$8;
               if (start > end) {tmp=start; start=end; end=tmp}
               print contig, start-1, end, gene
             }' "${BACMET_TXT}" > "${bed}"
    else
        echo "[INFO] BacMet BED already present for ${sample}, skipping."
    fi

    ########################################
    # 2.5) bedtools coverage -counts
    ########################################

    cov_tsv="${OUT_DIR}/coverage/${sample}_bacmet_counts.tsv"

    echo "[INFO] Calculating BacMet read counts for ${sample}"

    echo -e "Sample\tContig\tStart\tEnd\tBacMet_ID\tLength_bp\tRead_count" > "${cov_tsv}"

    bedtools coverage -a "${bed}" -b "${bam}" -counts \
      | awk -v S="${sample}" 'BEGIN{OFS="\t"}
                              {
                                len = $3 - $2;
                                print S, $1, $2+1, $3, $4, len, $5
                              }' >> "${cov_tsv}"

    echo "[INFO] Finished sample ${sample}"
done

echo
echo "[OK] BacMet read-count pipeline finished: $(date)"
