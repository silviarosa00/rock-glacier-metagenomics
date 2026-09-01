#!/usr/bin/env python3

import os
import glob
import csv
from collections import defaultdict

############################################
# CONFIGURATION
############################################

ARG_DIR = "/your/path/ARG_abundance/coverage"      # *_arg_counts.tsv
S16_DIR = "/your/path/16S_coverage"                # *_16S_counts.tsv
OUT_MATRIX = "/your/path/ARG_per16S_matrix.tsv"
OUT_16S_SUMMARY = "/your/path/16S_summary.tsv"

############################################
# 1) READ 16S COVERAGE
############################################

print(">>> Reading 16S coverage from:", S16_DIR)

tot_reads_16S = defaultdict(float)
tot_len_16S = defaultdict(float)

for path in glob.glob(os.path.join(S16_DIR, "*_16S_counts.tsv")):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:
            sample = row["Sample"]
            length_bp = float(row["Length_bp"])
            read_count = float(row["Read_count"])

            tot_reads_16S[sample] += read_count
            tot_len_16S[sample] += length_bp

if not tot_reads_16S:
    raise SystemExit(
        "No *_16S_counts.tsv files found or no valid rows detected."
    )

print("Samples with 16S:")
for sample in sorted(tot_reads_16S):
    print(
        f"  {sample}: "
        f"length={tot_len_16S[sample]:.0f} bp, "
        f"reads={tot_reads_16S[sample]:.0f}"
    )

# 16S coverage per kb
cov16S = {}

for sample in tot_reads_16S:
    cov16S[sample] = (
        tot_reads_16S[sample]
        / (tot_len_16S[sample] / 1000.0)
        if tot_len_16S[sample] > 0
        else 0.0
    )

# Write 16S summary
outdir = os.path.dirname(OUT_16S_SUMMARY)
if outdir:
    os.makedirs(outdir, exist_ok=True)

with open(OUT_16S_SUMMARY, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(
        [
            "Sample",
            "Total_16S_length_bp",
            "Total_16S_reads",
            "Coverage_16S_per_kb",
        ]
    )

    for sample in sorted(cov16S):
        writer.writerow(
            [
                sample,
                f"{tot_len_16S[sample]:.0f}",
                f"{tot_reads_16S[sample]:.0f}",
                f"{cov16S[sample]:.6g}",
            ]
        )

print(f">>> Written 16S summary to: {OUT_16S_SUMMARY}")

############################################
# 2) READ ARG COUNTS AND NORMALISE TO 16S
############################################

print(">>> Reading ARG coverage from:", ARG_DIR)

matrix = defaultdict(lambda: defaultdict(float))
samples_set = set()
genes_set = set()

for path in glob.glob(os.path.join(ARG_DIR, "*_arg_counts.tsv")):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:
            sample = row["Sample"]
            gene = row["ARG_ID"]
            length_bp = float(row["Length_bp"])
            read_count = float(row["Read_count"])

            if sample not in cov16S or cov16S[sample] <= 0.0:
                continue

            # ARG coverage per kb
            cov_arg = (
                read_count / (length_bp / 1000.0)
                if length_bp > 0
                else 0.0
            )

            # ARG / 16S ratio
            value = (
                cov_arg / cov16S[sample]
                if cov16S[sample] > 0
                else 0.0
            )

            matrix[gene][sample] += value
            samples_set.add(sample)
            genes_set.add(gene)

if not matrix:
    raise SystemExit(
        "No ARG_per16S values calculated. Check paths and input headers."
    )

samples = sorted(samples_set)
genes = sorted(genes_set)

print(f"Number of ARGs (ARG_ID): {len(genes)}")
print(f"Number of samples with ARG_per16S: {len(samples)}")

############################################
# 3) WRITE GENE x SAMPLE MATRIX
############################################

outdir = os.path.dirname(OUT_MATRIX)
if outdir:
    os.makedirs(outdir, exist_ok=True)

with open(OUT_MATRIX, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["ARG_ID"] + samples)

    for gene in genes:
        writer.writerow(
            [gene] + [
                f"{matrix[gene].get(sample, 0.0):.6g}"
                for sample in samples
            ]
        )

print(f">>> Written ARG_per16S matrix to: {OUT_MATRIX}")
print("Done.")
