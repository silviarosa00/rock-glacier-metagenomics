# Rock-glacier metagenomics

This repository contains the bioinformatic workflow and R scripts used for the analysis of shotgun metagenomic data from rock-glacier-fed alpine substrates.

The analyses were performed to investigate microbial taxonomic composition, metal resistance genes (MRGs), antimicrobial resistance genes (ARGs), their associations with hydrogeochemical gradients, and their potential co-localisation on putative plasmid-associated contigs.

## Repository structure

```text
rock-glacier-metagenomics/
├── bioinformatic_pipeline/   # Metagenomic processing and additional bioinformatic analyses
├── R_analysis/               # Statistical analyses and figure generation
└── README.md
```

## Bioinformatic workflow

The main analytical workflow consisted of:

1. Adapter removal
2. Quality trimming and quality control
3. Metagenome assembly
4. Taxonomic profiling
5. Functional profiling
6. ARG, MRG and 16S rRNA gene annotation using an external published workflow
7. Read-based quantification and 16S-normalisation of ARG and MRG abundance
8. Identification of putative plasmid-associated contigs and ARG/MRG co-localisation
9. Downstream statistical analyses and integration with hydrogeochemical data

Scripts and commands developed and used in this study are provided in the `bioinformatic_pipeline/` directory.

### Read processing and metagenomic analysis

The initial processing of shotgun metagenomic reads included adapter removal with Cutadapt, quality trimming with Sickle, post-trimming quality assessment with FastQC, metagenome assembly with MEGAHIT, taxonomic classification with Kraken2, and functional profiling with HUMAnN.

Detailed commands and parameters are provided in the corresponding scripts in `bioinformatic_pipeline/`.

### ARG, MRG and 16S annotation

ARG, MRG and 16S rRNA gene annotations were obtained using the previously developed **ARMonaut v1.0.1** workflow, available through Zenodo:

**DOI:** https://doi.org/10.5281/zenodo.18187986

The annotation workflow itself is not duplicated in this repository. Please refer to the original Zenodo repository for the corresponding scripts, software, databases and parameters.

### ARG and MRG quantification and normalisation

Following annotation, trimmed reads were mapped back to the assembled scaffolds. Read counts overlapping annotated ARG, MRG and 16S rRNA gene regions were subsequently obtained using Bowtie2, SAMtools and BEDTools.

ARG and MRG abundance was corrected for annotated feature length and subsequently normalised to 16S rRNA gene coverage.

Scripts used for read counting and ARG/16S and MRG/16S normalisation are provided in the `bioinformatic_pipeline/` directory.

The resulting normalised abundance matrices were used as input for downstream resistome analyses in R.

### Putative plasmid-associated contigs

Putative plasmid-associated contigs were identified and compared with ARG and MRG annotations to investigate potential ARG–MRG co-localisation on mobile genetic elements.

The corresponding scripts used for the downstream matching and analysis are provided in this repository.

## Downstream analyses in R

The `R_analysis/` directory contains scripts used for statistical analyses and visualisation, including:

- taxonomic composition and diversity analyses;
- ARG and MRG composition;
- beta-diversity analyses;
- associations between environmental variables and resistome composition;
- ARG–MRG network analysis;
- integration of hydrogeochemical and resistome data;
- path modelling;
- analysis of ARG/MRG co-localisation on putative plasmid-associated contigs;
- generation of figures reported in the manuscript.

## Software and databases

Major bioinformatic tools used in the workflow include:

- Cutadapt
- Sickle
- FastQC
- MEGAHIT
- Kraken2
- HUMAnN
- Bowtie2
- SAMtools
- BEDTools

Software versions, reference databases and relevant parameters are reported in the individual scripts and in the associated manuscript.

Software and databases used for ARG, MRG and 16S annotation are documented in the external ARMonaut workflow referenced above.

## Data availability

Raw shotgun metagenomic sequencing reads will be deposited in the European Nucleotide Archive (ENA).

**ENA accession:** [TO BE ADDED]

The ARMonaut annotation workflow is available through the Zenodo repository referenced above. Scripts for additional bioinformatic processing, abundance normalisation and downstream statistical analyses performed in this study are available in this GitHub repository.

## Citation

If you use scripts associated with this study, please cite the corresponding publication:

**[PUBLICATION REFERENCE TO BE ADDED]**

The ARMonaut workflow should be cited separately according to the citation information provided in its Zenodo repository.

## Contact

Silvia Rosa  
Free University of Bozen-Bolzano
