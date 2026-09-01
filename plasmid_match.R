# ============================================================
# Integrate geNomad plasmid-associated CARD/RGI ARGs,
# BACMET MRGs, and metadata1 sample names
# ============================================================
# Required input files in the working directory:
# - plasmid_scaffolds_with_CARD_ARG.tsv
# - plasmid_scaffolds_with_BACMET_MRG.tsv
# - metadata1.txt
#
# Main outputs:
# - plasmid_scaffolds_with_both_CARD_ARG_and_BACMET_MRG.tsv
# - plasmid_ARG_MRG_summary_by_sample.tsv
# - plasmid_ARG_MRG_summary_by_site_season.tsv
# - plasmid_ARG_MRG_overall_numbers.tsv
# - plasmid_ARG_MRG_readme.txt
# ============================================================

# -----------------------
# CONFIG
# -----------------------

setwd("C:/your/path")
# cambia solo questa cartella

card_file <- "plasmid_scaffolds_with_CARD_ARG.tsv"
bacmet_file <- "plasmid_scaffolds_with_BACMET_MRG.tsv"
metadata_file <- "metadata1.txt"

out_dir <- "ARG_MRG_integrated_outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# HELPER FUNCTIONS
# -------------------------

read_tsv_base <- function(file) {
  if (!file.exists(file)) {
    stop(paste("File not found:", file))
  }
  
  read.delim(
    file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

present_field <- function(x) {
  x <- as.character(x)
  !is.na(x) & x != "" & x != "NA" & x != "nan" & x != "<NA>" & x != "None"
}

collapse_unique <- function(x, sep = ";") {
  x <- as.character(x)
  x <- x[present_field(x)]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  values <- unlist(strsplit(paste(x, collapse = sep), split = sep, fixed = TRUE))
  values <- trimws(values)
  values <- values[present_field(values)]
  values <- unique(values)
  
  if (length(values) == 0) {
    return(NA_character_)
  }
  
  paste(values, collapse = sep)
}

count_unique_collapsed <- function(x, sep = ";") {
  x <- collapse_unique(x, sep = sep)
  
  if (!present_field(x)) {
    return(0)
  }
  
  values <- unlist(strsplit(x, split = sep, fixed = TRUE))
  values <- trimws(values)
  values <- values[present_field(values)]
  length(unique(values))
}

clean_bracken_id <- function(sample_ids) {
  cleaned <- gsub("\\.", "-", as.character(sample_ids))
  cleaned <- sub("_N.*$", "", cleaned)
  return(cleaned)
}

prepare_metadata_for_join <- function(metadata_file) {
  
  metadata_table <- read_tsv_base(metadata_file)
  
  required_cols <- c("Sample", "bracken_name")
  missing_cols <- setdiff(required_cols, colnames(metadata_table))
  
  if (length(missing_cols) > 0) {
    stop(paste(
      "Missing required metadata columns:",
      paste(missing_cols, collapse = ", ")
    ))
  }
  
  metadata_table$Sample_name <- as.character(metadata_table$Sample)
  metadata_table$bracken_name <- as.character(metadata_table$bracken_name)
  metadata_table$Sample_clean <- clean_bracken_id(metadata_table$bracken_name)
  
  metadata_table$Site <- substr(metadata_table$Sample_name, 1, 2)
  metadata_table$Location_code <- substr(metadata_table$Site, 1, 1)
  metadata_table$Type <- substr(metadata_table$Site, 2, 2)
  metadata_table$Site_type <- metadata_table$Type
  
  metadata_table$Location <- NA_character_
  metadata_table$Location[metadata_table$Location_code == "B"] <- "Bordolona"
  metadata_table$Location[metadata_table$Location_code == "P"] <- "Preghena"
  metadata_table$Location[metadata_table$Location_code == "V"] <- "Valbiolo"
  metadata_table$Location[metadata_table$Location_code == "S"] <- "Sadole"
  metadata_table$Location[metadata_table$Location_code == "C"] <- "Cavaion"
  
  metadata_table$Valley <- NA_character_
  metadata_table$Valley[metadata_table$Location %in% c("Bordolona", "Preghena")] <- "Val Bresimo"
  metadata_table$Valley[metadata_table$Location == "Valbiolo"] <- "Valbiolo"
  metadata_table$Valley[metadata_table$Location == "Sadole"] <- "Sadole_Lagorai"
  metadata_table$Valley[metadata_table$Location == "Cavaion"] <- "Val de La Mare"
  
  if ("Time" %in% colnames(metadata_table)) {
    metadata_table$Time <- as.character(metadata_table$Time)
    metadata_table$Season <- metadata_table$Time
    metadata_table$Season[metadata_table$Time == "T1"] <- "J"
    metadata_table$Season[metadata_table$Time == "T2"] <- "S"
  } else if ("Season" %in% colnames(metadata_table)) {
    metadata_table$Season <- as.character(metadata_table$Season)
  } else {
    stop("Metadata must contain either Time or Season.")
  }
  
  if ("Rep" %in% colnames(metadata_table)) {
    metadata_table$Rep <- as.character(metadata_table$Rep)
  } else {
    metadata_table$Rep <- sub(".*?([0-9]+)$", "\\1", metadata_table$Sample_name)
  }
  
  metadata_table$Group <- paste0(metadata_table$Site, metadata_table$Season)
  
  metadata_join <- metadata_table[
    ,
    c(
      "Sample_clean",
      "Sample_name",
      "bracken_name",
      "Site",
      "Season",
      "Group",
      "Rep",
      "Type",
      "Site_type",
      "Location",
      "Valley"
    )
  ]
  
  metadata_join <- metadata_join[!duplicated(metadata_join$Sample_clean), ]
  
  return(metadata_join)
}

# -------------------------
# READ INPUT DATA
# -------------------------

card <- read_tsv_base(card_file)
bacmet <- read_tsv_base(bacmet_file)
metadata_join <- prepare_metadata_for_join(metadata_file)

# -------------------------
# CHECK REQUIRED COLUMNS
# -------------------------

required_card <- c("Sample_clean", "seq_name", "CARD_ARGs")
required_bacmet <- c("Sample_clean", "seq_name", "MRG_genes")

missing_card <- setdiff(required_card, colnames(card))
missing_bacmet <- setdiff(required_bacmet, colnames(bacmet))

if (length(missing_card) > 0) {
  stop(paste(
    "Missing required columns in CARD/RGI file:",
    paste(missing_card, collapse = ", ")
  ))
}

if (length(missing_bacmet) > 0) {
  stop(paste(
    "Missing required columns in BACMET file:",
    paste(missing_bacmet, collapse = ", ")
  ))
}

# -------------------------
# STANDARDISE KEY COLUMNS
# -------------------------

card$Sample_clean <- as.character(card$Sample_clean)
card$seq_name <- as.character(card$seq_name)

bacmet$Sample_clean <- as.character(bacmet$Sample_clean)
bacmet$seq_name <- as.character(bacmet$seq_name)

metadata_join$Sample_clean <- as.character(metadata_join$Sample_clean)

# -------------------------
# ADD METADATA TO CARD AND BACMET TABLES
# -------------------------

card <- merge(
  card,
  metadata_join,
  by = "Sample_clean",
  all.x = TRUE,
  sort = FALSE
)

bacmet <- merge(
  bacmet,
  metadata_join,
  by = "Sample_clean",
  all.x = TRUE,
  sort = FALSE
)

unmatched_card <- unique(card$Sample_clean[is.na(card$Sample_name)])
unmatched_bacmet <- unique(bacmet$Sample_clean[is.na(bacmet$Sample_name)])

if (length(unmatched_card) > 0) {
  write.table(
    data.frame(Sample_clean = unmatched_card),
    file = file.path(out_dir, "unmatched_CARD_samples_to_metadata.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  warning("Some CARD/RGI samples did not match metadata. See unmatched_CARD_samples_to_metadata.tsv")
}

if (length(unmatched_bacmet) > 0) {
  write.table(
    data.frame(Sample_clean = unmatched_bacmet),
    file = file.path(out_dir, "unmatched_BACMET_samples_to_metadata.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  warning("Some BACMET samples did not match metadata. See unmatched_BACMET_samples_to_metadata.tsv")
}

# -------------------------
# ONE ROW PER PLASMID SCAFFOLD IN CARD TABLE
# -------------------------

card$match_id <- paste(card$Sample_clean, card$seq_name, sep = "__")
split_card <- split(card, card$match_id)

card_scaffold <- do.call(
  rbind,
  lapply(split_card, function(df) {
    
    data.frame(
      Sample = if ("Sample_name" %in% colnames(df) && present_field(df$Sample_name[1])) as.character(df$Sample_name[1]) else as.character(df$Sample_clean[1]),
      Sample_clean = as.character(df$Sample_clean[1]),
      Sample_raw_CARD = if ("Sample_raw" %in% colnames(df)) as.character(df$Sample_raw[1]) else NA,
      Timepoint = if ("Timepoint" %in% colnames(df)) as.character(df$Timepoint[1]) else NA,
      
      Site = if ("Site" %in% colnames(df)) as.character(df$Site[1]) else NA,
      Season = if ("Season" %in% colnames(df)) as.character(df$Season[1]) else NA,
      Group = if ("Group" %in% colnames(df)) as.character(df$Group[1]) else NA,
      Rep = if ("Rep" %in% colnames(df)) as.character(df$Rep[1]) else NA,
      Type = if ("Type" %in% colnames(df)) as.character(df$Type[1]) else NA,
      Site_type = if ("Site_type" %in% colnames(df)) as.character(df$Site_type[1]) else NA,
      Location = if ("Location" %in% colnames(df)) as.character(df$Location[1]) else NA,
      Valley = if ("Valley" %in% colnames(df)) as.character(df$Valley[1]) else NA,
      
      seq_name = as.character(df$seq_name[1]),
      
      length = if ("length" %in% colnames(df)) safe_numeric(df$length[1]) else NA,
      plasmid_score = if ("plasmid_score" %in% colnames(df)) safe_numeric(df$plasmid_score[1]) else NA,
      marker_enrichment = if ("marker_enrichment" %in% colnames(df)) safe_numeric(df$marker_enrichment[1]) else NA,
      conjugation_genes = if ("conjugation_genes" %in% colnames(df)) collapse_unique(df$conjugation_genes) else NA,
      geNomad_amr_genes = if ("amr_genes" %in% colnames(df)) collapse_unique(df$amr_genes) else NA,
      
      n_RGI_hits = if ("n_RGI_hits" %in% colnames(df)) sum(safe_numeric(df$n_RGI_hits), na.rm = TRUE) else nrow(df),
      n_unique_CARD_ARGs = count_unique_collapsed(df$CARD_ARGs),
      CARD_ARGs = collapse_unique(df$CARD_ARGs),
      ARO_IDs = if ("ARO_IDs" %in% colnames(df)) collapse_unique(df$ARO_IDs) else NA,
      Drug_classes = if ("Drug_classes" %in% colnames(df)) collapse_unique(df$Drug_classes) else NA,
      Resistance_mechanisms = if ("Resistance_mechanisms" %in% colnames(df)) collapse_unique(df$Resistance_mechanisms) else NA,
      AMR_gene_families = if ("AMR_gene_families" %in% colnames(df)) collapse_unique(df$AMR_gene_families) else NA,
      RGI_cutoffs = if ("RGI_cutoffs" %in% colnames(df)) collapse_unique(df$RGI_cutoffs) else NA,
      
      stringsAsFactors = FALSE
    )
  })
)

# -------------------------
# ONE ROW PER PLASMID SCAFFOLD IN BACMET TABLE
# -------------------------

bacmet$match_id <- paste(bacmet$Sample_clean, bacmet$seq_name, sep = "__")
split_bacmet <- split(bacmet, bacmet$match_id)

bacmet_scaffold <- do.call(
  rbind,
  lapply(split_bacmet, function(df) {
    
    data.frame(
      Sample_clean = as.character(df$Sample_clean[1]),
      Sample_raw_BACMET = if ("Sample_raw" %in% colnames(df)) as.character(df$Sample_raw[1]) else NA,
      seq_name = as.character(df$seq_name[1]),
      
      n_BACMET_hits = if ("n_BACMET_hits" %in% colnames(df)) sum(safe_numeric(df$n_BACMET_hits), na.rm = TRUE) else nrow(df),
      n_unique_MRG_genes = count_unique_collapsed(df$MRG_genes),
      MRG_genes = collapse_unique(df$MRG_genes),
      BACMET_IDs = if ("BACMET_IDs" %in% colnames(df)) collapse_unique(df$BACMET_IDs) else NA,
      BACMET_subjects = if ("BACMET_subjects" %in% colnames(df)) collapse_unique(df$BACMET_subjects) else NA,
      
      stringsAsFactors = FALSE
    )
  })
)

# -------------------------
# INNER JOIN: SAME SAMPLE + SAME SCAFFOLD
# -------------------------

both_ARG_MRG <- merge(
  card_scaffold,
  bacmet_scaffold,
  by = c("Sample_clean", "seq_name"),
  all = FALSE,
  sort = FALSE
)

both_ARG_MRG$Co_occurrence_type <- "Same predicted plasmidic scaffold carrying CARD/RGI ARG and BACMET MRG annotations"

# -------------------------
# REORDER COLUMNS
# -------------------------

preferred_cols <- c(
  "Sample",
  "Sample_clean",
  "Sample_raw_CARD",
  "Sample_raw_BACMET",
  "Timepoint",
  "Season",
  "Site",
  "Group",
  "Rep",
  "Type",
  "Site_type",
  "Location",
  "Valley",
  "seq_name",
  "length",
  "plasmid_score",
  "marker_enrichment",
  "conjugation_genes",
  "geNomad_amr_genes",
  "CARD_ARGs",
  "ARO_IDs",
  "Drug_classes",
  "Resistance_mechanisms",
  "AMR_gene_families",
  "RGI_cutoffs",
  "MRG_genes",
  "BACMET_IDs",
  "BACMET_subjects",
  "n_RGI_hits",
  "n_unique_CARD_ARGs",
  "n_BACMET_hits",
  "n_unique_MRG_genes",
  "Co_occurrence_type"
)

preferred_cols <- preferred_cols[preferred_cols %in% colnames(both_ARG_MRG)]
other_cols <- setdiff(colnames(both_ARG_MRG), preferred_cols)

both_ARG_MRG <- both_ARG_MRG[, c(preferred_cols, other_cols)]

write.table(
  both_ARG_MRG,
  file = file.path(out_dir, "plasmid_scaffolds_with_both_CARD_ARG_and_BACMET_MRG.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -------------------------
# SAMPLE-LEVEL SUMMARY
# -------------------------

if (nrow(both_ARG_MRG) > 0) {
  
  split_sample <- split(both_ARG_MRG, both_ARG_MRG$Sample)
  
  summary_by_sample <- do.call(
    rbind,
    lapply(split_sample, function(df) {
      
      all_args <- collapse_unique(df$CARD_ARGs)
      all_mrgs <- collapse_unique(df$MRG_genes)
      all_drugs <- collapse_unique(df$Drug_classes)
      all_mech <- collapse_unique(df$Resistance_mechanisms)
      
      data.frame(
        Sample = as.character(df$Sample[1]),
        Sample_clean = as.character(df$Sample_clean[1]),
        Timepoint = if ("Timepoint" %in% colnames(df)) as.character(df$Timepoint[1]) else NA,
        Season = if ("Season" %in% colnames(df)) as.character(df$Season[1]) else NA,
        Site = if ("Site" %in% colnames(df)) as.character(df$Site[1]) else NA,
        Group = if ("Group" %in% colnames(df)) as.character(df$Group[1]) else NA,
        Rep = if ("Rep" %in% colnames(df)) as.character(df$Rep[1]) else NA,
        Type = if ("Type" %in% colnames(df)) as.character(df$Type[1]) else NA,
        Site_type = if ("Site_type" %in% colnames(df)) as.character(df$Site_type[1]) else NA,
        Location = if ("Location" %in% colnames(df)) as.character(df$Location[1]) else NA,
        Valley = if ("Valley" %in% colnames(df)) as.character(df$Valley[1]) else NA,
        
        n_plasmid_scaffolds_with_both_ARG_and_MRG = nrow(df),
        total_length_ARG_MRG_plasmid_scaffolds_bp = sum(safe_numeric(df$length), na.rm = TRUE),
        total_length_ARG_MRG_plasmid_scaffolds_kb = sum(safe_numeric(df$length), na.rm = TRUE) / 1000,
        
        n_unique_CARD_ARGs_on_ARG_MRG_plasmids = count_unique_collapsed(df$CARD_ARGs),
        n_unique_MRG_genes_on_ARG_MRG_plasmids = count_unique_collapsed(df$MRG_genes),
        
        CARD_ARGs = all_args,
        MRG_genes = all_mrgs,
        Drug_classes = all_drugs,
        Resistance_mechanisms = all_mech,
        
        stringsAsFactors = FALSE
      )
    })
  )
  
} else {
  
  summary_by_sample <- data.frame(
    Sample = character(0),
    Sample_clean = character(0),
    Timepoint = character(0),
    Season = character(0),
    Site = character(0),
    Group = character(0),
    Rep = character(0),
    Type = character(0),
    Site_type = character(0),
    Location = character(0),
    Valley = character(0),
    n_plasmid_scaffolds_with_both_ARG_and_MRG = integer(0),
    total_length_ARG_MRG_plasmid_scaffolds_bp = numeric(0),
    total_length_ARG_MRG_plasmid_scaffolds_kb = numeric(0),
    n_unique_CARD_ARGs_on_ARG_MRG_plasmids = integer(0),
    n_unique_MRG_genes_on_ARG_MRG_plasmids = integer(0),
    CARD_ARGs = character(0),
    MRG_genes = character(0),
    Drug_classes = character(0),
    Resistance_mechanisms = character(0),
    stringsAsFactors = FALSE
  )
}

write.table(
  summary_by_sample,
  file = file.path(out_dir, "plasmid_ARG_MRG_summary_by_sample.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -------------------------
# SITE-SEASON SUMMARY
# -------------------------

if (nrow(summary_by_sample) > 0) {
  
  site_season_key <- paste(summary_by_sample$Site, summary_by_sample$Season, sep = "__")
  split_site_season <- split(summary_by_sample, site_season_key)
  
  summary_by_site_season <- do.call(
    rbind,
    lapply(split_site_season, function(df) {
      
      data.frame(
        Site = as.character(df$Site[1]),
        Season = as.character(df$Season[1]),
        Group = as.character(df$Group[1]),
        Type = as.character(df$Type[1]),
        Site_type = as.character(df$Site_type[1]),
        Location = as.character(df$Location[1]),
        Valley = as.character(df$Valley[1]),
        n_samples = nrow(df),
        
        total_ARG_MRG_plasmid_scaffolds = sum(df$n_plasmid_scaffolds_with_both_ARG_and_MRG, na.rm = TRUE),
        median_ARG_MRG_plasmid_scaffolds_per_sample = median(df$n_plasmid_scaffolds_with_both_ARG_and_MRG, na.rm = TRUE),
        mean_ARG_MRG_plasmid_scaffolds_per_sample = mean(df$n_plasmid_scaffolds_with_both_ARG_and_MRG, na.rm = TRUE),
        
        total_length_ARG_MRG_plasmid_scaffolds_kb = sum(df$total_length_ARG_MRG_plasmid_scaffolds_kb, na.rm = TRUE),
        median_length_ARG_MRG_plasmid_scaffolds_kb = median(df$total_length_ARG_MRG_plasmid_scaffolds_kb, na.rm = TRUE),
        
        CARD_ARGs = collapse_unique(df$CARD_ARGs),
        MRG_genes = collapse_unique(df$MRG_genes),
        Drug_classes = collapse_unique(df$Drug_classes),
        Resistance_mechanisms = collapse_unique(df$Resistance_mechanisms),
        
        stringsAsFactors = FALSE
      )
    })
  )
  
} else {
  summary_by_site_season <- data.frame()
}

write.table(
  summary_by_site_season,
  file = file.path(out_dir, "plasmid_ARG_MRG_summary_by_site_season.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -------------------------
# OVERALL NUMBERS
# -------------------------

n_samples_CARD_ARG <- length(unique(card_scaffold$Sample))
n_samples_BACMET_MRG <- length(unique(bacmet$Sample_name[present_field(bacmet$Sample_name)]))
n_samples_both <- ifelse(nrow(summary_by_sample) > 0, nrow(summary_by_sample), 0)

unique_args_both <- count_unique_collapsed(both_ARG_MRG$CARD_ARGs)
unique_mrgs_both <- count_unique_collapsed(both_ARG_MRG$MRG_genes)

overall_numbers <- data.frame(
  Metric = c(
    "Plasmid scaffolds with CARD/RGI ARGs",
    "Plasmid scaffolds with BACMET MRGs",
    "Plasmid scaffolds with both CARD/RGI ARGs and BACMET MRGs",
    "Samples with plasmid scaffolds carrying CARD/RGI ARGs",
    "Samples with plasmid scaffolds carrying BACMET MRGs",
    "Samples with plasmid scaffolds carrying both ARGs and MRGs",
    "Unique CARD/RGI ARGs on ARG-MRG plasmid scaffolds",
    "Unique BACMET MRGs on ARG-MRG plasmid scaffolds",
    "Total length of ARG-MRG plasmid scaffolds bp",
    "Total length of ARG-MRG plasmid scaffolds kb"
  ),
  Value = c(
    nrow(card_scaffold),
    nrow(bacmet_scaffold),
    nrow(both_ARG_MRG),
    n_samples_CARD_ARG,
    n_samples_BACMET_MRG,
    n_samples_both,
    unique_args_both,
    unique_mrgs_both,
    sum(safe_numeric(both_ARG_MRG$length), na.rm = TRUE),
    sum(safe_numeric(both_ARG_MRG$length), na.rm = TRUE) / 1000
  ),
  stringsAsFactors = FALSE
)

write.table(
  overall_numbers,
  file = file.path(out_dir, "plasmid_ARG_MRG_overall_numbers.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -------------------------
# README TEXT
# -------------------------

readme_lines <- c(
  "Plasmid-associated ARG-MRG co-occurrence summary",
  "",
  "Input files:",
  paste0("- CARD/RGI plasmid scaffold table: ", card_file),
  paste0("- BACMET plasmid scaffold table: ", bacmet_file),
  paste0("- Metadata file: ", metadata_file),
  "",
  "Matching criterion:",
  "- Same Sample_clean and same seq_name.",
  "- Sample_clean is the cleaned geNomad/RGI/BACMET sample ID, e.g. T1-1.",
  "- Sample is the metadata sample name, e.g. BFJ1.",
  "- This is a conservative scaffold-level match.",
  "",
  "Output files:",
  "- plasmid_scaffolds_with_both_CARD_ARG_and_BACMET_MRG.tsv",
  "- plasmid_ARG_MRG_summary_by_sample.tsv",
  "- plasmid_ARG_MRG_summary_by_site_season.tsv",
  "- plasmid_ARG_MRG_overall_numbers.tsv",
  "",
  "Overall numbers:",
  paste(overall_numbers$Metric, overall_numbers$Value, sep = ": "),
  "",
  "Interpretation:",
  "Scaffolds in plasmid_scaffolds_with_both_CARD_ARG_and_BACMET_MRG.tsv are predicted plasmidic scaffolds carrying both CARD/RGI antibiotic-resistance annotations and BACMET metal-resistance annotations.",
  "These scaffold-level co-occurrences provide stronger evidence for direct ARG-MRG co-localization than community-level correlations.",
  "However, because these are metagenomic scaffolds rather than experimentally validated complete plasmids, the wording should remain conservative: predicted plasmidic scaffolds or plasmid-associated scaffolds."
)

writeLines(
  readme_lines,
  con = file.path(out_dir, "plasmid_ARG_MRG_readme.txt")
)

# -------------------------
# PRINT QUICK SUMMARY
# -------------------------

message("Done.")
message("Outputs written to: ", out_dir)
message("")
message("Main file:")
message(file.path(out_dir, "plasmid_scaffolds_with_both_CARD_ARG_and_BACMET_MRG.tsv"))
message("")
message("Overall numbers:")
print(overall_numbers)


# ============================================================
# ADD-ON: selected-metal MRG scaffolds on plasmids
# and CARD ARG + selected-metal MRG co-occurrence on plasmids
# ============================================================

message("Creating selected-metal plasmid MRG and CARD+metal-MRG tables...")

# -------------------------
# CONFIG FOR BACMET MAP
# -------------------------
# Put the BacMet mapping file in the same folder.
# It must contain a BacMet ID column and a Compound/metal column.
# Common column names: BacMet_ID / BACMET_ID / bacmet_id and Compound / compound
bacmet_map_file <- "bacmet.txt"

if (!file.exists(bacmet_map_file)) {
  stop(paste("BacMet mapping file not found:", bacmet_map_file))
}

if (!exists("bacmet_scaffold")) {
  stop("Object bacmet_scaffold not found. Run the previous CARD/BACMET integration script first.")
}

if (!exists("both_ARG_MRG")) {
  stop("Object both_ARG_MRG not found. Run the previous CARD/BACMET integration script first.")
}

if (!exists("out_dir")) {
  out_dir <- "ARG_MRG_integrated_outputs"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# -------------------------
# SAFETY HELPERS IF NOT ALREADY DEFINED
# -------------------------

if (!exists("present_field")) {
  present_field <- function(x) {
    x <- as.character(x)
    !is.na(x) & x != "" & x != "NA" & x != "nan" & x != "<NA>" & x != "None"
  }
}

if (!exists("collapse_unique")) {
  collapse_unique <- function(x, sep = ";") {
    x <- as.character(x)
    x <- x[present_field(x)]
    if (length(x) == 0) return(NA_character_)
    
    values <- unlist(strsplit(paste(x, collapse = sep), split = sep, fixed = TRUE))
    values <- trimws(values)
    values <- values[present_field(values)]
    values <- unique(values)
    
    if (length(values) == 0) return(NA_character_)
    paste(values, collapse = sep)
  }
}

if (!exists("count_unique_collapsed")) {
  count_unique_collapsed <- function(x, sep = ";") {
    x <- collapse_unique(x, sep = sep)
    if (!present_field(x)) return(0)
    
    values <- unlist(strsplit(x, split = sep, fixed = TRUE))
    values <- trimws(values)
    values <- values[present_field(values)]
    length(unique(values))
  }
}

if (!exists("safe_numeric")) {
  safe_numeric <- function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  }
}

find_first_col <- function(df, possible_names) {
  hit <- possible_names[possible_names %in% colnames(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

expand_semicolon_column <- function(df, id_col, new_col = "BacMet_ID_clean") {
  
  if (!id_col %in% colnames(df)) {
    stop(paste("Column not found:", id_col))
  }
  
  out_list <- lapply(seq_len(nrow(df)), function(i) {
    
    ids <- as.character(df[[id_col]][i])
    
    if (!present_field(ids)) {
      return(NULL)
    }
    
    ids <- unlist(strsplit(ids, ";", fixed = TRUE))
    ids <- trimws(ids)
    ids <- ids[present_field(ids)]
    ids <- sub("\\|.*$", "", ids)
    ids <- regmatches(ids, regexpr("BAC[0-9]+", ids))
    ids <- ids[present_field(ids)]
    
    if (length(ids) == 0) {
      return(NULL)
    }
    
    tmp <- df[rep(i, length(ids)), , drop = FALSE]
    tmp[[new_col]] <- ids
    tmp
  })
  
  out <- do.call(rbind, out_list)
  
  if (is.null(out)) {
    out <- df[0, , drop = FALSE]
    out[[new_col]] <- character(0)
  }
  
  rownames(out) <- NULL
  out
}

# -------------------------
# READ BACMET MAP
# -------------------------

bacmet_map <- read.delim(
  bacmet_map_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

bacmet_id_col <- find_first_col(
  bacmet_map,
  c("BacMet_ID", "BACMET_ID", "bacmet_id", "BacMet ID", "BACMET ID", "ID")
)

compound_col <- find_first_col(
  bacmet_map,
  c("Compound", "compound", "Compounds", "Metal", "metal", "Substance", "substrate")
)

if (is.na(bacmet_id_col)) {
  stop("Could not find BacMet ID column in bacmet_map_file.")
}

if (is.na(compound_col)) {
  stop("Could not find Compound/Metal column in bacmet_map_file.")
}

bacmet_map$BacMet_ID_clean <- as.character(bacmet_map[[bacmet_id_col]])
bacmet_map$BacMet_ID_clean <- sub("\\|.*$", "", bacmet_map$BacMet_ID_clean)
bacmet_map$BacMet_ID_clean <- regmatches(
  bacmet_map$BacMet_ID_clean,
  regexpr("BAC[0-9]+", bacmet_map$BacMet_ID_clean)
)

bacmet_map$Compound_clean <- as.character(bacmet_map[[compound_col]])
bacmet_map$compound_low <- tolower(bacmet_map$Compound_clean)

# -------------------------
# SELECTED METAL DICTIONARY
# -------------------------

metal_terms_focus <- data.frame(
  Metal = c(
    "Cu", "Zn", "Ni", "Co", "Cd", "As", "Cr", "Mn", "Fe", "Pb",
    "Hg", "Ag", "Sn", "Sb", "Se", "Mo", "W", "Au", "Al", "V", "U"
  ),
  regex = c(
    "copper",
    "zinc",
    "nickel",
    "cobalt",
    "cadmium",
    "arsenic|arsenate|arsenite",
    "chromium|chromate",
    "manganese",
    "iron",
    "lead",
    "mercury|organomercury|organo-mercury",
    "silver",
    "tin|organotin|organo-tin",
    "antimony",
    "selenium|selenate|selenite",
    "molybdenum|molybdate",
    "tungsten|wolfram",
    "gold",
    "aluminium|aluminum",
    "vanadium|vanadate",
    "uranium|uranyl"
  ),
  stringsAsFactors = FALSE
)

map_list <- lapply(seq_len(nrow(metal_terms_focus)), function(i) {
  
  metal <- metal_terms_focus$Metal[i]
  pattern <- metal_terms_focus$regex[i]
  
  hit <- grepl(pattern, bacmet_map$compound_low, ignore.case = TRUE)
  
  if (!any(hit, na.rm = TRUE)) {
    return(NULL)
  }
  
  data.frame(
    BacMet_ID_clean = bacmet_map$BacMet_ID_clean[hit],
    Compound = bacmet_map$Compound_clean[hit],
    Metal = metal,
    stringsAsFactors = FALSE
  )
})

bacmet_selected_metal_map <- do.call(rbind, map_list)

if (is.null(bacmet_selected_metal_map) || nrow(bacmet_selected_metal_map) == 0) {
  stop("No selected-metal BacMet IDs found. Check bacmet_map_file and compound column.")
}

bacmet_selected_metal_map <- bacmet_selected_metal_map[
  present_field(bacmet_selected_metal_map$BacMet_ID_clean),
]

bacmet_selected_metal_map <- unique(bacmet_selected_metal_map)

write.table(
  bacmet_selected_metal_map,
  file = file.path(out_dir, "bacmet_selected_metal_ID_map.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 1) PLASMID + SELECTED-METAL MRG
# ============================================================

bacmet_scaffold_expanded <- expand_semicolon_column(
  bacmet_scaffold,
  id_col = "BACMET_IDs",
  new_col = "BacMet_ID_clean"
)

plasmid_metal_mrg_long <- merge(
  bacmet_scaffold_expanded,
  bacmet_selected_metal_map,
  by = "BacMet_ID_clean",
  all = FALSE,
  sort = FALSE
)

if (nrow(plasmid_metal_mrg_long) > 0) {
  
  plasmid_metal_mrg_key <- paste(
    plasmid_metal_mrg_long$Sample_clean,
    plasmid_metal_mrg_long$seq_name,
    sep = "__"
  )
  
  split_plasmid_metal_mrg <- split(plasmid_metal_mrg_long, plasmid_metal_mrg_key)
  
  plasmid_scaffolds_with_selected_metal_MRG <- do.call(
    rbind,
    lapply(split_plasmid_metal_mrg, function(df) {
      
      data.frame(
        Sample_clean = as.character(df$Sample_clean[1]),
        Sample_raw_BACMET = if ("Sample_raw_BACMET" %in% colnames(df)) as.character(df$Sample_raw_BACMET[1]) else NA,
        seq_name = as.character(df$seq_name[1]),
        
        n_selected_metal_BACMET_IDs = length(unique(df$BacMet_ID_clean)),
        n_selected_metals = length(unique(df$Metal)),
        selected_metals = collapse_unique(df$Metal),
        selected_metal_BACMET_IDs = collapse_unique(df$BacMet_ID_clean),
        selected_metal_compounds = collapse_unique(df$Compound),
        
        MRG_genes = if ("MRG_genes" %in% colnames(df)) collapse_unique(df$MRG_genes) else NA,
        BACMET_IDs = if ("BACMET_IDs" %in% colnames(df)) collapse_unique(df$BACMET_IDs) else NA,
        BACMET_subjects = if ("BACMET_subjects" %in% colnames(df)) collapse_unique(df$BACMET_subjects) else NA,
        
        stringsAsFactors = FALSE
      )
    })
  )
  
} else {
  
  plasmid_scaffolds_with_selected_metal_MRG <- data.frame(
    Sample_clean = character(0),
    Sample_raw_BACMET = character(0),
    seq_name = character(0),
    n_selected_metal_BACMET_IDs = integer(0),
    n_selected_metals = integer(0),
    selected_metals = character(0),
    selected_metal_BACMET_IDs = character(0),
    selected_metal_compounds = character(0),
    MRG_genes = character(0),
    BACMET_IDs = character(0),
    BACMET_subjects = character(0),
    stringsAsFactors = FALSE
  )
}

write.table(
  plasmid_scaffolds_with_selected_metal_MRG,
  file = file.path(out_dir, "plasmid_scaffolds_with_selected_metal_MRG.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 2) PLASMID + CARD ARG + SELECTED-METAL MRG
# ============================================================

if (nrow(both_ARG_MRG) > 0) {
  
  both_ARG_MRG_expanded <- expand_semicolon_column(
    both_ARG_MRG,
    id_col = "BACMET_IDs",
    new_col = "BacMet_ID_clean"
  )
  
  plasmid_card_metal_mrg_long <- merge(
    both_ARG_MRG_expanded,
    bacmet_selected_metal_map,
    by = "BacMet_ID_clean",
    all = FALSE,
    sort = FALSE
  )
  
} else {
  
  plasmid_card_metal_mrg_long <- both_ARG_MRG[0, , drop = FALSE]
  plasmid_card_metal_mrg_long$BacMet_ID_clean <- character(0)
  plasmid_card_metal_mrg_long$Compound <- character(0)
  plasmid_card_metal_mrg_long$Metal <- character(0)
}

if (nrow(plasmid_card_metal_mrg_long) > 0) {
  
  plasmid_card_metal_mrg_key <- paste(
    plasmid_card_metal_mrg_long$Sample_clean,
    plasmid_card_metal_mrg_long$seq_name,
    sep = "__"
  )
  
  split_plasmid_card_metal_mrg <- split(
    plasmid_card_metal_mrg_long,
    plasmid_card_metal_mrg_key
  )
  
  plasmid_scaffolds_with_CARD_ARG_and_selected_metal_MRG <- do.call(
    rbind,
    lapply(split_plasmid_card_metal_mrg, function(df) {
      
      data.frame(
        Sample = if ("Sample" %in% colnames(df)) as.character(df$Sample[1]) else NA,
        Sample_clean = as.character(df$Sample_clean[1]),
        Timepoint = if ("Timepoint" %in% colnames(df)) as.character(df$Timepoint[1]) else NA,
        Season = if ("Season" %in% colnames(df)) as.character(df$Season[1]) else NA,
        Site = if ("Site" %in% colnames(df)) as.character(df$Site[1]) else NA,
        Group = if ("Group" %in% colnames(df)) as.character(df$Group[1]) else NA,
        Rep = if ("Rep" %in% colnames(df)) as.character(df$Rep[1]) else NA,
        Type = if ("Type" %in% colnames(df)) as.character(df$Type[1]) else NA,
        Site_type = if ("Site_type" %in% colnames(df)) as.character(df$Site_type[1]) else NA,
        Location = if ("Location" %in% colnames(df)) as.character(df$Location[1]) else NA,
        Valley = if ("Valley" %in% colnames(df)) as.character(df$Valley[1]) else NA,
        
        seq_name = as.character(df$seq_name[1]),
        length = if ("length" %in% colnames(df)) safe_numeric(df$length[1]) else NA,
        plasmid_score = if ("plasmid_score" %in% colnames(df)) safe_numeric(df$plasmid_score[1]) else NA,
        marker_enrichment = if ("marker_enrichment" %in% colnames(df)) safe_numeric(df$marker_enrichment[1]) else NA,
        conjugation_genes = if ("conjugation_genes" %in% colnames(df)) collapse_unique(df$conjugation_genes) else NA,
        geNomad_amr_genes = if ("geNomad_amr_genes" %in% colnames(df)) collapse_unique(df$geNomad_amr_genes) else NA,
        
        CARD_ARGs = if ("CARD_ARGs" %in% colnames(df)) collapse_unique(df$CARD_ARGs) else NA,
        ARO_IDs = if ("ARO_IDs" %in% colnames(df)) collapse_unique(df$ARO_IDs) else NA,
        Drug_classes = if ("Drug_classes" %in% colnames(df)) collapse_unique(df$Drug_classes) else NA,
        Resistance_mechanisms = if ("Resistance_mechanisms" %in% colnames(df)) collapse_unique(df$Resistance_mechanisms) else NA,
        AMR_gene_families = if ("AMR_gene_families" %in% colnames(df)) collapse_unique(df$AMR_gene_families) else NA,
        
        MRG_genes = if ("MRG_genes" %in% colnames(df)) collapse_unique(df$MRG_genes) else NA,
        BACMET_IDs = if ("BACMET_IDs" %in% colnames(df)) collapse_unique(df$BACMET_IDs) else NA,
        BACMET_subjects = if ("BACMET_subjects" %in% colnames(df)) collapse_unique(df$BACMET_subjects) else NA,
        
        n_selected_metal_BACMET_IDs = length(unique(df$BacMet_ID_clean)),
        n_selected_metals = length(unique(df$Metal)),
        selected_metals = collapse_unique(df$Metal),
        selected_metal_BACMET_IDs = collapse_unique(df$BacMet_ID_clean),
        selected_metal_compounds = collapse_unique(df$Compound),
        
        Co_occurrence_type = "Same predicted plasmidic scaffold carrying CARD/RGI ARG and selected-metal BACMET MRG annotations",
        
        stringsAsFactors = FALSE
      )
    })
  )
  
} else {
  
  plasmid_scaffolds_with_CARD_ARG_and_selected_metal_MRG <- data.frame(
    Sample = character(0),
    Sample_clean = character(0),
    seq_name = character(0),
    CARD_ARGs = character(0),
    MRG_genes = character(0),
    selected_metals = character(0),
    selected_metal_BACMET_IDs = character(0),
    stringsAsFactors = FALSE
  )
}

write.table(
  plasmid_scaffolds_with_CARD_ARG_and_selected_metal_MRG,
  file = file.path(out_dir, "plasmid_scaffolds_with_CARD_ARG_and_selected_metal_MRG.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# SAMPLE-LEVEL SUMMARIES
# ============================================================

if (nrow(plasmid_scaffolds_with_selected_metal_MRG) > 0) {
  
  split_sample_metal <- split(
    plasmid_scaffolds_with_selected_metal_MRG,
    plasmid_scaffolds_with_selected_metal_MRG$Sample_clean
  )
  
  plasmid_selected_metal_MRG_summary_by_sample <- do.call(
    rbind,
    lapply(split_sample_metal, function(df) {
      data.frame(
        Sample_clean = as.character(df$Sample_clean[1]),
        n_plasmid_scaffolds_with_selected_metal_MRG = nrow(df),
        n_selected_metals = count_unique_collapsed(df$selected_metals),
        selected_metals = collapse_unique(df$selected_metals),
        MRG_genes = collapse_unique(df$MRG_genes),
        selected_metal_BACMET_IDs = collapse_unique(df$selected_metal_BACMET_IDs),
        stringsAsFactors = FALSE
      )
    })
  )
  
} else {
  plasmid_selected_metal_MRG_summary_by_sample <- data.frame()
}

write.table(
  plasmid_selected_metal_MRG_summary_by_sample,
  file = file.path(out_dir, "plasmid_selected_metal_MRG_summary_by_sample.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (nrow(plasmid_scaffolds_with_CARD_ARG_and_selected_metal_MRG) > 0) {
  
  split_sample_card_metal <- split(
    plasmid_scaffolds_with_CARD_ARG_and_selected_metal_MRG,
    plasmid_scaffolds_with_CARD_ARG_and_selected_metal_MRG$Sample
  )
  
  plasmid_CARD_selected_metal_MRG_summary_by_sample <- do.call(
    rbind,
    lapply(split_sample_card_metal, function(df) {
      data.frame(
        Sample = as.character(df$Sample[1]),
        Sample_clean = as.character(df$Sample_clean[1]),
        Season = if ("Season" %in% colnames(df)) as.character(df$Season[1]) else NA,
        Site = if ("Site" %in% colnames(df)) as.character(df$Site[1]) else NA,
        Group = if ("Group" %in% colnames(df)) as.character(df$Group[1]) else NA,
        Location = if ("Location" %in% colnames(df)) as.character(df$Location[1]) else NA,
        Valley = if ("Valley" %in% colnames(df)) as.character(df$Valley[1]) else NA,
        
        n_plasmid_scaffolds_with_CARD_ARG_and_selected_metal_MRG = nrow(df),
        total_length_CARD_selected_metal_MRG_plasmid_scaffolds_bp = sum(safe_numeric(df$length), na.rm = TRUE),
        total_length_CARD_selected_metal_MRG_plasmid_scaffolds_kb = sum(safe_numeric(df$length), na.rm = TRUE) / 1000,
        
        n_selected_metals = count_unique_collapsed(df$selected_metals),
        selected_metals = collapse_unique(df$selected_metals),
        CARD_ARGs = collapse_unique(df$CARD_ARGs),
        MRG_genes = collapse_unique(df$MRG_genes),
        selected_metal_BACMET_IDs = collapse_unique(df$selected_metal_BACMET_IDs),
        Drug_classes = collapse_unique(df$Drug_classes),
        Resistance_mechanisms = collapse_unique(df$Resistance_mechanisms),
        
        stringsAsFactors = FALSE
      )
    })
  )
  
} else {
  plasmid_CARD_selected_metal_MRG_summary_by_sample <- data.frame()
}

write.table(
  plasmid_CARD_selected_metal_MRG_summary_by_sample,
  file = file.path(out_dir, "plasmid_CARD_selected_metal_MRG_summary_by_sample.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# README / QUICK NUMBERS
# ============================================================

metal_mrg_overall <- data.frame(
  Metric = c(
    "Plasmid scaffolds with selected-metal BACMET MRGs",
    "Samples with selected-metal plasmid MRG scaffolds",
    "Plasmid scaffolds with CARD ARGs and selected-metal BACMET MRGs",
    "Samples with CARD ARG + selected-metal MRG plasmid scaffolds",
    "Selected metals on plasmid MRG scaffolds",
    "Selected metals on CARD+metal-MRG plasmid scaffolds"
  ),
  Value = c(
    nrow(plasmid_scaffolds_with_selected_metal_MRG),
    ifelse(nrow(plasmid_selected_metal_MRG_summary_by_sample) > 0,
           nrow(plasmid_selected_metal_MRG_summary_by_sample), 0),
    nrow(plasmid_scaffolds_with_CARD_ARG_and_selected_metal_MRG),
    ifelse(nrow(plasmid_CARD_selected_metal_MRG_summary_by_sample) > 0,
           nrow(plasmid_CARD_selected_metal_MRG_summary_by_sample), 0),
    collapse_unique(plasmid_scaffolds_with_selected_metal_MRG$selected_metals),
    collapse_unique(plasmid_scaffolds_with_CARD_ARG_and_selected_metal_MRG$selected_metals)
  ),
  stringsAsFactors = FALSE
)

write.table(
  metal_mrg_overall,
  file = file.path(out_dir, "plasmid_selected_metal_MRG_overall_numbers.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

readme_metal_lines <- c(
  "Selected-metal plasmid MRG add-on summary",
  "",
  "Generated files:",
  "- plasmid_scaffolds_with_selected_metal_MRG.tsv",
  "- plasmid_scaffolds_with_CARD_ARG_and_selected_metal_MRG.tsv",
  "- plasmid_selected_metal_MRG_summary_by_sample.tsv",
  "- plasmid_CARD_selected_metal_MRG_summary_by_sample.tsv",
  "- plasmid_selected_metal_MRG_overall_numbers.tsv",
  "- bacmet_selected_metal_ID_map.tsv",
  "",
  "Definition:",
  "Selected-metal MRG scaffolds are predicted plasmidic scaffolds with BACMET IDs mapping to selected trace-metal/metalloid compounds.",
  "The CARD+selected-metal MRG table is restricted to predicted plasmidic scaffolds carrying both CARD/RGI ARG annotations and selected-metal BACMET MRG annotations.",
  "",
  "Overall numbers:",
  paste(metal_mrg_overall$Metric, metal_mrg_overall$Value, sep = ": ")
)

writeLines(
  readme_metal_lines,
  con = file.path(out_dir, "plasmid_selected_metal_MRG_readme.txt")
)

message("Selected-metal plasmid MRG add-on completed.")
message("Main outputs:")
message(file.path(out_dir, "plasmid_scaffolds_with_selected_metal_MRG.tsv"))
message(file.path(out_dir, "plasmid_scaffolds_with_CARD_ARG_and_selected_metal_MRG.tsv"))
message("")
print(metal_mrg_overall)


