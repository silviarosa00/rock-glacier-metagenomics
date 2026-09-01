###############################################################################
# Bracken taxonomy pipeline: Families descriptive + Genus/Species statistics
# Input matrices must be the GitHub-style matrices produced from Bracken outputs:
#   name    taxonomy_id    taxonomy_lvl    sample1    sample2    sample3 ...
# Main analyses use Bracken relative abundances (_frac_matrix.tsv), not rarefied data.
###############################################################################
setwd("C:/your/path")
#### 0) Packages --------------------------------------------------------------
required_packages <- c(
  "readr", "dplyr", "tidyr", "stringr", "tibble",
  "ggplot2", "vegan", "ape", "pheatmap", "RColorBrewer"
)
packages_to_install <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(packages_to_install) > 0) {
  install.packages(packages_to_install, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_packages, library, character.only = TRUE))

set.seed(123)

#### 1) Input files: EDIT HERE ------------------------------------------------
# Family files are used for descriptive plots and alpha diversity.
family_frac_file  <- "bracken_combined_bacteria/bracken_F_T1_T2_Bacteria_frac_matrix.tsv"
family_num_file   <- "bracken_combined_bacteria/bracken_F_T1_T2_Bacteria_num_matrix.tsv"

# Genus and species files are used for PCoA and statistical tests.
genus_frac_file   <- "bracken_combined_bacteria/bracken_G_T1_T2_Bacteria_frac_matrix.tsv"
species_frac_file <- "bracken_combined_bacteria/bracken_S_T1_T2_Bacteria_frac_matrix.tsv"

# Metadata must contain at least:
#   Sample       = clean sample name used in plots, e.g. BFJ1, BFS2, CRJ1
#   bracken_name = sample/run name matching the Bracken matrix columns, e.g. T1-1_N2526
#   Time         = T1/T2 if Season is not already present
# Optional but useful: Rep, Valley, Site_type, Type, Location.
metadata_file <- "metadata1.txt"

#### 2) Output folders --------------------------------------------------------
main_output_dir <- "bracken_taxonomy_outputs"
descriptive_output_dir <- file.path(main_output_dir, "01_descriptive_family_level")
statistical_output_dir <- file.path(main_output_dir, "02_statistical_tests_genus_species")
figure_output_dir <- file.path(main_output_dir, "figures")
table_output_dir <- file.path(main_output_dir, "tables")

dir.create(descriptive_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(statistical_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_output_dir, showWarnings = FALSE, recursive = TRUE)

#### 3) Plot and analysis settings -------------------------------------------
minimum_mean_relative_abundance <- 0.001  # 0.1%
minimum_prevalence <- 0.10                # present in at least 10% of samples
family_top_n_barplot <- 15
family_top_n_heatmap <- 25
heatmap_pseudo_count <- 1e-6

site_levels <- c("BF", "BR", "CF", "CR", "PF", "PR", "SF", "SR", "VF", "VR")

site_cols_map <- c(
  BF = "#1b9e77", BR = "#a6dba0",
  CF = "#d95f02", CR = "#fdb863",
  PF = "#7570b3", PR = "#b2abd2",
  SF = "#1f9ac2", SR = "#a6dce7",
  VF = "#e7298a", VR = "#f2b2d4"
)

#### 4) Helper functions ------------------------------------------------------
clean_bracken_id <- function(sample_ids) {
  cleaned_ids <- sample_ids %>%
    stringr::str_replace_all("\\.", "-") %>%
    stringr::str_remove("_N.*$")
  return(cleaned_ids)
}

read_bracken_matrix <- function(input_file) {
  bracken_table <- readr::read_tsv(input_file, show_col_types = FALSE)

  required_columns <- c("name", "taxonomy_id", "taxonomy_lvl")
  missing_columns <- setdiff(required_columns, colnames(bracken_table))
  if (length(missing_columns) > 0) {
    stop("Missing required columns in ", input_file, ": ", paste(missing_columns, collapse = ", "))
  }

  taxon_metadata <- bracken_table %>%
    dplyr::select(name, taxonomy_id, taxonomy_lvl) %>%
    dplyr::mutate(
      taxon_label = ifelse(
        duplicated(name) | duplicated(name, fromLast = TRUE),
        paste0(name, "|taxid_", taxonomy_id),
        name
      )
    )

  abundance_matrix <- bracken_table %>%
    dplyr::select(-taxonomy_id, -taxonomy_lvl) %>%
    dplyr::mutate(name = taxon_metadata$taxon_label) %>%
    tibble::column_to_rownames("name") %>%
    as.data.frame()

  abundance_matrix[] <- lapply(abundance_matrix, function(column_values) {
    suppressWarnings(as.numeric(as.character(column_values)))
  })
  abundance_matrix[is.na(abundance_matrix)] <- 0

  community_table <- as.data.frame(t(abundance_matrix))

  return(list(
    community_table = community_table,
    taxon_metadata = taxon_metadata
  ))
}

prepare_metadata <- function(metadata_file) {
  
  metadata_table <- readr::read_tsv(metadata_file, show_col_types = FALSE)
  
  required_metadata_columns <- c("Sample", "bracken_name")
  missing_metadata_columns <- setdiff(required_metadata_columns, colnames(metadata_table))
  
  if (length(missing_metadata_columns) > 0) {
    stop(
      "Missing required metadata columns: ",
      paste(missing_metadata_columns, collapse = ", ")
    )
  }
  
  metadata_table <- metadata_table %>%
    dplyr::mutate(
      Sample = as.character(Sample),
      bracken_name = as.character(bracken_name),
      bracken_clean_id = clean_bracken_id(bracken_name),
      
      # Site code from sample name
      # Example: BFJ1 -> Site = BF
      Site = substr(Sample, 1, 2),
      
      # Location code and type from Site
      # Example: BF -> Location_code = B, Type = F
      Location_code = substr(Site, 1, 1),
      Type = substr(Site, 2, 2),
      Site_type = Type,
      
      # Extended location name
      Location_name = dplyr::case_when(
        Location_code == "B" ~ "Bordolona",
        Location_code == "P" ~ "Preghena",
        Location_code == "V" ~ "Valbiolo",
        Location_code == "S" ~ "Sadole",
        Location_code == "C" ~ "Cavaion",
        TRUE ~ NA_character_
      ),
      
      # Keep also a column called Location, because the rest of the script uses it
      Location = Location_name,
      
      # Valley assignment
      Valley = dplyr::case_when(
        Location_name %in% c("Bordolona", "Preghena") ~ "Val Bresimo",
        Location_name == "Valbiolo" ~ "Valbiolo",
        Location_name == "Sadole" ~ "Sadole_Lagorai",
        Location_name == "Cavaion" ~ "Val de La Mare",
        TRUE ~ NA_character_
      )
    )
  
  # Replicate column
  if ("Rep" %in% colnames(metadata_table)) {
    metadata_table <- metadata_table %>%
      dplyr::mutate(Rep = as.factor(Rep))
  } else {
    metadata_table <- metadata_table %>%
      dplyr::mutate(
        Rep = as.factor(stringr::str_extract(Sample, "[0-9]+$"))
      )
  }
  
  # Season column
  if ("Time" %in% colnames(metadata_table)) {
    metadata_table <- metadata_table %>%
      dplyr::mutate(
        Time = as.character(Time),
        Season = dplyr::case_when(
          Time == "T1" ~ "J",
          Time == "T2" ~ "S",
          TRUE ~ Time
        )
      )
  } else if ("Season" %in% colnames(metadata_table)) {
    metadata_table <- metadata_table %>%
      dplyr::mutate(
        Season = as.character(Season)
      )
  } else {
    stop("Metadata must contain either Time or Season.")
  }
  
  # Final formatting
  metadata_table <- metadata_table %>%
    dplyr::mutate(
      Season = factor(Season, levels = c("J", "S")),
      
      # Group = Site + Season
      # Example: Site BF + Season J -> BFJ
      Group = paste0(Site, Season),
      
      Site = factor(
        Site,
        levels = c("BF", "BR", "CF", "CR", "PF",
                   "PR", "SF", "SR", "VF", "VR")
      ),
      
      Location_code = factor(
        Location_code,
        levels = c("B", "P", "C", "S", "V")
      ),
      
      Location_name = factor(
        Location_name,
        levels = c("Bordolona", "Preghena", "Cavaion", "Sadole", "Valbiolo")
      ),
      
      Location = factor(
        Location,
        levels = c("Bordolona", "Preghena", "Cavaion", "Sadole", "Valbiolo")
      ),
      
      Valley = factor(
        Valley,
        levels = c("Val Bresimo", "Val de La Mare", "Sadole_Lagorai", "Valbiolo")
      ),
      
      Type = factor(Type, levels = c("F", "R")),
      Site_type = factor(Site_type, levels = c("F", "R")),
      Group = factor(Group)
    )
  
  return(metadata_table)
}
  # Harmonise optional names used in older scripts.
  if (!"Site_type" %in% colnames(metadata_table)) {
    metadata_table$Site_type <- metadata_table$Type
  }
  if (!"Valley" %in% colnames(metadata_table)) {
    metadata_table$Valley <- metadata_table$Location
  }

  return(metadata_table)
}

align_community_to_metadata <- function(community_table, metadata_table) {
  bracken_sample_ids <- rownames(community_table)
  cleaned_bracken_ids <- clean_bracken_id(bracken_sample_ids)
  metadata_match_index <- match(cleaned_bracken_ids, metadata_table$bracken_clean_id)

  if (any(is.na(metadata_match_index))) {
    unmatched_samples <- bracken_sample_ids[is.na(metadata_match_index)]
    stop(
      "These Bracken samples did not match metadata after clean_bracken_id():\n",
      paste(unmatched_samples, collapse = ", "),
      "\n\nCheck metadata$bracken_name and clean_bracken_id()."
    )
  }

  aligned_metadata <- metadata_table[metadata_match_index, , drop = FALSE]
  clean_sample_labels <- aligned_metadata$Sample

  if (any(duplicated(clean_sample_labels))) {
    warning("Duplicated Sample labels found. Making rownames unique.")
    clean_sample_labels <- make.unique(clean_sample_labels)
  }

  rownames(community_table) <- clean_sample_labels
  aligned_metadata$Sample <- clean_sample_labels

  stopifnot(all(rownames(community_table) == aligned_metadata$Sample))

  return(list(
    community_table = community_table,
    metadata_table = aligned_metadata
  ))
}

filter_taxa_for_beta <- function(relative_abundance_table) {
  taxa_prevalence <- colSums(relative_abundance_table > 0) / nrow(relative_abundance_table)
  taxa_mean <- colMeans(relative_abundance_table)

  taxa_to_keep <- taxa_mean >= minimum_mean_relative_abundance | taxa_prevalence >= minimum_prevalence
  filtered_table <- relative_abundance_table[, taxa_to_keep, drop = FALSE]

  filter_summary <- tibble::tibble(
    total_taxa_before_filtering = ncol(relative_abundance_table),
    total_taxa_after_filtering = ncol(filtered_table),
    minimum_mean_relative_abundance = minimum_mean_relative_abundance,
    minimum_prevalence = minimum_prevalence
  )

  return(list(
    filtered_table = filtered_table,
    filter_summary = filter_summary
  ))
}

save_taxa_matrix <- function(community_table, taxon_column_name, output_file) {
  readr::write_tsv(
    as.data.frame(t(community_table)) %>% tibble::rownames_to_column(taxon_column_name),
    output_file
  )
}

make_group_median_table <- function(relative_abundance_table, metadata_table) {
  metadata_aligned <- metadata_table[match(rownames(relative_abundance_table), metadata_table$Sample), , drop = FALSE]
  stopifnot(all(metadata_aligned$Sample == rownames(relative_abundance_table)))

  group_table <- as.data.frame(relative_abundance_table) %>%
    tibble::rownames_to_column("Sample") %>%
    dplyr::left_join(
      metadata_aligned %>% dplyr::select(Sample, Group, Site, Season, Location, Type),
      by = "Sample"
    ) %>%
    dplyr::group_by(Group, Site, Season, Location, Type) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), ~median(.x, na.rm = TRUE)), .groups = "drop")

  group_matrix <- group_table %>%
    dplyr::select(Group, where(is.numeric)) %>%
    tibble::column_to_rownames("Group") %>%
    as.matrix()
  storage.mode(group_matrix) <- "numeric"
  group_matrix[is.na(group_matrix)] <- 0

  group_metadata <- group_table %>%
    dplyr::select(Group, Site, Season, Location, Type) %>%
    dplyr::distinct()

  return(list(
    group_matrix = group_matrix,
    group_metadata = group_metadata
  ))
}

#### 5) Read and align all matrices ------------------------------------------
metadata_all <- prepare_metadata(metadata_file)

family_relative_raw <- read_bracken_matrix(family_frac_file)
family_counts_raw <- read_bracken_matrix(family_num_file)
genus_relative_raw <- read_bracken_matrix(genus_frac_file)
species_relative_raw <- read_bracken_matrix(species_frac_file)

family_relative_aligned <- align_community_to_metadata(family_relative_raw$community_table, metadata_all)
family_counts_aligned <- align_community_to_metadata(family_counts_raw$community_table, metadata_all)
genus_relative_aligned <- align_community_to_metadata(genus_relative_raw$community_table, metadata_all)
species_relative_aligned <- align_community_to_metadata(species_relative_raw$community_table, metadata_all)

family_relative_abundance <- family_relative_aligned$community_table
family_estimated_counts <- family_counts_aligned$community_table
genus_relative_abundance <- genus_relative_aligned$community_table
species_relative_abundance <- species_relative_aligned$community_table

sample_metadata <- family_relative_aligned$metadata_table

sample_metadata <- sample_metadata %>%
  dplyr::mutate(
    Location_code = substr(as.character(Site), 1, 1),
    Location = dplyr::case_when(
      Location_code == "B" ~ "Bordolona",
      Location_code == "P" ~ "Preghena",
      Location_code == "V" ~ "Valbiolo",
      Location_code == "S" ~ "Sadole",
      Location_code == "C" ~ "Cavaion",
      TRUE ~ NA_character_
    ),
    Location_name = Location,
    Valley = dplyr::case_when(
      Location %in% c("Bordolona", "Preghena") ~ "Val Bresimo",
      Location == "Valbiolo" ~ "Valbiolo",
      Location == "Sadole" ~ "Sadole_Lagorai",
      Location == "Cavaion" ~ "Val de La Mare",
      TRUE ~ NA_character_
    ),
    Type = substr(as.character(Site), 2, 2),
    Site_type = Type
  )

readr::write_tsv(sample_metadata, file.path(table_output_dir, "metadata_aligned_to_bracken_samples.tsv"))

#### 6) Family descriptive statistics ----------------------------------------
family_detection_summary <- tibble::tibble(
  Family = colnames(family_relative_abundance),
  Prevalence_fraction = colSums(family_relative_abundance > 0) / nrow(family_relative_abundance),
  Prevalence_n_samples = colSums(family_relative_abundance > 0),
  Mean_relative_abundance = colMeans(family_relative_abundance),
  Median_relative_abundance = apply(family_relative_abundance, 2, median),
  Max_relative_abundance = apply(family_relative_abundance, 2, max)
) %>%
  dplyr::arrange(dplyr::desc(Mean_relative_abundance))

readr::write_tsv(
  family_detection_summary,
  file.path(descriptive_output_dir, "family_detection_and_abundance_summary.tsv")
)

families_detected_50_percent <- family_detection_summary %>%
  dplyr::filter(Prevalence_fraction >= 0.50)

readr::write_tsv(
  families_detected_50_percent,
  file.path(descriptive_output_dir, "families_detected_in_at_least_50_percent_of_samples.tsv")
)

family_results_numbers <- tibble::tibble(
  n_samples = nrow(family_relative_abundance),
  n_families_total = ncol(family_relative_abundance),
  n_families_detected_in_at_least_50_percent_samples = nrow(families_detected_50_percent),
  top1_family = family_detection_summary$Family[1],
  top1_mean_relative_abundance = family_detection_summary$Mean_relative_abundance[1],
  top2_family = family_detection_summary$Family[2],
  top2_mean_relative_abundance = family_detection_summary$Mean_relative_abundance[2],
  top3_family = family_detection_summary$Family[3],
  top3_mean_relative_abundance = family_detection_summary$Mean_relative_abundance[3],
  top4_family = family_detection_summary$Family[4],
  top4_mean_relative_abundance = family_detection_summary$Mean_relative_abundance[4],
  top5_family = family_detection_summary$Family[5],
  top5_mean_relative_abundance = family_detection_summary$Mean_relative_abundance[5]
)

readr::write_tsv(
  family_results_numbers,
  file.path(descriptive_output_dir, "family_key_numbers_for_results_text.tsv")
)

#### 7) Fig. 2a-style family barplot: Top 10 + Other -------------------------

family_top_n_barplot <- 10

# Top 10 families based on mean relative abundance across all samples
family_top_taxa <- family_detection_summary %>%
  dplyr::slice_head(n = family_top_n_barplot) %>%
  dplyr::pull(Family)

# Long table: Sample x Family
family_long_for_barplot <- family_relative_abundance %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  tidyr::pivot_longer(
    cols = -Sample,
    names_to = "Family",
    values_to = "Family_relative_abundance"
  ) %>%
  dplyr::left_join(
    sample_metadata %>%
      dplyr::select(Sample, Site, Season, Location, Type),
    by = "Sample"
  ) %>%
  dplyr::mutate(
    Family_group = ifelse(Family %in% family_top_taxa, Family, "Other"),
    Site = factor(Site, levels = site_levels),
    Season = factor(Season, levels = c("J", "S"))
  )

# IMPORTANT:
# First sum all families belonging to "Other" within each sample.
# Then take the median across replicates for each Site x Season x Family_group.
family_barplot_table <- family_long_for_barplot %>%
  dplyr::group_by(Sample, Site, Season, Location, Type, Family_group) %>%
  dplyr::summarise(
    Family_relative_abundance = sum(Family_relative_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(Site, Season, Location, Type, Family_group) %>%
  dplyr::summarise(
    Family_relative_abundance = median(Family_relative_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(Site, Season) %>%
  dplyr::mutate(
    Family_relative_abundance =
      Family_relative_abundance / sum(Family_relative_abundance, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

# Other last factor level = Other on TOP of the stacked bar
# Legend breaks keep Other last in the legend too
family_barplot_table <- family_barplot_table %>%
  dplyr::mutate(
    Family_group = factor(
      Family_group,
      levels = c("Other", rev(family_top_taxa))
    )
  )

# Palette you used before: RColorBrewer Set3
# Top 10 means no repeated colours
family_colours <- RColorBrewer::brewer.pal(
  n = max(3, family_top_n_barplot),
  name = "Set3"
)[seq_along(family_top_taxa)]

family_colour_map <- stats::setNames(family_colours, family_top_taxa)
family_colour_map <- c(family_colour_map, Other = "grey70")

# Quick check: Other should exist and have non-zero abundance
print(
  family_barplot_table %>%
    dplyr::filter(Family_group == "Other") %>%
    dplyr::mutate(percent = 100 * Family_relative_abundance) %>%
    dplyr::select(Site, Season, percent)
)

family_barplot <- ggplot2::ggplot(
  family_barplot_table,
  ggplot2::aes(
    x = Site,
    y = Family_relative_abundance,
    fill = Family_group
  )
) +
  ggplot2::geom_col(
    color = "white",
    linewidth = 0.25,
    position = ggplot2::position_stack(reverse = FALSE)
  ) +
  ggplot2::facet_wrap(~ Season, nrow = 1) +
  ggplot2::scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    expand = ggplot2::expansion(mult = c(0, 0.005))
  ) +
  ggplot2::scale_fill_manual(
    values = family_colour_map,
    breaks = c(family_top_taxa, "Other"),
    drop = FALSE
  ) +
  ggplot2::labs(
    x = "Site",
    y = paste0(
      "Family composition (median across replicates)"
    ),
    fill = "Family"
  ) +
  ggplot2::theme_bw(base_size = 14) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 12
    ),
    axis.text.y = ggplot2::element_text(size = 12),
    axis.title = ggplot2::element_text(size = 15),
    strip.text = ggplot2::element_text(size = 14),
    legend.title = ggplot2::element_text(size = 15),
    legend.text = ggplot2::element_text(size = 14),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

family_barplot

# Save high-resolution PNG
ggplot2::ggsave(
  filename = file.path(
    descriptive_output_dir,
    paste0("Fig2a_family_composition_Top", family_top_n_barplot, "_bySiteSeason_HIGHRES.png")
  ),
  plot = family_barplot,
  width = 14,
  height = 6,
  dpi = 1200
)

# Save vector PDF
ggplot2::ggsave(
  filename = file.path(
    descriptive_output_dir,
    paste0("Fig2a_family_composition_Top", family_top_n_barplot, "_bySiteSeason_VECTOR.pdf")
  ),
  plot = family_barplot,
  width = 14,
  height = 6,
  device = cairo_pdf
)

# Save table used for the plot
readr::write_tsv(
  family_barplot_table,
  file.path(
    descriptive_output_dir,
    paste0("Fig2a_family_composition_Top", family_top_n_barplot, "_bySiteSeason_table.tsv")
  )
)
#### 8) Fig. 2b-style alpha diversity ---------------------------------------
# Main alpha diversity: no rarefaction.
# Observed richness is from non-zero Bracken estimated counts.
# Shannon and inverse Simpson are from Bracken relative abundances.
family_alpha_diversity <- sample_metadata %>%
  dplyr::mutate(
    Estimated_classified_reads = rowSums(family_estimated_counts),
    Observed_families = vegan::specnumber(family_estimated_counts > 0),
    Shannon_log2 = vegan::diversity(family_relative_abundance, index = "shannon") / log(2),
    Inverse_Simpson = vegan::diversity(family_relative_abundance, index = "invsimpson"),
    Pielou_evenness = Shannon_log2 / log2(pmax(Observed_families, 2)),
    Site = factor(Site, levels = site_levels),
    Season = factor(Season, levels = c("J", "S"))
  )

readr::write_tsv(
  family_alpha_diversity,
  file.path(descriptive_output_dir, "family_alpha_diversity_no_rarefaction.tsv")
)

family_alpha_summary <- family_alpha_diversity %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    classified_reads_min = min(Estimated_classified_reads, na.rm = TRUE),
    classified_reads_max = max(Estimated_classified_reads, na.rm = TRUE),
    observed_families_min = min(Observed_families, na.rm = TRUE),
    observed_families_max = max(Observed_families, na.rm = TRUE),
    shannon_log2_min = min(Shannon_log2, na.rm = TRUE),
    shannon_log2_max = max(Shannon_log2, na.rm = TRUE),
    inverse_simpson_min = min(Inverse_Simpson, na.rm = TRUE),
    inverse_simpson_max = max(Inverse_Simpson, na.rm = TRUE),
    pielou_evenness_min = min(Pielou_evenness, na.rm = TRUE),
    pielou_evenness_max = max(Pielou_evenness, na.rm = TRUE)
  )

readr::write_tsv(
  family_alpha_summary,
  file.path(descriptive_output_dir, "family_alpha_diversity_summary_for_results_text.tsv")
)

family_alpha_long <- family_alpha_diversity %>%
  dplyr::select(Sample, Site, Season, Observed_families, Shannon_log2, Inverse_Simpson) %>%
  tidyr::pivot_longer(
    cols = c(Observed_families, Shannon_log2, Inverse_Simpson),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    Metric = dplyr::recode(
      Metric,
      Observed_families = "Observed families",
      Shannon_log2 = "Shannon (log2)",
      Inverse_Simpson = "Inverse Simpson"
    ),
    Metric = factor(Metric, levels = c("Observed families", "Shannon (log2)", "Inverse Simpson"))
  )

readr::write_tsv(
  family_alpha_long,
  file.path(descriptive_output_dir, "family_alpha_diversity_long_table.tsv")
)

family_alpha_plot <- ggplot2::ggplot(
  family_alpha_long,
  ggplot2::aes(x = Site, y = Value, fill = Site)
) +
  ggplot2::geom_boxplot(outlier.shape = NA, color = "white", linewidth = 0.2) +
  ggplot2::geom_jitter(ggplot2::aes(shape = Season), width = 0.18, size = 2.2, alpha = 0.9) +
  ggplot2::facet_grid(Metric ~ Season, scales = "free_y") +
  ggplot2::scale_fill_manual(values = site_cols_map, drop = FALSE) +
  ggplot2::labs(x = NULL, y = NULL, fill = "Site", shape = "Season") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = ggplot2::element_blank(),
    legend.position = "right"
  )

ggplot2::ggsave(
  filename = file.path(figure_output_dir, "Fig2b_family_alpha_diversity_no_rarefaction_bySiteSeason.png"),
  plot = family_alpha_plot,
  width = 12,
  height = 7,
  dpi = 300
)
#### 8) Alpha diversity plot: family level, replicates separate --------------

alpha_df <- sample_metadata %>%
  dplyr::mutate(
    Estimated_classified_reads = rowSums(family_estimated_counts),
    Observed = vegan::specnumber(family_estimated_counts > 0),
    Shannon = vegan::diversity(family_relative_abundance, index = "shannon") / log(2),
    InvSimpson = vegan::diversity(family_relative_abundance, index = "invsimpson"),
    Site = factor(Site, levels = names(site_cols_map)),
    Season = factor(Season, levels = c("J", "S"))
  )

readr::write_tsv(
  alpha_df,
  file.path(descriptive_output_dir, "family_level_alpha_diversity_no_rarefaction.tsv")
)

alpha_summary <- alpha_df %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    classified_reads_min = min(Estimated_classified_reads, na.rm = TRUE),
    classified_reads_max = max(Estimated_classified_reads, na.rm = TRUE),
    observed_families_min = min(Observed, na.rm = TRUE),
    observed_families_max = max(Observed, na.rm = TRUE),
    shannon_log2_min = min(Shannon, na.rm = TRUE),
    shannon_log2_max = max(Shannon, na.rm = TRUE),
    inverse_simpson_min = min(InvSimpson, na.rm = TRUE),
    inverse_simpson_max = max(InvSimpson, na.rm = TRUE)
  )

readr::write_tsv(
  alpha_summary,
  file.path(descriptive_output_dir, "family_level_alpha_diversity_summary_for_results_text.tsv")
)

alpha_long <- alpha_df %>%
  dplyr::select(Sample, Site, Season, Observed, Shannon, InvSimpson) %>%
  tidyr::pivot_longer(
    cols = c(Observed, Shannon, InvSimpson),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    Metric = dplyr::recode(
      Metric,
      Observed = "Observed families",
      Shannon = "Shannon (log2)",
      InvSimpson = "Inverse Simpson"
    ),
    Metric = factor(
      Metric,
      levels = c("Observed families", "Shannon (log2)", "Inverse Simpson")
    ),
    Site = factor(Site, levels = names(site_cols_map)),
    Season = factor(Season, levels = c("J", "S"))
  )

readr::write_tsv(
  alpha_long,
  file.path(descriptive_output_dir, "family_level_alpha_diversity_long_table.tsv")
)

alpha_diversity_plot <- ggplot2::ggplot(
  alpha_long,
  ggplot2::aes(x = Site, y = Value, fill = Site)
) +
  ggplot2::geom_boxplot(
    outlier.shape = NA,
    color = "white",
    linewidth = 0.25,
    width = 0.65
  ) +
  ggplot2::geom_jitter(
    ggplot2::aes(shape = Season),
    width = 0.15,
    size = 2.2,
    alpha = 0.9,
    color = "black"
  ) +
  ggplot2::facet_grid(
    Metric ~ Season,
    scales = "free_y",
    switch = "y"
  ) +
  ggplot2::scale_fill_manual(values = site_cols_map, drop = FALSE) +
  ggplot2::scale_shape_manual(
    values = c("J" = 16, "S" = 17),
    drop = FALSE
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    fill = "Site",
    shape = "Season"
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 10
    ),
    axis.text.y = ggplot2::element_text(size = 10),
    strip.text.x = ggplot2::element_text(size = 12),
    strip.text.y = ggplot2::element_text(size = 11),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "right"
  )

alpha_diversity_plot

ggplot2::ggsave(
  filename = file.path(
    descriptive_output_dir,
    "Fig2b_family_alpha_diversity_boxplots_replicates_HIGHRES.png"
  ),
  plot = alpha_diversity_plot,
  width = 12,
  height = 7,
  dpi = 1200
)

ggplot2::ggsave(
  filename = file.path(
    descriptive_output_dir,
    "Fig2b_family_alpha_diversity_boxplots_replicates_VECTOR.pdf"
  ),
  plot = alpha_diversity_plot,
  width = 12,
  height = 7,
  device = cairo_pdf
)
#### Family prevalence across replicate samples ------------------------------

family_presence <- family_estimated_counts > 0

n_samples_family_prevalence <- nrow(family_presence)

family_prevalence_table <- tibble::tibble(
  Family = colnames(family_presence),
  n_samples_present = colSums(family_presence, na.rm = TRUE),
  prevalence = n_samples_present / n_samples_family_prevalence,
  prevalence_percent = 100 * prevalence
) %>%
  dplyr::arrange(dplyr::desc(prevalence), Family)

families_present_80 <- family_prevalence_table %>%
  dplyr::filter(prevalence >= 0.80)

families_present_90 <- family_prevalence_table %>%
  dplyr::filter(prevalence >= 0.90)

families_present_100 <- family_prevalence_table %>%
  dplyr::filter(prevalence == 1)

family_prevalence_summary <- tibble::tibble(
  threshold = c(">=80%", ">=90%", "100%"),
  n_samples_total = n_samples_family_prevalence,
  n_families = c(
    nrow(families_present_80),
    nrow(families_present_90),
    nrow(families_present_100)
  )
)

print(family_prevalence_summary)

cat("\nFamilies present in >=80% of samples:\n")
print(families_present_80)

cat("\nFamilies present in >=90% of samples:\n")
print(families_present_90)

cat("\nFamilies present in 100% of samples:\n")
print(families_present_100)

readr::write_tsv(
  family_prevalence_table,
  file.path(descriptive_output_dir, "family_prevalence_all_replicate_samples.tsv")
)

readr::write_tsv(
  family_prevalence_summary,
  file.path(descriptive_output_dir, "family_prevalence_threshold_summary.tsv")
)

readr::write_tsv(
  families_present_80,
  file.path(descriptive_output_dir, "families_present_in_80_percent_replicate_samples.tsv")
)

readr::write_tsv(
  families_present_90,
  file.path(descriptive_output_dir, "families_present_in_90_percent_replicate_samples.tsv")
)

readr::write_tsv(
  families_present_100,
  file.path(descriptive_output_dir, "families_present_in_100_percent_replicate_samples.tsv")
)
#### 9) Depth diagnostic instead of rarefaction ------------------------------
family_depth_control <- tibble::tibble(
  Sample = rownames(family_estimated_counts),
  Estimated_classified_reads = rowSums(family_estimated_counts),
  Detected_families_nonzero = vegan::specnumber(family_estimated_counts > 0)
) %>%
  dplyr::left_join(sample_metadata, by = "Sample") %>%
  dplyr::arrange(Estimated_classified_reads)

readr::write_tsv(
  family_depth_control,
  file.path(descriptive_output_dir, "CONTROL_family_depth_and_detected_families.tsv")
)

family_depth_plot <- ggplot2::ggplot(
  family_depth_control,
  ggplot2::aes(x = Estimated_classified_reads, y = Detected_families_nonzero, colour = Site, shape = Season)
) +
  ggplot2::geom_point(size = 3, alpha = 0.9) +
  ggplot2::scale_colour_manual(values = site_cols_map, drop = FALSE) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(
    title = "Family detection versus Bracken estimated read depth",
    x = "Estimated classified reads",
    y = "Detected families"
  )

ggplot2::ggsave(
  filename = file.path(figure_output_dir, "CONTROL_family_depth_vs_detected_families.png"),
  plot = family_depth_plot,
  width = 7,
  height = 5,
  dpi = 300
)

#### 10) Fig. 2c-style heatmap: family top 25, group medians -----------------
family_group_data <- make_group_median_table(family_relative_abundance, sample_metadata)
family_group_matrix <- family_group_data$group_matrix
family_group_metadata <- family_group_data$group_metadata

family_heatmap_top_taxa <- names(sort(colMeans(family_group_matrix), decreasing = TRUE))[1:min(family_top_n_heatmap, ncol(family_group_matrix))]
family_heatmap_relative <- family_group_matrix[, family_heatmap_top_taxa, drop = FALSE]

# Bray-Curtis clustering uses the relative abundance values.
family_heatmap_column_clustering <- hclust(vegan::vegdist(family_heatmap_relative, method = "bray"), method = "average")
family_heatmap_row_clustering <- hclust(vegan::vegdist(t(family_heatmap_relative), method = "bray"), method = "average")

# Display values are log10 relative abundance, clipped from 0.0001% to 10%.
# This makes the heatmap legend run from low percentage at the bottom to high percentage at the top.
family_heatmap_display <- log10(family_heatmap_relative + heatmap_pseudo_count)
family_heatmap_display[family_heatmap_display < -6] <- -6
family_heatmap_display[family_heatmap_display > -1] <- -1
family_heatmap_display <- t(family_heatmap_display)

family_heatmap_annotation <- family_group_metadata %>%
  dplyr::select(Group, Site) %>%
  tibble::column_to_rownames("Group")

used_sites <- intersect(names(site_cols_map), unique(as.character(family_heatmap_annotation$Site)))
family_heatmap_annotation_colours <- list(Site = site_cols_map[used_sites])

family_heatmap_colours <- grDevices::colorRampPalette(c("white", "#FDB863", "#E34A33", "#B30000"))(100)
family_heatmap_breaks <- seq(-6, -1, length.out = 101)
family_heatmap_legend_breaks <- c(-6, -5, -4, -3, -2, -1)
family_heatmap_legend_labels <- c("0.0001%", "0.001%", "0.01%", "0.1%", "1%", "10%")

readr::write_tsv(
  as.data.frame(family_heatmap_relative) %>% tibble::rownames_to_column("Group"),
  file.path(descriptive_output_dir, "family_heatmap_top25_group_median_relative_abundance.tsv")
)

pheatmap::pheatmap(
  family_heatmap_display,
  color = family_heatmap_colours,
  breaks = family_heatmap_breaks,
  cluster_rows = family_heatmap_row_clustering,
  cluster_cols = family_heatmap_column_clustering,
  annotation_col = family_heatmap_annotation,
  annotation_colors = family_heatmap_annotation_colours,
  border_color = NA,
  fontsize_row = 7,
  fontsize_col = 8,
  legend_breaks = family_heatmap_legend_breaks,
  legend_labels = family_heatmap_legend_labels,
  main = "Top 25 families (group medians; Bray-Curtis clustering; log10 relative abundance)",
  filename = file.path(figure_output_dir, "Fig2c_family_heatmap_Top25_groupMedian_BrayClustering_highPercentageTop.png"),
  width = 11,
  height = 8
)
#### Genus-level heatmap: Top 25, group-median, blue = high ------------------


# Settings
topN_genus_heatmap <- 25
heatmap_pseudo_count <- 1e-6

# Site colours, same as before
site_cols_map <- c(
  BF = "#1b9e77", BR = "#a6dba0",
  CF = "#d95f02", CR = "#fdb863",
  PF = "#7570b3", PR = "#b2abd2",
  SF = "#1f9ac2", SR = "#a6dce7",
  VF = "#e7298a", VR = "#f2b2d4"
)

# Make sure metadata are aligned to genus matrix
metadata_genus_heatmap <- sample_metadata[
  match(rownames(genus_relative_abundance), sample_metadata$Sample),
  ,
  drop = FALSE
]

stopifnot(all(metadata_genus_heatmap$Sample == rownames(genus_relative_abundance)))

metadata_genus_heatmap <- metadata_genus_heatmap %>%
  dplyr::mutate(
    Site = factor(Site, levels = names(site_cols_map)),
    Season = factor(Season, levels = c("J", "S")),
    Group = paste0(Site, Season)
  )

# Merge replicates by Group using median relative abundance
genus_group_table <- genus_relative_abundance %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::left_join(
    metadata_genus_heatmap %>%
      dplyr::select(Sample, Group, Site, Season),
    by = "Sample"
  ) %>%
  dplyr::group_by(Group, Site, Season) %>%
  dplyr::summarise(
    dplyr::across(where(is.numeric), ~ median(.x, na.rm = TRUE)),
    .groups = "drop"
  )

genus_group_metadata <- genus_group_table %>%
  dplyr::select(Group, Site, Season) %>%
  dplyr::distinct()

genus_group_matrix <- genus_group_table %>%
  dplyr::select(Group, where(is.numeric)) %>%
  tibble::column_to_rownames("Group") %>%
  as.matrix()

storage.mode(genus_group_matrix) <- "numeric"
genus_group_matrix[is.na(genus_group_matrix)] <- 0

# Select Top 25 genera by mean relative abundance
top_genus <- colMeans(genus_group_matrix, na.rm = TRUE) %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  head(topN_genus_heatmap)

genus_heatmap_relative <- genus_group_matrix[, top_genus, drop = FALSE]

# Transform to log10 relative abundance
# Higher abundance = less negative = blue
genus_heatmap_log <- log10(genus_heatmap_relative + heatmap_pseudo_count)

# Orient matrix: genera as rows, groups as columns
genus_heatmap_plot_matrix <- t(genus_heatmap_log)

# Bray-Curtis dendrogram for groups
dist_columns <- vegan::vegdist(genus_heatmap_relative, method = "bray")
hc_columns <- hclust(dist_columns, method = "average")

# Bray-Curtis dendrogram for genera
dist_rows <- vegan::vegdist(t(genus_heatmap_relative), method = "bray")
hc_rows <- hclust(dist_rows, method = "average")

# Column annotation
annotation_col <- genus_group_metadata %>%
  tibble::column_to_rownames("Group") %>%
  dplyr::select(Site)

annotation_col <- annotation_col[colnames(genus_heatmap_plot_matrix), , drop = FALSE]

annotation_colors <- list(
  Site = site_cols_map)

# Colour scale:
# low abundance = red
# intermediate = white/yellowish
# high abundance = blue
heatmap_breaks <- seq(
  min(genus_heatmap_plot_matrix, na.rm = TRUE),
  max(genus_heatmap_plot_matrix, na.rm = TRUE),
  length.out = 101
)

heatmap_colours <- colorRampPalette(
  c("#B2182B", "#FEE8C8", "#FFFFFF", "#D1E5F0", "#2166AC")
)(100)

# Legend labels for log10 relative abundance
legend_breaks <- log10(c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1))
legend_labels <- c("0.0001%", "0.001%", "0.01%", "0.1%", "1%", "10%")

# Save PNG
pheatmap::pheatmap(
  genus_heatmap_plot_matrix,
  color = heatmap_colours,
  breaks = heatmap_breaks,
  cluster_rows = hc_rows,
  cluster_cols = hc_columns,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  border_color = NA,
  fontsize_row = 12,
  fontsize_col = 8,
  treeheight_row = 45,
  treeheight_col = 45,
  legend_breaks = legend_breaks,
  legend_labels = legend_labels,
  main = "Top 25 genera - log10 relative abundance",
  filename = file.path(
    descriptive_output_dir,
    "Fig2c_genus_heatmap_Top25_groupMedian_BrayClustering_blueHigh_redLow_HIGHRES.png"
  ),
  width = 11,
  height = 8,
  dpi=1200
)

# Save PDF vector
pheatmap::pheatmap(
  genus_heatmap_plot_matrix,
  color = heatmap_colours,
  breaks = heatmap_breaks,
  cluster_rows = hc_rows,
  cluster_cols = hc_columns,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  border_color = NA,
  fontsize_row = 7,
  fontsize_col = 8,
  treeheight_row = 45,
  treeheight_col = 45,
  legend_breaks = legend_breaks,
  legend_labels = legend_labels,
  main = "Top 25 genera - log10 relative abundance",
  filename = file.path(
    descriptive_output_dir,
    "Fig_genus_heatmap_Top25_groupMedian_BrayClustering_blueHigh_redLow_VECTOR.pdf"
  ),
  width = 11,
  height = 8
)

# Save matrix used for the heatmap
readr::write_tsv(
  as.data.frame(genus_heatmap_relative) %>%
    tibble::rownames_to_column("Group"),
  file.path(
    descriptive_output_dir,
    "Fig_genus_heatmap_Top25_groupMedian_relative_abundance_matrix.tsv"
  )
)
#### 11) PCoA for Genus and Species, with and without rare taxa --------------
save_pcoa_for_taxonomic_level <- function(relative_abundance_table,
                                          metadata_table,
                                          taxonomic_level_label,
                                          analysis_label,
                                          output_subdir) {
  metadata_aligned <- metadata_table[match(rownames(relative_abundance_table), metadata_table$Sample), , drop = FALSE]
  stopifnot(all(metadata_aligned$Sample == rownames(relative_abundance_table)))

  bray_curtis_distance <- vegan::vegdist(relative_abundance_table, method = "bray")
  pcoa_result <- ape::pcoa(bray_curtis_distance)

  pcoa_scores <- as.data.frame(pcoa_result$vectors[, 1:2, drop = FALSE]) %>%
    tibble::rownames_to_column("Sample") %>%
    dplyr::rename(PCoA1 = Axis.1, PCoA2 = Axis.2) %>%
    dplyr::left_join(metadata_aligned, by = "Sample") %>%
    dplyr::mutate(
      Site = factor(Site, levels = site_levels),
      Season = factor(Season, levels = c("J", "S"))
    )

  variance_explained <- pcoa_result$values$Relative_eig[1:2] * 100

  pcoa_plot <- ggplot2::ggplot(
    pcoa_scores,
    ggplot2::aes(x = PCoA1, y = PCoA2, colour = Site, shape = Season)
  ) +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = site_cols_map, drop = FALSE) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(
      title = paste0("PCoA (Bray-Curtis, ", taxonomic_level_label, "; ", analysis_label, "; replicates not merged)"),
      x = paste0("PCoA1 (", round(variance_explained[1], 1), "%)"),
      y = paste0("PCoA2 (", round(variance_explained[2], 1), "%)"),
      colour = "Site",
      shape = "Season"
    )

  output_prefix <- paste0(taxonomic_level_label, "_", analysis_label)
  output_prefix <- stringr::str_replace_all(output_prefix, "[^A-Za-z0-9_]+", "_")

  readr::write_tsv(
    pcoa_scores,
    file.path(output_subdir, paste0(output_prefix, "_PCoA_coordinates.tsv"))
  )

  readr::write_tsv(
    tibble::tibble(
      Taxonomic_level = taxonomic_level_label,
      Analysis = analysis_label,
      PCoA1_percent = variance_explained[1],
      PCoA2_percent = variance_explained[2]
    ),
    file.path(output_subdir, paste0(output_prefix, "_PCoA_variance_explained.tsv"))
  )

  ggplot2::ggsave(
    filename = file.path(figure_output_dir, paste0(output_prefix, "_PCoA_noMergedReplicates.png")),
    plot = pcoa_plot,
    width = 8,
    height = 6,
    dpi = 1200
  )

  return(invisible(list(
    pcoa_scores = pcoa_scores,
    variance_explained = variance_explained,
    pcoa_plot = pcoa_plot
  )))
}

run_permanova_and_permdisp <- function(relative_abundance_table,
                                       metadata_table,
                                       taxonomic_level_label,
                                       analysis_label,
                                       output_subdir) {
  metadata_aligned <- metadata_table[match(rownames(relative_abundance_table), metadata_table$Sample), , drop = FALSE]
  stopifnot(all(metadata_aligned$Sample == rownames(relative_abundance_table)))

  metadata_aligned <- metadata_aligned %>%
    dplyr::mutate(
      Site = factor(Site),
      Location = factor(Location),
      Type = factor(Type),
      Season = factor(Season, levels = c("J", "S")),
      Valley = factor(Valley),
      Site_type = factor(Site_type)
    )

  output_prefix <- paste0(taxonomic_level_label, "_", analysis_label)
  output_prefix <- stringr::str_replace_all(output_prefix, "[^A-Za-z0-9_]+", "_")

  # Primary model: Location + Season + Type avoids testing Type after a Site factor that already contains Type.
  set.seed(123)
  permanova_location_type <- vegan::adonis2(
    relative_abundance_table ~ Location + Season + Type,
    data = metadata_aligned,
    method = "bray",
    permutations = 999,
    by = "margin"
  )

  permanova_location_type_table <- as.data.frame(permanova_location_type) %>%
    tibble::rownames_to_column("Term") %>%
    dplyr::mutate(Taxonomic_level = taxonomic_level_label, Analysis = analysis_label, Model = "Location + Season + Type")

  readr::write_tsv(
    permanova_location_type_table,
    file.path(output_subdir, paste0(output_prefix, "_PERMANOVA_Location_Season_Type_MARGIN.tsv"))
  )

  # Secondary model: Site + Season, useful if you want to report site-level heterogeneity.
  set.seed(123)
  permanova_site <- vegan::adonis2(
    relative_abundance_table ~ Site + Season,
    data = metadata_aligned,
    method = "bray",
    permutations = 999,
    by = "margin"
  )

  permanova_site_table <- as.data.frame(permanova_site) %>%
    tibble::rownames_to_column("Term") %>%
    dplyr::mutate(Taxonomic_level = taxonomic_level_label, Analysis = analysis_label, Model = "Site + Season")

  readr::write_tsv(
    permanova_site_table,
    file.path(output_subdir, paste0(output_prefix, "_PERMANOVA_Site_Season_MARGIN.tsv"))
  )

  bray_curtis_distance <- vegan::vegdist(relative_abundance_table, method = "bray")

  permdisp_one_factor <- function(factor_name) {
    grouping_vector <- metadata_aligned[[factor_name]]
    if (nlevels(droplevels(grouping_vector)) < 2) {
      return(tibble::tibble(Factor = factor_name, F_value = NA_real_, p_value = NA_real_))
    }
    betadisper_result <- vegan::betadisper(bray_curtis_distance, grouping_vector)
    permutation_test <- vegan::permutest(betadisper_result, permutations = 999)
    tibble::tibble(
      Factor = factor_name,
      F_value = permutation_test$tab$F[1],
      p_value = permutation_test$tab$`Pr(>F)`[1]
    )
  }

  permdisp_table <- dplyr::bind_rows(
    permdisp_one_factor("Location"),
    permdisp_one_factor("Season"),
    permdisp_one_factor("Type"),
    permdisp_one_factor("Site")
  ) %>%
    dplyr::mutate(
      p_adjusted_BH = p.adjust(p_value, method = "BH"),
      Taxonomic_level = taxonomic_level_label,
      Analysis = analysis_label
    )

  readr::write_tsv(
    permdisp_table,
    file.path(output_subdir, paste0(output_prefix, "_PERMDISP.tsv"))
  )

  return(list(
    permanova_location_type = permanova_location_type_table,
    permanova_site = permanova_site_table,
    permdisp = permdisp_table
  ))
}

# Genus and species: keep both unfiltered profiles and filtered profiles.
genus_filter <- filter_taxa_for_beta(genus_relative_abundance)
species_filter <- filter_taxa_for_beta(species_relative_abundance)

readr::write_tsv(genus_filter$filter_summary, file.path(statistical_output_dir, "genus_filter_summary.tsv"))
readr::write_tsv(species_filter$filter_summary, file.path(statistical_output_dir, "species_filter_summary.tsv"))

save_taxa_matrix(genus_relative_abundance, "Genus", file.path(table_output_dir, "genus_relative_abundance_all_taxa_taxa_x_sample.tsv"))
save_taxa_matrix(genus_filter$filtered_table, "Genus", file.path(table_output_dir, "genus_relative_abundance_filtered_taxa_x_sample.tsv"))
save_taxa_matrix(species_relative_abundance, "Species", file.path(table_output_dir, "species_relative_abundance_all_taxa_taxa_x_sample.tsv"))
save_taxa_matrix(species_filter$filtered_table, "Species", file.path(table_output_dir, "species_relative_abundance_filtered_taxa_x_sample.tsv"))

save_pcoa_for_taxonomic_level(genus_relative_abundance, sample_metadata, "Genus", "all_taxa_with_rare", statistical_output_dir)
save_pcoa_for_taxonomic_level(genus_filter$filtered_table, sample_metadata, "Genus", "filtered_rare_removed", statistical_output_dir)
save_pcoa_for_taxonomic_level(species_relative_abundance, sample_metadata, "Species", "all_taxa_with_rare", statistical_output_dir)
save_pcoa_for_taxonomic_level(species_filter$filtered_table, sample_metadata, "Species", "filtered_rare_removed", statistical_output_dir)

genus_all_stats <- run_permanova_and_permdisp(genus_relative_abundance, sample_metadata, "Genus", "all_taxa_with_rare", statistical_output_dir)
genus_filtered_stats <- run_permanova_and_permdisp(genus_filter$filtered_table, sample_metadata, "Genus", "filtered_rare_removed", statistical_output_dir)
species_all_stats <- run_permanova_and_permdisp(species_relative_abundance, sample_metadata, "Species", "all_taxa_with_rare", statistical_output_dir)
species_filtered_stats <- run_permanova_and_permdisp(species_filter$filtered_table, sample_metadata, "Species", "filtered_rare_removed", statistical_output_dir)

combined_permanova_results <- dplyr::bind_rows(
  genus_all_stats$permanova_location_type,
  genus_all_stats$permanova_site,
  genus_filtered_stats$permanova_location_type,
  genus_filtered_stats$permanova_site,
  species_all_stats$permanova_location_type,
  species_all_stats$permanova_site,
  species_filtered_stats$permanova_location_type,
  species_filtered_stats$permanova_site
)

combined_permdisp_results <- dplyr::bind_rows(
  genus_all_stats$permdisp,
  genus_filtered_stats$permdisp,
  species_all_stats$permdisp,
  species_filtered_stats$permdisp
)

readr::write_tsv(
  combined_permanova_results,
  file.path(statistical_output_dir, "combined_genus_species_PERMANOVA_results.tsv")
)

readr::write_tsv(
  combined_permdisp_results,
  file.path(statistical_output_dir, "combined_genus_species_PERMDISP_results.tsv")
)

#### 12) Copy-ready summary --------------------------------------------------
summary_lines <- c(
  "Bracken taxonomy output summary",
  "",
  paste0("Family-level descriptive outputs are in: ", descriptive_output_dir),
  paste0("Genus/species statistical outputs are in: ", statistical_output_dir),
  paste0("Figures are in: ", figure_output_dir),
  "",
  paste0("Family total taxa: ", ncol(family_relative_abundance)),
  paste0("Families detected in at least 50% of samples: ", nrow(families_detected_50_percent)),
  paste0("Top families by mean relative abundance: ", paste(head(family_detection_summary$Family, 5), collapse = ", ")),
  "",
  paste0("Alpha diversity range, observed families: ", min(family_alpha_diversity$Observed_families, na.rm = TRUE), "-", max(family_alpha_diversity$Observed_families, na.rm = TRUE)),
  paste0("Alpha diversity range, Shannon log2: ", round(min(family_alpha_diversity$Shannon_log2, na.rm = TRUE), 3), "-", round(max(family_alpha_diversity$Shannon_log2, na.rm = TRUE), 3)),
  paste0("Alpha diversity range, inverse Simpson: ", round(min(family_alpha_diversity$Inverse_Simpson, na.rm = TRUE), 3), "-", round(max(family_alpha_diversity$Inverse_Simpson, na.rm = TRUE), 3)),
  "",
  "Main beta-diversity analyses use Bracken relative abundances and no rarefaction.",
  "PCoA was saved for genus and species both with all taxa including rare taxa and after rare-taxon filtering.",
  "PERMANOVA and PERMDISP were saved for genus and species for both all-taxa and filtered profiles."
)

writeLines(summary_lines, file.path(main_output_dir, "README_results_summary.txt"))

cat("DONE. Outputs saved in: ", normalizePath(main_output_dir), "\n", sep = "")

#### Extra PERMANOVA: Valley + Season + Type ---------------------------------

# Output folder
extra_perm_dir <- file.path(statistical_output_dir, "extra_PERMANOVA_Valley_Season_Type")
dir.create(extra_perm_dir, showWarnings = FALSE, recursive = TRUE)

# Make sure metadata are aligned to community matrices
metadata_for_extra_perm <- sample_metadata %>%
  dplyr::mutate(
    Valley = factor(Valley),
    Season = factor(Season, levels = c("J", "S")),
    Type = factor(Type, levels = c("F", "R")),
    Site_type = factor(Site_type, levels = c("F", "R"))
  )

# Function for PERMANOVA + PERMDISP
run_extra_permanova <- function(relative_abundance_table,
                                metadata_table,
                                taxonomic_level,
                                output_dir) {
  
  metadata_aligned <- metadata_table[
    match(rownames(relative_abundance_table), metadata_table$Sample),
    ,
    drop = FALSE
  ]
  
  stopifnot(all(metadata_aligned$Sample == rownames(relative_abundance_table)))
  
  community_matrix <- as.matrix(relative_abundance_table)
  storage.mode(community_matrix) <- "numeric"
  community_matrix[is.na(community_matrix)] <- 0
  
  # Remove empty samples and empty taxa just in case
  community_matrix <- community_matrix[rowSums(community_matrix) > 0, , drop = FALSE]
  community_matrix <- community_matrix[, colSums(community_matrix) > 0, drop = FALSE]
  
  metadata_aligned <- metadata_aligned[
    match(rownames(community_matrix), metadata_aligned$Sample),
    ,
    drop = FALSE
  ]
  
  stopifnot(all(metadata_aligned$Sample == rownames(community_matrix)))
  
  # Drop unused levels
  metadata_aligned <- metadata_aligned %>%
    dplyr::mutate(
      Valley = droplevels(factor(Valley)),
      Season = droplevels(factor(Season)),
      Type = droplevels(factor(Type)),
      Site_type = droplevels(factor(Site_type))
    )
  
  # Bray-Curtis distance
  bray_distance <- vegan::vegdist(community_matrix, method = "bray")
  
  set.seed(123)
  
  # Main model requested: Valley + Season + Type
  permanova_valley_season_type <- vegan::adonis2(
    bray_distance ~ Valley + Season + Type,
    data = metadata_aligned,
    permutations = 999,
    by = "terms"
  )
  
  permanova_valley_season_type_margin <- vegan::adonis2(
    bray_distance ~ Valley + Season + Type,
    data = metadata_aligned,
    permutations = 999,
    by = "margin"
  )
  
  # Alternative name if you prefer Site_type instead of Type
  permanova_valley_season_sitetype <- vegan::adonis2(
    bray_distance ~ Valley + Season + Site_type,
    data = metadata_aligned,
    permutations = 999,
    by = "terms"
  )
  
  # Save PERMANOVA tables
  perm_terms_table <- as.data.frame(permanova_valley_season_type) %>%
    tibble::rownames_to_column("Term") %>%
    dplyr::mutate(
      taxonomic_level = taxonomic_level,
      model = "Valley + Season + Type",
      test = "PERMANOVA_by_terms"
    )
  
  perm_margin_table <- as.data.frame(permanova_valley_season_type_margin) %>%
    tibble::rownames_to_column("Term") %>%
    dplyr::mutate(
      taxonomic_level = taxonomic_level,
      model = "Valley + Season + Type",
      test = "PERMANOVA_by_margin"
    )
  
  perm_sitetype_table <- as.data.frame(permanova_valley_season_sitetype) %>%
    tibble::rownames_to_column("Term") %>%
    dplyr::mutate(
      taxonomic_level = taxonomic_level,
      model = "Valley + Season + Site_type",
      test = "PERMANOVA_by_terms"
    )
  
  readr::write_tsv(
    perm_terms_table,
    file.path(output_dir, paste0("PERMANOVA_", taxonomic_level, "_Valley_Season_Type_terms.tsv"))
  )
  
  readr::write_tsv(
    perm_margin_table,
    file.path(output_dir, paste0("PERMANOVA_", taxonomic_level, "_Valley_Season_Type_margin.tsv"))
  )
  
  readr::write_tsv(
    perm_sitetype_table,
    file.path(output_dir, paste0("PERMANOVA_", taxonomic_level, "_Valley_Season_SiteType_terms.tsv"))
  )
  
  # PERMDISP for each factor
  permdisp_one_factor <- function(factor_name) {
    
    factor_values <- metadata_aligned[[factor_name]]
    
    if (nlevels(droplevels(factor_values)) < 2) {
      return(NULL)
    }
    
    bd <- vegan::betadisper(bray_distance, factor_values)
    bd_perm <- vegan::permutest(bd, permutations = 999)
    
    tibble::tibble(
      taxonomic_level = taxonomic_level,
      Factor = factor_name,
      anova_F_value = as.numeric(anova(bd)$`F value`[1]),
      anova_p_value = as.numeric(anova(bd)$`Pr(>F)`[1]),
      permutation_F_value = as.numeric(bd_perm$tab$F[1]),
      permutation_p_value = as.numeric(bd_perm$tab$`Pr(>F)`[1])
    )
  }
  
  permdisp_table <- dplyr::bind_rows(
    permdisp_one_factor("Valley"),
    permdisp_one_factor("Season"),
    permdisp_one_factor("Type"),
    permdisp_one_factor("Site_type")
  ) %>%
    dplyr::mutate(
      permutation_p_adjusted_BH = p.adjust(permutation_p_value, method = "BH")
    )
  
  readr::write_tsv(
    permdisp_table,
    file.path(output_dir, paste0("PERMDISP_", taxonomic_level, "_Valley_Season_Type.tsv"))
  )
  
  return(list(
    permanova_terms = perm_terms_table,
    permanova_margin = perm_margin_table,
    permanova_sitetype = perm_sitetype_table,
    permdisp = permdisp_table
  ))
}

# Run on genus and species
extra_genus_results <- run_extra_permanova(
  relative_abundance_table = genus_relative_abundance,
  metadata_table = metadata_for_extra_perm,
  taxonomic_level = "genus",
  output_dir = extra_perm_dir
)

extra_species_results <- run_extra_permanova(
  relative_abundance_table = species_relative_abundance,
  metadata_table = metadata_for_extra_perm,
  taxonomic_level = "species",
  output_dir = extra_perm_dir
)

# Combined tables
combined_extra_permanova_terms <- dplyr::bind_rows(
  extra_genus_results$permanova_terms,
  extra_species_results$permanova_terms
) %>%
  dplyr::mutate(
    raw_p_value = `Pr(>F)`
  ) %>%
  dplyr::group_by(taxonomic_level) %>%
  dplyr::mutate(
    p_adjusted_BH = p.adjust(raw_p_value, method = "BH")
  ) %>%
  dplyr::ungroup()

combined_extra_permanova_margin <- dplyr::bind_rows(
  extra_genus_results$permanova_margin,
  extra_species_results$permanova_margin
) %>%
  dplyr::mutate(
    raw_p_value = `Pr(>F)`
  ) %>%
  dplyr::group_by(taxonomic_level) %>%
  dplyr::mutate(
    p_adjusted_BH = p.adjust(raw_p_value, method = "BH")
  ) %>%
  dplyr::ungroup()

combined_extra_permdisp <- dplyr::bind_rows(
  extra_genus_results$permdisp,
  extra_species_results$permdisp
) %>%
  dplyr::group_by(taxonomic_level) %>%
  dplyr::mutate(
    permutation_p_adjusted_BH_global = p.adjust(permutation_p_value, method = "BH")
  ) %>%
  dplyr::ungroup()

readr::write_tsv(
  combined_extra_permanova_terms,
  file.path(extra_perm_dir, "combined_PERMANOVA_genus_species_Valley_Season_Type_terms.tsv")
)

readr::write_tsv(
  combined_extra_permanova_margin,
  file.path(extra_perm_dir, "combined_PERMANOVA_genus_species_Valley_Season_Type_margin.tsv")
)

readr::write_tsv(
  combined_extra_permdisp,
  file.path(extra_perm_dir, "combined_PERMDISP_genus_species_Valley_Season_Type.tsv")
)

print(combined_extra_permanova_terms)
print(combined_extra_permdisp)

#### Extra PCoA genus level + envfit from metadata4 --------------------------
#### Genus-level PCoA + envfit metadata4, ARG/MRG style -----------------------

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(ggplot2)
library(vegan)
library(grid)

#### 0) Settings --------------------------------------------------------------

genus_envfit_dir <- file.path(statistical_output_dir, "extra_PCoA__genus_envfit_metadata4")
dir.create(genus_envfit_dir, showWarnings = FALSE, recursive = TRUE)

N_PERM <- 999
TOP_K_ARROWS <- 10
P_ADJ_CUTOFF <- 0.05

env_file <- "metadata4.csv"

shape_map <- c(
  "July" = 16,
  "September" = 17
)

#### 1) Genus community matrix ------------------------------------------------

genus_comm <- as.matrix(genus_relative_abundance)
storage.mode(genus_comm) <- "numeric"
genus_comm[is.na(genus_comm)] <- 0

# Remove empty samples/taxa, just in case
genus_comm <- genus_comm[rowSums(genus_comm) > 0, , drop = FALSE]
genus_comm <- genus_comm[, colSums(genus_comm) > 0, drop = FALSE]

meta_genus <- sample_metadata %>%
  dplyr::filter(Sample %in% rownames(genus_comm)) %>%
  dplyr::arrange(match(Sample, rownames(genus_comm))) %>%
  dplyr::mutate(
    Sample = as.character(Sample),
    Site = factor(Site, levels = names(site_cols_map)),
    Season = dplyr::case_when(
      as.character(Season) == "J" ~ "July",
      as.character(Season) == "S" ~ "September",
      TRUE ~ as.character(Season)
    ),
    Season = factor(Season, levels = c("July", "September")),
    Group = stringr::str_remove(Sample, "[0-9]+$")
  )

stopifnot(all(meta_genus$Sample == rownames(genus_comm)))

#### 2) Bray-Curtis + PCoA ----------------------------------------------------

genus_bray <- vegan::vegdist(genus_comm, method = "bray")

genus_pcoa <- cmdscale(genus_bray, k = 2, eig = TRUE)

genus_scores <- as.data.frame(genus_pcoa$points)
colnames(genus_scores) <- c("PCoA1", "PCoA2")
genus_scores$Sample <- rownames(genus_scores)

eig_pos <- genus_pcoa$eig[genus_pcoa$eig > 0]
var1 <- round(100 * genus_pcoa$eig[1] / sum(eig_pos), 1)
var2 <- round(100 * genus_pcoa$eig[2] / sum(eig_pos), 1)

plot_df <- genus_scores %>%
  dplyr::left_join(meta_genus, by = "Sample")

#### 3) Read metadata4.csv robustly ------------------------------------------

raw_lines <- readLines(env_file, encoding = "latin1", warn = FALSE)
raw_lines <- raw_lines[nchar(trimws(raw_lines)) > 0]

main_header <- strsplit(raw_lines[1], ",", fixed = TRUE)[[1]]
main_header <- gsub('"', "", main_header)
main_header <- trimws(main_header)

data_start <- grep("^C[0-9]+,", raw_lines)

if (length(data_start) == 0) {
  stop("Could not find data rows starting with C<number>, e.g. C49.")
}

first_data_line <- min(data_start)

extra_header_lines <- raw_lines[2:(first_data_line - 1)]

extra_names <- sapply(extra_header_lines, function(z) {
  zz <- strsplit(z, ",", fixed = TRUE)[[1]]
  zz <- gsub('"', "", zz)
  zz <- trimws(zz)
  zz[length(zz)]
})

extra_names <- extra_names[!is.na(extra_names) & extra_names != ""]

final_header <- c(main_header, extra_names)

data_text <- paste(raw_lines[first_data_line:length(raw_lines)], collapse = "\n")

env_raw <- readr::read_csv(
  I(data_text),
  col_names = final_header,
  show_col_types = FALSE,
  locale = readr::locale(encoding = "Latin1")
)

names(env_raw) <- names(env_raw) %>%
  iconv(from = "latin1", to = "UTF-8", sub = "") %>%
  stringr::str_replace_all("Â", "") %>%
  stringr::str_replace_all("µ", "u") %>%
  stringr::str_replace_all("μ", "u") %>%
  stringr::str_replace_all("°", "") %>%
  stringr::str_replace_all("\\s+", "_") %>%
  stringr::str_replace_all("\\(|\\)", "") %>%
  stringr::str_replace_all("/", "_") %>%
  stringr::str_replace_all("-", "_") %>%
  stringr::str_replace_all("_+", "_") %>%
  stringr::str_replace_all("^_|_$", "")

names(env_raw) <- make.unique(names(env_raw), sep = "_")

if (!("Sample" %in% names(env_raw))) {
  stop("Sample column not found in metadata4.csv after cleaning names.")
}

env_raw <- env_raw %>%
  dplyr::mutate(
    Sample = as.character(Sample),
    Group = stringr::str_remove(Sample, "[0-9]+$")
  )

#### 4) Clean environmental variable names -----------------------------------

clean_env_name <- function(x) {
  
  x <- as.character(x)
  
  x <- stringr::str_replace_all(x, "Â", "")
  x <- stringr::str_replace_all(x, "µ", "u")
  x <- stringr::str_replace_all(x, "μ", "u")
  x <- stringr::str_replace_all(x, "°", "")
  x <- stringr::str_replace_all(x, "\\.", "_")
  x <- stringr::str_replace_all(x, "\\s+", "_")
  x <- stringr::str_replace_all(x, "\\(|\\)", "")
  x <- stringr::str_replace_all(x, "/", "_")
  x <- stringr::str_replace_all(x, "-", "_")
  x <- stringr::str_replace_all(x, "_+", "_")
  x <- stringr::str_replace_all(x, "^_|_$", "")
  
  x <- dplyr::case_when(
    stringr::str_detect(x, stringr::regex("^EC|electrical", ignore_case = TRUE)) ~ "EC",
    stringr::str_detect(x, stringr::regex("^pH$", ignore_case = TRUE)) ~ "pH",
    stringr::str_detect(x, stringr::regex("Temp|Temperature|C$", ignore_case = TRUE)) ~ "Temp",
    stringr::str_detect(x, stringr::regex("Altitude|Elevation", ignore_case = TRUE)) ~ "Altitude",
    stringr::str_detect(x, stringr::regex("d2H|delta_2H|Mean_d2H", ignore_case = TRUE)) ~ "d2H",
    stringr::str_detect(x, stringr::regex("d18O|delta_18O|Mean_d18O", ignore_case = TRUE)) ~ "d18O",
    stringr::str_detect(x, stringr::regex("^F_", ignore_case = TRUE)) ~ "F",
    stringr::str_detect(x, stringr::regex("^Cl_", ignore_case = TRUE)) ~ "Cl",
    stringr::str_detect(x, stringr::regex("^SO4|Sulfate", ignore_case = TRUE)) ~ "SO4",
    stringr::str_detect(x, stringr::regex("^NO3|Nitrate", ignore_case = TRUE)) ~ "NO3",
    stringr::str_detect(x, stringr::regex("^NH4|Ammonium", ignore_case = TRUE)) ~ "NH4",
    TRUE ~ x
  )
  
  # Remove units from element/ion names
  x <- stringr::str_replace(x, "_?ug_L.*$", "")
  x <- stringr::str_replace(x, "_?mg_L.*$", "")
  x <- stringr::str_replace(x, "_?uS_cm.*$", "")
  x <- stringr::str_replace(x, "_?C$", "")
  
  x <- stringr::str_replace_all(x, "_+", "_")
  x <- stringr::str_replace_all(x, "^_|_$", "")
  
  return(x)
}

#### 5) Build chem_mat aligned to PCoA samples --------------------------------

exclude_env_columns <- c(
  "id_Lab", "idro", "Data", "Ora", "Area",
  "Codice_sorgente", "Sample", "Group", "Note"
)

env_numeric_raw <- env_raw %>%
  dplyr::select(-dplyr::any_of(exclude_env_columns)) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(),
      ~ {
        x <- as.character(.x)
        x <- gsub(",", ".", x)
        x <- gsub("<", "", x)
        x <- gsub(">", "", x)
        x <- gsub("n.d.", NA, x, fixed = TRUE)
        x <- gsub("nd", NA, x, ignore.case = TRUE)
        suppressWarnings(as.numeric(x))
      }
    )
  )

old_env_names <- colnames(env_numeric_raw)
new_env_names <- clean_env_name(old_env_names)
new_env_names <- make.unique(new_env_names, sep = "_")

rename_check <- tibble::tibble(
  old_name = old_env_names,
  clean_name = new_env_names
)

readr::write_tsv(
  rename_check,
  file.path(genus_envfit_dir, "metadata4_envfit_variable_name_cleaning.tsv")
)

colnames(env_numeric_raw) <- new_env_names

env_numeric_grouped <- dplyr::bind_cols(
  env_raw %>% dplyr::select(Group),
  env_numeric_raw
) %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::everything(),
      ~ median(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

env_genus <- plot_df %>%
  dplyr::select(Sample, Group) %>%
  dplyr::left_join(env_numeric_grouped, by = "Group") %>%
  dplyr::arrange(match(Sample, rownames(genus_comm)))

missing_env <- env_genus %>%
  dplyr::filter(if_all(-c(Sample, Group), is.na)) %>%
  dplyr::pull(Sample)

if (length(missing_env) > 0) {
  print(missing_env)
  stop("Some genus samples do not have matching environmental metadata in metadata4.csv.")
}

chem_mat <- env_genus %>%
  dplyr::select(-Sample, -Group) %>%
  as.data.frame()

# Keep only finite and variable columns
ok_cols <- sapply(chem_mat, function(x) {
  x2 <- x[is.finite(x)]
  length(x2) > 2 && stats::sd(x2, na.rm = TRUE) > 0
})

chem_mat <- chem_mat[, ok_cols, drop = FALSE]

if (ncol(chem_mat) == 0) {
  stop("No usable numeric environmental variables found for envfit.")
}

# Replace remaining NA/Inf with column median
chem_mat <- chem_mat %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(),
      ~ {
        x <- .x
        x[!is.finite(x)] <- NA
        if (any(is.na(x))) {
          x[is.na(x)] <- stats::median(x, na.rm = TRUE)
        }
        x
      }
    )
  ) %>%
  as.data.frame()

rownames(chem_mat) <- env_genus$Sample

stopifnot(all(rownames(chem_mat) == rownames(genus_comm)))

cat("Environmental variables used for envfit:\n")
print(colnames(chem_mat))

#### Keep only environmental variables used in ARG/MRG envfit -----------------

KEEP_ENV_VARS_ARG_MRG <- c(
  "T", "Cd", "Ce", "La", "alt", "Co", "U", "ph", "Ba", "Al",
  "Cr", "Cu", "Rb", "Fe", "Zn", "Mn", "meta_number", "Sr", "Pb",
  "Se", "As", "Mo", "Ti", "V", "Ni", "EC", "Sb", "Li", "Be", "E"
)

# Match genus-cleaned names to ARG/MRG-style names
colnames(chem_mat) <- dplyr::recode(
  colnames(chem_mat),
  "Temp" = "T",
  "Altitude" = "alt",
  "pH" = "ph",
  .default = colnames(chem_mat)
)

vars_present <- intersect(KEEP_ENV_VARS_ARG_MRG, colnames(chem_mat))
vars_missing <- setdiff(KEEP_ENV_VARS_ARG_MRG, colnames(chem_mat))

cat("Variables kept for genus envfit, matching ARG/MRG:\n")
print(vars_present)

cat("Variables requested but not found in genus chem_mat:\n")
print(vars_missing)

chem_mat <- chem_mat[, vars_present, drop = FALSE]

if (ncol(chem_mat) == 0) {
  stop("No ARG/MRG-matching environmental variables found in chem_mat.")
}

#### 6) envfit ---------------------------------------------------------------

set.seed(123)

fit <- vegan::envfit(
  genus_pcoa$points,
  chem_mat,
  permutations = N_PERM
)

vec <- as.data.frame(vegan::scores(fit, display = "vectors"))
vec$Variable <- rownames(vec)

r2 <- fit$vectors$r
p_value <- fit$vectors$pvals

env_tbl <- tibble::tibble(
  Variable = names(r2),
  r2 = as.numeric(r2),
  p_value = as.numeric(p_value)
) %>%
  dplyr::mutate(
    p_adj = p.adjust(p_value, method = "BH")
  ) %>%
  dplyr::arrange(p_adj, dplyr::desc(r2))

readr::write_tsv(
  env_tbl,
  file.path(genus_envfit_dir, "envfit_Genus_PCoA_metadata4_BH.tsv")
)

# Draw significant variables after BH.
# If none are significant, draw top variables by r2, like ARG/MRG scripts.
arrow_keep <- env_tbl %>%
  dplyr::filter(is.finite(p_adj), p_adj < P_ADJ_CUTOFF) %>%
  dplyr::slice_max(order_by = r2, n = TOP_K_ARROWS, with_ties = FALSE)

if (nrow(arrow_keep) == 0) {
  message("No envfit variables significant after BH correction. Drawing top variables by r2.")
  arrow_keep <- env_tbl %>%
    dplyr::slice_max(order_by = r2, n = TOP_K_ARROWS, with_ties = FALSE)
}

vec2 <- vec %>%
  dplyr::select(Variable, 1, 2)

names(vec2)[2:3] <- c("x", "y")

arrow_df <- vec2 %>%
  dplyr::inner_join(arrow_keep, by = "Variable")

#### 7) Scale arrows and place labels, ARG/MRG style -------------------------

xrange <- diff(range(plot_df$PCoA1, na.rm = TRUE))
yrange <- diff(range(plot_df$PCoA2, na.rm = TRUE))
mult <- 0.6 * min(xrange, yrange)

arrow_df <- arrow_df %>%
  dplyr::mutate(
    strength = sqrt(r2),
    xend = x * mult * strength,
    yend = y * mult * strength,
    angle = atan2(yend, xend)
  ) %>%
  dplyr::arrange(angle) %>%
  dplyr::mutate(
    off_id = rep(c(-1, 1, -2, 2, -3, 3, -4, 4, -5, 5), length.out = dplyr::n()),
    off = 0.04 * min(xrange, yrange) * off_id,
    perp_x = -sin(angle),
    perp_y = cos(angle),
    label_x = xend * 1.08 + off * perp_x,
    label_y = yend * 1.08 + off * perp_y
  )

# Optional manual label correction, same idea as ARG/MRG scripts.
# Edit dx/dy if labels overlap.
label_adjust <- tibble::tribble(
  ~Variable, ~dx,   ~dy,
  "Temp",     0.00,  0.015,
  "Ba",       0.00,  0.000,
  "Al",       0.00, -0.015,
  "Si",       0.015, 0.000
)

arrow_df <- arrow_df %>%
  dplyr::left_join(label_adjust, by = "Variable") %>%
  dplyr::mutate(
    dx = dplyr::coalesce(dx, 0),
    dy = dplyr::coalesce(dy, 0),
    label_x = label_x + dx,
    label_y = label_y + dy
  )

readr::write_tsv(
  arrow_df,
  file.path(genus_envfit_dir, "envfit_Genus_PCoA_arrows_drawn.tsv")
)

#### 8) Plot -----------------------------------------------------------------

plot_df <- plot_df %>%
  dplyr::mutate(
    Site = factor(Site, levels = names(site_cols_map)),
    Season = factor(Season, levels = c("July", "September"))
  )

missing_sites <- setdiff(unique(as.character(plot_df$Site)), names(site_cols_map))

if (length(missing_sites) > 0) {
  warning(
    "These Site levels are not in site_cols_map: ",
    paste(missing_sites, collapse = ", ")
  )
}

genus_pcoa_envfit_plot <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(x = PCoA1, y = PCoA2, color = Site, shape = Season)
) +
  ggplot2::geom_point(size = 3, alpha = 0.9) +
  ggplot2::scale_color_manual(values = site_cols_map, drop = FALSE) +
  ggplot2::scale_shape_manual(values = shape_map, drop = FALSE) +
  ggplot2::geom_segment(
    data = arrow_df,
    ggplot2::aes(x = 0, y = 0, xend = xend, yend = yend),
    inherit.aes = FALSE,
    linewidth = 0.6,
    arrow = grid::arrow(length = grid::unit(0.22, "cm"))
  ) +
  ggplot2::geom_text(
    data = arrow_df,
    ggplot2::aes(x = label_x, y = label_y, label = Variable),
    inherit.aes = FALSE,
    size = 3
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(10, 35, 10, 10)
  ) +
  ggplot2::labs(
    x = paste0("PCoA1 (", var1, "%)"),
    y = paste0("PCoA2 (", var2, "%)"),
    color = "Site",
    shape = "Season",
    title = "PCoA of genus-level community composition"
  )

genus_pcoa_envfit_plot

ggplot2::ggsave(
  filename = file.path(
    genus_envfit_dir,
    "PCoA_Genus_Bray_metadata4_envfit_ARG_MRG_style.png"
  ),
  plot = genus_pcoa_envfit_plot,
  width = 8,
  height = 6,
  dpi = 1200
)

ggplot2::ggsave(
  filename = file.path(
    genus_envfit_dir,
    "PCoA_Genus_Bray_metadata4_envfit_ARG_MRG_style.pdf"
  ),
  plot = genus_pcoa_envfit_plot,
  width = 8,
  height = 6,
  device = cairo_pdf
)

readr::write_tsv(
  plot_df,
  file.path(genus_envfit_dir, "PCoA_Genus_scores_metadata.tsv")
)

cat("Saved genus PCoA envfit outputs in:\n", normalizePath(genus_envfit_dir), "\n")

#### Check Ba correlations with other environmental variables -----------------

ba_cor_table <- chem_mat %>%
  as.data.frame() %>%
  dplyr::select(dplyr::any_of(c(
    "Ba", "T", "Temp", "Si", "Al", "Cd", "Pb", "Sb", "U", "Co",
    "SO4", "Cl", "pH", "ph", "EC", "alt"
  ))) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))

ba_cor <- cor(
  ba_cor_table,
  use = "pairwise.complete.obs",
  method = "spearman"
)

ba_cor_with_ba <- tibble::tibble(
  Variable = rownames(ba_cor),
  Spearman_rho_with_Ba = ba_cor[, "Ba"]
) %>%
  dplyr::filter(Variable != "Ba") %>%
  dplyr::arrange(dplyr::desc(abs(Spearman_rho_with_Ba)))

print(ba_cor_with_ba)
#### Ba distribution by Site --------------------------------------------------

plot_df_ba <- plot_df %>%
  dplyr::select(Sample, Site, Season, Group) %>%
  dplyr::left_join(
    env_genus %>% dplyr::select(Sample, Ba),
    by = "Sample"
  ) %>%
  dplyr::mutate(
    Site = factor(Site, levels = names(site_cols_map)),
    Season = factor(Season, levels = c("July", "September", "J", "S"))
  )

ba_site_plot <- ggplot2::ggplot(
  plot_df_ba,
  ggplot2::aes(x = Site, y = Ba, fill = Site)
) +
  ggplot2::geom_boxplot(
    outlier.shape = NA,
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  ggplot2::geom_jitter(
    ggplot2::aes(shape = Season),
    width = 0.15,
    size = 2.4,
    alpha = 0.9,
    color = "black"
  ) +
  ggplot2::scale_fill_manual(values = site_cols_map, drop = FALSE) +
  ggplot2::scale_shape_manual(
    values = c(
      "July" = 16,
      "September" = 17,
      "J" = 16,
      "S" = 17
    ),
    drop = FALSE
  ) +
  ggplot2::labs(
    x = "Site",
    y = "Ba",
    fill = "Site",
    shape = "Season"
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "right"
  )

ba_site_plot

readr::write_tsv(
  plot_df_ba,
  file.path(genus_envfit_dir, "Ba_values_by_site_for_genus_envfit.tsv")
)

ggplot2::ggsave(
  filename = file.path(
    genus_envfit_dir,
    "Ba_distribution_by_site_genus_envfit.png"
  ),
  plot = ba_site_plot,
  width = 8,
  height = 5,
  dpi = 600
)

ggplot2::ggsave(
  filename = file.path(
    genus_envfit_dir,
    "Ba_distribution_by_site_genus_envfit.pdf"
  ),
  plot = ba_site_plot,
  width = 8,
  height = 5,
  device = cairo_pdf
)
                 
                 
#### Genus-level Venn diagrams: locations and site-season groups --------------

library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggvenn)
library(patchwork)

venn_output_dir <- file.path(descriptive_output_dir, "genus_level_venn_diagrams")
dir.create(venn_output_dir, showWarnings = FALSE, recursive = TRUE)

ABUNDANCE_THRESHOLD <- 0.01
PREVALENCE_THRESHOLD <- 0.80

#### 1) Prepare genus matrix and metadata ------------------------------------

genus_mat <- as.matrix(genus_relative_abundance)
storage.mode(genus_mat) <- "numeric"
genus_mat[is.na(genus_mat)] <- 0

meta_venn <- sample_metadata %>%
  dplyr::filter(Sample %in% rownames(genus_mat)) %>%
  dplyr::arrange(match(Sample, rownames(genus_mat))) %>%
  dplyr::mutate(
    Sample = as.character(Sample),
    Site = as.character(Site),
    Season = as.character(Season),
    Location_code = substr(Site, 1, 1),
    Type = substr(Site, 2, 2),
    Group = paste0(Site, Season)
  )

stopifnot(all(meta_venn$Sample == rownames(genus_mat)))

#### 2) Helper: taxa present in >=80% samples with abundance >=1% -------------

get_core_taxa <- function(samples, abundance_threshold = 0.01, prevalence_threshold = 0.80) {
  
  samples <- intersect(samples, rownames(genus_mat))
  
  if (length(samples) == 0) {
    return(character(0))
  }
  
  submat <- genus_mat[samples, , drop = FALSE]
  
  presence <- submat >= abundance_threshold
  
  prevalence <- colSums(presence, na.rm = TRUE) / nrow(presence)
  
  taxa <- names(prevalence)[prevalence >= prevalence_threshold]
  
  return(taxa)
}

#### 3) Location-level sets: central Venn ------------------------------------

location_levels <- c("B", "P", "S", "V", "C")

location_labels <- c(
  B = "Bordolona",
  P = "Preghena",
  S = "Sadole",
  V = "Valbiolo",
  C = "Cavaion"
)

location_sets <- lapply(location_levels, function(loc) {
  samples_loc <- meta_venn %>%
    dplyr::filter(Location_code == loc) %>%
    dplyr::pull(Sample)
  
  get_core_taxa(
    samples = samples_loc,
    abundance_threshold = ABUNDANCE_THRESHOLD,
    prevalence_threshold = PREVALENCE_THRESHOLD
  )
})

names(location_sets) <- names(location_labels[location_levels])

#### 4) Group-level sets for each location -----------------------------------

make_location_group_sets <- function(location_code) {
  
  groups_here <- meta_venn %>%
    dplyr::filter(Location_code == location_code) %>%
    dplyr::distinct(Group, Site, Season, Type) %>%
    dplyr::arrange(Site, Season) %>%
    dplyr::pull(Group)
  
  group_sets <- lapply(groups_here, function(g) {
    samples_group <- meta_venn %>%
      dplyr::filter(Group == g) %>%
      dplyr::pull(Sample)
    
    get_core_taxa(
      samples = samples_group,
      abundance_threshold = ABUNDANCE_THRESHOLD,
      prevalence_threshold = PREVALENCE_THRESHOLD
    )
  })
  
  names(group_sets) <- groups_here
  
  return(group_sets)
}

group_sets_by_location <- list(
  B = make_location_group_sets("B"),
  P = make_location_group_sets("P"),
  S = make_location_group_sets("S"),
  V = make_location_group_sets("V"),
  C = make_location_group_sets("C")
)

#### 5) Colours: average between F and R colours for each location ------------

mix_two_colours <- function(col1, col2) {
  rgb1 <- grDevices::col2rgb(col1)
  rgb2 <- grDevices::col2rgb(col2)
  mixed <- rowMeans(cbind(rgb1, rgb2))
  grDevices::rgb(mixed[1], mixed[2], mixed[3], maxColorValue = 255)
}

location_cols <- c(
  B = mix_two_colours(site_cols_map["BF"], site_cols_map["BR"]),
  C = mix_two_colours(site_cols_map["CF"], site_cols_map["CR"]),
  P = mix_two_colours(site_cols_map["PF"], site_cols_map["PR"]),
  S = mix_two_colours(site_cols_map["SF"], site_cols_map["SR"]),
  V = mix_two_colours(site_cols_map["VF"], site_cols_map["VR"])
)

location_cols <- location_cols[location_levels]

#### 6) Create Venn plots -----------------------------------------------------

central_venn_plot <- ggvenn::ggvenn(
  location_sets,
  fill_color = unname(location_cols),
  stroke_size = 0.4,
  set_name_size = 4,
  text_size = 4,
  show_percentage = FALSE
) +
  ggplot2::labs(
    title = "Core genera shared among locations"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 13, face = "bold")
  )

make_group_venn_plot <- function(location_code) {
  
  group_sets <- group_sets_by_location[[location_code]]
  
  location_name <- location_labels[[location_code]]
  
  group_cols <- rep(location_cols[[location_code]], length(group_sets))
  
  ggvenn::ggvenn(
    group_sets,
    fill_color = group_cols,
    stroke_size = 0.35,
    set_name_size = 3.2,
    text_size = 3.2,
    show_percentage = FALSE
  ) +
    ggplot2::labs(
      title = location_name
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 11, face = "bold")
    )
}

venn_B <- make_group_venn_plot("B")
venn_P <- make_group_venn_plot("P")
venn_S <- make_group_venn_plot("S")
venn_V <- make_group_venn_plot("V")
venn_C <- make_group_venn_plot("C")

#### 7) Arrange: central Venn + five peripheral Venns -------------------------

empty_plot <- ggplot2::ggplot() + ggplot2::theme_void()

venn_layout <- "
ABC
DEF
GHI
"

combined_venn_plot <- 
  (venn_B | empty_plot | venn_P) /
  (venn_S | central_venn_plot | venn_V) /
  (empty_plot | venn_C | empty_plot) +
  patchwork::plot_annotation(
    title = "Genus-level core taxa across locations and site-season groups",
    subtitle = "Taxa included when relative abundance was >=1% in at least 80% of samples within each category",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 15, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 11)
    )
  )

combined_venn_plot

#### 8) Save figure -----------------------------------------------------------

ggplot2::ggsave(
  filename = file.path(
    venn_output_dir,
    "Genus_core_taxa_Venn_locations_and_siteSeasonGroups.png"
  ),
  plot = combined_venn_plot,
  width = 14,
  height = 12,
  dpi = 1200
)

ggplot2::ggsave(
  filename = file.path(
    venn_output_dir,
    "Genus_core_taxa_Venn_locations_and_siteSeasonGroups.pdf"
  ),
  plot = combined_venn_plot,
  width = 14,
  height = 12,
  device = cairo_pdf
)

#### 9) Save taxa lists -------------------------------------------------------

save_set_table <- function(set_list, output_file) {
  
  out <- tibble::tibble(
    Category = names(set_list),
    Taxa = lapply(set_list, paste, collapse = "; "),
    n_taxa = lengths(set_list)
  )
  
  readr::write_tsv(out, output_file)
}

save_set_table(
  location_sets,
  file.path(venn_output_dir, "Genus_core_taxa_location_sets.tsv")
)

for (loc in names(group_sets_by_location)) {
  save_set_table(
    group_sets_by_location[[loc]],
    file.path(
      venn_output_dir,
      paste0("Genus_core_taxa_", loc, "_siteSeason_group_sets.tsv")
    )
  )
}

#### 10) Print summary --------------------------------------------------------

cat("\nCentral location Venn taxa numbers:\n")
print(lengths(location_sets))

cat("\nPeripheral group Venn taxa numbers:\n")
print(lapply(group_sets_by_location, lengths))

cat("\nSaved Venn outputs in:\n")
cat(normalizePath(venn_output_dir), "\n")

#### Alternative figure: core genera presence/absence matrix ------------------

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)

#### 1) Location-level presence/absence table --------------------------------

location_presence_table <- tibble::tibble(
  Location = names(location_sets),
  Genus = location_sets
) %>%
  tidyr::unnest(Genus) %>%
  dplyr::distinct(Location, Genus) %>%
  dplyr::mutate(Present = 1)

all_location_genera <- sort(unique(location_presence_table$Genus))

location_presence_matrix <- tidyr::expand_grid(
  Location = names(location_sets),
  Genus = all_location_genera
) %>%
  dplyr::left_join(location_presence_table, by = c("Location", "Genus")) %>%
  dplyr::mutate(
    Present = ifelse(is.na(Present), 0, Present)
  )

genus_order_location <- location_presence_matrix %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarise(
    n_locations = sum(Present),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_locations), Genus) %>%
  dplyr::pull(Genus)

location_presence_matrix <- location_presence_matrix %>%
  dplyr::mutate(
    Genus = factor(Genus, levels = rev(genus_order_location)),
    Location = factor(Location, levels = names(location_sets))
  )

location_count_table <- location_presence_matrix %>%
  dplyr::group_by(Location) %>%
  dplyr::summarise(
    n_core_genera = sum(Present),
    .groups = "drop"
  )

location_presence_plot <- ggplot2::ggplot(
  location_presence_matrix,
  ggplot2::aes(x = Location, y = Genus)
) +
  ggplot2::geom_tile(
    ggplot2::aes(fill = Location, alpha = Present),
    color = "white",
    linewidth = 0.25
  ) +
  ggplot2::scale_fill_manual(
    values = location_cols,
    drop = FALSE
  ) +
  ggplot2::scale_alpha_continuous(
    range = c(0.08, 1),
    breaks = c(0, 1),
    labels = c("Absent", "Present"),
    name = "Core genus"
  ) +
  ggplot2::labs(
    title = "Core genera shared across locations",
    subtitle = "Genus included when relative abundance was ≥1% in at least 80% of samples within each location",
    x = "Location",
    y = "Genus",
    fill = "Location"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 11, face = "bold"),
    axis.text.y = ggplot2::element_text(size = 8, face = "italic"),
    panel.grid = ggplot2::element_blank(),
    legend.position = "right",
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = ggplot2::element_text(hjust = 0.5)
  )

#### 2) Within-location group-level presence/absence table --------------------

group_presence_table <- purrr::imap_dfr(
  group_sets_by_location,
  function(group_sets, location_code) {
    tibble::tibble(
      Location_code = location_code,
      Group = names(group_sets),
      Genus = group_sets
    ) %>%
      tidyr::unnest(Genus)
  }
) %>%
  dplyr::distinct(Location_code, Group, Genus) %>%
  dplyr::mutate(Present = 1)

all_group_genera <- sort(unique(group_presence_table$Genus))

group_presence_matrix <- group_presence_table %>%
  dplyr::group_by(Location_code) %>%
  tidyr::complete(
    Group,
    Genus = all_group_genera,
    fill = list(Present = 0)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Location = dplyr::recode(
      Location_code,
      B = "Bordolona",
      P = "Preghena",
      S = "Sadole",
      V = "Valbiolo",
      C = "Cavaion"
    )
  )

# keep only genera present in at least one group
group_presence_matrix <- group_presence_matrix %>%
  dplyr::group_by(Genus) %>%
  dplyr::filter(sum(Present) > 0) %>%
  dplyr::ungroup()

genus_order_group <- group_presence_matrix %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarise(
    n_groups = sum(Present),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_groups), Genus) %>%
  dplyr::pull(Genus)

group_presence_matrix <- group_presence_matrix %>%
  dplyr::mutate(
    Genus = factor(Genus, levels = rev(genus_order_group)),
    Group = factor(Group)
  )

group_fill_cols <- c(
  B = location_cols["B"],
  P = location_cols["P"],
  S = location_cols["S"],
  V = location_cols["V"],
  C = location_cols["C"]
)

group_presence_plot <- ggplot2::ggplot(
  group_presence_matrix,
  ggplot2::aes(x = Group, y = Genus)
) +
  ggplot2::geom_tile(
    ggplot2::aes(fill = Location_code, alpha = Present),
    color = "white",
    linewidth = 0.2
  ) +
  ggplot2::facet_wrap(~ Location, scales = "free_x", nrow = 1) +
  ggplot2::scale_fill_manual(values = group_fill_cols, drop = FALSE) +
  ggplot2::scale_alpha_continuous(
    range = c(0.06, 1),
    breaks = c(0, 1),
    labels = c("Absent", "Present"),
    name = "Core genus"
  ) +
  ggplot2::labs(
    title = "Core genera across site-season groups",
    subtitle = "Genus included when relative abundance was ≥1% in at least 80% of samples within each group",
    x = "Site-season group",
    y = "Genus",
    fill = "Location"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.text.y = ggplot2::element_text(size = 6, face = "italic"),
    strip.text = ggplot2::element_text(face = "bold", size = 10),
    panel.grid = ggplot2::element_blank(),
    legend.position = "right",
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = ggplot2::element_text(hjust = 0.5)
  )

#### 3) Save separately -------------------------------------------------------

location_presence_plot
group_presence_plot

ggplot2::ggsave(
  filename = file.path(
    venn_output_dir,
    "Genus_core_taxa_location_presence_absence_matrix.png"
  ),
  plot = location_presence_plot,
  width = 7,
  height = 8,
  dpi = 600
)

ggplot2::ggsave(
  filename = file.path(
    venn_output_dir,
    "Genus_core_taxa_location_presence_absence_matrix.pdf"
  ),
  plot = location_presence_plot,
  width = 7,
  height = 8,
  device = cairo_pdf
)

ggplot2::ggsave(
  filename = file.path(
    venn_output_dir,
    "Genus_core_taxa_siteSeason_presence_absence_matrix.png"
  ),
  plot = group_presence_plot,
  width = 13,
  height = 9,
  dpi = 600
)

ggplot2::ggsave(
  filename = file.path(
    venn_output_dir,
    "Genus_core_taxa_siteSeason_presence_absence_matrix.pdf"
  ),
  plot = group_presence_plot,
  width = 13,
  height = 9,
  device = cairo_pdf
)

#### 4) Save tables -----------------------------------------------------------

readr::write_tsv(
  location_presence_matrix,
  file.path(
    venn_output_dir,
    "Genus_core_taxa_location_presence_absence_matrix.tsv"
  )
)

readr::write_tsv(
  group_presence_matrix,
  file.path(
    venn_output_dir,
    "Genus_core_taxa_siteSeason_presence_absence_matrix.tsv"
  )
)

readr::write_tsv(
  location_count_table,
  file.path(
    venn_output_dir,
    "Genus_core_taxa_location_counts.tsv"
  )
)

#### Genus-level Venn diagram by Location, prevalence 75% --------------------

library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggvenn)

venn_output_dir <- file.path(descriptive_output_dir, "genus_level_venn_diagrams")
dir.create(venn_output_dir, showWarnings = FALSE, recursive = TRUE)

ABUNDANCE_THRESHOLD <- 0.001
PREVALENCE_THRESHOLD <- 0.75

#### 1) Prepare genus matrix and metadata ------------------------------------

genus_mat <- as.matrix(genus_relative_abundance)
storage.mode(genus_mat) <- "numeric"
genus_mat[is.na(genus_mat)] <- 0

meta_venn <- sample_metadata %>%
  dplyr::filter(Sample %in% rownames(genus_mat)) %>%
  dplyr::arrange(match(Sample, rownames(genus_mat))) %>%
  dplyr::mutate(
    Sample = as.character(Sample),
    Site = as.character(Site),
    Season = as.character(Season),
    Location_code = substr(Site, 1, 1)
  )

stopifnot(all(meta_venn$Sample == rownames(genus_mat)))

#### 2) Function: core genera per location -----------------------------------

get_core_genera_location <- function(location_code) {
  
  samples_loc <- meta_venn %>%
    dplyr::filter(Location_code == location_code) %>%
    dplyr::pull(Sample)
  
  submat <- genus_mat[samples_loc, , drop = FALSE]
  
  if (ABUNDANCE_THRESHOLD == 0) {
    presence <- submat > 0
  } else {
    presence <- submat >= ABUNDANCE_THRESHOLD
  }
  
  prevalence <- colSums(presence, na.rm = TRUE) / nrow(presence)
  
  core_genera <- names(prevalence)[prevalence >= PREVALENCE_THRESHOLD]
  
  return(core_genera)
}

#### 3) Build location sets ---------------------------------------------------

location_levels <- c("B", "P", "S", "V", "C")

location_labels <- c(
  B = "Bordolona",
  P = "Preghena",
  S = "Sadole",
  V = "Valbiolo",
  C = "Cavaion"
)

location_sets <- lapply(location_levels, get_core_genera_location)
names(location_sets) <- location_labels[location_levels]

#### 4) Colours: intermediate colour between F and R site colours -------------

mix_two_colours <- function(col1, col2) {
  rgb1 <- grDevices::col2rgb(col1)
  rgb2 <- grDevices::col2rgb(col2)
  mixed <- rowMeans(cbind(rgb1, rgb2))
  grDevices::rgb(mixed[1], mixed[2], mixed[3], maxColorValue = 255)
}

location_cols <- c(
  B = mix_two_colours(site_cols_map["BF"], site_cols_map["BR"]),
  P = mix_two_colours(site_cols_map["PF"], site_cols_map["PR"]),
  S = mix_two_colours(site_cols_map["SF"], site_cols_map["SR"]),
  V = mix_two_colours(site_cols_map["VF"], site_cols_map["VR"]),
  C = mix_two_colours(site_cols_map["CF"], site_cols_map["CR"])
)

location_cols <- location_cols[location_levels]

#### 5) Venn plot -------------------------------------------------------------

genus_location_venn_plot <- ggvenn::ggvenn(
  location_sets,
  fill_color = unname(location_cols),
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 4,
  show_percentage = FALSE
) +
  ggplot2::labs(
    title = "Core genera shared among locations",
    subtitle = "Genera included when relative abundance was ≥1% in at least 75% of samples within each location"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)
  )

genus_location_venn_plot

#### 6) Save figure -----------------------------------------------------------

ggplot2::ggsave(
  filename = file.path(
    venn_output_dir,
    "Genus_core_genera_Venn_locations_prev75_abund1.png"
  ),
  plot = genus_location_venn_plot,
  width = 8,
  height = 7,
  dpi = 600
)

ggplot2::ggsave(
  filename = file.path(
    venn_output_dir,
    "Genus_core_genera_Venn_locations_prev75_abund1.pdf"
  ),
  plot = genus_location_venn_plot,
  width = 8,
  height = 7,
  device = cairo_pdf
)

#### 7) Save taxa lists -------------------------------------------------------

location_sets_table <- tibble::tibble(
  Location = names(location_sets),
  n_core_genera = lengths(location_sets),
  Core_genera = sapply(location_sets, paste, collapse = "; ")
)

readr::write_tsv(
  location_sets_table,
  file.path(
    venn_output_dir,
    "Genus_core_genera_location_sets_prev75_abund1.tsv"
  )
)

location_presence_table <- tibble::tibble(
  Location = names(location_sets),
  Genus = location_sets
) %>%
  tidyr::unnest(Genus) %>%
  dplyr::distinct(Location, Genus)

readr::write_tsv(
  location_presence_table,
  file.path(
    venn_output_dir,
    "Genus_core_genera_location_presence_table_prev75_abund1.tsv"
  )
)

cat("\nNumber of core genera per location:\n")
print(location_sets_table %>% dplyr::select(Location, n_core_genera))

cat("\nSaved location-level Venn outputs in:\n")
cat(normalizePath(venn_output_dir), "\n")

#### Family-level Venn diagram by Location, prevalence 75% --------------------

library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggvenn)

venn_output_dir <- file.path(descriptive_output_dir, "family_level_venn_diagrams")
dir.create(venn_output_dir, showWarnings = FALSE, recursive = TRUE)

ABUNDANCE_THRESHOLD <- 0.01
PREVALENCE_THRESHOLD <- 0.75

#### 1) Prepare family matrix and metadata -----------------------------------

family_mat <- as.matrix(family_relative_abundance)
storage.mode(family_mat) <- "numeric"
family_mat[is.na(family_mat)] <- 0

meta_venn <- sample_metadata %>%
  dplyr::filter(Sample %in% rownames(family_mat)) %>%
  dplyr::arrange(match(Sample, rownames(family_mat))) %>%
  dplyr::mutate(
    Sample = as.character(Sample),
    Site = as.character(Site),
    Season = as.character(Season),
    Location_code = substr(Site, 1, 1)
  )

stopifnot(all(meta_venn$Sample == rownames(family_mat)))

#### 2) Function: core families per location ---------------------------------

get_core_families_location <- function(location_code) {
  
  samples_loc <- meta_venn %>%
    dplyr::filter(Location_code == location_code) %>%
    dplyr::pull(Sample)
  
  submat <- family_mat[samples_loc, , drop = FALSE]
  
  if (ABUNDANCE_THRESHOLD == 0) {
    presence <- submat > 0
  } else {
    presence <- submat >= ABUNDANCE_THRESHOLD
  }
  
  prevalence <- colSums(presence, na.rm = TRUE) / nrow(presence)
  
  core_families <- names(prevalence)[prevalence >= PREVALENCE_THRESHOLD]
  
  return(core_families)
}

#### 3) Build location sets ---------------------------------------------------

location_levels <- c("B", "P", "S", "V", "C")

location_labels <- c(
  B = "Bordolona",
  P = "Preghena",
  S = "Sadole",
  V = "Valbiolo",
  C = "Cavaion"
)

location_sets <- lapply(location_levels, get_core_families_location)
names(location_sets) <- location_labels[location_levels]

#### 4) Colours: intermediate colour between F and R site colours -------------

mix_two_colours <- function(col1, col2) {
  rgb1 <- grDevices::col2rgb(col1)
  rgb2 <- grDevices::col2rgb(col2)
  mixed <- rowMeans(cbind(rgb1, rgb2))
  grDevices::rgb(mixed[1], mixed[2], mixed[3], maxColorValue = 255)
}

location_cols <- c(
  B = mix_two_colours(site_cols_map["BF"], site_cols_map["BR"]),
  P = mix_two_colours(site_cols_map["PF"], site_cols_map["PR"]),
  S = mix_two_colours(site_cols_map["SF"], site_cols_map["SR"]),
  V = mix_two_colours(site_cols_map["VF"], site_cols_map["VR"]),
  C = mix_two_colours(site_cols_map["CF"], site_cols_map["CR"])
)

location_cols <- location_cols[location_levels]

#### 5) Venn plot -------------------------------------------------------------

family_location_venn_plot <- ggvenn::ggvenn(
  location_sets,
  fill_color = unname(location_cols),
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 4,
  show_percentage = FALSE
) +
  ggplot2::labs(
    title = "Core families shared among locations",
    subtitle = "Families included when relative abundance was ≥1% in at least 75% of samples within each location"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)
  )

family_location_venn_plot

#### 6) Save figure -----------------------------------------------------------

ggplot2::ggsave(
  filename = file.path(
    venn_output_dir,
    "Family_core_families_Venn_locations_prev75_abund1.png"
  ),
  plot = family_location_venn_plot,
  width = 8,
  height = 7,
  dpi = 600
)

ggplot2::ggsave(
  filename = file.path(
    venn_output_dir,
    "Family_core_families_Venn_locations_prev75_abund1.pdf"
  ),
  plot = family_location_venn_plot,
  width = 8,
  height = 7,
  device = cairo_pdf
)

#### 7) Save taxa lists -------------------------------------------------------

location_sets_table <- tibble::tibble(
  Location = names(location_sets),
  n_core_families = lengths(location_sets),
  Core_families = sapply(location_sets, paste, collapse = "; ")
)

readr::write_tsv(
  location_sets_table,
  file.path(
    venn_output_dir,
    "Family_core_families_location_sets_prev75_abund1.tsv"
  )
)

location_presence_table <- tibble::tibble(
  Location = names(location_sets),
  Family = location_sets
) %>%
  tidyr::unnest(Family) %>%
  dplyr::distinct(Location, Family)

readr::write_tsv(
  location_presence_table,
  file.path(
    venn_output_dir,
    "Family_core_families_location_presence_table_prev75_abund1.tsv"
  )
)

cat("\nNumber of core families per location:\n")
print(location_sets_table %>% dplyr::select(Location, n_core_families))

cat("\nFamilies shared by all locations:\n")
print(Reduce(intersect, location_sets))

cat("\nSaved family-level Venn outputs in:\n")
cat(normalizePath(venn_output_dir), "\n")

#### Family-level UpSet plots by Valley, Season and Site type -----------------
#### Family-level UpSet plot by Valley ---------------------------------------

library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ComplexUpset)

upset_output_dir <- file.path(descriptive_output_dir, "family_level_upset_plots")
dir.create(upset_output_dir, showWarnings = FALSE, recursive = TRUE)

ABUNDANCE_THRESHOLD <- 0.01
PREVALENCE_THRESHOLD <- 0.75

#### 1) Prepare matrix and metadata ------------------------------------------

family_mat <- as.matrix(family_relative_abundance)
storage.mode(family_mat) <- "numeric"
family_mat[is.na(family_mat)] <- 0

meta_upset <- sample_metadata %>%
  dplyr::filter(Sample %in% rownames(family_mat)) %>%
  dplyr::arrange(match(Sample, rownames(family_mat))) %>%
  dplyr::mutate(
    Sample = as.character(Sample),
    Site = as.character(Site),
    Season = dplyr::case_when(
      as.character(Season) == "J" ~ "July",
      as.character(Season) == "S" ~ "September",
      TRUE ~ as.character(Season)
    ),
    Site_type = dplyr::case_when(
      as.character(Site_type) == "F" ~ "Rock-glacier-fed",
      as.character(Site_type) == "R" ~ "Reference",
      TRUE ~ as.character(Site_type)
    ),
    Type = dplyr::case_when(
      as.character(Type) == "F" ~ "Rock-glacier-fed",
      as.character(Type) == "R" ~ "Reference",
      TRUE ~ as.character(Type)
    ),
    Valley = dplyr::case_when(
      as.character(Valley) == "Val Bresimo" ~ "Val Bresimo (B/P)",
      as.character(Valley) == "Val de La Mare" ~ "Val de La Mare (C)",
      as.character(Valley) == "Sadole_Lagorai" ~ "Val Sadole (S)",
      as.character(Valley) == "Valbiolo" ~ "Valbiolo (V)",
      TRUE ~ as.character(Valley)
    )
  )

stopifnot(all(meta_upset$Sample == rownames(family_mat)))

#### 2) Valley colours --------------------------------------------------------

mix_multiple_colours <- function(cols) {
  rgb_mat <- grDevices::col2rgb(cols)
  mixed <- rowMeans(rgb_mat)
  grDevices::rgb(mixed[1], mixed[2], mixed[3], maxColorValue = 255)
}

valley_cols <- c(
  "Val Bresimo (B/P)" = mix_multiple_colours(
    c(
      site_cols_map["BF"], site_cols_map["BR"],
      site_cols_map["PF"], site_cols_map["PR"]
    )
  ),
  "Val de La Mare (C)" = mix_multiple_colours(
    c(site_cols_map["CF"], site_cols_map["CR"])
  ),
  "Val Sadole (S)" = mix_multiple_colours(
    c(site_cols_map["SF"], site_cols_map["SR"])
  ),
  "Valbiolo (V)" = mix_multiple_colours(
    c(site_cols_map["VF"], site_cols_map["VR"])
  )
)

#### 3) Helper function -------------------------------------------------------

get_core_taxa_by_category <- function(category_column,
                                      category_levels = NULL,
                                      abundance_threshold = 0.01,
                                      prevalence_threshold = 0.75) {
  
  if (is.null(category_levels)) {
    category_levels <- sort(unique(meta_upset[[category_column]]))
  }
  
  category_levels <- category_levels[!is.na(category_levels)]
  
  sets <- lapply(category_levels, function(cat) {
    
    samples_cat <- meta_upset %>%
      dplyr::filter(.data[[category_column]] == cat) %>%
      dplyr::pull(Sample)
    
    submat <- family_mat[samples_cat, , drop = FALSE]
    
    if (nrow(submat) == 0) {
      return(character(0))
    }
    
    if (abundance_threshold == 0) {
      presence <- submat > 0
    } else {
      presence <- submat >= abundance_threshold
    }
    
    prevalence <- colSums(presence, na.rm = TRUE) / nrow(presence)
    
    names(prevalence)[prevalence >= prevalence_threshold]
  })
  
  names(sets) <- category_levels
  
  return(sets)
}

#### 4) Build Valley sets -----------------------------------------------------

valley_levels <- c(
  "Val Bresimo (B/P)",
  "Val de La Mare (C)",
  "Val Sadole (S)",
  "Valbiolo (V)"
)

valley_sets <- get_core_taxa_by_category(
  category_column = "Valley",
  category_levels = valley_levels,
  abundance_threshold = ABUNDANCE_THRESHOLD,
  prevalence_threshold = PREVALENCE_THRESHOLD
)

#### 5) Build UpSet input table ----------------------------------------------

all_families <- sort(unique(unlist(valley_sets)))

upset_df <- tibble::tibble(Family = all_families)

for (valley in names(valley_sets)) {
  upset_df[[valley]] <- upset_df$Family %in% valley_sets[[valley]]
}

sets_table <- tibble::tibble(
  Category = names(valley_sets),
  n_core_families = lengths(valley_sets),
  Core_families = sapply(valley_sets, paste, collapse = "; ")
)

readr::write_tsv(
  upset_df,
  file.path(
    upset_output_dir,
    "Family_core_families_UpSet_Valley_prev75_abund1_presence_absence.tsv"
  )
)

readr::write_tsv(
  sets_table,
  file.path(
    upset_output_dir,
    "Family_core_families_UpSet_Valley_prev75_abund1_sets.tsv"
  )
)

#### 6) Manual UpSet-like plot: aligned panels -------------------------------

library(patchwork)

valley_order <- c(
  "Val Bresimo (B/P)",
  "Val de La Mare (C)",
  "Val Sadole (S)",
  "Valbiolo (V)"
)

valley_order_plot <- rev(valley_order)

#### Build intersection table -------------------------------------------------

intersection_table <- upset_df %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    Intersection = paste(
      names(valley_sets)[as.logical(c_across(dplyr::all_of(names(valley_sets))))],
      collapse = " | "
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(Intersection != "") %>%
  dplyr::count(Intersection, name = "intersection_size") %>%
  dplyr::arrange(dplyr::desc(intersection_size), Intersection) %>%
  dplyr::mutate(
    Intersection_ID = paste0("I", dplyr::row_number()),
    Intersection_ID = factor(Intersection_ID, levels = Intersection_ID)
  )

intersection_matrix <- intersection_table %>%
  dplyr::select(Intersection, Intersection_ID) %>%
  tidyr::separate_rows(Intersection, sep = " \\| ") %>%
  dplyr::rename(Valley = Intersection) %>%
  dplyr::mutate(Present = TRUE)

matrix_plot_df <- tidyr::expand_grid(
  Intersection_ID = intersection_table$Intersection_ID,
  Valley = valley_order
) %>%
  dplyr::left_join(
    intersection_matrix,
    by = c("Intersection_ID", "Valley")
  ) %>%
  dplyr::mutate(
    Present = ifelse(is.na(Present), FALSE, Present),
    Intersection_ID = factor(
      Intersection_ID,
      levels = levels(intersection_table$Intersection_ID)
    ),
    Valley = factor(Valley, levels = valley_order_plot)
  )

segment_df <- matrix_plot_df %>%
  dplyr::filter(Present) %>%
  dplyr::mutate(
    Valley_num = as.numeric(Valley)
  ) %>%
  dplyr::group_by(Intersection_ID) %>%
  dplyr::summarise(
    y_min = min(Valley_num),
    y_max = max(Valley_num),
    n_valleys = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_valleys > 1)

#### Top barplot --------------------------------------------------------------

top_bar_df <- intersection_table

top_bar_plot <- ggplot2::ggplot(
  top_bar_df,
  ggplot2::aes(x = Intersection_ID, y = intersection_size)
) +
  ggplot2::geom_col(
    fill = "grey82",
    color = "grey45",
    linewidth = 0.25,
    width = 0.78
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = intersection_size),
    vjust = -0.45,
    size = 3.6
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, max(top_bar_df$intersection_size) * 1.28),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::labs(
    y = "Shared core families",
    x = NULL
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(5, 5, 0, 5)
  )

#### Matrix plot --------------------------------------------------------------

matrix_plot <- ggplot2::ggplot(
  matrix_plot_df,
  ggplot2::aes(x = Intersection_ID, y = Valley)
) +
  ggplot2::geom_tile(
    ggplot2::aes(fill = Valley),
    alpha = 0.10,
    width = 1,
    height = 0.75
  ) +
  ggplot2::geom_segment(
    data = segment_df,
    ggplot2::aes(
      x = Intersection_ID,
      xend = Intersection_ID,
      y = y_min,
      yend = y_max
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
  ggplot2::geom_point(
    data = matrix_plot_df %>% dplyr::filter(!Present),
    ggplot2::aes(x = Intersection_ID, y = Valley),
    color = "grey75",
    fill = "grey88",
    shape = 21,
    size = 3
  ) +
  ggplot2::geom_point(
    data = matrix_plot_df %>% dplyr::filter(Present),
    ggplot2::aes(x = Intersection_ID, y = Valley),
    color = "black",
    fill = "black",
    shape = 21,
    size = 3.2
  ) +
  ggplot2::scale_fill_manual(
    values = valley_cols,
    guide = "none"
  ) +
  ggplot2::labs(
    x = "Valley intersection",
    y = NULL
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 5, 5, 5)
  )

#### Valley label plot --------------------------------------------------------

label_df <- tibble::tibble(
  Valley = factor(valley_order_plot, levels = valley_order_plot),
  label = valley_order_plot
)

valley_label_plot <- ggplot2::ggplot(
  label_df,
  ggplot2::aes(x = 1, y = Valley, label = label)
) +
  ggplot2::geom_text(
    hjust = 1,
    size = 3.7
  ) +
  ggplot2::scale_x_continuous(limits = c(0, 1)) +
  ggplot2::theme_void() +
  ggplot2::theme(
    plot.margin = ggplot2::margin(0, 3, 5, 0)
  )

#### Left set-size barplot ----------------------------------------------------

set_size_df <- sets_table %>%
  dplyr::mutate(
    Category = factor(Category, levels = valley_order_plot)
  )

set_size_plot <- ggplot2::ggplot(
  set_size_df,
  ggplot2::aes(x = n_core_families, y = Category, fill = Category)
) +
  ggplot2::geom_col(
    color = "white",
    linewidth = 0.25,
    width = 0.65
  ) +
  ggplot2::scale_fill_manual(
    values = valley_cols,
    guide = "none"
  ) +
  ggplot2::scale_x_reverse(
    expand = ggplot2::expansion(mult = c(0.05, 0.05))
  ) +
  ggplot2::labs(
    x = "Core families per valley",
    y = NULL
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 5, 5)
  )

#### Empty spacers ------------------------------------------------------------

empty_plot <- ggplot2::ggplot() + ggplot2::theme_void()

#### Combine aligned panels ---------------------------------------------------

#### Combine aligned panels with spacer ---------------------------------------

spacer_plot <- ggplot2::ggplot() + ggplot2::theme_void()

family_valley_upset_plot <- (
  empty_plot + empty_plot + top_bar_plot +
    patchwork::plot_layout(widths = c(0.24, 0.18, 0.58))
) /
  (
    empty_plot + empty_plot + spacer_plot +
      patchwork::plot_layout(widths = c(0.24, 0.18, 0.58))
  ) /
  (
    set_size_plot + valley_label_plot + matrix_plot +
      patchwork::plot_layout(widths = c(0.24, 0.18, 0.58))
  ) +
  patchwork::plot_layout(
    heights = c(0.44, 0.06, 0.50)
  ) +
  patchwork::plot_annotation(
    title = "Core families shared among valleys",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      ),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )
  )

family_valley_upset_plot

#### 7) Save plot -------------------------------------------------------------

ggplot2::ggsave(
  filename = file.path(
    upset_output_dir,
    "Family_core_families_UpSet_Valley_prev75_abund1_manual.png"
  ),
  plot = family_valley_upset_plot,
  width = 14,
  height = 8,
  dpi = 1200
)

ggplot2::ggsave(
  filename = file.path(
    upset_output_dir,
    "Family_core_families_UpSet_Valley_prev75_abund1_manual.pdf"
  ),
  plot = family_valley_upset_plot,
  width = 14,
  height = 8,
  device = cairo_pdf
)

#### 8) Print useful summaries ------------------------------------------------

cat("\nCore families per Valley:\n")
print(sets_table %>% dplyr::select(Category, n_core_families))

cat("\nFamilies shared by all valleys:\n")
print(Reduce(intersect, valley_sets))

cat("\nSaved UpSet outputs in:\n")
cat(normalizePath(upset_output_dir), "\n")


