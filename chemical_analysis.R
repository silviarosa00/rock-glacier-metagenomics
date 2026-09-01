
# =========================================================
# CHEMICAL ANALYSIS FOR SECTION 3.1
# Environmental gradient across alpine spring catchments
# =========================================================

setwd("C:/Users/SiRosa/OneDrive - Scientific Network South Tyrol/Projects/Subsurface/unite")


# ---------------------------------------------------------
# 0) Libraries
# ---------------------------------------------------------

library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(vegan)
library(tibble)
library(FactoMineR)
library(factoextra)

# ---------------------------------------------------------
# 1) Output directory
# ---------------------------------------------------------

dir.create("chemical_analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("chemical_analysis/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("chemical_analysis/figures", showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------
# 2) Input file
# ---------------------------------------------------------
# Change this with your real file name.
# If your file is Excel, use read_excel().
# If your file is csv/tsv, use read_csv() or read_tsv().

chemical_file <- "metadata4.csv"

chemical_raw <- readr::read_delim(chemical_file, delim = ",", show_col_types = FALSE)



# ---------------------------------------------------------
# 3) Clean column names
# ---------------------------------------------------------
# Remove special characters, line breaks and units formatting.
# Example:
# "Cd\nug_L"  -> "Cd_ug_L"
# "Temp_°C"   -> "Temp_C"
# "HCO3_\nmg_L" -> "HCO3_mg_L"

chemical_data <- chemical_raw %>%
  rename_with(~ str_replace_all(.x, "\n", "_")) %>%
  rename_with(~ str_replace_all(.x, "µ", "u")) %>%
  rename_with(~ str_replace_all(.x, "°", "")) %>%
  rename_with(~ str_replace_all(.x, "δ", "d")) %>%
  rename_with(~ str_replace_all(.x, "[^A-Za-z0-9]+", "_")) %>%
  rename_with(~ str_replace_all(.x, "__+", "_")) %>%
  rename_with(~ str_replace_all(.x, "_$", "")) %>%
  rename_with(~ str_replace_all(.x, "^_", ""))

names(chemical_data)


# ---------------------------------------------------------
# 4) Create metadata from sample names
# ---------------------------------------------------------
# Sample code examples:
# SFJ = Sadole, reference spring, July
# SRS = Sadole, rock-glacier spring, September
#
# First letter = location
# Second letter = type
# Third letter = season

chemical_data <- chemical_data %>%
  mutate(
    Location_code = str_sub(Sample, 1, 1),
    Type_code = str_sub(Sample, 2, 2),
    Season_code = str_sub(Sample, 3, 3),
    
    Location = recode(
      Location_code,
      "B" = "Bordolona",
      "P" = "Preghena",
      "V" = "Valbiolo",
      "S" = "Sadole",
      "C" = "Cavaion"
    ),
    
    Valley = case_when(
      Location %in% c("Bordolona", "Preghena") ~ "Val Bresimo",
      Location == "Valbiolo" ~ "Valbiolo",
      Location == "Sadole" ~ "Sadole_Lagorai",
      Location == "Cavaion" ~ "Val de La Mare",
      TRUE ~ NA_character_
    ),
    
    Site_type = recode(
      Type_code,
      "R" = "RG_influenced",
      "F" = "Reference"
    ),
    
    Season = recode(
      Season_code,
      "J" = "July",
      "S" = "September"
    ),
    
    RG_status = case_when(
      Sample %in% c("CRJ", "CRS") ~ "Intact",
      Sample %in% c("BRJ", "BRS", "PRJ", "PRS", "SRJ", "SRS") ~ "Pseudo_relict",
      Sample %in% c("VRJ", "VRS") ~ "Relict",
      Site_type == "Reference" ~ "Reference",
      TRUE ~ NA_character_
    )
  )


# ---------------------------------------------------------
# 5) Convert chemical columns to numeric
# ---------------------------------------------------------
# "<LOD" and "NA" strings are converted to real NA.
# This is conservative.
# If you want to replace <LOD with half LOD, we can modify this.

chemical_data_numeric <- chemical_data %>%
  mutate(
    across(
      where(is.character),
      ~ na_if(.x, "<LOD")
    )
  ) %>%
  mutate(
    across(
      where(is.character),
      ~ na_if(.x, "NA")
    )
  )

# Convert all possible chemical columns to numeric, while keeping metadata as text.

metadata_columns <- c(
  "Sample", "Note", "Location_code", "Type_code", "Season_code",
  "Location", "Valley", "Site_type", "Season", "RG_status"
)

chemical_data_numeric <- chemical_data_numeric %>%
  mutate(
    across(
      -any_of(metadata_columns),
      ~ suppressWarnings(as.numeric(.x))
    )
  )


# ---------------------------------------------------------
# 6) Define variables for analysis
# ---------------------------------------------------------

basic_environmental_variables <- c(
  "EC_uS_cm",
  "pH",
  "Altitude",
  "Temp_C"
)

major_ion_variables <- c(
  "F_mg_L",
  "Cl_mg_L",
  "NO2_mg_L",
  "Br_mg_L",
  "NO3_mg_L",
  "PO43_mg_L",
  "SO42_mg_L",
  "HCO3_mg_L"
)

trace_element_variables <- c(
  "Al_ug_L",
  "As_ug_L",
  "B_ug_L",
  "Ba_ug_L",
  "Be_ug_L",
  "Ca_mg_L",
  "Cd_ug_L",
  "Ce_ug_L",
  "Co_ug_L",
  "Cr_ug_L",
  "Cu_ug_L",
  "Fe_ug_L",
  "K_mg_L",
  "La_ug_L",
  "Li_ug_L",
  "Mg_mg_L",
  "Mn_ug_L",
  "Mo_ug_L",
  "Na_mg_L",
  "Ni_ug_L",
  "Pb_ug_L",
  "Rb_ug_L",
  "Sb_ug_L",
  "Se_ug_L",
  "Si_mg_L",
  "Sr_ug_L",
  "Ti_ug_L",
  "U_ug_L",
  "V_ug_L",
  "Zn_ug_L"
)

chemical_variables_for_analysis <- c(
  basic_environmental_variables,
  major_ion_variables,
  trace_element_variables
)

chemical_variables_for_analysis <- chemical_variables_for_analysis[
  chemical_variables_for_analysis %in% names(chemical_data_numeric)
]

print(chemical_variables_for_analysis)

# ---------------------------------------------------------
# 4) Create metadata from sample names
# ---------------------------------------------------------
# Sample code examples:
# BFJ = Bordolona reference, July
# BRS = Bordolona RG-influenced, September
#
# First letter = location
# Second letter = type/site
# Third letter = season
# Base = first two letters, e.g. BF, BR, CF, CR

base_cols_map <- c(
  BF = "#1b9e77", BR = "#a6dba0",
  CF = "#d95f02", CR = "#fdb863",
  PF = "#7570b3", PR = "#b2abd2",
  SF = "#1f9ac2", SR = "#a6dce7",
  VF = "#e7298a", VR = "#f2b2d4"
)

season_shape_map <- c(
  July = 16,
  September = 17
)

chemical_data <- chemical_data %>%
  mutate(
    Base = str_sub(Sample, 1, 2),
    Location_code = str_sub(Sample, 1, 1),
    Type_code = str_sub(Sample, 2, 2),
    Season_code = str_sub(Sample, 3, 3),
    
    Location = recode(
      Location_code,
      "B" = "Bordolona",
      "P" = "Preghena",
      "V" = "Valbiolo",
      "S" = "Sadole",
      "C" = "Cavaion"
    ),
    
    Valley = case_when(
      Location %in% c("Bordolona", "Preghena") ~ "Val Bresimo",
      Location == "Valbiolo" ~ "Valbiolo",
      Location == "Sadole" ~ "Sadole_Lagorai",
      Location == "Cavaion" ~ "Val de La Mare",
      TRUE ~ NA_character_
    ),
    
    Site_type = recode(
      Type_code,
      "R" = "RG_influenced",
      "F" = "Reference"
    ),
    
    Season = recode(
      Season_code,
      "J" = "July",
      "S" = "September"
    ),
    
    RG_status = case_when(
      Base == "CR" ~ "Intact",
      Base %in% c("BR", "PR", "SR") ~ "Pseudo_relict",
      Base == "VR" ~ "Relict",
      Site_type == "Reference" ~ "Reference",
      TRUE ~ NA_character_
    )
  )

# ---------------------------------------------------------
# 7) Summary table by sampling site
# ---------------------------------------------------------

chemical_summary_by_sample <- chemical_data_numeric %>%
  select(
    Sample,
    Location,
    Valley,
    Site_type,
    RG_status,
    Season,
    all_of(chemical_variables_for_analysis)
  )

write_tsv(
  chemical_summary_by_sample,
  "chemical_analysis/tables/chemical_summary_by_sample.tsv"
)


# ---------------------------------------------------------
# 8) Summary table by Location and Site_type
# ---------------------------------------------------------
# This is useful for Supplementary Table.

chemical_summary_by_location_type <- chemical_data_numeric %>%
  group_by(Location, Valley, Site_type, RG_status) %>%
  summarise(
    across(
      all_of(chemical_variables_for_analysis),
      list(
        median = ~ median(.x, na.rm = TRUE),
        min = ~ min(.x, na.rm = TRUE),
        max = ~ max(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

write_tsv(
  chemical_summary_by_location_type,
  "chemical_analysis/tables/chemical_summary_by_location_type.tsv"
)


# ---------------------------------------------------------
# 9) Simple summary of key variables for main text
# ---------------------------------------------------------
# This table gives compact values for EC, pH, temperature and key metals.

key_variables_for_text <- c(
  "EC_uS_cm",
  "pH",
  "Temp_C",
  "Cd_ug_L",
  "Ce_ug_L",
  "Cu_ug_L",
  "Fe_ug_L",
  "Ni_ug_L",
  "Pb_ug_L",
  "Sb_ug_L",
  "Zn_ug_L"
)

key_variables_for_text <- key_variables_for_text[
  key_variables_for_text %in% names(chemical_data_numeric)
]

key_chemical_summary <- chemical_data_numeric %>%
  group_by(Location) %>%
  summarise(
    across(
      all_of(key_variables_for_text),
      list(
        median = ~ median(.x, na.rm = TRUE),
        min = ~ min(.x, na.rm = TRUE),
        max = ~ max(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

write_tsv(
  key_chemical_summary,
  "chemical_analysis/tables/key_chemical_summary_by_location.tsv"
)


# ---------------------------------------------------------
# 10) PCA of hydrochemical variables
# ---------------------------------------------------------
# Purpose:
# Check whether hydrochemistry separates mostly by location/valley,
# rather than simply by RG vs reference.

chemical_matrix_for_pca <- chemical_data_numeric %>%
  select(all_of(chemical_variables_for_analysis)) %>%
  mutate(across(everything(), ~ log1p(.x))) %>%
  mutate(across(everything(), ~ ifelse(is.nan(.x), NA, .x)))

# Remove variables with too many missing values or zero variance.

chemical_matrix_for_pca_clean <- chemical_matrix_for_pca %>%
  select(
    where(~ sum(!is.na(.x)) >= 0.7 * nrow(chemical_matrix_for_pca))
  ) %>%
  select(
    where(~ sd(.x, na.rm = TRUE) > 0)
  )

# Replace remaining missing values with variable median.

chemical_matrix_for_pca_imputed <- chemical_matrix_for_pca_clean %>%
  mutate(
    across(
      everything(),
      ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)
    )
  )

chemical_pca_result <- prcomp(
  chemical_matrix_for_pca_imputed,
  center = TRUE,
  scale. = TRUE
)

chemical_pca_scores <- as.data.frame(chemical_pca_result$x) %>%
  bind_cols(
    chemical_data_numeric %>%
      select(Sample, Base, Location, Valley, Site_type, Season, RG_status)
  )

chemical_pca_loadings <- as.data.frame(chemical_pca_result$rotation) %>%
  rownames_to_column("Variable")

chemical_pca_variance <- tibble(
  PC = paste0("PC", seq_along(chemical_pca_result$sdev)),
  Variance_percent = 100 * chemical_pca_result$sdev^2 / sum(chemical_pca_result$sdev^2)
)

write_tsv(
  chemical_pca_scores,
  "chemical_analysis/tables/chemical_PCA_scores.tsv"
)

write_tsv(
  chemical_pca_loadings,
  "chemical_analysis/tables/chemical_PCA_loadings.tsv"
)

write_tsv(
  chemical_pca_variance,
  "chemical_analysis/tables/chemical_PCA_variance.tsv"
)


# ---------------------------------------------------------
# ---------------------------------------------------------
# 11) PCA plot with PCoA-style colours
# ---------------------------------------------------------
# Colours follow the Base colour map used in the other ordination plots.
# Shapes indicate season.

pca_axis1_percent <- round(chemical_pca_variance$Variance_percent[1], 1)
pca_axis2_percent <- round(chemical_pca_variance$Variance_percent[2], 1)

chemical_pca_plot_base <- ggplot(
  chemical_pca_scores,
  aes(
    x = PC1,
    y = PC2,
    colour = Base,
    shape = Season
  )
) +
  geom_point(size = 3.5, alpha = 0.9) +
  scale_colour_manual(values = base_cols_map, drop = FALSE) +
  scale_shape_manual(values = season_shape_map, drop = FALSE) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = paste0("PC1 (", pca_axis1_percent, "%)"),
    y = paste0("PC2 (", pca_axis2_percent, "%)"),
    colour = "Base",
    shape = "Season",
    title = "Hydrochemical PCA"
  )



ggsave(
  "chemical_analysis/figures/chemical_PCA_by_base_season.png",
  chemical_pca_plot_base,
  width = 7,
  height = 5,
  dpi = 300
)



# ---------------------------------------------------------
# 12) PCA plot coloured by RG status
# ---------------------------------------------------------

chemical_pca_plot_rg_status <- ggplot(
  chemical_pca_scores,
  aes(
    x = PC1,
    y = PC2,
    colour = RG_status,
    shape = Season
  )
) +
  geom_point(size = 3, alpha = 0.9) +
  theme_bw(base_size = 12) +
  labs(
    x = paste0("PC1 (", pca_axis1_percent, "%)"),
    y = paste0("PC2 (", pca_axis2_percent, "%)"),
    colour = "RG status",
    shape = "Season",
    title = "Hydrochemical PCA"
  )

ggsave(
  "chemical_analysis/figures/chemical_PCA_by_RG_status.pdf",
  chemical_pca_plot_rg_status,
  width = 7,
  height = 5
)

ggsave(
  "chemical_analysis/figures/chemical_PCA_by_RG_status.png",
  chemical_pca_plot_rg_status,
  width = 7,
  height = 5,
  dpi = 300
)


# ---------------------------------------------------------
# 13) PERMANOVA on hydrochemical profiles
# ---------------------------------------------------------
# Purpose:
# Test whether chemical composition is structured by Location,
# Site_type and Season.

chemical_distance_euclidean <- dist(chemical_matrix_for_pca_imputed, method = "euclidean")

chemical_permanova_location <- adonis2(
  chemical_distance_euclidean ~ Location + Site_type + Season,
  data = chemical_data_numeric,
  permutations = 999,
  by = "margin"
)

chemical_permanova_valley <- adonis2(
  chemical_distance_euclidean ~ Valley + Site_type + Season,
  data = chemical_data_numeric,
  permutations = 999,
  by = "margin"
)

write_tsv(
  as.data.frame(chemical_permanova_location) %>%
    rownames_to_column("Term"),
  "chemical_analysis/tables/PERMANOVA_chemical_Location_SiteType_Season.tsv"
)

write_tsv(
  as.data.frame(chemical_permanova_valley) %>%
    rownames_to_column("Term"),
  "chemical_analysis/tables/PERMANOVA_chemical_Valley_SiteType_Season.tsv"
)


# ---------------------------------------------------------
# 14) PERMDISP for hydrochemical Location groups
# ---------------------------------------------------------
# Purpose:
# Check whether PERMANOVA may be influenced by different dispersion.

chemical_betadisper_location <- betadisper(
  chemical_distance_euclidean,
  chemical_data_numeric$Location
)

chemical_betadisper_location_anova <- anova(chemical_betadisper_location)
chemical_betadisper_location_permutation <- permutest(
  chemical_betadisper_location,
  permutations = 999
)

write_tsv(
  as.data.frame(chemical_betadisper_location_anova) %>%
    rownames_to_column("Term"),
  "chemical_analysis/tables/PERMDISP_chemical_Location_anova.tsv"
)

write_tsv(
  as.data.frame(chemical_betadisper_location_permutation$tab) %>%
    rownames_to_column("Term"),
  "chemical_analysis/tables/PERMDISP_chemical_Location_permutation.tsv"
)


# ---------------------------------------------------------
# 15) Boxplots of key environmental variables
# ---------------------------------------------------------


variables_for_boxplots <- c(
  "EC_uS_cm",
  "pH",
  "Temp_C",
  "SO42_mg_L",
  "Cd_ug_L",
  "Ce_ug_L",
  "Ni_ug_L",
  "Cu_ug_L",
  "Zn_ug_L",
  "Pb_ug_L",
  "Sb_ug_L",
  "As_ug_L",
  "Fe_ug_L",
  "Mn_ug_L",
  "U_ug_L",
  "V_ug_L"
)

variables_for_boxplots <- variables_for_boxplots[
  variables_for_boxplots %in% names(chemical_data_numeric)
]

print(variables_for_boxplots)


chemical_long_for_boxplots <- chemical_data_numeric %>%
  select(Sample, Location, Site_type, Season, all_of(variables_for_boxplots)) %>%
  pivot_longer(
    cols = all_of(variables_for_boxplots),
    names_to = "Variable",
    values_to = "Value"
  )
chemical_long_for_boxplots <- chemical_long_for_boxplots %>%
  mutate(
    Variable = str_replace_all(Variable, "_ug_L", " (ug/L)"),
    Variable = str_replace_all(Variable, "_mg_L", " (mg/L)"),
    Variable = str_replace_all(Variable, "_uS_cm", " (uS/cm)"),
    Variable = str_replace_all(Variable, "_C", " (°C)")
  )

chemical_long_for_boxplots <- chemical_long_for_boxplots %>%
  mutate(
    Variable = factor(
      Variable,
      levels = c(
        "Temp (°C)",
        "pH",
        "EC (uS/cm)",
        "SO42 (mg/L)",
        "Cd (ug/L)",
        "Ce (ug/L)",
        "Pb (ug/L)",
        "Sb (ug/L)",
        "Ni (ug/L)",
        "U (ug/L)",
        "Cu (ug/L)",
        "Zn (ug/L)",
        "As (ug/L)",
        "Fe (ug/L)",
        "Mn (ug/L)",
        "V (ug/L)"
      )
    )
  )
chemical_boxplot_by_location <- ggplot(
  chemical_long_for_boxplots,
  aes(
    x = Location,
    y = Value,
    fill = Site_type
  )
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(
    aes(shape = Season),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    size = 1.8,
    alpha = 0.8
  ) +
  facet_wrap(~ Variable, scales = "free_y", ncol=4) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Location",
    y = "Value",
    fill = "Site type",
    shape = "Season"
  )

ggsave(
  "chemical_analysis/figures/key_chemical_variables_by_location.pdf",
  chemical_boxplot_by_location,
  width = 10,
  height = 7
)

ggsave(
  "chemical_analysis/figures/key_chemical_variables_by_location.png",
  chemical_boxplot_by_location,
  width = 10,
  height = 7,
  dpi = 300
)


# ---------------------------------------------------------
# 16) Kruskal-Wallis tests for key variables by Location
# ---------------------------------------------------------

kruskal_location_results <- lapply(
  variables_for_boxplots,
  function(current_variable) {
    current_test <- kruskal.test(
      chemical_data_numeric[[current_variable]] ~ chemical_data_numeric$Location
    )
    
    tibble(
      Variable = current_variable,
      Test = "Kruskal_Wallis_by_Location",
      Statistic = as.numeric(current_test$statistic),
      df = as.numeric(current_test$parameter),
      p_value = current_test$p.value
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    p_adjusted_BH = p.adjust(p_value, method = "BH")
  )

write_tsv(
  kruskal_location_results,
  "chemical_analysis/tables/Kruskal_key_variables_by_Location.tsv"
)


# ---------------------------------------------------------
# 17) Wilcoxon tests for key variables by Site_type
# ---------------------------------------------------------

wilcoxon_site_type_results <- lapply(
  variables_for_boxplots,
  function(current_variable) {
    
    current_data <- chemical_data_numeric %>%
      select(Site_type, all_of(current_variable)) %>%
      filter(!is.na(.data[[current_variable]]))
    
    if (length(unique(current_data$Site_type)) == 2) {
      current_test <- wilcox.test(
        current_data[[current_variable]] ~ current_data$Site_type,
        exact = FALSE
      )
      
      tibble(
        Variable = current_variable,
        Test = "Wilcoxon_by_Site_type",
        Statistic = as.numeric(current_test$statistic),
        p_value = current_test$p.value
      )
    } else {
      tibble(
        Variable = current_variable,
        Test = "Wilcoxon_by_Site_type",
        Statistic = NA_real_,
        p_value = NA_real_
      )
    }
  }
) %>%
  bind_rows() %>%
  mutate(
    p_adjusted_BH = p.adjust(p_value, method = "BH")
  )

write_tsv(
  wilcoxon_site_type_results,
  "chemical_analysis/tables/Wilcoxon_key_variables_by_SiteType.tsv"
)


# ---------------------------------------------------------
# 18) Wilcoxon tests for key variables by Season
# ---------------------------------------------------------

wilcoxon_season_results <- lapply(
  variables_for_boxplots,
  function(current_variable) {
    
    current_data <- chemical_data_numeric %>%
      select(Season, all_of(current_variable)) %>%
      filter(!is.na(.data[[current_variable]]))
    
    if (length(unique(current_data$Season)) == 2) {
      current_test <- wilcox.test(
        current_data[[current_variable]] ~ current_data$Season,
        exact = FALSE
      )
      
      tibble(
        Variable = current_variable,
        Test = "Wilcoxon_by_Season",
        Statistic = as.numeric(current_test$statistic),
        p_value = current_test$p.value
      )
    } else {
      tibble(
        Variable = current_variable,
        Test = "Wilcoxon_by_Season",
        Statistic = NA_real_,
        p_value = NA_real_
      )
    }
  }
) %>%
  bind_rows() %>%
  mutate(
    p_adjusted_BH = p.adjust(p_value, method = "BH")
  )

write_tsv(
  wilcoxon_season_results,
  "chemical_analysis/tables/Wilcoxon_key_variables_by_Season.tsv"
)


# ---------------------------------------------------------
# 19) Optional: correlation matrix among selected chemical variables
# ---------------------------------------------------------

chemical_correlation_matrix <- chemical_data_numeric %>%
  select(all_of(variables_for_boxplots)) %>%
  cor(
    use = "pairwise.complete.obs",
    method = "spearman"
  )

write_tsv(
  as.data.frame(chemical_correlation_matrix) %>%
    rownames_to_column("Variable"),
  "chemical_analysis/tables/Spearman_correlation_key_chemical_variables.tsv"
)


# ---------------------------------------------------------
# 20) Save cleaned full table
# ---------------------------------------------------------

write_tsv(
  chemical_data_numeric,
  "chemical_analysis/tables/chemical_data_cleaned_with_metadata.tsv"
)


# ---------------------------------------------------------
# 21) Print main outputs in console
# ---------------------------------------------------------

print("PERMANOVA: chemical profiles by Location + Site_type + Season")
print(chemical_permanova_location)

print("PERMANOVA: chemical profiles by Valley + Site_type + Season")
print(chemical_permanova_valley)

print("Kruskal-Wallis tests by Location")
print(kruskal_location_results)

print("Wilcoxon tests by Site_type")
print(wilcoxon_site_type_results)

print("Wilcoxon tests by Season")
print(wilcoxon_season_results)


################################################################
# =========================================================
# 16S qPCR ANALYSIS
# =========================================================

library(readxl)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(tibble)

dir.create("chemical_analysis/16S_qPCR", showWarnings = FALSE, recursive = TRUE)
dir.create("chemical_analysis/16S_qPCR/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("chemical_analysis/16S_qPCR/figures", showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------
# 1) Read clean qPCR table
# ---------------------------------------------------------
qpcr_file <- "16s_log.xlsx"

qpcr_16S_data <- read_excel(
  qpcr_file,
  sheet = 1,
  col_names = FALSE
)
# ---------------------------------------------------------
# 2) Clean column names
# ---------------------------------------------------------

qpcr_16S_data <- qpcr_16S_data %>%
  rename(
    Sample = 1,
    Site = 2,
    R1_log10_16S_copies_ng_DNA = 3,
    R2_log10_16S_copies_ng_DNA = 4
  )

# ---------------------------------------------------------
# 3) Calculate mean log10 16S value
# ---------------------------------------------------------
# Since your table already contains log10 values, this gives
# one qPCR value per sample.

qpcr_16S_data <- qpcr_16S_data %>%
  mutate(
    Sample = as.character(Sample),
    Site = as.character(Site),
    R1_log10_16S_copies_ng_DNA = as.numeric(R1_log10_16S_copies_ng_DNA),
    R2_log10_16S_copies_ng_DNA = as.numeric(R2_log10_16S_copies_ng_DNA),
    
    mean_log10_16S_copies_ng_DNA = rowMeans(
      cbind(R1_log10_16S_copies_ng_DNA, R2_log10_16S_copies_ng_DNA),
      na.rm = TRUE
    )
  )

# ---------------------------------------------------------
# 4) Add metadata from Sample code
# ---------------------------------------------------------

qpcr_16S_data <- qpcr_16S_data %>%
  mutate(
    Base = str_sub(Sample, 1, 2),
    Location_code = str_sub(Sample, 1, 1),
    Type_code = str_sub(Sample, 2, 2),
    Season_code = str_sub(Sample, 3, 3),
    
    Location = recode(
      Location_code,
      "B" = "Bordolona",
      "P" = "Preghena",
      "V" = "Valbiolo",
      "S" = "Sadole",
      "C" = "Cavaion"
    ),
    
    Site_type = recode(
      Type_code,
      "R" = "RG_influenced",
      "F" = "Reference"
    ),
    
    Season = recode(
      Season_code,
      "J" = "July",
      "S" = "September"
    )
  )

# ---------------------------------------------------------
# 5) Save clean table
# ---------------------------------------------------------

write_tsv(
  qpcr_16S_data,
  "chemical_analysis/16S_qPCR/tables/16S_qPCR_clean_with_metadata.tsv"
)
# ---------------------------------------------------------
# 2) Add metadata from sample code
# ---------------------------------------------------------

base_cols_map <- c(
  BF = "#1b9e77", BR = "#a6dba0",
  CF = "#d95f02", CR = "#fdb863",
  PF = "#7570b3", PR = "#b2abd2",
  SF = "#1f9ac2", SR = "#a6dce7",
  VF = "#e7298a", VR = "#f2b2d4"
)

season_shape_map <- c(
  July = 16,
  September = 17
)

qpcr_16S_data <- qpcr_16S_data %>%
  mutate(
    Base = str_sub(Sample, 1, 2),
    Location_code = str_sub(Sample, 1, 1),
    Type_code = str_sub(Sample, 2, 2),
    Season_code = str_sub(Sample, 3, 3),
    
    Location = recode(
      Location_code,
      "B" = "Bordolona",
      "P" = "Preghena",
      "V" = "Valbiolo",
      "S" = "Sadole",
      "C" = "Cavaion"
    ),
    
    Valley = case_when(
      Location %in% c("Bordolona", "Preghena") ~ "Val Bresimo",
      Location == "Valbiolo" ~ "Valbiolo",
      Location == "Sadole" ~ "Sadole_Lagorai",
      Location == "Cavaion" ~ "Val de La Mare",
      TRUE ~ NA_character_
    ),
    
    Site_type = recode(
      Type_code,
      "R" = "RG_influenced",
      "F" = "Reference"
    ),
    
    Season = recode(
      Season_code,
      "J" = "July",
      "S" = "September"
    )
  )

# ---------------------------------------------------------
# 3) Calculate mean log10 16S value
# ---------------------------------------------------------
# The input table already contains log10-transformed qPCR values.
# Technical qPCR replicates are averaged to obtain one value per sample.

qpcr_16S_data <- qpcr_16S_data %>%
  mutate(
    R1_log10_16S_copies_ng_DNA = as.numeric(R1_log10_16S_copies_ng_DNA),
    R2_log10_16S_copies_ng_DNA = as.numeric(R2_log10_16S_copies_ng_DNA),
    
    log10_16S_copies_per_ng_DNA = rowMeans(
      cbind(R1_log10_16S_copies_ng_DNA, R2_log10_16S_copies_ng_DNA),
      na.rm = TRUE
    )
  )

# ---------------------------------------------------------
# 4) Overall summary
# ---------------------------------------------------------

qpcr_16S_overall_summary <- qpcr_16S_data %>%
  summarise(
    n_samples = n(),
    median_log10_16S = median(log10_16S_copies_per_ng_DNA, na.rm = TRUE),
    min_log10_16S = min(log10_16S_copies_per_ng_DNA, na.rm = TRUE),
    max_log10_16S = max(log10_16S_copies_per_ng_DNA, na.rm = TRUE)
  )

write_tsv(
  qpcr_16S_overall_summary,
  "chemical_analysis/16S_qPCR/tables/16S_qPCR_overall_summary.tsv"
)

# ---------------------------------------------------------
# 5) Summary by Base
# ---------------------------------------------------------

qpcr_16S_summary_by_base <- qpcr_16S_data %>%
  group_by(Base, Location, Site_type) %>%
  summarise(
    n_samples = n(),
    median_log10_16S = median(log10_16S_copies_per_ng_DNA, na.rm = TRUE),
    min_log10_16S = min(log10_16S_copies_per_ng_DNA, na.rm = TRUE),
    max_log10_16S = max(log10_16S_copies_per_ng_DNA, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(
  qpcr_16S_summary_by_base,
  "chemical_analysis/16S_qPCR/tables/16S_qPCR_summary_by_base.tsv"
)

# ---------------------------------------------------------
# 6) Statistical tests
# ---------------------------------------------------------

qpcr_16S_kruskal_location <- kruskal.test(
  log10_16S_copies_per_ng_DNA ~ Location,
  data = qpcr_16S_data
)

qpcr_16S_kruskal_valley <- kruskal.test(
  log10_16S_copies_per_ng_DNA ~ Valley,
  data = qpcr_16S_data
)

qpcr_16S_wilcoxon_site_type <- wilcox.test(
  log10_16S_copies_per_ng_DNA ~ Site_type,
  data = qpcr_16S_data,
  exact = FALSE
)

qpcr_16S_wilcoxon_season <- wilcox.test(
  log10_16S_copies_per_ng_DNA ~ Season,
  data = qpcr_16S_data,
  exact = FALSE
)

qpcr_16S_test_results <- tibble(
  Test = c(
    "Kruskal_Wallis_Location",
    "Kruskal_Wallis_Valley",
    "Wilcoxon_Site_type",
    "Wilcoxon_Season"
  ),
  Statistic = c(
    as.numeric(qpcr_16S_kruskal_location$statistic),
    as.numeric(qpcr_16S_kruskal_valley$statistic),
    as.numeric(qpcr_16S_wilcoxon_site_type$statistic),
    as.numeric(qpcr_16S_wilcoxon_season$statistic)
  ),
  p_value = c(
    qpcr_16S_kruskal_location$p.value,
    qpcr_16S_kruskal_valley$p.value,
    qpcr_16S_wilcoxon_site_type$p.value,
    qpcr_16S_wilcoxon_season$p.value
  )
) %>%
  mutate(p_adjusted_BH = p.adjust(p_value, method = "BH"))

write_tsv(
  qpcr_16S_test_results,
  "chemical_analysis/16S_qPCR/tables/16S_qPCR_tests.tsv"
)

# ---------------------------------------------------------
# 7) Plot by Base
# ---------------------------------------------------------

qpcr_16S_plot_by_base <- ggplot(
  qpcr_16S_data,
  aes(
    x = Base,
    y = log10_16S_copies_per_ng_DNA,
    colour = Base,
    shape = Season
  )
) +
  geom_boxplot(
    aes(group = Base),
    outlier.shape = NA,
    alpha = 0.25
  ) +
  geom_point(
    size = 3,
    position = position_jitter(width = 0.12, height = 0),
    alpha = 0.9
  ) +
  scale_colour_manual(values = base_cols_map, drop = FALSE) +
  scale_shape_manual(values = season_shape_map, drop = FALSE) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Base",
    y = expression(log[10]~"16S rRNA gene copies per ng DNA"),
    colour = "Base",
    shape = "Season"
  )

ggsave(
  "chemical_analysis/16S_qPCR/figures/16S_qPCR_log10_by_base.pdf",
  qpcr_16S_plot_by_base,
  width = 7,
  height = 5
)

ggsave(
  "chemical_analysis/16S_qPCR/figures/16S_qPCR_log10_by_base.png",
  qpcr_16S_plot_by_base,
  width = 7,
  height = 5,
  dpi = 300
)

# ---------------------------------------------------------
# 8) Print results
# ---------------------------------------------------------

print(qpcr_16S_overall_summary)
print(qpcr_16S_test_results)







