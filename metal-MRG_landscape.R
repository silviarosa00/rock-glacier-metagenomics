#### 0) Packages --------------------------------------------------------------
pkgs <- c("dplyr","readr","tidyr","stringr","ggplot2","tibble","effsize")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(to_install)) install.packages(to_install, repos="https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))
dir.create("tables", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

#### 1) Read MRG per16S -------------------------------------------------------
mrg_file <- "mrg_per16s_T1T2_renamed.tsv"
mrg_wide <- read_tsv(mrg_file, show_col_types = FALSE)

mrg_long <- mrg_wide %>%
  pivot_longer(cols = -BacMet_ID, names_to = "Sample", values_to = "Abundance") %>%
  mutate(
    Abundance = as.numeric(Abundance),
    Abundance = ifelse(is.na(Abundance), 0, Abundance),
    BacMet_ID_clean = str_extract(BacMet_ID, "^BAC\\d+")
  )

#### 2) Metadata from Sample name --------------------------------------------
meta <- mrg_long %>%
  distinct(Sample) %>%
  mutate(
    Base   = str_sub(Sample, 1, 2),
    Site   = str_sub(Sample, 1, 1),
    Season_code = str_sub(Sample, 3, 3),
    Season = recode(Season_code, "J"="July", "S"="September", .default = Season_code),
    Replicate = str_extract(Sample, "\\d+$"),
    Type = ifelse(str_sub(Base, 2, 2) == "R", "RG", "non-RG")
  )

mrg_long <- mrg_long %>% left_join(meta, by = "Sample")

#### 3) Read BacMet mapping + build gene -> metal group -----------------------
map_in <- "bacmet.txt"
bacmap <- read_tsv(map_in, show_col_types = FALSE) %>%
  transmute(
    BacMet_ID = as.character(BacMet_ID),
    Compound  = as.character(Compound)
  )

# Metalli che hai nel metadata (usati come dizionario)
metals <- c("Li","Be","Al","Ti","V","Cr","Mn","Fe","Co","Ni","Zn","Cu","As","Se",
            "Rb","Sr","Mo","Cd","Sb","Ba","La","Ce","Pb","U")

# Dizionario per match nel campo Compound (nomi + simboli)
metal_terms <- tibble::tribble(
  ~Metal, ~regex,
  "Li", "lithium|\\bli\\b",
  "Be", "beryllium|\\bbe\\b",
  "Al", "aluminium|aluminum|\\bal\\b",
  "Ti", "titanium|\\bti\\b",
  "V",  "vanadium|\\bv\\b",
  "Cr", "chromium|\\bcr\\b",
  "Mn", "manganese|\\bmn\\b",
  "Fe", "iron|\\bfe\\b",
  "Co", "cobalt|\\bco\\b",
  "Ni", "nickel|\\bni\\b",
  "Zn", "zinc|\\bzn\\b",
  "Cu", "copper|\\bcu\\b",
  "As", "arsenic|arsenate|arsenite|\\bas\\b",
  "Se", "selenium|selenite|selenate|\\bse\\b",
  "Rb", "rubidium|\\brb\\b",
  "Sr", "strontium|\\bsr\\b",
  "Mo", "molybdenum|\\bmo\\b",
  "Cd", "cadmium|\\bcd\\b",
  "Sb", "antimony|\\bsb\\b",
  "Ba", "barium|\\bba\\b",
  "La", "lanthanum|\\bla\\b",
  "Ce", "cerium|\\bce\\b",
  "Pb", "lead|\\bpb\\b",
  "U",  "uranium|\\bu\\b"
)

# Per ogni gene: quali metalli compaiono nel testo Compound?
gene_metal_long <- bacmap %>%
  mutate(comp_low = str_to_lower(ifelse(is.na(Compound), "", Compound))) %>%
  tidyr::crossing(metal_terms) %>%
  mutate(hit = str_detect(comp_low, regex)) %>%
  filter(hit) %>%
  select(BacMet_ID, Metal) %>%
  distinct()

# Se un gene matcha più metalli, lo contiamo in tutte le categorie (è ok per composizione)
# ma per evitare doppio conteggio nel totale metal, calcoleremo il totale a parte (is_metal).
gene_is_metal <- gene_metal_long %>%
  distinct(BacMet_ID) %>%
  mutate(is_metal = TRUE)

write_tsv(gene_metal_long, "tables/bacmet_gene_to_metal_long.tsv")
write_tsv(gene_is_metal,   "tables/bacmet_gene_is_metal.tsv")

#### 4) Total MRG + Metal MRG metrics per sample ------------------------------
alpha_tbl <- mrg_long %>%
  group_by(Sample, Base, Site, Type, Season, Replicate) %>%
  summarise(
    Total_MRG_per16S = sum(Abundance, na.rm = TRUE),
    MRG_richness     = sum(Abundance > 0, na.rm = TRUE),
    .groups = "drop"
  )

# Metal_MRG_per16S (senza doppio conteggio): somma abbondanze dei geni che sono "metal-associated"
mrg_long2 <- mrg_long %>%
  left_join(gene_is_metal, by = c("BacMet_ID_clean" = "BacMet_ID")) %>%
  mutate(is_metal = ifelse(is.na(is_metal), FALSE, is_metal))

metal_by_sample <- mrg_long2 %>%
  group_by(Sample, Base, Site, Type, Season, Replicate) %>%
  summarise(
    Metal_MRG_per16S = sum(Abundance[is_metal], na.rm = TRUE),
    .groups = "drop"
  )

alpha_tbl <- alpha_tbl %>%
  left_join(metal_by_sample, by = c("Sample","Base","Site","Type","Season","Replicate")) %>%
  mutate(
    Metal_MRG_per16S = ifelse(is.na(Metal_MRG_per16S), 0, Metal_MRG_per16S),
    Fraction_metal   = ifelse(Total_MRG_per16S > 0, Metal_MRG_per16S / Total_MRG_per16S, NA_real_)
  )

write_tsv(alpha_tbl, "tables/alpha_total_richness_metalmetrics.tsv")

#### 5) Stacked bar: composition by METAL (replicates aggregated) -------------
# Abundance per sample per metal (qui c'è doppio conteggio se gene matcha >1 metallo: va bene per "associazioni",
# ma per una composizione più conservativa, possiamo normalizzare sul totale metal-by-sample ottenuto sopra.)
metal_comp_sample <- mrg_long %>%
  inner_join(gene_metal_long, by = c("BacMet_ID_clean" = "BacMet_ID")) %>%
  group_by(Sample, Base, Site, Type, Season, Replicate, Metal) %>%
  summarise(Metal_abund = sum(Abundance, na.rm = TRUE), .groups="drop")

# Aggrega replicati: mediana per Base×Season×Metal
metal_comp_base <- metal_comp_sample %>%
  group_by(Base, Site, Type, Season, Metal) %>%
  summarise(Metal_abund_med = median(Metal_abund, na.rm=TRUE), .groups="drop")

# Converte in composizione relativa sul totale metal (per Base×Season)
metal_comp_base <- metal_comp_base %>%
  group_by(Base, Season) %>%
  mutate(Metal_rel = Metal_abund_med / sum(Metal_abund_med, na.rm=TRUE)) %>%
  ungroup()

# Top N metalli + Other (per leggibilità)
N <- 8
topN <- metal_comp_base %>%
  group_by(Metal) %>%
  summarise(mean_rel = mean(Metal_rel, na.rm=TRUE), .groups="drop") %>%
  arrange(desc(mean_rel)) %>%
  slice_head(n=N) %>%
  pull(Metal)

metal_plot <- metal_comp_base %>%
  mutate(Metal2 = ifelse(Metal %in% topN, Metal, "Other")) %>%
  group_by(Base, Site, Type, Season, Metal2) %>%
  summarise(Metal_rel = sum(Metal_rel, na.rm=TRUE), .groups="drop") %>%
  group_by(Base, Season) %>%
  mutate(Metal_rel = Metal_rel / sum(Metal_rel, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(Metal2 = factor(Metal2, levels = c(setdiff(sort(unique(Metal2)), "Other"), "Other")))

p <- ggplot(metal_plot, aes(x = Base, y = Metal_rel, fill = Metal2)) +
  geom_col(color="white", linewidth=0.2) +
  facet_wrap(~Season, nrow=1) +
  scale_y_continuous(labels = function(x) paste0(round(100*x), "%")) +
  labs(x="Base (site)", y="Metal-associated MRG composition (median across replicates)", fill="Element") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("figures/MRG_metal_compounds_TopN_byBaseSeason.png", p, width=11, height=4.5, dpi=300)

write_tsv(metal_plot, "tables/metal_composition_byBaseSeason_TopN.tsv")

#### 6) Season difference stats for Metal_MRG_per16S + Fraction_metal ----------
wilcox_es <- function(df, response, group_var) {
  x <- df[[response]]
  g <- df[[group_var]]
  lv <- unique(g)
  lv <- lv[!is.na(lv)]
  stopifnot(length(lv) == 2)
  
  g1 <- lv[1]; g2 <- lv[2]
  x1 <- x[g == g1]; x2 <- x[g == g2]
  
  wt <- wilcox.test(x1, x2, exact = FALSE)
  cd <- effsize::cliff.delta(x1, x2)$estimate
  
  tibble::tibble(
    metric = response,
    group = group_var,
    level1 = as.character(g1),
    level2 = as.character(g2),
    n1 = length(x1),
    n2 = length(x2),
    median1 = median(x1, na.rm=TRUE),
    median2 = median(x2, na.rm=TRUE),
    p_value = wt$p.value,
    cliffs_delta = as.numeric(cd)
  )
}

stats_metal_season <- dplyr::bind_rows(
  wilcox_es(alpha_tbl, "Metal_MRG_per16S", "Season"),
  wilcox_es(alpha_tbl, "Fraction_metal",   "Season")
) %>%
  mutate(p_adj = p.adjust(p_value, method="BH"))

write_tsv(stats_metal_season, "tables/stats_metalmetrics_Season.tsv")

#########################################################################
#########################################################################
#########################################################################
#########################################################################
################# PROVA CON PIU METALLI 
#########################################################################
#### 0) Packages --------------------------------------------------------------
pkgs <- c("dplyr","readr","tidyr","stringr","ggplot2","tibble","effsize","RColorBrewer")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(to_install)) install.packages(to_install, repos="https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))
dir.create("tables", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

#### 1) Read MRG per16S -------------------------------------------------------
mrg_file <- "mrg_per16s_T1T2_renamed.tsv"
mrg_wide <- read_tsv(mrg_file, show_col_types = FALSE)

mrg_long <- mrg_wide %>%
  pivot_longer(cols = -BacMet_ID, names_to = "Sample", values_to = "Abundance") %>%
  mutate(
    Abundance = as.numeric(Abundance),
    Abundance = ifelse(is.na(Abundance), 0, Abundance),
    BacMet_ID_clean = str_extract(BacMet_ID, "^BAC\\d+")
  )

#### 2) Metadata from Sample name --------------------------------------------
meta <- mrg_long %>%
  distinct(Sample) %>%
  mutate(
    Base   = str_sub(Sample, 1, 2),
    Site   = str_sub(Sample, 1, 1),
    Season_code = str_sub(Sample, 3, 3),
    Season = recode(Season_code, "J"="July", "S"="September", .default = Season_code),
    Replicate = str_extract(Sample, "\\d+$"),
    Type = ifelse(str_sub(Base, 2, 2) == "R", "RG", "non-RG")
  )

mrg_long <- mrg_long %>% left_join(meta, by = "Sample")

#### 3) Read BacMet mapping ---------------------------------------------------
map_in <- "bacmet.txt"   # <-- il tuo mapping con colonna Compound
bacmap <- read_tsv(map_in, show_col_types = FALSE) %>%
  transmute(
    BacMet_ID = as.character(BacMet_ID),
    Compound  = as.character(Compound),
    comp_low  = str_to_lower(ifelse(is.na(Compound), "", Compound))
  )

#### 4) BIG dictionary: metals/metalloids (names + symbols + common anions) ----
# Include molti elementi comuni in BacMet: Hg/Ag/Sn/Au/W ecc.
# Uso confini di parola per i simboli per evitare match casuali.
metal_terms_focus <- tibble::tribble(
  ~Metal, ~regex,
  "Cu", "copper",
  "Zn", "zinc",
  "Ni", "nickel",
  "Co", "cobalt",
  "Cd", "cadmium",
  "As", "arsenic|arsenate|arsenite",
  "Cr", "chromium|chromate",
  "Mn", "manganese",
  "Fe", "iron",
  "Pb", "lead",
  "Hg", "mercury|organomercury|organo-mercury",
  "Ag", "silver",
  "Sn", "tin|organotin|organo-tin",
  "Sb", "antimony",
  "Se", "selenium|selenate|selenite",
  "Mo", "molybdenum|molybdate",
  "W",  "tungsten|wolfram",
  "Au", "gold",
  "Al", "aluminium|aluminum",
  "V",  "vanadium|vanadate",
  "U",  "uranium|uranyl"
)

####To avoid inflating the “metal-associated” category with ubiquitous macronutrients and common ions (e.g., Na, K, Ca, Mg, Cl, S),
####metal-associated MRGs were defined using BacMet Compound annotations restricted to trace metals/metalloids typically 
####implicated in metal tolerance and co-selection (e.g., Cu, Zn, Ni, Co, Cd, As, Cr, Pb, Hg, Ag, Sn, Sb, Se, Mo).
.

#### 5) gene -> metal mapping (all metals) ------------------------------------
gene_metal_long_all <- bacmap %>%
  select(BacMet_ID, comp_low) %>%
  tidyr::crossing(metal_terms_focus)%>%
  mutate(hit = str_detect(comp_low, regex)) %>%
  filter(hit) %>%
  select(BacMet_ID, Metal) %>%
  distinct()

gene_is_metal_all <- gene_metal_long_all %>%
  distinct(BacMet_ID) %>%
  mutate(is_metal = TRUE)

write_tsv(gene_metal_long_all, "tables/bacmet_gene_to_metal_long_ALL.tsv")
write_tsv(gene_is_metal_all,   "tables/bacmet_gene_is_metal_ALL.tsv")

#### 6) Total + Metal metrics (ALL) -------------------------------------------
alpha_tbl <- mrg_long %>%
  group_by(Sample, Base, Site, Type, Season, Replicate) %>%
  summarise(
    Total_MRG_per16S = sum(Abundance, na.rm = TRUE),
    MRG_richness     = sum(Abundance > 0, na.rm = TRUE),
    .groups = "drop"
  )

mrg_long_all <- mrg_long %>%
  left_join(gene_is_metal_all, by = c("BacMet_ID_clean" = "BacMet_ID")) %>%
  mutate(is_metal = ifelse(is.na(is_metal), FALSE, is_metal))

metal_by_sample_all <- mrg_long_all %>%
  group_by(Sample, Base, Site, Type, Season, Replicate) %>%
  summarise(
    Metal_MRG_per16S_all = sum(Abundance[is_metal], na.rm = TRUE),
    .groups = "drop"
  )

alpha_tbl <- alpha_tbl %>%
  left_join(metal_by_sample_all, by = c("Sample","Base","Site","Type","Season","Replicate")) %>%
  mutate(
    Metal_MRG_per16S_all = ifelse(is.na(Metal_MRG_per16S_all), 0, Metal_MRG_per16S_all),
    Fraction_metal_all   = ifelse(Total_MRG_per16S > 0, Metal_MRG_per16S_all / Total_MRG_per16S, NA_real_)
  )

write_tsv(alpha_tbl, "tables/alpha_total_richness_metalmetrics_ALL.tsv")

#### 7) Stacked composition by metal (ALL), replicates aggregated -------------
metal_comp_sample_all <- mrg_long %>%
  inner_join(gene_metal_long_all, by = c("BacMet_ID_clean" = "BacMet_ID")) %>%
  group_by(Sample, Base, Site, Type, Season, Replicate, Metal) %>%
  summarise(Metal_abund = sum(Abundance, na.rm = TRUE), .groups="drop")

metal_comp_base_all <- metal_comp_sample_all %>%
  group_by(Base, Site, Type, Season, Metal) %>%
  summarise(Metal_abund_med = median(Metal_abund, na.rm=TRUE), .groups="drop") %>%
  group_by(Base, Season) %>%
  mutate(Metal_rel = Metal_abund_med / sum(Metal_abund_med, na.rm=TRUE)) %>%
  ungroup()

# Top N + Other metals
N <- 10
topN <- metal_comp_base_all %>%
  group_by(Metal) %>%
  summarise(mean_rel = mean(Metal_rel, na.rm=TRUE), .groups="drop") %>%
  arrange(desc(mean_rel)) %>%
  slice_head(n=N) %>%
  pull(Metal)

metal_plot_all <- metal_comp_base_all %>%
  mutate(Metal2 = ifelse(Metal %in% topN, Metal, "Other metals")) %>%
  group_by(Base, Site, Type, Season, Metal2) %>%
  summarise(Metal_rel = sum(Metal_rel, na.rm=TRUE), .groups="drop") %>%
  group_by(Base, Season) %>%
  mutate(Metal_rel = Metal_rel / sum(Metal_rel, na.rm=TRUE)) %>%
  ungroup()

# Metti Other metals in cima e in grigio
levs <- c(sort(setdiff(unique(metal_plot_all$Metal2), "Other metals")), "Other metals")
metal_plot_all <- metal_plot_all %>% mutate(Metal2 = factor(Metal2, levels = levs))


k <- length(levs)
pal <- RColorBrewer::brewer.pal(max(3, min(12, k)), "Set3")
names(pal) <- levs
pal["Other metals"] <- "grey70"

p <- ggplot(metal_plot_all, aes(x = Base, y = Metal_rel, fill = Metal2)) +
  geom_col(color="white", linewidth=0.2) +
  facet_wrap(~Season, nrow=1) +
  scale_y_continuous(labels = function(x) paste0(round(100*x), "%")) +
  scale_fill_manual(values = pal) +
  labs(x="Site", y="Metal-associated MRG composition (median across replicates)", fill="Metal") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("figures/MRG_metal_compounds_ALL_TopN_byBaseSeason1.png", p, width=12, height=5, dpi=600)
write_tsv(metal_plot_all, "tables/metal_composition_byBaseSeason_ALL_TopN.tsv")

#### 8) Season stats (ALL) ----------------------------------------------------
wilcox_es <- function(df, response, group_var) {
  x <- df[[response]]
  g <- df[[group_var]]
  lv <- unique(g); lv <- lv[!is.na(lv)]
  stopifnot(length(lv) == 2)
  g1 <- lv[1]; g2 <- lv[2]
  x1 <- x[g == g1]; x2 <- x[g == g2]
  wt <- wilcox.test(x1, x2, exact = FALSE)
  cd <- effsize::cliff.delta(x1, x2)$estimate
  tibble(
    metric = response, group = group_var,
    level1 = as.character(g1), level2 = as.character(g2),
    n1 = length(x1), n2 = length(x2),
    median1 = median(x1, na.rm=TRUE),
    median2 = median(x2, na.rm=TRUE),
    p_value = wt$p.value,
    cliffs_delta = as.numeric(cd)
  )
}

stats_metal_season_all <- bind_rows(
  wilcox_es(alpha_tbl, "Metal_MRG_per16S_all", "Season"),
  wilcox_es(alpha_tbl, "Fraction_metal_all",   "Season")
) %>% mutate(p_adj = p.adjust(p_value, method="BH"))

write_tsv(stats_metal_season_all, "tables/stats_metalmetrics_Season_ALL.tsv")

###################################################################
#################  TABELLA METALLI SELEZIONATI ###########################
##################################################################
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

dir.create("tables", showWarnings = FALSE)

# ---- input files ----
mrg_file <- "mrg_per16s_T1T2_renamed.tsv"
map_in   <- "bacmet.txt"

# ---- 1) leggi MRG per16S (wide: gene x sample) ----
mrg_wide <- read_tsv(mrg_file, show_col_types = FALSE)

# pulisci ID (BAC0001...) per join con BacMet
mrg_wide <- mrg_wide %>%
  mutate(BacMet_ID_clean = str_extract(BacMet_ID, "^BAC\\d+"))

# ---- 2) leggi BacMet mapping (Compound) ----
bacmap <- read_tsv(map_in, show_col_types = FALSE) %>%
  transmute(
    BacMet_ID = as.character(BacMet_ID),
    comp_low  = str_to_lower(ifelse(is.na(Compound), "", as.character(Compound)))
  )

# ---- 3) dizionario metals-of-interest (NO Na/K/Ca/Mg etc.) ----
metal_terms_focus <- tibble::tribble(
  ~Metal, ~regex,
  "Cu", "copper",
  "Zn", "zinc",
  "Ni", "nickel",
  "Co", "cobalt",
  "Cd", "cadmium",
  "As", "arsenic|arsenate|arsenite",
  "Cr", "chromium|chromate",
  "Mn", "manganese",
  "Fe", "iron",
  "Pb", "lead",
  "Hg", "mercury|organomercury|organo-mercury",
  "Ag", "silver",
  "Sn", "tin|organotin|organo-tin",
  "Sb", "antimony",
  "Se", "selenium|selenate|selenite",
  "Mo", "molybdenum|molybdate",
  "W",  "tungsten|wolfram",
  "Au", "gold",
  "Al", "aluminium|aluminum",
  "V",  "vanadium|vanadate",
  "U",  "uranium|uranyl"
)

# ---- 4) trova i geni BacMet che matchano almeno un metallo ----
gene_is_metal_focus <- bacmap %>%
  tidyr::crossing(metal_terms_focus) %>%
  mutate(hit = str_detect(comp_low, regex)) %>%
  filter(hit) %>%
  distinct(BacMet_ID) %>%
  mutate(is_metal = TRUE)

write_tsv(gene_is_metal_focus, "tables/bacmet_gene_is_metal_focus.tsv")

# ---- 5) filtra la matrice MRG per16S tenendo solo geni metal-associated ----
total_metal_mrg_per16s <- mrg_wide %>%
  semi_join(gene_is_metal_focus, by = c("BacMet_ID_clean" = "BacMet_ID")) %>%
  select(-BacMet_ID_clean)

# salva
write_tsv(total_metal_mrg_per16s, "tables/total_metal_mrg_per16s.tsv")

# (opzionale) stampa dimensioni
cat("Saved tables/total_metal_mrg_per16s.tsv with",
    nrow(total_metal_mrg_per16s), "genes x",
    (ncol(total_metal_mrg_per16s)-1), "samples\n")
##########################################################################
#######################################################################
##################barplot top geni #####################################

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)

dir.create("tables", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# --- inputs ---
mat_file <- "tables/total_metal_mrg_per16s.tsv"
map_in   <- "bacmet.txt"

# 1) Read gene x sample matrix (metal-only)
mat <- read_tsv(mat_file, show_col_types = FALSE)
stopifnot("BacMet_ID" %in% names(mat))

# Clean BacMet ID for join (BAC0001)
mat <- mat %>%
  mutate(BacMet_ID_clean = str_extract(BacMet_ID, "^BAC\\d+"))

# 2) Read BacMet mapping to get short gene names
# Assumes bacmet.txt has columns: BacMet_ID, Gene_name
bacmap <- read_tsv(map_in, show_col_types = FALSE) %>%
  transmute(
    BacMet_ID = as.character(BacMet_ID),
    Gene_name = as.character(Gene_name)
  ) %>%
  distinct()

# 3) Long format + metadata from Sample name
long <- mat %>%
  pivot_longer(cols = -c(BacMet_ID, BacMet_ID_clean),
               names_to = "Sample", values_to = "Abundance") %>%
  mutate(
    Abundance = as.numeric(Abundance),
    Abundance = ifelse(is.na(Abundance), 0, Abundance),
    Base   = str_sub(Sample, 1, 2),
    Site   = str_sub(Sample, 1, 1),
    Season_code = str_sub(Sample, 3, 3),
    Season = recode(Season_code, "J"="July", "S"="September", .default = Season_code),
    Replicate = str_extract(Sample, "\\d+$"),
    Type = ifelse(str_sub(Base, 2, 2) == "R", "RG", "non-RG")
  ) %>%
  left_join(bacmap, by = c("BacMet_ID_clean" = "BacMet_ID")) %>%
  mutate(
    # fallback se Gene_name manca
    Gene_label = ifelse(is.na(Gene_name) | Gene_name == "", BacMet_ID_clean, Gene_name)
  )

# 4) Aggregate replicates: median abundance per Base×Season×Gene
gene_base <- long %>%
  group_by(Base, Site, Type, Season, Gene_label) %>%
  summarise(gene_abund_med = median(Abundance, na.rm = TRUE), .groups = "drop")

# 5) Convert to relative composition within each Base×Season (sul totale metal-genes)
gene_base <- gene_base %>%
  group_by(Base, Season) %>%
  mutate(gene_rel = gene_abund_med / sum(gene_abund_med, na.rm = TRUE)) %>%
  ungroup()

# 6) Pick Top N genes by mean relative abundance across Base×Season
N <- 15
topN <- gene_base %>%
  group_by(Gene_label) %>%
  summarise(mean_rel = mean(gene_rel, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_rel)) %>%
  slice_head(n = N) %>%
  pull(Gene_label)

# 7) Collapse others -> "Other genes" and renormalize
gene_plot <- gene_base %>%
  mutate(Gene2 = ifelse(Gene_label %in% topN, Gene_label, "Other genes")) %>%
  group_by(Base, Site, Type, Season, Gene2) %>%
  summarise(gene_rel = sum(gene_rel, na.rm = TRUE), .groups = "drop") %>%
  group_by(Base, Season) %>%
  mutate(gene_rel = gene_rel / sum(gene_rel, na.rm = TRUE)) %>%
  ungroup()

# Put "Other genes" on top (last level)
levs <- c(sort(setdiff(unique(gene_plot$Gene2), "Other genes")), "Other genes")
gene_plot <- gene_plot %>% mutate(Gene2 = factor(Gene2, levels = levs))

# install.packages("colorspace")
# install.packages("Polychrome")
# install.packages("colorspace")
install.packages(c("randomcoloR","colorspace"), repos="https://cloud.r-project.org")
library(randomcoloR)
library(colorspace)

k <- length(levs)

# colori ben separati
pal_raw <- randomcoloR::distinctColorPalette(k)

# rendili pastello (più chiari) ma restano distinti
pal_soft <- colorspace::lighten(pal_raw, amount = 0.25)

pal <- setNames(pal_soft, levs)
pal["Other genes"] <- "grey70"



# 8) Plot stacked bar by Base×Season (replicates united)
p <- ggplot(gene_plot, aes(x = Base, y = gene_rel, fill = Gene2)) +
  geom_col(color = "white", linewidth = 0.2) +
  facet_wrap(~Season, nrow = 1) +
  scale_y_continuous(labels = function(x) paste0(round(100*x), "%")) +
  scale_fill_manual(values = pal) +
  labs(x = "Site",
       y = "Relative abundance (median across replicates)",
       fill = "MRG genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave("figures/MRG_metalGenes_Top15_byBaseSeason.png", p, width = 12, height = 5, dpi = 600)

write_tsv(gene_plot, "tables/metalGenes_rel_median_byBaseSeason_Top15.tsv")
write_tsv(tibble(Gene = topN), "tables/metalGenes_top15_list.tsv")

library(readr)
library(dplyr)
library(tidyr)
library(stringr)

in_file <- "tables/total_metal_mrg_per16s.tsv"   # <-- la tua matrice metal-MRG

# 1) Read matrix (genes x samples)
mat_wide <- read_tsv(in_file, show_col_types = FALSE)
stopifnot(ncol(mat_wide) > 2)

gene_col <- names(mat_wide)[1]   # prima colonna = gene ID
sample_cols <- setdiff(names(mat_wide), gene_col)

# 2) Basic overview on genes/features
n_genes_total <- nrow(mat_wide)

# gene present in >=1 sample
present_any <- mat_wide %>%
  mutate(any_present = rowSums(across(all_of(sample_cols), ~ as.numeric(.) > 0), na.rm = TRUE) > 0)

n_genes_present_any <- sum(present_any$any_present, na.rm = TRUE)

# 3) Long format to compute per-sample totals/richness
mat_long <- mat_wide %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "Sample", values_to = "Abundance") %>%
  mutate(Abundance = as.numeric(Abundance),
         Abundance = ifelse(is.na(Abundance), 0, Abundance))

sample_summary <- mat_long %>%
  group_by(Sample) %>%
  summarise(
    Total_metalMRG_per16S = sum(Abundance, na.rm = TRUE),
    Richness_metalMRG     = sum(Abundance > 0, na.rm = TRUE),
    .groups = "drop"
  )

# 4) Summary across samples
overall_summary <- sample_summary %>%
  summarise(
    n_samples = n(),
    mean_total = mean(Total_metalMRG_per16S),
    median_total = median(Total_metalMRG_per16S),
    sd_total = sd(Total_metalMRG_per16S),
    mean_richness = mean(Richness_metalMRG),
    median_richness = median(Richness_metalMRG),
    sd_richness = sd(Richness_metalMRG)
  )

# 5) (Optional) Top genes by mean abundance across samples
top_genes <- mat_long %>%
  group_by(!!sym(gene_col)) %>%
  summarise(mean_abund = mean(Abundance, na.rm = TRUE),
            prev = mean(Abundance > 0, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 10)

# 6) Save tables
write_tsv(sample_summary, "tables/metalMRG_sample_summary.tsv")
write_tsv(overall_summary, "tables/metalMRG_overall_summary.tsv")
write_tsv(top_genes, "tables/metalMRG_top10_genes_mean.tsv")

# 7) Print a ready-to-use overview sentence (fill numbers automatically)
cat("\nOVERVIEW:\n")
cat(sprintf(
  "The metal-associated MRG matrix contained %d gene features (%d detected in at least one sample; n=%d samples). Across samples, the total metal-MRG load (per 16S) had a median of %.1f (mean %.1f), and metal-MRG richness had a median of %.0f genes (mean %.1f).\n",
  n_genes_total,
  n_genes_present_any,
  overall_summary$n_samples,
  overall_summary$median_total,
  overall_summary$mean_total,
  overall_summary$median_richness,
  overall_summary$mean_richness
))

##################################################################
##################################################################
##################################################################
library(dplyr)
library(readr)
library(stringr)

infile <- "tables/metal_composition_byBaseSeason_ALL_TopN.tsv"
stopifnot(file.exists(infile))
df <- read_tsv(infile, show_col_types = FALSE)

# df atteso: Base, Site, Type, Season, Metal2, Metal_rel
df <- df %>%
  mutate(
    Metal2 = as.character(Metal2),
    Season = as.character(Season)
  )

# 1) Top metals usati nel grafico (tutti i livelli tranne "Other metals")
top_metals <- df %>%
  distinct(Metal2) %>%
  filter(Metal2 != "Other metals") %>%
  pull(Metal2) %>%
  sort()

# 2) Summary globale: mediana e range percentuale per ciascun metallo (su tutti i Base×Season)
sum_global <- df %>%
  filter(Metal2 != "Other metals") %>%
  group_by(Metal2) %>%
  summarise(
    median_pct = 100 * median(Metal_rel, na.rm = TRUE),
    min_pct    = 100 * min(Metal_rel, na.rm = TRUE),
    max_pct    = 100 * max(Metal_rel, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_pct))

# 3) Top metal per ciascun Base×Season (quello con % più alta)
top_by_group <- df %>%
  filter(Metal2 != "Other metals") %>%
  group_by(Base, Season) %>%
  slice_max(order_by = Metal_rel, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(pct = 100 * Metal_rel) %>%
  arrange(Season, Base)

# 4) Costruisci un testo "paper-ready"
txt_lines <- c(
  sprintf("Top metals in Figure 3 (as defined in the BacMet compound mapping): %s.",
          paste(top_metals, collapse = ", ")),
  "",
  "Across Base×season groups, the median (range) relative contribution of each top metal was:",
  paste0(
    "- ", sum_global$Metal2, ": ",
    sprintf("%.1f", sum_global$median_pct), "% (",
    sprintf("%.1f", sum_global$min_pct), "–",
    sprintf("%.1f", sum_global$max_pct), "%)"
  ),
  "",
  "Dominant metal by Base×season (highest within-group percentage):",
  paste0(
    "- ", top_by_group$Season, " / ", top_by_group$Base, ": ",
    top_by_group$Metal2, " (", sprintf("%.1f", top_by_group$pct), "%)"
  )
)

writeLines(txt_lines, "tables/Fig3_topMetals_summary.txt")
cat("Wrote: tables/Fig3_topMetals_summary.txt\n")

library(dplyr)
library(tidyr)
library(vegan)

# Parti da metal_comp_sample_all (che hai già creato)
# Sample, Base, Site, Type, Season, Replicate, Metal, Metal_abund

# 1) matrice Sample x Metal (abund) e poi relative within-sample
mat <- metal_comp_sample_all %>%
  group_by(Sample, Metal) %>%
  summarise(abund = sum(Metal_abund, na.rm=TRUE), .groups="drop") %>%
  group_by(Sample) %>%
  mutate(rel = abund / sum(abund, na.rm=TRUE)) %>%
  ungroup() %>%
  select(Sample, Metal, rel) %>%
  pivot_wider(names_from = Metal, values_from = rel, values_fill = 0)

X <- as.matrix(mat[,-1])
rownames(X) <- mat$Sample

meta_s <- meta %>% distinct(Sample, Base, Site, Type, Season)

# 2) PERMANOVA (Bray-Curtis)
set.seed(1)
adonis2(X ~ Season + Type, data = meta_s, permutations = 999, method = "bray", strata = meta_s$Site)
##################################################################
##################################################################
##################################################################
library(dplyr)
library(tidyr)
library(vegan)
library(readr)

# ---- (il tuo codice per costruire mat, X, meta_s) ----
# mat <- ...
# X <- ...
# meta_s <- ...

# 1) PERMANOVA
set.seed(1)
perm_met <- adonis2(X ~ Season + Type,
                    data = meta_s,
                    permutations = 999,
                    method = "bray",
                    strata = meta_s$Site)

# 2) Salva tabella "pulita" (tsv) + output testuale (txt)
perm_df <- as.data.frame(perm_met) %>%
  tibble::rownames_to_column("Term")

write_tsv(perm_df, "tables/PERMANOVA_metal_comp_Season+Type_strataSite.tsv")

# 3) Crea un TXT leggibile (paper-friendly)
txt <- c(
  "PERMANOVA (adonis2) on metal-associated MRG composition (metal categories)",
  "-----------------------------------------------------------------------",
  "Input: per-sample relative composition of metal categories (within-sample proportions).",
  "Distance: Bray-Curtis",
  "Model: Season + Type",
  "Permutation scheme: 999 permutations, constrained within Site (strata = Site).",
  "",
  "Link to Figure (Top-N stacked barplot):",
  "  The figure displays Top-N metals + 'Other metals' for visualization,",
  "  whereas this PERMANOVA is run on the full metal-category composition (recommended for inference).",
  "",
  sprintf("N samples used: %d", nrow(meta_s)),
  "",
  "Results table (Term / Df / SumOfSqs / R2 / F / Pr(>F)):",
  paste(capture.output(print(perm_met)), collapse = "\n")
)

writeLines(txt, "tables/PERMANOVA_metal_comp_Season+Type_strataSite.txt")
cat("Saved:\n",
    " - tables/PERMANOVA_metal_comp_Season+Type_strataSite.tsv\n",
    " - tables/PERMANOVA_metal_comp_Season+Type_strataSite.txt\n", sep="")

d <- vegdist(X, method="bray")
anova(betadisper(d, meta_s$Season))
anova(betadisper(d, meta_s$Type))
##################################################################
##################################################################
##################################################################
library(dplyr)
library(tidyr)
library(vegan)
library(readr)
library(tibble)

# matrice Sample x Metal (relative within-sample)
mat <- metal_comp_sample_all %>%
  group_by(Sample, Metal) %>%
  summarise(abund = sum(Metal_abund, na.rm=TRUE), .groups="drop") %>%
  group_by(Sample) %>%
  mutate(rel = abund / sum(abund, na.rm=TRUE)) %>%
  ungroup() %>%
  select(Sample, Metal, rel) %>%
  pivot_wider(names_from = Metal, values_from = rel, values_fill = 0)

X <- as.matrix(mat[,-1])
rownames(X) <- mat$Sample

meta_s <- meta %>% distinct(Sample, Base, Site, Type, Season) %>%
  filter(Sample %in% rownames(X))
X <- X[meta_s$Sample, , drop=FALSE]

pairwise_within_site_type <- function(site_letter) {
  sub_meta <- meta_s %>% filter(Site == site_letter)
  
  # serve che nel sito ci siano entrambi i tipi
  if (n_distinct(sub_meta$Type) < 2) return(NULL)
  
  sub_X <- X[sub_meta$Sample, , drop=FALSE]
  
  set.seed(1)
  fit <- adonis2(sub_X ~ Type,
                 data = sub_meta,
                 method = "bray",
                 permutations = 999,
                 strata = sub_meta$Season)  # controlla Season
  
  out <- as.data.frame(fit) %>% rownames_to_column("Term")
  type_row <- out %>% filter(Term == "Type")
  
  tibble(
    Site = site_letter,
    R2 = type_row$R2,
    F  = type_row$F,
    p  = type_row$`Pr(>F)`,
    n_samples = nrow(sub_meta),
    n_July = sum(sub_meta$Season == "July"),
    n_Sept = sum(sub_meta$Season == "September")
  )
}

sites <- sort(unique(meta_s$Site))
pw_tbl <- bind_rows(lapply(sites, pairwise_within_site_type)) %>%
  mutate(p_adj = p.adjust(p, method="BH"))

write_tsv(pw_tbl, "tables/pairwise_PERMANOVA_withinSite_Type_strataSeason.tsv")

txt <- c(
  "Pairwise PERMANOVA within Site: RG vs non-RG (Type)",
  "---------------------------------------------------",
  "Data: per-sample relative metal-category composition (Bray-Curtis).",
  "Model (within each Site): Type (RG vs non-RG).",
  "Permutations: 999, restricted within Season (strata = Season).",
  "",
  paste(capture.output(print(pw_tbl)), collapse="\n")
)
writeLines(txt, "tables/pairwise_PERMANOVA_withinSite_Type_strataSeason.txt")
cat("Saved pairwise PERMANOVA tables in tables/.\n")
##################################################################
##################################################################
##################################################################

library(vegan)

# Bray distance on samples
d <- vegdist(X, method="bray")

# centroide (media) per Base dentro ciascuna Season, poi distanza media tra centroidi
# (se vuoi replicare l'idea della figura che è per Season, questo è coerente)
centroids <- meta_s %>%
  select(Sample, Base, Site, Season) %>%
  group_by(Base, Site, Season) %>%
  summarise(Samples = list(Sample), .groups="drop")

centroid_vec <- function(samples) colMeans(X[samples, , drop=FALSE])
centroid_mat <- t(sapply(centroids$Samples, centroid_vec))
rownames(centroid_mat) <- paste(centroids$Base, centroids$Season, sep="__")

# distanze tra centroidi
dc <- as.matrix(vegdist(centroid_mat, method="bray"))

# costruisci lista di coppie centroidi (solo dentro la stessa Season)
pairs <- expand.grid(i = rownames(dc), j = rownames(dc), stringsAsFactors = FALSE) %>%
  filter(i < j) %>%
  mutate(
    Base_i = sub("__.*","", i),
    Base_j = sub("__.*","", j),
    Season_i = sub(".*__","", i),
    Season_j = sub(".*__","", j)
  ) %>%
  filter(Season_i == Season_j) %>%
  mutate(
    Site_i = substr(Base_i, 1, 1),
    Site_j = substr(Base_j, 1, 1),
    same_site = Site_i == Site_j,
    dist = dc[cbind(i, j)]
  )

# confronto: within-site (BF vs BR ecc.) vs between-site
w <- wilcox.test(dist ~ same_site, data = pairs, exact = FALSE)

write_tsv(pairs, "tables/centroidDistances_within_vs_between_site.tsv")
writeLines(
  c(
    "Centroid distance comparison (within-site vs between-site)",
    "----------------------------------------------------------",
    "Centroids computed per Base within each Season; distances are Bray-Curtis on metal-category composition.",
    sprintf("Wilcoxon test (within-site vs between-site): W=%.1f, p=%.4g", w$statistic, w$p.value)
  ),
  "tables/centroidDistances_within_vs_between_site.txt"
)
##################################################################
##################################################################
##################################################################
library(dplyr)
library(readr)

gp <- read_tsv("tables/metalGenes_rel_median_byBaseSeason_Top15.tsv", show_col_types = FALSE)

top1 <- gp %>%
  filter(Gene2 != "Other genes") %>%
  group_by(Base, Season) %>%
  slice_max(order_by = gene_rel, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(top_gene_pct = 100 * gene_rel) %>%
  arrange(Season, Base)

write_tsv(top1, "tables/top_gene_per_BaseSeason.tsv")
top1

dom <- top1 %>%
  summarise(
    median_top_gene_pct = median(top_gene_pct),
    max_top_gene_pct = max(top_gene_pct)
  )
dom

other <- gp %>%
  filter(Gene2 == "Other genes") %>%
  mutate(other_pct = 100 * gene_rel) %>%
  select(Base, Season, other_pct) %>%
  arrange(desc(other_pct))

write_tsv(other, "tables/other_genes_percent_by_BaseSeason.tsv")
other
##################################################################
###############Evenness + wtpC dominance##########################
##################################################################

library(dplyr)
library(readr)
library(tidyr)
library(stringr)

dir.create("tables", showWarnings = FALSE)

# ---- inputs ----
mat_file <- "tables/total_metal_mrg_per16s.tsv"  # genes x samples (metal-only)
map_in   <- "bacmet.txt"                         # has BacMet_ID, Gene_name

# ---- 1) read matrix + gene names ----
mat <- read_tsv(mat_file, show_col_types = FALSE)
stopifnot("BacMet_ID" %in% names(mat))

mat <- mat %>% mutate(BacMet_ID_clean = str_extract(BacMet_ID, "^BAC\\d+"))

bacmap <- read_tsv(map_in, show_col_types = FALSE) %>%
  transmute(
    BacMet_ID = as.character(BacMet_ID),
    Gene_name = as.character(Gene_name)
  ) %>%
  distinct()

# long format with gene labels
long <- mat %>%
  pivot_longer(cols = -c(BacMet_ID, BacMet_ID_clean),
               names_to = "Sample", values_to = "Abundance") %>%
  mutate(
    Abundance = as.numeric(Abundance),
    Abundance = ifelse(is.na(Abundance), 0, Abundance),
    Base = str_sub(Sample, 1, 2),
    Site = str_sub(Sample, 1, 1),
    Season_code = str_sub(Sample, 3, 3),
    Season = recode(Season_code, "J"="July", "S"="September", .default = Season_code),
    Replicate = str_extract(Sample, "\\d+$"),
    Type = ifelse(str_sub(Base, 2, 2) == "R", "RG", "non-RG")
  ) %>%
  left_join(bacmap, by = c("BacMet_ID_clean" = "BacMet_ID")) %>%
  mutate(Gene_label = ifelse(is.na(Gene_name) | Gene_name == "", BacMet_ID_clean, Gene_name))

# ---- 2) match wtpC (robust) ----
wtp_candidates <- long %>%
  distinct(Gene_label) %>%
  filter(str_detect(str_to_lower(Gene_label), "wtpc")) %>%
  pull(Gene_label)

if (length(wtp_candidates) == 0) {
  warning("No Gene_label matched 'wtpC'. Check spelling in Gene_label (bacmet Gene_name).")
}
if (length(wtp_candidates) > 1) {
  message("Multiple wtpC-like labels found. Using the first one. Candidates:\n",
          paste(wtp_candidates, collapse = ", "))
}
wtp_gene <- if (length(wtp_candidates) >= 1) wtp_candidates[1] else NA_character_

# ---- 3) replicate-median abundance per Base × Season × Gene (as in your Figure logic) ----
gene_base <- long %>%
  group_by(Base, Site, Type, Season, Gene_label) %>%
  summarise(gene_abund_med = median(Abundance, na.rm = TRUE), .groups = "drop")

# convert to within-group relative abundance (sum to 1 over ALL metal-associated genes)
gene_base <- gene_base %>%
  group_by(Base, Season) %>%
  mutate(gene_rel = ifelse(sum(gene_abund_med, na.rm=TRUE) > 0,
                           gene_abund_med / sum(gene_abund_med, na.rm=TRUE),
                           NA_real_)) %>%
  ungroup()

# ---- POINT (1): evenness metrics (Shannon + Pielou) per Base × Season ----
evenness_tbl <- gene_base %>%
  group_by(Base, Season) %>%
  summarise(
    S_richness = sum(gene_abund_med > 0, na.rm=TRUE),
    Shannon_H  = {
      p <- gene_rel[gene_rel > 0 & !is.na(gene_rel)]
      -sum(p * log(p))
    },
    Pielou_J   = ifelse(S_richness > 1, Shannon_H / log(S_richness), NA_real_),
    .groups = "drop"
  )

write_tsv(evenness_tbl, "tables/Fig3b_gene_evenness_Shannon_Pielou.tsv")

# optional: quick seasonal summary of evenness (median across Bases)
evenness_season_summary <- evenness_tbl %>%
  group_by(Season) %>%
  summarise(
    median_Shannon = median(Shannon_H, na.rm=TRUE),
    median_Pielou  = median(Pielou_J,  na.rm=TRUE),
    .groups = "drop"
  )
write_tsv(evenness_season_summary, "tables/Fig3b_evenness_bySeason_summary.tsv")

# ---- POINT (2): wtpC dominance (%), per Base × Season + overall summary ----
wtp_tbl <- gene_base %>%
  filter(!is.na(wtp_gene)) %>%
  filter(Gene_label == wtp_gene) %>%
  transmute(
    Base, Season, Site, Type,
    wtpC_pct = 100 * gene_rel
  ) %>%
  arrange(Season, Base)

write_tsv(wtp_tbl, "tables/Fig3b_wtpC_percent_byBaseSeason.tsv")

wtp_summary <- wtp_tbl %>%
  summarise(
    gene = wtp_gene,
    median_pct = median(wtpC_pct, na.rm=TRUE),
    min_pct    = min(wtpC_pct, na.rm=TRUE),
    max_pct    = max(wtpC_pct, na.rm=TRUE)
  )

write_tsv(wtp_summary, "tables/Fig3b_wtpC_percent_summary.tsv")

# ---- ALSO useful: top gene per Base × Season and its % (to check if wtpC always #1) ----
top_gene_tbl <- gene_base %>%
  group_by(Base, Season) %>%
  slice_max(order_by = gene_rel, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(Base, Season, top_gene = Gene_label, top_gene_pct = 100 * gene_rel) %>%
  arrange(Season, Base)

write_tsv(top_gene_tbl, "tables/Fig3b_top_gene_byBaseSeason.tsv")

# ---- print key outputs ----
cat("\nSaved:\n",
    "- tables/Fig3b_gene_evenness_Shannon_Pielou.tsv\n",
    "- tables/Fig3b_wtpC_percent_summary.tsv\n",
    "- tables/Fig3b_wtpC_percent_byBaseSeason.tsv\n",
    "- tables/Fig3b_top_gene_byBaseSeason.tsv\n\n", sep="")

print(wtp_summary)

#### 9) Spearman: metal concentration vs matching metal-MRG abundance ---------

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(readr)
library(tibble)

#----------------------------
# A) Read geochemistry
#----------------------------
geochem <- read_tsv("metadata3.txt", show_col_types = FALSE)

# Se hai "Sample" (consigliato), ricavo Base e Season da Sample
if ("Sample" %in% names(geochem)) {
  geochem <- geochem %>%
    mutate(
      Base = str_sub(Sample, 1, 2),
      Season_code = str_sub(Sample, 3, 3),
      Season = recode(Season_code, "J"="July", "S"="September", .default = Season_code)
    )
} else {
  stop("Nel geochem manca la colonna Sample: serve per ricavare Season (3a lettera).")
}

# Metalli disponibili nel geochem (colonne) e nel mapping MRG
metals_in_geo <- intersect(unique(metal_comp_sample_all$Metal), names(geochem))

# geochem in long: Base×Season×Metal
geochem_long <- geochem %>%
  pivot_longer(cols = all_of(metals_in_geo), names_to = "Metal", values_to = "Conc") %>%
  mutate(
    Metal = as.character(Metal),
    Conc  = as.numeric(Conc)
  ) %>%
  group_by(Base, Season, Metal) %>%                     # aggrega replicati se ci sono
  summarise(Conc = median(Conc, na.rm = TRUE), .groups="drop")

#----------------------------
# B) Aggregate MRG by Base×Season×Metal (median across replicates)
#----------------------------
metal_mrg_base <- metal_comp_sample_all %>%
  group_by(Base, Season, Metal) %>%
  summarise(MRG_metal = median(Metal_abund, na.rm = TRUE), .groups="drop")

#----------------------------
# C) Join + Spearman per Metal (metallo misurato vs MRG dello stesso metallo)
#----------------------------
df_corr_in <- metal_mrg_base %>%
  inner_join(geochem_long, by = c("Base","Season","Metal")) %>%
  filter(is.finite(Conc), is.finite(MRG_metal))

corr_tbl <- df_corr_in %>%
  group_by(Metal) %>%
  summarise(
    n = sum(complete.cases(Conc, MRG_metal)),
    rho = ifelse(n >= 3,
                 suppressWarnings(cor(Conc, MRG_metal, method="spearman", use="complete.obs")),
                 NA_real_),
    p_value = ifelse(n >= 3,
                     suppressWarnings(cor.test(Conc, MRG_metal, method="spearman", exact=FALSE)$p.value),
                     NA_real_),
    .groups="drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method="BH")) %>%
  arrange(p_adj)

write_tsv(corr_tbl, "tables/spearman_metalConc_vs_metalMRG_matched.tsv")

#----------------------------
# D) Plot stile Fig2c (barre di rho + stelline FDR)
#----------------------------
corr_plot <- corr_tbl %>%
  mutate(sig = case_when(
    is.na(p_adj) ~ "",
    p_adj < 0.001 ~ "***",
    p_adj < 0.01  ~ "**",
    p_adj < 0.05  ~ "*",
    TRUE ~ ""
  ))

p_corr <- ggplot(corr_plot, aes(x = reorder(Metal, rho), y = rho)) +
  geom_col() +
  geom_text(aes(label = sig), vjust = -0.2, size = 3) +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Metal",
    y = "Spearman rho (metal concentration vs matched metal-MRG)",
    title = "Matched metal correlations"
  )

ggsave("figures/spearman_metalConc_vs_metalMRG_matched.png", p_corr,
       width = 7, height = 4.5, dpi = 300)

library(dplyr)
library(ggplot2)
library(forcats)

corr_plot <- corr_tbl %>%
  filter(!is.na(rho)) %>%
  mutate(
    Metal = fct_reorder(Metal, rho),
    direction = ifelse(rho >= 0, "Positive", "Negative"),
    sig = case_when(
      is.na(p_adj) ~ "",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    ),
    star_hjust = ifelse(rho >= 0, -0.2, 1.2)
  )

p_corr_pretty2 <- ggplot(corr_plot, aes(x = rho, y = Metal, fill = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_col(width = 0.75, color = "white", linewidth = 0.2) +
  geom_text(aes(label = sig), hjust = corr_plot$star_hjust, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  scale_fill_manual(values = c(Positive = "#4E79A7", Negative = "#E15759"), guide = "none") +
  labs(
    title = "Matched metal correlations",
    x = "Spearman rho (metal concentration vs matched metal-MRG)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 11)
  )

ggsave("figures/spearman_metalConc_vs_metalMRG_matched_pretty_color.png",
       p_corr_pretty2, width = 8, height = 4.8, dpi = 300)


#############################################################
############################################################
library(dplyr)
library(ggplot2)
library(forcats)

corr_plot <- corr_tbl %>%
  filter(!is.na(rho)) %>%
  mutate(
    Metal = fct_reorder(Metal, rho),
    sig = case_when(
      is.na(p_adj) ~ "",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    ),
    # testo con valore rho
    rho_lab = sprintf("%.2f", rho),
    # posizione del numero (sopra se positivo, sotto se negativo)
    vjust_rho = ifelse(rho >= 0, -0.25, 1.25),
    # stelline un po' più distanti
    vjust_sig = ifelse(rho >= 0, -1.05, 2.05)
  )

max_abs <- max(abs(corr_plot$rho), na.rm = TRUE)
pad <- 0.22 * max_abs
y_lim <- c(-(max_abs + pad), (max_abs + pad))

p_vert <- ggplot(corr_plot, aes(x = Metal, y = rho)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey40") +
  geom_col(width = 0.85, fill = "grey35") +
  # numeri rho
  geom_text(aes(label = rho_lab), vjust = corr_plot$vjust_rho, size = 3.6, color = "black") +
  # stelline
  geom_text(aes(label = sig), vjust = corr_plot$vjust_sig, size = 4, color = "black") +
  coord_cartesian(ylim = y_lim) +
  labs(
    title = "Matched metal correlations",
    x = NULL,
    y = "Spearman rho"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

ggsave("figures/spearman_matched_vertical_WHITE_numbers.png",
       p_vert, width = 10, height = 6, dpi = 300, bg = "white")

## =========================
############################################################
## HEATMAP "like example": metal-MRG classes (merged reps, log10)
## - rows: metal categories (Cu, Zn, ...)
## - cols: Group (Base+Season = BFJ, BFS, ...)
## - reps merged: median across replicates
## - values: log10(relative abundance + 1e-6)  -> range ~ -6..0 like your example
## - top annotation: Base colours
############################################################
pkgs <- c("dplyr","readr","tidyr","stringr","tibble","pheatmap","RColorBrewer")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(to_install)) install.packages(to_install, repos="https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)

# ---- INPUTS ----
mat_file <- "tables/total_metal_mrg_per16s.tsv"  # metal-only genes x samples (per16S)
map_in   <- "bacmet.txt"                         # has BacMet_ID + Compound
stopifnot(file.exists(mat_file), file.exists(map_in))

# ---- SETTINGS ----
TOP_N_METALS <- 18        # cambia (es 12/15/20)
PSEUDO_REL   <- 1e-6      # per log10(relative + pseudo)
OUT_PREFIX   <- "HEATMAP_metalMRG_metalClasses_topN"

# ---- Base palette (your usual) ----
base_cols_map <- c(
  BF = "#1b9e77", BR = "#a6dba0",
  CF = "#d95f02", CR = "#fdb863",
  PF = "#7570b3", PR = "#b2abd2",
  SF = "#1f9ac2", SR = "#a6dce7",
  VF = "#e7298a", VR = "#f2b2d4"
)

# ---- metal dictionary (same as your "ALL") ----
metal_terms_focus <- tibble::tribble(
  ~Metal, ~regex,
  "Cu", "copper",
  "Zn", "zinc",
  "Ni", "nickel",
  "Co", "cobalt",
  "Cd", "cadmium",
  "As", "arsenic|arsenate|arsenite",
  "Cr", "chromium|chromate",
  "Mn", "manganese",
  "Fe", "iron",
  "Pb", "lead",
  "Hg", "mercury|organomercury|organo-mercury",
  "Ag", "silver",
  "Sn", "tin|organotin|organo-tin",
  "Sb", "antimony",
  "Se", "selenium|selenate|selenite",
  "Mo", "molybdenum|molybdate",
  "W",  "tungsten|wolfram",
  "Au", "gold",
  "Al", "aluminium|aluminum",
  "V",  "vanadium|vanadate",
  "U",  "uranium|uranyl"
)

# ---- 1) read metal-only matrix and go long ----
mat_wide <- read_tsv(mat_file, show_col_types = FALSE) %>%
  mutate(BacMet_ID_clean = str_extract(BacMet_ID, "^BAC\\d+"))

sample_cols <- setdiff(names(mat_wide), c("BacMet_ID","BacMet_ID_clean"))

metal_comp_sample <- mrg_long %>%
  inner_join(gene_metal_long,
             by = c("BacMet_ID_clean" = "BacMet_ID"),
             relationship = "many-to-many") %>%
  group_by(Sample, Group, Base, Season, Replicate, Metal) %>%
  summarise(Metal_abund = sum(Abundance, na.rm = TRUE), .groups="drop")

# ---- 2) gene -> metal mapping from BacMet Compound ----
bacmap <- read_tsv(map_in, show_col_types = FALSE) %>%
  transmute(
    BacMet_ID = as.character(BacMet_ID),
    comp_low  = str_to_lower(ifelse(is.na(Compound), "", as.character(Compound)))
  )

gene_metal_long <- bacmap %>%
  tidyr::crossing(metal_terms_focus) %>%
  mutate(hit = str_detect(comp_low, regex)) %>%
  filter(hit) %>%
  select(BacMet_ID, Metal) %>%
  distinct()

# ---- 3) per-sample metal category abundance (sum of genes that match that metal) ----
metal_comp_sample <- mrg_long %>%
  inner_join(gene_metal_long, by = c("BacMet_ID_clean" = "BacMet_ID")) %>%
  group_by(Sample, Group, Base, Season, Replicate, Metal) %>%
  summarise(Metal_abund = sum(Abundance, na.rm = TRUE), .groups="drop")

# ---- 4) merge reps: median per Group x Metal ----
metal_group <- metal_comp_sample %>%
  group_by(Group, Base, Season, Metal) %>%
  summarise(Metal_abund_med = median(Metal_abund, na.rm = TRUE), .groups="drop")

# ---- 5) build matrix (Metal x Group), then relative within Group + log10 ----
mat <- metal_group %>%
  select(Metal, Group, Metal_abund_med) %>%
  pivot_wider(names_from = Group, values_from = Metal_abund_med, values_fill = 0)

X <- as.matrix(mat[,-1])
rownames(X) <- mat$Metal

# relative within column (like composition), then log10 (gives negative range)
X_rel <- sweep(X, 2, pmax(colSums(X), 1), "/")
X_log <- log10(X_rel + PSEUDO_REL)

# ---- 6) keep Top N metals by mean relative abundance ----
metal_rank <- tibble(
  Metal = rownames(X_rel),
  mean_rel = rowMeans(X_rel, na.rm = TRUE)
) %>% arrange(desc(mean_rel))

n_keep <- min(TOP_N_METALS, nrow(metal_rank))
keep_metals <- metal_rank %>% slice_head(n = n_keep) %>% pull(Metal)

X_log_top <- X_log[keep_metals, , drop = FALSE]
X_rel_top <- X_rel[keep_metals, , drop = FALSE]

# save matrices used
write_tsv(as_tibble(X_rel_top, rownames="Metal"),
          sprintf("tables/%s_matrix_relative.tsv", OUT_PREFIX))
write_tsv(as_tibble(X_log_top, rownames="Metal"),
          sprintf("tables/%s_matrix_log10_relative.tsv", OUT_PREFIX))
write_tsv(metal_rank,
          sprintf("tables/%s_metals_rank.tsv", OUT_PREFIX))

# ---- 7) column annotation = Base (colours) ----
groups <- colnames(X_log_top)
base_lab <- substr(groups, 1, 2)

ann_col <- data.frame(Base = factor(base_lab, levels = sort(unique(base_lab))))
rownames(ann_col) <- groups

missing <- setdiff(levels(ann_col$Base), names(base_cols_map))
if (length(missing) > 0) {
  extra <- RColorBrewer::brewer.pal(max(3, min(8, length(missing))), "Set2")
  base_cols_map[missing] <- extra[seq_along(missing)]
}
ann_colors <- list(Base = base_cols_map[levels(ann_col$Base)])

# ---- 8) colors like your example (white -> orange -> red) ----
vals <- as.vector(X_log_top); vals <- vals[is.finite(vals)]
bk <- seq(min(vals), max(vals), length.out = 101)
cols <- colorRampPalette(c("white", "#FDB863", "#E34A33", "#B30000"))(100)

# sizes (avoid squished)
w <- max(10, 0.35 * ncol(X_log_top) + 3)
h <- max(7,  0.35 * nrow(X_log_top) + 3)

# ---- 9) plot ----
pheatmap::pheatmap(
  X_log_top,
  color = cols, breaks = bk,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  border_color = "grey80",
  fontsize_row = 10,
  fontsize_col = 10,
  treeheight_col = 45,
  treeheight_row = 45,
  main = "Top metal classes (merged reps, log10 relative abundance)",
  filename = sprintf("figures/%s.png", OUT_PREFIX),
  width = w, height = h
)

pheatmap::pheatmap(
  X_log_top,
  color = cols, breaks = bk,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  border_color = "grey80",
  fontsize_row = 10,
  fontsize_col = 10,
  treeheight_col = 45,
  treeheight_row = 45,
  main = "Top metal classes (merged reps, log10 relative abundance)",
  filename = sprintf("figures/%s.pdf", OUT_PREFIX),
  width = w, height = h
)

cat("Saved:\n",
    "- figures/", OUT_PREFIX, ".png\n",
    "- figures/", OUT_PREFIX, ".pdf\n",
    "- tables/", OUT_PREFIX, "_matrix_relative.tsv\n",
    "- tables/", OUT_PREFIX, "_matrix_log10_relative.tsv\n", sep="")


library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(vegan)
library(tibble)

dir.create("tables", showWarnings = FALSE)

# -----------------------
# INPUT
# -----------------------
mat_file  <- "tables/total_metal_mrg_per16s.tsv"  # metal-only genes x samples (per16S)
meta_file <- "metadata3.txt"                      # contiene Sample + Valley (e magari altro)

# scegli: PERMANOVA su abbondanze (Bray) oppure su composizione (relative)
USE_RELATIVE <- TRUE   # TRUE consigliato se vuoi confrontare composizione; FALSE = load assoluto

# -----------------------
# helper: read matrix robust
# -----------------------
read_matrix_tsv_robust <- function(path) {
  x <- read.delim(
    file = path, header = TRUE, sep = "\t",
    check.names = FALSE, fill = TRUE,
    quote = "", comment.char = "",
    stringsAsFactors = FALSE
  )
  idcol <- colnames(x)[1]
  ids <- x[[idcol]]
  keep <- !(is.na(ids) | trimws(ids) == "")
  x <- x[keep, , drop = FALSE]
  ids <- make.unique(as.character(x[[idcol]]))
  rownames(x) <- ids
  x[[idcol]] <- NULL
  x[] <- lapply(x, function(v) suppressWarnings(as.numeric(v)))
  m <- as.matrix(x); m[is.na(m)] <- 0
  m
}

# -----------------------
# 1) Read matrix + metadata
# -----------------------
m <- read_matrix_tsv_robust(mat_file)   # genes x samples
meta <- read_tsv(meta_file, show_col_types = FALSE)

# trova colonna Sample
if (!("Sample" %in% names(meta))) {
  cand <- names(meta)[tolower(names(meta)) %in% c("sample","id")]
  if (length(cand) == 0) stop("metadata3.txt: non trovo una colonna Sample (o sample/id).")
  meta <- meta %>% rename(Sample = all_of(cand[1]))
}
meta <- meta %>% mutate(Sample = as.character(Sample))

# check Valley
if (!("Valley" %in% names(meta))) stop("metadata3.txt: manca la colonna 'Valley'.")

# -----------------------
# 2) Allinea campioni
# -----------------------
common <- intersect(colnames(m), meta$Sample)
if (length(common) < 4) stop("Pochi campioni in comune tra matrice e metadata3.txt. Controlla i nomi.")

m <- m[, common, drop = FALSE]
meta2 <- meta %>% filter(Sample %in% common)

# ricostruisci Type/Season se non presenti (dal nome Sample)
if (!("Base" %in% names(meta2)))   meta2$Base   <- substr(meta2$Sample, 1, 2)
if (!("Season" %in% names(meta2))) meta2$Season <- substr(meta2$Sample, 3, 3)
if (!("Type" %in% names(meta2)))   meta2$Type   <- ifelse(substr(meta2$Base, 2, 2) == "R", "RG", "non-RG")

meta2 <- meta2 %>%
  mutate(
    Valley = factor(Valley),
    Type  = factor(Type),
    Season = factor(Season)   # J/S oppure July/September se già così
  )

# riordina meta come colonne matrice
meta2 <- meta2 %>% distinct(Sample, .keep_all = TRUE)
meta2 <- meta2[match(common, meta2$Sample), ]
stopifnot(all(meta2$Sample == common))

# -----------------------
# 3) Matrice samples x genes + standardizzazione
# -----------------------
X <- t(m)  # samples x genes

# opzionale: converti in composizione per campione
if (USE_RELATIVE) {
  rs <- rowSums(X)
  X <- sweep(X, 1, ifelse(rs > 0, rs, 1), "/")
}

# distanza Bray-Curtis
d <- vegdist(X, method = "bray")

# -----------------------
# 4) PERMANOVA
# -----------------------
set.seed(1)

# (A) effetti marginali (consigliato per interpretare ogni fattore “a parità degli altri”)
perm_marg <- adonis2(d ~ valley + Type + Season,
                     data = meta2,
                     permutations = 999,
                     by = "margin")

# (B) sequenziale (ordine conta): Valley poi Type poi Season
perm_seq <- adonis2(d ~ valley + Type + Season,
                    data = meta2,
                    permutations = 999,
                    by = "terms")

write_tsv(as.data.frame(perm_marg) %>% rownames_to_column("Term"),
          "tables/PERMANOVA_metalMRG_Valley_Type_Season_MARGIN.tsv")
write_tsv(as.data.frame(perm_seq) %>% rownames_to_column("Term"),
          "tables/PERMANOVA_metalMRG_Valley_Type_Season_TERMS.tsv")

# -----------------------
# 5) PERMDISP (dispersione) per i fattori
# -----------------------
bd_valley <- betadisper(d, meta2$valley); pd_valley <- anova(bd_valley)
bd_type   <- betadisper(d, meta2$Type);   pd_type   <- anova(bd_type)
bd_season <- betadisper(d, meta2$Season); pd_season <- anova(bd_season)

disp_tbl <- tibble(
  Factor = c("Valley","Type","Season"),
  F = c(pd_valley$`F value`[1], pd_type$`F value`[1], pd_season$`F value`[1]),
  p = c(pd_valley$`Pr(>F)`[1],  pd_type$`Pr(>F)`[1],  pd_season$`Pr(>F)`[1])
) %>% mutate(p_adj = p.adjust(p, method = "BH"))

write_tsv(disp_tbl, "tables/PERMDISP_metalMRG_Valley_Type_Season.tsv")

cat("\nSaved:\n",
    "- tables/PERMANOVA_metalMRG_Valley_Type_Season_MARGIN.tsv\n",
    "- tables/PERMANOVA_metalMRG_Valley_Type_Season_TERMS.tsv\n",
    "- tables/PERMDISP_metalMRG_Valley_Type_Season.tsv\n", sep = "")
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(vegan)
library(tibble)

dir.create("tables", showWarnings = FALSE)

# -----------------------
# INPUT
# -----------------------
mat_file  <- "tables/total_metal_mrg_per16s.tsv"  # metal-only genes x samples (per16S)
meta_file <- "metadata3.txt"                      # deve contenere Sample + Valley

USE_RELATIVE <- TRUE   # TRUE = composizione (consigliato); FALSE = load assoluto (per16S)

# -----------------------
# helper: robust read matrix
# -----------------------
read_matrix_tsv_robust <- function(path) {
  x <- read.delim(
    file = path, header = TRUE, sep = "\t",
    check.names = FALSE, fill = TRUE,
    quote = "", comment.char = "",
    stringsAsFactors = FALSE
  )
  idcol <- colnames(x)[1]
  ids <- x[[idcol]]
  keep <- !(is.na(ids) | trimws(ids) == "")
  x <- x[keep, , drop = FALSE]
  ids <- make.unique(as.character(x[[idcol]]))
  rownames(x) <- ids
  x[[idcol]] <- NULL
  x[] <- lapply(x, function(v) suppressWarnings(as.numeric(v)))
  m <- as.matrix(x); m[is.na(m)] <- 0
  m
}

# -----------------------
# 1) Read
# -----------------------
m <- read_matrix_tsv_robust(mat_file)   # genes x samples
meta <- read_tsv(meta_file, show_col_types = FALSE)

# find Sample column
if (!("Sample" %in% names(meta))) {
  cand <- names(meta)[tolower(names(meta)) %in% c("sample","id")]
  if (length(cand) == 0) stop("metadata3.txt: non trovo colonna Sample (o sample/id).")
  meta <- meta %>% rename(Sample = all_of(cand[1]))
}
meta <- meta %>% mutate(Sample = as.character(Sample))

if (!("valley" %in% names(meta))) stop("metadata3.txt: manca la colonna 'Valley'.")

# ensure core vars (fallback from Sample name)
if (!("Base" %in% names(meta)))   meta$Base   <- substr(meta$Sample, 1, 2)
if (!("Season" %in% names(meta))) meta$Season <- substr(meta$Sample, 3, 3)
if (!("Type" %in% names(meta)))   meta$Type   <- ifelse(substr(meta$Base, 2, 2) == "R", "RG", "non-RG")

# Group id for replicate-merge: Base + Season (BFJ, BFS, ...)
meta <- meta %>% mutate(Group = paste0(Base, Season))

# -----------------------
# 2) Align samples
# -----------------------
common <- intersect(colnames(m), meta$Sample)
if (length(common) < 4) stop("Pochi campioni in comune tra matrice e metadata. Controlla i nomi.")

m <- m[, common, drop = FALSE]
meta2 <- meta %>% filter(Sample %in% common)

# -----------------------
# 3) Merge replicates: median per Group (Base×Season) for each gene
# -----------------------
mat_long <- as.data.frame(m) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Abundance") %>%
  left_join(meta2 %>% select(Sample, Group), by = "Sample") %>%
  group_by(Gene, Group) %>%
  summarise(Abundance = median(Abundance, na.rm = TRUE), .groups = "drop")

m_group <- mat_long %>%
  pivot_wider(names_from = Group, values_from = Abundance, values_fill = 0)

# genes x groups matrix
G <- as.matrix(m_group[,-1])
rownames(G) <- m_group$Gene

# samples(groups) x genes
X <- t(G)

# optional: relative within group (composition)
if (USE_RELATIVE) {
  rs <- rowSums(X)
  X <- sweep(X, 1, ifelse(rs > 0, rs, 1), "/")
}

# -----------------------
# 4) Build group-level metadata + consistency check
# -----------------------
meta_g <- meta2 %>%
  group_by(Group) %>%
  summarise(
    n_reps = n(),
    valley_n = n_distinct(valley),
    Type_n  = n_distinct(Type),
    Season_n = n_distinct(Season),
    Base_n  = n_distinct(Base),
    valley = first(valley),
    Type   = first(Type),
    Season = first(Season),
    Base   = first(Base),
    .groups = "drop"
  )

# flag inconsistencies (should be 1)
incons <- meta_g %>% filter(Valley_n > 1 | Type_n > 1 | Season_n > 1 | Base_n > 1)
if (nrow(incons) > 0) {
  write_tsv(incons, "tables/WARN_inconsistent_metadata_within_Group.tsv")
  warning("Trovate inconsistenze in metadata dentro alcuni Group. Vedi tables/WARN_inconsistent_metadata_within_Group.tsv")
}

# align meta_g to X rows
meta_g <- meta_g %>% filter(Group %in% rownames(X))
meta_g <- meta_g[match(rownames(X), meta_g$Group), ]
stopifnot(all(meta_g$Group == rownames(X)))

meta_g <- meta_g %>%
  mutate(
    valley = factor(valley),
    Type   = factor(Type),
    Season = factor(Season)
  )

# -----------------------
# 5) PERMANOVA + PERMDISP
# -----------------------
d <- vegdist(X, method = "bray")
set.seed(1)

perm_marg <- adonis2(d ~ valley + Type + Season,
                     data = meta_g,
                     permutations = 999,
                     by = "margin")

perm_terms <- adonis2(d ~ valley + Type + Season,
                      data = meta_g,
                      permutations = 999,
                      by = "terms")

write_tsv(as.data.frame(perm_marg) %>% rownames_to_column("Term"),
          "tables/PERMANOVA_metalMRG_repMerged_Valley_Type_Season_MARGIN.tsv")
write_tsv(as.data.frame(perm_terms) %>% rownames_to_column("Term"),
          "tables/PERMANOVA_metalMRG_repMerged_Valley_Type_Season_TERMS.tsv")

# dispersion
bd_valley <- betadisper(d, meta_g$valley); pd_valley <- anova(bd_valley)
bd_type   <- betadisper(d, meta_g$Type);   pd_type   <- anova(bd_type)
bd_season <- betadisper(d, meta_g$Season); pd_season <- anova(bd_season)

disp_tbl <- tibble(
  Factor = c("valley","Type","Season"),
  F = c(pd_valley$`F value`[1], pd_type$`F value`[1], pd_season$`F value`[1]),
  p = c(pd_valley$`Pr(>F)`[1],  pd_type$`Pr(>F)`[1],  pd_season$`Pr(>F)`[1])
) %>% mutate(p_adj = p.adjust(p, method = "BH"))

write_tsv(disp_tbl, "tables/PERMDISP_metalMRG_repMerged_Valley_Type_Season.tsv")

# save group metadata used
write_tsv(meta_g, "tables/metalMRG_repMerged_group_metadata.tsv")

cat("\nSaved:\n",
    "- tables/PERMANOVA_metalMRG_repMerged_Valley_Type_Season_MARGIN.tsv\n",
    "- tables/PERMANOVA_metalMRG_repMerged_Valley_Type_Season_TERMS.tsv\n",
    "- tables/PERMDISP_metalMRG_repMerged_Valley_Type_Season.tsv\n",
    "- tables/metalMRG_repMerged_group_metadata.tsv\n", sep="")

