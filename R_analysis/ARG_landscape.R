###############################################################################
# ARG per 16S — replicate-united plots (like your MRG workflow)
# - Reads gene x sample matrix (ARG per 16S)
# - Extracts Base/Site/Season/Type/Replicate from sample names
# - OPTIONAL: set low bacterial Mb replicates to NA before aggregating
# - Aggregates replicates by median (Base × Season)
# - Computes Total & Richness, and produces plots + tables
###############################################################################

#### 0) Packages, folders -----------------------------------------------------
pkgs <- c("dplyr","readr","tidyr","stringr","ggplot2","tibble")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(to_install)) install.packages(to_install, repos="https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

dir.create("tables", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

#### 1) Read ARG per16S (gene x sample) --------------------------------------
# Replace with your file
arg_file <- "arg_per16s_T1T2_renamed.tsv"
arg_wide <- read_tsv(arg_file, show_col_types = FALSE)

# first column assumed = GeneID (change if needed)
gene_col <- names(arg_wide)[1]
sample_cols <- setdiff(names(arg_wide), gene_col)

arg_long <- arg_wide %>%
  pivot_longer(cols = all_of(sample_cols),
               names_to = "Sample", values_to = "Abundance") %>%
  mutate(
    Abundance = as.numeric(Abundance),
    Abundance = ifelse(is.na(Abundance), 0, Abundance),
    GeneID = .data[[gene_col]]
  ) %>%
  select(GeneID, Sample, Abundance)

#### 2) Metadata from Sample name --------------------------------------------
meta <- arg_long %>%
  distinct(Sample) %>%
  mutate(
    Base   = str_sub(Sample, 1, 2),
    Site   = str_sub(Sample, 1, 1),
    Season_code = str_sub(Sample, 3, 3),
    Season = recode(Season_code, "J"="July", "S"="September", .default = Season_code),
    Replicate = str_extract(Sample, "\\d+$"),
    Type = ifelse(str_sub(Base, 2, 2) == "R", "RG", "non-RG")
  )

arg_long <- arg_long %>% left_join(meta, by = "Sample")

#### 3) OPTIONAL: downweight/remove low bacterial Mb replicates ---------------
# If you have bacterial Mb per sample in a metadata file, load & join it here.
# Example metadata file with columns: Sample, bacterial_Mb
# meta_cov <- read_tsv("metadata_bacterialMb.tsv", show_col_types = FALSE)
# arg_long <- arg_long %>% left_join(meta_cov, by = "Sample")

# Then set Abundance to NA for low-Mb replicates (so median ignores them)
# Choose a rule that is not ad hoc: e.g., 5th percentile of bacterial_Mb
# if("bacterial_Mb" %in% names(arg_long)){
#   thr <- quantile(unique(arg_long$bacterial_Mb), probs = 0.05, na.rm = TRUE)
#   arg_long <- arg_long %>%
#     mutate(Abundance = ifelse(bacterial_Mb < thr, NA_real_, Abundance))
#   message("Low bacterial Mb threshold used: ", signif(thr, 3))
# }

#### 4) Aggregate replicates by median (Base × Season × Gene) -----------------
# This produces a replicate-united matrix for plotting composition etc.
arg_base_gene <- arg_long %>%
  group_by(Base, Site, Type, Season, GeneID) %>%
  summarise(Abund_med = median(Abundance, na.rm = TRUE), .groups = "drop")

# If ALL replicates were NA for a group, median(,na.rm=TRUE) becomes NA.
# Convert those to 0 (interpreted as undetected after QC)
arg_base_gene <- arg_base_gene %>%
  mutate(Abund_med = ifelse(is.na(Abund_med), 0, Abund_med))

write_tsv(arg_base_gene, "tables/ARG_byBaseSeason_gene_median.tsv")

#### 5) Alpha metrics (Total + Richness) on replicate-united data -------------
alpha_base <- arg_base_gene %>%
  group_by(Base, Site, Type, Season) %>%
  summarise(
    Total_ARG_per16S = sum(Abund_med, na.rm = TRUE),
    ARG_richness     = sum(Abund_med > 0, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(alpha_base, "tables/alpha_ARG_total_richness_byBaseSeason.tsv")

#### 6) Plot: totals and richness (like your MRG) -----------------------------
p_total <- ggplot(alpha_base, aes(x = Base, y = Total_ARG_per16S)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.12, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~Season, nrow = 1) +
  labs(x = "Base", y = "Total ARG (per 16S; median across replicates)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave("figures/ARG_total_per16S_byBaseSeason.png", p_total, width = 11, height = 4.5, dpi = 300)

p_rich <- ggplot(alpha_base, aes(x = Base, y = ARG_richness)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.12, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~Season, nrow = 1) +
  labs(x = "Base", y = "ARG richness (genes; median across replicates)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave("figures/ARG_richness_byBaseSeason.png", p_rich, width = 11, height = 4.5, dpi = 300)

#### 7) Stacked bar: composition by Top N GENES -------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)

# Assumo: arg_long ha colonne GeneID, Abundance, Base, Season, (Site/Type ok ma non obbligatori)
arg_long <- arg_long %>% mutate(GeneID = as.character(GeneID))

# 1) Unisci repliche (mediana) -> abbondanza per Base×Season×Gene
arg_base_gene <- arg_long %>%
  group_by(Base, Season, GeneID) %>%
  summarise(Abund_med = median(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(Abund_med = ifelse(is.na(Abund_med), 0, Abund_med))

# 2) Calcola TopN geni usando la "mean Abund_med" (o prevalenza, se preferisci)
N <- 20
topN <- arg_base_gene %>%
  group_by(GeneID) %>%
  summarise(mean_abund = mean(Abund_med, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = N) %>%
  pull(GeneID)

# 3) Collassa in TopN vs Other (SENZA perdere gli altri)
arg_plot <- arg_base_gene %>%
  mutate(Gene2 = ifelse(GeneID %in% topN, GeneID, "Other genes")) %>%
  group_by(Base, Season, Gene2) %>%
  summarise(Abund = sum(Abund_med, na.rm = TRUE), .groups = "drop")

# 4) Ordine livelli: Other in cima (ultima voce)
levs <- c(sort(setdiff(unique(arg_plot$Gene2), "Other genes")), "Other genes")
arg_plot <- arg_plot %>% mutate(Gene2 = factor(Gene2, levels = levs))



# 0) Se GeneID contiene "id|name", creiamo una label solo "name"
arg_base_gene2 <- arg_base_gene %>%
  mutate(
    GeneID = as.character(GeneID),
    GeneName = ifelse(grepl("\\|", GeneID),
                      sub("^.*\\|", "", GeneID),   # tutto dopo |
                      GeneID)                       # fallback se non c'è |
  )

# 1) Escludi voci che non vuoi (case-insensitive)
drop_terms <- c("streptomyces", "pseudomonas", "acinetobacter")
arg_base_gene2 <- arg_base_gene2 %>%
  filter(!tolower(GeneName) %in% drop_terms)

# 2) Seleziona TopN su questi dati filtrati
N <- 15
topN_names <- arg_base_gene2 %>%
  group_by(GeneName) %>%   # NB: scegliamo per nome gene, non per ID numerico
  summarise(mean_abund = mean(Abund_med, na.rm = TRUE), .groups="drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = N) %>%
  pull(GeneName)

# 3) Collassa in TopN vs Other
arg_plot2 <- arg_base_gene2 %>%
  mutate(Gene2 = ifelse(GeneName %in% topN_names, GeneName, "Other genes")) %>%
  group_by(Base, Season, Gene2) %>%
  summarise(Abund = sum(Abund_med, na.rm = TRUE), .groups="drop")

# livelli CORRETTI per arg_plot2
levs2 <- c(sort(setdiff(unique(arg_plot2$Gene2), "Other genes")), "Other genes")
arg_plot2 <- arg_plot2 %>% mutate(Gene2 = factor(Gene2, levels = levs2))

k2 <- length(levs2)

# palette distinta (no repeat) + "pastellizza"
pal_raw  <- randomcoloR::distinctColorPalette(k2)
pal_soft <- colorspace::lighten(pal_raw, amount = 0.25)

pal2 <- setNames(pal_soft, levs2)
pal2["Other genes"] <- "grey70"

p <- ggplot(arg_plot2, aes(x = Base, y = Abund, fill = Gene2)) +
  geom_col(position = "fill", color="white", linewidth=0.2) +
  facet_wrap(~Season, nrow = 1) +
  scale_y_continuous(labels = function(x) paste0(round(100*x), "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = pal2) +
  labs(x="Site", y="ARG composition (median across replicates)", fill="ARG genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("figures/ARG_TopGenes_Top15byBaseSeason_filteredNames.png", p, width=12, height=5, dpi=600)

ggsave("figures/ARG_TopGenes_Top15_byBaseSeason_filteredNames.png", p, width=12, height=5, dpi=300)
write_tsv(arg_plot2, "tables/ARG_gene_composition_byBaseSeason_Top20_filteredNames.tsv")
write_tsv(tibble(Gene=topN_names), "tables/ARG_top_genes_list_filteredNames.tsv")
###############################################################################
# ARG season report (replicates united) + top genes %
###############################################################################
pkgs <- c("dplyr","tidyr","stringr","readr","vegan","effsize")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(to_install)) install.packages(to_install, repos="https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

dir.create("tables", showWarnings = FALSE)

# --------------------------
# INPUT expected:
# arg_base_gene2: Base, Season, GeneName, Abund_med
# --------------------------

stopifnot(all(c("Base","Season","GeneName","Abund_med") %in% names(arg_base_gene2)))

# 1) Build Base×Season x Gene matrix -----------------------------------------
mat_long <- arg_base_gene2 %>%
  group_by(Base, Season, GeneName) %>%
  summarise(Abund = sum(Abund_med, na.rm=TRUE), .groups="drop")

mat_wide <- mat_long %>%
  unite(Group, Base, Season, sep="_") %>%
  pivot_wider(names_from = GeneName, values_from = Abund, values_fill = 0)

meta_g <- mat_wide %>%
  select(Group) %>%
  separate(Group, into=c("Base","Season"), sep="_", remove=FALSE)

X <- mat_wide %>% select(-Group) %>% as.data.frame()

# helper: relative within group
rel_long <- mat_long %>%
  group_by(Base, Season) %>%
  mutate(Total = sum(Abund, na.rm=TRUE),
         Rel   = ifelse(Total > 0, Abund/Total, 0)) %>%
  ungroup()

# 2) Alpha diversity on Base×Season ------------------------------------------
alpha_base <- mat_long %>%
  group_by(Base, Season) %>%
  summarise(
    Total_ARG_per16S = sum(Abund, na.rm=TRUE),
    Richness         = sum(Abund > 0, na.rm=TRUE),
    .groups="drop"
  )

# Shannon from relative (Pielou evenness = Shannon/log(Richness))
alpha_shannon <- rel_long %>%
  group_by(Base, Season) %>%
  summarise(
    Shannon = vegan::diversity(Rel, index="shannon"),
    .groups="drop"
  )

alpha_base <- alpha_base %>%
  left_join(alpha_shannon, by=c("Base","Season")) %>%
  mutate(
    Evenness_Pielou = ifelse(Richness > 1, Shannon / log(Richness), NA_real_)
  )

write_tsv(alpha_base, "tables/ARG_alpha_byBaseSeason.tsv")

# 3) Season tests (July vs September) ----------------------------------------
wilcox_cliff <- function(df, response){
  d <- df %>% filter(!is.na(.data[[response]]), !is.na(Season))
  xJ <- d %>% filter(Season=="July") %>% pull(.data[[response]])
  xS <- d %>% filter(Season=="September") %>% pull(.data[[response]])
  wt <- wilcox.test(xJ, xS, exact = FALSE)
  cd <- effsize::cliff.delta(xJ, xS)$estimate
  tibble::tibble(
    metric = response,
    n_July = length(xJ),
    n_Sept = length(xS),
    median_July = median(xJ, na.rm=TRUE),
    median_Sept = median(xS, na.rm=TRUE),
    p_value = wt$p.value,
    cliffs_delta = as.numeric(cd)
  )
}

season_stats <- dplyr::bind_rows(
  wilcox_cliff(alpha_base, "Total_ARG_per16S"),
  wilcox_cliff(alpha_base, "Richness"),
  wilcox_cliff(alpha_base, "Shannon"),
  wilcox_cliff(alpha_base, "Evenness_Pielou")
) %>%
  mutate(p_adj = p.adjust(p_value, method="BH"))

write_tsv(season_stats, "tables/ARG_season_tests.tsv")

# 4) Beta diversity: Bray-Curtis + PERMANOVA + PERMDISP -----------------------
bc <- vegan::vegdist(X, method="bray")

perm <- vegan::adonis2(bc ~ Season + Base, data = meta_g, permutations = 999, by="margin")

# dispersion by Season (important check for PERMANOVA)
bd_season <- vegan::betadisper(bc, group = meta_g$Season)
permdisp_aov  <- anova(bd_season)
permdisp_perm <- permutest(bd_season, permutations = 999)

# 5) Top genes: prevalence + mean% + max% ------------------------------------
gene_summary <- rel_long %>%
  group_by(GeneName) %>%
  summarise(
    prevalence_groups = sum(Rel > 0),
    prevalence_pct    = 100*mean(Rel > 0),
    mean_pct          = 100*mean(Rel, na.rm=TRUE),
    max_pct           = 100*max(Rel, na.rm=TRUE),
    .groups="drop"
  ) %>%
  arrange(desc(mean_pct))

write_tsv(gene_summary, "tables/ARG_gene_prevalence_meanMaxPct.tsv")

# top 20 by mean%
top20 <- gene_summary %>% slice_head(n=20)

# Top genes by season (mean% within season)
gene_by_season <- rel_long %>%
  group_by(Season, GeneName) %>%
  summarise(
    prevalence_groups = sum(Rel > 0),
    mean_pct = 100*mean(Rel, na.rm=TRUE),
    max_pct  = 100*max(Rel, na.rm=TRUE),
    .groups="drop"
  ) %>%
  arrange(Season, desc(mean_pct))

write_tsv(gene_by_season, "tables/ARG_gene_bySeason_meanMaxPct.tsv")

# 6) Compose TXT report -------------------------------------------------------
report_file <- "tables/ARG_season_report.txt"

# Helper to format PERMANOVA table nicely
perm_tbl <- as.data.frame(perm)
perm_tbl$Term <- rownames(perm_tbl)
perm_tbl <- perm_tbl %>% select(Term, Df, SumOfSqs, R2, F, `Pr(>F)`)

lines <- c()
lines <- c(lines, "ARG PROFILE COMPARISON: July vs September (replicates united by median)")
lines <- c(lines, "====================================================================")
lines <- c(lines, "")
lines <- c(lines, sprintf("N groups (Base×Season): %d (Bases: %d; Seasons: %s)",
                          nrow(alpha_base), dplyr::n_distinct(alpha_base$Base),
                          paste(sort(unique(alpha_base$Season)), collapse=", ")))
lines <- c(lines, "")

# Overall alpha summary by season
alpha_season_sum <- alpha_base %>%
  group_by(Season) %>%
  summarise(
    n = n(),
    Total_median = median(Total_ARG_per16S),
    Total_mean   = mean(Total_ARG_per16S),
    Rich_median  = median(Richness),
    Rich_mean    = mean(Richness),
    Shan_median  = median(Shannon, na.rm=TRUE),
    Shan_mean    = mean(Shannon, na.rm=TRUE),
    Even_median  = median(Evenness_Pielou, na.rm=TRUE),
    Even_mean    = mean(Evenness_Pielou, na.rm=TRUE),
    .groups="drop"
  )

lines <- c(lines, "ALPHA DIVERSITY (Base×Season)")
lines <- c(lines, "----------------------------")
lines <- c(lines, capture.output(print(alpha_season_sum, n=50)))
lines <- c(lines, "")
lines <- c(lines, "Season tests (Wilcoxon; BH-adjusted p; Cliff's delta as effect size)")
lines <- c(lines, capture.output(print(season_stats, n=50)))
lines <- c(lines, "")

lines <- c(lines, "BETA DIVERSITY (Bray–Curtis on Base×Season profiles)")
lines <- c(lines, "---------------------------------------------------")
lines <- c(lines, "PERMANOVA (adonis2, by='margin'):")
lines <- c(lines, capture.output(print(perm_tbl, row.names = FALSE)))
lines <- c(lines, "")
lines <- c(lines, "PERMDISP check (dispersion differences can inflate PERMANOVA):")
lines <- c(lines, "ANOVA on betadisper distances:")
lines <- c(lines, capture.output(print(permdisp_aov)))
lines <- c(lines, "Permutation test on betadisper distances:")
lines <- c(lines, capture.output(print(permdisp_perm)))
lines <- c(lines, "")

lines <- c(lines, "TOP GENES: prevalence and contribution (%)")
lines <- c(lines, "-----------------------------------------")
lines <- c(lines, "Top 20 genes by mean relative abundance across Base×Season groups:")
lines <- c(lines, capture.output(print(top20, n=50)))
lines <- c(lines, "")
lines <- c(lines, "Top genes by season (mean% and max% within season):")
lines <- c(lines, capture.output(print(gene_by_season %>% group_by(Season) %>% slice_head(n=10), n=50)))
lines <- c(lines, "")

writeLines(lines, report_file)
cat("Wrote:", report_file, "\n")
######################### 1) upgrade
library(dplyr)
library(tidyr)
library(effsize)

paired_wilcox <- function(df, metric){
  wide <- df %>%
    select(Base, Season, !!sym(metric)) %>%
    pivot_wider(names_from = Season, values_from = !!sym(metric))
  
  xJ <- wide$July
  xS <- wide$September
  
  wt <- wilcox.test(xJ, xS, paired = TRUE, exact = FALSE)
  cd <- effsize::cliff.delta(xJ, xS)$estimate  # effect size (non paired-specific ma ok come indicatore)
  
  tibble(
    metric = metric,
    n = sum(!is.na(xJ) & !is.na(xS)),
    median_July = median(xJ, na.rm=TRUE),
    median_Sept = median(xS, na.rm=TRUE),
    p_value = wt$p.value,
    cliffs_delta = as.numeric(cd)
  )
}

paired_stats <- bind_rows(
  paired_wilcox(alpha_base, "Total_ARG_per16S"),
  paired_wilcox(alpha_base, "Richness"),
  paired_wilcox(alpha_base, "Shannon"),
  paired_wilcox(alpha_base, "Evenness_Pielou")
) %>% mutate(p_adj = p.adjust(p_value, method="BH"))

paired_stats
write_tsv(paired_stats, "tables/ARG_season_tests_PAIRED.tsv")
################################## 2) upgrade
library(dplyr)
library(tidyr)

topk_share <- rel_long %>%
  group_by(Base, Season) %>%
  arrange(desc(Rel), .by_group = TRUE) %>%
  summarise(
    Top1_pct  = 100 * sum(Rel[1], na.rm=TRUE),
    Top5_pct  = 100 * sum(Rel[1:5], na.rm=TRUE),
    Top10_pct = 100 * sum(Rel[1:10], na.rm=TRUE),
    .groups="drop"
  )

topk_season_summary <- topk_share %>%
  group_by(Season) %>%
  summarise(
    Top1_median = median(Top1_pct),
    Top5_median = median(Top5_pct),
    Top10_median = median(Top10_pct),
    .groups="drop"
  )

topk_share
topk_season_summary

write_tsv(topk_share, "tables/ARG_topK_share_byBaseSeason.tsv")
write_tsv(topk_season_summary, "tables/ARG_topK_share_season_summary.tsv")

top1_gene <- rel_long %>%
  group_by(Base, Season) %>%
  slice_max(order_by = Rel, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Top1_pct = 100*Rel) %>%
  arrange(Season, desc(Top1_pct))

top1_gene
write_tsv(top1_gene, "tables/ARG_top1_gene_byBaseSeason.tsv")
################################################################################
#### meta_g deve avere anche Site (prima lettera) associato a Base##############
# se non ce l'hai:
meta_g$Site <- substr(meta_g$Base, 1, 1)

# July only
idxJ <- meta_g$Season == "July"
bcJ  <- vegdist(X[idxJ, , drop=FALSE], method="bray")
metaJ <- meta_g[idxJ, , drop=FALSE]
adonis2(bcJ ~ Site, data=metaJ, permutations=999, by="margin")
bdJ <- betadisper(bcJ, metaJ$Site); anova(bdJ); permutest(bdJ, permutations=999)

# September only
idxS <- meta_g$Season == "September"
bcS  <- vegdist(X[idxS, , drop=FALSE], method="bray")
metaS <- meta_g[idxS, , drop=FALSE]
adonis2(bcS ~ Site, data=metaS, permutations=999, by="margin")
bdS <- betadisper(bcS, metaS$Site); anova(bdS); permutest(bdS, permutations=999)

#### 8) OPTIONAL: stacked bar by ARG CLASS (if you have annotations) ----------
# If you have a mapping table, e.g. arg_annotation.tsv with columns: GeneID, Class
 ann <- read_tsv("arg_annotation.tsv", show_col_types = FALSE)
 arg_base_gene_class <- arg_base_gene %>%
   left_join(ann, by = "GeneID") %>%
   mutate(Class = ifelse(is.na(Class) | Class == "", "Unclassified", Class)) %>%
   group_by(Base, Site, Type, Season, Class) %>%
   summarise(Abund_med = sum(Abund_med, na.rm = TRUE), .groups = "drop") %>%
   group_by(Base, Season) %>%
   mutate(rel = ifelse(sum(Abund_med) > 0, Abund_med / sum(Abund_med), 0)) %>%
   ungroup()

 Nclass <- 10
 topC <- arg_base_gene_class %>%
   group_by(Class) %>%
   summarise(mean_rel = mean(rel, na.rm = TRUE), .groups="drop") %>%
   arrange(desc(mean_rel)) %>%
   slice_head(n=Nclass) %>%
   pull(Class)
#
 arg_class_plot <- arg_base_gene_class %>%
   mutate(Class2 = ifelse(Class %in% topC, Class, "Other classes")) %>%
   group_by(Base, Site, Type, Season, Class2) %>%
   summarise(rel = sum(rel, na.rm = TRUE), .groups="drop") %>%
   group_by(Base, Season) %>%
   mutate(rel = ifelse(sum(rel) > 0, rel/sum(rel), 0)) %>%
   ungroup()
#
 levsC <- c(sort(setdiff(unique(arg_class_plot$Class2), "Other classes")), "Other classes")
 arg_class_plot <- arg_class_plot %>% mutate(Class2 = factor(Class2, levels = levsC))
#
 p_class <- ggplot(arg_class_plot, aes(x=Base, y=rel, fill=Class2)) +
   geom_col(color="white", linewidth=0.2) +
   facet_wrap(~Season, nrow=1) +
   scale_y_continuous(labels=function(x) paste0(round(100*x), "%")) +
 labs(x="Base", y="ARG class composition (median across replicates)", fill="Class") +
   theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5))
#
# ggsave("figures/ARG_ClassComp_Top10_byBaseSeason.png", p_class, width=12, height=5, dpi=600)
# write_tsv(arg_class_plot, "tables/ARG_class_composition_byBaseSeason_Top10.tsv")
 ###############################################################################
 # ARG publication summary report (MRG-like) -> tables/ARG_publication_summary.txt
 ###############################################################################
 
 pkgs <- c("dplyr","tidyr","stringr","readr","vegan","effsize")
 to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
 if(length(to_install)) install.packages(to_install, repos="https://cloud.r-project.org")
 invisible(lapply(pkgs, library, character.only = TRUE))
 
 dir.create("tables", showWarnings = FALSE)
 
 # --------- INPUTS ---------
 arg_file <- "arg_per16s_T1T2_renamed.tsv"   # gene x sample (ARG/16S)
 out_txt  <- "tables/ARG_publication_summary.txt"
 
 # taxa to drop from gene labels (edit as you like)
 drop_terms <- c("streptomyces","pseudomonas","acinetobacter","staphylococcus")
 
 # --------- 1) Read wide matrix ---------
 arg_wide <- read_tsv(arg_file, show_col_types = FALSE)
 gene_col <- names(arg_wide)[1]
 sample_cols <- setdiff(names(arg_wide), gene_col)
 
 # force numeric
 arg_wide <- arg_wide %>%
   mutate(across(all_of(sample_cols), ~ as.numeric(.x)))
 
 # --------- 2) Per-sample alpha summary (headline numbers) ---------
 # N genes detected at least once
 n_ARG <- arg_wide %>%
   mutate(present = rowSums(across(all_of(sample_cols), ~ .x > 0), na.rm=TRUE) > 0) %>%
   summarise(n = sum(present, na.rm=TRUE)) %>% pull(n)
 
 # per-sample totals & richness
 per_sample <- arg_wide %>%
   pivot_longer(cols = all_of(sample_cols), names_to="Sample", values_to="Abund") %>%
   mutate(Abund = as.numeric(Abund)) %>%
   group_by(Sample) %>%
   summarise(
     Total_ARG_per16S = sum(Abund, na.rm=TRUE),
     Richness_ARG     = sum(Abund > 0, na.rm=TRUE),
     .groups="drop"
   )
 
 med_total <- median(per_sample$Total_ARG_per16S, na.rm=TRUE)
 med_rich  <- median(per_sample$Richness_ARG, na.rm=TRUE)
 
 write_tsv(per_sample, "tables/ARG_alpha_perSample.tsv")
 
 # --------- 3) Long format + metadata from Sample names ---------
 arg_long <- arg_wide %>%
   rename(GeneID = !!sym(gene_col)) %>%
   pivot_longer(cols = all_of(sample_cols), names_to="Sample", values_to="Abundance") %>%
   mutate(
     Abundance = as.numeric(Abundance),
     GeneID = as.character(GeneID),
     # label after "|" if present
     GeneName = ifelse(str_detect(GeneID, "\\|"), str_replace(GeneID, "^.*\\|", ""), GeneID),
     GeneName = as.character(GeneName),
     # parse metadata from sample code (same logic as before)
     Base = str_sub(Sample, 1, 2),
     Site = str_sub(Sample, 1, 1),
     Season_code = str_sub(Sample, 3, 3),
     Season = recode(Season_code, "J"="July", "S"="September", .default = Season_code),
     Replicate = str_extract(Sample, "\\d+$"),
     Type = ifelse(str_sub(Base, 2, 2) == "R", "RG", "non-RG")
   ) %>%
   # drop “taxa-like” gene labels
   filter(!tolower(GeneName) %in% drop_terms)
 
 # --------- 4) Unite replicates: median per Base×Season×Gene ---------
 arg_base_gene <- arg_long %>%
   group_by(Base, Site, Type, Season, GeneName) %>%
   summarise(Abund_med = median(Abundance, na.rm=TRUE), .groups="drop") %>%
   mutate(Abund_med = ifelse(is.na(Abund_med), 0, Abund_med))
 
 # build relative composition
 rel_long <- arg_base_gene %>%
   group_by(Base, Season) %>%
   mutate(
     Total = sum(Abund_med, na.rm=TRUE),
     Rel   = ifelse(Total > 0, Abund_med/Total, 0)
   ) %>%
   ungroup()
 
 # --------- 5) Dominance (Top1/Top5/Top10 shares) ---------
 topk_share <- rel_long %>%
   group_by(Base, Season) %>%
   arrange(desc(Rel), .by_group = TRUE) %>%
   summarise(
     Top1_pct  = 100*sum(Rel[1], na.rm=TRUE),
     Top5_pct  = 100*sum(Rel[1:5], na.rm=TRUE),
     Top10_pct = 100*sum(Rel[1:10], na.rm=TRUE),
     .groups="drop"
   )
 
 top1_summary <- topk_share %>%
   summarise(
     Top1_median = median(Top1_pct, na.rm=TRUE),
     Top1_min    = min(Top1_pct, na.rm=TRUE),
     Top1_max    = max(Top1_pct, na.rm=TRUE),
     Top5_median = median(Top5_pct, na.rm=TRUE),
     Top10_median= median(Top10_pct, na.rm=TRUE)
   )
 
 write_tsv(topk_share, "tables/ARG_topK_share_byBaseSeason.tsv")
 
 # top1 gene identity per Base×Season
 top1_gene <- rel_long %>%
   group_by(Base, Season) %>%
   slice_max(order_by = Rel, n = 1, with_ties = FALSE) %>%
   ungroup() %>%
   mutate(Top1_pct = 100*Rel) %>%
   arrange(Season, desc(Top1_pct))
 
 write_tsv(top1_gene, "tables/ARG_top1_gene_byBaseSeason.tsv")
 
 # --------- 6) Top genes table (prevalence + mean% + max% + median%) ---------
 gene_summary <- rel_long %>%
   group_by(GeneName) %>%
   summarise(
     prevalence_groups = sum(Rel > 0),
     prevalence_pct    = 100*mean(Rel > 0),
     mean_pct          = 100*mean(Rel, na.rm=TRUE),
     median_pct        = 100*median(Rel, na.rm=TRUE),
     max_pct           = 100*max(Rel, na.rm=TRUE),
     .groups="drop"
   ) %>%
   arrange(desc(mean_pct))
 
 top20_genes <- gene_summary %>% slice_head(n=20)
 
 write_tsv(gene_summary, "tables/ARG_gene_prevalence_meanMaxPct.tsv")
 write_tsv(top20_genes,  "tables/ARG_top20_genes_meanPct.tsv")
 
 # top gene by mean%
 top_gene <- gene_summary %>% slice_head(n=1)
 
 # --------- 7) Beta diversity (Bray–Curtis) + PERMANOVA ---------
 mat_wide <- rel_long %>%
   group_by(Base, Site, Type, Season, GeneName) %>%
   summarise(Abund = sum(Abund_med, na.rm=TRUE), .groups="drop") %>%
   unite(Group, Base, Season, sep="_") %>%
   pivot_wider(names_from = GeneName, values_from = Abund, values_fill = 0)
 
 meta_g <- mat_wide %>%
   select(Group) %>%
   separate(Group, into=c("Base","Season"), sep="_", remove=FALSE) %>%
   mutate(
     Site = substr(Base, 1, 1),
     Type = ifelse(substr(Base, 2, 2) == "R", "RG", "non-RG")
   )
 
 X <- mat_wide %>% select(-Group) %>% as.data.frame()
 bc <- vegdist(X, method="bray")
 
 # PERMANOVA: Base + Season + Type (Site is redundant with Base; report one of them)
 perm_base_season_type <- adonis2(bc ~ Base + Season + Type, data=meta_g, permutations=999, by="margin")
 
 # --------- 8) Within-site pairs vs between-site distances (MRG-like) ---------
 # Define within-site pairs by same first letter, different second letter (e.g., BF–BR)
 # We'll compute distances between Base profiles within each Season.
 dist_df <- as.data.frame(as.matrix(bc))
 dist_df$Group <- rownames(dist_df)
 
 pairs_long <- dist_df %>%
   pivot_longer(cols = -Group, names_to="Group2", values_to="BC") %>%
   filter(Group < Group2) %>%
   separate(Group, into=c("Base1","Season1"), sep="_", remove=FALSE) %>%
   separate(Group2, into=c("Base2","Season2"), sep="_", remove=FALSE) %>%
   filter(Season1 == Season2) %>%
   mutate(
     Season = Season1,
     Site1 = substr(Base1,1,1),
     Site2 = substr(Base2,1,1),
     within_site = (Site1 == Site2) & (Base1 != Base2)
   ) %>%
   select(Season, Base1, Base2, within_site, BC)
 
 # Wilcoxon within-site vs between-site, per season and overall
 wilcox_dist <- function(df){
   x <- df$BC[df$within_site]
   y <- df$BC[!df$within_site]
   if(length(x) < 1 || length(y) < 1) return(tibble(p_value=NA_real_, n_within=length(x), n_between=length(y)))
   tibble(
     n_within = length(x),
     n_between = length(y),
     median_within = median(x),
     median_between = median(y),
     p_value = wilcox.test(x, y, exact=FALSE)$p.value
   )
 }
 
 dist_overall <- wilcox_dist(pairs_long)
 dist_by_season <- pairs_long %>%
   group_by(Season) %>%
   group_modify(~ wilcox_dist(.x)) %>%
   ungroup()
 
 write_tsv(pairs_long, "tables/ARG_braycurtis_pairwise_within_between.tsv")
 
 # --------- 9) Write publication TXT ---------
 lines <- c()
 lines <- c(lines, "ARG LANDSCAPE SUMMARY (publication-ready numbers)")
 lines <- c(lines, "==============================================")
 lines <- c(lines, "")
 lines <- c(lines, sprintf("Distinct ARGs detected (present in ≥1 sample): %d", n_ARG))
 lines <- c(lines, sprintf("Per-sample (ARG/16S): median total = %.3f ; median richness = %.0f genes (n=%d samples)",
                           med_total, med_rich, nrow(per_sample)))
 lines <- c(lines, "")
 
 lines <- c(lines, "Dominance (Base×Season profiles; replicates united by median):")
 lines <- c(lines, sprintf("Top1%% median = %.2f%% (range %.2f–%.2f%%); Top5%% median = %.2f%%; Top10%% median = %.2f%%",
                           top1_summary$Top1_median, top1_summary$Top1_min, top1_summary$Top1_max,
                           top1_summary$Top5_median, top1_summary$Top10_median))
 lines <- c(lines, "")
 
 lines <- c(lines, "Most prevalent/dominant genes (Base×Season relative composition):")
 lines <- c(lines, sprintf("Top gene by mean%%: %s (mean %.2f%%; median %.2f%%; max %.2f%%; prevalence %.1f%% of groups)",
                           top_gene$GeneName, top_gene$mean_pct, top_gene$median_pct, top_gene$max_pct, top_gene$prevalence_pct))
 lines <- c(lines, "")
 lines <- c(lines, "Top 20 genes by mean% (mean%, max%, prevalence%):")
 lines <- c(lines, capture.output(print(top20_genes %>% select(GeneName, prevalence_pct, mean_pct, max_pct), n=50)))
 lines <- c(lines, "")
 
 lines <- c(lines, "Beta diversity (Bray–Curtis on Base×Season profiles) PERMANOVA by margin:")
 lines <- c(lines, capture.output(print(as.data.frame(perm_base_season_type))))
 lines <- c(lines, "")
 
 lines <- c(lines, "Centroid-distance style check: within-site pairs vs between-site comparisons (Bray–Curtis distances)")
 lines <- c(lines, "Overall (all seasons pooled):")
 lines <- c(lines, capture.output(print(dist_overall)))
 lines <- c(lines, "By season:")
 lines <- c(lines, capture.output(print(dist_by_season)))
 lines <- c(lines, "")
 
 writeLines(lines, out_txt)
 cat("Wrote:", out_txt, "\n")
 
 library(vegan)
 library(dplyr)
 library(tidyr)
 library(readr)
 
 arg_wide <- read_tsv("arg_per16s_T1T2_renamed.tsv", show_col_types = FALSE)
 gene_col <- names(arg_wide)[1]
 sample_cols <- setdiff(names(arg_wide), gene_col)
 
 X <- arg_wide %>%
   mutate(across(all_of(sample_cols), as.numeric)) %>%
   select(all_of(sample_cols)) %>%
   t() %>% as.data.frame()
 
 bc <- vegdist(X, method="bray")
 d <- as.vector(bc)
 
 cat(sprintf("Bray–Curtis dissimilarity across samples: median %.3f (range %.3f–%.3f)\n",
             median(d), min(d), max(d)))
 
 library(dplyr)
 
 # 1) Top1 gene per ciascun Base×Season
 top1_gene <- rel_long %>%
   group_by(Base, Season) %>%
   slice_max(order_by = Rel, n = 1, with_ties = FALSE) %>%
   ungroup() %>%
   mutate(Top1_pct = 100*Rel)
 
 # (a) Il gene più frequentemente Top1
 top1_most_frequent <- top1_gene %>%
   count(GeneName, sort = TRUE) %>%
   slice_head(n=10)
 
 # (b) Il Top1 più alto in assoluto (picco)
 top1_highest_peak <- top1_gene %>%
   arrange(desc(Top1_pct)) %>%
   slice_head(n=10)
 
 # (c) “Top1 median%” per gene (quando è Top1)
 top1_median_by_gene <- top1_gene %>%
   group_by(GeneName) %>%
   summarise(
     n_groups = n(),
     median_top1_pct = median(Top1_pct),
     max_top1_pct    = max(Top1_pct),
     .groups="drop"
   ) %>%
   arrange(desc(median_top1_pct))
 
 print(top1_most_frequent, n=10)
 print(top1_highest_peak, n=10)
 print(top1_median_by_gene, n=20)
 
 write_tsv(top1_gene, "tables/ARG_top1_gene_byBaseSeason.tsv")
 write_tsv(top1_most_frequent, "tables/ARG_top1_most_frequent.tsv")
 write_tsv(top1_highest_peak, "tables/ARG_top1_highest_peaks.tsv")
 write_tsv(top1_median_by_gene, "tables/ARG_top1_median_by_gene.tsv")
 
 ##########################################################################
 
 library(dplyr)
 library(tidyr)
 library(stringr)
 library(vegan)
 
 # 1) Ricostruisci matrice Base×Season x gene (da rel_long / arg_base_gene)
 # Usa Abund_med (non Rel) per Bray–Curtis
 mat_wide <- arg_base_gene %>%   # <-- deve esistere: Base, Season, GeneName, Abund_med
   group_by(Base, Season, GeneName) %>%
   summarise(Abund = sum(Abund_med, na.rm=TRUE), .groups="drop") %>%
   mutate(SeasonCode = ifelse(Season=="July","J",
                              ifelse(Season=="September","S", Season))) %>%
   mutate(Group = paste0(Base, SeasonCode)) %>%   # es: BFJ, BFS
   select(Group, GeneName, Abund) %>%
   pivot_wider(names_from = GeneName, values_from = Abund, values_fill = 0)
 
 # 2) Metti Group come rownames in modo esplicito
 X <- mat_wide %>% select(-Group) %>% as.data.frame()
 rownames(X) <- mat_wide$Group
 
 # 3) Controlli fondamentali
 stopifnot(nrow(X) > 2)
 print(head(rownames(X), 20))
 print(table(str_sub(rownames(X), 3, 3)))  # deve dare J e S
 
 # 4) Distanze
 bc <- vegdist(X, method="bray")
 D  <- as.matrix(bc)
 
 # 5) Pairwise table
 pairs <- expand.grid(i = seq_len(nrow(D)), j = seq_len(nrow(D))) %>%
   filter(i < j) %>%
   mutate(
     G1 = rownames(D)[i],
     G2 = colnames(D)[j],
     BC = D[cbind(i, j)],
     Base1 = str_sub(G1, 1, 2),
     Base2 = str_sub(G2, 1, 2),
     SeasonCode1 = str_sub(G1, 3, 3),
     SeasonCode2 = str_sub(G2, 3, 3),
     Season1 = recode(SeasonCode1, "J"="July", "S"="September", .default = NA_character_),
     Season2 = recode(SeasonCode2, "J"="July", "S"="September", .default = NA_character_),
     Site1 = str_sub(Base1, 1, 1),
     Site2 = str_sub(Base2, 1, 1)
   ) %>%
   filter(!is.na(Season1), Season1 == Season2) %>%
   mutate(
     Season = Season1,
     within_site = (Site1 == Site2) & (Base1 != Base2)
   )
 
 # 6) Check: ora NON deve essere vuoto
 print(dim(pairs))
 print(table(pairs$Season, pairs$within_site))
 
 # 7) Wilcoxon (solo se esistono within e between)
 wilcox_dist <- function(df){
   x <- df$BC[df$within_site]
   y <- df$BC[!df$within_site]
   if(length(x) < 2 || length(y) < 2){
     return(tibble(
       n_within = length(x),
       n_between = length(y),
       median_within = ifelse(length(x)>0, median(x), NA_real_),
       median_between = ifelse(length(y)>0, median(y), NA_real_),
       p_value = NA_real_
     ))
   }
   tibble(
     n_within = length(x),
     n_between = length(y),
     median_within = median(x),
     median_between = median(y),
     p_value = wilcox.test(x, y, exact=FALSE)$p.value
   )
 }
 
 dist_overall <- wilcox_dist(pairs)
 dist_by_season <- pairs %>% group_by(Season) %>% group_modify(~wilcox_dist(.x)) %>% ungroup()
 
 print(dist_overall)
 print(dist_by_season)
 
 # ==== Stacked bar: composition by Antibiotic Class (replicates aggregated) ====
 
 pkgs <- c("readr","dplyr","tibble","tidyr","stringr","ggplot2","RColorBrewer")
 to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
 if(length(to_install)) install.packages(to_install, repos="https://cloud.r-project.org")
 invisible(lapply(pkgs, library, character.only = TRUE))
 
 dir.create("figures", showWarnings = FALSE, recursive = TRUE)
 dir.create("tables",  showWarnings = FALSE, recursive = TRUE)
 
 # ---- INPUTS ----
 arg_matrix_file  <- "ARG_per16S_T1T2_renamed.tsv"      # genes x samples (prima col = ARG_ID)
 arg2class_file   <- "ARG_annotation_RGI_T1T2.tsv"      # colonne: ARG_ID, BestHit, DrugClass, ...
 metadata_file    <- "metadata2.txt"                    # opzionale ma consigliato
 
 # ---- helper: pulizia gene/ARG label ----
 clean_gene <- function(x){
   x <- as.character(x)
   x <- sub("^.*\\|", "", x, perl = TRUE)     # dopo ultimo |
   x <- sub("^[0-9]+", "", x, perl = TRUE)    # rimuove numeri iniziali
   x <- sub("^[-_ :]+", "", x, perl = TRUE)   # rimuove separatori iniziali
   stringr::str_trim(x)
 }
 
 # ---- 1) Read ARG matrix -> long ----
 arg_raw <- readr::read_tsv(arg_matrix_file, show_col_types = FALSE)
 arg_id_col <- names(arg_raw)[1]
 
 ARG_long <- arg_raw %>%
   tidyr::pivot_longer(-all_of(arg_id_col), names_to = "Sample", values_to = "Abundance") %>%
   dplyr::mutate(
     Abundance = as.numeric(Abundance),
     ARG_ID = as.character(.data[[arg_id_col]]),
     ARG_gene_clean = clean_gene(ARG_ID)
   ) %>%
   dplyr::select(Sample, ARG_ID, ARG_gene_clean, Abundance)
 
 # ---- 2) Read mapping ARG/BestHit -> DrugClass ----
 ab_raw <- readr::read_tsv(arg2class_file, show_col_types = FALSE)
 names(ab_raw) <- trimws(names(ab_raw))
 
 # attese: ARG_ID / BestHit / DrugClass
 stopifnot(all(c("ARG_ID","BestHit","DrugClass") %in% names(ab_raw)))
 
 ab_map <- ab_raw %>%
   dplyr::transmute(
     ARG_ID  = as.character(ARG_ID),
     BestHit = as.character(BestHit),
     DrugClass = as.character(DrugClass)
   ) %>%
   dplyr::mutate(
     DrugClass = stringr::str_trim(DrugClass),
     Gene_from_ARGID   = clean_gene(ARG_ID),
     Gene_from_BestHit = clean_gene(BestHit)
   ) %>%
   tidyr::pivot_longer(
     cols = c(Gene_from_ARGID, Gene_from_BestHit),
     names_to = "src", values_to = "Gene"
   ) %>%
   dplyr::mutate(Gene = stringr::str_trim(Gene)) %>%
   dplyr::filter(!is.na(Gene), Gene != "", !is.na(DrugClass), DrugClass != "") %>%
   tidyr::separate_rows(DrugClass, sep = "\\s*[,;/|]+\\s*") %>%  # se più classi per gene
   dplyr::mutate(DrugClass = stringr::str_trim(DrugClass)) %>%
   dplyr::filter(DrugClass != "") %>%
   dplyr::distinct(Gene, DrugClass)
 
 # NOTE: se un gene mappa a più classi, qui c'è "doppio conteggio" (come per metalli).
 # Se vuoi una composizione più conservativa, si può dividere Abundance per n_classi per gene.
 
 # ---- 3) Sample metadata (join) ----
 if(file.exists(metadata_file)){
   meta <- readr::read_tsv(metadata_file, show_col_types = FALSE)
   # prova a trovare una colonna Sample
   if(!("Sample" %in% names(meta))){
     # prova alternative
     cand <- names(meta)[tolower(names(meta)) %in% c("sample","id")]
     if(length(cand) == 0) stop("metadata_file non ha una colonna Sample (o sample/id).")
     meta <- meta %>% dplyr::rename(Sample = all_of(cand[1]))
   }
   meta <- meta %>% dplyr::mutate(Sample = as.character(Sample))
 } else {
   meta <- tibble(Sample = unique(ARG_long$Sample))
 }
 
 # fallback se mancano colonne (adatta se la codifica del Sample è diversa)
 if(!("Base" %in% names(meta)))   meta$Base   <- substr(meta$Sample, 1, 2)
 if(!("Season" %in% names(meta))) meta$Season <- substr(meta$Sample, 3, 3)
 if(!("Type" %in% names(meta)))   meta$Type   <- substr(meta$Sample, 4, 4)
 if(!("Site" %in% names(meta)))   meta$Site   <- meta$Base  # se non hai Site separato
 
 # ---- 4) Sum ARG abundances per sample per DrugClass ----
 abclass_sample <- ARG_long %>%
   dplyr::left_join(ab_map, by = c("ARG_gene_clean" = "Gene")) %>%
   dplyr::mutate(DrugClass = dplyr::if_else(is.na(DrugClass), "Unclassified", DrugClass)) %>%
   dplyr::left_join(meta, by = "Sample") %>%
   dplyr::group_by(Sample, Base, Site, Type, Season, DrugClass) %>%
   dplyr::summarise(Class_abund = sum(Abundance, na.rm = TRUE), .groups="drop")
 
 # ---- 5) Aggregate replicates: median per Base×Season×DrugClass (mantengo Site/Type come nel tuo) ----
 abclass_base <- abclass_sample %>%
   dplyr::group_by(Base, Site, Type, Season, DrugClass) %>%
   dplyr::summarise(Class_abund_med = median(Class_abund, na.rm=TRUE), .groups="drop")
 
 # ---- 6) Convert to relative composition within Base×Season ----
 abclass_base <- abclass_base %>%
   dplyr::group_by(Base, Season) %>%
   dplyr::mutate(Class_rel = Class_abund_med / sum(Class_abund_med, na.rm=TRUE)) %>%
   dplyr::ungroup()
 
 # ---- 7) Top N classes + Other ----
 N <- 10
 topN <- abclass_base %>%
   dplyr::group_by(DrugClass) %>%
   dplyr::summarise(mean_rel = mean(Class_rel, na.rm=TRUE), .groups="drop") %>%
   dplyr::arrange(dplyr::desc(mean_rel)) %>%
   dplyr::slice_head(n=N) %>%
   dplyr::pull(DrugClass)
 
 abclass_plot <- abclass_base %>%
   dplyr::mutate(DrugClass2 = ifelse(DrugClass %in% topN, DrugClass, "Other")) %>%
   dplyr::group_by(Base, Site, Type, Season, DrugClass2) %>%
   dplyr::summarise(Class_rel = sum(Class_rel, na.rm=TRUE), .groups="drop") %>%
   dplyr::group_by(Base, Season) %>%
   dplyr::mutate(Class_rel = Class_rel / sum(Class_rel, na.rm=TRUE)) %>%
   dplyr::ungroup()
 
 # ordine legenda: topN poi Other
 abclass_plot <- abclass_plot %>%
   dplyr::mutate(DrugClass2 = factor(DrugClass2, levels = c(setdiff(topN, "Other"), "Other")))
 
 # ---- 8) Plot ----
 # (in cima allo script, aggiungi il pacchetto)
 # pkgs <- c(..., "RColorBrewer")
 
 # ---- colori: pastello per classi, Other grigio ----
 levs <- levels(abclass_plot$DrugClass2)
 
 # prendi una palette pastello (Set3 è molto usata e morbida)
 pastel <- RColorBrewer::brewer.pal(max(3, min(length(levs), 12)), "Set3")
 
 # assegna colori a tutti tranne "Other"
 cls <- setdiff(levs, "Other")
 col_map <- setNames(rep(pastel, length.out = length(cls)), cls)
 
 # forza Other grigio
 col_map <- c(col_map, Other = "grey70")
 
 abclass_plot <- abclass_plot %>%
   mutate(Season = recode(Season, "J" = "July", "S" = "September"))
 
 p <- ggplot(abclass_plot, aes(x = Base, y = Class_rel, fill = DrugClass2)) +
   geom_col(color="white", linewidth=0.2) +
   facet_wrap(~Season, nrow=1) +
   scale_y_continuous(labels = function(x) paste0(round(100*x), "%")) +
   scale_fill_manual(values = col_map, drop = FALSE) +
   labs(
     x = "Site",
     y = "ARG drug-class composition (median across replicates)",
     fill = "Drug class"
   ) +
   theme_bw() +
   theme(axis.text.x = element_text(angle=90, vjust=0.5))
 
 ggsave("figures/ARG_drugclass_composition_TopN_byBaseSeason.png", p, width=12, height=5, dpi=600) 
 readr::write_tsv(abclass_plot, "tables/ARG_drugclass_composition_byBaseSeason_TopN.tsv")
 
 ############################################################
 ## HEATMAP like example: ARG Antibiotic classes (merged reps, log10)
 ## - rows: DrugClass (top N + optional "Other")
 ## - cols: Group = Base+Season (BFJ, BFS, ...)
 ## - reps merged: median across replicates
 ## - values: log10(relative composition + pseudo)  -> range ~ (-6..0)
 ## - top annotation: Base colours
 ############################################################
 pkgs <- c("readr","dplyr","tibble","tidyr","stringr","pheatmap","RColorBrewer")
 to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
 if(length(to_install)) install.packages(to_install, repos="https://cloud.r-project.org")
 invisible(lapply(pkgs, library, character.only = TRUE))
 
 dir.create("figures", showWarnings = FALSE, recursive = TRUE)
 dir.create("tables",  showWarnings = FALSE, recursive = TRUE)
 
 # ---- INPUTS ----
 arg_matrix_file  <- "ARG_per16S_T1T2_renamed.tsv"      # genes x samples (prima col = ARG_ID)
 arg2class_file   <- "ARG_annotation_RGI_T1T2.tsv"      # colonne: ARG_ID, BestHit, DrugClass
 metadata_file    <- "metadata2.txt"                    # opzionale
 
 # ---- SETTINGS ----
 TOP_N_CLASSES <- 18        # cambia (es 12/15/20)
 USE_OTHER     <- TRUE      # TRUE = collassa non-top in "Other"
 PSEUDO_REL    <- 1e-6      # per log10(relative + pseudo)
 OUT_PREFIX    <- "HEATMAP_ARG_drugClasses_topN"
 
 # ---- palette Base (prime 2 lettere) ----
 base_cols_map <- c(
   BF = "#1b9e77", BR = "#a6dba0",
   CF = "#d95f02", CR = "#fdb863",
   PF = "#7570b3", PR = "#b2abd2",
   SF = "#1f9ac2", SR = "#a6dce7",
   VF = "#e7298a", VR = "#f2b2d4"
 )
 
 # ---- helper: pulizia gene/ARG label ----
 clean_gene <- function(x){
   x <- as.character(x)
   x <- sub("^.*\\|", "", x, perl = TRUE)     # dopo ultimo |
   x <- sub("^[0-9]+", "", x, perl = TRUE)    # rimuove numeri iniziali
   x <- sub("^[-_ :]+", "", x, perl = TRUE)   # rimuove separatori iniziali
   stringr::str_trim(x)
 }
 
 # ---- 1) Read ARG matrix -> long ----
 arg_raw <- readr::read_tsv(arg_matrix_file, show_col_types = FALSE)
 arg_id_col <- names(arg_raw)[1]
 
 ARG_long <- arg_raw %>%
   tidyr::pivot_longer(-all_of(arg_id_col), names_to = "Sample", values_to = "Abundance") %>%
   dplyr::mutate(
     Abundance = suppressWarnings(as.numeric(Abundance)),
     Abundance = ifelse(is.na(Abundance), 0, Abundance),
     ARG_ID = as.character(.data[[arg_id_col]]),
     ARG_gene_clean = clean_gene(ARG_ID)
   ) %>%
   dplyr::select(Sample, ARG_ID, ARG_gene_clean, Abundance)
 
 # ---- 2) Read mapping -> DrugClass ----
 ab_raw <- readr::read_tsv(arg2class_file, show_col_types = FALSE)
 names(ab_raw) <- trimws(names(ab_raw))
 stopifnot(all(c("ARG_ID","BestHit","DrugClass") %in% names(ab_raw)))
 
 ab_map <- ab_raw %>%
   dplyr::transmute(
     ARG_ID  = as.character(ARG_ID),
     BestHit = as.character(BestHit),
     DrugClass = as.character(DrugClass)
   ) %>%
   dplyr::mutate(
     DrugClass = stringr::str_trim(DrugClass),
     Gene_from_ARGID   = clean_gene(ARG_ID),
     Gene_from_BestHit = clean_gene(BestHit)
   ) %>%
   tidyr::pivot_longer(
     cols = c(Gene_from_ARGID, Gene_from_BestHit),
     names_to = "src", values_to = "Gene"
   ) %>%
   dplyr::mutate(Gene = stringr::str_trim(Gene)) %>%
   dplyr::filter(!is.na(Gene), Gene != "", !is.na(DrugClass), DrugClass != "") %>%
   tidyr::separate_rows(DrugClass, sep = "\\s*[,;/|]+\\s*") %>%  # se più classi per gene
   dplyr::mutate(DrugClass = stringr::str_trim(DrugClass)) %>%
   dplyr::filter(DrugClass != "") %>%
   dplyr::distinct(Gene, DrugClass)
 
 # ---- 3) Metadata (optional) ----
 if (file.exists(metadata_file)) {
   meta <- readr::read_tsv(metadata_file, show_col_types = FALSE)
   if(!("Sample" %in% names(meta))){
     cand <- names(meta)[tolower(names(meta)) %in% c("sample","id")]
     if(length(cand) == 0) stop("metadata_file non ha una colonna Sample (o sample/id).")
     meta <- meta %>% dplyr::rename(Sample = all_of(cand[1]))
   }
   meta <- meta %>% mutate(Sample = as.character(Sample))
 } else {
   meta <- tibble(Sample = unique(ARG_long$Sample))
 }
 
 # fallback parse dal nome campione
 if(!("Base"   %in% names(meta))) meta$Base   <- substr(meta$Sample, 1, 2)
 if(!("Season" %in% names(meta))) meta$Season <- substr(meta$Sample, 3, 3)
 
 meta <- meta %>%
   mutate(
     Base   = as.character(Base),
     Season = as.character(Season),
     Group  = paste0(Base, Season)   # BFJ, BFS...
   )
 
 # ---- 4) Sum ARG abundances per sample per DrugClass ----
 abclass_sample <- ARG_long %>%
   left_join(ab_map,
             by = c("ARG_gene_clean" = "Gene"),
             relationship = "many-to-many") %>%  # silenzia warning (atteso)
   mutate(DrugClass = if_else(is.na(DrugClass), "Unclassified", DrugClass)) %>%
   left_join(meta %>% select(Sample, Base, Season, Group), by = "Sample") %>%
   group_by(Sample, Base, Season, Group, DrugClass) %>%
   summarise(Class_abund = sum(Abundance, na.rm = TRUE), .groups="drop")
 
 # ---- 5) Merge replicates: median per Group x DrugClass ----
 # (qui Group = Base+Season, quindi mediana sulle repliche dentro BFJ ecc.)
 abclass_group <- abclass_sample %>%
   group_by(Group, Base, Season, DrugClass) %>%
   summarise(Class_abund_med = median(Class_abund, na.rm=TRUE), .groups="drop")
 
 # ---- 6) Relative composition within Group ----
 abclass_group <- abclass_group %>%
   group_by(Group) %>%
   mutate(rel = Class_abund_med / sum(Class_abund_med, na.rm=TRUE)) %>%
   ungroup()
 
 # ---- 7) Top N classes (+ optional Other) ----
 class_rank <- abclass_group %>%
   group_by(DrugClass) %>%
   summarise(mean_rel = mean(rel, na.rm=TRUE), .groups="drop") %>%
   arrange(desc(mean_rel))
 
 n_keep <- min(TOP_N_CLASSES, nrow(class_rank))
 top_classes <- class_rank %>% slice_head(n = n_keep) %>% pull(DrugClass)
 
 if (USE_OTHER) {
   abclass_group2 <- abclass_group %>%
     mutate(DrugClass2 = ifelse(DrugClass %in% top_classes, DrugClass, "Other")) %>%
     group_by(Group, Base, Season, DrugClass2) %>%
     summarise(rel = sum(rel, na.rm=TRUE), .groups="drop") %>%
     group_by(Group) %>%
     mutate(rel = rel / sum(rel, na.rm=TRUE)) %>%
     ungroup()
 } else {
   abclass_group2 <- abclass_group %>%
     filter(DrugClass %in% top_classes) %>%
     mutate(DrugClass2 = DrugClass)
 }
 
 # ---- 8) Matrix rows=DrugClass, cols=Group ----
 mat <- abclass_group2 %>%
   select(DrugClass = DrugClass2, Group, rel) %>%
   pivot_wider(names_from = Group, values_from = rel, values_fill = 0)
 
 X_rel <- as.matrix(mat[,-1])
 rownames(X_rel) <- mat$DrugClass
 
 # log10(relative + pseudo) -> valori negativi tipo heatmap esempio
 X_pct <- 100 * X_rel   # X_rel è la composizione
 FIX_MIN <- 0
 FIX_MAX <- 30          # oppure 50
 
 bk <- seq(FIX_MIN, FIX_MAX, length.out = 101)
 cols <- viridisLite::viridis(100)
 X_cap <- pmin(pmax(X_pct, FIX_MIN), FIX_MAX)
 
 # salva matrici
 write_tsv(as_tibble(X_rel, rownames="DrugClass"),
           sprintf("tables/%s_matrix_relative.tsv", OUT_PREFIX))
 write_tsv(class_rank,
           sprintf("tables/%s_class_rank.tsv", OUT_PREFIX))
 
 # ---- 9) Column annotation: Base colours ----
 groups <- colnames(X_log)
 base_lab <- substr(groups, 1, 2)
 
 ann_col <- data.frame(Base = factor(base_lab, levels = sort(unique(base_lab))))
 rownames(ann_col) <- groups
 
 missing <- setdiff(levels(ann_col$Base), names(base_cols_map))
 if (length(missing) > 0) {
   extra <- RColorBrewer::brewer.pal(max(3, min(8, length(missing))), "Set2")
   base_cols_map[missing] <- extra[seq_along(missing)]
 }
 ann_colors <- list(Base = base_cols_map[levels(ann_col$Base)])
 
 # (opzionale) ordine colonne: dentro Base, J poi S
 col_order <- tibble(Group = colnames(X_log)) %>%
   mutate(Base = substr(Group, 1, 2),
          Season = substr(Group, 3, 3),
          Season = factor(Season, levels = c("J","S"))) %>%
   arrange(Base, Season) %>%
   pull(Group)
 X_log <- X_log[, col_order, drop = FALSE]
 
 # ---- 10) Fixed color scale like other figures (e.g. -6 to -3) ----
 FIX_MIN <- 0
 FIX_MAX <- 30  # prova 30 o 40 o 50
 
 bk <- seq(FIX_MIN, FIX_MAX, length.out = 101)
 cols <- colorRampPalette(c("white", "#FEE5D9", "#FC9272", "#DE2D26", "#A50F15"))(100)
 
 # cap SOLO per visualizzazione (opzionale ma utile)
 X_show <- pmin(pmax(X_pct, FIX_MIN), FIX_MAX)
 
 legend_breaks <- c(0,5,10,20,30)
 legend_labels <- paste0(legend_breaks, "%")
 
 pheatmap::pheatmap(
   X_show,
   color = cols, breaks = bk,
   cluster_rows = TRUE,
   cluster_cols = TRUE,
   annotation_col = ann_col[col_order, , drop=FALSE],
   annotation_colors = ann_colors,
   border_color = "grey80",
   legend_breaks = legend_breaks,
   legend_labels = legend_labels,
   main = "Top antibiotic classes (merged reps, % composition)",
   filename = sprintf("figures/%s_PERCENT.png", OUT_PREFIX),
   width = w, height = h
 )
 
 # sizes (avoid squished)
 w <- max(10, 0.35 * ncol(X_log) + 3)
 h <- max(7,  0.35 * nrow(X_log) + 3)
 
 # ---- 11) Plot heatmap ----
 pheatmap::pheatmap(
   X_log_cap,
   color = cols, breaks = bk,
   cluster_rows = TRUE,
   cluster_cols = TRUE,
   annotation_col = ann_col[col_order, , drop=FALSE],
   annotation_colors = ann_colors,
   border_color = "grey80",
   fontsize_row = 10,
   fontsize_col = 10,
   treeheight_col = 45,
   treeheight_row = 45,
   main = "Top antibiotic classes (merged reps, log10 relative abundance)",
   filename = sprintf("figures/%s.png", OUT_PREFIX),
   width = w, height = h
 )

 
 
 cat("Saved:\n",
     "- figures/", OUT_PREFIX, ".png\n",
     "- tables/", OUT_PREFIX, "_matrix_relative.tsv\n",
     "- tables/", OUT_PREFIX, "_matrix_log10_relative.tsv\n", sep="")
 
 library(dplyr)
 library(tibble)
 library(readr)
 library(tidyr)
 library(stringr)
 library(vegan)
 
 # -----------------------
 # INPUT
 # -----------------------
 mat_file  <- "ARG_per16S_T1T2_renamed.tsv"   # genes x samples (prima col = ARG_ID)
 meta_file <- "metadata3.txt"                 # replicate-level; deve contenere Sample + Valley
 
 out_dir <- "tables"
 dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
 
 USE_RELATIVE <- TRUE   # TRUE consigliato (composizione); FALSE = load assoluto per16S
 
 # -----------------------
 # read matrix robust (genes x samples)
 # -----------------------
 read_matrix_tsv_robust <- function(path) {
   x <- read.delim(path, header=TRUE, sep="\t",
                   check.names=FALSE, fill=TRUE,
                   quote="", comment.char="",
                   stringsAsFactors=FALSE)
   idcol <- names(x)[1]
   ids <- x[[idcol]]
   keep <- !(is.na(ids) | trimws(ids)=="")
   x <- x[keep,,drop=FALSE]
   rownames(x) <- make.unique(as.character(x[[idcol]]))
   x[[idcol]] <- NULL
   x[] <- lapply(x, function(v) suppressWarnings(as.numeric(v)))
   m <- as.matrix(x); m[is.na(m)] <- 0
   storage.mode(m) <- "numeric"
   m
 }
 
 m <- read_matrix_tsv_robust(mat_file)   # genes x samples
 meta <- read_tsv(meta_file, show_col_types = FALSE)
 
 # -----------------------
 # 1) Find Sample column + create Group
 # -----------------------
 if(!("Sample" %in% names(meta))){
   cand <- names(meta)[tolower(names(meta)) %in% c("sample","id","group")]
   if(length(cand)==0) stop("metadata3.txt: non trovo Sample/Group.")
   meta <- meta %>% rename(Sample = all_of(cand[1]))
 }
 
 meta <- meta %>% mutate(
   Sample = as.character(Sample),
   # BFJ1 -> BFJ; se già BFJ, resta BFJ
   Group  = ifelse(grepl("\\d+$", Sample), sub("\\d+$","", Sample), Sample)
 )
 
 # Valley column (case-insensitive)
 if(!("Valley" %in% names(meta))){
   candv <- names(meta)[tolower(names(meta)) %in% c("valley","vallata","valle")]
   if(length(candv)==0) stop("metadata3.txt: manca colonna Valley (o valley/valle).")
   meta <- meta %>% rename(Valley = all_of(candv[1]))
 }
 
 # ensure Season/Type exist (derive from Group if missing)
 if(!("Season" %in% names(meta))) meta$Season <- substr(meta$Group, 3, 3)  # J/S
 if(!("Base"   %in% names(meta))) meta$Base   <- substr(meta$Group, 1, 2)  # BF/BR...
 if(!("Type"   %in% names(meta))) meta$Type   <- ifelse(substr(meta$Base, 2, 2) == "R", "RG", "non-RG")
 
 # -----------------------
 # 2) Collapse metadata to Group (1 row per Group)
 # -----------------------
 meta_g <- meta %>%
   group_by(Group) %>%
   summarise(
     Valley = first(Valley),
     Season = first(Season),
     Type   = first(Type),
     Base   = first(Base),
     n_reps = n(),
     .groups="drop"
   ) %>%
   mutate(
     Group  = as.character(Group),
     Valley = factor(Valley),
     Season = factor(Season, levels=c("J","S")),
     Type   = factor(Type),
     Base   = factor(Base)
   )
 
 # -----------------------
 # 3) Align samples + merge replicates to Group (median)
 # -----------------------
 common_samples <- intersect(colnames(m), meta$Sample)
 if(length(common_samples) < 4) stop("Pochi campioni in comune tra matrice ARG e metadata3.txt.")
 
 m <- m[, common_samples, drop=FALSE]
 meta_s <- meta %>% filter(Sample %in% common_samples) %>% select(Sample, Group)
 
 # long -> median per Gene x Group -> wide
 arg_long <- as.data.frame(m) %>%
   rownames_to_column("Gene") %>%
   pivot_longer(-Gene, names_to="Sample", values_to="Abundance") %>%
   left_join(meta_s, by="Sample") %>%
   group_by(Gene, Group) %>%
   summarise(Abundance = median(Abundance, na.rm=TRUE), .groups="drop")
 
 arg_group_wide <- arg_long %>%
   pivot_wider(names_from = Gene, values_from = Abundance, values_fill = 0)
 
 X <- as.matrix(arg_group_wide[,-1])           # rows = Group, cols = genes
 rownames(X) <- arg_group_wide$Group
 storage.mode(X) <- "numeric"
 
 # align groups with meta_g
 common_groups <- intersect(rownames(X), meta_g$Group)
 if(length(common_groups) < 4) stop("Pochi gruppi in comune dopo group-merge. Controlla Group naming.")
 X <- X[common_groups, , drop=FALSE]
 meta2 <- meta_g %>% filter(Group %in% common_groups)
 meta2 <- meta2[match(common_groups, meta2$Group), ]
 stopifnot(all(meta2$Group == common_groups))
 
 # optional: relative within group
 if(USE_RELATIVE){
   rs <- rowSums(X)
   X <- sweep(X, 1, ifelse(rs>0, rs, 1), "/")
 }
 
 # Bray distance
 d <- vegdist(X, method="bray")
 
 # -----------------------
 # 4) Drop factors with <2 levels (avoid contrasts error)
 # -----------------------
 nlev <- c(
   Valley = nlevels(droplevels(meta2$Valley)),
   Type   = nlevels(droplevels(meta2$Type)),
   Season = nlevels(droplevels(meta2$Season))
 )
 print(nlev)
 
 keep_terms <- names(nlev)[nlev >= 2]
 if(length(keep_terms) == 0) stop("Nessun fattore ha >=2 livelli dopo allineamento. Controlla Valley/Season/Type.")
 
 form <- as.formula(paste("d ~", paste(keep_terms, collapse = " + ")))
 
 set.seed(1)
 perm_marg  <- adonis2(form, data=meta2, permutations=999, by="margin")
 perm_terms <- adonis2(form, data=meta2, permutations=999, by="terms")
 
 write_tsv(as.data.frame(perm_marg)  %>% rownames_to_column("Term"),
           file.path(out_dir, "PERMANOVA_ARG_GroupMerged_Valley_Type_Season_MARGIN.tsv"))
 write_tsv(as.data.frame(perm_terms) %>% rownames_to_column("Term"),
           file.path(out_dir, "PERMANOVA_ARG_GroupMerged_Valley_Type_Season_TERMS.tsv"))
 
 # -----------------------
 # 5) PERMDISP for factors with >=2 levels
 # -----------------------
 disp_list <- list()
 
 if(nlev["Valley"] >= 2){
   bd <- betadisper(d, meta2$Valley)
   a  <- anova(bd)
   disp_list$Valley <- tibble(Factor="Valley", F=a$`F value`[1], p=a$`Pr(>F)`[1])
 }
 if(nlev["Type"] >= 2){
   bd <- betadisper(d, meta2$Type)
   a  <- anova(bd)
   disp_list$Type <- tibble(Factor="Type", F=a$`F value`[1], p=a$`Pr(>F)`[1])
 }
 if(nlev["Season"] >= 2){
   bd <- betadisper(d, meta2$Season)
   a  <- anova(bd)
   disp_list$Season <- tibble(Factor="Season", F=a$`F value`[1], p=a$`Pr(>F)`[1])
 }
 
 disp_tbl <- bind_rows(disp_list) %>%
   mutate(p_adj = p.adjust(p, "BH"))
 
 write_tsv(disp_tbl, file.path(out_dir, "PERMDISP_ARG_GroupMerged_Valley_Type_Season.tsv"))
 
 cat("\nSaved:\n",
     "- PERMANOVA_ARG_GroupMerged_Valley_Type_Season_MARGIN.tsv\n",
     "- PERMANOVA_ARG_GroupMerged_Valley_Type_Season_TERMS.tsv\n",
     "- PERMDISP_ARG_GroupMerged_Valley_Type_Season.tsv\n", sep="")
 
 library(readr)
 library(dplyr)
 library(tidyr)
 library(stringr)
 library(tibble)
 
 arg_matrix_file <- "ARG_per16S_T1T2_renamed.tsv"
 arg2class_file  <- "ARG_annotation_RGI_T1T2.tsv"
 
 clean_gene <- function(x){
   x <- as.character(x)
   x <- sub("^.*\\|", "", x, perl = TRUE)     # dopo ultimo |
   x <- sub("^[0-9]+", "", x, perl = TRUE)    # rimuove numeri iniziali
   x <- sub("^[-_ :]+", "", x, perl = TRUE)   # rimuove separatori iniziali
   stringr::str_trim(x)
 }
 
 # --- 1) read mapping (RGI) ---
 ann <- readr::read_tsv(arg2class_file, show_col_types = FALSE) %>%
   mutate(
     ARG_ID  = as.character(ARG_ID),
     BestHit = as.character(BestHit),
     DrugClass = stringr::str_trim(as.character(DrugClass))
   ) %>%
   tidyr::separate_rows(DrugClass, sep="\\s*[,;/|]+\\s*") %>%  # se più classi in una cella
   filter(!is.na(DrugClass), DrugClass != "") %>%
   mutate(
     Gene_from_ARGID   = clean_gene(ARG_ID),
     Gene_from_BestHit = clean_gene(BestHit)
   )
 
 n_classes_mapping <- n_distinct(ann$DrugClass)
 
 # --- 2) read ARG matrix and list observed genes (presenti >0 in almeno un campione) ---
 arg_raw <- readr::read_tsv(arg_matrix_file, show_col_types = FALSE)
 arg_id_col <- names(arg_raw)[1]
 
 ARG_mat <- arg_raw %>%
   tibble::column_to_rownames(arg_id_col) %>%
   as.data.frame()
 
 ARG_mat[] <- lapply(ARG_mat, as.numeric)
 ARG_mat <- as.matrix(ARG_mat)
 
 # geni osservati (almeno un campione con abbondanza >0)
 genes_observed <- rownames(ARG_mat)[rowSums(ARG_mat > 0, na.rm = TRUE) > 0] %>%
   clean_gene() %>%
   unique()
 
 # mapping Gene -> DrugClass usando sia ARG_ID che BestHit
 map_gene_class <- ann %>%
   pivot_longer(cols = c(Gene_from_ARGID, Gene_from_BestHit),
                names_to = "src", values_to = "Gene") %>%
   filter(!is.na(Gene), Gene != "") %>%
   distinct(Gene, DrugClass)
 
 # classi effettivamente osservate nei campioni
 classes_observed <- map_gene_class %>%
   filter(Gene %in% genes_observed) %>%
   distinct(DrugClass) %>%
   pull(DrugClass)
 
 cat("Unique antibiotic classes in RGI annotation file:", n_classes_mapping, "\n")
 cat("Unique antibiotic classes observed in samples:", length(classes_observed), "\n\n")
 cat("Observed classes:\n")
 print(sort(classes_observed))
 
 # --- 3) opzionale: quante classi per campione (presenza/assenza) ---
 ARG_long <- arg_raw %>%
   pivot_longer(-all_of(arg_id_col), names_to="Sample", values_to="Abund") %>%
   mutate(
     Abund = as.numeric(Abund),
     Gene = clean_gene(.data[[arg_id_col]])
   ) %>%
   filter(Abund > 0)
 
 classes_per_sample <- ARG_long %>%
   left_join(map_gene_class, by = c("Gene" = "Gene")) %>%
   filter(!is.na(DrugClass)) %>%
   distinct(Sample, DrugClass) %>%
   count(Sample, name="n_drug_classes") %>%
   arrange(desc(n_drug_classes))
 
 cat("\nClasses per sample (top 10):\n")
 print(head(classes_per_sample, 10))
 
 cat("\nSummary classes per sample:\n")
 print(summary(classes_per_sample$n_drug_classes))
 cat("Median:", median(classes_per_sample$n_drug_classes, na.rm=TRUE), "\n")
 
 ###############################################################################
 # PCoA of ARG genes + envfit
 # - Input ARG matrix: genes x samples (ARG per 16S)
 # - Colors = Base
 # - Shapes = Season
 # - envfit on numeric variables from metadata3.txt
 # - Default: replicates NOT merged
 #   Set MERGE_REPS <- TRUE to merge replicates by median within Base x Season
 ###############################################################################
 
 #### 0) Packages --------------------------------------------------------------
 pkgs <- c("dplyr","readr","tidyr","stringr","tibble","ggplot2","vegan")
 to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
 if(length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
 invisible(lapply(pkgs, library, character.only = TRUE))
 
 dir.create("figures", showWarnings = FALSE, recursive = TRUE)
 dir.create("tables",  showWarnings = FALSE, recursive = TRUE)
 
 #### 1) Settings --------------------------------------------------------------
 mat_file   <- "ARG_per16S_T1T2_renamed.tsv"   # genes x samples
 env_file   <- "metadata3.txt"                 # must contain Sample column
 
 MERGE_REPS    <- FALSE   # TRUE = merge replicates by median within Base x Season
 TOP_K_ARROWS  <- 10
 P_ADJ_CUTOFF  <- 0.05
 N_PERM        <- 999
 
 # if you want to force specific env vars, write them here, e.g. c("T","ph","Cd","Pb")
 KEEP_ENV_VARS <- NULL
 
 # numeric columns to exclude automatically from envfit if present
 EXCLUDE_NUMERIC <- c("Replicate", "Rep", "n_reps")
 
 base_cols_map <- c(
   BF = "#1b9e77", BR = "#a6dba0",
   CF = "#d95f02", CR = "#fdb863",
   PF = "#7570b3", PR = "#b2abd2",
   SF = "#1f9ac2", SR = "#a6dce7",
   VF = "#e7298a", VR = "#f2b2d4"
 )
 
 shape_map <- c(July = 16, September = 17)
 
 #### 2) Read ARG matrix and build long table ----------------------------------
 arg_wide <- read_tsv(mat_file, show_col_types = FALSE)
 gene_col <- names(arg_wide)[1]
 sample_cols <- setdiff(names(arg_wide), gene_col)
 
 arg_long <- arg_wide %>%
   pivot_longer(cols = all_of(sample_cols), names_to = "Sample", values_to = "Abundance") %>%
   mutate(
     Abundance = as.numeric(Abundance),
     Abundance = ifelse(is.na(Abundance), 0, Abundance),
     GeneID    = as.character(.data[[gene_col]])
   ) %>%
   select(GeneID, Sample, Abundance)
 
 #### 3) Metadata from sample names --------------------------------------------
 meta_parse <- arg_long %>%
   distinct(Sample) %>%
   mutate(
     Base        = str_sub(Sample, 1, 2),
     Site        = str_sub(Sample, 1, 1),
     Season_code = str_sub(Sample, 3, 3),
     Season      = recode(Season_code, "J" = "July", "S" = "September", .default = Season_code),
     Replicate   = str_extract(Sample, "\\d+$"),
     Type        = ifelse(str_sub(Base, 2, 2) == "R", "RG", "non-RG"),
     Group       = paste0(Base, Season_code)
   )
 
 arg_long <- arg_long %>% left_join(meta_parse, by = "Sample")
 
 #### 4) Build ordination matrix ------------------------------------------------
 if (!MERGE_REPS) {
   
   # sample-level matrix (replicates separate)
   X_wide <- arg_long %>%
     select(Sample, GeneID, Abundance) %>%
     pivot_wider(names_from = GeneID, values_from = Abundance, values_fill = 0)
   
   X <- X_wide %>% select(-Sample) %>% as.data.frame()
   rownames(X) <- X_wide$Sample
   X <- as.matrix(X)
   storage.mode(X) <- "numeric"
   
   meta_ord <- meta_parse %>%
     filter(Sample %in% rownames(X)) %>%
     arrange(match(Sample, rownames(X)))
   
   stopifnot(all(meta_ord$Sample == rownames(X)))
   
 } else {
   
   # merge replicates by median within Base x Season
   arg_group <- arg_long %>%
     group_by(Group, Base, Site, Season, Season_code, Type, GeneID) %>%
     summarise(Abundance = median(Abundance, na.rm = TRUE), .groups = "drop") %>%
     mutate(Abundance = ifelse(is.na(Abundance), 0, Abundance))
   
   X_wide <- arg_group %>%
     select(Group, GeneID, Abundance) %>%
     pivot_wider(names_from = GeneID, values_from = Abundance, values_fill = 0)
   
   X <- X_wide %>% select(-Group) %>% as.data.frame()
   rownames(X) <- X_wide$Group
   X <- as.matrix(X)
   storage.mode(X) <- "numeric"
   
   meta_ord <- arg_group %>%
     distinct(Group, Base, Site, Season, Season_code, Type) %>%
     arrange(match(Group, rownames(X)))
   
   stopifnot(all(meta_ord$Group == rownames(X)))
 }
 
 # remove all-zero rows
 keep <- rowSums(X) > 0
 X <- X[keep, , drop = FALSE]
 
 if (!MERGE_REPS) {
   meta_ord <- meta_ord[keep, , drop = FALSE]
   row_id <- meta_ord$Sample
 } else {
   meta_ord <- meta_ord[keep, , drop = FALSE]
   row_id <- meta_ord$Group
 }
 
 stopifnot(all(row_id == rownames(X)))
 
 #### 5) Relative abundance + Bray-Curtis + PCoA -------------------------------
 X_rel  <- vegan::decostand(X, method = "total")
 d_bray <- vegan::vegdist(X_rel, method = "bray")
 
 pcoa <- cmdscale(d_bray, k = 2, eig = TRUE)
 
 coords <- as.data.frame(pcoa$points)
 colnames(coords) <- c("PCoA1", "PCoA2")
 coords$ID <- rownames(coords)
 
 eig_pos <- pcoa$eig[pcoa$eig > 0]
 var_expl <- round(100 * eig_pos / sum(eig_pos), 1)
 var1 <- ifelse(length(var_expl) >= 1, var_expl[1], NA)
 var2 <- ifelse(length(var_expl) >= 2, var_expl[2], NA)
 
 if (!MERGE_REPS) {
   plot_df <- coords %>%
     left_join(meta_ord %>% rename(ID = Sample), by = "ID")
 } else {
   plot_df <- coords %>%
     left_join(meta_ord %>% rename(ID = Group), by = "ID")
 }
 
 #### 6) Read env data ---------------------------------------------------------
 env_raw <- read_tsv(env_file, show_col_types = FALSE)
 
 if (!("Sample" %in% names(env_raw))) {
   cand <- names(env_raw)[tolower(names(env_raw)) %in% c("sample", "id")]
   if (length(cand) == 0) stop("env_file must contain a Sample column.")
   env_raw <- env_raw %>% rename(Sample = all_of(cand[1]))
 }
 
 env_raw <- env_raw %>%
   mutate(Sample = as.character(Sample))
 
 # derive grouping info from Sample if needed
 env_raw <- env_raw %>%
   mutate(
     Base        = ifelse("Base" %in% names(.), as.character(.data$Base), str_sub(Sample, 1, 2)),
     Site        = ifelse("Site" %in% names(.), as.character(.data$Site), str_sub(Sample, 1, 1)),
     Season_code = ifelse("Season_code" %in% names(.), as.character(.data$Season_code), str_sub(Sample, 3, 3)),
     Season      = ifelse("Season" %in% names(.), as.character(.data$Season),
                          recode(Season_code, "J" = "July", "S" = "September", .default = Season_code)),
     Type        = ifelse("Type" %in% names(.), as.character(.data$Type),
                          ifelse(str_sub(Base, 2, 2) == "R", "RG", "non-RG")),
     Group       = paste0(Base, Season_code)
   )
 
 if (!MERGE_REPS) {
   
   env2 <- env_raw %>%
     filter(Sample %in% row_id) %>%
     arrange(match(Sample, row_id))
   
   stopifnot(all(env2$Sample == row_id))
   
   if (is.null(KEEP_ENV_VARS)) {
     num_cols <- names(env2)[sapply(env2, is.numeric)]
     env_vars <- setdiff(num_cols, EXCLUDE_NUMERIC)
   } else {
     env_vars <- KEEP_ENV_VARS
   }
   
   chem_mat <- env2 %>%
     select(all_of(env_vars)) %>%
     mutate(across(everything(), as.numeric)) %>%
     as.data.frame()
   
 } else {
   
   # merge numeric env variables by Group (median)
   if (is.null(KEEP_ENV_VARS)) {
     num_cols <- names(env_raw)[sapply(env_raw, is.numeric)]
     env_vars <- setdiff(num_cols, EXCLUDE_NUMERIC)
   } else {
     env_vars <- KEEP_ENV_VARS
   }
   
   env2 <- env_raw %>%
     filter(Group %in% row_id) %>%
     group_by(Group, Base, Site, Season, Type) %>%
     summarise(across(all_of(env_vars), ~ median(as.numeric(.x), na.rm = TRUE)), .groups = "drop") %>%
     arrange(match(Group, row_id))
   
   stopifnot(all(env2$Group == row_id))
   
   chem_mat <- env2 %>%
     select(all_of(env_vars)) %>%
     mutate(across(everything(), as.numeric)) %>%
     as.data.frame()
 }
 
 # keep only finite, variable columns
 ok_cols <- sapply(chem_mat, function(x) {
   x2 <- x[is.finite(x)]
   length(x2) > 2 && sd(x2, na.rm = TRUE) > 0
 })
 
 chem_mat <- chem_mat[, ok_cols, drop = FALSE]
 
 if (ncol(chem_mat) == 0) {
   stop("No usable numeric environmental variables found for envfit.")
 }
 
 #### 7) envfit ---------------------------------------------------------------
 fit <- vegan::envfit(pcoa$points, chem_mat, permutations = N_PERM)
 
 vec <- as.data.frame(scores(fit, display = "vectors"))
 vec$Variable <- rownames(vec)
 
 r2 <- fit$vectors$r
 p  <- fit$vectors$pvals
 
 env_tbl <- tibble(
   Variable = names(r2),
   r2       = as.numeric(r2),
   p_value  = as.numeric(p)
 ) %>%
   mutate(
     p_adj = p.adjust(p_value, method = "BH")
   ) %>%
   arrange(p_adj, desc(r2))
 
 write_tsv(env_tbl,
           ifelse(MERGE_REPS,
                  "tables/envfit_ARGgenes_PCoA_GroupMerged.tsv",
                  "tables/envfit_ARGgenes_PCoA_sampleLevel.tsv"))
 
 # select arrows to draw
 arrow_keep <- env_tbl %>%
   filter(is.finite(p_adj), p_adj < P_ADJ_CUTOFF) %>%
   slice_max(order_by = r2, n = TOP_K_ARROWS, with_ties = FALSE)
 
 if (nrow(arrow_keep) == 0) {
   arrow_keep <- env_tbl %>%
     slice_max(order_by = r2, n = TOP_K_ARROWS, with_ties = FALSE)
 }
 
 vec2 <- vec %>%
   select(Variable, 1, 2)
 
 names(vec2)[2:3] <- c("x", "y")
 
 arrow_df <- vec2 %>%
   inner_join(arrow_keep, by = "Variable")
 
 # scale arrows to ordination space
 xrange <- diff(range(plot_df$PCoA1, na.rm = TRUE))
 yrange <- diff(range(plot_df$PCoA2, na.rm = TRUE))
 mult <- 0.6 * min(xrange, yrange)
 
 arrow_df <- arrow_df %>%
   mutate(
     strength = sqrt(r2),
     xend = x * mult * strength,
     yend = y * mult * strength
   )
 
 # label placement
 arrow_df <- arrow_df %>%
   mutate(
     alen = sqrt(xend^2 + yend^2),
     off = pmin(pmax(0.015 * min(xrange, yrange), 0.12 * alen), 0.04 * min(xrange, yrange)),
     perp_x = -yend,
     perp_y =  xend,
     perp_len = sqrt(perp_x^2 + perp_y^2),
     perp_xu = ifelse(perp_len > 0, perp_x / perp_len, 0),
     perp_yu = ifelse(perp_len > 0, perp_y / perp_len, 0),
     label_x = xend + off * perp_xu,
     label_y = yend + off * perp_yu,
     hjust = ifelse(label_x >= 0, -0.1, 1.1),
     vjust = ifelse(label_y >= 0, -0.2, 1.2)
   )
 
 #### 8) Optional PERMANOVA / PERMDISP ----------------------------------------
 # sample-level: can use strata = Site
 if (!MERGE_REPS) {
   perm1 <- vegan::adonis2(d_bray ~ Season + Type, data = meta_ord, permutations = N_PERM, strata = meta_ord$Site)
   perm2 <- vegan::adonis2(d_bray ~ Site + Season + Type, data = meta_ord, permutations = N_PERM)
   
   write.table(perm1,
               "tables/PERMANOVA_ARGgenes_Season_Type_strataSite.tsv",
               sep = "\t", quote = FALSE, col.names = NA)
   write.table(perm2,
               "tables/PERMANOVA_ARGgenes_Site_Season_Type.tsv",
               sep = "\t", quote = FALSE, col.names = NA)
   
   bd_season <- vegan::betadisper(d_bray, meta_ord$Season)
   bd_type   <- vegan::betadisper(d_bray, meta_ord$Type)
   
   write.table(anova(bd_season),
               "tables/betadisper_ARGgenes_Season.tsv",
               sep = "\t", quote = FALSE, col.names = NA)
   write.table(anova(bd_type),
               "tables/betadisper_ARGgenes_Type.tsv",
               sep = "\t", quote = FALSE, col.names = NA)
   
 } else {
   perm2 <- vegan::adonis2(d_bray ~ Site + Season + Type, data = meta_ord, permutations = N_PERM)
   
   write.table(perm2,
               "tables/PERMANOVA_ARGgenes_GroupMerged_Site_Season_Type.tsv",
               sep = "\t", quote = FALSE, col.names = NA)
   
   bd_season <- vegan::betadisper(d_bray, meta_ord$Season)
   bd_type   <- vegan::betadisper(d_bray, meta_ord$Type)
   
   write.table(anova(bd_season),
               "tables/betadisper_ARGgenes_GroupMerged_Season.tsv",
               sep = "\t", quote = FALSE, col.names = NA)
   write.table(anova(bd_type),
               "tables/betadisper_ARGgenes_GroupMerged_Type.tsv",
               sep = "\t", quote = FALSE, col.names = NA)
 }
 
 #### 9) Plot ------------------------------------------------------------------
 pkgs <- c("dplyr","readr","tidyr","stringr","tibble","ggplot2","vegan","ggrepel")
  
 plot_df$Base   <- factor(plot_df$Base, levels = names(base_cols_map))
 plot_df$Season <- factor(plot_df$Season, levels = c("July", "September"))
 labs(
   x = paste0("PCoA1 (", var1, "%)"),
   y = paste0("PCoA2 (", var2, "%)"),
   color = "Site",
   shape = "Season",
   title = ifelse(MERGE_REPS,
                  "PCoA of ARG genes (merged replicates, Bray-Curtis)",
                  "PCoA of ARG genes (sample level, Bray-Curtis)")
 )

 missing_base <- setdiff(unique(as.character(plot_df$Base)), names(base_cols_map))
 if (length(missing_base) > 0) {
   warning("These Base levels are not in base_cols_map: ",
           paste(missing_base, collapse = ", "))
 }
 # separa le label senza disegnare segmenti aggiuntivi
 arrow_df <- arrow_df %>%
   mutate(
     strength = sqrt(r2),
     xend = x * mult * strength,
     yend = y * mult * strength,
     angle = atan2(yend, xend)
   ) %>%
   arrange(angle) %>%
   mutate(
     # offset alternato per frecce vicine
     off_id = c(-1, 1, -2, 2, -3, 3, -4, 4)[seq_len(n())],
     off = 0.04 * min(xrange, yrange) * off_id,
     
     # direzione perpendicolare alla freccia
     perp_x = -sin(angle),
     perp_y =  cos(angle),
     
     # label poco oltre la punta + offset perpendicolare
     label_x = xend * 1.08 + off * perp_x
     label_y = yend * 1.08 + off * perp_y
   )
 
 arrow_df <- arrow_df %>%
   mutate(
     dx = 0,
     dy = 0
   )
 label_adjust <- tibble::tribble(
   ~Variable, ~dx,   ~dy,
   "Cd",      -0.02,  0.01,
   "Ce",       0.03, -0.01,
   "T",        0.01,-0.01
 )
 
 arrow_df <- arrow_df %>%
   select(-any_of(c("dx","dy","dx.x","dx.y","dy.x","dy.y","label_x","label_y"))) %>%
   left_join(label_adjust, by = "Variable") %>%
   mutate(
     dx = dplyr::coalesce(dx, 0),
     dy = dplyr::coalesce(dy, 0),
     label_x = xend + dx,
     label_y = yend + dy
   )
 
 plot_df$Base   <- factor(plot_df$Base, levels = names(base_cols_map))
 plot_df$Season <- factor(plot_df$Season, levels = c("July", "September"))
 
 p <- ggplot(plot_df, aes(PCoA1, PCoA2, color = Base, shape = Season)) +
   geom_point(size = 3, alpha = 0.9) +
   scale_color_manual(values = base_cols_map, drop = FALSE) +
   scale_shape_manual(values = shape_map, drop = FALSE) +
   geom_segment(
     data = arrow_df,
     aes(x = 0, y = 0, xend = xend, yend = yend),
     inherit.aes = FALSE,
     linewidth = 0.6,
     arrow = arrow(length = grid::unit(0.22, "cm"))
   ) +
   geom_text(
     data = arrow_df,
     aes(x = label_x, y = label_y, label = Variable),
     inherit.aes = FALSE,
     size = 3
   ) +
   coord_cartesian(clip = "off") +
   theme_bw() +
   labs(
     x = paste0("PCoA1 (", var1, "%)"),
     y = paste0("PCoA2 (", var2, "%)"),
     color = "Site",
     shape = "Season",
     title = ifelse(MERGE_REPS,
                    "PCoA of ARG genes (merged replicates, Bray-Curtis)",
                    "PCoA of ARG genes (sample level, Bray-Curtis)")
   )
 
 outfile_plot <- ifelse(MERGE_REPS,
                        "figures/PCoA_ARGgenes_Bray_colorBase_shapeSeason_envfit_GroupMerged.png",
                        "figures/PCoA_ARGgenes_Bray_colorBase_shapeSeason_envfit_sampleLevel.png")
 
 ggsave(outfile_plot, p, width = 8, height = 6, dpi = 300)
 
 cat("Saved:\n")
 cat("-", outfile_plot, "\n")
 cat("-", ifelse(MERGE_REPS,
                 "tables/envfit_ARGgenes_PCoA_GroupMerged.tsv",
                 "tables/envfit_ARGgenes_PCoA_sampleLevel.tsv"), "\n") 
 
 
 ###############################################################################
 # ARG drug-class summary on replicate-merged profiles
 # - Merge replicates by median within sampling site and season
 # - Summarise drug classes:
 #     * number of observed classes
 #     * mean / median / max relative abundance
 #     * prevalence across merged groups
 #     * number of classes per merged group
 # - Save human-readable TXT report
 # - Save matrices/tables ready for Bray-Curtis later
 ###############################################################################
 
 #### 0) Packages --------------------------------------------------------------
 pkgs <- c("readr","dplyr","tibble","tidyr","stringr")
 to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
 if(length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
 invisible(lapply(pkgs, library, character.only = TRUE))
 
 dir.create("tables",  showWarnings = FALSE, recursive = TRUE)
 dir.create("figures", showWarnings = FALSE, recursive = TRUE)
 
 #### 1) Inputs ----------------------------------------------------------------
 arg_matrix_file <- "ARG_per16S_T1T2_renamed.tsv"   # genes x samples
 arg2class_file  <- "ARG_annotation_RGI_T1T2.tsv"   # must contain ARG_ID, BestHit, DrugClass
 
 # If TRUE: if one gene maps to multiple drug classes, abundance is split equally
 # If FALSE: abundance is counted in every mapped class (same logic as before)
 SPLIT_MULTI_CLASS <- FALSE
 
 report_file <- "tables/ARG_drugclass_summary_report_merged.txt"
 
 #### 2) Helper ----------------------------------------------------------------
 clean_gene <- function(x){
   x <- as.character(x)
   x <- sub("^.*\\|", "", x, perl = TRUE)   # keep text after last |
   x <- sub("^[0-9]+", "", x, perl = TRUE)  # remove leading numbers
   x <- sub("^[-_ :]+", "", x, perl = TRUE) # remove leading separators
   stringr::str_trim(x)
 }
 
 #### 3) Read ARG matrix -> long -----------------------------------------------
 arg_raw <- readr::read_tsv(arg_matrix_file, show_col_types = FALSE)
 arg_id_col <- names(arg_raw)[1]
 
 ARG_long <- arg_raw %>%
   pivot_longer(-all_of(arg_id_col), names_to = "Sample", values_to = "Abundance") %>%
   mutate(
     Abundance       = suppressWarnings(as.numeric(Abundance)),
     Abundance       = ifelse(is.na(Abundance), 0, Abundance),
     ARG_ID          = as.character(.data[[arg_id_col]]),
     ARG_gene_clean  = clean_gene(ARG_ID)
   ) %>%
   select(Sample, ARG_ID, ARG_gene_clean, Abundance)
 
 #### 4) Read ARG -> DrugClass mapping -----------------------------------------
 ann <- readr::read_tsv(arg2class_file, show_col_types = FALSE)
 names(ann) <- trimws(names(ann))
 
 stopifnot(all(c("ARG_ID","BestHit","DrugClass") %in% names(ann)))
 
 ab_map <- ann %>%
   transmute(
     ARG_ID    = as.character(ARG_ID),
     BestHit   = as.character(BestHit),
     DrugClass = stringr::str_trim(as.character(DrugClass))
   ) %>%
   separate_rows(DrugClass, sep = "\\s*[,;/|]+\\s*") %>%
   filter(!is.na(DrugClass), DrugClass != "") %>%
   mutate(
     Gene_from_ARGID   = clean_gene(ARG_ID),
     Gene_from_BestHit = clean_gene(BestHit)
   ) %>%
   pivot_longer(
     cols = c(Gene_from_ARGID, Gene_from_BestHit),
     names_to = "src", values_to = "Gene"
   ) %>%
   mutate(Gene = stringr::str_trim(Gene)) %>%
   filter(!is.na(Gene), Gene != "") %>%
   distinct(Gene, DrugClass)
 
 # how many classes per gene? (needed only if SPLIT_MULTI_CLASS = TRUE)
 gene_nclass <- ab_map %>%
   distinct(Gene, DrugClass) %>%
   count(Gene, name = "n_classes")
 
 #### 5) Metadata from sample names --------------------------------------------
 meta <- tibble(Sample = unique(ARG_long$Sample)) %>%
   mutate(
     SamplingSite = str_sub(Sample, 1, 2),   # BF, BR, CF, CR...
     Location     = str_sub(Sample, 1, 1),   # B, C, P, S, V
     Season_code  = str_sub(Sample, 3, 3),   # J, S
     Season       = recode(Season_code, "J" = "July", "S" = "September", .default = Season_code),
     Replicate    = str_extract(Sample, "\\d+$"),
     Type         = ifelse(str_sub(SamplingSite, 2, 2) == "R", "RG", "non-RG"),
     Group        = paste(SamplingSite, Season, sep = "_")
   )
 
 #### 6) Join mapping + optionally split multi-class abundance ------------------
 ARG_class_long <- ARG_long %>%
   left_join(ab_map, by = c("ARG_gene_clean" = "Gene"), relationship = "many-to-many") %>%
   mutate(DrugClass = ifelse(is.na(DrugClass), "Unclassified", DrugClass)) %>%
   left_join(meta, by = "Sample")
 
 if (SPLIT_MULTI_CLASS) {
   ARG_class_long <- ARG_class_long %>%
     left_join(gene_nclass, by = c("ARG_gene_clean" = "Gene")) %>%
     mutate(
       n_classes = ifelse(is.na(n_classes) | n_classes < 1, 1, n_classes),
       Abundance_class = Abundance / n_classes
     )
 } else {
   ARG_class_long <- ARG_class_long %>%
     mutate(Abundance_class = Abundance)
 }
 
 #### 7) Sum per sample and merge replicates by median --------------------------
 # class abundance per original sample
 class_sample <- ARG_class_long %>%
   group_by(Sample, SamplingSite, Location, Season, Type, Group, DrugClass) %>%
   summarise(Class_abund = sum(Abundance_class, na.rm = TRUE), .groups = "drop")
 
 # replicate-merged by median within SamplingSite + Season + DrugClass
 class_group <- class_sample %>%
   group_by(Group, SamplingSite, Location, Season, Type, DrugClass) %>%
   summarise(Class_abund_med = median(Class_abund, na.rm = TRUE), .groups = "drop") %>%
   mutate(Class_abund_med = ifelse(is.na(Class_abund_med), 0, Class_abund_med))
 
 #### 8) Relative composition within merged group -------------------------------
 class_group <- class_group %>%
   group_by(Group) %>%
   mutate(
     Total_group = sum(Class_abund_med, na.rm = TRUE),
     Class_rel   = ifelse(Total_group > 0, Class_abund_med / Total_group, 0)
   ) %>%
   ungroup()
 
 #### 9) Summary per drug class ------------------------------------------------
 drugclass_summary <- class_group %>%
   group_by(DrugClass) %>%
   summarise(
     prevalence_groups_n = sum(Class_rel > 0, na.rm = TRUE),
     prevalence_groups_pct = 100 * mean(Class_rel > 0, na.rm = TRUE),
     mean_rel   = mean(Class_rel, na.rm = TRUE),
     median_rel = median(Class_rel, na.rm = TRUE),
     max_rel    = max(Class_rel, na.rm = TRUE),
     mean_pct   = 100 * mean_rel,
     median_pct = 100 * median_rel,
     max_pct    = 100 * max_rel,
     .groups = "drop"
   ) %>%
   arrange(desc(mean_rel))
 
 #### 10) Number of classes per merged group -----------------------------------
 classes_per_group <- class_group %>%
   group_by(Group, SamplingSite, Location, Season, Type) %>%
   summarise(
     n_drug_classes = sum(Class_rel > 0, na.rm = TRUE),
     total_class_abundance = sum(Class_abund_med, na.rm = TRUE),
     .groups = "drop"
   ) %>%
   arrange(desc(n_drug_classes))
 
 classes_per_season <- classes_per_group %>%
   group_by(Season) %>%
   summarise(
     n_groups = n(),
     median_n_classes = median(n_drug_classes, na.rm = TRUE),
     mean_n_classes   = mean(n_drug_classes, na.rm = TRUE),
     min_n_classes    = min(n_drug_classes, na.rm = TRUE),
     max_n_classes    = max(n_drug_classes, na.rm = TRUE),
     .groups = "drop"
   )
 
 #### 11) Overall counts -------------------------------------------------------
 n_classes_annotated <- ann %>%
   transmute(DrugClass = stringr::str_trim(as.character(DrugClass))) %>%
   separate_rows(DrugClass, sep = "\\s*[,;/|]+\\s*") %>%
   filter(!is.na(DrugClass), DrugClass != "") %>%
   summarise(n = n_distinct(DrugClass)) %>%
   pull(n)
 
 n_classes_observed_groups <- drugclass_summary %>%
   filter(prevalence_groups_n > 0, DrugClass != "Unclassified") %>%
   summarise(n = n()) %>%
   pull(n)
 
 n_groups_total <- n_distinct(class_group$Group)
 
 #### 12) Matrix ready for Bray-Curtis later -----------------------------------
 class_matrix_rel <- class_group %>%
   select(Group, DrugClass, Class_rel) %>%
   pivot_wider(names_from = DrugClass, values_from = Class_rel, values_fill = 0) %>%
   arrange(Group)
 
 #### 13) Save tables ----------------------------------------------------------
 write_tsv(drugclass_summary,
           "tables/ARG_drugclass_summary_mergedGroups.tsv")
 
 write_tsv(classes_per_group,
           "tables/ARG_nDrugClasses_perMergedGroup.tsv")
 
 write_tsv(classes_per_season,
           "tables/ARG_nDrugClasses_perMergedGroup_bySeason.tsv")
 
 write_tsv(class_group,
           "tables/ARG_drugclass_composition_mergedGroups_long.tsv")
 
 write_tsv(class_matrix_rel,
           "tables/ARG_drugclass_matrix_relative_mergedGroups.tsv")
 
 #### 14) Human-readable TXT report --------------------------------------------
 top10_classes <- drugclass_summary %>% slice_head(n = 10)
 
 lines <- c()
 lines <- c(lines, "ARG DRUG-CLASS SUMMARY (replicate-merged profiles)")
 lines <- c(lines, "=================================================")
 lines <- c(lines, "")
 lines <- c(lines, sprintf("Input ARG matrix: %s", arg_matrix_file))
 lines <- c(lines, sprintf("Input annotation file: %s", arg2class_file))
 lines <- c(lines, sprintf("Replicates merged by median within sampling site and season."))
 lines <- c(lines, sprintf("Multi-class handling: %s",
                           ifelse(SPLIT_MULTI_CLASS,
                                  "abundance split equally across classes",
                                  "abundance counted in each mapped class")))
 lines <- c(lines, "")
 lines <- c(lines, "OVERVIEW")
 lines <- c(lines, "--------")
 lines <- c(lines, sprintf("Number of annotated drug classes in the RGI file: %d", n_classes_annotated))
 lines <- c(lines, sprintf("Number of drug classes observed across replicate-merged groups: %d", n_classes_observed_groups))
 lines <- c(lines, sprintf("Number of replicate-merged groups (sampling site + season): %d", n_groups_total))
 lines <- c(lines, "")
 
 lines <- c(lines, "NUMBER OF DRUG CLASSES PER MERGED GROUP")
 lines <- c(lines, "--------------------------------------")
 lines <- c(lines, capture.output(print(classes_per_group, n = nrow(classes_per_group))))
 lines <- c(lines, "")
 lines <- c(lines, "Summary by season:")
 lines <- c(lines, capture.output(print(classes_per_season, n = nrow(classes_per_season))))
 lines <- c(lines, "")
 
 lines <- c(lines, "TOP DRUG CLASSES BY MEAN RELATIVE ABUNDANCE")
 lines <- c(lines, "-------------------------------------------")
 lines <- c(lines, capture.output(print(
   top10_classes %>%
     select(DrugClass, prevalence_groups_n, prevalence_groups_pct, mean_pct, median_pct, max_pct),
   n = 10
 )))
 lines <- c(lines, "")
 
 lines <- c(lines, "FULL DRUG-CLASS SUMMARY")
 lines <- c(lines, "-----------------------")
 lines <- c(lines, capture.output(print(
   drugclass_summary %>%
     select(DrugClass, prevalence_groups_n, prevalence_groups_pct, mean_pct, median_pct, max_pct),
   n = nrow(drugclass_summary)
 )))
 lines <- c(lines, "")
 
 lines <- c(lines, "FILES WRITTEN")
 lines <- c(lines, "-------------")
 lines <- c(lines, "tables/ARG_drugclass_summary_mergedGroups.tsv")
 lines <- c(lines, "tables/ARG_nDrugClasses_perMergedGroup.tsv")
 lines <- c(lines, "tables/ARG_nDrugClasses_perMergedGroup_bySeason.tsv")
 lines <- c(lines, "tables/ARG_drugclass_composition_mergedGroups_long.tsv")
 lines <- c(lines, "tables/ARG_drugclass_matrix_relative_mergedGroups.tsv")
 lines <- c(lines, report_file)
 
 writeLines(lines, report_file)
 cat("Wrote:", report_file, "\n") 
 
 classes_per_sample <- class_sample %>%
   group_by(Sample) %>%
   summarise(
     n_drug_classes = sum(Class_abund > 0, na.rm = TRUE),
     .groups = "drop"
   )
 
 mean_classes_per_sample   <- mean(classes_per_sample$n_drug_classes, na.rm = TRUE)
 median_classes_per_sample <- median(classes_per_sample$n_drug_classes, na.rm = TRUE)
 
 mean_classes_per_sample
 median_classes_per_sample
 
 
 ###############################################################################
 # ARG drug classes BETADIVERSITY
 # - PCoA (Bray-Curtis) on replicate-level profiles
 # - PERMANOVA / PERMDISP on replicate-merged profiles
 # - Input ARG matrix: genes x samples (ARG per 16S)
 # - Input annotation: ARG_ID / BestHit / DrugClass
 ###############################################################################
 
 #### 0) Packages --------------------------------------------------------------
 pkgs <- c("dplyr","readr","tidyr","stringr","tibble","ggplot2","vegan")
 to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
 if(length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
 invisible(lapply(pkgs, library, character.only = TRUE))
 
 dir.create("figures", showWarnings = FALSE, recursive = TRUE)
 dir.create("tables",  showWarnings = FALSE, recursive = TRUE)
 
 #### 1) Inputs ----------------------------------------------------------------
 arg_matrix_file <- "ARG_per16S_T1T2_renamed.tsv"   # genes x samples
 arg2class_file  <- "ARG_annotation_RGI_T1T2.tsv"   # columns: ARG_ID, BestHit, DrugClass
 meta_file       <- "metadata3.txt"                 # optional but recommended
 
 # If TRUE: split abundance equally if one gene maps to multiple drug classes
 # If FALSE: abundance is counted in each mapped class (same logic as before)
 SPLIT_MULTI_CLASS <- FALSE
 
 # plot settings
 base_cols_map <- c(
   BF = "#1b9e77", BR = "#a6dba0",
   CF = "#d95f02", CR = "#fdb863",
   PF = "#7570b3", PR = "#b2abd2",
   SF = "#1f9ac2", SR = "#a6dce7",
   VF = "#e7298a", VR = "#f2b2d4"
 )
 
 shape_map <- c(July = 16, September = 17)
 
 #### 2) Helper ----------------------------------------------------------------
 clean_gene <- function(x){
   x <- as.character(x)
   x <- sub("^.*\\|", "", x, perl = TRUE)   # keep text after last |
   x <- sub("^[0-9]+", "", x, perl = TRUE)  # remove leading numbers
   x <- sub("^[-_ :]+", "", x, perl = TRUE) # remove leading separators
   stringr::str_trim(x)
 }
 
 read_matrix_tsv_robust <- function(path){
   x <- read.delim(path, header = TRUE, sep = "\t",
                   check.names = FALSE, fill = TRUE,
                   quote = "", comment.char = "",
                   stringsAsFactors = FALSE)
   idcol <- names(x)[1]
   ids <- x[[idcol]]
   keep <- !(is.na(ids) | trimws(ids) == "")
   x <- x[keep, , drop = FALSE]
   rownames(x) <- make.unique(as.character(x[[idcol]]))
   x[[idcol]] <- NULL
   x[] <- lapply(x, function(v) suppressWarnings(as.numeric(v)))
   m <- as.matrix(x)
   m[is.na(m)] <- 0
   storage.mode(m) <- "numeric"
   m
 }
 
 #### 3) Read ARG matrix -> long -----------------------------------------------
 arg_raw <- readr::read_tsv(arg_matrix_file, show_col_types = FALSE)
 arg_id_col <- names(arg_raw)[1]
 
 ARG_long <- arg_raw %>%
   pivot_longer(-all_of(arg_id_col), names_to = "Sample", values_to = "Abundance") %>%
   mutate(
     Abundance      = suppressWarnings(as.numeric(Abundance)),
     Abundance      = ifelse(is.na(Abundance), 0, Abundance),
     ARG_ID         = as.character(.data[[arg_id_col]]),
     ARG_gene_clean = clean_gene(ARG_ID)
   ) %>%
   select(Sample, ARG_ID, ARG_gene_clean, Abundance)
 
 #### 4) Read mapping ARG -> DrugClass -----------------------------------------
 ann <- readr::read_tsv(arg2class_file, show_col_types = FALSE)
 names(ann) <- trimws(names(ann))
 
 stopifnot(all(c("ARG_ID","BestHit","DrugClass") %in% names(ann)))
 
 ab_map <- ann %>%
   transmute(
     ARG_ID    = as.character(ARG_ID),
     BestHit   = as.character(BestHit),
     DrugClass = stringr::str_trim(as.character(DrugClass))
   ) %>%
   separate_rows(DrugClass, sep = "\\s*[,;/|]+\\s*") %>%
   filter(!is.na(DrugClass), DrugClass != "") %>%
   mutate(
     Gene_from_ARGID   = clean_gene(ARG_ID),
     Gene_from_BestHit = clean_gene(BestHit)
   ) %>%
   pivot_longer(
     cols = c(Gene_from_ARGID, Gene_from_BestHit),
     names_to = "src", values_to = "Gene"
   ) %>%
   mutate(Gene = stringr::str_trim(Gene)) %>%
   filter(!is.na(Gene), Gene != "") %>%
   distinct(Gene, DrugClass)
 
 gene_nclass <- ab_map %>%
   distinct(Gene, DrugClass) %>%
   count(Gene, name = "n_classes")
 
 #### 5) Metadata ---------------------------------------------------------------
 # Use metadata file if available, otherwise parse from sample names
 if (file.exists(meta_file)) {
   meta <- readr::read_tsv(meta_file, show_col_types = FALSE)
   if (!("Sample" %in% names(meta))) {
     cand <- names(meta)[tolower(names(meta)) %in% c("sample","id")]
     if (length(cand) == 0) stop("metadata3.txt must contain a Sample column.")
     meta <- meta %>% rename(Sample = all_of(cand[1]))
   }
   meta <- meta %>% mutate(Sample = as.character(Sample))
 } else {
   meta <- tibble(Sample = unique(ARG_long$Sample))
 }
 
 # fallback parsing if needed
 if (!("SamplingSite" %in% names(meta))) meta$SamplingSite <- substr(meta$Sample, 1, 2)  # BF, BR...
 if (!("Location" %in% names(meta)))     meta$Location     <- substr(meta$Sample, 1, 1)  # B, C...
 if (!("Season" %in% names(meta))) {
   scode <- substr(meta$Sample, 3, 3)
   meta$Season <- dplyr::recode(scode, "J" = "July", "S" = "September", .default = scode)
 }
 if (!("Type" %in% names(meta))) {
   meta$Type <- ifelse(substr(meta$SamplingSite, 2, 2) == "R", "RG", "non-RG")
 }
 if (!("Valley" %in% names(meta))) {
   candv <- names(meta)[tolower(names(meta)) %in% c("valley","vallata","valle")]
   if (length(candv) > 0) meta <- meta %>% rename(Valley = all_of(candv[1]))
 }
 if (!("Valley" %in% names(meta))) {
   warning("Valley column not found in metadata; valley-based PERMANOVA will not run.")
   meta$Valley <- NA_character_
 }
 
 meta <- meta %>%
   mutate(
     Sample       = as.character(Sample),
     SamplingSite = as.character(SamplingSite),
     Location     = as.character(Location),
     Season       = as.character(Season),
     Type         = as.character(Type),
     Valley       = as.character(Valley),
     Group        = paste(SamplingSite, Season, sep = "_")
   )
 
 #### 6) Build class-level abundance table -------------------------------------
 ARG_class_long <- ARG_long %>%
   left_join(ab_map, by = c("ARG_gene_clean" = "Gene"), relationship = "many-to-many") %>%
   mutate(DrugClass = ifelse(is.na(DrugClass), "Unclassified", DrugClass)) %>%
   left_join(meta, by = "Sample")
 
 if (SPLIT_MULTI_CLASS) {
   ARG_class_long <- ARG_class_long %>%
     left_join(gene_nclass, by = c("ARG_gene_clean" = "Gene")) %>%
     mutate(
       n_classes = ifelse(is.na(n_classes) | n_classes < 1, 1, n_classes),
       Abundance_class = Abundance / n_classes
     )
 } else {
   ARG_class_long <- ARG_class_long %>%
     mutate(Abundance_class = Abundance)
 }
 
 # Class abundance per sample
 class_sample <- ARG_class_long %>%
   group_by(Sample, SamplingSite, Location, Season, Type, Valley, Group, DrugClass) %>%
   summarise(Class_abund = sum(Abundance_class, na.rm = TRUE), .groups = "drop")
 
 ###############################################################################
 # PART A. PCoA on replicate-level profiles ------------------------------------
 ###############################################################################
 
 #### 7A) Sample-level matrix: samples x drug classes ---------------------------
 class_sample_wide <- class_sample %>%
   select(Sample, DrugClass, Class_abund) %>%
   pivot_wider(names_from = DrugClass, values_from = Class_abund, values_fill = 0)
 
 X_sample <- class_sample_wide %>% select(-Sample) %>% as.data.frame()
 rownames(X_sample) <- class_sample_wide$Sample
 X_sample <- as.matrix(X_sample)
 storage.mode(X_sample) <- "numeric"
 
 meta_sample <- meta %>%
   filter(Sample %in% rownames(X_sample)) %>%
   arrange(match(Sample, rownames(X_sample)))
 
 stopifnot(all(meta_sample$Sample == rownames(X_sample)))
 
 # remove all-zero samples (e.g. one sample with no ARGs)
 keep_sample <- rowSums(X_sample) > 0
 X_sample <- X_sample[keep_sample, , drop = FALSE]
 meta_sample <- meta_sample[keep_sample, , drop = FALSE]
 
 # relative composition within sample
 X_sample_rel <- vegan::decostand(X_sample, method = "total")
 d_sample <- vegan::vegdist(X_sample_rel, method = "bray")
 
 # PCoA
 pcoa_sample <- cmdscale(d_sample, k = 2, eig = TRUE)
 coords_sample <- as.data.frame(pcoa_sample$points)
 colnames(coords_sample) <- c("PCoA1","PCoA2")
 coords_sample$Sample <- rownames(coords_sample)
 
 eig_pos <- pcoa_sample$eig[pcoa_sample$eig > 0]
 var_expl <- round(100 * eig_pos / sum(eig_pos), 1)
 var1 <- ifelse(length(var_expl) >= 1, var_expl[1], NA)
 var2 <- ifelse(length(var_expl) >= 2, var_expl[2], NA)
 
 plot_df <- coords_sample %>%
   left_join(meta_sample, by = "Sample")
 
 plot_df$SamplingSite <- factor(plot_df$SamplingSite, levels = names(base_cols_map))
 plot_df$Season <- factor(plot_df$Season, levels = c("July","September"))
 
 p_class_pcoa <- ggplot(plot_df, aes(PCoA1, PCoA2, color = SamplingSite, shape = Season)) +
   geom_point(size = 3, alpha = 0.9) +
   scale_color_manual(values = base_cols_map, drop = FALSE, name = "Site") +
   scale_shape_manual(values = shape_map, drop = FALSE) +
   theme_bw() +
   theme(axis.text = element_text(color = "black")) +
   labs(
     x = paste0("PCoA1 (", var1, "%)"),
     y = paste0("PCoA2 (", var2, "%)"),
     shape = "Season",
     title = "PCoA of ARG drug classes (sample level, Bray-Curtis)"
   )
 
 ggsave("figures/PCoA_ARGdrugClasses_Bray_sampleLevel.png",
        p_class_pcoa, width = 8, height = 6, dpi = 300)
 
 write_tsv(plot_df, "tables/PCoA_ARGdrugClasses_sampleLevel_scores.tsv")
 
 ###############################################################################
 # PART B. PERMANOVA on replicate-merged profiles ------------------------------
 ###############################################################################
 
 #### 7B) Merge replicates by median within site + season + class --------------
 class_group <- class_sample %>%
   group_by(Group, SamplingSite, Location, Season, Type, Valley, DrugClass) %>%
   summarise(Class_abund_med = median(Class_abund, na.rm = TRUE), .groups = "drop") %>%
   mutate(Class_abund_med = ifelse(is.na(Class_abund_med), 0, Class_abund_med))
 
 # relative within group
 class_group <- class_group %>%
   group_by(Group) %>%
   mutate(
     Total_group = sum(Class_abund_med, na.rm = TRUE),
     Class_rel   = ifelse(Total_group > 0, Class_abund_med / Total_group, 0)
   ) %>%
   ungroup()
 
 # matrix groups x classes
 class_group_wide <- class_group %>%
   select(Group, DrugClass, Class_rel) %>%
   pivot_wider(names_from = DrugClass, values_from = Class_rel, values_fill = 0)
 
 X_group <- class_group_wide %>% select(-Group) %>% as.data.frame()
 rownames(X_group) <- class_group_wide$Group
 X_group <- as.matrix(X_group)
 storage.mode(X_group) <- "numeric"
 
 meta_group <- class_group %>%
   distinct(Group, SamplingSite, Location, Season, Type, Valley) %>%
   arrange(match(Group, rownames(X_group)))
 
 stopifnot(all(meta_group$Group == rownames(X_group)))
 
 # remove all-zero groups, if any
 keep_group <- rowSums(X_group) > 0
 X_group <- X_group[keep_group, , drop = FALSE]
 meta_group <- meta_group[keep_group, , drop = FALSE]
 
 d_group <- vegdist(X_group, method = "bray")
 
 #### 8B) PERMANOVA -------------------------------------------------------------
 # only keep terms with >= 2 levels and not all NA
 nlev <- c(
   Valley = nlevels(droplevels(factor(meta_group$Valley))),
   Type   = nlevels(droplevels(factor(meta_group$Type))),
   Season = nlevels(droplevels(factor(meta_group$Season)))
 )
 
 keep_terms <- names(nlev)[nlev >= 2]
 
 if (length(keep_terms) == 0) {
   stop("No factor with >=2 levels available for PERMANOVA.")
 }
 
 meta_group <- meta_group %>%
   mutate(
     Valley = factor(Valley),
     Type   = factor(Type),
     Season = factor(Season, levels = c("July","September"))
   )
 
 form <- as.formula(paste("d_group ~", paste(keep_terms, collapse = " + ")))
 
 set.seed(1)
 perm_marg  <- adonis2(form, data = meta_group, permutations = 999, by = "margin")
 perm_terms <- adonis2(form, data = meta_group, permutations = 999, by = "terms")
 
 write_tsv(as.data.frame(perm_marg) %>% rownames_to_column("Term"),
           "tables/PERMANOVA_ARGdrugClasses_GroupMerged_MARGIN.tsv")
 write_tsv(as.data.frame(perm_terms) %>% rownames_to_column("Term"),
           "tables/PERMANOVA_ARGdrugClasses_GroupMerged_TERMS.tsv")
 
 #### 9B) PERMDISP --------------------------------------------------------------
 disp_list <- list()
 
 if (nlev["Valley"] >= 2) {
   bd <- betadisper(d_group, meta_group$Valley)
   a  <- anova(bd)
   disp_list$Valley <- tibble(Factor = "Valley", F = a$`F value`[1], p = a$`Pr(>F)`[1])
 }
 
 if (nlev["Type"] >= 2) {
   bd <- betadisper(d_group, meta_group$Type)
   a  <- anova(bd)
   disp_list$Type <- tibble(Factor = "Type", F = a$`F value`[1], p = a$`Pr(>F)`[1])
 }
 
 if (nlev["Season"] >= 2) {
   bd <- betadisper(d_group, meta_group$Season)
   a  <- anova(bd)
   disp_list$Season <- tibble(Factor = "Season", F = a$`F value`[1], p = a$`Pr(>F)`[1])
 }
 
 disp_tbl <- bind_rows(disp_list) %>%
   mutate(p_adj = p.adjust(p, method = "BH"))
 
 write_tsv(disp_tbl, "tables/PERMDISP_ARGdrugClasses_GroupMerged.tsv")
 
 #### 10) Save matrix for later use --------------------------------------------
 write_tsv(class_group_wide, "tables/ARGdrugClasses_matrix_relative_GroupMerged.tsv")
 
 cat("Saved:\n")
 cat("- figures/PCoA_ARGdrugClasses_Bray_sampleLevel.png\n")
 cat("- tables/PCoA_ARGdrugClasses_sampleLevel_scores.tsv\n")
 cat("- tables/PERMANOVA_ARGdrugClasses_GroupMerged_MARGIN.tsv\n")
 cat("- tables/PERMANOVA_ARGdrugClasses_GroupMerged_TERMS.tsv\n")
 cat("- tables/PERMDISP_ARGdrugClasses_GroupMerged.tsv\n")
 cat("- tables/ARGdrugClasses_matrix_relative_GroupMerged.tsv\n")
 
 
