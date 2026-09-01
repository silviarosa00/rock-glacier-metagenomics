#### Beta diversity of metal-associated MRG genes (Bray–Curtis) + metadata integration
# This chunk performs:
# 1) PCoA (Bray–Curtis) on metal-associated gene abundances (tables/total_metal_mrg_per16s.tsv)
#    - points coloured by Base (BF, BR, ...)
#    - point shapes by Season (July vs September)
# 2) PERMANOVA on Bray–Curtis distances
#    - model: Season + Type with permutations constrained within Site (strata = Site)
#    - model: Site + Season + Type (no strata) to quantify spatial effect
#    - outputs saved to tables/
# 3) Homogeneity of dispersion checks (betadisper) for Season and Type
#    - outputs saved to tables/
# 4) Environmental fitting (envfit) of water chemistry metals (µg/L) onto the ordination
#    - reads chemistry table (user-specified path) with a 'Sample' column matching IDs (e.g., BFJ1, BFS1, ...)
#    - tests each metal by permutation, applies FDR correction (BH)
#    - automatically draws arrows for metals with FDR-adjusted p < 0.05 (up to top_k by R²);
#      if none pass FDR, draws top_k by R² (to show strongest correlations)

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(vegan)

dir.create("tables", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# -----------------------------
# INPUTS
# -----------------------------
mat_file  <- "tables/total_metal_mrg_per16s.tsv"

# <-- cambia questo con il TUO file di chimica (campioni x metalli in µg/L)
chem_file <- "metadata3.txt"   # deve avere colonna Sample

# -----------------------------
# 1) Read matrix gene x sample (metal-only genes) -> sample x gene
# -----------------------------
mat <- read_tsv(mat_file, show_col_types = FALSE)
stopifnot("BacMet_ID" %in% names(mat))

X <- mat %>%
  pivot_longer(cols = -BacMet_ID, names_to = "Sample", values_to = "Abundance") %>%
  mutate(Abundance = as.numeric(Abundance),
         Abundance = ifelse(is.na(Abundance), 0, Abundance)) %>%
  pivot_wider(names_from = BacMet_ID, values_from = Abundance, values_fill = 0)

samples <- X$Sample
X_mat <- X %>% select(-Sample) %>% as.data.frame()
rownames(X_mat) <- samples

# -----------------------------
# 2) Metadata from sample names
# -----------------------------
meta <- tibble(Sample = samples) %>%
  mutate(
    Base   = str_sub(Sample, 1, 2),
    Site   = str_sub(Sample, 1, 1),
    Season_code = str_sub(Sample, 3, 3),
    Season = recode(Season_code, "J"="July", "S"="September", .default = Season_code),
    Replicate = str_extract(Sample, "\\d+$"),
    Type = ifelse(str_sub(Base, 2, 2) == "R", "RG", "non-RG")
  )

# remove all-zero samples (safety)
keep <- rowSums(X_mat) > 0
X_mat <- X_mat[keep, , drop = FALSE]
meta  <- meta[keep, , drop = FALSE]

# -----------------------------
# 3) Bray–Curtis + PCoA
# -----------------------------
X_rel  <- vegan::decostand(X_mat, method = "total")   # relative abundance per sample
d_bray <- vegan::vegdist(X_rel, method = "bray") # Bray-Curtis on relative abundance

pcoa <- cmdscale(d_bray, k = 2, eig = TRUE)
coords <- as.data.frame(pcoa$points)
colnames(coords) <- c("PCoA1","PCoA2")
coords$Sample <- rownames(coords)

var_expl <- round(100 * pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]), 1)
var1 <- var_expl[1]; var2 <- var_expl[2]

plot_df <- coords %>% left_join(meta, by = "Sample")

# -----------------------------
# 4) PERMANOVA + betadisper (save to files)
# -----------------------------
# Test "within Site" (recommended): strata = Site
perm_season_type_strata <- vegan::adonis2(d_bray ~ Season + Type, data = meta, permutations = 999, strata = meta$Site)
write.table(perm_season_type_strata,
            "tables/PERMANOVA_metalGenes_Season_Type_strataSite.tsv",
            sep = "\t", quote = FALSE, col.names = NA)

# Test including Site as predictor (no strata)
perm_site_season_type <- vegan::adonis2(d_bray ~ Site + Season + Type, data = meta, permutations = 999)
write.table(perm_site_season_type,
            "tables/PERMANOVA_metalGenes_Site_Season_Type.tsv",
            sep = "\t", quote = FALSE, col.names = NA)

# Dispersion checks (important!)
bd_season <- vegan::betadisper(d_bray, meta$Season)
bd_type   <- vegan::betadisper(d_bray, meta$Type)

write.table(anova(bd_season),
            "tables/betadisper_metalGenes_Season.tsv",
            sep = "\t", quote = FALSE, col.names = NA)
write.table(anova(bd_type),
            "tables/betadisper_metalGenes_Type.tsv",
            sep = "\t", quote = FALSE, col.names = NA)

# -----------------------------
# 5) envfit: add arrows for metals (choose automatically)
# -----------------------------
# Read chemistry (must have Sample column matching)
chem <- read_tsv(chem_file, show_col_types = FALSE)

if (!("Sample" %in% names(chem))) stop("chem_file must contain a 'Sample' column matching your sample IDs (e.g., BFJ1).")

# align chemistry rows to meta order
chem2 <- chem %>%
  semi_join(meta, by = "Sample") %>%
  arrange(match(Sample, meta$Sample))

# keep only numeric columns (metals)
chem_mat <- chem2 %>%
  select(-Sample) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()

# remove columns all NA or zero variance
ok_cols <- sapply(chem_mat, function(x) {
  x2 <- x[is.finite(x)]
  length(x2) > 2 && sd(x2, na.rm = TRUE) > 0
})
chem_mat <- chem_mat[, ok_cols, drop = FALSE]

# envfit on the ordination (use pcoa points)
fit <- vegan::envfit(pcoa$points, chem_mat, permutations = 999)

# -----------------------------
# extract envfit vectors + stats (ALL variables in chem_mat)
# -----------------------------
vec <- as.data.frame(scores(fit, display = "vectors"))
vec$Metal <- rownames(vec)   # 'Metal' here means "variable name" (metals + pH/EC/etc.)

# stats from envfit
r2 <- fit$vectors$r
p  <- fit$vectors$pvals

env_tbl <- tibble::tibble(
  Metal   = names(r2),
  r2      = as.numeric(r2),
  p_value = as.numeric(p)
) %>%
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  dplyr::arrange(p_adj, dplyr::desc(r2))

readr::write_tsv(env_tbl, "tables/envfit_metals_on_PCoA.tsv")

# -----------------------------
# choose arrows to draw
# -----------------------------
top_k <- 10

sig_metals <- env_tbl %>%
  dplyr::filter(is.finite(p_adj), p_adj < 0.05) %>%
  dplyr::slice_max(order_by = r2, n = top_k, with_ties = FALSE)

# fallback: if none significant, draw top_k by r2
if (nrow(sig_metals) == 0) {
  sig_metals <- env_tbl %>%
    dplyr::slice_max(order_by = r2, n = top_k, with_ties = FALSE)
}

# -----------------------------
# build arrow dataframe (robust to Dim1/Dim2 naming)
# -----------------------------
# Take first two columns of vec as the arrow coordinates (Dim1/Dim2, or whatever they are)
vec2 <- vec %>%
  dplyr::select(Metal, 1, 2)

names(vec2)[2:3] <- c("x", "y")

arrow_df <- vec2 %>%
  dplyr::inner_join(sig_metals, by = "Metal")


# Scale arrows to the plotting space
# (simple: scale by max range of axes)
xrange <- diff(range(plot_df$PCoA1, na.rm = TRUE))
yrange <- diff(range(plot_df$PCoA2, na.rm = TRUE))
mult <- 0.6 * min(xrange, yrange)

arrow_df <- arrow_df %>%
  mutate(
    strength = sqrt(r2),          # 0..1
    xend = x * mult * strength,
    yend = y * mult * strength
  )


# -----------------------------
# 6) Plot: color by Base, shape by Season, add arrows
# -----------------------------
#################################################
# --- Base colour map (as you want) ---
base_cols_map <- c(
  BF = "#1b9e77", BR = "#a6dba0",
  CF = "#d95f02", CR = "#fdb863",
  PF = "#7570b3", PR = "#b2abd2",
  SF = "#1f9ac2", SR = "#a6dce7",
  VF = "#e7298a", VR = "#f2b2d4"
)

# offset perpendicolare alla freccia (in coordinate del plot)
xrange <- diff(range(plot_df$PCoA1, na.rm = TRUE))
yrange <- diff(range(plot_df$PCoA2, na.rm = TRUE))

# --- label placement: small offset that scales with arrow length ---
xrange <- diff(range(plot_df$PCoA1, na.rm = TRUE))
yrange <- diff(range(plot_df$PCoA2, na.rm = TRUE))

arrow_df <- arrow_df %>%
  mutate(
    # arrow length in plot units
    alen = sqrt(xend^2 + yend^2),
    
    # offset size: proportional to arrow length (min/max to avoid extremes)
    off = pmin(pmax(0.015 * min(xrange, yrange), 0.12 * alen), 0.04 * min(xrange, yrange)),
    
    # unit perpendicular direction
    perp_x = -yend,
    perp_y =  xend,
    perp_len = sqrt(perp_x^2 + perp_y^2),
    perp_xu = ifelse(perp_len > 0, perp_x / perp_len, 0),
    perp_yu = ifelse(perp_len > 0, perp_y / perp_len, 0),
    
    # label near arrow tip + scaled perpendicular offset
    label_x = xend + off * perp_xu,
    label_y = yend + off * perp_yu,
    
    # auto text alignment: push label outward from the origin
    hjust = ifelse(label_x >= 0, -0.1, 1.1),
    vjust = ifelse(label_y >= 0, -0.2, 1.2)
  )

# forza Base come factor con i livelli della palette (così i colori matchano al 100%)
plot_df$Base <- factor(plot_df$Base, levels = names(base_cols_map))

# (opzionale) controlla se ci sono Base non riconosciuti
setdiff(unique(as.character(plot_df$Base)), names(base_cols_map))


p <- ggplot(plot_df, aes(PCoA1, PCoA2, color = Base, shape = Season)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = base_cols_map, drop = FALSE) +
  geom_segment(
    data = arrow_df,
    aes(x = 0, y = 0, xend = xend, yend = yend),
    inherit.aes = FALSE,
    linewidth = 0.6,
    arrow = arrow(length = unit(0.22, "cm"))
  ) +
  geom_text(
    data = arrow_df,
    aes(x = label_x, y = label_y, label = Metal, hjust = hjust, vjust = vjust),
    inherit.aes = FALSE,
    size = 3
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(10, 25, 10, 10)) +
  theme_bw() +
  labs(
    x = paste0("PCoA1 (", var1, "%)"),
    y = paste0("PCoA2 (", var2, "%)"),
    title = "PCoA of metal-associated MRG genes (Bray–Curtis)"
  )
p <- p + labs(color = "Site")


ggsave("figures/PCoA_metalGenes_Bray_colorBase_shapeSeason_envfit.png",
       p, width = 8, height = 6, dpi = 300)


############################################################################
###########################################################################
########################### NMDS 

set.seed(1)

# NMDS (2D) on Bray-Curtis
nmds <- metaMDS(X_rel, distance = "bray", k = 2, trymax = 200, autotransform = FALSE)

# stress (da riportare)
nmds$stress

# coordinates
nmds_df <- as.data.frame(scores(nmds, display = "sites"))
nmds_df$Sample <- rownames(nmds_df)
colnames(nmds_df)[1:2] <- c("NMDS1", "NMDS2")

plot_df2 <- nmds_df %>% left_join(meta, by = "Sample")

# envfit on NMDS
fit2 <- envfit(nmds, chem_mat, permutations = 999)

vec <- as.data.frame(scores(fit2, display = "vectors"))
vec$Metal <- rownames(vec)

r2 <- fit2$vectors$r
p  <- fit2$vectors$pvals

env_tbl <- tibble(
  Metal = names(r2),
  r2 = as.numeric(r2),
  p_value = as.numeric(p)
) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj, desc(r2))

top_k <- 6
sig_metals <- env_tbl %>%
  filter(is.finite(p_adj), p_adj < 0.05) %>%
  slice_max(order_by = r2, n = top_k, with_ties = FALSE)

if (nrow(sig_metals) == 0) {
  sig_metals <- env_tbl %>% slice_max(order_by = r2, n = top_k, with_ties = FALSE)
}

# robust: first two columns are NMDS1/NMDS2 directions
vec2 <- vec %>% select(Metal, 1, 2)
names(vec2)[2:3] <- c("x", "y")

arrow_df <- vec2 %>% inner_join(sig_metals, by = "Metal")

# scale arrows; optionally proportional to sqrt(r2)
xrange <- diff(range(plot_df2$NMDS1, na.rm = TRUE))
yrange <- diff(range(plot_df2$NMDS2, na.rm = TRUE))
mult <- 0.6 * min(xrange, yrange)

arrow_df <- arrow_df %>%
  mutate(
    strength = sqrt(r2),
    xend = x * mult * strength,
    yend = y * mult * strength
  )

# label offset perpendicular (like your PCoA version)
off <- 0.03 * min(xrange, yrange)
arrow_df <- arrow_df %>%
  mutate(
    perp_x = -y, perp_y = x,
    perp_len = sqrt(perp_x^2 + perp_y^2),
    perp_xu = ifelse(perp_len > 0, perp_x / perp_len, 0),
    perp_yu = ifelse(perp_len > 0, perp_y / perp_len, 0),
    label_x = xend + off * perp_xu,
    label_y = yend + off * perp_yu
  )

p_nmds <- ggplot(plot_df2, aes(NMDS1, NMDS2, color = Base, shape = Season)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = base_cols_map, name = "Site") +  # se vuoi solo titolo
  geom_segment(
    data = arrow_df,
    aes(x = 0, y = 0, xend = xend, yend = yend),
    inherit.aes = FALSE,
    linewidth = 0.6,
    arrow = arrow(length = unit(0.22, "cm"))
  ) +
  geom_text(
    data = arrow_df,
    aes(x = label_x, y = label_y, label = Metal),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_bw() +
  labs(
    title = paste0("NMDS of metal-associated MRG genes (Bray–Curtis), stress = ", round(nmds$stress, 3))
  )

ggsave("figures/NMDS_metalGenes_Bray_colorBase_shapeSeason_envfit.png",
       p_nmds, width = 8, height = 6, dpi = 300)

write_tsv(env_tbl, "tables/envfit_metals_on_NMDS.tsv")


#####################################
#####################################
#####################################

perm_terms <- vegan::adonis2(d_bray ~ Site + Season + Type, data = meta,
                             permutations = 999, by = "terms")
print(perm_terms)
perm_margin <- vegan::adonis2(d_bray ~ Site + Season + Type, data = meta,
                              permutations = 999, by = "margin")
print(perm_margin)


perm_terms_tbl <- as.data.frame(perm_terms) %>%
  rownames_to_column("Term")

write_tsv(perm_terms_tbl, "tables/PERMANOVA_metalGenes_Site_Season_Type_byTerms.tsv")



# ==== SOLO METAL-MRG ====
mrg_file <- "tables/total_metal_mrg_per16s.tsv"

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
  m <- as.matrix(x)
  m[is.na(m)] <- 0
  m
}

m <- read_matrix_tsv_robust(mrg_file)  # genes x samples

# (1) geni distinti presenti in >=1 campione
n_genes_detected <- sum(rowSums(m > 0, na.rm = TRUE) > 0)

# (2) total load per 16S per campione
total_load <- colSums(m, na.rm = TRUE)

# (3) richness per campione
richness <- colSums(m > 0, na.rm = TRUE)

cat("\n=== METAL-MRG SUMMARY ===\n")
cat("File:", mrg_file, "\n")
cat("N samples:", ncol(m), "\n")
cat("Detected metal-MRG genes (>=1 sample):", n_genes_detected, "\n")
cat("Median total metal-MRG load per 16S:", round(median(total_load, na.rm = TRUE)), "\n")
cat("Median metal-MRG richness (genes):", round(median(richness, na.rm = TRUE)), "\n\n")

cat(sprintf(
  "Paper sentence:\nA total of %d distinct metal-associated MRGs were identified across all samples; the median total metal-MRG load per 16S was %.0f, and the median per-sample metal-MRG richness was %.0f genes.\n",
  n_genes_detected,
  median(total_load, na.rm = TRUE),
  median(richness, na.rm = TRUE)
))

###################################################
#### BRAY CURTIS ACROSS SAMPLES ##################


mrg_wide <- read_tsv("tables/total_metal_mrg_per16s.tsv", show_col_types = FALSE)
gene_col <- names(mrg_wide)[1]
sample_cols <- setdiff(names(mrg_wide), gene_col)

X <- mrg_wide %>%
  mutate(across(all_of(sample_cols), as.numeric)) %>%
  select(all_of(sample_cols)) %>%
  t() %>% as.data.frame()

# opzionale ma prudente: rimuovi eventuali profili all-zero
X <- X[rowSums(X, na.rm = TRUE) > 0, , drop = FALSE]

bc <- vegdist(X, method = "bray")
d <- as.vector(bc)

cat(sprintf(
  "Bray-Curtis dissimilarity across individual replicate profiles: median %.3f (range %.3f–%.3f)\n",
  median(d), min(d), max(d)
))


# 1) leggi matrice metal-MRG (genes x samples)
mrg_wide <- read_tsv("tables/total_metal_mrg_per16s.tsv", show_col_types = FALSE)
gene_col <- names(mrg_wide)[1]
sample_cols <- setdiff(names(mrg_wide), gene_col)

# 2) long format + metadata dai nomi campione
mrg_long <- mrg_wide %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "Sample", values_to = "Abundance") %>%
  mutate(
    Abundance = as.numeric(Abundance),
    Abundance = ifelse(is.na(Abundance), 0, Abundance),
    sampling_site = str_sub(Sample, 1, 2),   # BF, BR, CF...
    location      = str_sub(Sample, 1, 1),   # B, C, P, S, V
    SeasonCode    = str_sub(Sample, 3, 3),
    Season        = recode(SeasonCode, "J" = "July", "S" = "September", .default = SeasonCode)
  )

# 3) ricostruisci matrice replicate-merged: sampling_site x season x gene
mat_wide <- mrg_long %>%
  group_by(sampling_site, Season, !!sym(gene_col)) %>%
  summarise(Abund_med = median(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    SeasonCode = ifelse(Season == "July", "J",
                        ifelse(Season == "September", "S", Season)),
    Group = paste0(sampling_site, SeasonCode)   # BFJ, BFS, ...
  ) %>%
  select(Group, !!sym(gene_col), Abund_med) %>%
  pivot_wider(names_from = !!sym(gene_col), values_from = Abund_med, values_fill = 0)

# 4) matrix con rownames
X <- mat_wide %>% select(-Group) %>% as.data.frame()
rownames(X) <- mat_wide$Group

# opzionale ma prudente
X <- X[rowSums(X, na.rm = TRUE) > 0, , drop = FALSE]

# 5) Bray-Curtis
bc <- vegdist(X, method = "bray")
D  <- as.matrix(bc)

# 6) pairwise table
pairs <- expand.grid(i = seq_len(nrow(D)), j = seq_len(nrow(D))) %>%
  filter(i < j) %>%
  mutate(
    G1 = rownames(D)[i],
    G2 = colnames(D)[j],
    BC = D[cbind(i, j)],
    sampling_site1 = str_sub(G1, 1, 2),
    sampling_site2 = str_sub(G2, 1, 2),
    SeasonCode1 = str_sub(G1, 3, 3),
    SeasonCode2 = str_sub(G2, 3, 3),
    Season1 = recode(SeasonCode1, "J" = "July", "S" = "September", .default = NA_character_),
    Season2 = recode(SeasonCode2, "J" = "July", "S" = "September", .default = NA_character_),
    location1 = str_sub(sampling_site1, 1, 1),
    location2 = str_sub(sampling_site2, 1, 1)
  ) %>%
  filter(!is.na(Season1), Season1 == Season2) %>%
  mutate(
    Season = Season1,
    within_location = (location1 == location2) & (sampling_site1 != sampling_site2)
  )

# 7) Wilcoxon within-location vs between-location
wilcox_dist <- function(df) {
  x <- df$BC[df$within_location]
  y <- df$BC[!df$within_location]
  
  if (length(x) < 2 || length(y) < 2) {
    return(tibble(
      n_within = length(x),
      n_between = length(y),
      median_within = ifelse(length(x) > 0, median(x), NA_real_),
      median_between = ifelse(length(y) > 0, median(y), NA_real_),
      p_value = NA_real_
    ))
  }
  
  tibble(
    n_within = length(x),
    n_between = length(y),
    median_within = median(x),
    median_between = median(y),
    p_value = wilcox.test(x, y, exact = FALSE)$p.value
  )
}

dist_overall <- wilcox_dist(pairs)

dist_by_season <- pairs %>%
  group_by(Season) %>%
  group_modify(~ wilcox_dist(.x)) %>%
  ungroup()

print(dist_overall)
print(dist_by_season)

# opzionale: salva anche i risultati
write_tsv(pairs, "tables/metalMRG_pairwise_within_vs_between_location.tsv")
write_tsv(dist_overall, "tables/metalMRG_within_vs_between_location_overall.tsv")
write_tsv(dist_by_season, "tables/metalMRG_within_vs_between_location_bySeason.tsv")

