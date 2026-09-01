###############################################################################
# Multi-omics vs Geochemistry (Group-merged) - ONE SCRIPT
###############################################################################

pkgs <- c("dplyr","tibble","tidyr","readr","stringr","vegan","ape","ggplot2")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(to_install)) install.packages(to_install, repos="https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# --------------------------- PATHS (EDIT) -----------------------------------
out_dir <- "integrated_multiomics_out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "matrices"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "models"),   showWarnings = FALSE, recursive = TRUE)

# 1) Families: usa comm_filt se ce l'hai già in workspace.
#    Se invece vuoi leggere da file, metti il path qui (taxa x sample o sample x taxa)
families_file <-"output_beta_F/Families_filtered_rel_taxa_x_sample.tsv" 

# 2) HUMAnN pathways
pathway_rel_file <- "T12_pathabundance_rel_renamed.tsv"     # <-- cambia


# 3) ARG /16S
arg_file <- "ARG_per16S_T1T2_renamed.tsv"          # <-- cambia

# 4) MRG /16S (BacMet)
mrg_file <- "total_metal_mrg_per16s.tsv"           # <-- cambia

# 5) BacMet annotation (metal only)
bacmet_is_metal_file <- "tables/bacmet_gene_is_metal_ALL.tsv"  # <-- cambia

# 6) Geochemistry metadata
geo_file <- "metadata4.csv"                         # <-- cambia
# -----------------------"-----------------------------------------------------

# --------------------------- HELPERS ----------------------------------------

read_feature_table <- function(path, id_col = 1){
  x <- readr::read_delim(path, delim = "\t", show_col_types = FALSE)
  if(is.numeric(id_col)) id_col <- names(x)[id_col]
  stopifnot(id_col %in% names(x))
  ids <- x[[id_col]]
  mat <- x %>% select(-all_of(id_col)) %>% as.data.frame()
  mat <- as.matrix(sapply(mat, function(v) suppressWarnings(as.numeric(v))))
  rownames(mat) <- make.unique(as.character(ids))
  mat[is.na(mat)] <- 0
  mat
}

# accepts either feature x sample OR sample x feature; returns sample x feature
as_sample_x_feature <- function(mat){
  # heuristic: samples are columns if column names look like BFJ1/BFJ etc
  # If rownames look like samples and colnames like features, keep as is.
  rn <- rownames(mat); cn <- colnames(mat)
  is_sample_name <- function(x) grepl("^[A-Z]{2}[JS]\\d*$", x) | grepl("^[A-Z]{3}\\d*$", x) | grepl("^[A-Z]{2}[A-Z][JS]\\d*$", x)
  rn_samples <- !is.null(rn) && mean(is_sample_name(rn)) > 0.6
  cn_samples <- !is.null(cn) && mean(is_sample_name(cn)) > 0.6
  if(cn_samples && !rn_samples) {
    return(t(mat))  # feature x sample -> sample x feature
  } else {
    return(mat)     # already sample x feature
  }
}

merge_replicates_to_group <- function(mat_sxf){
  stopifnot(!is.null(rownames(mat_sxf)))
  samp <- rownames(mat_sxf)
  grp <- sub("\\d+$", "", samp)  # BFJ1->BFJ, BFJ2->BFJ
  out <- sapply(split(seq_len(nrow(mat_sxf)), grp), function(ix){
    colMedians <- apply(mat_sxf[ix, , drop=FALSE], 2, median, na.rm=TRUE)
    colMedians
  })
  out <- t(out)  # groups x features
  rownames(out) <- names(split(seq_len(nrow(mat_sxf)), grp))
  storage.mode(out) <- "numeric"
  out[is.na(out)] <- 0
  out
}

# Extract BacMet_ID prefix (BAC0003) from strings like "BAC0003|acr3|..."
extract_bacmet_prefix <- function(x){
  sub("^([A-Z]{3}\\d+).*", "\\1", x)
}

# Clean geochem column names -> T, ph, EC, Al, Cd, ...
clean_geo_names <- function(nm){
  nm2 <- nm
  nm2 <- gsub("\n", " ", nm2)
  nm2 <- gsub("\"", "", nm2)
  
  # standard names
  nm2 <- gsub("^EC.*", "EC", nm2)
  nm2 <- gsub("^pH.*", "ph", nm2)
  nm2 <- gsub("^Temp.*", "T", nm2)
  nm2 <- gsub("^alt.*", "alt", nm2)
  
  # elements: take first token (letters/digits) at beginning
  nm2 <- ifelse(nm2 %in% c("EC","ph","T","alt","Sample"), nm2,
                sub("^\\s*([A-Za-z]{1,2})\\b.*", "\\1", nm2))
  nm2
}

to_num <- function(x){
  if(is.character(x)){
    x <- gsub(",", ".", x)
    x <- gsub("\\s*\\[.*?\\]", "", x)  # remove [..]
    x <- gsub("<LOD", NA, x)
    x <- gsub("<LOD\\s*", NA, x)
  }
  suppressWarnings(as.numeric(x))
}

median_impute <- function(df){
  for(nm in names(df)){
    if(is.numeric(df[[nm]])){
      med <- median(df[[nm]], na.rm=TRUE)
      df[[nm]][is.na(df[[nm]])] <- med
    }
  }
  df
}

# Run dbRDA with VIF pruning
run_capscale_vif <- function(Y_sxf, X_df, dist="bray", vif_thr=10){
  # remove near-zero variance response cols
  keepY <- apply(Y_sxf, 2, function(x) var(x, na.rm=TRUE) > 0)
  Y <- Y_sxf[, keepY, drop=FALSE]
  
  X <- X_df
  mod <- vegan::capscale(Y ~ ., data=X, distance=dist)
  v <- vegan::vif.cca(mod)
  
  while(length(v) > 1 && max(v, na.rm=TRUE) > vif_thr){
    drop <- names(which.max(v))
    X <- X[, setdiff(names(X), drop), drop=FALSE]
    mod <- vegan::capscale(Y ~ ., data=X, distance=dist)
    v <- vegan::vif.cca(mod)
  }
  
  list(model=mod, X=X, vif=v)
}

extract_capscale_tables <- function(mod){
  glob <- as.data.frame(anova.cca(mod, permutations=999)) %>% tibble::rownames_to_column("term")
  byterm <- as.data.frame(anova.cca(mod, by="term", permutations=999)) %>% tibble::rownames_to_column("term")
  byterm <- byterm %>% mutate(p_adj = p.adjust(`Pr(>F)`, "BH"))
  rsq <- vegan::RsquareAdj(mod)
  list(global=glob, term=byterm, rsq=rsq)
}

# PCoA scores (first 2 axes) from Bray–Curtis on sample x feature
pcoa2 <- function(Y_sxf){
  d <- vegan::vegdist(Y_sxf, method="bray")
  p <- ape::pcoa(d)
  sc <- as.data.frame(p$vectors[,1:2])
  colnames(sc) <- c("Axis1","Axis2")
  sc$Sample <- rownames(sc)
  list(scores=sc, var = p$values$Relative_eig[1:2]*100)
}

# --------------------------- 1) READ / BUILD MATRICES ------------------------

# Families
if(exists("comm_filt")){
  fam_sxf <- comm_filt  # samples x families (rel)
} else if(!is.na(families_file)){
  fam_raw <- read_feature_table(families_file, 1)
  fam_sxf <- as_sample_x_feature(fam_raw)
} else {
  stop("Non ho comm_filt e families_file è NA. Fornisci comm_filt oppure un file.")
}
stopifnot(!is.null(rownames(fam_sxf)))

# Pathways (rel + optional abs)
path_rel_raw <- read_feature_table(pathway_rel_file, 1)
path_rel_sxf <- as_sample_x_feature(path_rel_raw)

path_abs_sxf <- NULL
if(!is.na(pathway_abs_file) && file.exists(pathway_abs_file)){
  path_abs_raw <- read_feature_table(pathway_abs_file, 1)
  path_abs_sxf <- as_sample_x_feature(path_abs_raw)
}

# ARG /16S
arg_raw <- read_feature_table(arg_file, 1)
arg_sxf <- as_sample_x_feature(arg_raw)

# MRG /16S (BacMet)
mrg_raw <- read_feature_table(mrg_file, 1)
mrg_sxf <- as_sample_x_feature(mrg_raw)

# Filter MRG to metal-only using bacmet_is_metal_ALL
bacmet_map <- readr::read_delim(bacmet_is_metal_file, delim="\t", show_col_types = FALSE)
stopifnot(all(c("BacMet_ID","is_metal") %in% names(bacmet_map)))
metal_ids <- bacmet_map %>% filter(is_metal == TRUE) %>% pull(BacMet_ID) %>% unique()

mrg_prefix <- extract_bacmet_prefix(colnames(t(mrg_sxf))) # dummy to keep function loaded
# filter on feature names (columns)
mrg_feat <- colnames(mrg_sxf)
mrg_feat_prefix <- extract_bacmet_prefix(mrg_feat)
keep_mrg <- mrg_feat_prefix %in% metal_ids
mrg_sxf_metal <- mrg_sxf[, keep_mrg, drop=FALSE]

# --------------------------- 2) MERGE REPLICATES -> GROUP -------------------

fam_grp  <- merge_replicates_to_group(fam_sxf)
path_grp <- merge_replicates_to_group(path_rel_sxf)
arg_grp  <- merge_replicates_to_group(arg_sxf)
mrg_grp  <- merge_replicates_to_group(mrg_sxf_metal)

# Save matrices
write_tsv(as.data.frame(fam_grp)  %>% rownames_to_column("Group"), file.path(out_dir, "matrices/Families_groupMerged.tsv"))
write_tsv(as.data.frame(path_grp) %>% rownames_to_column("Group"), file.path(out_dir, "matrices/PathwaysRel_groupMerged.tsv"))
write_tsv(as.data.frame(arg_grp)  %>% rownames_to_column("Group"), file.path(out_dir, "matrices/ARG_per16S_groupMerged.tsv"))
write_tsv(as.data.frame(mrg_grp)  %>% rownames_to_column("Group"), file.path(out_dir, "matrices/MRG_metalOnly_per16S_groupMerged.tsv"))

if(!is.null(path_abs_sxf)){
  path_abs_grp <- merge_replicates_to_group(path_abs_sxf)
  write_tsv(as.data.frame(path_abs_grp) %>% rownames_to_column("Group"), file.path(out_dir, "matrices/PathwaysAbs_groupMerged.tsv"))
}

# --------------------------- 3) READ + CLEAN GEOCHEM ------------------------

geo <- readr::read_delim(geo_file, delim = ",", show_col_types = FALSE)
# keep Sample and numeric predictors; rename columns
names(geo) <- clean_geo_names(names(geo))

stopifnot("Sample" %in% names(geo))
geo <- geo %>% rename(Group = Sample)

# Keep only env vars present (EC/ph/T + elements + alt if available)
env_keep <- intersect(env_cols, names(geo))
# NOTE: your file uses T/ph/EC naming after clean_geo_names; elements become symbols.
geo_X <- geo %>% select(Group, all_of(env_keep))

# numeric conversion
for(nm in env_keep){
  geo_X[[nm]] <- to_num(geo_X[[nm]])
}

# remove vars with too many NA
na_prop <- sapply(geo_X[env_keep], function(x) mean(is.na(x)))
env_keep2 <- env_keep[na_prop <= 0.3]  # keep vars with <=30% NA
geo_X <- geo_X %>% select(Group, all_of(env_keep2))

# impute + scale
X <- geo_X %>%
  select(-Group) %>%
  as.data.frame()
X <- median_impute(X)
X <- as.data.frame(scale(X))
rownames(X) <- geo_X$Group

# --------------------------- 4) ALIGN SAMPLES -------------------------------

common_groups <- Reduce(intersect, list(rownames(fam_grp), rownames(path_grp), rownames(arg_grp), rownames(mrg_grp), rownames(X)))
if(length(common_groups) < 8) stop("Troppi pochi campioni comuni dopo l'allineamento. Controlla nomi Group e metadata4.csv.")
message("Common groups used: ", length(common_groups))

fam_grp  <- fam_grp[common_groups, , drop=FALSE]
path_grp <- path_grp[common_groups, , drop=FALSE]
arg_grp  <- arg_grp[common_groups, , drop=FALSE]
mrg_grp  <- mrg_grp[common_groups, , drop=FALSE]
X        <- X[common_groups, , drop=FALSE]

# optional transformations (helpful)
fam_grp_h <- vegan::decostand(fam_grp, "hellinger")
path_grp_h <- vegan::decostand(path_grp, "hellinger")
arg_grp_l <- log1p(arg_grp)
mrg_grp_l <- log1p(mrg_grp)

# --------------------------- 5) STRADA A: dbRDA per layer -------------------

layers <- list(
  Families = fam_grp_h,
  PathwaysRel = path_grp_h,
  ARG_per16S = arg_grp_l,
  MRG_metal_per16S = mrg_grp_l
)

A_results <- list()

for(nm in names(layers)){
  Y <- layers[[nm]]
  fit <- run_capscale_vif(Y, X, dist="bray", vif_thr=10)
  tabs <- extract_capscale_tables(fit$model)
  
  write_tsv(tabs$global, file.path(out_dir, paste0("models/A_", nm, "_capscale_global.tsv")))
  write_tsv(tabs$term,   file.path(out_dir, paste0("models/A_", nm, "_capscale_terms_BH.tsv")))
  writeLines(c(
    paste0("Layer: ", nm),
    paste0("Predictors kept (after VIF): ", paste(colnames(fit$X), collapse=", ")),
    paste0("R2 adj: ", round(tabs$rsq$r.squared, 3), " ; adj: ", round(tabs$rsq$adj.r.squared, 3))
  ), file.path(out_dir, paste0("models/A_", nm, "_summary.txt")))
  
  A_results[[nm]] <- list(model=fit$model, predictors=colnames(fit$X), rsq=tabs$rsq)
}

saveRDS(A_results, file.path(out_dir, "models/StradaA_models.rds"))

# --------------------------- 6) STRADA B: integra (PCoA axes) + RDA ----------

get_axes <- function(Y, prefix){
  pc <- pcoa2(Y)
  sc <- pc$scores %>% select(Sample, Axis1, Axis2)
  colnames(sc) <- c("Group", paste0(prefix,"_PCoA1"), paste0(prefix,"_PCoA2"))
  list(scores=sc, var=pc$var)
}

ax_tax <- get_axes(fam_grp_h, "TAX")
ax_fun <- get_axes(path_grp_h, "FUN")
ax_arg <- get_axes(arg_grp_l, "ARG")
ax_mrg <- get_axes(mrg_grp_l, "MRG")

Y_axes <- ax_tax$scores %>%
  left_join(ax_fun$scores, by="Group") %>%
  left_join(ax_arg$scores, by="Group") %>%
  left_join(ax_mrg$scores, by="Group") %>%
  column_to_rownames("Group") %>%
  as.matrix()

# RDA on axes (euclidean)
rda_mod <- vegan::rda(Y_axes ~ ., data = X)
rda_glob <- as.data.frame(anova.cca(rda_mod, permutations=999)) %>% rownames_to_column("term")
rda_term <- as.data.frame(anova.cca(rda_mod, by="term", permutations=999)) %>% rownames_to_column("term") %>%
  mutate(p_adj = p.adjust(`Pr(>F)`, "BH"))
rda_rsq <- vegan::RsquareAdj(rda_mod)

write_tsv(rda_glob, file.path(out_dir, "models/B_integratedAxes_RDA_global.tsv"))
write_tsv(rda_term, file.path(out_dir, "models/B_integratedAxes_RDA_terms_BH.tsv"))
writeLines(c(
  "Integrated (PCoA axes) RDA",
  paste0("R2: ", round(rda_rsq$r.squared, 3), " ; adj: ", round(rda_rsq$adj.r.squared, 3)),
  paste0("Tax PCoA1-2 var (%): ", paste(round(ax_tax$var,1), collapse=", ")),
  paste0("Fun PCoA1-2 var (%): ", paste(round(ax_fun$var,1), collapse=", ")),
  paste0("ARG PCoA1-2 var (%): ", paste(round(ax_arg$var,1), collapse=", ")),
  paste0("MRG PCoA1-2 var (%): ", paste(round(ax_mrg$var,1), collapse=", "))
), file.path(out_dir, "models/B_integratedAxes_summary.txt"))

saveRDS(list(rda=rda_mod, rsq=rda_rsq), file.path(out_dir, "models/StradaB_model.rds"))

message("DONE. Outputs in: ", normalizePath(out_dir))

# =========================== FIGURES (add at end) ===========================
dir.create(file.path(out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

# metadata per plotting (da Group code)
meta_plot <- tibble(Group = common_groups) %>%
  mutate(
    Base   = substr(Group, 1, 2),
    Season = substr(Group, 3, 3),
    Type   = substr(Group, 2, 2),
    Site   = substr(Group, 1, 1)
  ) %>%
  mutate(
    Season = factor(Season, levels = c("J","S")),
    Base   = factor(Base)
  )

# colori Base "soliti" (metti i tuoi se vuoi)
base_cols_map <- c(
  BF = "#1b9e77", BR = "#a6dba0",
  CF = "#d95f02", CR = "#fdb863",
  PF = "#7570b3", PR = "#b2abd2",
  SF = "#1f9ac2", SR = "#a6dce7",
  VF = "#e7298a", VR = "#f2b2d4"
)

# helper: plot dbRDA/capscale con frecce
plot_capscale <- function(mod, meta_df, title, fn_png, show_top_if_none = 5){
  sc_sites <- as.data.frame(vegan::scores(mod, display = "sites", choices = 1:2))
  sc_sites$Group <- rownames(sc_sites)
  
  dfp <- sc_sites %>% left_join(meta_df, by = "Group")
  
  # biplot arrows (constraints)
  bp <- try(as.data.frame(vegan::scores(mod, display = "bp", choices = 1:2)), silent = TRUE)
  arrow_df <- NULL
  
  if(!inherits(bp, "try-error") && nrow(bp) > 0){
    bp$var <- rownames(bp)
    colnames(bp)[1:2] <- c("x","y")
    
    # significance per term (BH) dal tuo output già calcolato
    # rifaccio qui in modo indipendente
    tt <- as.data.frame(vegan::anova.cca(mod, by="term", permutations=999)) %>%
      tibble::rownames_to_column("var") %>%
      filter(!var %in% c("Residual")) %>%
      mutate(p_adj = p.adjust(`Pr(>F)`, "BH"))
    
    arrow_df <- bp %>% left_join(tt, by="var")
    
    # tieni solo significative BH<0.05
    arrow_sig <- arrow_df %>% filter(!is.na(p_adj) & p_adj < 0.05)
    
    # se zero frecce BH, mostra top N per F (o r2 proxy) solo per orientamento
    if(nrow(arrow_sig) == 0 && show_top_if_none > 0){
      arrow_sig <- arrow_df %>%
        filter(!is.na(F)) %>%
        arrange(desc(F)) %>%
        head(show_top_if_none)
    }
    
    # scala frecce per essere visibili
    xr <- diff(range(dfp$CAP1, na.rm=TRUE))
    yr <- diff(range(dfp$CAP2, na.rm=TRUE))
    mult <- 0.7 * min(xr, yr)
    
    arrow_sig <- arrow_sig %>%
      mutate(x0=0, y0=0,
             x1 = x * mult,
             y1 = y * mult)
    
    arrow_df <- arrow_sig
  }
  
  p <- ggplot(dfp, aes(x = CAP1, y = CAP2, colour = Base, shape = Season)) +
    geom_point(size = 4, alpha = 0.95) +
    scale_colour_manual(values = base_cols_map, drop = FALSE) +
    theme_bw(base_size = 13) +
    labs(title = title, x = "dbRDA1", y = "dbRDA2", colour = "Base", shape = "Season")
  
  if(!is.null(arrow_df) && nrow(arrow_df) > 0){
    p <- p +
      geom_segment(
        data = arrow_df,
        aes(x = x0, y = y0, xend = x1, yend = y1),
        inherit.aes = FALSE,
        linewidth = 0.7,
        arrow = grid::arrow(length = grid::unit(0.02, "npc"))
      ) +
      ggrepel::geom_text_repel(
        data = arrow_df,
        aes(x = x1, y = y1, label = var),
        inherit.aes = FALSE,
        size = 3,
        max.overlaps = Inf
      )
  }
  
  ggsave(fn_png, p, width = 8, height = 6, dpi = 300)
  p
}

# helper: scree plot (base) per capscale/rda
save_scree <- function(mod, fn_png, main){
  png(fn_png, width=1200, height=800, res=150)
  try({
    vegan::screeplot(mod, main = main)
  }, silent=TRUE)
  dev.off()
}

# ------------------- Strada A: 4 dbRDA plots + scree ------------------------
for(nm in names(A_results)){
  mod <- A_results[[nm]]$model
  
  p <- plot_capscale(
    mod,
    meta_plot,
    title = paste0("dbRDA (Bray–Curtis) — ", nm, " (Group-merged)"),
    fn_png = file.path(out_dir, "figures", paste0("A_dbRDA_", nm, ".png")),
    show_top_if_none = 5
  )
  
  save_scree(mod,
             fn_png = file.path(out_dir, "figures", paste0("A_scree_", nm, ".png")),
             main = paste0("Scree — ", nm))
}

# ------------------- Strada B: RDA integrata plot + scree -------------------
# Punteggi siti + frecce vincoli
sc_sites <- as.data.frame(vegan::scores(rda_mod, display="sites", choices=1:2))
sc_sites$Group <- rownames(sc_sites)
dfB <- sc_sites %>% left_join(meta_plot, by="Group")

bpB <- as.data.frame(vegan::scores(rda_mod, display="bp", choices=1:2))
bpB$var <- rownames(bpB)
colnames(bpB)[1:2] <- c("x","y")

ttB <- as.data.frame(vegan::anova.cca(rda_mod, by="term", permutations=999)) %>%
  tibble::rownames_to_column("var") %>%
  filter(!var %in% c("Residual")) %>%
  mutate(p_adj = p.adjust(`Pr(>F)`, "BH"))

arB <- bpB %>% left_join(ttB, by="var")
arB_sig <- arB %>% filter(!is.na(p_adj) & p_adj < 0.05)

# se zero frecce BH, prendi top 5 per F
if(nrow(arB_sig) == 0){
  arB_sig <- arB %>% filter(!is.na(F)) %>% arrange(desc(F)) %>% head(5)
}

xr <- diff(range(dfB$RDA1, na.rm=TRUE))
yr <- diff(range(dfB$RDA2, na.rm=TRUE))
mult <- 0.7 * min(xr, yr)
arB_sig <- arB_sig %>% mutate(x0=0,y0=0, x1=x*mult, y1=y*mult)

pB <- ggplot(dfB, aes(RDA1, RDA2, colour=Base, shape=Season)) +
  geom_point(size=4, alpha=0.95) +
  scale_colour_manual(values=base_cols_map, drop=FALSE) +
  theme_bw(base_size=13) +
  labs(title="Integrated (PCoA axes) RDA — Group-merged", x="RDA1", y="RDA2")

if(nrow(arB_sig) > 0){
  pB <- pB +
    geom_segment(
      data=arB_sig,
      aes(x=x0,y=y0,xend=x1,yend=y1),
      inherit.aes=FALSE,
      linewidth=0.7,
      arrow=grid::arrow(length=grid::unit(0.02,"npc"))
    ) +
    ggrepel::geom_text_repel(
      data=arB_sig,
      aes(x=x1,y=y1,label=var),
      inherit.aes=FALSE,
      size=3,
      max.overlaps=Inf
    )
}

ggsave(file.path(out_dir, "figures", "B_integratedAxes_RDA.png"), pB, width=8, height=6, dpi=300)
save_scree(rda_mod, file.path(out_dir, "figures", "B_integratedAxes_RDA_scree.png"), "Scree — Integrated RDA")

message("FIGURES saved in: ", normalizePath(file.path(out_dir, "figures")))
# =========================================================================== 

library(dplyr)
library(tibble)
library(ggplot2)
library(vegan)
library(ggrepel)

out_dir <- "integrated_multiomics_out"
dir.create(file.path(out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

base_cols_map <- c(
  BF = "#1b9e77", BR = "#a6dba0",
  CF = "#d95f02", CR = "#fdb863",
  PF = "#7570b3", PR = "#b2abd2",
  SF = "#1f9ac2", SR = "#a6dce7",
  VF = "#e7298a", VR = "#f2b2d4"
)

# -------- helper: points + envfit arrows (works even if model is overfit) ----
plot_sites_envfit <- function(ord_scores, env_df, meta_df, title, fn_png, topN=6){
  dfp <- ord_scores %>% left_join(meta_df, by="Group")
  
  # envfit on ordination
  set.seed(123)
  ef <- envfit(as.matrix(dfp[,c("Axis1","Axis2")]), env_df, permutations=999)
  
  ef_tbl <- tibble(
    var = names(ef$vectors$pvals),
    r2  = as.numeric(ef$vectors$r),
    p   = as.numeric(ef$vectors$pvals),
    p_adj = p.adjust(as.numeric(ef$vectors$pvals), "BH")
  ) %>% arrange(p_adj)
  
  # arrows: BH<0.05 else topN by r2
  keep <- ef_tbl %>% filter(p_adj < 0.05) %>% pull(var)
  if(length(keep) == 0) keep <- ef_tbl %>% arrange(desc(r2)) %>% head(topN) %>% pull(var)
  
  ar <- as.data.frame(scores(ef, display="vectors"))
  ar$var <- rownames(ar)
  ar <- ar %>% filter(var %in% keep) %>% left_join(ef_tbl, by="var")
  
  # scale arrows
  xr <- diff(range(dfp$Axis1, na.rm=TRUE))
  yr <- diff(range(dfp$Axis2, na.rm=TRUE))
  mult <- 0.6 * min(xr, yr)
  ar <- ar %>% mutate(x0=0,y0=0, x1=Axis1*sqrt(r2)*mult, y1=Axis2*sqrt(r2)*mult)
  
  p <- ggplot(dfp, aes(Axis1, Axis2, colour=Base, shape=Season)) +
    geom_point(size=4, alpha=0.95) +
    scale_colour_manual(values=base_cols_map, drop=FALSE) +
    theme_bw(base_size=13) +
    labs(title=title, x="Axis 1", y="Axis 2")
  
  if(nrow(ar) > 0){
    p <- p +
      geom_segment(data=ar, aes(x=x0,y=y0,xend=x1,yend=y1),
                   inherit.aes=FALSE, linewidth=0.7,
                   arrow=grid::arrow(length=grid::unit(0.02,"npc"))) +
      ggrepel::geom_text_repel(data=ar, aes(x=x1,y=y1,label=var),
                               inherit.aes=FALSE, size=3, max.overlaps=Inf)
  }
  
  ggsave(fn_png, p, width=8, height=6, dpi=300)
  p
}

# -------- reconstruct meta_plot from model sample names ----------------------
make_meta_plot <- function(groups){
  tibble(Group = groups) %>%
    mutate(
      Base   = substr(Group, 1, 2),
      Season = substr(Group, 3, 3)
    ) %>%
    mutate(
      Season = factor(Season, levels=c("J","S")),
      Base   = factor(Base, levels=intersect(names(base_cols_map), unique(Base)))
    )
}

# -------- load matrices used in modeling (needed for envfit predictors) ------
# these files were saved by your big script
X_scaled <- readRDS(file.path(out_dir, "models", "X_scaled.rds"))  # if you didn't save it, see note below

# If you DON'T have X_scaled.rds, create it now from your geochem table:
# (tell me and I’ll give you a 10-line block)

# ===================== STRADA B PLOT (integrated RDA) =======================
objB <- readRDS(file.path(out_dir, "models/StradaB_model.rds"))
rda_mod <- objB$rda

groupsB <- rownames(scores(rda_mod, display="sites"))
metaB <- make_meta_plot(groupsB)

scB <- as.data.frame(scores(rda_mod, display="sites", choices=1:2))
colnames(scB)[1:2] <- c("Axis1","Axis2")
scB$Group <- rownames(scB)

# envfit uses X_scaled aligned
envB <- X_scaled[groupsB, , drop=FALSE]

plot_sites_envfit(scB, envB, metaB,
                  title="Integrated (PCoA axes) RDA — Group-merged (envfit arrows)",
                  fn_png=file.path(out_dir,"figures","B_integratedAxes_envfit.png"),
                  topN=6
)

# ===================== STRADA A PLOTS (4 layers) ============================
A_results <- readRDS(file.path(out_dir, "models/StradaA_models.rds"))

for(nm in names(A_results)){
  mod <- A_results[[nm]]$model
  groups <- rownames(scores(mod, display="sites"))
  metaA <- make_meta_plot(groups)
  
  sc <- as.data.frame(scores(mod, display="sites", choices=1:2))
  colnames(sc)[1:2] <- c("Axis1","Axis2")
  sc$Group <- rownames(sc)
  
  envA <- X_scaled[groups, , drop=FALSE]
  
  plot_sites_envfit(sc, envA, metaA,
                    title=paste0("dbRDA — ", nm, " (Group-merged; envfit arrows)"),
                    fn_png=file.path(out_dir,"figures",paste0("A_dbRDA_",nm,"_envfit.png")),
                    topN=6
  )
}

message("Saved plots in: ", normalizePath(file.path(out_dir, "figures")))

####################################
####################################
####################################

###############################################################################
# RUN 1: physchem (T, ph, EC)
# RUN 2: trace elements (auto-select top K by envfit r2 to avoid overfitting)
###############################################################################

pkgs <- c("dplyr","tibble","tidyr","readr","stringr","vegan","ape","ggplot2","ggrepel")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(to_install)) install.packages(to_install, repos="https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# ------------------------- EDIT PATHS ----------------------------------------
base_dir <- "integrated_multiomics_out"     # dove hai già salvato matrices/
geo_file <- "metadata4.csv"                # tuo metadata geochimico

mat_dir <- file.path(base_dir, "matrices")
out_dir <- file.path(base_dir, "redo_physchem_vs_trace")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# input matrices (Group x features)
fam_file  <- file.path(mat_dir, "Families_groupMerged.tsv")
path_file <- file.path(mat_dir, "PathwaysRel_groupMerged.tsv")
arg_file  <- file.path(mat_dir, "ARG_per16S_groupMerged.tsv")
mrg_file  <- file.path(mat_dir, "MRG_metalOnly_per16S_groupMerged.tsv")
# ---------------------------------------------------------------------------

# ------------------------- COLORS (Base) ------------------------------------
base_cols_map <- c(
  BF="#1b9e77", BR="#a6dba0", CF="#d95f02", CR="#fdb863",
  PF="#7570b3", PR="#b2abd2", SF="#1f9ac2", SR="#a6dce7",
  VF="#e7298a", VR="#f2b2d4"
)
# ---------------------------------------------------------------------------

# ------------------------- HELPERS ------------------------------------------
read_group_matrix <- function(path){
  x <- readr::read_tsv(path, show_col_types = FALSE)
  stopifnot(names(x)[1] %in% c("Group","group","Sample","sample"))
  rn <- x[[1]]
  mat <- x %>% select(-1) %>% as.data.frame()
  mat <- as.matrix(sapply(mat, function(v) suppressWarnings(as.numeric(v))))
  rownames(mat) <- rn
  mat[is.na(mat)] <- 0
  storage.mode(mat) <- "numeric"
  mat
}

to_num <- function(x){
  if(is.character(x)){
    x <- gsub(",", ".", x)
    x <- gsub("\\s*\\[.*?\\]", "", x)
    x <- gsub("<LOD", NA, x)
    x <- gsub("<LOD\\s*", NA, x)
  }
  suppressWarnings(as.numeric(x))
}

clean_geo_names <- function(nm){
  nm2 <- nm
  nm2 <- gsub("\n", " ", nm2)
  nm2 <- gsub("\"", "", nm2)
  
  nm2 <- gsub("^EC.*", "EC", nm2)
  nm2 <- gsub("^pH.*", "ph", nm2)
  nm2 <- gsub("^Temp.*", "T", nm2)
  nm2 <- gsub("^alt.*", "alt", nm2)
  
  # elements: take first token
  nm2 <- ifelse(nm2 %in% c("EC","ph","T","alt","Sample"), nm2,
                sub("^\\s*([A-Za-z]{1,2})\\b.*", "\\1", nm2))
  nm2
}

median_impute_scale <- function(df){
  for(nm in names(df)){
    if(is.numeric(df[[nm]])){
      med <- median(df[[nm]], na.rm=TRUE)
      df[[nm]][is.na(df[[nm]])] <- med
    }
  }
  as.data.frame(scale(df))
}

meta_from_group <- function(groups){
  tibble(Group = groups) %>%
    mutate(
      Base   = substr(Group, 1, 2),
      Season = substr(Group, 3, 3),
      Type   = substr(Group, 2, 2),
      Site   = substr(Group, 1, 1)
    ) %>%
    mutate(
      Season = factor(Season, levels = c("J","S")),
      Base   = factor(Base, levels = intersect(names(base_cols_map), unique(Base)))
    )
}

# envfit-based selection of top K predictors (by r2), optionally BH filter
select_topK_by_envfit <- function(Y_sxf, X_df, K=6){
  d <- vegan::vegdist(Y_sxf, method="bray")
  p <- ape::pcoa(d)
  ord <- as.data.frame(p$vectors[,1:2])
  colnames(ord) <- c("Axis1","Axis2")
  
  set.seed(123)
  ef <- vegan::envfit(as.matrix(ord), X_df, permutations=999)
  
  tbl <- tibble(
    var = names(ef$vectors$pvals),
    r2  = as.numeric(ef$vectors$r),
    p   = as.numeric(ef$vectors$pvals),
    p_adj = p.adjust(as.numeric(ef$vectors$pvals), "BH")
  ) %>% arrange(desc(r2))
  
  # take top K by r2 (stable, avoids overfitting)
  keep <- head(tbl$var, min(K, nrow(tbl)))
  list(keep=keep, table=tbl)
}

# capscale model + plot (arrows from envfit on CAP scores)
fit_and_plot_capscale <- function(Y_sxf, X_df, meta_plot, tag, layer, K_arrows=6){
  mod <- vegan::capscale(Y_sxf ~ ., data = X_df, distance = "bray")
  
  # tests
  glob <- as.data.frame(anova.cca(mod, permutations=999)) %>% rownames_to_column("term")
  term <- as.data.frame(anova.cca(mod, by="term", permutations=999)) %>% rownames_to_column("term") %>%
    mutate(p_adj = p.adjust(`Pr(>F)`, "BH"))
  
  # save tables
  dir.create(file.path(out_dir, tag, "tables"), recursive=TRUE, showWarnings=FALSE)
  write_tsv(glob, file.path(out_dir, tag, "tables", paste0("A_",layer,"_global.tsv")))
  write_tsv(term, file.path(out_dir, tag, "tables", paste0("A_",layer,"_terms_BH.tsv")))
  
  # sites scores
  sc <- as.data.frame(scores(mod, display="sites", choices=1:2))
  colnames(sc)[1:2] <- c("Axis1","Axis2")
  sc$Group <- rownames(sc)
  dfp <- sc %>% left_join(meta_plot, by="Group")
  
  # envfit arrows on CAP axes (more stable than bp if aliased)
  set.seed(123)
  ef <- envfit(as.matrix(dfp[,c("Axis1","Axis2")]), X_df, permutations=999)
  ef_tbl <- tibble(
    var = names(ef$vectors$pvals),
    r2  = as.numeric(ef$vectors$r),
    p   = as.numeric(ef$vectors$pvals),
    p_adj = p.adjust(as.numeric(ef$vectors$pvals), "BH")
  ) %>% arrange(p_adj)
  
  keep <- ef_tbl %>% filter(p_adj < 0.05) %>% pull(var)
  if(length(keep)==0) keep <- ef_tbl %>% arrange(desc(r2)) %>% head(K_arrows) %>% pull(var)
  
  ar <- as.data.frame(scores(ef, display="vectors"))
  ar$var <- rownames(ar)
  ar <- ar %>% filter(var %in% keep) %>% left_join(ef_tbl, by="var")
  
  # scale arrows
  xr <- diff(range(dfp$Axis1, na.rm=TRUE))
  yr <- diff(range(dfp$Axis2, na.rm=TRUE))
  mult <- 0.6 * min(xr, yr)
  ar <- ar %>% mutate(x0=0,y0=0, x1=Axis1*sqrt(r2)*mult, y1=Axis2*sqrt(r2)*mult)
  
  # plot
  dir.create(file.path(out_dir, tag, "figures"), recursive=TRUE, showWarnings=FALSE)
  
  p <- ggplot(dfp, aes(Axis1, Axis2, colour=Base, shape=Season)) +
    geom_point(size=4, alpha=0.95) +
    scale_colour_manual(values=base_cols_map, drop=FALSE) +
    theme_bw(base_size=13) +
    labs(title=paste0(tag," — dbRDA ",layer), x="dbRDA1", y="dbRDA2")
  
  if(nrow(ar)>0){
    p <- p +
      geom_segment(data=ar, aes(x=x0,y=y0,xend=x1,yend=y1),
                   inherit.aes=FALSE, linewidth=0.7,
                   arrow=grid::arrow(length=grid::unit(0.02,"npc"))) +
      ggrepel::geom_text_repel(data=ar, aes(x=x1,y=y1,label=var),
                               inherit.aes=FALSE, size=3, max.overlaps=Inf)
  }
  
  ggsave(file.path(out_dir, tag, "figures", paste0("A_dbRDA_",layer,".png")),
         p, width=8, height=6, dpi=300)
  
  # save envfit table
  write_tsv(ef_tbl, file.path(out_dir, tag, "tables", paste0("A_",layer,"_envfit_BH.tsv")))
  
  mod
}

# Strada B: build integrated axes + RDA + plot with envfit arrows
fit_and_plot_integrated_RDA <- function(fam, path, arg, mrg, X_df, meta_plot, tag, K_arrows=6){
  pcoa2 <- function(Y){
    d <- vegdist(Y, method="bray")
    p <- ape::pcoa(d)
    sc <- as.data.frame(p$vectors[,1:2])
    colnames(sc) <- c("PCoA1","PCoA2")
    sc$Group <- rownames(sc)
    sc
  }
  
  Y_axes <- pcoa2(fam) %>%
    rename(TAX1=PCoA1, TAX2=PCoA2) %>%
    left_join(pcoa2(path) %>% rename(FUN1=PCoA1, FUN2=PCoA2), by="Group") %>%
    left_join(pcoa2(arg)  %>% rename(ARG1=PCoA1, ARG2=PCoA2), by="Group") %>%
    left_join(pcoa2(mrg)  %>% rename(MRG1=PCoA1, MRG2=PCoA2), by="Group") %>%
    column_to_rownames("Group") %>%
    as.matrix()
  
  # RDA
  rda_mod <- rda(Y_axes ~ ., data=X_df)
  
  dir.create(file.path(out_dir, tag, "tables"), recursive=TRUE, showWarnings=FALSE)
  dir.create(file.path(out_dir, tag, "figures"), recursive=TRUE, showWarnings=FALSE)
  
  glob <- as.data.frame(anova.cca(rda_mod, permutations=999)) %>% rownames_to_column("term")
  term <- as.data.frame(anova.cca(rda_mod, by="term", permutations=999)) %>% rownames_to_column("term") %>%
    mutate(p_adj = p.adjust(`Pr(>F)`, "BH"))
  
  write_tsv(glob, file.path(out_dir, tag, "tables", "B_integrated_global.tsv"))
  write_tsv(term, file.path(out_dir, tag, "tables", "B_integrated_terms_BH.tsv"))
  
  sc <- as.data.frame(scores(rda_mod, display="sites", choices=1:2))
  colnames(sc)[1:2] <- c("Axis1","Axis2")
  sc$Group <- rownames(sc)
  dfp <- sc %>% left_join(meta_plot, by="Group")
  
  # envfit arrows (same rule)
  set.seed(123)
  ef <- envfit(as.matrix(dfp[,c("Axis1","Axis2")]), X_df, permutations=999)
  ef_tbl <- tibble(
    var = names(ef$vectors$pvals),
    r2  = as.numeric(ef$vectors$r),
    p   = as.numeric(ef$vectors$pvals),
    p_adj = p.adjust(as.numeric(ef$vectors$pvals), "BH")
  ) %>% arrange(p_adj)
  
  keep <- ef_tbl %>% filter(p_adj < 0.05) %>% pull(var)
  if(length(keep)==0) keep <- ef_tbl %>% arrange(desc(r2)) %>% head(K_arrows) %>% pull(var)
  
  ar <- as.data.frame(scores(ef, display="vectors"))
  ar$var <- rownames(ar)
  ar <- ar %>% filter(var %in% keep) %>% left_join(ef_tbl, by="var")
  
  xr <- diff(range(dfp$Axis1, na.rm=TRUE))
  yr <- diff(range(dfp$Axis2, na.rm=TRUE))
  mult <- 0.6 * min(xr, yr)
  ar <- ar %>% mutate(x0=0,y0=0, x1=Axis1*sqrt(r2)*mult, y1=Axis2*sqrt(r2)*mult)
  
  p <- ggplot(dfp, aes(Axis1, Axis2, colour=Base, shape=Season)) +
    geom_point(size=4, alpha=0.95) +
    scale_colour_manual(values=base_cols_map, drop=FALSE) +
    theme_bw(base_size=13) +
    labs(title=paste0(tag," — Integrated axes RDA"), x="RDA1", y="RDA2")
  
  if(nrow(ar)>0){
    p <- p +
      geom_segment(data=ar, aes(x=x0,y=y0,xend=x1,yend=y1),
                   inherit.aes=FALSE, linewidth=0.7,
                   arrow=grid::arrow(length=grid::unit(0.02,"npc"))) +
      ggrepel::geom_text_repel(data=ar, aes(x=x1,y=y1,label=var),
                               inherit.aes=FALSE, size=3, max.overlaps=Inf)
  }
  
  ggsave(file.path(out_dir, tag, "figures", "B_integrated_RDA.png"),
         p, width=8, height=6, dpi=300)
  
  write_tsv(ef_tbl, file.path(out_dir, tag, "tables", "B_integrated_envfit_BH.tsv"))
  
  rda_mod
}
# ---------------------------------------------------------------------------

# ------------------------- LOAD MATRICES ------------------------------------
stopifnot(file.exists(fam_file), file.exists(path_file), file.exists(arg_file), file.exists(mrg_file))

fam  <- read_group_matrix(fam_file)
path <- read_group_matrix(path_file)
arg  <- read_group_matrix(arg_file)
mrg  <- read_group_matrix(mrg_file)

# transforms
fam_h  <- decostand(fam,  "hellinger")
path_h <- decostand(path, "hellinger")
arg_l  <- log1p(arg)
mrg_l  <- log1p(mrg)

# ------------------------- LOAD GEOCHEM -------------------------------------
geo <- read_delim(geo_file, delim=",", show_col_types = FALSE)
names(geo) <- clean_geo_names(names(geo))
stopifnot("Sample" %in% names(geo))
geo <- geo %>% rename(Group = Sample)

# numeric conversion for all potential predictors
for(nm in setdiff(names(geo), "Group")) geo[[nm]] <- to_num(geo[[nm]])

# align samples
common <- Reduce(intersect, list(rownames(fam_h), rownames(path_h), rownames(arg_l), rownames(mrg_l), geo$Group))
if(length(common) < 8) stop("Troppi pochi Group comuni tra matrici e geochimica. Controlla i nomi (es. SFJ vs SFJ1).")

fam_h  <- fam_h[common, , drop=FALSE]
path_h <- path_h[common, , drop=FALSE]
arg_l  <- arg_l[common, , drop=FALSE]
mrg_l  <- mrg_l[common, , drop=FALSE]

meta_plot <- meta_from_group(common)

# ------------------------- DEFINE PREDICTOR SETS -----------------------------
physchem_vars <- c("T","ph","EC", "Altitude", "Mean δ2H" , "dev st δ2H" , "Mean δ18O" , "dev st δ18O", "d-excess")

trace_vars <- c("Li","Be","Al","Ti","V","Cr","Mn","Fe","Co","Ni","Zn","Cu","As","Se","Rb","Sr",
                "Mo","Cd","Sb","Ba","La","Ce","Pb","U","alt")

# keep only available columns
physchem_vars <- intersect(physchem_vars, names(geo))
trace_vars    <- intersect(trace_vars,    names(geo))

# ------------------------- BUILD X (impute+scale) ----------------------------
build_X <- function(vars){
  X <- geo %>% filter(Group %in% common) %>% arrange(match(Group, common)) %>% select(all_of(vars))
  X <- as.data.frame(X)
  Xs <- median_impute_scale(X)
  rownames(Xs) <- common
  Xs
}

# ------------------------- RUN ANALYSES --------------------------------------
run_all <- function(tag, Xs, Kselect=6){
  dir.create(file.path(out_dir, tag), recursive=TRUE, showWarnings=FALSE)
  
  # optional: if too many predictors, select top K by envfit r2 on families
  if(ncol(Xs) > Kselect){
    sel <- select_topK_by_envfit(fam_h, Xs, K=Kselect)
    Xuse <- Xs[, sel$keep, drop=FALSE]
    write_tsv(sel$table, file.path(out_dir, tag, "tables_envfit_rank_allPredictors.tsv"))
  } else {
    Xuse <- Xs
  }
  
  # Strada A (4 layers)
  A_mod_fam  <- fit_and_plot_capscale(fam_h,  Xuse, meta_plot, tag, "Families",   K_arrows=6)
  A_mod_path <- fit_and_plot_capscale(path_h, Xuse, meta_plot, tag, "Pathways",   K_arrows=6)
  A_mod_arg  <- fit_and_plot_capscale(arg_l,  Xuse, meta_plot, tag, "ARG_per16S", K_arrows=6)
  A_mod_mrg  <- fit_and_plot_capscale(mrg_l,  Xuse, meta_plot, tag, "MRG_per16S", K_arrows=6)
  
  saveRDS(list(Families=A_mod_fam, Pathways=A_mod_path, ARG=A_mod_arg, MRG=A_mod_mrg),
          file.path(out_dir, tag, "models_StradaA.rds"))
  
  # Strada B (integrated)
  B_mod <- fit_and_plot_integrated_RDA(fam_h, path_h, arg_l, mrg_l, Xuse, meta_plot, tag, K_arrows=6)
  saveRDS(B_mod, file.path(out_dir, tag, "model_StradaB.rds"))
  
  write_tsv(as.data.frame(Xuse) %>% rownames_to_column("Group"),
            file.path(out_dir, tag, "X_used_scaled.tsv"))
  
  message("DONE: ", tag, " -> ", normalizePath(file.path(out_dir, tag)))
}

# RUN 1: physchem only
X_phys <- build_X(physchem_vars)
run_all("RUN_physchem_T_ph_EC", X_phys, Kselect = 3)

# RUN 2: trace only (auto-select top 6 predictors if too many)
X_trace <- build_X(trace_vars)
run_all("RUN_traceElements_top6", X_trace, Kselect = 6)

message("All outputs in: ", normalizePath(out_dir))


###############################################################################
# Hydrologist package: summary table + 2 plots (Families + Geochemistry)
# Outputs:
#  - out_hydro/Summary_for_hydrologists.xlsx
#  - out_hydro/Fig1_scatter_trace_vs_family.png
#  - out_hydro/Fig2_S07_heatmap_features.png
###############################################################################

# ---------------------- PACKAGES ----------------------
pkgs <- c("dplyr","tibble","tidyr","readr","stringr","ggplot2","ggrepel","openxlsx","vegan")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(to_install)) install.packages(to_install, repos="https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# ---------------------- PATHS (EDIT) -------------------
# 1) Families matrix GROUP-MERGED (Group x families, relative abundance)
fam_file <- "integrated_multiomics_out/matrices/Families_groupMerged.tsv"

# 2) Geochemistry table
geo_file <- "metadata4.csv"

# 3) Output folder
out_dir <- "out_hydro"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# TARGET: how S07_22 is written in the geochem file ("Codice sorgente")
target_code <- "S07_22"   # <-- EDIT if needed
# -------------------------------------------------------

# ---------------------- HELPERS ------------------------
to_num <- function(x){
  if(is.character(x)){
    x <- gsub(",", ".", x)
    x <- gsub("\\s*\\[.*?\\]", "", x)  # remove brackets like "4.5 [4.3]"
    x <- gsub("<LOD", NA, x)
    x <- gsub("<LOD\\s*", NA, x)
    x <- trimws(x)
  }
  suppressWarnings(as.numeric(x))
}

clean_geo_names <- function(nm){
  nm2 <- nm
  nm2 <- gsub("\n", " ", nm2)
  nm2 <- gsub("\"", "", nm2)
  
  # keep these as-is (for mapping)
  nm2 <- gsub("^Codice sorgente.*", "Code", nm2)
  nm2 <- gsub("^Sample.*", "Group", nm2)
  nm2 <- gsub("^Area.*", "Area", nm2)
  nm2 <- gsub("^Note.*", "Note", nm2)
  
  # phys-chem
  nm2 <- gsub("^EC.*", "EC", nm2)
  nm2 <- gsub("^pH.*", "ph", nm2)
  nm2 <- gsub("^Temp.*", "T", nm2)
  nm2 <- gsub("^alt.*", "alt", nm2)
  
  # isotopes
  nm2 <- gsub("Mean.*δ2H.*", "d2H", nm2)
  nm2 <- gsub("dev st.*δ2H.*", "sd_d2H", nm2)
  nm2 <- gsub("Mean.*δ18O.*", "d18O", nm2)
  nm2 <- gsub("dev st.*δ18O.*", "sd_d18O", nm2)
  nm2 <- gsub("^d-?excess.*", "d_excess", nm2)
  
  # major ions
  nm2 <- gsub("^F-.*", "F", nm2)
  nm2 <- gsub("^Cl-.*", "Cl", nm2)
  nm2 <- gsub("^NO2-.*", "NO2", nm2)
  nm2 <- gsub("^Br-.*", "Br", nm2)
  nm2 <- gsub("^NO3-.*", "NO3", nm2)
  nm2 <- gsub("^PO43-.*", "PO4", nm2)
  nm2 <- gsub("^SO42-.*", "SO4", nm2)
  nm2 <- gsub("^HCO3-.*", "HCO3", nm2)
  
  # trace elements: take first token (Al, Cd, Pb...)
  nm2 <- ifelse(nm2 %in% c("Group","Code","Area","Note","EC","ph","T","alt",
                           "d2H","sd_d2H","d18O","sd_d18O","d_excess",
                           "F","Cl","NO2","Br","NO3","PO4","SO4","HCO3"),
                nm2,
                sub("^\\s*([A-Za-z]{1,2})\\b.*", "\\1", nm2))
  nm2
}

# robust Group metadata from code: BFJ / BFS etc.
add_group_meta <- function(df){
  df %>%
    mutate(
      Base   = substr(Group, 1, 2),
      Type   = substr(Group, 2, 2),
      Season = substr(Group, 3, 3),
      Site   = substr(Group, 1, 1)
    )
}

# z-score per feature (row)
z_row <- function(m){
  t(scale(t(m)))
}
# -------------------------------------------------------

# ---------------------- READ FAMILIES -------------------
fam <- readr::read_tsv(fam_file, show_col_types = FALSE)
stopifnot("Group" %in% names(fam))
fam_mat <- fam %>% column_to_rownames("Group") %>% as.matrix()
storage.mode(fam_mat) <- "numeric"
fam_mat[is.na(fam_mat)] <- 0

# alpha (Group-merged)
alpha_tbl <- tibble(Group = rownames(fam_mat)) %>%
  mutate(
    Observed   = rowSums(fam_mat > 0),
    Shannon    = vegan::diversity(fam_mat, index="shannon", base=2),
    InvSimpson = vegan::diversity(fam_mat, index="invsimpson"),
    Pielou     = vegan::diversity(fam_mat, index="shannon", base=2) /
      log2(pmax(rowSums(fam_mat > 0), 2))
  ) %>%
  add_group_meta()

# ---------------------- READ GEOCHEM --------------------
geo_raw <- readr::read_csv(geo_file, show_col_types = FALSE)
names(geo_raw) <- clean_geo_names(names(geo_raw))

# need Group column
stopifnot("Group" %in% names(geo_raw))

geo <- geo_raw %>%
  mutate(Group = as.character(Group))

# convert numeric columns
num_cols <- setdiff(names(geo), c("Group","Code","Area","Note"))
for(nm in intersect(num_cols, names(geo))){
  geo[[nm]] <- to_num(geo[[nm]])
}

# ---------------------- SELECT VARIABLES ----------------
# Phys-chem + isotopes + major ions (non-trace)
physchem_vars <- c("EC","ph","T","alt","d2H","sd_d2H","d18O","sd_d18O","d_excess",
                   "F","Cl","NO2","Br","NO3","PO4","SO4","HCO3")
physchem_vars <- intersect(physchem_vars, names(geo))

# Trace elements: keep those available
trace_vars <- c("Li","Be","Al","Ti","V","Cr","Mn","Fe","Co","Ni","Zn","Cu",
                "As","Se","Rb","Sr","Mo","Cd","Sb","Ba","La","Ce","Pb","U")
trace_vars <- intersect(trace_vars, names(geo))

# If trace list is huge, for the hydrologist summary we’ll pick top 6 by variance
trace_top <- trace_vars
if(length(trace_vars) > 6){
  vv <- sapply(geo[trace_vars], function(x) var(x, na.rm=TRUE))
  trace_top <- names(sort(vv, decreasing = TRUE))[1:6]
}

# ---------------------- JOIN SUMMARY --------------------
# Focus table: 1 row per Group, columns: meta + selected chemistry + alpha
chem_keep <- c(physchem_vars, trace_top)
sum_tbl <- alpha_tbl %>%
  left_join(geo %>% select(Group, any_of(c("Code","Area")), any_of(chem_keep)), by="Group")

# Top 10 families per Group (for Excel)
topN <- 10
topfam_long <- as.data.frame(fam_mat) %>%
  rownames_to_column("Group") %>%
  pivot_longer(-Group, names_to="Family", values_to="RelAbund") %>%
  group_by(Group) %>%
  arrange(desc(RelAbund), .by_group = TRUE) %>%
  slice_head(n = topN) %>%
  ungroup()

topfam_wide <- topfam_long %>%
  group_by(Group) %>%
  summarise(
    TopFamilies = paste0(Family, " (", round(RelAbund,3), ")", collapse = "; "),
    .groups="drop"
  )

sum_tbl2 <- sum_tbl %>% left_join(topfam_wide, by="Group")

# ---------------------- TARGET GROUPS (S07_22) ----------
# try to map S07_22 in geochem "Code" -> Group codes
target_groups <- geo %>%
  filter(!is.na(Code) & Code == target_code) %>%
  pull(Group) %>% unique()

# fallback: if not found, user can manually set it
if(length(target_groups) == 0){
  message("WARNING: target_code not found in geo$Code. ",
          "Set target_groups manually (e.g., c('SFJ','SFS')).")
  target_groups <- character(0)
}

# ---------------------- FIG 1: Scatter (trace vs family) --------------------
# choose 1 trace element (highest variance) + 2 families (most enriched in target)
# if target unknown, choose top variance trace + top abundant families.
choose_trace <- if(length(trace_top) >= 1) trace_top[1] else NA_character_

# compute enrichment if we have target
fam_means <- colMeans(fam_mat, na.rm=TRUE)
top_fams_global <- names(sort(fam_means, decreasing=TRUE))[1:15]

if(length(target_groups) > 0){
  tg <- intersect(target_groups, rownames(fam_mat))
  other <- setdiff(rownames(fam_mat), tg)
  if(length(tg) > 0 && length(other) > 0){
    enr <- colMeans(fam_mat[tg, , drop=FALSE]) - colMeans(fam_mat[other, , drop=FALSE])
    top_fams <- names(sort(enr, decreasing=TRUE))[1:2]
  } else {
    top_fams <- top_fams_global[1:2]
  }
} else {
  top_fams <- top_fams_global[1:2]
}

plot_df <- sum_tbl %>%
  select(Group, Base, Type, Season, any_of(c(choose_trace))) %>%
  left_join(
    as.data.frame(fam_mat[, top_fams, drop=FALSE]) %>% rownames_to_column("Group"),
    by="Group"
  ) %>%
  pivot_longer(cols = all_of(top_fams), names_to="Family", values_to="RelAbund") %>%
  mutate(is_target = Group %in% target_groups)

if(!is.na(choose_trace)){
  p1 <- ggplot(plot_df, aes_string(x=choose_trace, y="RelAbund")) +
    geom_point(aes(shape=Season, colour=Type), size=3, alpha=0.9) +
    facet_wrap(~Family, scales="free_y") +
    ggrepel::geom_text_repel(
      data = plot_df %>% filter(is_target),
      aes(label = Group),
      size = 3,
      max.overlaps = Inf
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = paste0("Trace element vs selected families (", target_code, " highlighted)"),
      x = choose_trace,
      y = "Relative abundance (family)"
    )
  
  ggsave(file.path(out_dir, "Fig1_scatter_trace_vs_family.png"),
         p1, width=10, height=5, dpi=300)
}

# ---------------------- FIG 2: Heatmap (target vs all) ----------------------
# Build feature matrix: selected trace + top 15 families
hm_fams <- top_fams_global[1:15]
hm_feats <- c(trace_top, hm_fams)

# feature matrix: rows = features, cols = groups
hm_mat <- matrix(NA_real_, nrow=length(hm_feats), ncol=nrow(sum_tbl),
                 dimnames=list(hm_feats, sum_tbl$Group))

# fill trace from geo
for(v in trace_top){
  hm_mat[v, sum_tbl$Group] <- sum_tbl[[v]]
}
# fill families from fam_mat
for(f in hm_fams){
  hm_mat[f, sum_tbl$Group] <- fam_mat[sum_tbl$Group, f]
}

# z-score per row (feature)
hm_z <- z_row(hm_mat)
hm_z[is.na(hm_z)] <- 0

hm_long <- as.data.frame(hm_z) %>%
  rownames_to_column("Feature") %>%
  pivot_longer(-Feature, names_to="Group", values_to="Z") %>%
  left_join(add_group_meta(tibble(Group=unique(sum_tbl$Group))), by="Group") %>%
  mutate(is_target = Group %in% target_groups)

# order columns: target first if present
ord_groups <- unique(c(target_groups, setdiff(unique(hm_long$Group), target_groups)))
hm_long$Group <- factor(hm_long$Group, levels = ord_groups)

p2 <- ggplot(hm_long, aes(x=Group, y=Feature, fill=Z)) +
  geom_tile() +
  theme_bw(base_size=11) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(
    title = paste0("S07_22-focused heatmap (z-score by feature)"),
    x = NULL, y = NULL
  )

ggsave(file.path(out_dir, "Fig2_S07_heatmap_features.png"),
       p2, width=11, height=6, dpi=300)

# ---------------------- EXCEL SUMMARY -------------------
wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Overview")
openxlsx::writeData(wb, "Overview", sum_tbl2)

openxlsx::addWorksheet(wb, "TopFamilies_long")
openxlsx::writeData(wb, "TopFamilies_long", topfam_long)

openxlsx::addWorksheet(wb, "Chem_vars_used")
openxlsx::writeData(wb, "Chem_vars_used",
                    tibble(
                      physchem_vars = physchem_vars,
                      trace_vars_available = trace_vars,
                      trace_vars_used_in_summary = trace_top,
                      target_code = target_code,
                      target_groups = paste(target_groups, collapse=", ")
                    ))

openxlsx::saveWorkbook(wb, file.path(out_dir, "Summary_for_hydrologists.xlsx"), overwrite = TRUE)

cat("DONE.\nSaved in: ", normalizePath(out_dir), "\n", sep="")


##################################################################
###################################################################
########### PLSM ###############################################
library(readr); library(dplyr); library(vegan); library(tibble)

fam <- read_tsv("integrated_multiomics_out/matrices/Families_groupMerged.tsv", show_col_types = FALSE)
X <- as.matrix(fam %>% select(-Group))
rownames(X) <- fam$Group

d <- vegdist(X, method = "bray")
sc <- cmdscale(d, k = 5, eig = TRUE)

coord <- as.data.frame(sc$points) %>%
  rownames_to_column("Group") %>%
  rename(PCoA1=V1, PCoA2=V2, PCoA3=V3, PCoA4=V4, PCoA5=V5)

var <- tibble(
  Axis = paste0("PCoA", seq_along(sc$eig)),
  Eigenvalue = sc$eig,
  Variance_explained = sc$eig / sum(sc$eig),
  Variance_explained_pct = 100 * sc$eig / sum(sc$eig)
)

write_tsv(coord, "PCoA_Families_groupMerged.tsv")
write_tsv(var,   "PCoA_Families_groupMerged_variance.tsv")
write.table(as.matrix(d), "BrayCurtis_Families_groupMerged.tsv", sep="\t", quote=FALSE)

###############################################################################
# PLS-PM: Environment + Taxonomy + ARG + MRG (Group-level)
# - metadata3.txt: Sample-level env/metals (replicates)
# - Families_groupMerged.tsv: family composition by Group (rel. abundance)
# - PCoA_Families_groupMerged.tsv: ordination coords by Group
# - MRG_metalOnly_per16S_groupMerged.tsv: metal MRG by Group
# - ARG matrix (NOT groupMerged): features x Sample OR Sample x features
###############################################################################

#### 0) Packages --------------------------------------------------------------
pkgs <- c("dplyr","readr","tidyr","stringr","tibble","vegan","plspm",
          "DiagrammeR","DiagrammeRsvg","rsvg")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(to_install)) install.packages(to_install, repos="https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))


dir.create("integrated_multiomics_out/tables", showWarnings = FALSE)
dir.create("integrated_multiomics_out/figures", showWarnings = FALSE)

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(tibble)
library(vegan)
library(plspm)

# ---- FILES ----
meta_file <- "metadata3.txt"
fam_file  <- "integrated_multiomics_out/matrices/Families_groupMerged.tsv"
pcoa_file <- "PCoA_Families_groupMerged.tsv"
mrg_file  <- "integrated_multiomics_out/matrices/MRG_metalOnly_per16S_groupMerged.tsv"
mrg  <- "integrated_multiomics_out/matrices/MRG_metalOnly_per16S_groupMerged.tsv"
arg_file  <- "ARG_per16S_T1T2_renamed.tsv"   # <<--- metti qui quello giusto

# se per qualche motivo la prima colonna non si chiama esattamente Group (es. "roup"), forzala:
names(mrg)[1] <- "Group"

mrg <- mrg %>%
  mutate(Group = str_trim(as.character(Group))) %>%
  mutate(across(-Group, ~ as.numeric(str_replace_all(as.character(.x), ",", "."))))

mrg_metrics <- mrg %>%
  rowwise() %>%
  mutate(
    MRG_total = sum(c_across(-Group), na.rm = TRUE),
    MRG_rich  = sum(c_across(-Group) > 0, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(Group, MRG_total, MRG_rich)

# check: deve essere 20 righe e zero NA
print(dim(mrg_metrics))
print(colSums(is.na(mrg_metrics)))




# ---- helper numerico robusto ----
to_num <- function(x) as.numeric(str_replace_all(as.character(x), ",", "."))

# ---- 1) metadata3 -> env_group (Group median) ----
meta3 <- read_tsv(meta_file, show_col_types = FALSE) %>%
  mutate(
    Sample = as.character(Sample),
    Group  = str_trim(str_sub(Sample, 1, 3))
  )

env_vars   <- intersect(c("T","ph","EC","alt"), names(meta3))
metal_vars <- intersect(c("Cd","Sb","Pb"), names(meta3))

env_group <- meta3 %>%
  mutate(across(all_of(c(env_vars, metal_vars)), to_num)) %>%
  group_by(Group) %>%
  summarise(
    across(all_of(env_vars),   ~ median(.x, na.rm = TRUE)),
    across(all_of(metal_vars), ~ median(.x, na.rm = TRUE)),
    .groups="drop"
  )

# ---- 2) taxonomy diversity from Families_groupMerged ----
fam <- read_tsv(fam_file, show_col_types = FALSE) %>%
  mutate(Group = str_trim(as.character(Group)))

Xfam <- fam %>% select(-Group) %>% mutate(across(everything(), to_num)) %>% as.data.frame()
rownames(Xfam) <- fam$Group

tax_div <- tibble(
  Group    = fam$Group,
  Shannon  = vegan::diversity(Xfam, index="shannon"),
  Richness = vegan::specnumber(Xfam)
)

# ---- 3) PCoA coords ----
pcoa <- read_tsv(pcoa_file, show_col_types = FALSE) %>%
  transmute(
    Group = str_trim(as.character(Group)),
    PCoA1 = to_num(PCoA1),
    PCoA2 = to_num(PCoA2)
  )



# ---- 5) ARG metrics from NOT-merged matrix (replicate-level) ----
arg_raw <- read_tsv(arg_file, show_col_types = FALSE)

if ("Sample" %in% names(arg_raw)) {
  # Sample x features
  arg_long <- arg_raw %>%
    mutate(
      Sample = as.character(Sample),
      Group  = str_trim(str_sub(Sample, 1, 3))
    ) %>%
    pivot_longer(cols = -c(Sample, Group), names_to="Feature", values_to="Abundance") %>%
    mutate(Abundance = to_num(Abundance), Abundance = ifelse(is.na(Abundance), 0, Abundance))
} else {
  # features x Sample
  feature_col <- names(arg_raw)[1]
  sample_cols <- setdiff(names(arg_raw), feature_col)
  
  arg_long <- arg_raw %>%
    pivot_longer(cols = all_of(sample_cols), names_to="Sample", values_to="Abundance") %>%
    mutate(
      Sample = as.character(Sample),
      Group  = str_trim(str_sub(Sample, 1, 3)),
      Abundance = to_num(Abundance),
      Abundance = ifelse(is.na(Abundance), 0, Abundance)
    )
}

arg_sample <- arg_long %>%
  group_by(Sample, Group) %>%
  summarise(
    ARG_total = sum(Abundance, na.rm=TRUE),
    ARG_rich  = sum(Abundance > 0, na.rm=TRUE),
    .groups="drop"
  )

arg_group <- arg_sample %>%
  group_by(Group) %>%
  summarise(
    ARG_total = median(ARG_total, na.rm=TRUE),
    ARG_rich  = median(ARG_rich,  na.rm=TRUE),
    .groups="drop"
  )

# ---- 6) keep only groups present everywhere (EVITA colonne tutte NA) ----
groups_ok <- Reduce(intersect, list(env_group$Group, tax_div$Group, pcoa$Group, mrg_metrics$Group, arg_group$Group))
cat("Groups kept for PLS:", length(groups_ok), "\n")

env_group   <- env_group   %>% filter(Group %in% groups_ok)
tax_div     <- tax_div     %>% filter(Group %in% groups_ok)
pcoa        <- pcoa        %>% filter(Group %in% groups_ok)
mrg_metrics <- mrg_metrics %>% filter(Group %in% groups_ok)
arg_group   <- arg_group   %>% filter(Group %in% groups_ok)

df_sem <- env_group %>%
  inner_join(tax_div, by="Group") %>%
  inner_join(pcoa, by="Group") %>%
  inner_join(arg_group, by="Group") %>%
  inner_join(mrg_metrics, by="Group")

# ---- 7) transforms + numeric sanity ----
df_sem <- df_sem %>%
  mutate(
    across(all_of(metal_vars), ~ log1p(pmax(.x, 0))),  # sicurezza: no negativi
    ARG_total = log1p(pmax(ARG_total, 0)),
    MRG_total = log1p(pmax(MRG_total, 0)),
    EC = if ("EC" %in% names(.)) log1p(pmax(EC, 0)) else EC
  )

# ---- 8) blocks (>=2 indicatori se possibile) ----
blocks <- list(
  Physical  = intersect(c("T","alt"), names(df_sem)),
  WaterChem = intersect(c("ph","EC"), names(df_sem)),
  Metals    = intersect(c("Cd","Sb","Pb","Zn","Cu","Ni"), names(df_sem)),
  MicroComp = c("PCoA1","PCoA2"),
  MicroDiv  = c("Shannon","Richness"),
  ARGs      = c("ARG_total","ARG_rich"),
  MRGs      = c("MRG_total","MRG_rich")
)

# rimuovi blocchi vuoti
blocks2 <- blocks[sapply(blocks, length) > 0]

# ---- 9) remove bad manifest vars (NA all / sd=0) then scale ----
raw <- df_sem %>% select(-Group) %>% mutate(across(everything(), to_num))

diag <- tibble(
  var = names(raw),
  n_na = sapply(raw, function(x) sum(is.na(x))),
  sd   = sapply(raw, function(x) sd(x, na.rm=TRUE)),
  n_unique = sapply(raw, function(x) n_distinct(x, na.rm=TRUE))
)

keep_vars <- diag %>%
  filter(n_na < nrow(raw), is.finite(sd), sd > 0, n_unique >= 2) %>%
  pull(var)

raw2 <- raw %>% select(all_of(keep_vars))
dat2 <- as.data.frame(scale(raw2))
stopifnot(all(is.finite(as.matrix(dat2))))

# aggiorna blocks2
blocks2 <- lapply(blocks2, function(v) intersect(v, keep_vars))
blocks2 <- blocks2[sapply(blocks2, length) > 0]
print(sapply(blocks2, length))

# ---- 10) inner model (lower triangular) and trim to remaining LVs ----
# ordine
LV_order <- c("Physical","WaterChem","Metals","MicroComp","MicroDiv","ARGs","MRGs")
LVs <- intersect(LV_order, names(blocks2))
blocks2 <- blocks2[LVs]

inner <- matrix(0, nrow=length(LVs), ncol=length(LVs),
                dimnames=list(LVs, LVs))

add_path <- function(from, to) if(all(c(from,to) %in% LVs)) inner[to, from] <<- 1

# env -> microbiome
add_path("Physical",  "MicroComp")
add_path("WaterChem", "MicroComp")
add_path("Physical",  "MicroDiv")
add_path("WaterChem", "MicroDiv")

# microbiome -> ARG
add_path("MicroComp", "ARGs")
add_path("MicroDiv",  "ARGs")

# metals + ARG -> MRG (co-selection)
add_path("Metals", "MRGs")
add_path("ARGs",   "MRGs")

# opzionale (se vuoi): composition -> MRG
add_path("MicroComp","MRGs")

stopifnot(all(inner[upper.tri(inner)] == 0))

# ---- 11) fit ----
modes <- rep("A", length(blocks2))

set.seed(1)
pls0 <- plspm(dat2, inner, blocks2, modes=modes, boot.val=FALSE)

# se gira:
set.seed(1)
pls <- plspm(dat2, inner, blocks2, modes=modes, boot.val=TRUE, br=2000)
# quick outputs
write.table(pls$path_coefs, "tables1/PLSPM_path_coefs.txt", sep="\t", quote=FALSE)
write.table(pls$outer_model, "tables1/PLSPM_outer_loadings.txt", sep="\t", quote=FALSE, row.names=FALSE)
writeLines(sprintf("Goodness of Fit: %.3f", pls$gof), "tables/PLSPM_GOF.txt")

# salva dataset finale usato
df_used <- df_sem %>% select(Group, all_of(keep_vars))
write_tsv(df_used, "tables1/PLSPM_input_USED.tsv")
write.table(df_used, "tables1/PLSPM_input_USED.txt", sep="\t", quote=FALSE, row.names=FALSE)

##################Estrai edges + bootstrap (CI e significatività 95%)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)

# (A) path coefficients original (non-zero)
pc <- as.data.frame(pls$path_coefs) %>%
  rownames_to_column("from") %>%
  pivot_longer(-from, names_to = "to", values_to = "beta") %>%
  filter(beta != 0)

# (B) bootstrap paths: from rownames "A -> B"
boot_df <- as.data.frame(pls$boot$paths) %>%
  rownames_to_column("path") %>%
  separate(path, into = c("from","to"), sep = "\\s*->\\s*", remove = FALSE)

# (C) merge + sig
edges <- pc %>%
  left_join(boot_df, by = c("from","to")) %>%
  mutate(
    sig   = !(perc.025 <= 0 & perc.975 >= 0),
    style = ifelse(sig, "solid", "dashed"),
    col   = ifelse(beta >= 0, "#D55E00", "#0072B2"),
    label = sprintf("%.3f", beta)
  )

write_tsv(edges, "tables1/PLSPM_edges_with_bootstrap.tsv")


library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(tibble)
library(dplyr)

LVs <- names(blocks2)

node_pos <- tibble::tribble(
  ~id,         ~label,                          ~x, ~y,
  "Physical",   "Physical\nparameters",          1,  3,
  "WaterChem",  "Water\nchemistry",              3,  3,
  "Metals",     "Metals\n(PC1–PC2)",             5,  3,
  "MicroComp",  "Microbial\ncomposition",        3,  2,
  "MicroDiv",   "Microbial\ndiversity",          3,  1,
  "ARGs",       "ARGs",                          2,  0,
  "MRGs",       "MRGs",                          4,  0
) %>% filter(id %in% LVs)

node_lines <- apply(node_pos, 1, function(r){
  sprintf('"%s" [label="%s", shape=box, style="rounded,filled", fillcolor="white", pos="%s,%s!"];',
          r[["id"]], r[["label"]], r[["x"]], r[["y"]])
})

edge_lines <- apply(edges, 1, function(r){
  sprintf('"%s" -> "%s" [label="%s", color="%s", fontcolor="%s", style="%s", penwidth=2];',
          r[["from"]], r[["to"]], r[["label"]], r[["col"]], r[["col"]], r[["style"]])
})

dot <- paste0(
  "digraph G {\n",
  'graph [layout=neato, overlap=false, splines=true];\n',
  'node [fontname=Helvetica, fontsize=16];\n',
  'edge [fontname=Helvetica, fontsize=14];\n',
  paste(node_lines, collapse="\n"), "\n",
  paste(edge_lines, collapse="\n"), "\n",
  sprintf('gof [label="Goodness of Fit: %.3f", shape=plaintext, fontsize=14, pos="3,-0.8!"];\n', pls$gof),
  "}\n"
)

writeLines(dot, "tables1/PLSPM_graph.dot")

gr <- DiagrammeR::grViz(dot)

DiagrammeRsvg::export_svg(gr) %>%
  charToRaw() %>%
  rsvg::rsvg_png("figures1/PLSPM_path_model.png", width = 1600)



################## change direction 
# Order (must be causal order)
LV_order <- c("Physical","WaterChem","Metals","MicroComp","MicroDiv","MRGs","ARGs")
LVs <- intersect(LV_order, names(blocks2))
blocks2 <- blocks2[LVs]

inner <- matrix(0, nrow=length(LVs), ncol=length(LVs), dimnames=list(LVs, LVs))
add_path <- function(from, to) if(all(c(from,to) %in% LVs)) inner[to, from] <<- 1

# Physical/Water/Metals -> everything central + resistomes
add_path("Physical","MicroComp")
add_path("Physical","MicroDiv")
add_path("Physical","ARGs")
add_path("Physical","MRGs")

add_path("WaterChem","MicroComp")
add_path("WaterChem","MicroDiv")
add_path("WaterChem","ARGs")
add_path("WaterChem","MRGs")

add_path("Metals","MicroComp")
add_path("Metals","MicroDiv")
add_path("Metals","ARGs")
add_path("Metals","MRGs")

# MRGs -> ARGs (your requested direction)
add_path("MRGs","ARGs")

# MUST be lower triangular
stopifnot(all(inner[upper.tri(inner)] == 0))

modes <- rep("A", length(blocks2))

set.seed(1)
pls0 <- plspm(dat2, inner, blocks2, modes=modes, boot.val=FALSE)

set.seed(1)
pls <- plspm(dat2, inner, blocks2, modes=modes, boot.val=TRUE, br=2000)



# quick outputs
write.table(pls$path_coefs, "tables1/PLSPM_path_coefs1.txt", sep="\t", quote=FALSE)
write.table(pls$outer_model, "tables1/PLSPM_outer_loadings1.txt", sep="\t", quote=FALSE, row.names=FALSE)
writeLines(sprintf("Goodness of Fit: %.3f", pls$gof), "tables1/PLSPM_GOF1.txt")

# salva dataset finale usato
df_used <- df_sem %>% select(Group, all_of(keep_vars))
write_tsv(df_used, "tables/PLSPM_input_USED1.tsv")
write.table(df_used, "tables/PLSPM_input_USED1.txt", sep="\t", quote=FALSE, row.names=FALSE)

library(dplyr)
library(tidyr)
library(tibble)
library(readr)

# path coefs original (non-zero)
pc <- as.data.frame(pls$path_coefs) %>%
  rownames_to_column("from") %>%
  pivot_longer(-from, names_to="to", values_to="beta") %>%
  filter(beta != 0)

# bootstrap table: rownames like "A -> B"
boot_df <- as.data.frame(pls$boot$paths) %>%
  rownames_to_column("path") %>%
  separate(path, into=c("from","to"), sep="\\s*->\\s*", remove=FALSE)

edges <- pc %>%
  left_join(boot_df, by=c("from","to")) %>%
  mutate(
    # significance from bootstrap CI
    sig = !(perc.025 <= 0 & perc.975 >= 0),
    
    # approximate p-value from boot SE (useful for stars)
    z = abs(Original) / Std.Error,
    p_approx = 2 * (1 - pnorm(z)),
    
    stars = case_when(
      is.na(p_approx) ~ "",
      p_approx < 0.001 ~ "***",
      p_approx < 0.01  ~ "**",
      TRUE ~ ""
    ),
    
    style = ifelse(sig, "solid", "dashed"),   # dashed = not significant
    col   = ifelse(beta >= 0, "#D55E00", "#0072B2"),
    # label = coef + stars
    label = paste0(sprintf("%.3f", beta), stars)
  )

write_tsv(edges, "tables1/PLSPM_edges_with_stars.tsv")


# install.packages(c("tidygraph","ggraph","ggrepel","dplyr","tibble","grid"))
library(dplyr)
library(tibble)
library(tidygraph)
library(ggraph)
library(ggrepel)
library(grid)

# ---- 1) Prepara edges per il plot ----
# Atteso: edges con colonne from, to, beta, sig (TRUE/FALSE), p_approx, label, col, style
# Se nel tuo edges hai "style" = solid/dashed, lo traduciamo in linetype.

edges_plot <- edges %>%
  mutate(
    linetype = ifelse(sig, "solid", "dotted"),
    edge_col = ifelse(sig, col, "grey70"),
    # etichetta: coeff + stars (già in label) — se preferisci togliere i non-sig:
    lab = label
  )

# ---- 2) Nodi con posizioni (modifica se vuoi) ----
nodes <- tibble::tribble(
  ~name,        ~label,                          ~x, ~y,
  "Physical",    "Physical\nparameters",          1,  3,
  "WaterChem",   "Water\nchemistry",              3,  3,
  "Metals",      "Metals\n(PC1–PC2)",             5,  3,
  "MicroComp",   "Microbial\ncomposition",        3,  2,
  "MicroDiv",    "Microbial\ndiversity",          3,  1,
  "ARGs",        "ARGs",                          2,  0,
  "MRGs",        "MRGs",                          4,  0
)

# tieni solo nodi usati davvero
nodes <- nodes %>% filter(name %in% unique(c(edges_plot$from, edges_plot$to)))

# ---- 3) Costruisci grafo ----
g <- tbl_graph(nodes = nodes, edges = edges_plot, directed = TRUE)

# ---- 4) Per i "dot" all'origine: coord del nodo sorgente di ogni edge ----
dot_df <- edges_plot %>%
  left_join(nodes %>% select(name, x, y), by = c("from" = "name")) %>%
  mutate(dx = x, dy = y)

# ---- 5) Posizione label (midpoint) + repel ----
lab_df <- edges_plot %>%
  left_join(nodes %>% select(name, x, y), by = c("from" = "name")) %>%
  rename(x_from = x, y_from = y) %>%
  left_join(nodes %>% select(name, x, y), by = c("to" = "name")) %>%
  rename(x_to = x, y_to = y) %>%
  mutate(
    lx = (x_from + x_to) / 2,
    ly = (y_from + y_to) / 2
  )

# ---- 6) Plot ----
library(dplyr)
library(tidyr)
library(tibble)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# -----------------------------
# 1) EDGES: direction + dotted + stars
# -----------------------------
boot_df <- as.data.frame(pls$boot$paths) %>%
  rownames_to_column("path") %>%
  separate(path, into = c("from","to"), sep = "\\s*->\\s*", remove = FALSE)

if ("Std. Error" %in% names(boot_df)) boot_df <- boot_df %>% rename(Std.Error = `Std. Error`)

# p approx (per le stelline)
boot_df <- boot_df %>%
  mutate(
    z = abs(Original) / Std.Error,
    p_approx = 2 * (1 - pnorm(z)),
    sig = p_approx < 0.05,
    stars = case_when(
      p_approx < 0.01 ~ "**",
      p_approx < 0.05 ~ "*",
      TRUE ~ ""
    ),
    # colore in base al segno del coefficiente
    col = ifelse(Original >= 0, "#D55E00", "#0072B2"),
    # dotted se non significativo
    style = ifelse(sig, "solid", "dotted"),
    penw  = ifelse(sig, 2.6, 1.2),
    # label SOLO per significativi (molto più leggibile)
    label = ifelse(sig, paste0(sprintf("%.3f", Original), stars), "")
  ) %>%
  select(from, to, Original, p_approx, sig, stars, col, style, penw, label)

edges <- boot_df

write_tsv(edges, "tables1/PLSPM_edges_for_plot.tsv")

# -----------------------------
# 2) NODES (stesso layout tuo)
# -----------------------------
LVs <- names(blocks2)

node_pos <- tibble::tribble(
  ~id,         ~label,                          ~x, ~y,
  "Physical",   "Physical\nparameters",          0,  6,
  "WaterChem",  "Water\nchemistry",              6,  6,
  "Metals",     "Metals\n(PC1–PC2)",            12,  6,
  
  "MicroComp",  "Microbial\ncomposition",        3,  4,
  "MicroDiv",   "Microbial\ndiversity",          9,  4,
  
  "ARGs",       "ARGs",                          3,  1,
  "MRGs",       "MRGs",                          9,  1
) %>% filter(id %in% LVs)

node_lines <- apply(node_pos, 1, function(r){
  sprintf('"%s" [label="%s", shape=box, style="rounded,filled", fillcolor="white", pos="%s,%s!"];',
          r[["id"]], r[["label"]], r[["x"]], r[["y"]])
})

# tieni solo non-sig con |beta| >= 0.2 (esempio)
edges2 <- edges2 %>% filter(sig | abs(Original) >= 0.2)

# -----------------------------
# 3) EDGES lines: dotted if non-sig, stars in label if sig
# -----------------------------
# ricrea edges2 dal dataframe corretto
edges2 <- edges %>% filter(sig | abs(Original) >= 0.2)

edge_lines <- apply(edges2, 1, function(r){
  sprintf(
    '"%s" -> "%s" [label="%s", color="%s", fontcolor="%s", style="%s", penwidth=%s, arrowsize=0.8];',
    r[["from"]], r[["to"]], r[["label"]],
    r[["col"]], r[["col"]], r[["style"]],
    r[["penw"]]
  )
})

dot <- paste0(
  "digraph G {\n",
  'graph [layout=neato, overlap=false, splines=true];\n',
  'node [fontname=Helvetica, fontsize=16];\n',
  'edge [fontname=Helvetica, fontsize=14];\n',
  paste(node_lines, collapse="\n"), "\n",
  paste(edge_lines, collapse="\n"), "\n",
  sprintf('gof [label="Goodness of Fit: %.3f", shape=plaintext, fontsize=14, pos="3,-0.8!"];\n', pls$gof),
  "}\n"
)

writeLines(dot, "tables1/PLSPM_graph.dot")

gr <- DiagrammeR::grViz(dot)

DiagrammeRsvg::export_svg(gr) %>%
  charToRaw() %>%
  rsvg::rsvg_png("figures1/PLSPM_path_model_clean_dotted_stars2.png", width = 1600)

# --- EDGES (Graphviz) ---
edge_lines <- apply(edges2, 1, function(r){
  
  # ports: defaults
  tailport <- "s"
  headport <- "n"
  
  # make top -> middle enter from top
  if (r[["from"]] %in% c("Physical","WaterChem","Metals")) tailport <- "s"
  if (r[["to"]]   %in% c("MicroComp","MicroDiv"))          headport <- "n"
  
  # make middle -> bottom enter from top of bottom boxes
  if (r[["to"]]   %in% c("ARGs","MRGs"))                   headport <- "n"
  
  # labels: push away from edge
  ld <- 2.7
  la <- 18
  
  # special case: MRGs -> ARGs label (often overlaps)
  if (r[["from"]] == "MRGs" && r[["to"]] == "ARGs") { ld <- 3.4; la <- -25 }
  
  sprintf(
    '"%s" -> "%s" [label="%s", color="%s", fontcolor="%s", style="%s", penwidth=%s, arrowsize=%s, tailport=%s, headport=%s, labeldistance=%s, labelangle=%s];',
    r[["from"]], r[["to"]], r[["label"]],
    r[["col"]],  r[["col"]], r[["style"]],
    r[["penw"]], r[["arrowsize"]],
    tailport, headport, ld, la
  )
})

dot <- paste0(
  "digraph G {\n",
  'graph [layout=neato, overlap=false, splines=curved, sep="+30", pad="0.35"];\n',
  'node [fontname=Helvetica, fontsize=18, margin="0.2,0.2"];\n',
  'edge [fontname=Helvetica, fontsize=10, labelfontname=Helvetica, labelfontsize=10];\n',
  paste(node_lines, collapse="\n"), "\n",
  paste(edge_lines, collapse="\n"), "\n",
  sprintf('gof [label="Goodness of Fit: %.3f", shape=plaintext, fontsize=30, pos="6,6"];\n', pls$gof),
  "}\n"
)

writeLines(dot, "tables1/PLSPM_graph.dot")

gr <- DiagrammeR::grViz(dot)

DiagrammeRsvg::export_svg(gr) %>%
  charToRaw() %>%
  rsvg::rsvg_png("figures1/PLSPM_path_model_clean.png", width = 1600)


library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

dir.create("tables", showWarnings = FALSE)

# ==========
# INPUT
# ==========
infile <- "tables/PLSPM_input_USED1.tsv"   # <-- cambia solo questo se serve
df <- read_tsv(infile, show_col_types = FALSE)

# helper numerico robusto (virgole -> punto)
to_num <- function(x) as.numeric(str_replace_all(as.character(x), ",", "."))

# ==========
# TABLE S8: lista variabili del PLS-PM
# ==========

# metalli usati per PCA (metti qui quelli che vuoi dichiarare come input PCA)
metals_declared <- c("Cd","Sb","Pb","Zn","Cu","Ni")
metals_used <- intersect(metals_declared, names(df))

# definisci i blocchi come nel tuo modello
blocks_s8 <- list(
  "Physical parameters"   = intersect(c("T","alt"), names(df)),
  "Water chemistry"       = intersect(c("ph","EC"), names(df)),
  "Metals"                = intersect(c("Metals_PC1","Metals_PC2"), names(df)),
  "Microbial composition" = intersect(c("PCoA1","PCoA2"), names(df)),
  "Microbial diversity"   = intersect(c("Shannon","Richness"), names(df)),
  "ARGs"                  = intersect(c("ARG_total","ARG_rich"), names(df)),
  "MRGs"                  = intersect(c("MRG_total","MRG_rich"), names(df))
)

source_of <- function(v){
  if (v %in% c("T","alt","ph","EC") || v %in% c("Metals_PC1","Metals_PC2") || v %in% metals_used)
    return("metadata3.txt")
  if (v %in% c("PCoA1","PCoA2"))
    return("PCoA_Families_groupMerged.tsv (from Families_groupMerged.tsv)")
  if (v %in% c("Shannon","Richness"))
    return("Families_groupMerged.tsv (vegan)")
  if (str_detect(v, "^ARG_"))
    return("ARG per16S matrix (replicates) → median per Group")
  if (str_detect(v, "^MRG_"))
    return("MRG_metalOnly_per16S_groupMerged.tsv")
  return("")
}

preproc_of <- function(v){
  if (v == "EC") return("Median per Group; log1p; z-score")
  if (v %in% c("T","alt","ph")) return("Median per Group; z-score")
  if (v %in% c("Metals_PC1","Metals_PC2")) {
    if (length(metals_used) > 0) {
      return(paste0("Median per Group; log1p metals; PCA on {", paste(metals_used, collapse=", "), "}; z-score"))
    } else {
      return("Median per Group; log1p metals; PCA; z-score")
    }
  }
  if (v %in% c("PCoA1","PCoA2")) return("Bray–Curtis PCoA (families); z-score")
  if (v %in% c("Shannon","Richness")) return("Computed on families; z-score")
  if (v %in% c("ARG_total","MRG_total")) return("Group median; log1p; z-score")
  if (v %in% c("ARG_rich","MRG_rich")) return("Group median; z-score")
  return("z-score")
}

notes_of <- function(v){
  if (v=="T") return("Temperature (°C)")
  if (v=="alt") return("Altitude (m)")
  if (v=="ph") return("pH")
  if (v=="EC") return("Electrical conductivity (µS/cm)")
  if (v %in% c("Metals_PC1","Metals_PC2")) return("Trace-metal gradient (PCA score)")
  if (v %in% c("PCoA1","PCoA2")) return("Taxonomic composition (PCoA axis)")
  if (v=="Shannon") return("Shannon diversity (families)")
  if (v=="Richness") return("Observed family richness")
  if (v=="ARG_total") return("Total ARG abundance (per16S)")
  if (v=="ARG_rich") return("ARG richness (features >0)")
  if (v=="MRG_total") return("Total MRG abundance (per16S)")
  if (v=="MRG_rich") return("MRG richness (features >0)")
  return("")
}

table_s8 <- enframe(blocks_s8, name="Block", value="Indicator") %>%
  unnest(Indicator) %>%
  mutate(
    Source_file = sapply(Indicator, source_of),
    Preprocessing = sapply(Indicator, preproc_of),
    Notes_units = sapply(Indicator, notes_of)
  )

write_tsv(table_s8, "tables/TableS8_PLS_PM_input_variables.tsv")
write.csv(table_s8, "tables/TableS8_PLS_PM_input_variables.csv", row.names = FALSE)

# ==========
# TABLE S9: PCA loadings dei metalli + varianza spiegata
# ==========
if (length(metals_used) >= 2) {
  M <- df %>%
    select(all_of(metals_used)) %>%
    mutate(across(everything(), to_num)) %>%
    # se nel file non sono già log1p, applicalo qui (sicuro)
    mutate(across(everything(), ~ log1p(pmax(.x, 0)))) %>%
    as.data.frame()
  
  pca_met <- prcomp(M, center = TRUE, scale. = TRUE)
  
  loadings <- as.data.frame(pca_met$rotation[,1:2]) %>%
    rownames_to_column("Metal") %>%
    rename(PC1_loading = PC1, PC2_loading = PC2)
  
  var_exp <- tibble(
    Axis = c("PC1","PC2"),
    Variance_explained = (pca_met$sdev[1:2]^2) / sum(pca_met$sdev^2),
    Variance_explained_pct = 100 * Variance_explained
  )
  
  write_tsv(loadings, "tables/TableS9_MetalsPCA_loadings_PC1_PC2.tsv")
  write_tsv(var_exp,  "tables/TableS9_MetalsPCA_variance_PC1_PC2.tsv")
} else {
  warning("Non trovo abbastanza metalli nel file per fare Table S9 (serve >=2 metalli). Controlla metals_declared e i nomi colonne.")
}


