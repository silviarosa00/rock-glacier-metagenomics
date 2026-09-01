### 1.- Working Directory

setwd("C:/Users/SiRosa/OneDrive - Scientific Network South Tyrol/Projects/Subsurface/unite")
dir.create("SpiecEasi", showWarnings = FALSE, recursive = TRUE)
dir.create("SpiecEasi/plots", showWarnings = FALSE, recursive = TRUE)
dir.create("SpiecEasi/tables", showWarnings = FALSE, recursive = TRUE)
################################################################################
### 2.- required libraries.
# Required packages

install.packages("devtools")
install.packages("BiocManager")

# Since two of NetCoMi's dependencies are only available on GitHub, 
# it is recommended to install them first:
devtools::install_github("zdk123/SpiecEasi")
devtools::install_github("GraceYoon/SPRING")

# Install NetCoMi

devtools::install_github("stefpeschel/NetCoMi", 
                         repos = c("https://cloud.r-project.org/",
                                   BiocManager::repositories()))
####################

library(readr)
library(NetCoMi)
library(igraph)
library(tidygraph)
library(ggraph)
library(dplyr)

####################################################
#### 3.- DATA
save.image("data_myProject1.RData")

############################################################
######### MAKING MATRIX #####################
###########################################################
# ============================================================
# CREATE ONE COMBINED ARG + MRG MATRIX FOR SPIEC-EASI / NetCoMi
# - reads original matrices
# - keeps ALL genes
# - aligns common samples
# - cleans ARG names (remove numeric IDs, keep gene after "|")
# - cleans MRG names (keep gene after "|" when present)
# - adds ARG__/MRG__ prefixes
# - makes names unique
# - creates X_spiec = samples x features
# - saves optional mapping tables and matrix file
# ============================================================

# -----------------------
# 1) PACKAGES
# -----------------------
pkgs <- c("readr", "dplyr", "tibble")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install)) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}
invisible(lapply(pkgs, library, character.only = TRUE))

# -----------------------
# 2) WORKING DIRECTORY / INPUTS
# -----------------------
setwd("C:/Users/SiRosa/OneDrive - Scientific Network South Tyrol/Projects/Subsurface/unite")

arg_matrix_file <- "ARG_per16S_T1T2_renamed.tsv"   # genes x samples
mrg_matrix_file <- "total_metal_mrg_per16s.tsv"    # genes x samples

outdir <- "spiec_input_all_genes"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------
# 3) HELPERS
# -----------------------
clean_arg_name <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  
  # keep only part after the last "|", if present
  x <- ifelse(grepl("\\|", x), sub("^.*\\|", "", x), x)
  
  # remove leading numeric IDs if still present
  x <- sub("^[0-9]+[-_ :]*", "", x)
  
  x <- trimws(x)
  x
}

clean_mrg_name <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  
  # if format is BACxxxx|gene, keep gene
  x <- ifelse(grepl("\\|", x), sub("^.*\\|", "", x), x)
  
  # remove leading numeric IDs if present
  x <- sub("^[0-9]+[-_ :]*", "", x)
  
  x <- trimws(x)
  x
}

# -----------------------
# 4) READ ORIGINAL MATRICES
# -----------------------
arg_raw <- readr::read_tsv(arg_matrix_file, show_col_types = FALSE)
mrg_raw <- readr::read_tsv(mrg_matrix_file, show_col_types = FALSE)

arg_id_col <- names(arg_raw)[1]
mrg_id_col <- names(mrg_raw)[1]

ARG <- arg_raw %>%
  tibble::column_to_rownames(arg_id_col) %>%
  as.data.frame(check.names = FALSE)

MRG <- mrg_raw %>%
  tibble::column_to_rownames(mrg_id_col) %>%
  as.data.frame(check.names = FALSE)

# force numeric
ARG[] <- lapply(ARG, as.numeric)
MRG[] <- lapply(MRG, as.numeric)

ARG <- as.matrix(ARG)
MRG <- as.matrix(MRG)

# -----------------------
# 5) KEEP COMMON SAMPLES
# -----------------------
common_samp <- intersect(colnames(ARG), colnames(MRG))
if (length(common_samp) < 2) {
  stop("Too few common samples between ARG and MRG matrices.")
}

ARG <- ARG[, common_samp, drop = FALSE]
MRG <- MRG[, common_samp, drop = FALSE]

# optional: remove samples with zero total signal across ARG + MRG
keep_samp <- (colSums(ARG, na.rm = TRUE) + colSums(MRG, na.rm = TRUE)) > 0
ARG <- ARG[, keep_samp, drop = FALSE]
MRG <- MRG[, keep_samp, drop = FALSE]

# -----------------------
# 6) CLEAN FEATURE NAMES
# -----------------------
old_arg_names <- rownames(ARG)
old_mrg_names <- rownames(MRG)

new_arg_names <- make.unique(paste0("ARG__", clean_arg_name(old_arg_names)))
new_mrg_names <- make.unique(paste0("MRG__", clean_mrg_name(old_mrg_names)))

rownames(ARG) <- new_arg_names
rownames(MRG) <- new_mrg_names

# save name mappings
arg_map <- data.frame(
  original_name = old_arg_names,
  cleaned_name  = new_arg_names,
  stringsAsFactors = FALSE
)

mrg_map <- data.frame(
  original_name = old_mrg_names,
  cleaned_name  = new_mrg_names,
  stringsAsFactors = FALSE
)

readr::write_tsv(arg_map, file.path(outdir, "ARG_name_mapping.tsv"))
readr::write_tsv(mrg_map, file.path(outdir, "MRG_name_mapping.tsv"))

# -----------------------
# 7) COMBINE ARG + MRG
# -----------------------
# X = features x samples
X <- rbind(ARG, MRG)

# X_spiec = samples x features  --> ready for SPIEC-EASI / NetCoMi
X_spiec <- t(X)
X_spiec <- as.matrix(X_spiec)
storage.mode(X_spiec) <- "numeric"

# -----------------------
# 8) OPTIONAL SAVE MATRIX
# -----------------------
# save as TSV with samples in rows and features in columns
X_spiec_df <- as.data.frame(X_spiec, check.names = FALSE)
X_spiec_df <- tibble::rownames_to_column(X_spiec_df, var = "Sample")

readr::write_tsv(
  X_spiec_df,
  file.path(outdir, "ARG_MRG_all_genes_samples_x_features.tsv")
)

# also save R object
saveRDS(X_spiec, file.path(outdir, "X_spiec_all_genes.rds"))

# -----------------------
# 9) QUICK CHECKS
# -----------------------
cat("Saved files in:", outdir, "\n\n")
cat("N samples:", nrow(X_spiec), "\n")
cat("N total features:", ncol(X_spiec), "\n")
cat("N ARG features:", sum(grepl("^ARG__", colnames(X_spiec))), "\n")
cat("N MRG features:", sum(grepl("^MRG__", colnames(X_spiec))), "\n\n")

cat("First sample names:\n")
print(head(rownames(X_spiec)))

cat("\nFirst feature names:\n")
print(head(colnames(X_spiec)))

### Idealy, read pre-processed asv data as specific ASVs dataframes as a RDS file:                                           
#NetInput  <- readRDS("../data.rds" )  # 

### or read specific as dataframes:  
#NetInput <- read_delim("NetInput.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)


####################################################
### 4.- Filter and prepare the InPut for NetCoMi
# Data has to be a numeric Matrix
# as Matrix
pre_NetInput <- as.matrix(X_spiec)

#  as "numeric"
pre_NetInput <- as.matrix(pre_NetInput)
storage.mode(pre_NetInput) <- "numeric"
# -----------------------
# 1) basic checks
# -----------------------
cat("Class:", class(pre_NetInput), "\n")
cat("Mode:", mode(pre_NetInput), "\n")
cat("Dimensions (samples x features):", dim(pre_NetInput)[1], "x", dim(pre_NetInput)[2], "\n\n")

### Remove zeros
###    That is, only summ ASVs that are in at least on 60% of the sampling dates

numNA <- colSums(is.na(pre_NetInput))   # number of NA in each column
summary(numNA)
sum(numNA)                              # total number of NA in the matrix
sum(numNA > 0)                          # how many columns contain at least one NA
##no columns

numSamp <- rep(nrow(pre_NetInput), ncol(pre_NetInput))   # total sample size for each column
summary(numSamp)

number_ofZeros <- colSums(pre_NetInput == 0, na.rm = TRUE)   # number of zeros in each column
summary(number_ofZeros)
head(number_ofZeros)

# percentage of zeros in each column
propZeros <- number_ofZeros / nrow(pre_NetInput)
summary(propZeros)

# how many columns are very sparse
sum(propZeros > 0.5)    # more than 50% zeros
sum(propZeros > 0.8)    # more than 80% zeros
sum(propZeros > 0.9)    # more than 90% zeros

#filter per prevalence 5%
NetInput <- pre_NetInput[, colSums(pre_NetInput > 0, na.rm = TRUE) >= ceiling(0.05 * nrow(pre_NetInput))]

# check how many features remain
ncol(pre_NetInput) #969
ncol(NetInput) #605

### Check
# average number of reads per taxa
summary(colMeans(NetInput))

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.001282  0.041774  0.515507  4.832437  2.857696 97.760911


write.table(NetInput, "NetInput.txt" , quote = F , row.names = TRUE , col.names = TRUE , na = " ", sep = "\t" , eol = "\n") 


####################################################################  
### 5.- run NetCoMi
## by SpiecEasi -> all taxa (but set prevalence higher than 50%) 
# Safiya here I keep your values to 8 out of 12, check this
netConstructed <- netConstruct(NetInput,
                               filtTax = "none",
                               #filtTaxPar = ,  ## originale list(numbSamp = 8)
                               measure = "spieceasi",                # this should be at least 40 for publication. Use 10 in exploratory fase
                               measurePar = list(
                                 method = "mb",
                                 pulsar.params = list(rep.num = 50),
                                 nlambda = 20,
                                 symBetaMode = "ave"
                               ),  
                               #dissFunc = "signed",     # Transformation associations into dissimilarities
                               weighted = TRUE,
                               normMethod = "none",
                               zeroMethod = "none",
                               sparsMethod = "none",    # Sparsification: "none", or "t-test"
                               #alpha = 0.05,              # p-value for the Student's t-test or bootstrap
                               #adjust = "adaptBH",
                               #trueNullMethod = "convest",  # "convest"(default), "lfdr", "mean", and "hist"
                               #cores = 4L,
                               #softThreshCut = 0.8,
                               verbose = 3, seed = 88)

saveRDS(netConstructed, file= "./SpiecEasi/netConstructed_SpiecEasi.rds" ) 


########

netAnalyzd <- netAnalyze(netConstructed,
                         clustMethod = "cluster_fast_greedy",
                         hubPar = "degree", 
                         centrLCC = TRUE, # if TRUE, then compute centralities only for the largest connected component
                         hubQuant = 0.90,  # a node is identified as hub if for degree", the node’s centrality value is above the 90% quantile of the fitted log-normal distribution.
                         lnormFit = FALSE, #normal distribution is fitted to the centrality values to identify nodes with “highest” centrality values.
                         weightDeg = FALSE, normDeg = FALSE,
                         graphlet= FALSE, )     

# Weighted degree used (unweighted degree not meaningful for a fully connected network).

summary(netAnalyzd, numbNodes = 4L)
saveRDS(netAnalyzd, file= "./SpiecEasi/netAnalyzd_SpiecEasi.rds" )

save.image("data_myProject.RData")

###############################################
# I OBTAINED FEW ARG SO I TRY THIS 
####################################################
### Extract ARG-MRG subnetwork from NetCoMi object
####################################################

# 1) take matrices from netConstructed
A <- netConstructed$adjaMat1   # adjacency matrix
W <- netConstructed$assoMat1   # association matrix (signed weights)

# 2) identify ARG and MRG nodes
arg_nodes <- grep("^ARG__", colnames(A), value = TRUE)
mrg_nodes <- grep("^MRG__", colnames(A), value = TRUE)

length(arg_nodes) #170
length(mrg_nodes) #435

# 3) keep only ARG x MRG block
A_cross <- A[arg_nodes, mrg_nodes, drop = FALSE]
W_cross <- W[arg_nodes, mrg_nodes, drop = FALSE]

# 4) extract only non-zero edges
ix <- which(A_cross != 0, arr.ind = TRUE)

edges_ARG_MRG <- data.frame(
  ARG = rownames(A_cross)[ix[,1]],
  MRG = colnames(A_cross)[ix[,2]],
  weight = W_cross[ix],
  sign = ifelse(W_cross[ix] > 0, "pos", "neg"),
  stringsAsFactors = FALSE
)

# 5) order edges by absolute weight
edges_ARG_MRG <- edges_ARG_MRG[order(-abs(edges_ARG_MRG$weight)), ]

# save edge list
write.table(edges_ARG_MRG,
            file = "edges_ARG_MRG_from_NetCoMi.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 6) simple hub tables = number of cross-domain connections
arg_degree <- rowSums(A_cross != 0)
mrg_degree <- colSums(A_cross != 0)

arg_hubs <- data.frame(
  ARG = names(arg_degree),
  degree = as.numeric(arg_degree),
  stringsAsFactors = FALSE
)
arg_hubs <- arg_hubs[order(-arg_hubs$degree), ]

mrg_hubs <- data.frame(
  MRG = names(mrg_degree),
  degree = as.numeric(mrg_degree),
  stringsAsFactors = FALSE
)
mrg_hubs <- mrg_hubs[order(-mrg_hubs$degree), ]

write.table(arg_hubs,
            file = "ARG_hubs_from_NetCoMi_cross.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(mrg_hubs,
            file = "MRG_hubs_from_NetCoMi_cross.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 7) quick checks
cat("Number of ARG nodes in full network:", length(arg_nodes), "\n")   #170
cat("Number of MRG nodes in full network:", length(mrg_nodes), "\n") #435
cat("Number of ARG-MRG edges:", nrow(edges_ARG_MRG), "\n") #740
cat("Positive ARG-MRG edges:", sum(edges_ARG_MRG$sign == "pos"), "\n") #478
cat("Negative ARG-MRG edges:", sum(edges_ARG_MRG$sign == "neg"), "\n\n") #262

cat("Top ARG hubs:\n") 
print(head(arg_hubs, 20)) #

cat("\nTop MRG hubs:\n")
print(head(mrg_hubs, 20))

cat("\nTop ARG-MRG edges:\n")
print(head(edges_ARG_MRG, 20))

##                           ARG               MRG    weight sign
#472                  ARG__YajC   MRG__NCCX_ALCXX 0.6424745  pos
#237                  ARG__vanB MRG__B2I419_MYCMM 0.5376524  pos
#188                ARG__vanY.2   MRG__YBBM_ECOLI 0.5154565  pos
#633                 ARG__FEZ-1 MRG__H6WCN2_9FLAO 0.4572244  pos
#225                ARG__vanH.1   MRG__MERC_THIFE 0.4443969  pos
#76                 ARG__ADC-43   MRG__COPZ_ENTHA 0.4146274  pos
#138                  ARG__mecD MRG__D3UHP9_HELM1 0.4093344  pos
#429 ARG__AAC(3)-Ib/AAC(6')-Ib3    MRG__AIS_SALTY 0.4062509  pos
#598        ARG__Staphylococcus   MRG__MERB_SERMA 0.4007495  pos
#595                   ARG__bmr   MRG__MERB_STAAU 0.3754170  pos
#330                  ARG__H-NS MRG__Q9ZHD1_SALTM 0.3605039  pos
#71                   ARG__basS   MRG__COPK_RALME 0.3551231  pos
#189                ARG__vanY.4   MRG__YBBM_ECOLI 0.3533605  pos
#628                  ARG__vanW MRG__Q9WWL1_BACSR 0.3390614  pos
#399         ARG__Mycobacterium  MRG__Y1248_HAEIN 0.3353441  pos
#468                  ARG__MexT   MRG__NCCX_ALCXX 0.3266489  pos
#155                  ARG__vanB    MRG__DPS_STRSU 0.3206932  pos
#69                 ARG__CMY-26   MRG__COPK_RALME 0.3195364  pos
#564                  ARG__vanB MRG__F4ZBX9_XANCI 0.3160119  pos
#410                  ARG__rsmA   MRG__ZNUA_ECOLI 0.3040517  pos

#### ho taxa in arg come gestirli tolgo prima o dopo dell'analisi? una figura pulita → puoi toglierli dopo
# una rete finale da interpretare biologicamente → meglio toglierli prima e rifare netConstruct()


#--------------------------------------------------------------------------------
# Delete taxa before the analisis 
#-------------------------------------------------------------------------------
####################################################
### 4.- Filter and prepare input for NetCoMi
####################################################

pre_NetInput <- as.matrix(X_spiec)
storage.mode(pre_NetInput) <- "numeric"

cat("Class:", class(pre_NetInput), "\n")
cat("Mode:", mode(pre_NetInput), "\n")
cat("Dimensions (samples x features):", dim(pre_NetInput)[1], "x", dim(pre_NetInput)[2], "\n\n")

# ---- check NAs / zeros ----
numNA <- colSums(is.na(pre_NetInput))
summary(numNA)
sum(numNA)        # total NAs
sum(numNA > 0)    # number of columns with at least one NA

number_ofZeros <- colSums(pre_NetInput == 0, na.rm = TRUE)
summary(number_ofZeros)

propZeros <- number_ofZeros / nrow(pre_NetInput)
summary(propZeros)

sum(propZeros > 0.5)
sum(propZeros > 0.8)
sum(propZeros > 0.9)

# ---- light prevalence filter: keep features present in at least 5% of samples ----
NetInput <- pre_NetInput[, colSums(pre_NetInput > 0, na.rm = TRUE) >= ceiling(0.05 * nrow(pre_NetInput))]

cat("Original features:", ncol(pre_NetInput), "\n")
cat("Features after 5% prevalence filter:", ncol(NetInput), "\n")
cat("ARG after 5% filter:", sum(grepl("^ARG__", colnames(NetInput))), "\n")
cat("MRG after 5% filter:", sum(grepl("^MRG__", colnames(NetInput))), "\n\n")

# ---- remove ARG labels that are clearly taxa names BEFORE network inference ----
taxa_regex <- paste(
  c("Listeria",
    "Streptomyces", "Staphylococcus", "Pseudomonas", "Mycobacterium",
    "Acinetobacter", "Bacillus", "Burkholderia", "Enterococcus",
    "Escherichia", "Klebsiella", "Salmonella", "Listeria"),
  collapse = "|"
)
arg_cols <- grepl("^ARG__", colnames(NetInput))
arg_taxa_cols <- arg_cols & grepl(taxa_regex, colnames(NetInput), ignore.case = TRUE)

# also remove empty/broken labels
bad_empty <- colnames(NetInput) %in% c("ARG__", "MRG__", "ARG__NA", "MRG__NA", "")

# optional check: which features will be removed
removed_features <- colnames(NetInput)[arg_taxa_cols | bad_empty]
print(removed_features)

NetInput_clean <- NetInput[, !(arg_taxa_cols | bad_empty)]

cat("Features after removing taxon-like ARG labels:", ncol(NetInput_clean), "\n")
cat("ARG after cleaning:", sum(grepl("^ARG__", colnames(NetInput_clean))), "\n")
cat("MRG after cleaning:", sum(grepl("^MRG__", colnames(NetInput_clean))), "\n\n")

write.table(NetInput_clean,
            "NetInput_clean.txt",
            quote = FALSE, row.names = TRUE, col.names = TRUE,
            na = " ", sep = "\t", eol = "\n")


####################################################
### 5.- Run NetCoMi / SPIEC-EASI
####################################################

netConstructed <- netConstruct(
  NetInput_clean,
  filtTax = "none",
  filtSamp = "none",
  measure = "spieceasi",
  measurePar = list(
    method = "mb",
    pulsar.params = list(rep.num = 50),
    nlambda = 20,
    symBetaMode = "ave"
  ),
  weighted = TRUE,
  normMethod = "none",
  zeroMethod = "none",
  sparsMethod = "none",
  verbose = 3,
  seed = 88
)

saveRDS(netConstructed, file = "./SpiecEasi/netConstructed_SpiecEasi_clean.rds")


####################################################
### 6.- Analyze network
####################################################

netAnalyzd <- netAnalyze(
  netConstructed,
  clustMethod = "cluster_fast_greedy",
  hubPar = "degree",
  centrLCC = TRUE,
  hubQuant = 0.90,
  lnormFit = FALSE,
  weightDeg = FALSE,
  normDeg = FALSE,
  graphlet = FALSE
)

summary(netAnalyzd, numbNodes = 4L)
saveRDS(netAnalyzd, file = "./SpiecEasi/netAnalyzd_SpiecEasi_clean.rds")


####################################################
### 7.- Extract ARG-MRG subnetwork only
####################################################

A <- netConstructed$adjaMat1
W <- netConstructed$assoMat1

arg_nodes <- grep("^ARG__", colnames(A), value = TRUE)
mrg_nodes <- grep("^MRG__", colnames(A), value = TRUE)

cat("ARG nodes in full network:", length(arg_nodes), "\n")
cat("MRG nodes in full network:", length(mrg_nodes), "\n")

A_cross <- A[arg_nodes, mrg_nodes, drop = FALSE]
W_cross <- W[arg_nodes, mrg_nodes, drop = FALSE]

ix <- which(A_cross != 0, arr.ind = TRUE)

edges_ARG_MRG <- data.frame(
  ARG = rownames(A_cross)[ix[,1]],
  MRG = colnames(A_cross)[ix[,2]],
  weight = W_cross[ix],
  sign = ifelse(W_cross[ix] > 0, "pos", "neg"),
  stringsAsFactors = FALSE
)

edges_ARG_MRG <- edges_ARG_MRG[order(-abs(edges_ARG_MRG$weight)), ]

write.table(edges_ARG_MRG,
            file = "edges_ARG_MRG_from_NetCoMi_clean.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# simple cross-domain degree
arg_degree <- rowSums(A_cross != 0)
mrg_degree <- colSums(A_cross != 0)

arg_hubs <- data.frame(
  ARG = names(arg_degree),
  degree = as.numeric(arg_degree),
  stringsAsFactors = FALSE
)
arg_hubs <- arg_hubs[order(-arg_hubs$degree), ]

mrg_hubs <- data.frame(
  MRG = names(mrg_degree),
  degree = as.numeric(mrg_degree),
  stringsAsFactors = FALSE
)
mrg_hubs <- mrg_hubs[order(-mrg_hubs$degree), ]

write.table(arg_hubs,
            file = "ARG_hubs_from_NetCoMi_cross_clean.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(mrg_hubs,
            file = "MRG_hubs_from_NetCoMi_cross_clean.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Number of ARG-MRG edges:", nrow(edges_ARG_MRG), "\n")
cat("Positive ARG-MRG edges:", sum(edges_ARG_MRG$sign == "pos"), "\n")
cat("Negative ARG-MRG edges:", sum(edges_ARG_MRG$sign == "neg"), "\n\n")

cat("Top ARG hubs:\n")
print(head(arg_hubs, 20))

cat("\nTop MRG hubs:\n")
print(head(mrg_hubs, 20))

cat("\nTop ARG-MRG edges:\n")
print(head(edges_ARG_MRG, 20))


##################################################
# PLOT
########################################
####################################################
### CLEANER ARG–MRG NETWORK PLOT FROM NetCoMi OUTPUT
### Uses:
### - edges_ARG_MRG
### - arg_hubs
### - mrg_hubs
####################################################

library(dplyr)
library(igraph)
library(scales)



# -----------------------
# 1) optional: remove suspicious ARG labels again for plotting
# -----------------------
taxa_regex <- "Streptomyces|Staphylococcus|Pseudomonas|Mycobacterium|Acinetobacter|Bacillus|Burkholderia|Enterococcus|Escherichia|Klebsiella|Salmonella|Listeria"

edges_plot <- edges_ARG_MRG %>%
  filter(!is.na(ARG), !is.na(MRG), ARG != "", MRG != "") %>%
  filter(!grepl(taxa_regex, ARG, ignore.case = TRUE)) %>%
  filter(ARG != "ARG__", MRG != "MRG__")

# -----------------------
# 2) keep only strongest edges for readability
# change n = 40 / 60 / 80 depending on clutter
# -----------------------
top_edges <- edges_plot %>%
  slice_max(order_by = abs(weight), n = 60, with_ties = FALSE)

# -----------------------
# 3) build nodes table
# -----------------------
node_names <- unique(c(top_edges$ARG, top_edges$MRG))

nodes <- data.frame(
  name  = node_names,
  group = ifelse(grepl("^ARG__", node_names), "ARG", "MRG"),
  stringsAsFactors = FALSE
)

# -----------------------
# 4) build graph
# -----------------------
g <- graph_from_data_frame(
  d = top_edges[, c("ARG", "MRG", "weight", "sign")],
  directed = FALSE,
  vertices = nodes
)

# -----------------------
# 5) pastel palette
# -----------------------
# Nodes
col_ARG_node <- "#FFD166"   # ARG: giallo caldo
col_MRG_node <- "#CDB4DB"   # MRG: viola/lilla, molto distinto

# Edges
col_pos_edge <- "#4E79A7"   # positivo: verde
col_neg_edge <- "#D73027"   # negativo: rosso

# -----------------------
# 6) node styling
# -----------------------
V(g)$color <- ifelse(V(g)$group == "ARG", col_ARG_node, col_MRG_node)
V(g)$frame.color <- NA

deg <- degree(g)
V(g)$size <- 3 + 6 * (deg - min(deg)) / (max(deg) - min(deg) + 1e-12)

# -----------------------
# 7) edge styling
# -----------------------
E(g)$base_col <- ifelse(E(g)$sign == "pos", col_pos_edge, col_neg_edge)

w <- abs(E(g)$weight)
w01 <- (w - min(w)) / (max(w) - min(w) + 1e-12)
E(g)$width <- 1.2 + 3.5 * w01
E(g)$color <- adjustcolor(E(g)$base_col, alpha.f = 0.35)

# -----------------------
# 8) labels: top ARG hubs + top MRG hubs only
# -----------------------
# ---- labels from original mapping ----

# function: for MRG original names, keep only the part after the first "|" and before the second "|"
extract_mrg_gene <- function(x) {
  x <- as.character(x)
  out <- ifelse(grepl("\\|", x), sub("^[^|]*\\|([^|]+)\\|.*$", "\\1", x), x)
  out[is.na(out) | out == ""] <- x[is.na(out) | out == ""]
  out
}

# mapping: cleaned node name -> short MRG gene label
mrg_label_map <- setNames(
  extract_mrg_gene(mrg_map$original_name),
  mrg_map$cleaned_name
)

# ARG labels
lab_arg <- sub("^ARG__", "", V(g)$name)

# MRG labels from original mapping, fallback to cleaned name without prefix
lab_mrg <- unname(mrg_label_map[V(g)$name])
lab_mrg[is.na(lab_mrg)] <- sub("^MRG__", "", V(g)$name[is.na(lab_mrg)])

# final labels
V(g)$label <- ifelse(V(g)$group == "ARG", lab_arg, lab_mrg)

# plot settings
V(g)$label.cex <- 0.80
V(g)$label.color <- "black"
V(g)$label.dist <- 0.08
V(g)$label.degree <- 0

# optional manual nudges for specific labels if they exist
for (nm in c("vanY.2", "PA0320", "ACC(3)-lb/ACC(6')-lb", "")) {
  i <- which(V(g)$name == nm)
  if (length(i) == 1) {
    V(g)$label.dist[i] <- 0.30
  }
}
i <- which(V(g)$label == "AAC(3)-Ib/AAC(6')-Ib3")
V(g)$label.dist[i] <- 0.55
V(g)$label.degree[i] <- V(g)$label.degree[i] + 0.8
# -----------------------
# 9) layout: FR with positive weights only
# -----------------------
set.seed(1)
lay_fr <- layout_with_fr(g, weights = abs(E(g)$weight) + 1e-6)
lay_fr <- layout_with_fr(g, weights = abs(E(g)$weight) + 1e-6)

# stretch the layout a bit
lay_fr[,1] <- lay_fr[,1] * 1.25
lay_fr[,2] <- lay_fr[,2] * 1.15
# -----------------------
# 10) save PDF
# -----------------------
pdf(file.path(outdir, "FIG_network_ARG_MRG_NetCoMi_FR_clean.pdf"),
    width = 11, height = 7)

par(mar = c(0.2, 0.2, 0.2, 0.2), xpd = NA)

plot(
  g,
  layout = lay_fr,
  vertex.label = V(g)$label,
  vertex.label.cex = V(g)$label.cex,
  vertex.label.color = V(g)$label.color,
  vertex.label.dist = V(g)$label.dist,
  vertex.label.degree = V(g)$label.degree,
  vertex.size = V(g)$size,
  vertex.color = V(g)$color,
  vertex.frame.color = V(g)$frame.color,
  edge.color = E(g)$color,
  edge.width = E(g)$width,
  main = "ARG–MRG network (NetCoMi / SPIEC-EASI)"
)

legend(
  "topleft", bty = "n",
  legend = c("ARG node", "MRG node", "Positive edge", "Negative edge"),
  pch = c(16, 16, NA, NA),
  col = c(col_ARG_node, col_MRG_node, col_pos_edge, col_neg_edge),
  lty = c(NA, NA, 1, 1),
  lwd = c(NA, NA, 2, 2)
)

dev.off()

# -----------------------
# 11) save PNG
# -----------------------
png(file.path(outdir, "FIG_network_ARG_MRG_NetCoMi_FR_clean.png"),
    width = 6000, height = 4300, res = 600)

par(mar = c(0.5, 0.5, 2.5, 0.5))

plot(
  g,
  layout = lay_fr,
  vertex.label = V(g)$label,
  vertex.label.cex = V(g)$label.cex,
  vertex.label.color = V(g)$label.color,
  vertex.label.dist = V(g)$label.dist,
  vertex.label.degree = V(g)$label.degree,
  vertex.size = V(g)$size,
  vertex.color = V(g)$color,
  vertex.frame.color = V(g)$frame.color,
  edge.color = E(g)$color,
  edge.width = E(g)$width,
  main = "ARG–MRG network (NetCoMi / SPIEC-EASI)" , cex= 0.9
)

legend(
  "topleft", bty = "n",
  legend = c("ARG node", "MRG node", "Positive edge", "Negative edge"),
  pch = c(16, 16, NA, NA),
  col = c(col_ARG_node, col_MRG_node, col_pos_edge, col_neg_edge),
  cex = 0.8, 
  lty = c(NA, NA, 1, 1),
  lwd = c(NA, NA, 2, 2)
)

dev.off()



deg_table <- data.frame(
  gene = V(g)$label,
  node = V(g)$name,
  type = V(g)$group,
  degree = degree(g)
) %>%
  arrange(desc(degree))

print(deg_table)
summary(deg_table$degree)

hist(deg_table$degree, breaks = 20,
     main = "Degree distribution",
     xlab = "Degree")


print(hub_frequency)

top_ARG_hubs <- deg_table %>%
  filter(type == "ARG") %>%
  arrange(desc(degree)) %>%
  head(10)

print(top_ARG_hubs)
top_MRG_hubs <- deg_table %>%
  filter(type == "MRG") %>%
  arrange(desc(degree)) %>%
  head(15)

print(top_MRG_hubs)

# =========================================================
# NETWORK SUMMARY TABLE - COMI ARG-MRG NETWORK
# =========================================================

library(igraph)
library(dplyr)
library(readr)
library(tibble)

node_number <- vcount(g)
edge_number <- ecount(g)

positive_edge_percentage <- mean(E(g)$sign == "pos") * 100
negative_edge_percentage <- mean(E(g)$sign == "neg") * 100

average_path_length <- mean_distance(
  g,
  directed = FALSE,
  weights = NA,
  unconnected = TRUE
)

deg <- degree(g, mode = "all")

mean_degree <- mean(deg, na.rm = TRUE)
median_degree <- median(deg, na.rm = TRUE)
max_degree <- max(deg, na.rm = TRUE)

btw <- betweenness(
  g,
  directed = FALSE,
  weights = NA,
  normalized = TRUE
)

mean_betweenness <- mean(btw, na.rm = TRUE)
median_betweenness <- median(btw, na.rm = TRUE)
max_betweenness <- max(btw, na.rm = TRUE)

comm <- cluster_fast_greedy(
  g,
  weights = abs(E(g)$weight)
)

modularity_value <- modularity(
  comm,
  weights = abs(E(g)$weight)
)

clustering_coefficient <- transitivity(
  g,
  type = "global",
  isolates = "zero"
)

average_clustering_coefficient <- transitivity(
  g,
  type = "average",
  isolates = "zero"
)

network_summary_COMI <- tibble(
  metric = c(
    "Node number",
    "Edge number",
    "Average path length",
    "Mean betweenness",
    "Median betweenness",
    "Max betweenness",
    "Mean degree",
    "Median degree",
    "Max degree",
    "Modularity",
    "Positive edge percentage",
    "Negative edge percentage",
    "Global clustering coefficient",
    "Average clustering coefficient"
  ),
  value = c(
    node_number,
    edge_number,
    average_path_length,
    mean_betweenness,
    median_betweenness,
    max_betweenness,
    mean_degree,
    median_degree,
    max_degree,
    modularity_value,
    positive_edge_percentage,
    negative_edge_percentage,
    clustering_coefficient,
    average_clustering_coefficient
  )
)

print(network_summary_COMI)

write_tsv(
  network_summary_COMI,
  file.path(outdir, "network_summary_COMI_ARG_MRG.tsv")
)




