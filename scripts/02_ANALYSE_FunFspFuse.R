#-------------------------------------------------------------------------------
# A.1 — Computation of FUn and FSp in each functional space
#       (strict alignment with Pimiento et al. 2020, Sci. Adv.)
#
# Author  : A. Toussaint
# Output  : species_df augmented with FUn_<dim>, FSp_<dim>, FUn_std_<dim>,
#           FSp_std_<dim> and FUSE_<dim> for all three dimensions, saved as RDS.
#
# Definitions (following Pimiento et al. 2020, Materials & Methods, eq. 1-3):
#   * FUn = mean Euclidean distance to the k nearest neighbours in
#           the retained PCoA space (k = 5 by default, following Pimiento)
#   * FSp = Euclidean distance to the global centroid of the space
#   * Standardisation [0,1] then multiplication by 4
#   * GE = ordinal IUCN score: LC=0, NT=1, VU=2, EN=3, CR=4
#   * FUSE = ln(1 + FUn_std * GE) + ln(1 + FSp_std * GE)
#
# Prerequisites:
#   - coords : named list of data.frames, one per dimension
#              ("locomotion", "diet", "reproduction"), with columns
#              species, PC1, PC2 (...). The script uses by default
#              the first two PCoA axes, adjustable via PARAMS$pcoa_axes.
#   - species_df : data.frame with columns species and iucn (5 IUCN levels)
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c("dplyr", "tibble", "purrr", "scales", "RANN")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))


# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(
  # PCoA axes used to compute distances
  # default: first 2 (consistent with Methods: 2D PCoA)
  pcoa_axes        = c("PC1", "PC2"),
  
  # nombre de voisins pour FUn (Pimiento 2020 : k = 5)
  k_neighbours     = 5,
  
  # multiplicative factor after [0,1] standardisation
  # (Pimiento 2020, Materials & Methods : "FUn and FSp standardised, multiplied by 4")
  fuse_scale       = 4,
  
  # IUCN -> ordinal score mapping (Pimiento 2020; same as Mooers et al. 2008)
  iucn_to_GE       = c(LC = 0, NT = 1, VU = 2, EN = 3, CR = 4),
  
  # dimensions to process (output column suffixes)
  dim_suffixes     = c(locomotion = "loco",
                       diet       = "diet",
                       reproduction = "repro"),
  
  # sortie
  out_dir          = "data/processed",
  out_file         = "species_metrics_FUn_FSp_FUSE.rds"
)
dir.create(PARAMS$out_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================================
# 2. USER INPUTS — adjust paths -------------------------------
# ============================================================================
coords     <- readRDS("data/processed/pcoa_coords.rds")
species_df <- readRDS("data/processed/species_table.rds")

stopifnot(
  exists("coords"),
  is.list(coords),
  all(names(PARAMS$dim_suffixes) %in% names(coords)),
  exists("species_df"),
  all(c("species", "iucn") %in% names(species_df))
)

# Check: all spaces share the same species list
species_lists <- lapply(coords, function(d) sort(unique(d$species)))
if (!all(sapply(species_lists[-1], identical, species_lists[[1]]))) {
  warning("Species lists differ across dimensions; ",
          "le calcul utilisera l'intersection.")
  common <- Reduce(intersect, species_lists)
  coords <- lapply(coords, function(d) d[d$species %in% common, , drop = FALSE])
  species_df <- species_df[species_df$species %in% common, , drop = FALSE]
}


# ============================================================================
# 3. Fonctions de calcul ------------------------------------------------------
# ============================================================================

#' Functional Uniqueness (Pimiento 2020) : distance euclidienne moyenne
#' aux k plus proches voisins dans un espace euclidien.
#'
#' kd-tree implementation (RANN::nn2) — scalable to n > 10^4 without
#' allocating the full distance matrix.
#'
#' @param mat  numeric matrix n x p (rows = species, columns = PCoA axes)
#' @param k    number of neighbours (default 5)
#' @return     numeric vector of length n
compute_FUn <- function(mat, k = PARAMS$k_neighbours) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  n <- nrow(mat)
  if (n <= k) stop("n (", n, ") doit être > k (", k, ").")
  # nn2 returns k+1 neighbours (1st column = the species itself)
  nn <- RANN::nn2(data = mat, query = mat, k = k + 1)
  # on retire la colonne 1 (auto-distance = 0) et on moyenne
  rowMeans(nn$nn.dists[, -1, drop = FALSE])
}

#' Functional Specialisation (Pimiento 2020) : distance euclidienne au
#' global centroid of the space.
#'
#' @param mat  numeric matrix n x p
#' @return     numeric vector of length n
compute_FSp <- function(mat) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  centroid <- colMeans(mat, na.rm = TRUE)
  sqrt(rowSums(sweep(mat, 2, centroid, FUN = "-")^2))
}

#' Standardisation [0,1] puis multiplication par un facteur (×4 chez Pimiento).
rescale_x4 <- function(x, factor = PARAMS$fuse_scale) {
  factor * scales::rescale(x, to = c(0, 1))
}

#' FUSE selon Pimiento et al. (2020), eq. 1-3 :
#'   FUSE = ln(1 + FUn_std * GE) + ln(1 + FSp_std * GE)
compute_FUSE <- function(FUn_std, FSp_std, GE) {
  log(1 + FUn_std * GE) + log(1 + FSp_std * GE)
}


# ============================================================================
# 4. Mapping IUCN -> GE -------------------------------------------------------
# ============================================================================
# Coverage check: DD/NE species must be imputed
# BEFORE this step (see pipeline: imputation by sampling from the
# distribution empirique des statuts par ordre, suivant Pimiento 2020).
if (any(!species_df$iucn %in% names(PARAMS$iucn_to_GE))) {
  stop("Statuts IUCN non reconnus : ",
       paste(unique(species_df$iucn[!species_df$iucn %in%
                                      names(PARAMS$iucn_to_GE)]),
             collapse = ", "),
       ". Imputez DD/NE avant d'appeler ce script.")
}
species_df$GE <- PARAMS$iucn_to_GE[as.character(species_df$iucn)]


# ============================================================================
# 5. Compute metrics per dimension -------- -----------------------------------------------------
# ============================================================================

results <- list()

for (dim_name in names(PARAMS$dim_suffixes)) {
  
  message("--- ", dim_name, " ---")
  
  suf <- PARAMS$dim_suffixes[[dim_name]]
  cdf <- coords[[dim_name]]
  
  # axis check
  if (!all(PARAMS$pcoa_axes %in% names(cdf)))
    stop("Axes ", paste(PARAMS$pcoa_axes, collapse = ", "),
         " absents pour la dimension ", dim_name, ".")
  
  # align to species_df
  cdf <- cdf[match(species_df$species, cdf$species), , drop = FALSE]
  if (any(is.na(cdf$species)))
    stop("Espèces absentes dans coords[[", dim_name, "]].")
  
  mat <- as.matrix(cdf[, PARAMS$pcoa_axes, drop = FALSE])
  
  # 5a. raw metrics
  FUn_raw <- compute_FUn(mat, k = PARAMS$k_neighbours)
  FSp_raw <- compute_FSp(mat)
  
  # 5b. standardisation Pimiento ([0,1] x 4)
  FUn_std <- rescale_x4(FUn_raw)
  FSp_std <- rescale_x4(FSp_raw)
  
  # 5c. FUSE
  FUSE_val <- compute_FUSE(FUn_std, FSp_std, species_df$GE)
  
  # 5d. FUS sans poids d'extinction (Griffin et al. 2020) :
  #     FUS = (FUn_std + FSp_std) / 2 — utile pour les analyses descriptives
  FUS_val <- (FUn_std + FSp_std) / 2
  
  # 5e. stockage
  results[[paste0("FUn_",     suf)]] <- FUn_raw
  results[[paste0("FSp_",     suf)]] <- FSp_raw
  results[[paste0("FUn_std_", suf)]] <- FUn_std
  results[[paste0("FSp_std_", suf)]] <- FSp_std
  results[[paste0("FUS_",     suf)]] <- FUS_val
  results[[paste0("FUSE_",    suf)]] <- FUSE_val
  
  # diagnostic
  message(sprintf(
    "  n = %d | FUn (raw): mean=%.3f sd=%.3f | FSp (raw): mean=%.3f sd=%.3f",
    nrow(mat),
    mean(FUn_raw), sd(FUn_raw),
    mean(FSp_raw), sd(FSp_raw)
  ))
  message(sprintf(
    "  Top-3 FUSE: %s",
    paste(species_df$species[order(-FUSE_val)][1:3], collapse = ", ")
  ))
}


# ============================================================================
# 6. Assemblage et sauvegarde -------------------------------------------------
# ============================================================================
metrics_df <- as_tibble(c(list(species = species_df$species,
                               iucn   = species_df$iucn,
                               GE     = species_df$GE),
                          results))

# sanity: pas de NA inattendus
n_na <- sapply(metrics_df, function(x) sum(is.na(x)))
if (any(n_na > 0)) {
  message("Colonnes contenant des NA :")
  print(n_na[n_na > 0])
}

out_path <- file.path(PARAMS$out_dir, PARAMS$out_file)
saveRDS(metrics_df, out_path)
message("\nMétriques sauvegardées dans : ", out_path)


# ============================================================================
# 7. Cross-checks -------- ---------------------------------------------------
# ============================================================================

# 7a. Spearman correlations between FUn per pair of dimensions
#     (direct equivalent of the ρ values reported in the Results)
fun_cor <- cor(
  metrics_df[, paste0("FUn_", PARAMS$dim_suffixes)],
  method = "spearman", use = "pairwise.complete.obs"
)
message("\nCorrélations Spearman entre FUn par dimension :")
print(round(fun_cor, 3))

# 7b. FUSE correlations per pair of dimensions
fuse_cor <- cor(
  metrics_df[, paste0("FUSE_", PARAMS$dim_suffixes)],
  method = "spearman", use = "pairwise.complete.obs"
)
message("\nCorrélations Spearman entre FUSE par dimension :")
print(round(fuse_cor, 3))

# 7c. Distribution of top 5% FUSE per dimension
top_5pct <- lapply(PARAMS$dim_suffixes, function(suf) {
  v <- metrics_df[[paste0("FUSE_", suf)]]
  thr <- quantile(v, 0.95, na.rm = TRUE)
  metrics_df$species[v >= thr]
})
message("
Number of species in top 5% FUSE per dimension:")
print(sapply(top_5pct, length))


# ============================================================================
# 8. Reproducibility -------- ---------------------------------------------------------
# ============================================================================
sessionInfo()
