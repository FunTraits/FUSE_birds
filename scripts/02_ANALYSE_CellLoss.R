#-------------------------------------------------------------------------------
# Biogeographical extinction-scenario analysis — cell-by-cell.
#
# For each DGGS cell (res 7) and each of the three functional
# dimensions (locomotion, diet, reproduction), computes under two
# extinction scenarios (stochastic IUCN 100 and deterministic IUCN AT):
#   * % change in species richness (Species change)
#   * % de changement de FRic    (FRic change)
#   * % de changement de FUn     (FUn change)
#
# Values are expressed as signed percentage changes,
# i.e. negative for losses (consistent with Fig. 2).
#
# Author  : A. Toussaint
# Prerequisites (objects in memory or reloaded via `readRDS`):
#   - tpd_spaces : named list of TPDs (locomotion / diet / reproduction)
#   - coords     : named list of data.frames (species, PC1, PC2, ...)
#   - species_df : data.frame with columns species, iucn (LC..CR)
#   - sitesdggs7 : named list cell_id -> species vector
#   - dggs_grid  : sf object of DGGS cells (with column 'seqnum')
#                  (can be omitted: used only if
#                   exporting results as spatial object)
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c("dplyr", "tibble", "tidyr", "purrr", "RANN",
                   "future.apply", "progressr")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))


# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(
  
  # Number of stochastic iterations for IUCN 100 PER CELL
  # (default 200 — more efficient than 1000; sufficient to estimate the
  #  mean and 95% CI of loss; see stability test)
  n_iter           = 200,
  
  # IUCN 100: 100-year extinction probabilities
  # (Mooers et al. 2008 / Davis et al. 2018)
  p_extinction_100yr = c(
    LC = 0.0009,
    NT = 0.0071,
    VU = 0.10,
    EN = 0.667,
    CR = 0.999
  ),
  
  # Categories considered threatened (IUCN AT scenario)
  threatened_cats  = c("VU", "EN", "CR"),
  
  # Axes PCoA pour le calcul de FUn
  pcoa_axes        = c("PC1", "PC2"),
  k_neighbours     = 5,
  
  # Minimum number of species in a cell to compute FUn
  # (sinon FUn = NA pour cette cellule)
  min_n_for_FUn    = 8,
  
  # Parallelisation
  parallel         = TRUE,
  n_workers        = max(1, parallel::detectCores() - 1),
  
  # Reproducibility
  seed             = 20251028,
  
  # Output
  out_dir          = "data/processed",
  out_file         = "extinction_scenarios_cellwise.rds"
)
dir.create(PARAMS$out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(PARAMS$seed)
options(future.globals.maxSize = 2 * 1024^3)


# ============================================================================
# 2. USER INPUTS ------------------------------------------------------
# ============================================================================
tpd_spaces <- readRDS("data/processed/tpd_lists.rds")
coords     <- readRDS("data/processed/pcoa_coords.rds")
species_df <- readRDS("data/processed/species_table.rds")
sitesdggs7 <- readRDS("data/raw/sitesdggs7.RDS")
library(dggridR)
dggs <- dgconstruct(res = 7)
grid <- dgearthgrid(dggs)
dggs_grid <- grid   # alias

stopifnot(
  exists("tpd_spaces"),
  all(c("locomotion", "diet", "reproduction") %in% names(tpd_spaces)),
  exists("coords"),
  all(c("locomotion", "diet", "reproduction") %in% names(coords)),
  exists("species_df"),
  all(c("species", "iucn") %in% names(species_df)),
  exists("sitesdggs7"),
  is.list(sitesdggs7),
  all(species_df$iucn %in% names(PARAMS$p_extinction_100yr))
)

# Clean assemblages: retain only species present in all three spaces
# le pool global (intersection avec species_df$species)
pool_species <- species_df$species
sitesdggs7_clean <- lapply(sitesdggs7, function(sp) intersect(sp, pool_species))
# Remove cells that are empty after filtering
n_per_cell <- vapply(sitesdggs7_clean, length, integer(1))
sitesdggs7_clean <- sitesdggs7_clean[n_per_cell > 0]
message(sprintf("Cellules conservées : %d / %d (avec >= 1 espèce du pool)",
                length(sitesdggs7_clean), length(sitesdggs7)))


# ============================================================================
# 3. Fonctions de calcul ------------------------------------------------------
# ============================================================================

#' Build compact representation of kernel supports.
#' Each species -> integer indices of non-zero TPD grid cells.
build_compact_supports <- function(tpd_obj) {
  kernels <- tpd_obj$TPDs
  if (is.null(kernels))
    stop("tpd_obj$TPDs est NULL : vérifier l'objet TPD::TPDs.")
  n_cells <- length(kernels[[1]])
  supports <- lapply(kernels, function(k) which(k > 0))
  list(supports = supports, n_cells = n_cells)
}

#' FRic of an assemblage from compact supports.
#' Returns the fraction of the global trait-space occupied by the assemblage.
compute_FRic_compact <- function(compact, subset_sp) {
  subset_sp <- intersect(subset_sp, names(compact$supports))
  if (length(subset_sp) == 0) return(0)
  occupied <- unique(unlist(compact$supports[subset_sp], use.names = FALSE))
  length(occupied) / compact$n_cells
}

#' FUn moyen sur un sous-ensemble dans un espace PCoA.
#' Renvoie NA si n_subset <= k.
compute_FUn_mean <- function(mat, sp_names, subset_sp, k) {
  idx <- which(sp_names %in% subset_sp)
  if (length(idx) <= k) return(NA_real_)
  m_sub <- mat[idx, , drop = FALSE]
  nn <- RANN::nn2(data = m_sub, query = m_sub, k = k + 1)
  mean(rowMeans(nn$nn.dists[, -1, drop = FALSE]))
}


# ============================================================================
# 4. Prepare dimensional objects -------- ------------------------------------
# ============================================================================
prep_dim <- function(dim_name) {
  
  tpd_obj <- tpd_spaces[[dim_name]]
  compact <- build_compact_supports(tpd_obj)
  rm(tpd_obj); gc(verbose = FALSE)
  
  cdf <- coords[[dim_name]]
  cdf <- cdf[match(species_df$species, cdf$species), , drop = FALSE]
  if (any(is.na(cdf$species)))
    stop("Espèces absentes dans coords[[", dim_name, "]].")
  mat <- as.matrix(cdf[, PARAMS$pcoa_axes, drop = FALSE])
  
  size_mb <- as.numeric(object.size(compact)) / 1024^2
  message(sprintf("  Compact supports: %.1f MiB", size_mb))
  
  list(
    compact   = compact,
    mat       = mat,
    sp_names  = species_df$species
  )
}


# ============================================================================
# 5. Coeur de l'analyse cellule ----------------------------------------------
# ============================================================================

#' For a given cell, compute all relative changes
#' under IUCN 100 (averaged over n_iter iterations) and IUCN AT (deterministic).
#'
#' Returns a named vector with 9 values:
#'   n_sp_baseline,
#'   sp_change_100, FRic_change_100, FUn_change_100,
#'   sp_change_AT,  FRic_change_AT,  FUn_change_AT,
#'   sp_change_100_sd, FRic_change_100_sd
#'  (les SD pour IUCN 100 — utiles pour cartographier l'incertitude)
process_cell <- function(cell_species, prep, p_ext_lookup,
                         threatened_lookup, n_iter, k_min) {
  
  cell_species <- unique(cell_species)
  n_baseline <- length(cell_species)
  
  # Baseline FRic et FUn pour la cellule
  FRic_base <- compute_FRic_compact(prep$compact, cell_species)
  FUn_base  <- compute_FUn_mean(prep$mat, prep$sp_names, cell_species,
                                PARAMS$k_neighbours)
  
  if (FRic_base == 0 || is.na(FUn_base) || n_baseline < k_min) {
    # Cellule trop pauvre : on retourne NA pour tous les changements
    return(c(
      n_sp_baseline       = n_baseline,
      sp_change_100       = NA_real_,
      FRic_change_100     = NA_real_,
      FUn_change_100      = NA_real_,
      sp_change_AT        = NA_real_,
      FRic_change_AT      = NA_real_,
      FUn_change_AT       = NA_real_,
      sp_change_100_sd    = NA_real_,
      FRic_change_100_sd  = NA_real_
    ))
  }
  
  # ----- IUCN AT (deterministic: 1 computation) -----
  threatened_in_cell <- threatened_lookup[cell_species]
  survivors_AT <- cell_species[!threatened_in_cell]
  n_extinct_AT <- sum(threatened_in_cell)
  
  FRic_AT <- compute_FRic_compact(prep$compact, survivors_AT)
  FUn_AT  <- compute_FUn_mean(prep$mat, prep$sp_names, survivors_AT,
                              PARAMS$k_neighbours)
  
  sp_change_AT   <- -100 * n_extinct_AT / n_baseline
  FRic_change_AT <- -100 * (1 - FRic_AT / FRic_base)
  FUn_change_AT  <- if (is.na(FUn_AT)) NA_real_
  else 100 * (FUn_AT - FUn_base) / FUn_base
  
  # ----- IUCN 100 (stochastic: n_iter iterations) -----
  p_ext_in_cell <- p_ext_lookup[cell_species]
  
  sp_loss_100   <- numeric(n_iter)
  FRic_100      <- numeric(n_iter)
  FUn_100       <- numeric(n_iter)
  
  for (i in seq_len(n_iter)) {
    extinct_i <- rbinom(n_baseline, 1, p_ext_in_cell) == 1
    survivors_i <- cell_species[!extinct_i]
    sp_loss_100[i] <- sum(extinct_i)
    FRic_100[i]    <- compute_FRic_compact(prep$compact, survivors_i)
    FUn_100[i]     <- compute_FUn_mean(prep$mat, prep$sp_names, survivors_i,
                                       PARAMS$k_neighbours)
  }
  
  # Average over iterations then convert to % change
  mean_sp_loss_100  <- mean(sp_loss_100)
  sd_sp_loss_100    <- sd(sp_loss_100)
  mean_FRic_100     <- mean(FRic_100)
  sd_FRic_100       <- sd(FRic_100)
  mean_FUn_100      <- mean(FUn_100, na.rm = TRUE)
  
  sp_change_100     <- -100 * mean_sp_loss_100 / n_baseline
  FRic_change_100   <- -100 * (1 - mean_FRic_100 / FRic_base)
  FUn_change_100    <- if (is.na(mean_FUn_100)) NA_real_
  else 100 * (mean_FUn_100 - FUn_base) / FUn_base
  
  # SD de la perte en pourcentage (utile pour barres d'erreur cartographiques)
  sp_change_100_sd   <- 100 * sd_sp_loss_100 / n_baseline
  FRic_change_100_sd <- 100 * (sd_FRic_100 / FRic_base)
  
  c(
    n_sp_baseline       = n_baseline,
    sp_change_100       = sp_change_100,
    FRic_change_100     = FRic_change_100,
    FUn_change_100      = FUn_change_100,
    sp_change_AT        = sp_change_AT,
    FRic_change_AT      = FRic_change_AT,
    FUn_change_AT       = FUn_change_AT,
    sp_change_100_sd    = sp_change_100_sd,
    FRic_change_100_sd  = FRic_change_100_sd
  )
}


# ============================================================================
# 6. Lookup tables (fast access inside the parallelised loop) -------- --------
# ============================================================================
p_ext_lookup     <- setNames(PARAMS$p_extinction_100yr[as.character(species_df$iucn)],
                             species_df$species)
threatened_lookup <- setNames(species_df$iucn %in% PARAMS$threatened_cats,
                              species_df$species)


# ============================================================================
# 7. Boucle principale : par dimension, par cellule --------------------------
# ============================================================================
if (PARAMS$parallel) {
  future::plan(future::multisession, workers = PARAMS$n_workers)
} else {
  future::plan(future::sequential)
}
progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")

cell_results <- list()

for (dim_name in c("locomotion", "diet", "reproduction")) {
  
  message("\n========== ", dim_name, " ==========")
  prep <- prep_dim(dim_name)
  
  cells <- names(sitesdggs7_clean)
  message(sprintf("  Traitement de %d cellules x %d itérations IUCN 100",
                  length(cells), PARAMS$n_iter))
  
  progressr::with_progress({
    p <- progressr::progressor(steps = length(cells))
    res_mat <- future.apply::future_sapply(
      cells,
      function(cid) {
        out <- process_cell(
          cell_species      = sitesdggs7_clean[[cid]],
          prep              = prep,
          p_ext_lookup      = p_ext_lookup,
          threatened_lookup = threatened_lookup,
          n_iter            = PARAMS$n_iter,
          k_min             = PARAMS$min_n_for_FUn
        )
        p()
        out
      },
      future.seed = TRUE
    )
  })
  
  # Reshape: cells in rows, metrics in columns
  res_df <- as.data.frame(t(res_mat))
  res_df <- tibble::tibble(cell = cells, !!!res_df)
  
  cell_results[[dim_name]] <- res_df
  
  # Diagnostics rapides
  msg <- sprintf(
    paste0("  IUCN 100 : sp = %+0.1f%% (sd %.1f) | FRic = %+0.1f%% (sd %.1f) | FUn = %+0.1f%%\n",
           "  IUCN AT  : sp = %+0.1f%%        | FRic = %+0.1f%%        | FUn = %+0.1f%%"),
    mean(res_df$sp_change_100,    na.rm = TRUE),
    mean(res_df$sp_change_100_sd, na.rm = TRUE),
    mean(res_df$FRic_change_100,   na.rm = TRUE),
    mean(res_df$FRic_change_100_sd, na.rm = TRUE),
    mean(res_df$FUn_change_100,    na.rm = TRUE),
    mean(res_df$sp_change_AT,      na.rm = TRUE),
    mean(res_df$FRic_change_AT,    na.rm = TRUE),
    mean(res_df$FUn_change_AT,     na.rm = TRUE)
  )
  message(msg)
}

future::plan(future::sequential)


# ============================================================================
# 8. Assemblage en table large unique ----------------------------------------
# ============================================================================
# Each cell has 3 dimensions x (n_sp + 6 metrics + 2 SD) = 27 columns.
# Prefix by dimension to disambiguate.

cell_long <- bind_rows(lapply(names(cell_results), function(d) {
  cell_results[[d]] %>%
    mutate(dimension = d)
}))

cell_wide <- cell_long %>%
  pivot_wider(
    id_cols = cell,
    names_from = dimension,
    values_from = c(n_sp_baseline,
                    sp_change_100, FRic_change_100, FUn_change_100,
                    sp_change_AT,  FRic_change_AT,  FUn_change_AT,
                    sp_change_100_sd, FRic_change_100_sd),
    names_glue = "{.value}__{dimension}"
  )


# ============================================================================
# 9. Sauvegarde ---------------------------------------------------------------
# ============================================================================
out_path <- file.path(PARAMS$out_dir, PARAMS$out_file)
saveRDS(list(cell_long = cell_long, cell_wide = cell_wide,
             params = PARAMS), out_path)
message("\nRésultats cellule sauvegardés dans : ", out_path)


# ============================================================================
# 10. Optional: enrich with DGGS geometry for mapping -------- -------------
# ============================================================================
# If dggs_grid is available, join results to cell geometries.
# Join expects dggs_grid to have a numeric column 'seqnum'.
# and that cell_wide$cell be comparable (also numeric).

if (exists("dggs_grid")) {
  cell_wide$seqnum <- as.numeric(cell_wide$cell)
  cell_sf <- dggs_grid %>%
    dplyr::left_join(cell_wide, by = "seqnum")
  saveRDS(cell_sf,
          file.path(PARAMS$out_dir, "extinction_scenarios_cellwise_sf.rds"))
  message("Version sf sauvegardée pour cartographie.")
}


# ============================================================================
# 11. Reproducibility -------- --------------------------------------------------------
# ============================================================================
sessionInfo()