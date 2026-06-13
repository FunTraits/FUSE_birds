#-------------------------------------------------------------------------------
# B — Extinction scenarios (core analysis) — SES version
#
#   For each functional dimension (locomotion, diet, reproduction):
#   - simulates IUCN 100 and IUCN AT scenarios (Pimiento et al. 2020)
#   - computes FRic loss and FUn change per iteration
#   - runs null randomisations preserving the number of extinctions
#   - quantifies deviation from null via Standardized Effect Size (SES)
#     with permutation p-value (Davison & Hinkley 1997; +1/+1 rule)
#
#   Changes from previous version:
#   - Replaces wilcox_paired() and wilcox_one_sample() with compute_SES().
#   - $ses_stats replaces $pvals (alias $pvals kept for backward compatibility).
#   - Global summary table now includes SES, null_mean, null_sd.
#   - Panel annotations in Fig. 2 display "SES = X.XX, p = Y.YY"
#
# Author  : A. Toussaint
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
  
  n_iter           = 1000,
  
  p_extinction_100yr = c(
    LC = 0.0009,
    NT = 0.0071,
    VU = 0.10,
    EN = 0.667,
    CR = 0.999
  ),
  
  threatened_cats  = c("VU", "EN", "CR"),
  
  pcoa_axes        = c("PC1", "PC2"),
  k_neighbours     = 5,
  
  parallel         = TRUE,
  n_workers        = max(1, parallel::detectCores() - 1),
  
  seed             = 20251028,
  
  # Direction des tests SES
  #   FRic loss : alternative "greater" (perte empirique > null attendu)
  #   FUn change : alternative "greater" (gain empirique > null attendu)
  #   Use "two.sided" for a two-tailed test.
  ses_alternative_FRic = "greater",
  ses_alternative_FUn  = "greater",
  
  out_dir          = "data/processed",
  out_file         = "extinction_scenarios_results.rds"
)
dir.create(PARAMS$out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(PARAMS$seed)

options(future.globals.maxSize = 2 * 1024^3)


# ============================================================================
# 2. USER INPUTS — adjust paths -------------------------------
# ============================================================================
tpd_spaces <- readRDS("data/processed/tpd_lists.rds")
coords     <- readRDS("data/processed/pcoa_coords.rds")
species_df <- readRDS("data/processed/species_table.rds")

stopifnot(
  exists("tpd_spaces"),
  all(c("locomotion","diet","reproduction") %in% names(tpd_spaces)),
  exists("species_df"),
  all(c("species","iucn") %in% names(species_df)),
  exists("coords"),
  all(c("locomotion","diet","reproduction") %in% names(coords)),
  all(species_df$iucn %in% names(PARAMS$p_extinction_100yr))
)


# ============================================================================
# 3. Fonctions de calcul ------------------------------------------------------
# ============================================================================

build_compact_supports <- function(tpd_obj) {
  kernels <- tpd_obj$TPDs
  if (is.null(kernels))
    stop("tpd_obj$TPDs est NULL : vérifier la structure de l'objet TPD::TPDs.")
  n_cells <- length(kernels[[1]])
  supports <- lapply(kernels, function(k) which(k > 0))
  list(supports = supports, n_cells = n_cells)
}

compute_FRic_compact <- function(compact, subset_sp) {
  subset_sp <- intersect(subset_sp, names(compact$supports))
  if (length(subset_sp) == 0) return(0)
  occupied <- unique(unlist(compact$supports[subset_sp], use.names = FALSE))
  length(occupied) / compact$n_cells
}

compute_FUn_mean <- function(mat, sp_names, subset_sp, k) {
  idx <- which(sp_names %in% subset_sp)
  if (length(idx) <= k) return(NA_real_)
  m_sub <- mat[idx, , drop = FALSE]
  nn <- RANN::nn2(data = m_sub, query = m_sub, k = k + 1)
  mean(rowMeans(nn$nn.dists[, -1, drop = FALSE]))
}


# ============================================================================
# 3 bis. Standardized Effect Size (SES) --------------------------------------
# ============================================================================
#' SES + p-value de permutation
#'
#' Compare an empirical distribution (obs) to a null distribution (null_vec)
#' via le SES, avec p-value de permutation suivant Davison & Hinkley (1997 :
#' "+1/+1" rule that avoids p = 0 and remains conservative for permutation
#' tests. Convention shared by vegan::permutest, picante::ses.mpd,
#' and most macroecology packages.
#'
#' For the deterministic IUCN AT scenario, `obs` can be a scalar (the single
#' empirique constante) : la fonction le traite alors comme un singleton.
#'
#' SES = (mean(obs) - mean(null)) / sd(null)
#' p_two_sided = 2 * min(p_left, p_right) avec
#'   p_left  = (#{null <= obs_mean} + 1) / (n_null + 1)
#'   p_right = (#{null >= obs_mean} + 1) / (n_null + 1)
#'
#' @param obs        empirical vector (or scalar for the deterministic scenario)
#' @param null_vec   vecteur null (n_iter randomisations)
#' @param alternative "two.sided", "less" ou "greater"
#' @return data.frame 1 ligne avec obs_mean, null_mean, null_sd, SES, p_value,
#'         and n_obs, n_null for traceability
compute_SES <- function(obs, null_vec, alternative = "two.sided") {
  
  obs      <- obs[is.finite(obs)]
  null_vec <- null_vec[is.finite(null_vec)]
  
  if (length(obs) == 0L || length(null_vec) < 2L) {
    return(data.frame(
      obs_mean    = if (length(obs)) mean(obs) else NA_real_,
      obs_median  = if (length(obs)) median(obs) else NA_real_,
      null_mean   = if (length(null_vec)) mean(null_vec) else NA_real_,
      null_sd     = if (length(null_vec)) sd(null_vec)   else NA_real_,
      n_obs       = length(obs),
      n_null      = length(null_vec),
      SES         = NA_real_,
      p_value     = NA_real_,
      alternative = alternative,
      stringsAsFactors = FALSE
    ))
  }
  
  obs_mean  <- mean(obs)
  null_mean <- mean(null_vec)
  null_sd   <- sd(null_vec)
  n_null    <- length(null_vec)
  
  SES <- if (null_sd > 0) (obs_mean - null_mean) / null_sd else NA_real_
  
  p_left  <- (sum(null_vec <= obs_mean) + 1L) / (n_null + 1L)
  p_right <- (sum(null_vec >= obs_mean) + 1L) / (n_null + 1L)
  
  p_value <- switch(
    alternative,
    "two.sided" = min(2 * min(p_left, p_right), 1),
    "less"      = p_left,
    "greater"   = p_right,
    stop("alternative must be 'two.sided', 'less' or 'greater'")
  )
  
  data.frame(
    obs_mean    = obs_mean,
    obs_median  = median(obs),
    null_mean   = null_mean,
    null_sd     = null_sd,
    n_obs       = length(obs),
    n_null      = n_null,
    SES         = SES,
    p_value     = p_value,
    alternative = alternative,
    stringsAsFactors = FALSE
  )
}


# ============================================================================
# 4. Prepare per-dimension compact objects -------- ----------------------------
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
  
  FRic_baseline <- compute_FRic_compact(compact, species_df$species)
  FUn_baseline  <- compute_FUn_mean(mat, species_df$species,
                                    species_df$species, PARAMS$k_neighbours)
  
  size_mb <- as.numeric(object.size(compact)) / 1024^2
  message(sprintf("  Compact supports: %.1f MiB", size_mb))
  
  list(
    compact       = compact,
    mat           = mat,
    sp_names      = species_df$species,
    FRic_baseline = FRic_baseline,
    FUn_baseline  = FUn_baseline
  )
}


# ============================================================================
# 5. One simulation iteration -------- ---------------------------------------------
# ============================================================================
one_iter <- function(i, prep, p_ext, threatened, n_sp, all_sp,
                     precomputed_AT) {
  
  extinct_100  <- rbinom(n_sp, 1, p_ext) == 1
  survivors_100 <- all_sp[!extinct_100]
  n_extinct_100 <- sum(extinct_100)
  
  FRic_AT  <- precomputed_AT$FRic_AT
  FUn_AT   <- precomputed_AT$FUn_AT
  n_extinct_AT <- precomputed_AT$n_extinct_AT
  
  FRic_100 <- compute_FRic_compact(prep$compact, survivors_100)
  FUn_100  <- compute_FUn_mean(prep$mat, prep$sp_names, survivors_100,
                               PARAMS$k_neighbours)
  
  rand_100 <- sample(all_sp, n_sp - n_extinct_100)
  rand_AT  <- sample(all_sp, n_sp - n_extinct_AT)
  FRic_rand_100 <- compute_FRic_compact(prep$compact, rand_100)
  FRic_rand_AT  <- compute_FRic_compact(prep$compact, rand_AT)
  FUn_rand_100  <- compute_FUn_mean(prep$mat, prep$sp_names, rand_100,
                                    PARAMS$k_neighbours)
  FUn_rand_AT   <- compute_FUn_mean(prep$mat, prep$sp_names, rand_AT,
                                    PARAMS$k_neighbours)
  
  base_FRic <- prep$FRic_baseline
  base_FUn  <- prep$FUn_baseline
  
  c(
    n_extinct_100   = n_extinct_100,
    n_extinct_AT    = n_extinct_AT,
    
    loss_FRic_100        = 1 - FRic_100      / base_FRic,
    loss_FRic_AT         = 1 - FRic_AT       / base_FRic,
    loss_FRic_rand_100   = 1 - FRic_rand_100 / base_FRic,
    loss_FRic_rand_AT    = 1 - FRic_rand_AT  / base_FRic,
    
    delta_FUn_100        = (FUn_100      - base_FUn) / base_FUn,
    delta_FUn_AT         = (FUn_AT       - base_FUn) / base_FUn,
    delta_FUn_rand_100   = (FUn_rand_100 - base_FUn) / base_FUn,
    delta_FUn_rand_AT    = (FUn_rand_AT  - base_FUn) / base_FUn
  )
}


# ============================================================================
# 6. Boucle principale --------------------------------------------------------
# ============================================================================
p_ext_vec <- PARAMS$p_extinction_100yr[as.character(species_df$iucn)]
threatened_mask <- species_df$iucn %in% PARAMS$threatened_cats
n_sp_total <- nrow(species_df)
all_species <- species_df$species

if (PARAMS$parallel) {
  future::plan(future::multisession, workers = PARAMS$n_workers)
} else {
  future::plan(future::sequential)
}
progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")

extinction_results <- list()

for (dim_name in c("locomotion", "diet", "reproduction")) {
  
  message("\n========== ", dim_name, " ==========")
  prep <- prep_dim(dim_name)
  message(sprintf("  Baseline FRic = %.4f | Baseline FUn = %.4f",
                  prep$FRic_baseline, prep$FUn_baseline))
  
  survivors_AT <- all_species[!threatened_mask]
  precomputed_AT <- list(
    FRic_AT      = compute_FRic_compact(prep$compact, survivors_AT),
    FUn_AT       = compute_FUn_mean(prep$mat, prep$sp_names, survivors_AT,
                                    PARAMS$k_neighbours),
    n_extinct_AT = sum(threatened_mask)
  )
  message(sprintf("  IUCN AT (déterministe) : n_extinct = %d, FRic_AT = %.4f, FUn_AT = %.4f",
                  precomputed_AT$n_extinct_AT,
                  precomputed_AT$FRic_AT,
                  precomputed_AT$FUn_AT))
  
  progressr::with_progress({
    p <- progressr::progressor(steps = PARAMS$n_iter)
    iter_mat <- future.apply::future_sapply(
      seq_len(PARAMS$n_iter),
      function(i) {
        p()
        one_iter(i, prep,
                 p_ext = p_ext_vec,
                 threatened = threatened_mask,
                 n_sp = n_sp_total,
                 all_sp = all_species,
                 precomputed_AT = precomputed_AT)
      },
      future.seed = TRUE
    )
  })
  
  iter_df <- as.data.frame(t(iter_mat))
  
  
  # ====================================================================
  # 6 bis. SES + p-values empirique vs. randomisation
  # --------------------------------------------------------------------
  # Remplace les wilcox.test (paired pour IUCN 100, one-sample pour
  # deterministic IUCN AT). compute_SES handles both cases
  # uniform: a deterministic scenario simply has obs = scalar.
  # ====================================================================
  ses_FRic_100 <- compute_SES(
    obs         = iter_df$loss_FRic_100,
    null_vec    = iter_df$loss_FRic_rand_100,
    alternative = PARAMS$ses_alternative_FRic
  )
  ses_FRic_AT <- compute_SES(
    obs         = unique(iter_df$loss_FRic_AT),       # scalaire
    null_vec    = iter_df$loss_FRic_rand_AT,
    alternative = PARAMS$ses_alternative_FRic
  )
  ses_FUn_100 <- compute_SES(
    obs         = iter_df$delta_FUn_100,
    null_vec    = iter_df$delta_FUn_rand_100,
    alternative = PARAMS$ses_alternative_FUn
  )
  ses_FUn_AT <- compute_SES(
    obs         = unique(iter_df$delta_FUn_AT),       # scalaire
    null_vec    = iter_df$delta_FUn_rand_AT,
    alternative = PARAMS$ses_alternative_FUn
  )
  
  # Table SES par dimension (4 lignes : FRic/FUn × IUCN 100/AT)
  ses_stats <- dplyr::bind_rows(
    cbind(metric = "FRic loss", scenario = "IUCN 100", ses_FRic_100),
    cbind(metric = "FRic loss", scenario = "IUCN AT",  ses_FRic_AT),
    cbind(metric = "FUn change", scenario = "IUCN 100", ses_FUn_100),
    cbind(metric = "FUn change", scenario = "IUCN AT",  ses_FUn_AT)
  )
  
  # Backward-compatible alias: $pvals remains accessible for downstream scripts
  # qui lirait extinction_results[[dim]]$pvals$p_FRic_100 etc.
  pvals <- list(
    p_FRic_100 = ses_FRic_100$p_value,
    p_FRic_AT  = ses_FRic_AT$p_value,
    p_FUn_100  = ses_FUn_100$p_value,
    p_FUn_AT   = ses_FUn_AT$p_value
  )
  ses_values <- list(
    SES_FRic_100 = ses_FRic_100$SES,
    SES_FRic_AT  = ses_FRic_AT$SES,
    SES_FUn_100  = ses_FUn_100$SES,
    SES_FUn_AT   = ses_FUn_AT$SES
  )
  
  
  # ---- summaries (means and 95% bootstrap CI) ----
  bootstrap_ci <- function(x, n_boot = 1000) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(c(NA, NA))
    bs <- replicate(n_boot, mean(sample(x, replace = TRUE)))
    quantile(bs, c(0.025, 0.975))
  }
  
  summarise_metric <- function(v) {
    c(mean = mean(v, na.rm = TRUE),
      sd   = sd(v, na.rm = TRUE),
      median = median(v, na.rm = TRUE),
      ci_lo = bootstrap_ci(v)[1],
      ci_hi = bootstrap_ci(v)[2])
  }
  
  summary_tbl <- tibble::tibble(
    metric = c("loss_FRic_100", "loss_FRic_AT",
               "loss_FRic_rand_100", "loss_FRic_rand_AT",
               "delta_FUn_100", "delta_FUn_AT",
               "delta_FUn_rand_100", "delta_FUn_rand_AT")
  ) %>%
    rowwise() %>%
    mutate(s = list(summarise_metric(iter_df[[metric]]))) %>%
    tidyr::unnest_wider(s) %>%
    ungroup()
  
  extinction_results[[dim_name]] <- list(
    baseline    = list(FRic = prep$FRic_baseline,
                       FUn  = prep$FUn_baseline),
    iter_df     = iter_df,
    ses_stats   = ses_stats,    # full SES table (4 rows)
    ses_values  = ses_values,   # NOUVEAU : alias scalaire
    pvals       = pvals,        # rétro-compatibilité : juste les p-values
    summary     = summary_tbl
  )
  
  # Console display of SES (useful to monitor during execution)
  message("  --- SES results ---")
  print(format(ses_stats[, c("metric","scenario","obs_mean","null_mean",
                             "null_sd","SES","p_value")],
               digits = 3, scientific = FALSE), row.names = FALSE)
}

future::plan(future::sequential)


# ============================================================================
# 7. Sauvegarde ---------------------------------------------------------------
# ============================================================================
out_path <- file.path(PARAMS$out_dir, PARAMS$out_file)
saveRDS(extinction_results, out_path)
message("\nRésultats sauvegardés dans : ", out_path)


# ============================================================================
# 8. Global summary table (with SES) -------- -----------------------------------
# ============================================================================
synth <- dplyr::bind_rows(lapply(names(extinction_results), function(d) {
  r <- extinction_results[[d]]
  s <- r$ses_stats
  s$dimension <- d
  s
}))

# Reorder columns for readability
synth <- synth[, c("dimension", "scenario", "metric",
                   "obs_mean", "obs_median",
                   "null_mean", "null_sd",
                   "n_obs", "n_null",
                   "SES", "p_value", "alternative")]

message("\n=== Tableau de synthèse global (SES + p-values) ===")
print(format(synth, digits = 3, scientific = FALSE), row.names = FALSE)

write.csv(synth,
          file.path(PARAMS$out_dir, "extinction_scenarios_summary.csv"),
          row.names = FALSE)
# ============================================================================
# F7. Reproducibility -------- --------------------------------------------------------
# ============================================================================
sessionInfo()