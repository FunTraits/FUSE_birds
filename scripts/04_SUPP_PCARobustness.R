#-------------------------------------------------------------------------------
# SUPP — Robustness of functional spaces (PCoA) to trait subsampling

#
# For each functional space (Locomotion, Life history, Diet,
# Combined LMD), subsamples traits from 3 to n_total and evaluates
# ordination stability via:
#   * Procrustes test (statistic t0: correlation between ordinations)
#   * Mantel correlation between distance matrices
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c("cluster", "ade4", "vegan", "progress",
                   "dplyr", "ggplot2", "readr")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))


# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(

  morpho_traits = c("Tarsus.Length", "Wing.Length", "Kipps.Distance",
                    "Secondary1", "Hand.Wing.Index", "Tail.Length",
                    "Mass", "adult_svl_cm"),

  lht_traits    = c("litter_or_clutch_size_n", "incubation_d", "longevity_y",
                    "fledging_age_d", "litters_or_clutches_per_y"),

  diet_traits   = c("Diet.Inv", "Diet.Vend", "Diet.Vect", "Diet.Vfish",
                    "Diet.Vunk", "Diet.Scav", "Diet.Fruit", "Diet.Nect",
                    "Diet.Seed", "Diet.PlantO"),

  n_iter        = 30,     # iterations per number of traits
  min_traits    = 3,      # minimum number of traits to subsample
  seed          = 123,

  procrustes_threshold = 0.8,   # robustness threshold shown as dashed line

  out_dir       = "data/processed",
  fig_dir       = "results/figures"
)
dir.create(PARAMS$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(PARAMS$fig_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(PARAMS$seed)


# ============================================================================
# 2. USER INPUTS — adjust paths --------------------------------
# ============================================================================
phenoBird <- read_csv("data/processed/phenoBirdsImputedREADY.csv")
rownames(phenoBird) <- phenoBird$scientificNameStd

stopifnot(
  all(PARAMS$morpho_traits %in% names(phenoBird)),
  all(PARAMS$lht_traits    %in% names(phenoBird)),
  all(PARAMS$diet_traits   %in% names(phenoBird))
)


# ============================================================================
# 3. Helper functions -------- ----------------------------------------------------
# ============================================================================

rescale_fuzzy <- function(df) {
  row_sums <- rowSums(df)
  df[row_sums > 0, ] <- df[row_sums > 0, ] / row_sums[row_sums > 0]
  df
}

fix_all_zero <- function(df) {
  zero_rows <- which(rowSums(df) == 0)
  if (length(zero_rows) > 0) {
    df[zero_rows, ] <- 1e-6
    df <- rescale_fuzzy(df)
  }
  df
}

# PCoA from a trait table + distance type(s)
run_pcoa <- function(traits, groups) {
  if (length(groups) > 1) {
    diet_f <- prep.fuzzy(
      traits[, colnames(traits) %in% PARAMS$diet_traits],
      col.blocks = sum(colnames(traits) %in% PARAMS$diet_traits),
      label = "diet"
    )
    diet_f[diet_f < 0] <- 0
    ktabList <- ktab.list.df(list(
      morpho = traits[, colnames(traits) %in% PARAMS$morpho_traits],
      lh     = traits[, colnames(traits) %in% PARAMS$lht_traits],
      diet   = diet_f
    ))
  } else if (groups == "F") {
    diet_f <- prep.fuzzy(
      traits[, colnames(traits) %in% PARAMS$diet_traits],
      col.blocks = sum(colnames(traits) %in% PARAMS$diet_traits),
      label = "diet"
    )
    diet_f[diet_f < 0] <- 0
    ktabList <- ktab.list.df(list(diet_f))
  } else {
    ktabList <- ktab.list.df(list(traits))
  }

  dis    <- dist.ktab(ktabList, groups, scan = FALSE, option = "scaledBYrange")
  pcoa_r <- cmdscale(dis, k = nrow(traits) - 1, eig = TRUE)
  list(scores = pcoa_r$points, dist = dis)
}


# ============================================================================
# 4. Subsampling — single space (M, L, or D) -------- ------------------------
# ============================================================================

pcoa_trait_subsampling <- function(traits, groups, space_name,
                                   n_iter = 50, min_traits = 3) {
  n_total     <- ncol(traits)
  if (n_total < min_traits) stop(paste("Not enough traits in", space_name))
  full        <- run_pcoa(traits, groups)
  trait_counts <- seq(min_traits, n_total)
  pb <- progress_bar$new(
    format = paste0("⏳ ", space_name, " [:bar] :percent | ETA: :eta"),
    total  = length(trait_counts) * n_iter, clear = FALSE, width = 60
  )
  results <- list()
  for (n in trait_counts) {
    for (i in seq_len(n_iter)) {
      chosen    <- sample(colnames(traits), n)
      sub       <- traits[, chosen, drop = FALSE]
      if (groups == "F") sub <- fix_all_zero(rescale_fuzzy(sub))
      sub_pcoa  <- run_pcoa(sub, groups)
      proc      <- protest(full$scores, sub_pcoa$scores, permutations = 0)
      mantel_r  <- mantel(full$dist, sub_pcoa$dist, permutations = 0)
      results[[length(results) + 1]] <- data.frame(
        space    = space_name, n_traits = n, iteration = i,
        proc_stat = proc$t0, mantel_r = mantel_r$statistic
      )
      pb$tick()
    }
  }
  bind_rows(results)
}


# ============================================================================
# 5. Subsampling — combined LMD space -------- --------------------------------
# ============================================================================

pcoa_trait_subsampling_combined <- function(data, morpho_traits, lht_traits,
                                             diet_traits, n_iter = 50,
                                             min_traits = 3) {
  all_traits <- c(morpho_traits, lht_traits, diet_traits)
  if (length(all_traits) < min_traits) stop("Not enough traits for combined space")

  diet_full <- fix_all_zero(rescale_fuzzy(data[, diet_traits]))
  all_mat   <- cbind(data[, morpho_traits], data[, lht_traits], diet_full)
  full      <- run_pcoa(all_mat, groups = c("Q", "Q", "F"))

  trait_counts <- seq(min_traits, length(all_traits))
  pb <- progress_bar$new(
    format = "⏳ LMD [:bar] :percent | ETA: :eta",
    total  = length(trait_counts) * n_iter, clear = FALSE, width = 60
  )
  results <- list()
  for (n in trait_counts) {
    for (i in seq_len(n_iter)) {
      pick_m    <- sample(morpho_traits, 1)
      pick_l    <- sample(lht_traits, 1)
      pick_d    <- sample(diet_traits, 2)
      remaining <- setdiff(all_traits, c(pick_m, pick_l, pick_d))
      extra     <- if (n > 4) sample(remaining, n - 4) else c()
      chosen    <- c(pick_m, pick_l, pick_d, extra)
      sub       <- data[, chosen, drop = FALSE]
      d_cols    <- intersect(colnames(sub), diet_traits)
      if (length(d_cols) > 0) sub[, d_cols] <- fix_all_zero(rescale_fuzzy(sub[, d_cols]))
      sub_pcoa  <- run_pcoa(sub, groups = c("Q", "Q", "F"))
      proc      <- protest(full$scores, sub_pcoa$scores, permutations = 0)
      mantel_r  <- mantel(full$dist, sub_pcoa$dist, permutations = 0)
      results[[length(results) + 1]] <- data.frame(
        space    = "LMD", n_traits = n, iteration = i,
        proc_stat = proc$t0, mantel_r = mantel_r$statistic
      )
      pb$tick()
    }
  }
  bind_rows(results)
}


# ============================================================================
# 6. Run analyses -------- ---------------------------------------------------
# ============================================================================

# Option A: compute from scratch
results_raw <- bind_rows(
  pcoa_trait_subsampling(phenoBird[, PARAMS$morpho_traits], "Q", "Locomotion",
                         PARAMS$n_iter, PARAMS$min_traits),
  pcoa_trait_subsampling(phenoBird[, PARAMS$lht_traits],    "Q", "Reproduction",
                         PARAMS$n_iter, PARAMS$min_traits),
  pcoa_trait_subsampling(phenoBird[, PARAMS$diet_traits],   "F", "Diet",
                         PARAMS$n_iter, PARAMS$min_traits),
  pcoa_trait_subsampling_combined(phenoBird,
                                  PARAMS$morpho_traits, PARAMS$lht_traits,
                                  PARAMS$diet_traits,
                                  PARAMS$n_iter, min_traits = 4)
)
saveRDS(results_raw,
        file.path(PARAMS$out_dir, "pcoa_trait_subsampling_results.rds"))

# Option B: load pre-computed results (uncomment if needed)
# results <- bind_rows(
#   read_csv(file.path(PARAMS$out_dir, "pcoa_robustness_results_M.csv"))   %>% mutate(space = "Locomotion"),
#   read_csv(file.path(PARAMS$out_dir, "pcoa_robustness_results_LHT.csv")) %>% mutate(space = "Reproduction"),
#   read_csv(file.path(PARAMS$out_dir, "pcoa_robustness_results_D.csv"))   %>% mutate(space = "Diet"),
#   read_csv(file.path(PARAMS$out_dir, "pcoa_robustness_results_LMD.csv")) %>% mutate(space = "Combined")
# )


# ============================================================================
# 7. Figure — Procrustes robustness vs. number of traits -------- ----------
# ============================================================================

p_rob <- ggplot(results_raw,
                aes(x = factor(n_traits), y = proc_stat, fill = space)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_hline(yintercept = PARAMS$procrustes_threshold,
             linetype = "dashed", color = "red", linewidth = 1) +
  labs(x     = "Number of traits",
       y     = "Procrustes correlation",
       fill  = "Trait space",
       title = "") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")

ggsave(file.path(PARAMS$fig_dir, "PCoA_robustness.png"),
       p_rob, width = 8, height = 6, dpi = 300)
