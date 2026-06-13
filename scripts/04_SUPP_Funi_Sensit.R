#-------------------------------------------------------------------------------
# SUPP — Sensitivity of functional uniqueness to trait subsampling

#
# For each functional space (Locomotion, Life history, Diet,
# Combined), randomly resamples a trait subset and
# recomputes the most functionally unique species (top 10%).
# Compares Jaccard overlap between observed and simulated species sets
# to test the robustness of priority sets.
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c("tidyverse", "ade4", "funrar", "ggvenn", "progressr")
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

  top_frac      = 0.10,   # fraction of most unique species retained
  reps          = 100,    # number of resampling iterations

  # minimum number of traits per space (= size of the smallest space)
  n_traits_M    = 5,
  n_traits_L    = 5,
  n_traits_D    = 7,
  n_traits_C    = 14,

  out_dir       = "data/processed",
  fig_dir       = "figures"
)
dir.create(PARAMS$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(PARAMS$fig_dir, showWarnings = FALSE, recursive = TRUE)


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

# Functional uniqueness = mean distance to all other species
compute_uniqueness <- function(dist_matrix) {
  dist_mat <- as.matrix(dist_matrix)
  diag(dist_mat) <- NA
  rowMeans(dist_mat, na.rm = TRUE)
}

# Rescale rows of a fuzzy table to sum to 1
rescale_fuzzy <- function(df) {
  row_sums <- rowSums(df)
  df[row_sums > 0, ] <- df[row_sums > 0, ] / row_sums[row_sums > 0]
  df
}

# Fix all-zero rows (add small constant ε)
fix_all_zero <- function(df) {
  zero_rows <- which(rowSums(df) == 0)
  if (length(zero_rows) > 0) {
    df[zero_rows, ] <- 1e-6
    df <- rescale_fuzzy(df)
  }
  df
}

# Jaccard index between two sets
jaccard_overlap <- function(set1, set2) {
  length(intersect(set1, set2)) / length(union(set1, set2))
}


# ============================================================================
# 4. Observed functional uniqueness -------- ---------------------------------------------
# ============================================================================

diet_fuzzy <- prep.fuzzy(phenoBird[, PARAMS$diet_traits],
                         col.blocks = ncol(phenoBird[, PARAMS$diet_traits]),
                         label = "diet")
diet_fuzzy[diet_fuzzy < 0] <- 0
rownames(diet_fuzzy) <- rownames(phenoBird)

dist_morpho <- dist.ktab(ktab.list.df(list(morpho = phenoBird[, PARAMS$morpho_traits])),
                         type = "Q")
dist_lh     <- dist.ktab(ktab.list.df(list(lifehistory = phenoBird[, PARAMS$lht_traits])),
                         type = "Q")
dist_diet   <- dist.ktab(ktab.list.df(list(diet = diet_fuzzy)),
                         type = "F")
dist_all    <- dist.ktab(
  ktab.list.df(list(morpho = phenoBird[, PARAMS$morpho_traits],
                    lifehistory = phenoBird[, PARAMS$lht_traits],
                    diet = diet_fuzzy)),
  type = c("Q", "Q", "F")
)

fd_df <- tibble(
  species     = rownames(phenoBird),
  morpho      = compute_uniqueness(dist_morpho),
  lifehistory = compute_uniqueness(dist_lh),
  diet        = compute_uniqueness(dist_diet),
  aggregated  = compute_uniqueness(dist_all)
)

top_n <- ceiling(nrow(fd_df) * PARAMS$top_frac)

set_rare_obs <- list(
  M       = fd_df %>% slice_max(morpho,      n = top_n) %>% pull(species),
  L       = fd_df %>% slice_max(lifehistory, n = top_n) %>% pull(species),
  D       = fd_df %>% slice_max(diet,        n = top_n) %>% pull(species),
  Combine = fd_df %>% slice_max(aggregated,  n = top_n) %>% pull(species)
)
saveRDS(set_rare_obs, file.path(PARAMS$out_dir, "set_rare_obs.rds"))


# ============================================================================
# 5. Trait resampling (sensitivity analysis) -------- -------------------
# ============================================================================

resample_uniqueness <- function(data, nn_sp, type, n_traits, reps = 100,
                                top_frac = 0.1,
                                morpho = NULL, lifehistory = NULL, diet = NULL) {

  diet_cols_all <- c("Diet.Inv", "Diet.Vend", "Diet.Vect", "Diet.Vfish",
                     "Diet.Vunk", "Diet.Scav", "Diet.Fruit", "Diet.Nect",
                     "Diet.Seed", "Diet.PlantO")
  res       <- vector("list", reps)
  all_traits <- colnames(data)

  handlers(global = TRUE)
  with_progress({
    p <- progressor(steps = reps)

    for (i in seq_len(reps)) {

      # --- Combined LMD space ---
      if (!is.null(morpho) && !is.null(lifehistory) && !is.null(diet)) {
        pick_m <- sample(morpho, 1)
        pick_l <- sample(lifehistory, 1)
        pick_d <- sample(diet, 2)
        remaining <- setdiff(all_traits, c(pick_m, pick_l, pick_d))
        extra  <- if (n_traits > 4) sample(remaining, n_traits - 4) else c()
        chosen <- c(pick_m, pick_l, pick_d, extra)

        sub    <- data[, chosen, drop = FALSE]
        d_cols <- intersect(colnames(sub), diet_cols_all)
        if (length(d_cols) > 0) {
          sub[, d_cols] <- fix_all_zero(rescale_fuzzy(sub[, d_cols]))
        }
        m_sub <- sub[, intersect(colnames(sub), morpho), drop = FALSE]
        l_sub <- sub[, intersect(colnames(sub), lifehistory), drop = FALSE]
        d_sub <- prep.fuzzy(sub[, intersect(colnames(sub), diet_cols_all)],
                            col.blocks = length(intersect(colnames(sub), diet_cols_all)),
                            label = "diet")
        dist_tmp <- dist.ktab(ktab.list.df(list(morpho = m_sub,
                                                lifehistory = l_sub,
                                                diet = d_sub)),
                              type = c("Q", "Q", "F"))

      } else {
        # --- Single space (M, L, or D) ---
        sub <- data[, sample(ncol(data), n_traits), drop = FALSE]
        if (type == "Q") {
          dist_tmp <- dist.ktab(ktab.list.df(list(sub = sub)), type = "Q")
        } else {
          d_cols <- intersect(colnames(sub), diet_cols_all)
          if (length(d_cols) > 0) {
            sub[, d_cols] <- fix_all_zero(rescale_fuzzy(sub[, d_cols]))
          }
          df_fuzzy <- prep.fuzzy(sub, col.blocks = ncol(sub), label = "diet")
          df_fuzzy[df_fuzzy < 0] <- 0
          dist_tmp <- dist.ktab(ktab.list.df(list(sub = df_fuzzy)), type = "F")
        }
      }

      uniq   <- compute_uniqueness(dist_tmp)
      names(uniq) <- nn_sp
      cutoff <- quantile(uniq, probs = 1 - top_frac)
      res[[i]] <- names(uniq[uniq >= cutoff])
      p()
    }
  })
  res
}

set.seed(42)
res_M <- resample_uniqueness(
  data = phenoBird[, PARAMS$morpho_traits],
  nn_sp = phenoBird$scientificNameStd, type = "Q",
  n_traits = PARAMS$n_traits_M, reps = PARAMS$reps, top_frac = PARAMS$top_frac,
  morpho = PARAMS$morpho_traits, lifehistory = NULL, diet = NULL
)
res_L <- resample_uniqueness(
  data = phenoBird[, PARAMS$lht_traits],
  nn_sp = phenoBird$scientificNameStd, type = "Q",
  n_traits = PARAMS$n_traits_L, reps = PARAMS$reps, top_frac = PARAMS$top_frac,
  morpho = NULL, lifehistory = PARAMS$lht_traits, diet = NULL
)
res_D <- resample_uniqueness(
  data = phenoBird[, PARAMS$diet_traits],
  nn_sp = phenoBird$scientificNameStd, type = "F",
  n_traits = PARAMS$n_traits_D, reps = PARAMS$reps, top_frac = PARAMS$top_frac,
  morpho = NULL, lifehistory = NULL, diet = PARAMS$diet_traits
)
res_C <- resample_uniqueness(
  data = phenoBird[, c(PARAMS$morpho_traits, PARAMS$lht_traits, PARAMS$diet_traits)],
  nn_sp = phenoBird$scientificNameStd, type = c("Q", "Q", "F"),
  n_traits = PARAMS$n_traits_C, reps = PARAMS$reps, top_frac = PARAMS$top_frac,
  morpho = PARAMS$morpho_traits, lifehistory = PARAMS$lht_traits, diet = PARAMS$diet_traits
)

saveRDS(list(M = res_M, L = res_L, D = res_D, C = res_C),
        file.path(PARAMS$out_dir, "resampling_output.rds"))


# ============================================================================
# 6. Observed vs. resampled comparison -------- -----------------------------------
# ============================================================================

compare_overlap <- function(obs1, obs2, res1, res2) {
  obs <- jaccard_overlap(obs1, obs2)
  sim <- map2_dbl(res1, res2, jaccard_overlap)
  tibble(
    observed      = obs,
    mean_resampled = mean(sim),
    lower95       = quantile(sim, 0.025),
    upper95       = quantile(sim, 0.975)
  )
}

pairs <- combn(names(set_rare_obs), 2, simplify = FALSE)

results <- map_dfr(pairs, function(x) {
  compare_overlap(
    set_rare_obs[[x[1]]], set_rare_obs[[x[2]]],
    get(paste0("res_", substr(x[1], 1, 1))),
    get(paste0("res_", substr(x[2], 1, 1)))
  ) %>%
    mutate(pair = paste(x, collapse = "-"))
})
print(results)


# ============================================================================
# 7. Figure — simulated vs. observed Jaccard distributions -------- -------------------
# ============================================================================

resampled_df <- bind_rows(
  tibble(pair = "M-L",       value = map2_dbl(res_M, res_L, jaccard_overlap)),
  tibble(pair = "M-D",       value = map2_dbl(res_M, res_D, jaccard_overlap)),
  tibble(pair = "M-Combine", value = map2_dbl(res_M, res_C, jaccard_overlap)),
  tibble(pair = "L-D",       value = map2_dbl(res_L, res_D, jaccard_overlap)),
  tibble(pair = "L-Combine", value = map2_dbl(res_L, res_C, jaccard_overlap)),
  tibble(pair = "D-Combine", value = map2_dbl(res_D, res_C, jaccard_overlap))
)

obs_df <- results %>% select(pair, observed)

p_uni <- ggplot(resampled_df, aes(x = pair, y = value)) +
  geom_boxplot(fill = "grey82", color = "black", outlier.shape = NA) +
  geom_point(data = obs_df, aes(y = observed), color = "red", size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Trait space comparison",
    y = "Jaccard overlap (rarest 10%)",
    title = "Sensitivity of functional uniqueness overlap to trait resampling"
  ) +
  theme_minimal(base_size = 14)

ggsave(file.path(PARAMS$fig_dir, "Uniqueness_overlap_sensitivity.png"),
       p_uni, width = 8, height = 6, dpi = 300)
