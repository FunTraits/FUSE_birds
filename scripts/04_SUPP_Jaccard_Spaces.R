#-------------------------------------------------------------------------------
# SUPP — Cross-space concordance of functionally distinct species
#
# Tests whether species identified as functionally distinct are consistent
# across trait spaces (Locomotion [M], Life History [L], Diet [D],
# Combined [C]).
#
# 1. Pairwise Jaccard similarity between spaces based on the top 10%
#    most functionally unique species (computed from full trait sets).
# 2. Classification of distinct species as exclusive (identified in
#    exactly one trait space) or shared (identified in two or more).
# 3. Taxonomic representation: families with exclusive vs. shared
#    distinct species to identify space-dependent distinctiveness.
# 4. Resampling sensitivity (n = 1 000): for each pair of spaces,
#    random trait subsets are drawn independently within each space,
#    functional uniqueness is recomputed, and cross-space Jaccard is
#    calculated. The resulting null distribution is compared to the
#    observed Jaccard to test whether concordance/discordance is
#    driven by trait choice alone.
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c("dplyr", "tibble", "purrr", "tidyr", "forcats",
                   "ggplot2", "patchwork",
                   "ade4", "future", "future.apply", "readr")
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

  top_frac      = 0.10,   # fraction defining functionally distinct species
  reps          = 50,   # resampling iterations per pair
  seed          = 42,

  # Fixed subsample sizes per space (resampling sensitivity)
  n_sub_M       = 5,      # traits drawn from Locomotion
  n_sub_L       = 3,      # traits drawn from Life History
  n_sub_D       = 6,      # traits drawn from Diet
  # Combined: at least 1 M + 1 L + 2 D, remainder drawn randomly to reach n_sub_C
  n_sub_C       = 10,

  out_dir       = "tables",
  fig_dir       = "figures"
)
dir.create(PARAMS$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(PARAMS$fig_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(PARAMS$seed)


# ============================================================================
# 2. USER INPUTS — adjust paths -----------------------------------------------
# ============================================================================
phenoBird <- read_csv("data/processed/phenoBirdsImputedREADY.csv")
rownames(phenoBird) <- phenoBird$scientificNameStd

taxo <- read.csv("data/processed/Taxo_Birds.csv")   # must contain: species, family
names(taxo)[1] = "species"

stopifnot(
  all(PARAMS$morpho_traits %in% names(phenoBird)),
  all(PARAMS$lht_traits    %in% names(phenoBird)),
  all(PARAMS$diet_traits   %in% names(phenoBird)),
  all(c("species", "family") %in% names(taxo))
)


# ============================================================================
# 3. Helper functions ---------------------------------------------------------
# ============================================================================

# Mean distance to all other species (functional uniqueness sensu Mouillot 2013)
compute_uniqueness <- function(dist_matrix) {
  m <- as.matrix(dist_matrix)
  diag(m) <- NA
  rowMeans(m, na.rm = TRUE)
}

# Row-normalise a fuzzy-coded matrix (rows sum to 1)
rescale_fuzzy <- function(df) {
  rs <- rowSums(df)
  df[rs > 0, ] <- df[rs > 0, ] / rs[rs > 0]
  df
}

# Replace all-zero rows with uniform epsilon, then renormalise
fix_all_zero <- function(df) {
  zero <- which(rowSums(df) == 0)
  if (length(zero) > 0) { df[zero, ] <- 1e-6; df <- rescale_fuzzy(df) }
  df
}

# Jaccard similarity between two species sets
jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

# Return names of the top-frac species by functional uniqueness
top_distinct <- function(uniq_vec, top_frac = 0.10) {
  n <- ceiling(length(uniq_vec) * top_frac)
  names(sort(uniq_vec, decreasing = TRUE))[seq_len(n)]
}

# Build a Gower distance matrix from a trait data frame.
#   groups      : "Q", "F", or c("Q","Q","F") for combined
#   diet_cols   : column names belonging to the fuzzy (diet) block
#   morpho_cols : column names in the locomotion block (combined only)
#   lht_cols    : column names in the life-history block (combined only)
build_dist <- function(trait_mat,
                       groups,
                       diet_cols   = NULL,
                       morpho_cols = NULL,
                       lht_cols    = NULL) {
  trait_mat <- as.data.frame(trait_mat)

  if (length(groups) > 1) {
    # Combined space: three blocks M + L + D
    dc <- intersect(colnames(trait_mat), diet_cols)
    mc <- intersect(colnames(trait_mat), morpho_cols)
    lc <- intersect(colnames(trait_mat), lht_cols)

    d_sub <- fix_all_zero(rescale_fuzzy(trait_mat[, dc, drop = FALSE]))
    d_fuz <- prep.fuzzy(d_sub, col.blocks = ncol(d_sub), label = "diet")
    d_fuz[d_fuz < 0] <- 0

    kt <- ktab.list.df(list(
      morpho = trait_mat[, mc, drop = FALSE],
      lht    = trait_mat[, lc, drop = FALSE],
      diet   = d_fuz
    ))
    dist.ktab(kt, type = groups, scan = FALSE, option = "scaledBYrange")

  } else if (groups == "F") {
    d_sub <- fix_all_zero(rescale_fuzzy(trait_mat))
    d_fuz <- prep.fuzzy(d_sub, col.blocks = ncol(d_sub), label = "diet")
    d_fuz[d_fuz < 0] <- 0
    kt <- ktab.list.df(list(diet = d_fuz))
    dist.ktab(kt, type = "F", scan = FALSE, option = "scaledBYrange")

  } else {
    kt <- ktab.list.df(list(traits = trait_mat))
    dist.ktab(kt, type = "Q", scan = FALSE, option = "scaledBYrange")
  }
}


# ============================================================================
# 4. Observed functional uniqueness and top-distinct sets ---------------------
# ============================================================================

pheno_df  <- as.data.frame(phenoBird)
diet_full <- fix_all_zero(rescale_fuzzy(pheno_df[, PARAMS$diet_traits]))
rownames(diet_full) <- rownames(pheno_df)

data_M <- pheno_df[, PARAMS$morpho_traits]
data_L <- pheno_df[, PARAMS$lht_traits]
data_D <- diet_full
data_C <- cbind(data_M, data_L, data_D)

dist_M <- build_dist(data_M, "Q")
dist_L <- build_dist(data_L, "Q")
dist_D <- build_dist(data_D, "F", diet_cols = PARAMS$diet_traits)
dist_C <- build_dist(data_C, c("Q", "Q", "F"),
                     diet_cols   = PARAMS$diet_traits,
                     morpho_cols = PARAMS$morpho_traits,
                     lht_cols    = PARAMS$lht_traits)

uniq_M <- compute_uniqueness(dist_M)
uniq_L <- compute_uniqueness(dist_L)
uniq_D <- compute_uniqueness(dist_D)
uniq_C <- compute_uniqueness(dist_C)

sets_obs <- list(
  M = top_distinct(uniq_M, PARAMS$top_frac),
  L = top_distinct(uniq_L, PARAMS$top_frac),
  D = top_distinct(uniq_D, PARAMS$top_frac),
  C = top_distinct(uniq_C, PARAMS$top_frac)
)
cat(sprintf(
  "Distinct species per space:  M = %d  |  L = %d  |  D = %d  |  C = %d\n",
  length(sets_obs$M), length(sets_obs$L),
  length(sets_obs$D), length(sets_obs$C)
))


# ============================================================================
# 5. Cross-space Jaccard — observed -------------------------------------------
# ============================================================================

space_pairs <- combn(names(sets_obs), 2, simplify = FALSE)

jac_obs <- purrr::map_dfr(space_pairs, function(pair) {
  tibble(
    space_A = pair[1], space_B = pair[2],
    jaccard = jaccard(sets_obs[[pair[1]]], sets_obs[[pair[2]]]),
    n_inter = length(intersect(sets_obs[[pair[1]]], sets_obs[[pair[2]]])),
    n_union = length(union(sets_obs[[pair[1]]], sets_obs[[pair[2]]])),
    pair    = paste(pair, collapse = "-")
  )
})
print(jac_obs)
write.csv(jac_obs,
          file.path(PARAMS$out_dir, "jaccard_cross_space_observed.csv"),
          row.names = FALSE)


# ============================================================================
# 6. Exclusive vs. shared species ---------------------------------------------
# ============================================================================

all_distinct <- Reduce(union, sets_obs)

membership <- tibble(species = all_distinct) %>%
  mutate(
    in_M     = species %in% sets_obs$M,
    in_L     = species %in% sets_obs$L,
    in_D     = species %in% sets_obs$D,
    in_C     = species %in% sets_obs$C,
    n_spaces = in_M + in_L + in_D + in_C,
    category = if_else(n_spaces == 1, "exclusive", "shared"),
    excl_space = case_when(
      n_spaces == 1 & in_M ~ "M",
      n_spaces == 1 & in_L ~ "L",
      n_spaces == 1 & in_D ~ "D",
      n_spaces == 1 & in_C ~ "C",
      TRUE ~ NA_character_
    )
  )

cat(sprintf(
  "Total distinct species: %d\n  Exclusive: %d (%.1f%%)\n  Shared:    %d (%.1f%%)\n",
  nrow(membership),
  sum(membership$category == "exclusive"),
  100 * mean(membership$category == "exclusive"),
  sum(membership$category == "shared"),
  100 * mean(membership$category == "shared")
))
write.csv(membership,
          file.path(PARAMS$out_dir, "distinct_species_membership.csv"),
          row.names = FALSE)


# ============================================================================
# 7. Taxonomic representation -------------------------------------------------
# ============================================================================

membership_tax <- membership %>%
  left_join(taxo %>% select(species, family), by = "species")

# Per-family summary
family_summary <- membership_tax %>%
  group_by(family) %>%
  summarise(
    n_distinct  = n(),
    n_exclusive = sum(category == "exclusive"),
    n_shared    = sum(category == "shared"),
    pct_excl    = round(100 * n_exclusive / n_distinct, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_distinct))

# Exclusive species count per family × space
excl_by_space_family <- membership_tax %>%
  filter(category == "exclusive") %>%
  count(family, excl_space, name = "n") %>%
  pivot_wider(names_from = excl_space, values_from = n, values_fill = 0L)

write.csv(family_summary,
          file.path(PARAMS$out_dir, "family_distinctiveness_summary.csv"),
          row.names = FALSE)
write.csv(excl_by_space_family,
          file.path(PARAMS$out_dir, "family_exclusive_by_space.csv"),
          row.names = FALSE)


# ============================================================================
# 8. Resampling sensitivity (n = 1 000 iterations) ---------------------------
# ============================================================================
# For each pair of spaces (A, B):
#   Draw random trait subsets independently within A and within B
#   (preserving trait type for the combined space).
#   Recompute functional uniqueness → identify top 10% → compute Jaccard.
#   The resulting distribution is the null expectation for cross-space
#   Jaccard driven solely by trait choice, against which observed values
#   are compared (outside 95% CI = robust concordance/discordance).

# Internal subsampler — draws a random trait subset for one space
draw_subset <- function(data, groups, n_sub,
                        morpho_cols = NULL, lht_cols = NULL, diet_cols = NULL) {
  if (length(groups) > 1) {
    # Combined: guarantee at least 1 M + 1 L + 2 D
    mc <- intersect(colnames(data), morpho_cols)
    lc <- intersect(colnames(data), lht_cols)
    dc <- intersect(colnames(data), diet_cols)
    fixed  <- c(sample(mc, 1), sample(lc, 1), sample(dc, 2))
    rest   <- setdiff(colnames(data), fixed)
    extra  <- if (n_sub > 4 && length(rest) > 0)
                sample(rest, min(n_sub - 4, length(rest)))
              else character(0)
    data[, c(fixed, extra), drop = FALSE]
  } else {
    data[, sample(ncol(data), min(n_sub, ncol(data))), drop = FALSE]
  }
}

# Jaccard from one pair of random subsamples
one_resample <- function(data_a, groups_a, n_sub_a,
                         data_b, groups_b, n_sub_b,
                         morpho_cols, lht_cols, diet_cols, top_frac) {
  tryCatch({
    sub_a <- draw_subset(data_a, groups_a, n_sub_a, morpho_cols, lht_cols, diet_cols)
    sub_b <- draw_subset(data_b, groups_b, n_sub_b, morpho_cols, lht_cols, diet_cols)

    d_a <- build_dist(sub_a, groups_a, diet_cols, morpho_cols, lht_cols)
    d_b <- build_dist(sub_b, groups_b, diet_cols, morpho_cols, lht_cols)

    set_a <- top_distinct(compute_uniqueness(d_a), top_frac)
    set_b <- top_distinct(compute_uniqueness(d_b), top_frac)

    jaccard(set_a, set_b)
  }, error = function(e) NA_real_)
}

plan(multisession)   # parallel execution — adjust workers with workers = N

run_resampling <- function(data_a, groups_a, n_sub_a,
                           data_b, groups_b, n_sub_b,
                           morpho_cols, lht_cols, diet_cols,
                           top_frac, reps, seed) {
  future_sapply(
    seq_len(reps),
    function(i) one_resample(data_a, groups_a, n_sub_a,
                             data_b, groups_b, n_sub_b,
                             morpho_cols, lht_cols, diet_cols, top_frac),
    future.seed = seed
  )
}

mc  <- PARAMS$morpho_traits
lc  <- PARAMS$lht_traits
dc  <- PARAMS$diet_traits
tfr <- PARAMS$top_frac
rps <- PARAMS$reps
sd  <- PARAMS$seed

message("Running resampling (1 000 iterations × 6 pairs) ...")

resamp_results <- list(
  "M-L" = run_resampling(data_M, "Q",             PARAMS$n_sub_M,
                          data_L, "Q",             PARAMS$n_sub_L,
                          mc, lc, dc, tfr, rps, sd),
  "M-D" = run_resampling(data_M, "Q",             PARAMS$n_sub_M,
                          data_D, "F",             PARAMS$n_sub_D,
                          mc, lc, dc, tfr, rps, sd),
  "M-C" = run_resampling(data_M, "Q",             PARAMS$n_sub_M,
                          data_C, c("Q","Q","F"),  PARAMS$n_sub_C,
                          mc, lc, dc, tfr, rps, sd),
  "L-D" = run_resampling(data_L, "Q",             PARAMS$n_sub_L,
                          data_D, "F",             PARAMS$n_sub_D,
                          mc, lc, dc, tfr, rps, sd),
  "L-C" = run_resampling(data_L, "Q",             PARAMS$n_sub_L,
                          data_C, c("Q","Q","F"),  PARAMS$n_sub_C,
                          mc, lc, dc, tfr, rps, sd),
  "D-C" = run_resampling(data_D, "F",             PARAMS$n_sub_D,
                          data_C, c("Q","Q","F"),  PARAMS$n_sub_C,
                          mc, lc, dc, tfr, rps, sd)
)

plan(sequential)

# Summarise null distributions and compare to observed
resamp_summary <- purrr::imap_dfr(resamp_results, function(vals, pair) {
  v <- vals[!is.na(vals)]
  tibble(
    pair     = pair,
    n_valid  = length(v),
    mean_jac = mean(v),
    lower95  = quantile(v, 0.025),
    upper95  = quantile(v, 0.975)
  )
}) %>%
  left_join(jac_obs %>% select(pair, obs_jaccard = jaccard), by = "pair") %>%
  mutate(
    direction = case_when(
      obs_jaccard > upper95 ~ "higher than null",
      obs_jaccard < lower95 ~ "lower than null",
      TRUE                  ~ "within null"
    )
  )

print(resamp_summary %>% select(pair, obs_jaccard, lower95, upper95, direction))

saveRDS(resamp_results,
        file.path(PARAMS$out_dir, "resampling_cross_space_jaccard.rds"))
write.csv(resamp_summary %>% select(-n_valid),
          file.path(PARAMS$out_dir, "resampling_cross_space_summary.csv"),
          row.names = FALSE)


# ============================================================================
# 9. Figures ------------------------------------------------------------------
# ============================================================================

space_labels <- c(M = "Locomotion", L = "Life History", D = "Diet", C = "Combined")
space_cols   <- c(M = "#4393C3",    L = "#F4A582",      D = "#74C476", C = "#9970AB")

# ── Fig 1: Jaccard matrix (observed) ─────────────────────────────────────────
jac_grid <- bind_rows(
  jac_obs %>% select(x = space_A, y = space_B, jaccard),
  jac_obs %>% select(x = space_B, y = space_A, jaccard),
  tibble(x = names(sets_obs), y = names(sets_obs), jaccard = 1)
) %>%
  mutate(across(c(x, y), ~ factor(., levels = c("M", "L", "D", "C"),
                                    labels = space_labels)))

p_matrix <- ggplot(jac_grid, aes(x = x, y = y, fill = jaccard)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", jaccard)), size = 4.5, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#D6604D",
                       midpoint = 0.5, limits = c(0, 1), name = "Jaccard") +
  labs(
    title = "Pairwise Jaccard similarity between trait spaces",
    subtitle = sprintf("Top %.0f%% most functionally unique species", PARAMS$top_frac * 100),
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        plot.title = element_text(face = "bold"))

ggsave(file.path(PARAMS$fig_dir, "jaccard_matrix_observed.png"),
       p_matrix, width = 6, height = 5.5, dpi = 300)


# ── Fig 2: Exclusive vs. shared bar chart ────────────────────────────────────
memb_plot <- membership %>%
  mutate(
    label = if_else(category == "shared", "Shared (≥2 spaces)",
                    paste0("Exclusive — ", space_labels[excl_space])),
    label = factor(label, levels = c(
      "Exclusive — Locomotion", "Exclusive — Life History",
      "Exclusive — Diet",       "Exclusive — Combined",
      "Shared (≥2 spaces)"
    ))
  )

fill_vals <- c(
  "Exclusive — Locomotion"   = space_cols["M"],
  "Exclusive — Life History" = space_cols["L"],
  "Exclusive — Diet"         = space_cols["D"],
  "Exclusive — Combined"     = space_cols["C"],
  "Shared (≥2 spaces)"       = "grey60"
)

p_bar <- ggplot(memb_plot, aes(x = label, fill = label)) +
  geom_bar(color = "white", linewidth = 0.4) +
  scale_fill_manual(values = fill_vals, guide = "none") +
  labs(
    title    = "Exclusive vs. shared functionally distinct species",
    subtitle = sprintf("n = %d distinct species across all spaces", nrow(membership)),
    x = NULL, y = "Number of species"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold"))

ggsave(file.path(PARAMS$fig_dir, "distinct_species_exclusive_shared.png"),
       p_bar, width = 6.5, height = 5, dpi = 300)


# ── Fig 3: Top families with exclusive distinct species ──────────────────────
top20_fam <- family_summary %>%
  filter(n_exclusive > 0) %>%
  slice_max(n_exclusive, n = 20) %>%
  mutate(family = fct_reorder(family, n_exclusive))

p_fam <- ggplot(top20_fam, aes(x = family, y = n_exclusive, fill = pct_excl)) +
  geom_col(color = NA) +
  coord_flip() +
  scale_fill_gradient(low = "#C7E9C0", high = "#005A32",
                      name = "% exclusive\nspecies") +
  labs(
    title = "Families with the most exclusive distinct species",
    x = NULL, y = "Number of exclusive distinct species"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(PARAMS$fig_dir, "family_exclusive_distinct.png"),
       p_fam, width = 7.5, height = 6.5, dpi = 300)


# ── Fig 4: Resampling null distributions vs. observed Jaccard ─────────────────
null_long <- purrr::imap_dfr(resamp_results, function(vals, pair) {
  tibble(pair = pair, jaccard = vals[!is.na(vals)])
}) %>%
  left_join(resamp_summary %>% select(pair, obs_jaccard, direction), by = "pair")

p_resamp <- ggplot(null_long, aes(x = pair, y = jaccard)) +
  geom_boxplot(fill = "grey88", color = "grey40",
               outlier.shape = NA, alpha = 0.9, width = 0.6) +
  geom_point(
    data    = resamp_summary,
    mapping = aes(x = pair, y = obs_jaccard, color = direction),
    size = 4, shape = 18
  ) +
  geom_errorbar(
    data    = resamp_summary,
    mapping = aes(x = pair, ymin = lower95, ymax = upper95),
    width = 0.2, color = "grey40", linewidth = 0.6
  ) +
  scale_color_manual(
    values = c("higher than null" = "#D6604D",
               "lower than null"  = "#2166AC",
               "within null"      = "black"),
    name = "Observed vs. null"
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title    = "Cross-space Jaccard: observed vs. resampling null",
    subtitle = sprintf(
      "Boxplot = null distribution (n = %d); diamond = observed (full trait set)\nError bars = 95%% CI of null; colour indicates position relative to CI",
      PARAMS$reps
    ),
    x = "Trait space pair",
    y = sprintf("Jaccard similarity (top %.0f%% distinct species)", PARAMS$top_frac * 100)
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(PARAMS$fig_dir, "resampling_cross_space_jaccard.pdf"),
       p_resamp, width = 9, height = 6)

p_resamp
