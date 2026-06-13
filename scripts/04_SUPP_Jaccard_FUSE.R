#-------------------------------------------------------------------------------
# SUPP — Jaccard curve: FUS vs. FUSE priorities
#
# For species retention proportions from 1% to 100%, computes the Jaccard
# between priority sets defined by FUS_sum3 (aggregated functional uniqueness)
# and FUSE_sum3 (uniqueness weighted by extinction risk).
# Both indices are built as the sum of 3 normalised components
# [0,1] (locomotion, life history, diet).
#
# Prerequisites — fd_df must exist with columns:
#   species, morpho_FNo, lifehistory_FNo, diet_FNo,
#             morpho_FUSE, lifehistory_FUSE, diet_FUSE
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c("dplyr", "ggplot2", "readr", "scales")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))


# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(

  # Proportions tested (1% to 100%, step 1%)
  props   = seq(0.01, 1, by = 0.01),

  out_dir = "tables",
  fig_dir = "figures"
)
dir.create(PARAMS$out_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================================
# 2. USER INPUTS — adjust paths --------------------------------
# ============================================================================
# fd_df must be loaded upstream (from 02_ANALYSE_FunFspFuse.R)
stopifnot(
  exists("fd_df"),
  all(c("species",
        "morpho_FNo", "lifehistory_FNo", "diet_FNo",
        "morpho_FUSE", "lifehistory_FUSE", "diet_FUSE") %in% names(fd_df))
)


# ============================================================================
# 3. Helper functions -------- ----------------------------------------------------
# ============================================================================

# Min-max normalisation to [0, 1]
minmax01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2])) return(rep(NA_real_, length(x)))
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

# Jaccard index between two sets
jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

# Top-p species by score (high score = high priority)
top_set <- function(df, score_col, prop) {
  top_n <- ceiling(nrow(df) * prop)
  df %>%
    slice_max(order_by = .data[[score_col]], n = top_n, with_ties = FALSE) %>%
    pull(species)
}


# ============================================================================
# 4. Build aggregated scores FUS_sum3 and FUSE_sum3 -------- -------------------
# ============================================================================

sum_df <- fd_df %>%
  transmute(
    species,
    morpho_FNo_s       = minmax01(morpho_FNo),
    lifehistory_FNo_s  = minmax01(lifehistory_FNo),
    diet_FNo_s         = minmax01(diet_FNo),
    morpho_FUSE_s      = minmax01(morpho_FUSE),
    lifehistory_FUSE_s = minmax01(lifehistory_FUSE),
    diet_FUSE_s        = minmax01(diet_FUSE)
  ) %>%
  mutate(
    FUS_sum3  = morpho_FNo_s  + lifehistory_FNo_s  + diet_FNo_s,
    FUSE_sum3 = morpho_FUSE_s + lifehistory_FUSE_s + diet_FUSE_s,
    n_nonNA_FUS  = rowSums(!is.na(cbind(morpho_FNo_s, lifehistory_FNo_s, diet_FNo_s))),
    n_nonNA_FUSE = rowSums(!is.na(cbind(morpho_FUSE_s, lifehistory_FUSE_s, diet_FUSE_s)))
  ) %>%
  filter(n_nonNA_FUS == 3, n_nonNA_FUSE == 3) %>%
  select(species, FUS_sum3, FUSE_sum3)


# ============================================================================
# 5. Compute Jaccard curve -------- ------------------------------------------
# ============================================================================

jac_df <- lapply(PARAMS$props, function(p) {
  set_fus  <- top_set(sum_df, "FUS_sum3",  p)
  set_fuse <- top_set(sum_df, "FUSE_sum3", p)
  data.frame(
    top_percent    = p * 100,
    top_n          = ceiling(nrow(sum_df) * p),
    jaccard        = jaccard(set_fus, set_fuse),
    intersection_n = length(intersect(set_fus, set_fuse)),
    union_n        = length(union(set_fus, set_fuse))
  )
}) %>% bind_rows()

write_csv(jac_df,
          file.path(PARAMS$out_dir, "jaccard_FUSsum3_vs_FUSEsum3_top1_to_100.csv"))


# ============================================================================
# 6. Figure — Jaccard curve -------- -----------------------------------------------
# ============================================================================

p_jac <- ggplot(jac_df, aes(x = top_percent, y = jaccard)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(
    title    = "Jaccard overlap between FUS_sum3 and FUSE_sum3 priorities",
    subtitle = "Top sets defined by summed (0–3) indices; components scaled to [0,1]",
    x        = "Top proportion of species retained (%)",
    y        = "Jaccard similarity"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(PARAMS$fig_dir, "jaccard_curve_FUSsum3_vs_FUSEsum3.png"),
       p_jac, width = 8.5, height = 5.5, dpi = 320)
ggsave(file.path(PARAMS$fig_dir, "jaccard_curve_FUSsum3_vs_FUSEsum3.pdf"),
       p_jac, width = 8.5, height = 5.5)

p_jac
