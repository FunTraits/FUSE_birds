#-------------------------------------------------------------------------------
# Figure S3. Aggregated vs. disaggregated FUSE.
#
# Shows that aggregating traits into a single functional space shifts the
# identity of priority species identified by FUSE, compared to the
# multi-space disaggregated approach (FUSE_global = sum of FUSE per dimension).
#
# Three panels:
#   (a) Scatter plot of FUSE_aggregated vs. FUSE_disaggregated_global ranks.
#       Coloured by IUCN status, 1:1 line, top-25 highlighted.
#   (b) 2-set Venn diagram: top-25 species in each approach.
#   (c) Identity of top-25 species in each approach, side by side.
#
# Author : A. Toussaint (CRBE/CNRS, Toulouse)
#
# Prerequisites:
#   - metrics_df : tibble produced by A1_recompute_FUn_FSp_FUSE.R, with
#                  FUSE_loco, FUSE_diet, FUSE_repro, iucn, GE
#   - PCA_Birds_LMD : RDS object produced by PCoAGraph() upstream, containing
#                  the aggregated functional space across the three trait sets
#                  (LHT + morpho + diet) via ade4::dist.ktab.
#                  Loaded via PARAMS$pca_birds_path.
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c(
  "ggplot2", "dplyr", "tidyr", "tibble", "purrr", "scales",
  "RANN", "patchwork"
)
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# Polyfill %||%
if (!exists("%||%", envir = baseenv())) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  assign("%||%", `%||%`, envir = .GlobalEnv)
}


# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(
  
  # Aggregation metric for global disaggregated FUSE (sum or max)
  disaggregated_combine = "sum",   # "sum" ou "max"
  
  # Number of PCoA axes for the aggregated space
  # (consistent with dimensional spaces in the manuscript: 2 axes)
  n_axes_aggregated = 2,
  
  # k pour FUn (Pimiento 2020)
  k_neighbours      = 5,
  
  # Path to the pre-built aggregated PCA object (PCoAGraph)
  pca_birds_path    = "data/processed/PCA_Birds.rds",
  
  # Mapping IUCN -> GE
  iucn_to_GE        = c(LC = 0, NT = 1, VU = 2, EN = 3, CR = 4),
  iucn_levels       = c("LC", "NT", "VU", "EN", "CR"),
  iucn_colors       = c(LC = "#60C659", NT = "#CCE226",
                        VU = "#F9E814", EN = "#FC7F3F", CR = "#D81E05"),
  fuse_scale        = 4,
  
  # Top-N
  top_n             = 25,
  
  # Apparence
  base_font_size    = 8,
  
  # Output
  out_dir           = "figures",
  fig_basename      = "FigS3_FUSE_aggregated_vs_disaggregated",
  width_mm          = 180,
  height_mm         = 220,
  dpi               = 600,
  
  data_out_dir      = "data/processed",
  data_out_file     = "FigS3_aggregated_vs_disaggregated.rds"
)
dir.create(PARAMS$out_dir,      showWarnings = FALSE, recursive = TRUE)
dir.create(PARAMS$data_out_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================================
# 2. USER INPUTS ------------------------------------------------------
# ============================================================================
metrics_df = readRDS("data/processed/species_metrics_FUn_FSp_FUSE.rds")
traits_full <- read.csv("data/processed/phenoBirdsImputedREADY.csv")
colnames(traits_full)[1] = 'species'
colnames(traits_full)[3] = 'order'
colnames(traits_full)[67] = 'iucn'

stopifnot(
  exists("metrics_df"),
  all(c("species", "iucn", "GE",
        "FUSE_loco", "FUSE_diet", "FUSE_repro") %in% names(metrics_df))
)

metrics_df <- metrics_df %>% arrange(species)


# ============================================================================
# 3. Load aggregated functional space (pre-built) -------- ---------------
# ============================================================================
PCA_Birds_LMD <- readRDS(PARAMS$pca_birds_path)$LMD

stopifnot(
  is.list(PCA_Birds_LMD),
  "PCoA" %in% names(PCA_Birds_LMD),
  "vectors" %in% names(PCA_Birds_LMD$PCoA),
  "values"  %in% names(PCA_Birds_LMD$PCoA)
)

# Extract PCoA coordinates and explained variance
agg_vectors_full <- PCA_Birds_LMD$PCoA$vectors
agg_eigenvalues  <- PCA_Birds_LMD$PCoA$values$Eigenvalues
eig_pos <- agg_eigenvalues[agg_eigenvalues > 0]
var_explained <- eig_pos / sum(eig_pos)

n_axes_avail <- ncol(agg_vectors_full)
n_axes <- min(PARAMS$n_axes_aggregated, n_axes_avail)
if (n_axes < PARAMS$n_axes_aggregated) {
  message(sprintf(
    "ATTENTION : seulement %d axes disponibles dans PCA_Birds_LMD (%d demandés).",
    n_axes_avail, PARAMS$n_axes_aggregated
  ))
}

message(sprintf(
  "Espace agrégé chargé : %d espèces x %d axes utilisés.",
  nrow(agg_vectors_full), n_axes
))
message(sprintf(
  "Variance expliquée par les %d premiers axes : %s ; cumul = %.1f %%",
  n_axes,
  paste(sprintf("%.1f%%",
                100 * var_explained[1:n_axes]),
        collapse = ", "),
  100 * sum(var_explained[1:n_axes])
))

# Conversion en data.frame avec colonne species
agg_coords <- as.data.frame(agg_vectors_full[, 1:n_axes, drop = FALSE])
colnames(agg_coords) <- paste0("PC", 1:n_axes)
agg_coords$species <- rownames(agg_vectors_full)

# Alignement avec metrics_df via species
common_sp <- intersect(metrics_df$species, agg_coords$species)
if (length(common_sp) < nrow(metrics_df)) {
  message(sprintf(
    "Filtrage à %d espèces communes (sur %d en metrics_df, %d en PCoA agrégé).",
    length(common_sp), nrow(metrics_df), nrow(agg_coords)
  ))
  metrics_df <- metrics_df %>% filter(species %in% common_sp)
}
agg_coords <- agg_coords[match(metrics_df$species, agg_coords$species), ,
                         drop = FALSE]
stopifnot(identical(agg_coords$species, metrics_df$species))


# ============================================================================
# 4. Compute aggregated FUSE -------- ----------------------------------------------------
# ============================================================================

mat_agg <- as.matrix(agg_coords[, paste0("PC", 1:n_axes), drop = FALSE])

# FUn (kNN in aggregated space)
nn <- RANN::nn2(data = mat_agg, query = mat_agg, k = PARAMS$k_neighbours + 1)
FUn_agg_raw <- rowMeans(nn$nn.dists[, -1, drop = FALSE])

# FSp (distance to centroid)
centroid <- colMeans(mat_agg, na.rm = TRUE)
FSp_agg_raw <- sqrt(rowSums(sweep(mat_agg, 2, centroid, FUN = "-")^2))

# Standardisation Pimiento ([0,1] x 4)
FUn_agg_std <- PARAMS$fuse_scale * scales::rescale(FUn_agg_raw, to = c(0, 1))
FSp_agg_std <- PARAMS$fuse_scale * scales::rescale(FSp_agg_raw, to = c(0, 1))

# Aggregated FUSE
FUSE_agg <- log(1 + FUn_agg_std * metrics_df$GE) +
  log(1 + FSp_agg_std * metrics_df$GE)


# ============================================================================
# 5. Disaggregated global FUSE ------ ---------------------------------------------------
# ============================================================================

if (PARAMS$disaggregated_combine == "sum") {
  FUSE_disagg <- metrics_df$FUSE_loco +
    metrics_df$FUSE_diet +
    metrics_df$FUSE_repro
} else if (PARAMS$disaggregated_combine == "max") {
  FUSE_disagg <- pmax(metrics_df$FUSE_loco,
                      metrics_df$FUSE_diet,
                      metrics_df$FUSE_repro)
} else {
  stop("disaggregated_combine doit être 'sum' ou 'max'.")
}


# ============================================================================
# 6. Tableau comparatif -------------------------------------------------------
# ============================================================================
comp_df <- tibble(
  species         = metrics_df$species,
  iucn            = metrics_df$iucn,
  FUSE_aggregated   = FUSE_agg,
  FUSE_disaggregated = FUSE_disagg
) %>%
  mutate(
    rank_aggregated   = rank(-FUSE_aggregated, ties.method = "min"),
    rank_disaggregated = rank(-FUSE_disaggregated, ties.method = "min"),
    iucn = factor(iucn, levels = PARAMS$iucn_levels)
  )

# top-25 par approche
top_agg    <- comp_df %>% slice_min(rank_aggregated,
                                    n = PARAMS$top_n) %>% pull(species)
top_disagg <- comp_df %>% slice_min(rank_disaggregated,
                                    n = PARAMS$top_n) %>% pull(species)

# Membership categorisation
comp_df <- comp_df %>%
  mutate(
    in_top_agg    = species %in% top_agg,
    in_top_disagg = species %in% top_disagg,
    overlap_class = case_when(
      in_top_agg &  in_top_disagg ~ "Both",
      in_top_agg & !in_top_disagg ~ "Aggregated only",
      !in_top_agg &  in_top_disagg ~ "Disaggregated only",
      TRUE ~ "Neither"
    ),
    overlap_class = factor(overlap_class,
                           levels = c("Both",
                                      "Aggregated only",
                                      "Disaggregated only",
                                      "Neither"))
  )

# Sauvegarde
saveRDS(list(comp_df = comp_df,
             FUn_agg_raw = FUn_agg_raw,
             FSp_agg_raw = FSp_agg_raw,
             var_explained = var_explained,
             params = PARAMS),
        file.path(PARAMS$data_out_dir, PARAMS$data_out_file))

# Spearman correlation entre les deux approches
rho_full <- cor(comp_df$FUSE_aggregated, comp_df$FUSE_disaggregated,
                method = "spearman", use = "pairwise.complete.obs")
n_overlap <- sum(comp_df$in_top_agg & comp_df$in_top_disagg)
n_agg_only    <- sum(comp_df$in_top_agg & !comp_df$in_top_disagg)
n_disagg_only <- sum(!comp_df$in_top_agg & comp_df$in_top_disagg)

message(sprintf("\nCorrélation Spearman FUSE_agg vs FUSE_disagg : %.3f",
                rho_full))
message(sprintf("Top-%d: %d shared, %d aggregated only, %d disaggregated only",
                PARAMS$top_n, n_overlap, n_agg_only, n_disagg_only))


# ============================================================================
# 7. (a) Scatter plot des rangs ----------------------------------------------
# ============================================================================
# Ranks: low = higher priority (rank 1 = top FUSE).
# On inverse l'axe pour que le coin haut-droit corresponde aux meilleurs.
n_total <- nrow(comp_df)

panel_a <- ggplot(comp_df,
                  aes(x = rank_aggregated, y = rank_disaggregated)) +
  geom_abline(intercept = 0, slope = 1,
              color = "grey60", linetype = "dashed", linewidth = 0.3) +
  # majority of species (low FUSE): scattered points
  geom_point(data = . %>% filter(overlap_class == "Neither"),
             color = "grey80", size = 0.3, alpha = 0.4) +
  # top-25 species: highlighted
  geom_point(data = . %>% filter(overlap_class != "Neither"),
             aes(color = iucn, shape = overlap_class),
             size = 1.4, alpha = 0.9, stroke = 0.4) +
  scale_color_manual(values = PARAMS$iucn_colors,
                     drop = FALSE,
                     name = "IUCN") +
  scale_shape_manual(values = c(Both = 16,
                                `Aggregated only` = 4,
                                `Disaggregated only` = 5,
                                Neither = 1),
                     drop = FALSE,
                     name = "Top-25 membership") +
  scale_x_reverse() +   # rang 1 à droite
  scale_y_reverse() +
  annotate("text",
           x = n_total * 0.02, y = n_total * 0.02,
           label = sprintf("Spearman ρ = %.3f\nN = %d",
                           rho_full, n_total),
           hjust = 1, vjust = 0,
           size = PARAMS$base_font_size / 3, color = "grey20") +
  labs(x = "Rank under FUSE_aggregated  (1 = top)",
       y = "Rank under FUSE_disaggregated  (1 = top)",
       title = "(a) Species rankings: aggregated vs. disaggregated FUSE") +
  theme_minimal(base_size = PARAMS$base_font_size) +
  theme(
    plot.title          = element_text(size = PARAMS$base_font_size + 1,
                                       face = "bold", hjust = 0,
                                       margin = margin(b = 3)),
    plot.title.position = "plot",
    panel.grid.minor    = element_blank(),
    panel.grid.major    = element_line(color = "grey92", linewidth = 0.25),
    axis.text           = element_text(size = PARAMS$base_font_size - 1.5,
                                       color = "grey20"),
    axis.title          = element_text(size = PARAMS$base_font_size - 0.5),
    axis.ticks          = element_line(color = "grey60", linewidth = 0.25),
    axis.ticks.length   = unit(1, "mm"),
    panel.background    = element_rect(fill = "white", color = NA),
    plot.background     = element_rect(fill = "white", color = NA),
    panel.border        = element_rect(color = "grey60", fill = NA,
                                       linewidth = 0.3),
    legend.position     = "right",
    legend.title        = element_text(size = PARAMS$base_font_size - 0.5),
    legend.text         = element_text(size = PARAMS$base_font_size - 1.5),
    legend.key.size     = unit(3, "mm")
  )


# ============================================================================
# 8. (b) Venn 2-ensembles top-25 ---------------------------------------------
# ============================================================================
# Built in pure ggplot2, no external dependency

cx_A <- -0.6;  cy_A <- 0
cx_D <-  0.6;  cy_D <- 0
r    <- 1.0

make_circle <- function(x0, y0, r, n = 200) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(x = x0 + r * cos(theta), y = y0 + r * sin(theta))
}
circle_A <- make_circle(cx_A, cy_A, r) %>% mutate(set = "Aggregated")
circle_D <- make_circle(cx_D, cy_D, r) %>% mutate(set = "Disaggregated")

panel_b <- ggplot() +
  geom_polygon(data = circle_A, aes(x, y, fill = set),
               alpha = 0.25, color = NA) +
  geom_polygon(data = circle_D, aes(x, y, fill = set),
               alpha = 0.25, color = NA) +
  geom_path(data = circle_A, aes(x, y), color = "#7B3294", linewidth = 0.5) +
  geom_path(data = circle_D, aes(x, y), color = "#008837", linewidth = 0.5) +
  # Comptes
  annotate("text", x = -1.0, y = 0,
           label = n_agg_only,
           size = PARAMS$base_font_size / 2,
           fontface = "bold", color = "grey15") +
  annotate("text", x = 1.0, y = 0,
           label = n_disagg_only,
           size = PARAMS$base_font_size / 2,
           fontface = "bold", color = "grey15") +
  annotate("text", x = 0, y = 0,
           label = n_overlap,
           size = PARAMS$base_font_size / 2,
           fontface = "bold", color = "grey15") +
  # Étiquettes de sets
  annotate("text", x = -1.0, y = -1.25,
           label = "Aggregated\nFUSE top-25",
           size = PARAMS$base_font_size / 2.5,
           fontface = "bold", color = "#7B3294") +
  annotate("text", x = 1.0, y = -1.25,
           label = "Disaggregated\nFUSE top-25",
           size = PARAMS$base_font_size / 2.5,
           fontface = "bold", color = "#008837") +
  scale_fill_manual(values = c(Aggregated = "#7B3294",
                               Disaggregated = "#008837"),
                    guide = "none") +
  coord_equal(xlim = c(-2.2, 2.2), ylim = c(-1.8, 1.5)) +
  labs(title = sprintf("(b) Top-%d overlap", PARAMS$top_n)) +
  theme_void(base_size = PARAMS$base_font_size) +
  theme(
    plot.title          = element_text(size = PARAMS$base_font_size + 1,
                                       face = "bold", hjust = 0,
                                       margin = margin(b = 3)),
    plot.title.position = "plot",
    plot.background     = element_rect(fill = "white", color = NA),
    plot.margin         = margin(2, 2, 2, 2)
  )


# ============================================================================
# 9. (c) Top-25 species lists side by side -------- ------------------------------------------
# ============================================================================
# Top-25 species per approach, flagged if they appear in the other approach
# dans les DEUX listes (point plein) ou seulement dans une seule (point vide).

build_top_table <- function(comp_df, rank_col, n) {
  comp_df %>%
    arrange(.data[[rank_col]]) %>%
    slice_head(n = n) %>%
    mutate(rank_position = row_number(),
           species_short = sub("^([A-Za-z])[a-z]+_", "\\1. ", species),
           species_label = factor(species_short,
                                  levels = rev(species_short)),
           shared = species %in% intersect(top_agg, top_disagg))
}

top_agg_df    <- build_top_table(comp_df, "rank_aggregated",    PARAMS$top_n)
top_disagg_df <- build_top_table(comp_df, "rank_disaggregated", PARAMS$top_n)

make_top_list_plot <- function(df, title, set_color) {
  ggplot(df,
         aes(x = 1, y = species_label, fill = iucn)) +
    geom_point(aes(shape = shared), size = 2.5,
               stroke = 0.4, color = "grey20") +
    scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 23),
                       labels = c(`TRUE` = "Shared with other approach",
                                  `FALSE` = "Unique to this approach"),
                       name = NULL,
                       drop = FALSE) +
    scale_fill_manual(values = PARAMS$iucn_colors,
                      drop = FALSE,
                      name = "IUCN") +
    scale_x_continuous(limits = c(0.7, 1.3), breaks = NULL) +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = PARAMS$base_font_size) +
    theme(
      plot.title          = element_text(size = PARAMS$base_font_size + 1,
                                         face = "bold", color = set_color,
                                         hjust = 0,
                                         margin = margin(b = 2)),
      plot.title.position = "plot",
      panel.grid.major.y  = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.minor    = element_blank(),
      panel.grid.major.x  = element_blank(),
      axis.text.y         = element_text(size = PARAMS$base_font_size - 2,
                                         face = "italic", color = "grey20"),
      axis.text.x         = element_blank(),
      axis.ticks          = element_blank(),
      panel.background    = element_rect(fill = "white", color = NA),
      plot.background     = element_rect(fill = "white", color = NA),
      panel.border        = element_blank(),
      legend.position     = "none",
      plot.margin         = margin(2, 4, 2, 2)
    )
}

panel_c_left  <- make_top_list_plot(top_agg_df,
                                    "(c) Top-25 FUSE_aggregated",
                                    "#7B3294")
panel_c_right <- make_top_list_plot(top_disagg_df,
                                    "    Top-25 FUSE_disaggregated",
                                    "#008837")


# Global legend for panel (c)
make_iucn_legend <- function() {
  df <- tibble(x = 1, y = 1,
               iucn = factor(PARAMS$iucn_levels, levels = PARAMS$iucn_levels))
  ggplot(df, aes(x, y, fill = iucn)) +
    geom_point(shape = 21, size = 3) +
    scale_fill_manual(values = PARAMS$iucn_colors, drop = FALSE,
                      name = "IUCN") +
    coord_cartesian(xlim = c(2, 3), ylim = c(2, 3)) +
    theme_void(base_size = PARAMS$base_font_size) +
    theme(legend.position = "bottom",
          legend.title    = element_text(size = PARAMS$base_font_size,
                                         face = "bold"),
          legend.text     = element_text(size = PARAMS$base_font_size - 1),
          legend.key.size = unit(3, "mm"))
}

make_shape_legend <- function() {
  df <- tibble(x = 1, y = 1, shared = factor(c(TRUE, FALSE)))
  ggplot(df, aes(x, y, shape = shared)) +
    geom_point(size = 2.5, fill = "grey60", color = "grey20", stroke = 0.4) +
    scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 23),
                       labels = c(`TRUE` = "Shared",
                                  `FALSE` = "Unique to one approach"),
                       name = NULL) +
    coord_cartesian(xlim = c(2, 3), ylim = c(2, 3)) +
    theme_void(base_size = PARAMS$base_font_size) +
    theme(legend.position = "bottom",
          legend.text     = element_text(size = PARAMS$base_font_size - 1),
          legend.key.size = unit(3, "mm"))
}

iucn_legend  <- make_iucn_legend()
shape_legend <- make_shape_legend()


# ============================================================================
# 10. Assemblage final --------------------------------------------------------
# ============================================================================
row_ab <- patchwork::wrap_plots(
  panel_a, panel_b,
  nrow = 1,
  widths = c(0.6, 0.4)
)

row_c <- patchwork::wrap_plots(
  panel_c_left, panel_c_right,
  nrow = 1
)

row_legends <- patchwork::wrap_plots(
  iucn_legend, shape_legend,
  nrow = 1
)

final_fig <- (row_ab / row_c / row_legends) +
  patchwork::plot_layout(heights = c(1.4, 2.5, 0.18)) 
  # patchwork::plot_annotation(
  #   caption = sprintf(
  #     "Aggregated FUSE = single multidimensional space built from all traits combined (Gower distance + PCoA, %d axes). \nDisaggregated FUSE = sum of dimension-specific FUSE scores (locomotion + diet + reproduction). \nThe disaggregation reframes the priority list and identifies species that aggregated approaches systematically overlook.",
  #     PARAMS$n_axes_aggregated
  #   ),
  #   theme = theme(
  #     plot.caption    = element_text(size = PARAMS$base_font_size - 1.5,
  #                                    color = "grey25", hjust = 0,
  #                                    margin = margin(t = 4)),
  #     plot.background = element_rect(fill = "white", color = NA)
  #   )
  # )


# ============================================================================
# 11. Export ------------------------------------------------------------------
# ============================================================================
out_pdf <- file.path(PARAMS$out_dir, paste0(PARAMS$fig_basename, ".pdf"))

ggsave(out_pdf, final_fig,
       width = PARAMS$width_mm, height = PARAMS$height_mm,
       units = "mm", device = cairo_pdf)


# ============================================================================
# 12. Summary tables for the manuscript -------- --------------------------------
# ============================================================================
message("\n=== Espèces top-",
        PARAMS$top_n, " AGGREGATED only (absent from disaggregated) ===")
print(comp_df %>% filter(in_top_agg & !in_top_disagg) %>%
        select(species, iucn, FUSE_aggregated, FUSE_disaggregated,
               rank_aggregated, rank_disaggregated))

message("\n=== Espèces top-",
        PARAMS$top_n, " DISAGGREGATED only (absent from aggregated) ===")
print(comp_df %>% filter(!in_top_agg & in_top_disagg) %>%
        select(species, iucn, FUSE_aggregated, FUSE_disaggregated,
               rank_aggregated, rank_disaggregated))

message("\n=== Espèces communes aux deux top-", PARAMS$top_n, " ===")
print(comp_df %>% filter(in_top_agg & in_top_disagg) %>%
        select(species, iucn, FUSE_aggregated, FUSE_disaggregated))


# ============================================================================
# 13. Reproducibility -------- --------------------------------------------------------
# ============================================================================
sessionInfo()