#-------------------------------------------------------------------------------
# Figure S5. Venn diagrams of top-5% FUS and FUSE overlaps across dimensions.
#
# Retained as supplement: original version of the paper (before the
# transition to disaggregated FUSE in main figures). Shows that:
#  - FUS (no extinction weighting): low overlap among the 3
#    dimensions, reflecting ecological independence of trait spaces.
#  - FUSE (extinction-weighted): higher overlap, consequence
#    of the shared GE factor across spaces, though asymmetry persists.
#
# Layout: 2 rows (FUS, FUSE) x 1 column (3 sets: Locomotion, Diet,
#          Reproduction). Each panel is a 3-set Venn diagram built
#          in pure ggplot (no eulerr/VennDiagram dependency).
#
# Format: 180 mm x 130 mm.
#
# Author  : A. Toussaint
#
# Prerequisites:
#   - metrics_df : tibble produced by script A1, with columns:
#       species, iucn,
#       FUS_loco, FUS_diet, FUS_repro,
#       FUSE_loco, FUSE_diet, FUSE_repro
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c(
  "ggplot2", "dplyr", "tidyr", "tibble", "purrr", "scales",
  "patchwork"
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
  
  # Top-X% pour les Venn (5% comme dans la version originale du papier)
  top_pct          = 0.05,
  
  # Suffixes des colonnes par dimension dans metrics_df
  fus_cols  = c(locomotion = "FUS_loco",
                diet       = "FUS_diet",
                reproduction = "FUS_repro"),
  fuse_cols = c(locomotion = "FUSE_loco",
                diet       = "FUSE_diet",
                reproduction = "FUSE_repro"),
  
  # Apparence
  base_font_size = 8,
  dim_order      = c("locomotion", "diet", "reproduction"),
  dim_labels     = c(locomotion   = "Locomotion",
                     diet         = "Diet",
                     reproduction = "Reproduction"),
  dim_colors     = c(locomotion   = "#2E7D32",
                     diet         = "#C62828",
                     reproduction = "#1565C0"),
  
  # Output
  out_dir       = "figures",
  fig_basename  = "FigS5_Venn_top5pct_FUS_FUSE",
  width_mm      = 180,
  height_mm     = 130,
  dpi           = 600
)
dir.create(PARAMS$out_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================================
# 2. USER INPUTS ------------------------------------------------------
# ============================================================================
metrics_df <- readRDS("data/processed/species_metrics_FUn_FSp_FUSE.rds")

stopifnot(
  exists("metrics_df"),
  all(c("species", PARAMS$fus_cols, PARAMS$fuse_cols) %in% names(metrics_df))
)


# ============================================================================
# 3. Select top-5% per dimension and index -------- --------------------------
# ============================================================================
n_total <- nrow(metrics_df)
n_top   <- ceiling(n_total * PARAMS$top_pct)
message(sprintf("Top-%.0f%%: %d species out of %d.",
                100 * PARAMS$top_pct, n_top, n_total))

#' For an index (FUS or FUSE) and a dimension, returns the top-X% species
get_top_species <- function(metrics_df, score_col, n_top) {
  metrics_df %>%
    arrange(desc(.data[[score_col]])) %>%
    slice_head(n = n_top) %>%
    pull(species)
}

# Top-5% sets
top_FUS <- list(
  Locomotion   = get_top_species(metrics_df, PARAMS$fus_cols[["locomotion"]],   n_top),
  Diet         = get_top_species(metrics_df, PARAMS$fus_cols[["diet"]],         n_top),
  Reproduction = get_top_species(metrics_df, PARAMS$fus_cols[["reproduction"]], n_top)
)
top_FUSE <- list(
  Locomotion   = get_top_species(metrics_df, PARAMS$fuse_cols[["locomotion"]],   n_top),
  Diet         = get_top_species(metrics_df, PARAMS$fuse_cols[["diet"]],         n_top),
  Reproduction = get_top_species(metrics_df, PARAMS$fuse_cols[["reproduction"]], n_top)
)


# ============================================================================
# 4. Fonction Venn 3-ensembles native ggplot ---------------------------------
# ============================================================================

#' Build a 3-set Venn diagram from three lists.
#' @param sets       named list of 3 species identifier vectors
#'                   (noms : "Locomotion", "Diet", "Reproduction")
#' @param title      titre du panneau
#' @param dim_colors palette des trois dimensions
make_venn3 <- function(sets, title, dim_colors) {
  
  L <- sets[["Locomotion"]]
  D <- sets[["Diet"]]
  R <- sets[["Reproduction"]]
  
  # Counts per region
  union_all <- unique(c(L, D, R))
  in_L <- union_all %in% L
  in_D <- union_all %in% D
  in_R <- union_all %in% R
  
  n_only_L  <- sum( in_L & !in_D & !in_R)
  n_only_D  <- sum(!in_L &  in_D & !in_R)
  n_only_R  <- sum(!in_L & !in_D &  in_R)
  n_LD      <- sum( in_L &  in_D & !in_R)
  n_LR      <- sum( in_L & !in_D &  in_R)
  n_DR      <- sum(!in_L &  in_D &  in_R)
  n_LDR     <- sum( in_L &  in_D &  in_R)
  
  # Circle centres (equilateral triangle)
  r <- 1.4
  cx_L <- -0.85; cy_L <- -0.5
  cx_D <-  0.85; cy_D <- -0.5
  cx_R <-  0.0;  cy_R <-  0.85
  
  circles <- tibble(
    x0 = c(cx_L, cx_D, cx_R),
    y0 = c(cy_L, cy_D, cy_R),
    r  = r,
    set = c("Locomotion", "Diet", "Reproduction"),
    fill_color = c(dim_colors[["locomotion"]],
                   dim_colors[["diet"]],
                   dim_colors[["reproduction"]])
  )
  
  make_circle <- function(x0, y0, r, n = 200) {
    theta <- seq(0, 2 * pi, length.out = n)
    data.frame(x = x0 + r * cos(theta), y = y0 + r * sin(theta))
  }
  circle_points <- bind_rows(lapply(seq_len(nrow(circles)), function(i) {
    make_circle(circles$x0[i], circles$y0[i], circles$r[i]) %>%
      mutate(set = circles$set[i],
             fill_color = circles$fill_color[i])
  }))
  
  # Étiquettes de comptes
  labels_df <- tibble(
    x = c(-1.6, 1.6, 0.0, 0.0, -0.85, 0.85, 0.0),
    y = c(-0.85, -0.85, 1.45, -0.95, 0.45, 0.45, 0.0),
    count = c(n_only_L, n_only_D, n_only_R, n_LD, n_LR, n_DR, n_LDR)
  )
  
  # Étiquettes de set
  set_labels <- tibble(
    x = c(cx_L - 0.7, cx_D + 0.7, cx_R),
    y = c(cy_L - r - 0.15, cy_D - r - 0.15, cy_R + r + 0.15),
    set = c("Locomotion", "Diet", "Reproduction"),
    color = c(dim_colors[["locomotion"]],
              dim_colors[["diet"]],
              dim_colors[["reproduction"]])
  )
  
  # Total (pour le sous-titre)
  n_union <- length(union_all)
  
  ggplot() +
    geom_polygon(data = circle_points,
                 aes(x = x, y = y, group = set, fill = set),
                 alpha = 0.25, color = NA) +
    geom_path(data = circle_points,
              aes(x = x, y = y, group = set, color = set),
              linewidth = 0.5) +
    geom_text(data = labels_df,
              aes(x = x, y = y, label = count),
              size = PARAMS$base_font_size / 2,
              fontface = "bold", color = "grey15") +
    geom_text(data = set_labels,
              aes(x = x, y = y, label = set, color = set),
              size = PARAMS$base_font_size / 2.5,
              fontface = "bold", show.legend = FALSE) +
    scale_fill_manual(values = setNames(circles$fill_color, circles$set),
                      guide = "none") +
    scale_color_manual(values = setNames(circles$fill_color, circles$set),
                       guide = "none") +
    coord_equal(xlim = c(-2.6, 2.6), ylim = c(-2.4, 2.6)) +
    labs(title = title,
         subtitle = sprintf("union = %d species", n_union)) +
    theme_void(base_size = PARAMS$base_font_size) +
    theme(
      plot.title          = element_text(size = PARAMS$base_font_size + 1,
                                         face = "bold", hjust = 0,
                                         margin = margin(b = 1)),
      plot.title.position = "plot",
      plot.subtitle       = element_text(size = PARAMS$base_font_size - 0.5,
                                         color = "grey40", hjust = 0,
                                         margin = margin(b = 4)),
      plot.background     = element_rect(fill = "white", color = NA),
      plot.margin         = margin(2, 4, 2, 2)
    )
}


# ============================================================================
# 5. Construction des deux panneaux ------------------------------------------
# ============================================================================
panel_FUS <- make_venn3(
  sets       = top_FUS,
  title      = sprintf("(a) FUS — top %.0f%% per dimension",
                       100 * PARAMS$top_pct),
  dim_colors = PARAMS$dim_colors
)

panel_FUSE <- make_venn3(
  sets       = top_FUSE,
  title      = sprintf("(b) FUSE — top %.0f%% per dimension",
                       100 * PARAMS$top_pct),
  dim_colors = PARAMS$dim_colors
)


# ============================================================================
# 6. Assemblage final --------------------------------------------------------
# ============================================================================
final_fig <- patchwork::wrap_plots(
  panel_FUS, panel_FUSE,
  nrow = 1
) 
  # patchwork::plot_annotation(
  #   caption = sprintf(
  #     "Overlap among the top %.0f%% of bird species ranked by Functional Uniqueness–Specialisation (FUS, panel a) and threat-weighted FUSE (panel b), \nacross three independent ecological trait spaces.\nFUS captures intrinsic functional distinctiveness; FUSE additionally weights by IUCN-derived extinction probability. \nLimited overlap under FUS demonstrates that functional importance is dimension-dependent. \nIncreased overlap under FUSE reflects the integrating effect of extinction risk but does not eliminate trait-space specificity.",
  #     100 * PARAMS$top_pct
  #   ),
  #   theme = theme(
  #     plot.caption    = element_text(size = PARAMS$base_font_size - 1.5,
  #                                    color = "grey25",
  #                                    hjust = 0,
  #                                    margin = margin(t = 4)),
  #     plot.background = element_rect(fill = "white", color = NA)
  #   )
  # )


# ============================================================================
# 7. Export -------------------------------------------------------------------
# ============================================================================
out_pdf <- file.path(PARAMS$out_dir, paste0(PARAMS$fig_basename, ".pdf"))

ggsave(out_pdf, final_fig,
       width = PARAMS$width_mm, height = PARAMS$height_mm,
       units = "mm", device = cairo_pdf)
message("Figure S5 written to:\n  - ", out_pdf,
        "\n  - ", out_png, "\n  - ", out_svg)


# ============================================================================
# 8. Summary table for the manuscript -------- -----------------------------------
# ============================================================================
summary_table <- function(sets, label) {
  L <- sets[["Locomotion"]]
  D <- sets[["Diet"]]
  R <- sets[["Reproduction"]]
  union_all <- unique(c(L, D, R))
  in_L <- union_all %in% L
  in_D <- union_all %in% D
  in_R <- union_all %in% R
  
  tibble(
    index            = label,
    n_total_per_dim  = length(L),
    n_union          = length(union_all),
    n_loco_only      = sum( in_L & !in_D & !in_R),
    n_diet_only      = sum(!in_L &  in_D & !in_R),
    n_repro_only     = sum(!in_L & !in_D &  in_R),
    n_loco_diet      = sum( in_L &  in_D & !in_R),
    n_loco_repro     = sum( in_L & !in_D &  in_R),
    n_diet_repro     = sum(!in_L &  in_D &  in_R),
    n_all_three      = sum( in_L &  in_D &  in_R),
    pct_all_three    = 100 * sum(in_L & in_D & in_R) / length(union_all)
  )
}

synth <- bind_rows(
  summary_table(top_FUS,  "FUS"),
  summary_table(top_FUSE, "FUSE")
)

message("\n=== Chevauchement top-",
        100 * PARAMS$top_pct, "% par dimension ===")
print(synth)


# ============================================================================
# 9. Reproducibility -------- ---------------------------------------------------------
# ============================================================================
sessionInfo()