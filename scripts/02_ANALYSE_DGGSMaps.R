#-------------------------------------------------------------------------------
# Figure S4. DGGS maps of community FRic loss under IUCN AT.
#
# Analogous to Pimiento et al. (2020), Fig. 4: for each functional dimension
# (locomotion, diet, reproduction), a global map of the percentage
# change in community FRic under the IUCN AT scenario (deterministic
# extinction of all VU + EN + CR species), computed per DGGS cell.
#
# Differences from manuscript Fig. 3:
#   * single row (FRic only, not FUn)
#   * SHARED colour scale across the three dimensions, allowing
#     direct visual comparison of loss magnitudes across
#     locomotion / diet / reproduction (unchanged)
#   * supplement format (180 x 80 mm)
#
# Author  : A. Toussaint
#
# Prerequisites:
#   - extinction_scenarios_cellwise_sf.rds : produced by script C, with
#     columns FRic_change_AT__{locomotion,diet,reproduction} per cell
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c("ggplot2", "dplyr", "tidyr", "sf",
                   "patchwork", "scales",
                   "rnaturalearth", "rnaturalearthdata")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))


# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(
  # Input
  cellwise_sf_path = "data/processed/extinction_scenarios_cellwise_sf.rds",
  
  # Output
  out_dir       = "figures",
  fig_basename  = "FigS4_FRic_loss_maps_IUCNAT",
  width_mm      = 180,
  height_mm     = 80,
  dpi           = 600,
  
  # Apparence
  base_font_size = 8,
  
  # Palette divergente : rouge profond pour pertes intenses,
  # gris clair pour pas de changement, bleu profond pour gains.
  pal_low       = "#67001f",
  pal_mid       = "#f7f7f7",
  pal_high      = "#053061",
  
  # Common limits (computed over the union of the three dimensions)
  # ; leave NULL for automatic quantile-based limits
  fric_limits   = NULL,
  clip_quantile = 0.99,
  
  # Projection
  crs_target    = "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs",
  
  # Dimension order and labels
  dim_order     = c("locomotion", "diet", "reproduction"),
  dim_labels    = c(locomotion   = "Locomotion",
                    diet         = "Diet",
                    reproduction = "Reproduction"),
  dim_colors    = c(locomotion   = "#2E7D32",
                    diet         = "#C62828",
                    reproduction = "#1565C0")
)
dir.create(PARAMS$out_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================================
# 2. Chargement ---------------------------------------------------------------
# ============================================================================
cell_sf <- readRDS(PARAMS$cellwise_sf_path)

# Check IUCN AT columns are present for the three dimensions
needed_cols <- paste0("FRic_change_AT__", PARAMS$dim_order)
missing_cols <- setdiff(needed_cols, names(cell_sf))
if (length(missing_cols) > 0) {
  stop("Colonnes manquantes dans cell_sf : ",
       paste(missing_cols, collapse = ", "))
}

# Wrap dateline to avoid stretched polygons in projection
if (is.na(sf::st_crs(cell_sf))) cell_sf <- sf::st_set_crs(cell_sf, 4326)
cell_sf <- sf::st_wrap_dateline(cell_sf,
                                options = c("WRAPDATELINE=YES",
                                            "DATELINEOFFSET=180"))
cell_sf_proj <- sf::st_transform(cell_sf, PARAMS$crs_target)

world <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf") %>%
  sf::st_wrap_dateline(options = c("WRAPDATELINE=YES",
                                   "DATELINEOFFSET=180")) %>%
  sf::st_transform(PARAMS$crs_target)


# ============================================================================
# 3. SHARED palette limits (key for Fig. S4 message) -------- -----------------
# ============================================================================
# Pimiento Fig. 4 uses a single colour scale to allow direct
# visual comparison between regions; here extended across dimensions.

if (is.null(PARAMS$fric_limits)) {
  fric_vals <- unlist(lapply(PARAMS$dim_order, function(d) {
    cell_sf_proj[[paste0("FRic_change_AT__", d)]]
  }))
  fric_vals <- fric_vals[!is.na(fric_vals)]
  ext <- quantile(abs(fric_vals), PARAMS$clip_quantile, na.rm = TRUE)
  fric_lims <- c(-ext, ext)
} else {
  fric_lims <- PARAMS$fric_limits
}

message(sprintf("Limites de palette FRic_change (symétriques) : [%.1f, %.1f]",
                fric_lims[1], fric_lims[2]))


# ============================================================================
# 4. Constructeur de carte ---------------------------------------------------
# ============================================================================
theme_map <- function(base_size = PARAMS$base_font_size) {
  theme_void(base_size = base_size) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin      = margin(1, 1, 1, 1),
      legend.position  = "none",
      plot.title       = element_text(size = base_size + 1,
                                      face = "bold", hjust = 0.5,
                                      margin = margin(b = 1))
    )
}

make_fric_map <- function(sf_data, value_col, lims, dim_name, dim_color) {
  ggplot() +
    geom_sf(data = world, fill = "grey92",
            color = "grey75", linewidth = 0.1) +
    geom_sf(data = sf_data,
            aes(fill = .data[[value_col]]),
            color = NA) +
    scale_fill_gradient2(
      low      = PARAMS$pal_low,
      mid      = PARAMS$pal_mid,
      high     = PARAMS$pal_high,
      midpoint = 0,
      limits   = lims,
      oob      = scales::squish,
      na.value = "grey95"
    ) +
    coord_sf(crs = PARAMS$crs_target, expand = FALSE) +
    labs(title = dim_name) +
    theme_map() +
    theme(plot.title = element_text(color = dim_color))
}


# ============================================================================
# 5. Construction des trois cartes -------------------------------------------
# ============================================================================
maps_FRic <- lapply(PARAMS$dim_order, function(d) {
  col <- paste0("FRic_change_AT__", d)
  make_fric_map(cell_sf_proj, col, fric_lims,
                dim_name  = PARAMS$dim_labels[d],
                dim_color = PARAMS$dim_colors[d])
})


# ============================================================================
# 6. Shared legend (bottom) -------- --------------------------------------
# ============================================================================
make_legend_only <- function(lims, label) {
  df <- data.frame(x = c(0, 1), y = c(0, 1), v = lims)
  ggplot(df, aes(x, y, fill = v)) +
    geom_tile() +
    scale_fill_gradient2(
      low      = PARAMS$pal_low,
      mid      = PARAMS$pal_mid,
      high     = PARAMS$pal_high,
      midpoint = 0,
      limits   = lims,
      oob      = scales::squish,
      name     = label,
      guide    = guide_colorbar(
        title.position = "top",
        title.hjust    = 0.5,
        barwidth       = unit(60, "mm"),
        barheight      = unit(2.5, "mm"),
        ticks.colour   = "grey30",
        frame.colour   = "grey30"
      )
    ) +
    coord_cartesian(xlim = c(2, 3), ylim = c(2, 3)) +
    theme_void(base_size = PARAMS$base_font_size) +
    theme(
      legend.position    = "bottom",
      legend.title       = element_text(size = PARAMS$base_font_size,
                                        face = "bold"),
      legend.text        = element_text(size = PARAMS$base_font_size - 1.5),
      plot.background    = element_rect(fill = "white", color = NA),
      plot.margin        = margin(0, 0, 0, 0)
    )
}

legend_FRic <- make_legend_only(fric_lims, "FRic change (%) under IUCN AT")


# ============================================================================
# 7. Assemblage final --------------------------------------------------------
# ============================================================================
maps_row <- patchwork::wrap_plots(
  maps_FRic[[1]], maps_FRic[[2]], maps_FRic[[3]],
  nrow = 1
)

final_fig <- (maps_row / legend_FRic) +
  patchwork::plot_layout(heights = c(1, 0.18)) 
  # patchwork::plot_annotation(
  #   caption = "Cell-level percent change in functional richness under the IUCN All-Threatened scenario, computed independently for each ecological dimension. \nNegative values indicate functional contraction. \nCommon colour scale across panels facilitates direct comparison of the amplitude of losses between dimensions.",
  #   theme = theme(
  #     plot.caption = element_text(size = PARAMS$base_font_size - 1.5,
  #                                 color = "grey25",
  #                                 hjust = 0,
  #                                 margin = margin(t = 3)),
  #     plot.background = element_rect(fill = "white", color = NA)
  #   )
  # )


# ============================================================================
# 8. Export -------------------------------------------------------------------
# ============================================================================
out_pdf <- file.path(PARAMS$out_dir, paste0(PARAMS$fig_basename, ".pdf"))

ggsave(out_pdf, final_fig,
       width = PARAMS$width_mm, height = PARAMS$height_mm,
       units = "mm", device = cairo_pdf)

message("Figure S4 written to:\n  - ", out_pdf,
        "\n  - ", out_png, "\n  - ", out_svg)


# ============================================================================
# 9. Summary statistics for the manuscript -------- ------------------------------
# ============================================================================
fric_summary <- bind_rows(lapply(PARAMS$dim_order, function(d) {
  v <- cell_sf_proj[[paste0("FRic_change_AT__", d)]]
  v <- v[!is.na(v)]
  tibble(
    dimension     = PARAMS$dim_labels[d],
    n_cells       = length(v),
    mean_change   = mean(v),
    median_change = median(v),
    q05           = quantile(v, 0.05),
    q95           = quantile(v, 0.95),
    pct_below_minus10 = 100 * mean(v < -10),
    pct_below_minus25 = 100 * mean(v < -25)
  )
}))

message("\n=== Statistiques par dimension (FRic change sous IUCN AT) ===")
print(fric_summary)


# ============================================================================
# 10. Reproducibility -------- --------------------------------------------------------
# ============================================================================
sessionInfo()