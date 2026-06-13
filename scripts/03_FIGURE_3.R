#-------------------------------------------------------------------------------
# Figure 3. Cell-level functional change under extinction scenarios.
#
# 2-row (FRic change, FUn change) x 3-column (Locomotion, Diet,
# Reproduction) grid of global DGGS maps. The default scenario shown
# is IUCN AT (loss if all VU+EN+CR go extinct); the script
# allows switching to IUCN 100 via PARAMS$scenario_to_plot.
#
# Author  : A. Toussaint
#
# Prerequisites:
#   - extinction_scenarios_cellwise_sf.rds : sf object produced by
#     script C (cellwise extinction scenarios), with DGGS geometry
#     and columns:
#       FRic_change_AT__locomotion, FRic_change_AT__diet,
#       FRic_change_AT__reproduction, idem pour _100 et pour FUn_change_*
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
  fig_basename  = "Fig3_cellwise_maps",
  width_mm      = 180,
  height_mm     = 110,
  dpi           = 600,
  
  # Scenario to display: "AT" or "100"
  scenario_to_plot = "AT",
  
  # Apparence
  base_font_size = 8,
  
  # Diverging palette — IMPORTANT: centred on 0, light for neutral,
  # deep red for extreme negative values (large FRic loss),
  # deep blue for extreme positive values.
  # For FRic: mostly negative or zero values
  # For FUn : mostly positive or zero values
  pal_low       = "#08221CFF",   # rouge profond (perte intense)
  pal_mid       = "#72972FFF",   # gris clair (pas de changement)
  pal_high      = "#D4D5BFFF",   # bleu profond (gain intense)
  
  
  # different; limits computed separately)
  # If NULL, derived from data (with clipping at a robust quantile
  # to prevent one extreme cell from crushing the palette).
  fric_limits   = NULL,
  fun_limits    = NULL,
  clip_quantile = 0.99,        # écrêtage pour limites automatiques
  
  # Projection pour les cartes
  # Robinson (compromis classique pour cartes globales)
  crs_target    = "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs",
  
  # Dimension order and labels
  dim_order     = c("locomotion", "diet", "reproduction"),
  dim_labels    = c(locomotion   = "Locomotion",
                    diet         = "Diet",
                    reproduction = "Reproduction"),
  
  # Étiquettes des couleurs des titres de colonne
  dim_colors    = c(locomotion   = "#2E7D32",
                    diet         = "#C62828",
                    reproduction = "#1565C0")
)
dir.create(PARAMS$out_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================================
# 2. Load and prepare data -------- -----------------------------------------------
# ============================================================================
cell_sf <- readRDS(PARAMS$cellwise_sf_path)

# Check expected columns are present
needed_cols <- as.vector(outer(
  c("FRic_change", "FUn_change"),
  PARAMS$dim_order,
  function(m, d) paste0(m, "_", PARAMS$scenario_to_plot, "__", d)
))
missing_cols <- setdiff(needed_cols, names(cell_sf))
if (length(missing_cols) > 0) {
  stop("Colonnes manquantes dans cell_sf : ",
       paste(missing_cols, collapse = ", "))
}

# Reprojection
# IMPORTANT: before reprojecting, clip DGGS polygons that
# straddle the antimeridian (180°/-180°). Without this step, those hexagons
# turn into horizontal bands spanning the whole globe, which
# masks nearly the entire map. st_wrap_dateline() splits those
# polygons along the date meridian.
# Prerequisites: cell_sf must be in lon/lat (EPSG:4326). CRS is enforced
# if missing, which is the default when coming from dgearthgrid().
if (is.na(sf::st_crs(cell_sf))) {
  cell_sf <- sf::st_set_crs(cell_sf, 4326)
}
cell_sf_wrapped <- sf::st_wrap_dateline(
  cell_sf,
  options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")
)
cell_sf_proj <- sf::st_transform(cell_sf_wrapped, PARAMS$crs_target)

# Basemap: world countries (Natural Earth) - simplified to thin lines
world <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf") %>%
  sf::st_wrap_dateline(options = c("WRAPDATELINE=YES",
                                   "DATELINEOFFSET=180")) %>%
  sf::st_transform(PARAMS$crs_target)


# ============================================================================
# 3. Determine palette limits per row -------- --------------------------
# ============================================================================
get_symmetric_limits <- function(values, q = 0.99) {
  v <- values[!is.na(values)]
  if (length(v) == 0) return(c(-1, 1))
  ext <- quantile(abs(v), q, na.rm = TRUE)
  c(-ext, ext)
}

if (is.null(PARAMS$fric_limits)) {
  fric_vals <- unlist(lapply(PARAMS$dim_order, function(d) {
    cell_sf[[paste0("FRic_change_", PARAMS$scenario_to_plot, "__", d)]]
  }))
  fric_lims <- get_symmetric_limits(fric_vals, q = PARAMS$clip_quantile)
  fric_lims[2] <- 0
} else fric_lims <- PARAMS$fric_limits

if (is.null(PARAMS$fun_limits)) {
  fun_vals <- unlist(lapply(PARAMS$dim_order, function(d) {
    cell_sf[[paste0("FUn_change_", PARAMS$scenario_to_plot, "__", d)]]
  }))
  fun_lims <- get_symmetric_limits(fun_vals, q = PARAMS$clip_quantile)
  fun_lims[1] <- 0
  fun_lims[2] <- 10
} else fun_lims <- PARAMS$fun_limits

message("FRic limits (symmetric): [", round(fric_lims[1], 2),
        ", ", round(fric_lims[2], 2), "]")
message("FUn  limits (symmetric): [", round(fun_lims[1], 2),
        ", ", round(fun_lims[2], 2), "]")


# ============================================================================
# 4. Common theme for maps -------- --------------------------------------------
# ============================================================================
theme_map <- function(base_size = PARAMS$base_font_size) {
  theme_void(base_size = base_size) +
    theme(
      plot.background      = element_rect(fill = "white", color = NA),
      panel.background     = element_rect(fill = "white", color = NA),
      plot.margin          = margin(1, 1, 1, 1),
      legend.position      = "none",  # géré globalement
      plot.title           = element_blank(),
      plot.subtitle        = element_blank()
    )
}


# ============================================================================
# 5. Constructeur de carte ---------------------------------------------------
# ============================================================================
make_map_FRic <- function(sf_data, value_col, lims, world_bg) {
  
  # To prevent extreme values from crushing the palette,
  # clip to defined limits (oob = squish)
  ggplot() +
    geom_sf(data = world_bg, fill = "grey92",
            color = "grey75", linewidth = 0.1) +
    geom_sf(data = sf_data,
            aes(fill = .data[[value_col]]),
            color = NA) +
    scale_fill_gradient2(
      low      = PARAMS$pal_low,
      mid      = PARAMS$pal_mid,
      high     = PARAMS$pal_high,
      midpoint = (0+fric_lims[1])/2,
      limits   = lims,
      oob      = scales::squish,
      na.value = "grey95"
    ) +
    coord_sf(crs = PARAMS$crs_target, expand = FALSE) +
    theme_map()
}

make_map_FUn <- function(sf_data, value_col, lims, world_bg) {
  
  # To prevent extreme values from crushing the palette,
  # clip to defined limits (oob = squish)
  ggplot() +
    geom_sf(data = world_bg, fill = "grey92",
            color = "grey75", linewidth = 0.1) +
    geom_sf(data = sf_data,
            aes(fill = .data[[value_col]]),
            color = NA) +
    scale_fill_gradient2(
      low      = PARAMS$pal_high,
      mid      = PARAMS$pal_mid,
      high     = PARAMS$pal_low,
      midpoint = (0+fun_lims[2])/2,
      limits   = lims,
      oob      = scales::squish,
      na.value = "grey95"
    ) +
    coord_sf(crs = PARAMS$crs_target, expand = FALSE) +
    theme_map()
}


# ============================================================================
# 6. Build the 6 panels + legends -------- ----------------------------------
# ============================================================================

# One reference map per row to extract the row-specific legend
maps_FRic <- lapply(PARAMS$dim_order, function(d) {
  col <- paste0("FRic_change_", PARAMS$scenario_to_plot, "__", d)
  make_map_FRic(cell_sf_proj, col, fric_lims, world)
})
names(maps_FRic) <- PARAMS$dim_order

maps_FUn <- lapply(PARAMS$dim_order, function(d) {
  col <- paste0("FUn_change_", PARAMS$scenario_to_plot, "__", d)
  make_map_FUn(cell_sf_proj, col, fun_lims, world)
})
names(maps_FUn) <- PARAMS$dim_order


# Build legend plots (one per row) — see make_legend_only below

#' Build a ggplot whose sole purpose is to display a legend
#' horizontal colourbar for a given value range. More robust than
#' l'extraction par cowplot, qui peut renvoyer plusieurs guide-boxes selon la
#' version de ggplot2 et casser le layout patchwork.
make_legend_only <- function(lims, label) {
  # dummy data: 2 points covering the limits to force the scale
  df <- data.frame(x = c(0, 1), y = c(0, 1), v = lims)
  ggplot(df, aes(x, y, fill = v)) +
    geom_tile() +
    scale_fill_gradient2(
      low      = PARAMS$pal_low,
      mid      = PARAMS$pal_mid,
      high     = PARAMS$pal_high,
      midpoint = (0+lims[1])/2,
      limits   = lims,
      oob      = scales::squish,
      name     = label,
      guide    = guide_colorbar(
        title.position = "top",
        title.hjust    = 0.5,
        barwidth       = unit(40, "mm"),
        barheight      = unit(2.5, "mm"),
        ticks.colour   = "grey30",
        frame.colour   = "grey30"
      )
    ) +
    theme_void(base_size = PARAMS$base_font_size) +
    theme(
      legend.position    = "bottom",
      legend.title       = element_text(size = PARAMS$base_font_size,
                                        face = "bold"),
      legend.text        = element_text(size = PARAMS$base_font_size - 1.5),
      legend.box.margin  = margin(0, 0, 0, 0),
      plot.margin        = margin(0, 0, 0, 0)
    ) +
    # hide the tile layer and axes: keep only the legend
    guides(fill = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(40, "mm"),
      barheight      = unit(2.5, "mm")
    )) +
    coord_cartesian(xlim = c(2, 3), ylim = c(2, 3))   # hors-cadre, rien ne se voit
}
legend_FRic <- make_legend_only(fric_lims, "FRic change (%)")

make_legend_only <- function(lims, label) {
  # dummy data: 2 points covering the limits to force the scale
  df <- data.frame(x = c(0, 1), y = c(0, 1), v = lims)
  ggplot(df, aes(x, y, fill = v)) +
    geom_tile() +
    scale_fill_gradient2(
      low      = PARAMS$pal_high,
      mid      = PARAMS$pal_mid,
      high     = PARAMS$pal_low,
      midpoint = (0+lims[2])/2,
      limits   = lims,
      oob      = scales::squish,
      name     = label,
      guide    = guide_colorbar(
        title.position = "top",
        title.hjust    = 0.5,
        barwidth       = unit(40, "mm"),
        barheight      = unit(2.5, "mm"),
        ticks.colour   = "grey30",
        frame.colour   = "grey30"
      )
    ) +
    theme_void(base_size = PARAMS$base_font_size) +
    theme(
      legend.position    = "bottom",
      legend.title       = element_text(size = PARAMS$base_font_size,
                                        face = "bold"),
      legend.text        = element_text(size = PARAMS$base_font_size - 1.5),
      legend.box.margin  = margin(0, 0, 0, 0),
      plot.margin        = margin(0, 0, 0, 0)
    ) +
    # hide the tile layer and axes: keep only the legend
    guides(fill = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(40, "mm"),
      barheight      = unit(2.5, "mm")
    )) +
    coord_cartesian(xlim = c(2, 3), ylim = c(2, 3))   # hors-cadre, rien ne se voit
}


legend_FUn  <- make_legend_only(fun_lims,  "FUn change (%)")


# ============================================================================
# 7. Ajout des titres de colonne et de ligne via patchwork ------------------
# ============================================================================

# Bandeau de titres de colonne (une fois en haut)
make_col_title <- function(label, color) {
  ggplot() +
    annotate("text", x = 0, y = 0,
             label = label,
             size  = PARAMS$base_font_size / 2.5,
             color = color, fontface = "bold") +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))
}

col_titles <- lapply(PARAMS$dim_order, function(d) {
  make_col_title(PARAMS$dim_labels[d], PARAMS$dim_colors[d])
})

# Étiquettes de ligne (rotation 90°)
make_row_label <- function(label) {
  ggplot() +
    annotate("text", x = 0, y = 0,
             label = label,
             size  = PARAMS$base_font_size / 2.5,
             angle = 90, fontface = "bold") +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))
}

row_label_FRic <- make_row_label("FRic change (%)")
row_label_FUn  <- make_row_label("FUn change (%)")


# ============================================================================
# 8. Assemblage final --------------------------------------------------------
# ============================================================================
# Layout :
#
#   [vide]   [Loco]   [Diet]   [Repro]      <- titres de colonne
#   [FRic]   m1       m2       m3            <- ligne FRic + cartes
#   [FUn ]   m4       m5       m6            <- ligne FUn  + cartes
#   [empty]           legends                <- legends at bottom
#
# patchwork avec wrap_plots et plot_layout

# Bandeau de titres
header_row <- patchwork::wrap_plots(
  c(list(patchwork::plot_spacer()), col_titles),
  nrow = 1,
  widths = c(0.06, 1, 1, 1)
)

# Ligne FRic
fric_row <- patchwork::wrap_plots(
  c(list(row_label_FRic),
    maps_FRic[PARAMS$dim_order]),
  nrow = 1,
  widths = c(0.06, 1, 1, 1)
)

# Ligne FUn
fun_row <- patchwork::wrap_plots(
  c(list(row_label_FUn),
    maps_FUn[PARAMS$dim_order]),
  nrow = 1,
  widths = c(0.06, 1, 1, 1)
)

# Empilement
maps_block <- header_row / fric_row / fun_row +
  patchwork::plot_layout(heights = c(0.06, 1, 1))

# Legends at the bottom, side by side
legends_block <- patchwork::wrap_plots(
  legend_FRic, legend_FUn,
  nrow = 1
)

final_fig <- header_row / fric_row / fun_row / legends_block +
  patchwork::plot_layout(heights = c(0.06, 1, 1, 0.18))  


# ============================================================================
# 9. Export -------------------------------------------------------------------
# ============================================================================
suffix  <- paste0("_IUCN", PARAMS$scenario_to_plot)
out_pdf <- file.path(PARAMS$out_dir,
                     paste0(PARAMS$fig_basename, suffix, ".pdf"))
ggsave(out_pdf, final_fig,
       width = PARAMS$width_mm, height = PARAMS$height_mm,
       units = "mm", device = cairo_pdf)
# ============================================================================
# 10. Reproducibility -------- --------------------------------------------------------
# ============================================================================
sessionInfo()