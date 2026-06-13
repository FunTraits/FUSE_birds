#-------------------------------------------------------------------------------
# Figure 4. Top FUSE species and biogeographical patterns of community FUSE.
#
# Layout (180 mm x 200 mm) :
#   Ligne 1 (a-c) : top 25 FUSE par dimension, barplots horizontaux
#                   coloured by IUCN category
#   Ligne 2 (d)   : UpSet plot du chevauchement entre les 3 listes top-25
#   Ligne 3 (e)   : trois cartes globales DGGS (top 10% mean) du FUSE
#                   FUSE, shared colour scale
#
# Author  : A. Toussaint
#
# Prerequisites:
#   - metrics_df : RDS object produced by A1_recompute_FUn_FSp_FUSE.R,
#                  with columns species, iucn, FUSE_loco, FUSE_diet,
#                  FUSE_repro
#   - sitesdggs7 : named list cell_id -> species vector
#   - dggs_grid  : sf of DGGS cells (with column 'seqnum')
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c(
  "ggplot2", "dplyr", "tidyr", "tibble", "purrr", "scales",
  "patchwork", "sf",
  "rnaturalearth", "rnaturalearthdata"
)
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# Polyfill %||% (cf. discussion script D)
if (!exists("%||%", envir = baseenv())) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  assign("%||%", `%||%`, envir = .GlobalEnv)
}


# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(
  
  # Input
  metrics_path   = "data/processed/species_metrics_FUn_FSp_FUSE.rds",
  
  # Suffixes des colonnes FUSE par dimension dans metrics_df
  fuse_cols      = c(locomotion = "FUSE_loco",
                     diet       = "FUSE_diet",
                     reproduction = "FUSE_repro"),
  
  # Top-N pour les barplots et l'UpSet
  top_n          = 25,
  
  # Top-X% of species in a cell for the community score
  top_pct        = 0.10,
  
  # Apparence
  base_font_size = 8,
  dim_order      = c("locomotion", "diet", "reproduction"),
  dim_labels     = c(locomotion   = "Locomotion",
                     diet         = "Diet",
                     reproduction = "Reproduction"),
  dim_colors     = c(locomotion   = "#2E7D32",
                     diet         = "#C62828",
                     reproduction = "#1565C0"),
  
  # Couleurs IUCN (palette IUCN officielle ; LC bleu-gris)
  iucn_levels    = c("LC", "NT", "VU", "EN", "CR"),
  iucn_colors    = c(LC = "#60C659", NT = "#CCE226",
                     VU = "#F9E814", EN = "#FC7F3F", CR = "#D81E05"),
  
  # Cartes
  crs_target     = "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs",
  map_palette    = "magma",     # viridisLite option
  clip_quantile  = 0.99,        # écrêtage pour limites de palette
  
  # Format figure
  out_dir        = "figures",
  fig_basename   = "Fig4_top_FUSE_overlap_maps",
  width_mm       = 180,
  height_mm      = 200,
  dpi            = 600
)
dir.create(PARAMS$out_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================================
# 2. USER INPUTS — adjust paths -------------------------------
# ============================================================================
metrics_df <- readRDS(PARAMS$metrics_path)
sitesdggs7 <- readRDS("data/raw/sitesdggs7.RDS")
library(dggridR)
dggs <- dgconstruct(res = 7)
grid <- dgearthgrid(dggs)
dggs_grid <- grid   # alias

stopifnot(
  exists("metrics_df"),
  all(c("species", "iucn", PARAMS$fuse_cols) %in% names(metrics_df)),
  exists("sitesdggs7"),
  exists("dggs_grid"),
  inherits(dggs_grid, "sf"),
  "seqnum" %in% names(dggs_grid)
)

# Force IUCN factor canonique
metrics_df <- metrics_df %>%
  mutate(iucn = factor(iucn, levels = PARAMS$iucn_levels))


# ============================================================================
# 3. (a-c) Top-25 FUSE par dimension : tables et barplots --------------------
# ============================================================================

#' For a given dimension, return a tibble of the top-N species
#' (species, FUSE, iucn) sorted by decreasing FUSE.
#' Names are abbreviated (G. species) to reduce Y-axis footprint
#' et permettre un meilleur alignement avec les autres panneaux.
get_top_n <- function(metrics_df, fuse_col, n = PARAMS$top_n) {
  metrics_df %>%
    select(species, iucn, !!fuse_col) %>%
    rename(FUSE = !!fuse_col) %>%
    arrange(desc(FUSE)) %>%
    slice_head(n = n) %>%
    mutate(
      # "Genus_species" -> "G. species"
      species_short = sub("^([A-Za-z])[a-z]+_", "\\1. ", species),
      species_label = factor(species_short, levels = rev(species_short))
    )
}

top_lists <- lapply(setNames(names(PARAMS$fuse_cols), names(PARAMS$fuse_cols)),
                    function(d) get_top_n(metrics_df, PARAMS$fuse_cols[[d]]))


#' Build the horizontal barplot for one dimension.
make_top_barplot <- function(top_df, dim_name, dim_color) {
  ggplot(top_df,
         aes(x = FUSE, y = species_label, fill = iucn)) +
    geom_col(color = "grey20", linewidth = 0.15, width = 0.75) +
    scale_fill_manual(values = PARAMS$iucn_colors,
                      drop = FALSE,
                      name = "IUCN") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = "FUSE", y = NULL, title = dim_name) +
    theme_minimal(base_size = PARAMS$base_font_size) +
    theme(
      plot.title         = element_text(size = PARAMS$base_font_size + 1,
                                        face = "bold", color = dim_color,
                                        hjust = 0,
                                        margin = margin(b = 2)),
      plot.title.position = "plot",   # aligne sur bord du dispositif
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.25),
      panel.background   = element_rect(fill = "white", color = NA),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.border       = element_blank(),
      axis.text.y        = element_text(size = PARAMS$base_font_size - 2,
                                        face = "italic", color = "grey20"),
      axis.text.x        = element_text(size = PARAMS$base_font_size - 1.5,
                                        color = "grey20"),
      axis.title.x       = element_text(size = PARAMS$base_font_size - 0.5),
      axis.ticks.x       = element_line(color = "grey60", linewidth = 0.25),
      legend.position    = "none",   # collectée en bas via patchwork
      plot.margin        = margin(2, 4, 2, 2)
    )
}

panels_top <- lapply(names(top_lists), function(d) {
  make_top_barplot(top_lists[[d]],
                   dim_name  = PARAMS$dim_labels[d],
                   dim_color = PARAMS$dim_colors[d])
})
names(panels_top) <- names(top_lists)


# ============================================================================
# 4. (d) Diagramme de Venn 3-ensembles du chevauchement ----------------------
# ============================================================================
# Note : le UpSet de ComplexUpset n'est pas compatible avec ggplot2 4.0
# (theme validation S7). For 3 sets, a Venn is in any case
# more readable. Built in pure ggplot2 to ensure compatibility.

# Construire une matrice membership : species x dimension (TRUE/FALSE)
top_species_all <- unique(unlist(lapply(top_lists, function(t) t$species)))
membership_df <- tibble(species = top_species_all)
for (d in names(top_lists)) {
  membership_df[[PARAMS$dim_labels[d]]] <- membership_df$species %in%
    top_lists[[d]]$species
}
membership_df <- membership_df %>%
  left_join(metrics_df %>% select(species, iucn), by = "species")

# Count species in each of the 7 regions of the 3-set Venn diagram
loco_lab  <- PARAMS$dim_labels[["locomotion"]]
diet_lab  <- PARAMS$dim_labels[["diet"]]
repro_lab <- PARAMS$dim_labels[["reproduction"]]

L <- membership_df[[loco_lab]]
D <- membership_df[[diet_lab]]
R <- membership_df[[repro_lab]]

n_only_L  <- sum( L & !D & !R)
n_only_D  <- sum(!L &  D & !R)
n_only_R  <- sum(!L & !D &  R)
n_LD      <- sum( L &  D & !R)
n_LR      <- sum( L & !D &  R)
n_DR      <- sum(!L &  D &  R)
n_LDR     <- sum( L &  D &  R)

# Centres of the three circles (equilateral triangle)
r <- 1.4
cx_L <- -0.85; cy_L <- -0.5
cx_D <-  0.85; cy_D <- -0.5
cx_R <-  0.0;  cy_R <-  0.85

circles <- tibble(
  x0 = c(cx_L, cx_D, cx_R),
  y0 = c(cy_L, cy_D, cy_R),
  r  = r,
  set = c("Locomotion", "Diet", "Reproduction"),
  fill_color = c(PARAMS$dim_colors[["locomotion"]],
                 PARAMS$dim_colors[["diet"]],
                 PARAMS$dim_colors[["reproduction"]])
)

# Helper: generate points of a circle
make_circle <- function(x0, y0, r, n = 200) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(x = x0 + r * cos(theta), y = y0 + r * sin(theta))
}
circle_points <- bind_rows(lapply(seq_len(nrow(circles)), function(i) {
  make_circle(circles$x0[i], circles$y0[i], circles$r[i]) %>%
    mutate(set = circles$set[i],
           fill_color = circles$fill_color[i])
}))

# Position of count labels in each region
# (determined geometrically by trial-and-error for 3 circles r=1.4)
labels_df <- tibble(
  x = c(-1.6, 1.6, 0.0, 0.0, -0.85, 0.85, 0.0),
  y = c(-0.85, -0.85, 1.45, -0.95, 0.45, 0.45, 0.0),
  count = c(n_only_L, n_only_D, n_only_R, n_LD, n_LR, n_DR, n_LDR),
  region = c("Locomotion only", "Diet only", "Reproduction only",
             "Loco ∩ Diet", "Loco ∩ Repro", "Diet ∩ Repro",
             "All three")
)

# Étiquettes des sets (au-dessus de chaque cercle)
set_labels <- tibble(
  x = c(cx_L - 0.7, cx_D + 0.7, cx_R),
  y = c(cy_L - r - 0.15, cy_D - r - 0.15, cy_R + r + 0.15),
  set = c("Locomotion", "Diet", "Reproduction"),
  color = c(PARAMS$dim_colors[["locomotion"]],
            PARAMS$dim_colors[["diet"]],
            PARAMS$dim_colors[["reproduction"]])
)

panel_venn <- ggplot() +
  # Cercles transparents par dimension
  geom_polygon(data = circle_points,
               aes(x = x, y = y, group = set, fill = set),
               alpha = 0.25, color = NA) +
  geom_path(data = circle_points,
            aes(x = x, y = y, group = set, color = set),
            linewidth = 0.5) +
  # Comptes
  geom_text(data = labels_df,
            aes(x = x, y = y, label = count),
            size = PARAMS$base_font_size / 2,
            color = "grey15") +
  # Étiquettes de set
  geom_text(data = set_labels,
            aes(x = x, y = y, label = set, color = set),
            size = PARAMS$base_font_size / 3,
            show.legend = FALSE) +
  scale_fill_manual(values = setNames(circles$fill_color, circles$set),
                    guide = "none") +
  scale_color_manual(values = setNames(circles$fill_color, circles$set),
                     guide = "none") +
  coord_equal(xlim = c(-2.6, 2.6), ylim = c(-2.4, 2.6)) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(2, 2, 2, 2))


# ============================================================================
# 5. (e) Cartes communautaires FUSE par dimension ----------------------------
# ============================================================================

# Preparation: for each cell and dimension, compute the mean
# of the top 10% FUSE of present species.

#' Score communautaire : moyenne du top X% des FUSE dans une cellule.
community_score <- function(species_in_cell, fuse_lookup,
                            top_pct = PARAMS$top_pct) {
  fuse_vals <- fuse_lookup[species_in_cell]
  fuse_vals <- fuse_vals[!is.na(fuse_vals)]
  if (length(fuse_vals) == 0) return(NA_real_)
  k <- max(1, ceiling(length(fuse_vals) * top_pct))
  mean(sort(fuse_vals, decreasing = TRUE)[1:k])
}

# Filter assemblages to species with defined FUSE
pool_species <- metrics_df$species
sitesdggs7_clean <- lapply(sitesdggs7, function(sp) intersect(sp, pool_species))
sitesdggs7_clean <- sitesdggs7_clean[lengths(sitesdggs7_clean) > 0]

# Compute per dimension
cell_fuse <- tibble(cell = names(sitesdggs7_clean))

for (d in names(PARAMS$fuse_cols)) {
  fuse_lookup <- setNames(metrics_df[[PARAMS$fuse_cols[[d]]]],
                          metrics_df$species)
  cell_fuse[[paste0("FUSE_", d)]] <- vapply(
    sitesdggs7_clean,
    function(sp) community_score(sp, fuse_lookup),
    numeric(1)
  )
}

# Join DGGS geometry
cell_fuse$seqnum <- as.numeric(cell_fuse$cell)
cell_sf <- dggs_grid %>%
  dplyr::left_join(cell_fuse, by = "seqnum")

# Wrap dateline et reprojection
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
# 5b. Boxplot per continent (secondary panel beside the Venn) -------- --------------
# ============================================================================
# Approach: join DGGS cells to their continent IN LON/LAT (WGS84),
# avant la reprojection Robinson — la combinaison st_wrap_dateline +
# st_transform(Robinson) produces invalid geometries that cause
# GEOS lors d'un st_join post-projection ("TopologyException: side location
# conflict").
#
# Use cell_sf (still in WGS84 at this stage) and the raw Natural Earth
# world (unprojected).

# Cells in WGS84 (cell_sf is still in lon/lat at this script stage)
cell_centroids_ll <- sf::st_centroid(cell_sf)

# Continents NaturalEarth en WGS84
world_ll <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")
continents_ll <- world_ll %>%
  dplyr::select(continent, geometry) %>%
  dplyr::filter(!is.na(continent), continent != "Seven seas (open ocean)")

# Safety: if WGS84 geometries are invalid (rare at this resolution
# mais possible), on les corrige avant join.
if (any(!sf::st_is_valid(continents_ll))) {
  continents_ll <- sf::st_make_valid(continents_ll)
}

# Disable S2 if active (extra safety for st_join in lon/lat)
old_s2 <- sf::sf_use_s2()
suppressMessages(sf::sf_use_s2(FALSE))

cell_with_continent <- sf::st_join(cell_centroids_ll, continents_ll,
                                   join = sf::st_within, left = FALSE)

# Restore S2 setting
suppressMessages(sf::sf_use_s2(old_s2))

# Ordre des continents : ranking par richesse en oiseaux (approximatif)
continent_order <- c("South America", "Africa", "Asia",
                     "North America", "Oceania", "Europe", "Antarctica")

# Format long pour ggplot : une ligne par cellule × dimension
boxplot_df <- cell_with_continent %>%
  sf::st_drop_geometry() %>%
  dplyr::select(continent,
                FUSE_locomotion, FUSE_diet, FUSE_reproduction) %>%
  tidyr::pivot_longer(starts_with("FUSE_"),
                      names_to = "dimension", values_to = "FUSE") %>%
  dplyr::mutate(
    dimension = sub("FUSE_", "", dimension),
    dimension = factor(dimension,
                       levels = PARAMS$dim_order,
                       labels = unname(PARAMS$dim_labels[PARAMS$dim_order])),
    continent = factor(continent,
                       levels = intersect(continent_order, unique(continent)))
  ) %>%
  dplyr::filter(!is.na(continent), !is.na(FUSE))

# --- Summary statistics per continent x dimension: mean and SE -------- -------
boxplot_summary <- boxplot_df %>%
  dplyr::group_by(continent, dimension) %>%
  dplyr::summarise(
    mean_FUSE = mean(FUSE, na.rm = TRUE),
    n         = sum(!is.na(FUSE)),
    se_FUSE   = sd(FUSE, na.rm = TRUE) / sqrt(pmax(n, 1)),
    .groups   = "drop"
  )

# Dimension colours (named mapping for ggplot)
dim_colors_named <- setNames(
  unname(PARAMS$dim_colors[PARAMS$dim_order]),
  unname(PARAMS$dim_labels[PARAMS$dim_order])
)

panel_boxplot <- ggplot(
  boxplot_df,
  aes(x = continent, y = FUSE,
      color = dimension, fill = dimension,
      group = interaction(continent, dimension))
) +
  # --- Jitter de points individuels (1 point = 1 cellule DGGS) ---
  geom_point(
    position = position_jitterdodge(
      jitter.width  = 0.18,
      jitter.height = 0,
      dodge.width   = 0.85
    ),
    size = 0.35, alpha = 0.25, shape = 16
  ) +
  # --- Hollow boxplot (no fill), coloured border ---
  geom_boxplot(
    fill = NA,
    outlier.shape = NA,           # outliers déjà visibles via les points
    linewidth = 0.45,
    width = 0.7,
    position = position_dodge(width = 0.85),
    alpha = 1
  ) +
  # --- Point central : moyenne par continent x dimension ---
  geom_point(
    data = boxplot_summary,
    aes(x = continent, y = mean_FUSE,
        fill = dimension, group = dimension),
    color = "black",
    shape = 21,
    size = 2,
    stroke = 0.4,
    position = position_dodge(width = 0.85),
    inherit.aes = FALSE
  ) +
  # --- Barres d'erreur sur la moyenne (SE) ---
  geom_errorbar(
    data = boxplot_summary,
    aes(x = continent,
        ymin = mean_FUSE - se_FUSE,
        ymax = mean_FUSE + se_FUSE,
        group = dimension),
    color = "black",
    width = 0,
    linewidth = 0.4,
    position = position_dodge(width = 0.85),
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = dim_colors_named, name = NULL) +
  scale_fill_manual( values = dim_colors_named, name = NULL,
                     guide = "none") +
  labs(x = NULL, y = "Community FUSE") +
  guides(color = guide_legend(override.aes = list(
    size = 2.5, alpha = 1, shape = 15
  ))) +
  theme_minimal(base_size = PARAMS$base_font_size) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.border       = element_rect(color = "grey60", fill = NA,
                                      linewidth = 0.3),
    axis.text.x        = element_text(size = PARAMS$base_font_size - 1.5,
                                      angle = 30, hjust = 1, color = "grey20"),
    axis.text.y        = element_text(size = PARAMS$base_font_size - 1.5,
                                      color = "grey20"),
    axis.title.y       = element_text(size = PARAMS$base_font_size - 0.5),
    axis.ticks         = element_line(color = "grey60", linewidth = 0.25),
    axis.ticks.length  = unit(1, "mm"),
    legend.position    = "top",
    legend.text        = element_text(size = PARAMS$base_font_size - 1),
    legend.key.size    = unit(3, "mm"),
    legend.box.margin  = margin(-2, 0, -2, 0),
    plot.margin        = margin(2, 4, 2, 2)
  )

# Échelle commune sur les trois dimensions
fuse_vals_all <- unlist(lapply(names(PARAMS$fuse_cols),
                               function(d) cell_sf_proj[[paste0("FUSE_", d)]]))
fuse_lim_hi <- quantile(fuse_vals_all, PARAMS$clip_quantile, na.rm = TRUE)
fuse_lim_lo <- min(fuse_vals_all, na.rm = TRUE)

#' Build a map for a given dimension.
make_fuse_map <- function(sf_data, fuse_col, lims, dim_name, dim_color) {
  ggplot() +
    geom_sf(data = world, fill = "grey92", color = "grey75",
            linewidth = 0.1) +
    geom_sf(data = sf_data,
            aes(fill = .data[[fuse_col]]),
            color = NA) +
    scale_fill_viridis_c(option = PARAMS$map_palette,
                         limits = lims, oob = scales::squish,
                         na.value = "grey95",
                         name = "Community FUSE",
                         guide = guide_colorbar(
                           title.position = "bottom",
                           title.hjust    = 0.5,
                           barwidth       = unit(40, "mm"),
                           barheight      = unit(2.5, "mm")
                         )) +
    coord_sf(crs = PARAMS$crs_target, expand = FALSE) +
    labs(caption = dim_name) +
    theme_void(base_size = PARAMS$base_font_size) +
    theme(
      plot.caption       = element_text(size = PARAMS$base_font_size + 3,
                                        face = "bold", color = dim_color,
                                        hjust = 0.5, vjust = -5,
                                        margin = margin(b = 1)),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position  = "none",     # collectée séparément
      plot.margin      = margin(1, 1, 1, 1)
    )
}

panels_maps <- lapply(names(PARAMS$fuse_cols), function(d) {
  make_fuse_map(
    cell_sf_proj,
    fuse_col  = paste0("FUSE_", d),
    lims      = c(fuse_lim_lo, fuse_lim_hi),
    dim_name  = PARAMS$dim_labels[d],
    dim_color = PARAMS$dim_colors[d]
  )
})
names(panels_maps) <- names(PARAMS$fuse_cols)


# ============================================================================
# 6. Build global legends -------- --------------------------------------
# ============================================================================

#' Standalone IUCN legend (used for barplots and UpSet)
make_iucn_legend <- function() {
  df <- tibble(x = 1, y = 1,
               iucn = factor(PARAMS$iucn_levels, levels = PARAMS$iucn_levels))
  ggplot(df, aes(x, y, fill = iucn)) +
    geom_tile() +
    scale_fill_manual(values = PARAMS$iucn_colors, drop = FALSE,
                      name = "IUCN status") +
    coord_cartesian(xlim = c(2, 3), ylim = c(2, 3)) +
    theme_void(base_size = PARAMS$base_font_size) +
    theme(legend.position = "bottom",
          legend.title    = element_text(size = PARAMS$base_font_size,
                                         face = "bold"),
          legend.text     = element_text(size = PARAMS$base_font_size - 1),
          legend.key.size = unit(3, "mm"),
          plot.margin     = margin(0, 0, 0, 0))
}

#' Community FUSE legend (continuous, viridis magma)
make_fuse_map_legend <- function() {
  df <- data.frame(x = c(0, 1), y = c(0, 1),
                   v = c(fuse_lim_lo, fuse_lim_hi))
  ggplot(df, aes(x, y, fill = v)) +
    geom_tile() +
    scale_fill_viridis_c(option = PARAMS$map_palette,
                         limits = c(fuse_lim_lo, fuse_lim_hi),
                         name = "Community FUSE",
                         guide = guide_colorbar(
                           title.position = "top",
                           title.hjust    = 0.5,
                           barwidth       = unit(40, "mm"),
                           barheight      = unit(2.5, "mm")
                         )) +
    coord_cartesian(xlim = c(2, 3), ylim = c(2, 3)) +
    theme_void(base_size = PARAMS$base_font_size) +
    theme(legend.position = "bottom",
          legend.title    = element_text(size = PARAMS$base_font_size,
                                         face = "bold"),
          legend.text     = element_text(size = PARAMS$base_font_size - 1.5),
          plot.margin     = margin(0, 0, 0, 0))
}

iucn_legend     <- make_iucn_legend()
fuse_map_legend <- make_fuse_map_legend()


# ============================================================================
# 7. Assemblage final --------------------------------------------------------
# ============================================================================

# Row 1: barplots a-b-c (titles already coloured by dimension in make_top_barplot)
row_top <- patchwork::wrap_plots(
  panels_top$locomotion + labs(title = "(a) Locomotion"),
  panels_top$diet         + labs(title = "(b) Diet"),
  panels_top$reproduction + labs(title = "(c) Reproduction"),
  nrow = 1
)

# Ligne 2 : Venn (gauche) + Boxplot (droite), avec titres natifs
# On attache les titres directement aux panneaux via labs() + theme(plot.title)
# pour que l'alignement soit relatif au panneau (pas au dispositif global).

panel_venn_titled <- panel_venn +
  labs(title = "(d) Top-25 overlap") +
  theme(plot.title = element_text(size = PARAMS$base_font_size + 1,
                                  face = "bold", hjust = 0,
                                  margin = margin(b = 2)),
        plot.title.position = "plot")   # aligne sur le bord du PLOT, pas du panel

panel_boxplot_titled <- panel_boxplot +
  labs(title = "(e) Community FUSE by continent") +
  theme(plot.title = element_text(size = PARAMS$base_font_size + 1,
                                  face = "bold", hjust = 0,
                                  margin = margin(b = 2)),
        plot.title.position = "plot")

row_venn_box <- patchwork::wrap_plots(
  panel_venn_titled, panel_boxplot_titled,
  nrow = 1,
  widths = c(0.42, 0.58)
)

# Ligne 3 : cartes avec titre natif sur la 1re carte qui couvre toute la ligne.
# Astuce : on attache un titre composite au row entier via plot_annotation,
# but this raises the same alignment issue. Simpler: title as
# subtitle on the leftmost panel, which aligns with the
# bord gauche du panneau, donc visuellement OK car la carte de gauche n'a
# no Y axis that offsets it.
panel_loco_map_titled <- panels_maps$locomotion +
  labs(subtitle = "(f) Community-level FUSE (mean of top 10% species per cell)") +
  theme(plot.subtitle = element_text(size = PARAMS$base_font_size + 1,
                                     face = "bold", hjust = 0,
                                     color = "black",
                                     margin = margin(t = 12, b = 0)),
        plot.title.position = "plot")

row_maps <- patchwork::wrap_plots(
  panel_loco_map_titled,
  panels_maps$diet,
  panels_maps$reproduction,
  nrow = 1
)

# Row 4: legends side by side
row_legends <- patchwork::wrap_plots(
  iucn_legend, fuse_map_legend,
  nrow = 1
)

# Stack (no more interleaved title bands: native titles on panels)
final_fig <- (row_top /
                row_venn_box /
                row_maps /
                row_legends) +
  patchwork::plot_layout(
    heights = c(1.8, 1.4, 1.5, 0.2)
  ) +
  patchwork::plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  )


# ============================================================================
# 8. Export -------------------------------------------------------------------
# ============================================================================
out_pdf <- file.path(PARAMS$out_dir, paste0(PARAMS$fig_basename, ".pdf"))

ggsave(out_pdf, final_fig,
       width = PARAMS$width_mm, height = PARAMS$height_mm,
       units = "mm", device = cairo_pdf)


# ============================================================================
# 9. Summary table: top-25 and overlaps -------- --------------------------
# ============================================================================
n_overlap <- list(
  loco_only        = sum(membership_df$Locomotion &
                           !membership_df$Diet & !membership_df$Reproduction),
  diet_only        = sum(!membership_df$Locomotion &
                           membership_df$Diet & !membership_df$Reproduction),
  repro_only       = sum(!membership_df$Locomotion &
                           !membership_df$Diet &  membership_df$Reproduction),
  loco_diet        = sum(membership_df$Locomotion &
                           membership_df$Diet & !membership_df$Reproduction),
  loco_repro       = sum(membership_df$Locomotion &
                           !membership_df$Diet &  membership_df$Reproduction),
  diet_repro       = sum(!membership_df$Locomotion &
                           membership_df$Diet &  membership_df$Reproduction),
  all_three        = sum(membership_df$Locomotion &
                           membership_df$Diet &  membership_df$Reproduction)
)
message("\nChevauchement des top-", PARAMS$top_n, " FUSE :")
print(unlist(n_overlap))
message(sprintf("Total d'espèces uniques dans l'union des top-%d : %d",
                PARAMS$top_n, nrow(membership_df)))


# ============================================================================
# 10. Reproducibility -------- --------------------------------------------------------
# ============================================================================
sessionInfo()