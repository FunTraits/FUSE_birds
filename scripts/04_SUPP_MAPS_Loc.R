#-------------------------------------------------------------------------------
# SUPP — Global DGGS maps of community priority indices
#         (FUS_sum3 and FUSE_sum3)
#
# Aggregates species-level scores (FUS / FUSE) to DGGS cells,
# produces two global maps on a shared colour scale, and a
# cell-level scatter plot for comparison.
#
# Prerequisites — the following objects must exist:
#   comm_long : named list (cell -> species vector)
#   site      : data.frame with columns long, lat, order, hole, piece, group, cell
#   fd_df     : species data.frame with morpho_FNo, lifehistory_FNo, diet_FNo,
#               morpho_FUSE, lifehistory_FUSE, diet_FUSE
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c("dplyr", "tidyr", "purrr", "sf",
                   "ggplot2", "patchwork", "readr", "rnaturalearth")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))


# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(

  # Variable to map ("FUS_mean_top10" or "FUS_mean_all", etc.)
  map_var_FUS  = "FUS_mean_top10",
  map_var_FUSE = "FUSE_mean_top10",

  # Fraction of top species for the cell-level mean
  top_prop     = 0.10,

  # Palette: green -> yellow -> red
  palette      = c("#006837", "#FFFFBF", "#A50026"),

  out_dir      = "results",
  fig_dir      = "results"
)
dir.create(PARAMS$out_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================================
# 2. USER INPUTS — adjust paths --------------------------------
# ============================================================================
# comm_long, site, and fd_df must be loaded upstream
stopifnot(
  exists("comm_long"), is.list(comm_long),
  exists("site"),      is.data.frame(site),
  exists("fd_df"),
  all(c("species",
        "morpho_FNo", "lifehistory_FNo", "diet_FNo",
        "morpho_FUSE", "lifehistory_FUSE", "diet_FUSE") %in% names(fd_df))
)


# ============================================================================
# 3. Helper functions -------- ----------------------------------------------------
# ============================================================================

minmax01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2])) return(rep(NA_real_, length(x)))
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

rescale01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (!is.finite(r[1]) || !is.finite(r[2])) return(rep(NA_real_, length(x)))
  if (r[1] == r[2]) return(rep(0.5, length(x)))
  (x - r[1]) / (r[2] - r[1])
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else max(x)
}

mean_top_prop <- function(x, prop = 0.10) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  k <- max(1, ceiling(length(x) * prop))
  mean(sort(x, decreasing = TRUE)[seq_len(k)])
}

# Build a coordinate ring (DGGS polygon)
make_ring <- function(df_ring) {
  df_ring <- df_ring %>% arrange(order)
  coords  <- as.matrix(df_ring[, c("long", "lat")])
  if (!all(coords[1, ] == coords[nrow(coords), ])) {
    coords <- rbind(coords, coords[1, ])
  }
  coords
}

# Build the sf geometry of a DGGS cell
build_cell_geom <- function(df_cell) {
  units <- df_cell %>% group_by(piece, group) %>% group_split()

  polys <- purrr::map(units, function(u) {
    outers <- u %>% filter(!hole)
    if (nrow(outers) == 0) outers <- u
    outer_ring <- make_ring(outers)
    holes      <- u %>% filter(hole)
    hole_rings <- if (nrow(holes) > 0) list(make_ring(holes)) else list()
    c(list(outer_ring), hole_rings)
  })

  sf::st_multipolygon(polys)
}


# ============================================================================
# 4. FUS / FUSE scores per species -------- --------------------------------------------
# ============================================================================

scores <- fd_df %>%
  transmute(
    species,
    FUS_sum3  = minmax01(morpho_FNo)  + minmax01(lifehistory_FNo)  + minmax01(diet_FNo),
    FUSE_sum3 = minmax01(morpho_FUSE) + minmax01(lifehistory_FUSE) + minmax01(diet_FUSE)
  )


# ============================================================================
# 5. Aggregation to DGGS cell level -------- --------------------------------
# ============================================================================

comm_tbl <- tibble(
  cell    = names(comm_long),
  species = comm_long
) %>%
  unnest(species) %>%
  mutate(cell = as.character(cell))

cell_scores <- comm_tbl %>%
  left_join(scores, by = "species") %>%
  group_by(cell) %>%
  summarise(
    S               = n_distinct(species),
    FUS_mean_all    = mean(FUS_sum3,  na.rm = TRUE),
    FUSE_mean_all   = mean(FUSE_sum3, na.rm = TRUE),
    FUS_max         = safe_max(FUS_sum3),
    FUSE_max        = safe_max(FUSE_sum3),
    FUS_mean_top10  = mean_top_prop(FUS_sum3,  prop = PARAMS$top_prop),
    FUSE_mean_top10 = mean_top_prop(FUSE_sum3, prop = PARAMS$top_prop),
    n_scored        = sum(is.finite(FUS_sum3) & is.finite(FUSE_sum3)),
    .groups         = "drop"
  )

write_csv(cell_scores, file.path(PARAMS$out_dir, "cell_scores_FUS_FUSE_sum3.csv"))


# ============================================================================
# 6. Build DGGS sf polygons -------- --------------------------------------
# ============================================================================

site2 <- site %>%
  mutate(cell  = as.character(cell),
         piece = as.character(piece),
         group = as.character(group),
         hole  = as.logical(hole))

cell_sf <- site2 %>%
  group_by(cell) %>%
  group_split() %>%
  purrr::map_dfr(function(df_cell) {
    sf::st_sf(
      cell     = df_cell$cell[1],
      geometry = sf::st_sfc(build_cell_geom(df_cell), crs = 4326)
    )
  })

map_sf <- cell_sf %>%
  left_join(cell_scores, by = "cell") %>%
  mutate(
    FUS_mean_top10  = rescale01(.data[[PARAMS$map_var_FUS]]),
    FUSE_mean_top10 = rescale01(.data[[PARAMS$map_var_FUSE]])
  )


# ============================================================================
# 7. Figures — maps + scatter plot -------- --------------------------------------
# ============================================================================

land <- ne_countries(scale = "medium", returnclass = "sf")

lims <- range(c(map_sf[[PARAMS$map_var_FUS]], map_sf[[PARAMS$map_var_FUSE]]),
              na.rm = TRUE)

base_theme <- theme_void(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0),
    plot.subtitle    = element_text(hjust = 0),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.position  = "right"
  )

p_fus <- ggplot() +
  geom_sf(data = land, fill = "grey85", colour = NA) +
  geom_sf(data = map_sf, aes(fill = .data[[PARAMS$map_var_FUS]]), colour = NA) +
  coord_sf(expand = FALSE, datum = NA) +
  scale_fill_gradientn(colours = PARAMS$palette, limits = lims,
                       na.value = "white", name = NULL) +
  labs(title = "FUS (community index)", subtitle = PARAMS$map_var_FUS) +
  base_theme

p_fuse <- ggplot() +
  geom_sf(data = land, fill = "grey85", colour = NA) +
  geom_sf(data = map_sf, aes(fill = .data[[PARAMS$map_var_FUSE]]), colour = NA) +
  coord_sf(expand = FALSE, datum = NA) +
  scale_fill_gradientn(colours = PARAMS$palette, limits = lims,
                       na.value = "white", name = NULL) +
  labs(title = "FUSE (community index)", subtitle = PARAMS$map_var_FUSE) +
  base_theme

# Scatter plot: FUSE vs. FUS per cell
cell_xy <- map_sf %>%
  st_drop_geometry() %>%
  transmute(cell,
            FUS  = .data[[PARAMS$map_var_FUS]],
            FUSE = .data[[PARAMS$map_var_FUSE]]) %>%
  filter(is.finite(FUS), is.finite(FUSE))

sp  <- suppressWarnings(cor.test(cell_xy$FUS, cell_xy$FUSE,
                                 method = "spearman", exact = FALSE))
ann <- paste0("Spearman ρ = ", round(unname(sp$estimate), 2),
              "   p = ", format.pval(sp$p.value, digits = 2, eps = 1e-4),
              "\nN cells = ", nrow(cell_xy))

p_scatter <- ggplot(cell_xy, aes(x = FUS, y = FUSE)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.6) +
  geom_vline(xintercept = 0.5, linewidth = 0.4) +
  geom_hline(yintercept = 0.5, linewidth = 0.4) +
  geom_point(size = 1.4, alpha = 0.5) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(title    = "Cell-level comparison",
       subtitle = "Each point is a DGGS cell (community indices scaled to [0,1])",
       x = "FUS (scaled)", y = "FUSE (scaled)") +
  annotate("label", x = 0.02, y = 0.98, hjust = 0, vjust = 1,
           label = ann, size = 3) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

final_fig <- (p_fus | p_fuse) / p_scatter +
  plot_layout(heights = c(1, 2.5)) +
  plot_annotation(
    title    = "Global patterns of functional priority indices",
    subtitle = "Top: DGGS maps · Bottom: cell-level biplot (FUSE vs FUS)"
  )

ggsave(file.path(PARAMS$fig_dir, "maps_plus_biplot_FUS_FUSE.png"),
       final_fig, width = 14, height = 14, dpi = 320, bg = "white")

final_fig
