#-------------------------------------------------------------------------------
# Figure 1. Functional structure of birds across three independent ecological
# dimensions (locomotion, diet, reproduction).
#
# Author : A. Toussaint
#
# Output:
#   * three-panel main figure (TPD heat map + IUCN contours + key species labels)
#   * inset 3 x 3 Spearman correlation matrix of FUn across dimensions
#   * single composed figure exported in PDF, PNG and SVG
#
# Required inputs (see "USER INPUTS" block below):
#   - tpd_lists   : named list of three TPDs objects (one per dimension)
#                   produced by TPD::TPDs() on the 200 x 200 (= 40,000 cells) grid
#   - coords      : named list of three data.frames with columns
#                     species, PC1, PC2  (the 2D PCoA scores used to build the TPDs)
#   - species_df  : a data.frame with at least
#                     species, iucn (factor: "LC","NT","VU","EN","CR"),
#                     FUn_loco, FUn_diet, FUn_repro
#   - key_species : a data.frame with columns
#                     species, dimension ("locomotion"/"diet"/"reproduction"),
#                     label (display name)
#                   listing 4-6 emblematic species per dimension to annotate.
#
# All other choices (palette, axis ranges, contour levels) are parameterised
# at the top of the script.
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c(
  "ggplot2", "dplyr", "tidyr", "tibble", "purrr",
  "patchwork", "scales", "ggrepel",
  "MASS",          # kde2d for IUCN contours
  "viridisLite",   # heatmap palette
  "ks",            # bandwidth selection (Hpi)
  "ggpubr",        # ggarrange
  "rphylopic"      # bird silhouettes from phylopic.org
)
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))


# ============================================================================
# 1. Global parameters --------------------------------------------------------
# ============================================================================
PARAMS <- list(
  # output
  out_dir         = "figures",
  fig_basename    = "Fig1_functional_structure",
  width_mm        = 360,    # Science Advances full-width
  height_mm       = 140,     # short panel
  dpi             = 600,
  
  # heatmap
  heat_palette    = "mako",      # viridisLite option; alternative: "rocket"
  heat_alpha      = 1,
  n_grid_breaks   = 6,           # for axis ticks
  
  # IUCN contours
  iucn_levels     = c("LC", "NT", "VU", "EN", "CR"),
  iucn_colors     = c(LC = "#3B9AB2", NT = "#78B7C5",
                      VU = "#EBCC2A", EN = "#E1AF00", CR = "#F21A00"),
  contour_quantile = 0.5,        # contour drawn at this density quantile
  min_n_for_contour = 30,        # below this n, plot points instead of contour
  
  # labels
  label_size      = 2.4,
  label_box_padding = 0.4,
  point_size      = 0.6,
  point_alpha     = 0.7,
  
  # general
  base_font_size  = 8,
  panel_titles    = c(locomotion  = "Locomotion (morphology)",
                      diet        = "Diet (foraging)",
                      reproduction = "Reproduction (life history)"),
  dim_colors    = c(locomotion   = "#2E7D32",
                    diet         = "#C62828",
                    reproduction = "#1565C0"),
  
  # ----- Phylopic illustrations -----
  # For each dimension, emblematic species with:
  #   - name    : PhyloPic taxon (genus or species, to retrieve the silhouette)
  #   - species : binomial name "Genus_species" present in the dataset
  #               (to retrieve the TRUE position in the PCoA and draw
  #                a coloured dot + segment to the silhouette)
  #   - h : hauteur de la silhouette en fraction du span vertical (0.05-0.10)
  #
  # Silhouettes are automatically distributed along an elliptical ring
  # around the density cloud, each silhouette angle being determined
  # by the species position relative to the cloud centroid.
  phylopic_species = list(
    locomotion = list(
      list(name = "Struthio camelus",       species = "Struthio_camelus",       h = 0.07),
      list(name = "Rhea americana",         species = "Rhea_americana",         h = 0.07),
      list(name = "Aquila chrysaetos",      species = "Aquila_chrysaetos",      h = 0.06),
      list(name = "Spheniscus demersus",    species = "Spheniscus_demersus",    h = 0.06),
      list(name = "Diomedea exulans",       species = "Diomedea_exulans",       h = 0.07),
      list(name = "Apus apus",              species = "Apus_apus",              h = 0.05),
      list(name = "Passer domesticus",      species = "Passer_domesticus",      h = 0.05),
      list(name = "Anas platyrhynchos",     species = "Anas_platyrhynchos",     h = 0.05),
      list(name = "Gallus gallus",          species = "Gallus_gallus",          h = 0.06)
    ),
    diet = list(
      list(name = "Gyps fulvus",            species = "Gyps_fulvus",            h = 0.07),  # DSv (scavenger)
      list(name = "Aquila chrysaetos",      species = "Aquila_chrysaetos",      h = 0.06),  # DVend (vertebrate)
      list(name = "Ramphastos toco",        species = "Ramphastos_toco",        h = 0.07),  # DF (frugivore)
      list(name = "Trochilidae",            species = "Phaethornis_superciliosus", h = 0.05),  # DN (nectarivore)
      list(name = "Anas platyrhynchos",     species = "Anas_platyrhynchos",     h = 0.06),  # DP (plant)
      list(name = "Pelecanus onocrotalus",  species = "Pelecanus_onocrotalus",  h = 0.06),  # DVfish
      list(name = "Hirundo rustica",        species = "Hirundo_rustica",        h = 0.05),  # DI (insectivore aerial)
      list(name = "Loxia curvirostra",      species = "Loxia_curvirostra",      h = 0.05),  # DSd (granivore)
      list(name = "Turdus merula",          species = "Turdus_merula",          h = 0.05)   # généraliste
    ),
    reproduction = list(
      list(name = "Diomedea exulans",       species = "Diomedea_exulans",       h = 0.07),  # longévité extrême
      list(name = "Macrocephalon maleo",    species = "Macrocephalon_maleo",    h = 0.07),  # ponte unique
      list(name = "Strigops habroptila",    species = "Strigops_habroptila",    h = 0.07),  # slow LH
      list(name = "Passer domesticus",      species = "Passer_domesticus",      h = 0.05),  # fast LH
      list(name = "Gallus gallus",          species = "Gallus_gallus",          h = 0.06),  # large clutch
      list(name = "Aquila chrysaetos",      species = "Aquila_chrysaetos",      h = 0.06),  # K-strategist
      list(name = "Spheniscus demersus",    species = "Spheniscus_demersus",    h = 0.06),  # seabird
      list(name = "Apteryx",                species = "Apteryx_australis",      h = 0.06)   # kiwi (oeuf énorme)
    )
  ),
  # Silhouette colour (dark grey for sober rendering)
  silhouette_color = "grey20",
  silhouette_alpha = 0.85,
  # Dot + connecting segment to the true species position
  position_point_size  = 1.6,
  position_point_alpha = 1,
  segment_color        = "grey50",
  segment_size         = 0.2,
  segment_alpha        = 0.5,
  # Anneau de placement des silhouettes (en fraction du demi-span depuis le centre)
  # 1.0 = panel edge; > 1 = outside; < 1 = inside
  ring_radius_x        = 1.20,
  ring_radius_y        = 1.20,
  # Minimum angular spacing between adjacent silhouettes (degrees)
  min_angular_sep_deg  = 25,
  phylopic_cache_path  = "data/processed/phylopic_uuid_cache.rds"
)
dir.create(PARAMS$out_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================================
# 2. USER INPUTS — replace with paths to your saved RDS files -----------------
#    (these objects are produced by your existing pipeline)
# ============================================================================

# Example loaders (uncomment & adjust paths):
tpd_lists  <- readRDS("data/processed/tpd_lists.rds")
coords     <- readRDS("data/processed/pcoa_coords.rds")
species_df <- readRDS("data/processed/species_metrics_FUn_FSp_FUSE.rds")
# key_species <- read.csv("data/key_species_to_label.csv")

key_species = species_df[1:12,1]
key_species$dimension = rep(c("locomotion","diet","reproduction"),4)
key_species$label = unlist(lapply(strsplit(key_species$species,"_"),function(x){paste0(substr(x[1],1,1),". ",substr(x[2],1,3))}))

PCA <- readRDS("data/processed/PCA_Birds.rds")
shortNames <- read.csv("data/processed/Shortnames_Birds.csv")
type <- names(PCA)
title_space = c('Combined','Locomotion','Reproduction','Diet')

# Sanity checks
stopifnot(
  exists("tpd_lists"),
  all(c("locomotion", "diet", "reproduction") %in% names(tpd_lists)),
  exists("coords"),
  all(c("locomotion", "diet", "reproduction") %in% names(coords)),
  exists("species_df"),
  all(c("species", "iucn",
        "FUn_loco", "FUn_diet", "FUn_repro") %in% names(species_df)),
  exists("key_species"),
  all(c("species", "dimension", "label") %in% names(key_species))
)

# Force IUCN factor with the canonical ordering (drops DD/NE)
species_df <- species_df %>%
  filter(iucn %in% PARAMS$iucn_levels) %>%
  mutate(iucn = factor(iucn, levels = PARAMS$iucn_levels))


# ============================================================================
# 3. Helper functions ---------------------------------------------------------
# ============================================================================

PCA_plot <- function(PCoAPlot,PCoACorPlot,multAx1,multAx,title,legend){
  Plan12=ggplot(data = data.frame(PCoAPlot$vectors)) +
    geom_point(aes(x = Axis.1, y = Axis.2),col="grey89") +
    theme_classic()+
    ggtitle(title)+
    xlab(paste0("PCoA 1 (",round(PCoAPlot$values[1,2]*100,2),"%)"))+
    ylab(paste0("PCoA 2 (",round(PCoAPlot$values[2,2]*100,2),"%)"))+
    #coord_fixed() + ## need aspect ratio of 1!
    geom_segment(data = data.frame(PCoACorPlot),
                 aes(x = 0, xend = Axis.1*multAx, y = 0, yend = Axis.2*multAx),
                 arrow = arrow(length = unit(0.25, "cm")), colour = PCoACorPlot$color) +
    geom_text(data = data.frame(PCoACorPlot), 
              aes(x =  Axis.1*multAx1, y = Axis.2*multAx1, 
                  label =  names), colour = PCoACorPlot$color,
              size = 5) +
    theme(
      plot.title = element_text(
        size = 18,       # Increase size (adjust as needed)
        face = "bold"    # Make it bold
      )
    )
  
  if(!is.null(legend)){
    Plan12 = Plan12 + annotate(geom = "text",label = legend, x = Inf, y = -Inf, hjust = 1, vjust = -0.5)      # Make text bold)
  }
  
  
  
  Plan12
}

# ---------- utilities ----------
.sanitize_axes <- function(df) {
  stopifnot(all(c("Axis.1","Axis.2") %in% names(df)))
  df <- df[is.finite(df$Axis.1) & is.finite(df$Axis.2), , drop = FALSE]
  df$Axis.1 <- as.numeric(df$Axis.1)
  df$Axis.2 <- as.numeric(df$Axis.2)
  df
}

.safe_limits <- function(x, pad_factor = 0.08) {
  xr <- range(x, finite = TRUE)
  if (!is.finite(xr[1]) || !is.finite(xr[2])) return(c(-1, 1))
  if (xr[1] == xr[2]) {
    eps <- ifelse(abs(xr[1]) > 0, abs(xr[1]) * pad_factor, 1e-6)
    return(c(xr[1] - eps, xr[2] + eps))
  }
  width <- diff(xr)
  c(xr[1] - width * pad_factor, xr[2] + width * pad_factor)
}

# KDE grid using {ks}
.kde_grid <- function(df, gridsize = 181, H = NULL) {
  if (!requireNamespace("ks", quietly = TRUE)) {
    stop("Please install.packages('ks') to use the KDE-based density background.")
  }
  df <- .sanitize_axes(df)
  if (nrow(df) < 3) return(list(grid = NULL, x = NULL, y = NULL, z = NULL))
  
  xlim <- .safe_limits(df$Axis.1); ylim <- .safe_limits(df$Axis.2)
  
  X <- cbind(
    pmin(pmax(df$Axis.1, xlim[1]), xlim[2]),
    pmin(pmax(df$Axis.2, ylim[1]), ylim[2])
  )
  
  if (is.null(H)) {
    H <- try(ks::Hpi(x = X), silent = TRUE)
    if (inherits(H, "try-error") || any(!is.finite(as.vector(H)))) {
      H <- try(ks::Hscv(x = X), silent = TRUE)
    }
    if (inherits(H, "try-error") || any(!is.finite(as.vector(H)))) {
      v <- apply(X, 2, stats::var, na.rm = TRUE); v[v <= 0 | !is.finite(v)] <- 1e-6
      H <- diag(v)
    }
  }
  
  kde <- ks::kde(
    x = X, H = H,
    xmin = c(xlim[1], ylim[1]), xmax = c(xlim[2], ylim[2]),
    gridsize = c(gridsize, gridsize)
  )
  
  gx <- kde$eval.points[[1]]
  gy <- kde$eval.points[[2]]
  gz <- kde$estimate  # matrix [length(gx) x length(gy)]
  
  grid <- expand.grid(x = gx, y = gy)
  grid$z <- as.vector(gz)
  
  list(grid = grid, x = gx, y = gy, z = gz)
}

# Compute density cutoffs (levels) that enclose given probabilities on the grid
.hdr_levels_grid <- function(gx, gy, gz, probs = c(0.25, 0.5, 0.99)) {
  # assume regular grid (true for ks::kde)
  dx <- mean(diff(gx)); dy <- mean(diff(gy))
  w <- as.vector(gz)
  ord <- order(w, decreasing = TRUE)
  cum_mass <- cumsum(w[ord]) * dx * dy
  # normalize to 1 (numerical integration may be slightly off)
  cum_mass <- cum_mass / max(cum_mass, na.rm = TRUE)
  
  sapply(probs, function(p) {
    idx <- which(cum_mass >= p)[1]
    if (is.na(idx)) min(w, na.rm = TRUE) else w[ord][idx]
  })
}

# Build contour paths + label positions from a raster grid
.contour_data <- function(gx, gy, gz, levels, prob_labels) {
  cls <- contourLines(x = gx, y = gy, z = gz, levels = levels)  # note t()
  if (length(cls) == 0) return(list(lines = NULL, labels = NULL))
  lines_df <- do.call(rbind, lapply(seq_along(cls), function(i) {
    data.frame(x = cls[[i]]$x, y = cls[[i]]$y, level = cls[[i]]$level, id = i)
  }))
  labels_df <- do.call(rbind, lapply(seq_along(cls), function(i) {
    n <- length(cls[[i]]$x); j <- floor(n/2)
    data.frame(x = cls[[i]]$x[j], y = cls[[i]]$y[j], level = cls[[i]]$level, id = i)
  }))
  # map level -> probability text
  lvl_map <- setNames(prob_labels, levels)
  labels_df$prob <- lvl_map[as.character(labels_df$level)]
  list(lines = lines_df, labels = labels_df)
}

# ---------- funspace-style plot ----------
PCA_plot_funspace <- function(PCoAPlot, PCoACorPlot, multAx1, multAx, title, legend = NULL, 
                              colLeg = NULL,
                              probs = c(0.25, 0.5, 0.99),
                              gridsize = 181, H = NULL, bins_filled = 30,pts = T) {
  
  df <- data.frame(PCoAPlot$vectors)
  df <- .sanitize_axes(df)
  
  # KDE grid
  KG <- .kde_grid(df, gridsize = gridsize, H = H)
  
  # density background
  # continuous raster background — gradient palette white -> orange -> dark red
  # Inspired by RColorBrewer::OrRd (5 stops) but with a more saturated centre
  p <- ggplot() +
    geom_raster(
      data = KG$grid,
      aes(x = x, y = y, fill = z)
    ) +
    scale_fill_gradientn(
      colours = c(
        "#FFFFFF",  # white          (bord externe)
        "#FFEDA0",  # pale yellow    (transition douce)
        "#FEB24C",  # warm orange    (zone intermédiaire)
        "#FC8D59",  # medium orange-red
        "#EF6548",  # strong red-orange
        "#D7301F",  # deep red
        "#990000",  # dark brownish red
        "#7F0000"   # darkest center (cœur saturé)
      ),
      values = scales::rescale(c(0, 0.05, 0.20, 0.40, 0.60, 0.80, 0.95, 1)),
      guide = "none"
    )
  if (!is.null(pts)) {
    p <- p +
      geom_point(data = df, aes(Axis.1, Axis.2),
                 colour = "grey36", size = 0.3,alpha = 1/15)
  }
  
  
  # HDR probability contours + labels (0.25, 0.5, 0.95 by default)
  # Dark-red contours (instead of black) to harmonise with the palette
  levs <- .hdr_levels_grid(KG$x, KG$y, KG$z, probs = probs)
  cd <- .contour_data(KG$x, KG$y, KG$z, levels = as.numeric(levs),
                      prob_labels = paste0(probs))
  if (!is.null(cd$lines)) {
    p <- p +
      geom_path(data = cd$lines, aes(x, y, group = id),
                colour = "#7F0000", linewidth = 0.5)
  }
  if (!is.null(cd$labels)) {
    p <- p +
      geom_text(data = cd$labels, aes(x, y, label = prob),
                colour = "#7F0000", size = 3)
  }
  

  # trait loading arrows + labels
  p <- p +
    geom_segment(
      data = data.frame(PCoACorPlot),
      aes(x = 0, xend = Axis.1 * multAx, y = 0, yend = Axis.2 * multAx),
      arrow = arrow(length = unit(0.25, "cm")),
      colour = PCoACorPlot$color, linewidth = 0.5
    ) +
    geom_text(
      data = data.frame(PCoACorPlot),
      aes(x = Axis.1 * multAx1, y = Axis.2 * multAx1, label = names),
      colour = PCoACorPlot$color, size = 5
    ) +
    theme_classic() +
    ggtitle(title) +
    xlab(paste0("PCoA 1 (", round(PCoAPlot$values[1, 2] * 100, 2), "%)")) +
    ylab(paste0("PCoA 2 (", round(PCoAPlot$values[2, 2] * 100, 2), "%)")) +
    coord_fixed() +
    theme(plot.title = element_text(size = 18, face = "bold",colour  = colLeg))
  
  if (!is.null(legend)) {
    p <- p + annotate("text", label = legend, x = Inf, y = -Inf, hjust = 1, vjust = -0.5)
  }
  
  p
}


## 2.1 Generate a PCA correlation plot with density legend
make_pca_plot <- function(pca_object, correlation_data, multAx1, multAx, title) {
  legend_text <- paste0(
    "Space occupied by:\n",
    "50% = ", round(pca_object$ALLDensity[1, "0.5"], 2), "\n",
    "99% = ", round(pca_object$ALLDensity[1, "0.99"], 2)
  )
  PCA_plot(pca_object$PCoA, correlation_data, multAx1, multAx, title, legend_text)
}
make_pca_plot_34 <- function(pca_object, correlation_data, multAx1, multAx, title) {
  legend_text <- paste0(
    "Space occupied by:\n",
    "50% = ", round(pca_object$ALLDensity[1, "0.5"], 2), "\n",
    "99% = ", round(pca_object$ALLDensity[1, "0.99"], 2)
  )
  PCA_plot_34(pca_object$PCoA, correlation_data, multAx1, multAx, title, legend_text)
}

## 2.2 Format PCA correlation table with colors and labels
format_correlation_table <- function(PCoACorPlot, shortNames) {
  cbind.data.frame(
    PCoACorPlot,
    color = shortNames[match(rownames(PCoACorPlot), shortNames$original), "color"],
    names = shortNames[match(rownames(PCoACorPlot), shortNames$original), "short"]
  )
}


## 2.3 Phylopic silhouette helpers --------------------------------------------
#' Retrieve and cache PhyloPic UUIDs for a list of names
#'
#' La fonction interroge l'API phylopic via rphylopic::get_uuid().
#' If a lookup fails (taxon not found), the UUID is NA and a warning
#' is raised. The cache is saved to disk to avoid re-querying
#' on every run.
fetch_phylopic_uuids <- function(species_lists, cache_path) {
  
  # Charger le cache existant
  cache <- if (file.exists(cache_path)) readRDS(cache_path) else list()
  
  # Flatten the list to retrieve all unique names
  all_names <- unique(unlist(lapply(species_lists, function(lst) {
    sapply(lst, `[[`, "name")
  })))
  
  for (nm in all_names) {
    if (!is.null(cache[[nm]]) && !is.na(cache[[nm]])) next  # déjà en cache
    uuid <- tryCatch(
      rphylopic::get_uuid(name = nm, n = 1L),
      error = function(e) {
        warning("Phylopic : pas d'UUID trouvé pour '", nm, "'. ",
                "Vérifier l'orthographe ou choisir un taxon parent. ",
                "(Erreur : ", conditionMessage(e), ")", call. = FALSE)
        NA_character_
      }
    )
    cache[[nm]] <- as.character(uuid)
    Sys.sleep(0.3)  # politesse vis-à-vis de l'API
  }
  
  # Sauvegarder le cache
  dir.create(dirname(cache_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(cache, cache_path)
  
  cache
}

#' Add PhyloPic silhouettes distributed around the functional space
#'
#' Algorithme de placement automatique :
#'   1) For each species, compute angle theta of its true position
#'      (Axis.1, Axis.2) relative to the density cloud centroid.
#'   2) Sort species by theta (i.e. by angular position around
#'      du centre).
#'   3) On force un espacement angulaire minimum entre voisins (par
#'      soft iteration, "force-directed" style) to avoid
#'      chevauchements de silhouettes.
#'   4) Each silhouette is placed on an elliptical ring defined
#'      par `ring_radius_x` et `ring_radius_y` (en fraction du demi-span).
#'
#' Result: silhouettes evenly distributed around the functional space,
#' straight segments from the inner point (true position) to the
#' corresponding outer silhouette, without crossings.
#'
#' @param p           ggplot dont on veut ajouter des silhouettes
#' @param species_list liste de listes (name, species, h)
#' @param uuid_cache  cache des UUID
#' @param coords_df   data.frame of coordinates (Axis.1, Axis.2),
#'                    rownames = species names
#' @param dim_color   couleur des points (couleur de la dimension)
#' @param sil_color, sil_alpha  apparence des silhouettes
#' @param point_size, point_alpha appearance of the dot at the true position
#' @param seg_color, seg_size, seg_alpha apparence du segment
#' @param ring_radius_x, ring_radius_y  rayons elliptiques de l'anneau
#'                    (as fraction of half-span; > 1 = beyond the panel edge)
#' @param min_angular_sep_deg  espacement angulaire minimum entre 2 silhouettes
add_phylopic_layer <- function(p, species_list, uuid_cache, coords_df,
                               dim_color = "grey20",
                               sil_color = "grey20",
                               sil_alpha = 0.85,
                               point_size = 1.6, point_alpha = 1,
                               seg_color = "grey50", seg_size = 0.2,
                               seg_alpha = 0.5,
                               ring_radius_x = 1.20,
                               ring_radius_y = 1.20,
                               min_angular_sep_deg = 25) {
  
  if (length(species_list) == 0L) return(p)
  
  # Panel bounds and cloud centroid
  x_range  <- range(coords_df$Axis.1, na.rm = TRUE)
  y_range  <- range(coords_df$Axis.2, na.rm = TRUE)
  x_span   <- diff(x_range)
  y_span   <- diff(y_range)
  cx       <- mean(x_range)
  cy       <- mean(y_range)
  half_x   <- x_span / 2
  half_y   <- y_span / 2
  
  # Étape 1 : extraction des vraies positions et calcul des angles
  spec_info <- lapply(species_list, function(sp) {
    uuid <- uuid_cache[[sp$name]]
    if (is.null(uuid) || is.na(uuid)) return(NULL)
    
    species_name <- sp$species %||% NA_character_
    if (is.na(species_name) || !species_name %in% rownames(coords_df))
      return(NULL)
    
    x_true <- coords_df[species_name, "Axis.1"]
    y_true <- coords_df[species_name, "Axis.2"]
    
    # Angle (in radians) from the centroid, normalised by span
    # to remain invariant to asymmetric scales between PC1 and PC2.
    dx_norm <- (x_true - cx) / half_x
    dy_norm <- (y_true - cy) / half_y
    theta   <- atan2(dy_norm, dx_norm)  # [-pi, pi]
    
    list(
      uuid    = uuid,
      species = species_name,
      h       = sp$h,
      x_true  = x_true,
      y_true  = y_true,
      theta   = theta
    )
  })
  spec_info <- spec_info[!vapply(spec_info, is.null, logical(1L))]
  
  if (length(spec_info) == 0L) return(p)
  
  # Étape 2 : tri par angle croissant
  spec_info <- spec_info[order(sapply(spec_info, `[[`, "theta"))]
  
  # Step 3: minimum angular spacing (greedy iterative algorithm)
  thetas <- sapply(spec_info, `[[`, "theta")
  min_sep <- min_angular_sep_deg * pi / 180
  
  # Force-directed iteration: push overly close neighbours apart
  for (iter in seq_len(100)) {
    moved <- FALSE
    for (i in seq_along(thetas)) {
      j <- if (i == length(thetas)) 1L else i + 1L
      diff_ang <- thetas[j] - thetas[i]
      if (j == 1L) diff_ang <- diff_ang + 2 * pi  # boucle circulaire
      if (diff_ang < min_sep) {
        push <- (min_sep - diff_ang) / 2
        thetas[i] <- thetas[i] - push
        thetas[j] <- thetas[j] + push
        moved <- TRUE
      }
    }
    if (!moved) break
    # Renormaliser dans [-pi, pi]
    thetas <- ((thetas + pi) %% (2 * pi)) - pi
    # Re-sort to preserve angular order
    spec_info <- spec_info[order(thetas)]
    thetas <- sort(thetas)
  }
  
  # Étape 4 : position des silhouettes sur l'anneau elliptique
  for (i in seq_along(spec_info)) {
    sp <- spec_info[[i]]
    th <- thetas[i]
    
    x_sil <- cx + ring_radius_x * half_x * cos(th)
    y_sil <- cy + ring_radius_y * half_y * sin(th)
    h_abs <- sp$h * y_span
    
    spec_info[[i]]$x_sil <- x_sil
    spec_info[[i]]$y_sil <- y_sil
    spec_info[[i]]$h_abs <- h_abs
  }
  
  # Étape 5 : construction des couches
  # 5a. Segments (couche du fond)
  seg_df <- do.call(rbind, lapply(spec_info, function(sp) {
    data.frame(x = sp$x_true, y = sp$y_true,
               xend = sp$x_sil, yend = sp$y_sil)
  }))
  if (nrow(seg_df) > 0L) {
    p <- p + geom_segment(
      data = seg_df,
      aes(x = x, y = y, xend = xend, yend = yend),
      color = seg_color, linewidth = seg_size, alpha = seg_alpha,
      inherit.aes = FALSE
    )
  }
  
  # 5b. Dots at the true position
  point_df <- do.call(rbind, lapply(spec_info, function(sp) {
    data.frame(x = sp$x_true, y = sp$y_true)
  }))
  if (nrow(point_df) > 0L) {
    p <- p + geom_point(
      data = point_df,
      aes(x = x, y = y),
      color = dim_color, size = point_size, alpha = point_alpha,
      inherit.aes = FALSE
    )
  }
  
  # 5c. Silhouettes (couche du dessus)
  for (sp in spec_info) {
    p <- p + rphylopic::geom_phylopic(
      data = data.frame(x = sp$x_sil, y = sp$y_sil),
      aes(x = x, y = y),
      uuid = sp$uuid,
      height = sp$h_abs,
      color = sil_color,
      alpha = sil_alpha,
      inherit.aes = FALSE
    )
  }
  
  p
}

# Helper : coalescing nul (R < 4.4 compatibility)
`%||%` <- function(a, b) if (is.null(a)) b else a


## 2.4 Perform Procrustes analysis between PCA spaces
run_procrustes <- function(PCA_list) {
  combNames <- combn(names(PCA_list), 2)
  combNames <- apply(combNames, 2, function(x) paste0(x[1], "_", x[2]))
  procrustes_table <- matrix(NA, ncol = 1, nrow = length(combNames),
                             dimnames = list(combNames, 'Birds'))
  
  for (j in 1:length(PCA_list)) {
    x <- PCA_list[[j]]$PCoA$vectors
    for (i in 1:length(PCA_list)) {
      if (i > j) {
        y <- PCA_list[[i]]$PCoA$vectors
        rownames(y) = rownames(x)
        prcTest <- ade4::procuste.rtest(as.data.frame(x), as.data.frame(y), nrepet = 999)
        cor_coef <- prcTest$obs
        p_val <- prcTest$pvalue
        signif <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", "ns")))
        procrustes_table[paste0(names(PCA_list)[j], "_", names(PCA_list)[i]), 1] <- 
          paste0(round(cor_coef, 3), " ", signif)
      }
    }
  }
  return(procrustes_table)
}

# ============================================================================
# 5. Fetch phylopic UUIDs (with disk cache) -----------------------------------
# ============================================================================
message("Récupération des UUID Phylopic (cache : ", PARAMS$phylopic_cache_path, ")")
phylopic_uuids <- fetch_phylopic_uuids(
  PARAMS$phylopic_species,
  cache_path = PARAMS$phylopic_cache_path
)

# Diagnostic: how many species have a valid UUID?
n_total  <- length(phylopic_uuids)
n_found  <- sum(!is.na(unlist(phylopic_uuids)))
message(sprintf("  -> %d/%d silhouettes disponibles dans le cache.",
                n_found, n_total))


# ============================================================================
# 6. Build the three main panels (with phylopic silhouettes) -----------------
# ============================================================================
multAx <- 0.25
multAx1 <- 0.29

# PARAMS$panel_titles[c(1,3,2)] matches PCA order (loco, repro, diet)
# From the previous version: keep the same indexation and associate
# each panel to its ecological dimension.
dim_order_in_PCA <- c("locomotion", "diet", "reproduction")  # ordre des panneaux

main_panels <- lapply(seq_len(3), function(d) {
  
  dim_name  <- dim_order_in_PCA[d]
  pca_obj   <- PCA[[names(PCA)[-1][d]]]
  cor_table <- format_correlation_table(pca_obj$PCoACor, shortNames)
  
  # 1. Construction du panneau de base (heatmap + contours + arrows)
  p <- PCA_plot_funspace(
    pca_obj$PCoA, cor_table, multAx1, multAx,
    paste0(PARAMS$panel_titles[d]),
    colLeg = PARAMS$dim_colors[d]
  )
  
  # 2. Ajout des silhouettes phylopic pour cette dimension
  #    coords_df must have species names in rownames so that
  #    add_phylopic_layer can retrieve the true position.
  coords_df <- data.frame(pca_obj$PCoA$vectors)
  coords_df <- coords_df[, 1:2]
  colnames(coords_df) <- c("Axis.1", "Axis.2")
  # rownames are already species names in pca_obj$PCoA$vectors
  
  p <- add_phylopic_layer(
    p,
    species_list = PARAMS$phylopic_species[[dim_name]],
    uuid_cache   = phylopic_uuids,
    coords_df    = coords_df,
    dim_color    = PARAMS$dim_colors[[dim_name]],
    sil_color    = PARAMS$silhouette_color,
    sil_alpha    = PARAMS$silhouette_alpha,
    point_size   = PARAMS$position_point_size,
    point_alpha  = PARAMS$position_point_alpha,
    seg_color    = PARAMS$segment_color,
    seg_size     = PARAMS$segment_size,
    seg_alpha    = PARAMS$segment_alpha,
    ring_radius_x       = PARAMS$ring_radius_x,
    ring_radius_y       = PARAMS$ring_radius_y,
    min_angular_sep_deg = PARAMS$min_angular_sep_deg
  )
  
  p
})

# ============================================================================
# 6. Inset: 3 x 3 Spearman correlation matrix of FUn --------------------------
# ============================================================================
fun_mat <- species_df %>%
  select(FUn_loco, FUn_diet, FUn_repro) %>%
  rename(Locomotion = FUn_loco,
         Diet       = FUn_diet,
         Reproduction = FUn_repro)

cor_mat <- cor(fun_mat, method = "spearman", use = "pairwise.complete.obs")

cor_long <- as.data.frame(as.table(cor_mat)) %>%
  rename(Var1 = Var1, Var2 = Var2, rho = Freq) %>%
  mutate(
    Var1 = factor(Var1, levels = rev(c("Locomotion", "Diet", "Reproduction"))),
    Var2 = factor(Var2, levels = c("Locomotion", "Diet", "Reproduction")),
    label = sprintf("%.2f", rho)
  )

inset_cor <- ggplot(cor_long, aes(x = Var2, y = Var1, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 2.4, color = "black") +
  scale_fill_gradient2(
    low = "#3B9AB2", mid = "white", high = "#F21A00",
    midpoint = 0, limits = c(-1, 1),
    name = expression(paste("Spearman ", rho))
  ) +
  coord_fixed() +
  theme_void(base_size = PARAMS$base_font_size - 1) +
  theme(
    plot.background    = element_rect(fill = "white", color = "grey50",
                                      linewidth = 0.3),
    axis.text.x        = element_text(angle = 30, hjust = 1,
                                      size = PARAMS$base_font_size - 2),
    axis.text.y        = element_text(size = PARAMS$base_font_size - 2),
    legend.position    = "none",
    plot.title         = element_text(size = PARAMS$base_font_size - 1,
                                      hjust = 0.5,
                                      margin = margin(b = 1)),
    plot.margin        = margin(1, 1, 1, 1)
  ) +
  ggtitle("FUn correlations")


# ============================================================================
# 7. Compose final figure (panels + inset overlay) ----------------------------
# ============================================================================
# Inset overlay top-right of the leftmost panel via patchwork::inset_element
final_fig <- ggarrange(
  main_panels[[1]],main_panels[[3]],main_panels[[2]],
  hjust = 1, align = "h", ncol = 3, nrow = 1
)
# ============================================================================
# 8. Export -------------------------------------------------------------------
# ============================================================================
out_pdf <- file.path(PARAMS$out_dir, paste0(PARAMS$fig_basename, ".pdf"))

ggsave(out_pdf, final_fig,
       width = PARAMS$width_mm, height = PARAMS$height_mm,
       units = "mm", device = cairo_pdf)


# ============================================================================
# 9. Reproducibility ----------------------------------------------------------
# ============================================================================
sessionInfo()