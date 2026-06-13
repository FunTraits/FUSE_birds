#-------------------------------------------------------------------------------
# Procrustes analyses — congruence between functional spaces
#
#   Compares the multivariate configurations of the three functional spaces
#   (locomotion, diet, reproduction) as a complement to Spearman correlations
#   on FUn and FUSE.
#
#   Outputs:
#     1. "Main text" table: pairwise Procrustes (3 pairs)
#        - Procrustes correlation (proc_corr) and m12² (sum of squared residuals)
#        - PROTEST p-value (999 permutations, Jackson 1995)
#        - conservative test stratified by taxonomic order
#     2. "Main text" table: generalised Procrustes analysis (GPA, 3 spaces)
#        - consensus score per species (useful to identify outliers)
#     3. Combined supplementary figure (4 panels):
#        (a-c) Pairwise Procrustes biplots with residual vectors per species
#        (d)   Distribution of consensus residuals per species
#        Labels: abbreviations "G. spe" (5 per panel, ggrepel)
#        Caption: abbreviation-to-full-name correspondence table
#
# Methodological reference:
#   Jackson DA (1995) PROTEST : a Procrustean randomization test of community
#   environment concordance. Écoscience 2:297-303.
#   Peres-Neto PR & Jackson DA (2001) How well do multivariate data sets match?
#   The advantages of a Procrustean superimposition approach over the Mantel
#   test. Oecologia 129:169-178.
#
# Author : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c("vegan", "dplyr", "tibble", "tidyr", "ggplot2",
                   "patchwork", "purrr", "shapes", "ggrepel")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))


# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(
  # Inputs
  pcoa_paths = c(
    locomotion   = "data/processed/PCA_Birds_L.rds",
    diet         = "data/processed/PCA_Birds_D.rds",
    reproduction = "data/processed/PCA_Birds_M.rds"  # à adapter selon vos noms
  ),
  species_table_path = "data/processed/species_table.rds",  # contient $species, $iucn, $order
  
  # Analytical parameters
  k_axes        = 2,        # nombre d'axes PCoA pour la comparaison
  n_perm        = 999,      # nombre de permutations pour PROTEST
  symmetric     = TRUE,     # symmetric Procrustes (Peres-Neto & Jackson 2001)
  
  # Outputs
  out_dir       = "data/processed",
  fig_dir       = "figures",
  fig_basename  = "FigSX_procrustes",
  width_mm      = 180,
  height_mm     = 180,
  dpi           = 600,
  
  # Apparence
  base_font_size = 8,
  dim_colors     = c(locomotion   = "#2E7D32",
                     diet         = "#C62828",
                     reproduction = "#1565C0"),
  dim_labels     = c(locomotion   = "Locomotion",
                     diet         = "Diet",
                     reproduction = "Reproduction"),
  point_size     = 0.4,
  point_alpha    = 0.4,
  arrow_alpha    = 0.25,
  arrow_size     = 0.15,
  
  # Labelling: number of "outlier" species to label per panel
  highlight_top  = 5,
  label_size     = 2.3
)
dir.create(PARAMS$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(PARAMS$fig_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(20251029)


# ============================================================================
# 2. Load and align PCoA coordinates -------- --------------------------------
# ============================================================================
load_pcoa_axes <- function(path, k = PARAMS$k_axes) {
  if (!file.exists(path))
    stop("Fichier PCoA introuvable : ", path)
  
  obj <- readRDS(path)
  
  if (!is.null(obj$PCoA) && !is.null(obj$PCoA$vectors)) {
    mat <- obj$PCoA$vectors
  } else if (!is.null(obj$vectors)) {
    mat <- obj$vectors
  } else if (is.matrix(obj)) {
    mat <- obj
  } else {
    stop("Structure inconnue pour ", path,
         ". Attendu : objet avec $PCoA$vectors ou $vectors.")
  }
  
  if (ncol(mat) < k)
    stop(sprintf("Le fichier %s contient seulement %d axes (k = %d demandé).",
                 path, ncol(mat), k))
  
  mat <- mat[, seq_len(k), drop = FALSE]
  colnames(mat) <- paste0("PC", seq_len(k))
  mat
}

coords_list <- lapply(PARAMS$pcoa_paths, load_pcoa_axes)
names(coords_list) <- names(PARAMS$pcoa_paths)

# Common species subset
common_sp <- Reduce(intersect, lapply(coords_list, rownames))
n_common  <- length(common_sp)

if (n_common < 100)
  stop(sprintf("Too few species shared across the three spaces (%d).",
               n_common))

message(sprintf("Procrustes: %d species shared across %d functional spaces.",
                n_common, length(coords_list)))

coords_list <- lapply(coords_list, function(m) m[common_sp, , drop = FALSE])

# Species metadata
species_meta <- if (file.exists(PARAMS$species_table_path)) {
  sp_tab <- readRDS(PARAMS$species_table_path)
  sp_tab[match(common_sp, sp_tab$species), , drop = FALSE]
} else {
  warning("species_table.rds introuvable : la stratification par ordre ne sera pas faite.")
  NULL
}


# ============================================================================
# 3. Helper: species name abbreviation -------- --------------------------------
# ============================================================================
#' Abbreviation in format "G. spe" from "Genus_species"
#' Exemple : "Rhea_americana" -> "R. ame"
abbreviate_species <- function(sp_name) {
  parts <- strsplit(as.character(sp_name), "_", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) < 2L) return(as.character(sp_name)[1])
    sprintf("%s. %s", substr(p[1], 1, 1), substr(p[2], 1, 3))
  }, character(1L))
}

# Shared buffer to collect abbreviation <-> full name correspondences
abbrev_lookup <- list()


# ============================================================================
# 4. Procrustes pairwise (3 paires) ------------------------------------------
# ============================================================================
pairs_to_test <- combn(names(coords_list), 2, simplify = FALSE)

run_protest <- function(pair, strata = NULL) {
  
  X <- coords_list[[pair[1]]]
  Y <- coords_list[[pair[2]]]
  
  prot <- vegan::protest(
    X = X, Y = Y,
    scores       = "sites",
    permutations = if (is.null(strata)) PARAMS$n_perm else
      vegan::how(nperm = PARAMS$n_perm, blocks = strata),
    symmetric    = PARAMS$symmetric
  )
  
  data.frame(
    dim_1        = pair[1],
    dim_2        = pair[2],
    n_species    = nrow(X),
    proc_corr    = prot$t0,
    m12_squared  = prot$ss,
    p_value      = prot$signif,
    n_perm       = prot$permutations,
    stratified   = !is.null(strata),
    stringsAsFactors = FALSE
  )
}

procrustes_table_unstrat <- do.call(
  rbind,
  lapply(pairs_to_test, run_protest, strata = NULL)
)

if (!is.null(species_meta) && "order" %in% names(species_meta)) {
  strata_vec <- factor(species_meta$order)
  procrustes_table_strat <- do.call(
    rbind,
    lapply(pairs_to_test, run_protest, strata = strata_vec)
  )
} else {
  procrustes_table_strat <- NULL
}

procrustes_main <- procrustes_table_unstrat %>%
  dplyr::transmute(
    Comparison    = paste(PARAMS$dim_labels[dim_1], "vs",
                          PARAMS$dim_labels[dim_2]),
    n_species     = n_species,
    proc_corr     = round(proc_corr, 3),
    m12_squared   = round(m12_squared, 3),
    p_value       = ifelse(p_value < 1 / (n_perm + 1),
                           sprintf("< %.3f", 1 / (n_perm + 1)),
                           sprintf("%.3f", p_value))
  )

if (!is.null(procrustes_table_strat)) {
  procrustes_main$p_value_strata_order <- ifelse(
    procrustes_table_strat$p_value < 1 / (PARAMS$n_perm + 1),
    sprintf("< %.3f", 1 / (PARAMS$n_perm + 1)),
    sprintf("%.3f", procrustes_table_strat$p_value)
  )
}

message("\n=== Procrustes pairwise (main text Table 1) ===\n")
print(procrustes_main, row.names = FALSE)

write.csv(procrustes_main,
          file.path(PARAMS$out_dir, "procrustes_pairwise_main.csv"),
          row.names = FALSE)


# ============================================================================
# 5. Conservation des objets `protest` complets pour la figure ---------------
# ============================================================================
protest_objects <- lapply(pairs_to_test, function(pair) {
  vegan::protest(coords_list[[pair[1]]], coords_list[[pair[2]]],
                 scores = "sites",
                 permutations = PARAMS$n_perm,
                 symmetric = PARAMS$symmetric)
})
names(protest_objects) <- sapply(pairs_to_test, paste, collapse = "_vs_")

saveRDS(protest_objects,
        file.path(PARAMS$out_dir, "procrustes_objects.rds"))


# ============================================================================
# 6. Generalised Procrustes analysis (GPA, 3 simultaneous spaces) -------- ----------------------
# ============================================================================
gpa_array <- array(NA_real_,
                   dim = c(n_common, PARAMS$k_axes, length(coords_list)),
                   dimnames = list(common_sp, paste0("PC", 1:PARAMS$k_axes),
                                   names(coords_list)))
for (i in seq_along(coords_list)) {
  gpa_array[, , i] <- coords_list[[i]]
}

gpa_fit <- shapes::procGPA(
  gpa_array,
  scale       = TRUE,
  reflect     = FALSE,
  proc.output = TRUE
)

# Consensus residual per species
consensus_resid <- rowSums(
  sapply(seq_len(dim(gpa_fit$rotated)[3]), function(i) {
    sqrt(rowSums((gpa_fit$rotated[, , i] - gpa_fit$mshape)^2))
  })
)

# Helper pour forcer une extraction scalaire stricte
as_scalar <- function(x) {
  x <- as.numeric(unname(x))
  if (length(x) == 0L) return(NA_real_)
  x[1]
}

# Build GPA table — robust to vectorial rho
gpa_main <- rbind(
  data.frame(Statistic = "Number of species (n)",
             Value = as_scalar(n_common)),
  data.frame(Statistic = "Number of axes per space (k)",
             Value = as_scalar(PARAMS$k_axes)),
  data.frame(Statistic = "Number of spaces (m)",
             Value = as_scalar(length(coords_list))),
  if (length(gpa_fit$rho) > 1L) {
    data.frame(
      Statistic = paste("Riemannian shape distance (rho) -",
                        names(coords_list)),
      Value = round(as.numeric(gpa_fit$rho), 4)
    )
  } else {
    data.frame(Statistic = "Riemannian shape distance (rho)",
               Value = round(as_scalar(gpa_fit$rho), 4))
  },
  data.frame(Statistic = "Sum of squares (consensus)",
             Value = round(sum(consensus_resid^2, na.rm = TRUE), 3)),
  data.frame(Statistic = "Mean per-species residual",
             Value = round(mean(consensus_resid, na.rm = TRUE), 4)),
  data.frame(Statistic = "Median per-species residual",
             Value = round(median(consensus_resid, na.rm = TRUE), 4)),
  data.frame(Statistic = "95th percentile per-species residual",
             Value = round(quantile(consensus_resid, 0.95,
                                    na.rm = TRUE, names = FALSE), 4))
)

message("\n=== Generalized Procrustes Analysis (main text Table 2) ===\n")
print(gpa_main, row.names = FALSE)

write.csv(gpa_main,
          file.path(PARAMS$out_dir, "procrustes_gpa_main.csv"),
          row.names = FALSE)

consensus_df <- tibble::tibble(
  species          = common_sp,
  consensus_resid  = consensus_resid
) %>%
  dplyr::arrange(dplyr::desc(consensus_resid))

if (!is.null(species_meta) && "iucn" %in% names(species_meta)) {
  consensus_df$iucn <- species_meta$iucn[match(consensus_df$species,
                                               species_meta$species)]
}

write.csv(consensus_df,
          file.path(PARAMS$out_dir, "procrustes_consensus_residuals.csv"),
          row.names = FALSE)

saveRDS(list(gpa_fit = gpa_fit, consensus_df = consensus_df),
        file.path(PARAMS$out_dir, "procrustes_gpa_objects.rds"))


# ============================================================================
# 7. Supplementary figure — 4 combined panels -------- -----------------------------
# ============================================================================

# ---- 7a. Biplots Procrustes pairwise --------------------------------------
make_procrustes_biplot <- function(pair_idx) {
  
  pair    <- pairs_to_test[[pair_idx]]
  prot    <- protest_objects[[pair_idx]]
  proc_r  <- prot$t0
  pv      <- prot$signif
  
  X    <- prot$X
  Yrot <- prot$Yrot
  
  df_seg <- data.frame(
    species = rownames(X),
    x_start = X[, 1],    y_start = X[, 2],
    x_end   = Yrot[, 1], y_end   = Yrot[, 2],
    resid   = sqrt((X[, 1] - Yrot[, 1])^2 + (X[, 2] - Yrot[, 2])^2)
  )
  
  df_label <- df_seg %>%
    dplyr::arrange(dplyr::desc(resid)) %>%
    dplyr::slice_head(n = PARAMS$highlight_top) %>%
    dplyr::mutate(
      species_abbr = abbreviate_species(species),
      species_full = sub("_", " ", species),
      x_mid        = (x_start + x_end) / 2,
      y_mid        = (y_start + y_end) / 2
    )
  
  # Stocker les correspondances pour la caption finale
  abbrev_lookup[[paste(pair, collapse = "_")]] <<- df_label %>%
    dplyr::select(species_abbr, species_full) %>%
    dplyr::distinct()
  
  pv_label <- if (pv < 1 / (PARAMS$n_perm + 1))
    sprintf("p < %.3f", 1 / (PARAMS$n_perm + 1)) else
      sprintf("p = %.3f", pv)
  
  ggplot() +
    geom_segment(
      data = df_seg,
      aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
      color = "grey50",
      linewidth = PARAMS$arrow_size,
      alpha = PARAMS$arrow_alpha
    ) +
    geom_point(
      data = df_seg,
      aes(x = x_start, y = y_start),
      color = PARAMS$dim_colors[[pair[1]]],
      size = PARAMS$point_size,
      alpha = PARAMS$point_alpha
    ) +
    geom_point(
      data = df_seg,
      aes(x = x_end, y = y_end),
      color = PARAMS$dim_colors[[pair[2]]],
      size = PARAMS$point_size,
      alpha = PARAMS$point_alpha
    ) +
    ggrepel::geom_text_repel(
      data = df_label,
      aes(x = x_mid, y = y_mid, label = species_abbr),
      size = PARAMS$label_size, fontface = "italic", color = "grey15",
      box.padding   = 0.4,
      point.padding = 0.15,
      segment.color = "grey50",
      segment.size  = 0.2,
      min.segment.length = 0.1,
      max.overlaps = Inf,
      seed = 42
    ) +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = sprintf("%s vs %s\nProcrustes r = %.3f\n%s",
                      PARAMS$dim_labels[[pair[1]]],
                      PARAMS$dim_labels[[pair[2]]],
                      proc_r, pv_label),
      hjust = -0.05, vjust = 1.2,
      size = 2.4, fontface = "bold", lineheight = 0.95,
      color = "grey15"
    ) +
    coord_fixed() +
    labs(x = "Procrustes axis 1", y = "Procrustes axis 2") +
    theme_minimal(base_size = PARAMS$base_font_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(color = "grey94", linewidth = 0.2),
      panel.background   = element_rect(fill = "white", color = NA),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.border       = element_rect(color = "grey60", fill = NA,
                                        linewidth = 0.3),
      axis.text          = element_text(size = PARAMS$base_font_size - 1.5,
                                        color = "grey20"),
      axis.title         = element_text(size = PARAMS$base_font_size - 0.5),
      axis.ticks         = element_line(color = "grey60", linewidth = 0.25),
      axis.ticks.length  = unit(1, "mm"),
      plot.margin        = margin(2, 2, 2, 2)
    )
}

biplot_panels <- lapply(seq_along(pairs_to_test), make_procrustes_biplot)


# ---- 7b. Distribution of consensus residuals per species -------- ------------------
iucn_colors <- c(LC = "#7BC97A", NT = "#E8E552", VU = "#F7C82E",
                 EN = "#F08C2E", CR = "#D32F2F",
                 DD = "grey70", `NA` = "grey70")

# Prepare top-N consensus labels
top_consensus_labels <- consensus_df %>%
  dplyr::slice_head(n = PARAMS$highlight_top) %>%
  dplyr::transmute(
    species_abbr   = abbreviate_species(species),
    species_full   = sub("_", " ", species),
    consensus_resid = consensus_resid
  )

# Stocker les correspondances pour la caption finale
abbrev_lookup[["consensus_top"]] <- top_consensus_labels %>%
  dplyr::select(species_abbr, species_full)

panel_consensus <- if ("iucn" %in% names(consensus_df)) {
  
  consensus_df_plot <- consensus_df %>%
    dplyr::mutate(iucn_plot = ifelse(is.na(iucn) | iucn == "DD",
                                     "DD", as.character(iucn)),
                  iucn_plot = factor(iucn_plot,
                                     levels = c("LC", "NT", "VU", "EN",
                                                "CR", "DD")))
  
  ggplot(consensus_df_plot, aes(x = consensus_resid)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 60, fill = "grey85", color = "grey60",
                   linewidth = 0.2) +
    geom_density(color = "grey20", linewidth = 0.4) +
    geom_rug(
      data = consensus_df_plot %>%
        dplyr::filter(iucn_plot %in% c("VU", "EN", "CR")),
      aes(x = consensus_resid, color = iucn_plot),
      sides = "b", length = unit(2, "mm"),
      alpha = 0.6, linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = top_consensus_labels,
      aes(x = consensus_resid, y = 0, label = species_abbr),
      size = PARAMS$label_size, fontface = "italic", color = "grey15",
      direction = "x",
      nudge_y = 0.4,
      segment.size = 0.2, segment.color = "grey50",
      box.padding   = 0.3,
      point.padding = 0.1,
      max.overlaps = Inf,
      inherit.aes = FALSE,
      seed = 42
    ) +
    scale_color_manual(
      values = iucn_colors,
      breaks = c("VU", "EN", "CR"),
      name = "IUCN status"
    ) +
    labs(
      x = "Per-species residual to GPA consensus",
      y = "Density"
    ) +
    theme_minimal(base_size = PARAMS$base_font_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(color = "grey94", linewidth = 0.2),
      panel.background   = element_rect(fill = "white", color = NA),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.border       = element_rect(color = "grey60", fill = NA,
                                        linewidth = 0.3),
      axis.text          = element_text(size = PARAMS$base_font_size - 1.5,
                                        color = "grey20"),
      axis.title         = element_text(size = PARAMS$base_font_size - 0.5),
      axis.ticks         = element_line(color = "grey60", linewidth = 0.25),
      axis.ticks.length  = unit(1, "mm"),
      legend.position    = "top",
      legend.text        = element_text(size = PARAMS$base_font_size - 1.5),
      legend.title       = element_text(size = PARAMS$base_font_size - 1,
                                        face = "bold"),
      legend.key.size    = unit(3, "mm"),
      legend.box.margin  = margin(-2, 0, -2, 0),
      plot.margin        = margin(2, 4, 4, 4)
    )
  
} else {
  ggplot(consensus_df, aes(x = consensus_resid)) +
    geom_histogram(bins = 60, fill = "grey85", color = "grey60",
                   linewidth = 0.2) +
    ggrepel::geom_text_repel(
      data = top_consensus_labels,
      aes(x = consensus_resid, y = 0, label = species_abbr),
      size = PARAMS$label_size, fontface = "italic", color = "grey15",
      direction = "x", nudge_y = 5,
      segment.size = 0.2, segment.color = "grey50",
      box.padding = 0.3, max.overlaps = Inf,
      inherit.aes = FALSE, seed = 42
    ) +
    labs(x = "Per-species residual to GPA consensus", y = "Count") +
    theme_minimal(base_size = PARAMS$base_font_size)
}


# ---- 7c. Build abbreviation -> full name table -------- --------------
abbrev_full_lookup <- dplyr::bind_rows(abbrev_lookup) %>%
  dplyr::distinct(species_abbr, species_full) %>%
  dplyr::arrange(species_abbr)

abbrev_caption <- paste0(
  "Species abbreviations: ",
  paste(sprintf("%s = %s",
                abbrev_full_lookup$species_abbr,
                abbrev_full_lookup$species_full),
        collapse = "; "),
  "."
)

# Also save as CSV for downstream use
write.csv(abbrev_full_lookup,
          file.path(PARAMS$out_dir, "procrustes_species_abbreviations.csv"),
          row.names = FALSE)


# ---- 7d. Assemblage des 4 panneaux via patchwork ---------------------------
top_row <- patchwork::wrap_plots(biplot_panels, nrow = 1, ncol = 3)

final_fig <- (top_row / panel_consensus) +
  patchwork::plot_layout(heights = c(1, 0.9)) 
  # patchwork::plot_annotation(
  #   tag_levels = "a",
  #   tag_prefix = "(", tag_suffix = ")",
  #   caption = paste(
  #     "Procrustes analyses comparing functional spaces.",
  #     "(a-c) Pairwise Procrustes biplots: each segment links a species'",
  #     "position in the first space (coloured dot at one end) to its",
  #     "position in the second space (coloured dot at the other end) after",
  #     "optimal Procrustes alignment; longer segments indicate greater",
  #     "discordance between spaces.",
  #     sprintf("The %d species with the largest Procrustes residual per panel are",
  #             PARAMS$highlight_top),
  #     "labelled with abbreviations of the form 'G. spe' (first letter of",
  #     "the genus + first three letters of the species epithet).",
  #     "(d) Distribution of per-species residuals to the consensus configuration",
  #     "from a Generalized Procrustes Analysis combining all three spaces.",
  #     "Coloured rugs indicate threatened species (VU, EN, CR).",
  #     sprintf("PROTEST p-values from %d permutations.", PARAMS$n_perm),
  #     "\n\n",
  #     abbrev_caption
  #   ),
  #   theme = theme(
  #     plot.tag = element_text(size = PARAMS$base_font_size + 1,
  #                             face = "bold"),
  #     plot.caption = element_text(size = PARAMS$base_font_size - 1.5,
  #                                 color = "grey25", hjust = 0,
  #                                 lineheight = 1.15,
  #                                 margin = margin(t = 4))
  #   )
  # )


# ============================================================================
# 8. Export -------------------------------------------------------------------
# ============================================================================
out_pdf <- file.path(PARAMS$fig_dir, paste0(PARAMS$fig_basename, ".pdf"))

ggsave(out_pdf, final_fig,
       width = PARAMS$width_mm, height = PARAMS$height_mm,
       units = "mm", device = cairo_pdf)

# ============================================================================
# 9. Reproducibility -------- --------------------------------------------------------
# ============================================================================
message("\n=== Output: tables for main text ===")
message("  - ", file.path(PARAMS$out_dir, "procrustes_pairwise_main.csv"))
message("  - ", file.path(PARAMS$out_dir, "procrustes_gpa_main.csv"))
message("  - ", file.path(PARAMS$out_dir, "procrustes_consensus_residuals.csv"))
message("  - ", file.path(PARAMS$out_dir, "procrustes_species_abbreviations.csv"))
message("\n=== Objets sauvegardés ===")
message("  - ", file.path(PARAMS$out_dir, "procrustes_objects.rds"))
message("  - ", file.path(PARAMS$out_dir, "procrustes_gpa_objects.rds"))

sessionInfo()