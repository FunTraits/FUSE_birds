#-------------------------------------------------------------------------------
# Figure S2. Sensitivity to the choice of k (number of neighbours for FUn).
#
# For each functional dimension (locomotion / diet / reproduction)
# and each value of k in {1, 3, 5, 10}, recomputes FUn and FUSE.
# Then assesses result stability across k via:
#
#   (a) Spearman correlations between FUn computed at different k.
#       4x4 heatmaps per dimension.
#   (b) Distribution of FUn values per k and dimension (boxplots).
#   (c) Stability of top-25 FUSE: fraction of top-25 species shared
#       between k = 5 (Pimiento 2020 reference) and each other k.
#
# Author  : A. Toussaint
#
# Prerequisites:
#   - coords     : named list of data.frames (species, PC1, PC2, ...)
#   - species_df : data.frame with columns species, iucn (LC..CR)
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
  
  # Values of k to test (Pimiento 2020 uses k = 5)
  k_values         = c(1, 3, 5, 10),
  k_reference      = 5,
  
  # PCoA axes for distances (consistent with script A1)
  pcoa_axes        = c("PC1", "PC2"),
  
  # Mapping IUCN -> GE (ordinal)
  iucn_to_GE       = c(LC = 0, NT = 1, VU = 2, EN = 3, CR = 4),
  fuse_scale       = 4,
  
  # Top-N for stability analysis (consistent with Fig. 4)
  top_n            = 25,
  
  # Apparence
  base_font_size   = 8,
  dim_order        = c("locomotion", "diet", "reproduction"),
  dim_labels       = c(locomotion = "Locomotion",
                       diet       = "Diet",
                       reproduction = "Reproduction"),
  dim_colors       = c(locomotion   = "#2E7D32",
                       diet         = "#C62828",
                       reproduction = "#1565C0"),
  
  # Output
  out_dir          = "figures",
  fig_basename     = "FigS2_FUn_k_sensitivity",
  width_mm         = 180,
  height_mm        = 220,
  dpi              = 600,
  
  # RDS file for saving numerical results
  data_out_dir     = "data/processed",
  data_out_file    = "FigS2_k_sensitivity_results.rds"
)
dir.create(PARAMS$out_dir,      showWarnings = FALSE, recursive = TRUE)
dir.create(PARAMS$data_out_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================================
# 2. USER INPUTS ------------------------------------------------------
# ============================================================================
coords     <- readRDS("data/processed/pcoa_coords.rds")
species_df <- readRDS("data/processed/species_table.rds")
species_df = na.omit(species_df[!species_df$iucn %in% c("DD","EW","EX","RE"),])

stopifnot(
  exists("coords"),
  is.list(coords),
  all(PARAMS$dim_order %in% names(coords)),
  exists("species_df"),
  all(c("species", "iucn") %in% names(species_df)),
  all(species_df$iucn %in% names(PARAMS$iucn_to_GE))
)

species_df$GE <- PARAMS$iucn_to_GE[as.character(species_df$iucn)]


# ============================================================================
# 3. Computation functions (aligned with script A1) -------- ----------------------------
# ============================================================================

compute_FUn <- function(mat, k) {
  n <- nrow(mat)
  if (n <= k) stop("n (", n, ") <= k (", k, ")")
  nn <- RANN::nn2(data = mat, query = mat, k = k + 1)
  rowMeans(nn$nn.dists[, -1, drop = FALSE])
}

compute_FSp <- function(mat) {
  centroid <- colMeans(mat, na.rm = TRUE)
  sqrt(rowSums(sweep(mat, 2, centroid, FUN = "-")^2))
}

rescale_x4 <- function(x) PARAMS$fuse_scale * scales::rescale(x, to = c(0, 1))

compute_FUSE <- function(FUn_std, FSp_std, GE) {
  log(1 + FUn_std * GE) + log(1 + FSp_std * GE)
}


# ============================================================================
# 4. Computations: FUn and FUSE per dimension x k -------- -------------------------
# ============================================================================

# Structure : results[[dim_name]] = tibble (species, FUn_k1, FUn_k3, ...,
#                                            FUSE_k1, FUSE_k3, ...)
results <- list()

for (dim_name in PARAMS$dim_order) {
  
  message("\n--- ", dim_name, " ---")
  
  cdf <- coords[[dim_name]]
  cdf <- cdf[match(species_df$species, cdf$species), , drop = FALSE]
  if (any(is.na(cdf$species)))
    stop("Espèces absentes dans coords[[", dim_name, "]].")
  
  mat <- as.matrix(cdf[, PARAMS$pcoa_axes, drop = FALSE])
  FSp_raw <- compute_FSp(mat)
  FSp_std <- rescale_x4(FSp_raw)
  
  out <- tibble(species = species_df$species,
                iucn    = species_df$iucn,
                GE      = species_df$GE)
  
  for (k in PARAMS$k_values) {
    FUn_raw <- compute_FUn(mat, k = k)
    FUn_std <- rescale_x4(FUn_raw)
    FUSE_v  <- compute_FUSE(FUn_std, FSp_std, species_df$GE)
    
    out[[paste0("FUn_k", k)]]  <- FUn_raw
    out[[paste0("FUSE_k", k)]] <- FUSE_v
    
    message(sprintf("  k=%d : FUn mean=%.4f sd=%.4f | FUSE top-3 : %s",
                    k, mean(FUn_raw), sd(FUn_raw),
                    paste(out$species[order(-FUSE_v)][1:3], collapse = ", ")))
  }
  
  results[[dim_name]] <- out
}

saveRDS(results,
        file.path(PARAMS$data_out_dir, PARAMS$data_out_file))
message("\nDonnées sauvegardées : ",
        file.path(PARAMS$data_out_dir, PARAMS$data_out_file))


# ============================================================================
# 5. (a) Spearman correlation heatmaps of FUn per dimension -------- ------------
# ============================================================================

#' FUn x FUn correlation matrix across different k values, for one dimension.
build_cor_long <- function(res_dim, dim_name) {
  fun_mat <- res_dim %>%
    select(starts_with("FUn_k"))
  cmat <- cor(fun_mat, method = "spearman", use = "pairwise.complete.obs")
  rn <- sub("FUn_k", "k=", rownames(cmat))
  colnames(cmat) <- rn
  rownames(cmat) <- rn
  as.data.frame(as.table(cmat)) %>%
    rename(k1 = Var1, k2 = Var2, rho = Freq) %>%
    mutate(dimension = dim_name,
           k1 = factor(k1, levels = rn),
           k2 = factor(k2, levels = rev(rn)))   # Y reversed for triangular display
}

cor_long <- bind_rows(lapply(PARAMS$dim_order, function(d) {
  build_cor_long(results[[d]], PARAMS$dim_labels[d])
})) %>%
  mutate(dimension = factor(dimension,
                            levels = unname(PARAMS$dim_labels[PARAMS$dim_order])))

panel_a <- ggplot(cor_long, aes(x = k1, y = k2, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", rho)),
            size = PARAMS$base_font_size / 2.8, color = "black") +
  scale_fill_gradient2(low = "#3B9AB2", mid = "white", high = "#F21A00",
                       midpoint = 0, limits = c(-1, 1),
                       name = expression(paste("Spearman ", rho))) +
  facet_wrap(~ dimension, nrow = 1) +
  coord_fixed() +
  labs(x = NULL, y = NULL,
       title = "(a) Spearman correlation between FUn computed at different k") +
  theme_minimal(base_size = PARAMS$base_font_size) +
  theme(
    plot.title          = element_text(size = PARAMS$base_font_size + 1,
                                       face = "bold", hjust = 0,
                                       margin = margin(b = 3)),
    plot.title.position = "plot",
    panel.grid          = element_blank(),
    strip.text          = element_text(face = "bold",
                                       size = PARAMS$base_font_size),
    strip.background    = element_rect(fill = "grey95", color = NA),
    axis.text           = element_text(size = PARAMS$base_font_size - 1.5,
                                       color = "grey20"),
    legend.position     = "right",
    legend.title        = element_text(size = PARAMS$base_font_size - 1),
    legend.text         = element_text(size = PARAMS$base_font_size - 1.5),
    legend.key.height   = unit(15, "mm"),
    legend.key.width    = unit(2.5, "mm"),
    plot.background     = element_rect(fill = "white", color = NA)
  )


# ============================================================================
# 6. (b) Distribution des valeurs FUn par k ----------------------------------
# ============================================================================

dist_long <- bind_rows(lapply(PARAMS$dim_order, function(d) {
  results[[d]] %>%
    select(species, starts_with("FUn_k")) %>%
    pivot_longer(starts_with("FUn_k"),
                 names_to = "k", values_to = "FUn") %>%
    mutate(k = as.integer(sub("FUn_k", "", k)),
           dimension = PARAMS$dim_labels[d])
})) %>%
  mutate(dimension = factor(dimension,
                            levels = unname(PARAMS$dim_labels[PARAMS$dim_order])),
         k = factor(k, levels = PARAMS$k_values))

panel_b <- ggplot(dist_long,
                  aes(x = k, y = FUn, fill = dimension)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3,
               linewidth = 0.3, width = 0.7) +
  scale_fill_manual(values = setNames(unname(PARAMS$dim_colors[PARAMS$dim_order]),
                                      unname(PARAMS$dim_labels[PARAMS$dim_order])),
                    guide = "none") +
  facet_wrap(~ dimension, nrow = 1, scales = "free_y") +
  labs(x = "k (number of nearest neighbours)",
       y = "FUn (raw, before standardisation)",
       title = "(b) Distribution of FUn values across k, per dimension") +
  theme_minimal(base_size = PARAMS$base_font_size) +
  theme(
    plot.title          = element_text(size = PARAMS$base_font_size + 1,
                                       face = "bold", hjust = 0,
                                       margin = margin(b = 3)),
    plot.title.position = "plot",
    panel.grid.major.x  = element_blank(),
    panel.grid.minor    = element_blank(),
    panel.grid.major.y  = element_line(color = "grey92", linewidth = 0.25),
    strip.text          = element_text(face = "bold",
                                       size = PARAMS$base_font_size),
    strip.background    = element_rect(fill = "grey95", color = NA),
    axis.text           = element_text(size = PARAMS$base_font_size - 1.5,
                                       color = "grey20"),
    axis.title          = element_text(size = PARAMS$base_font_size - 0.5),
    axis.ticks          = element_line(color = "grey60", linewidth = 0.25),
    axis.ticks.length   = unit(1, "mm"),
    panel.background    = element_rect(fill = "white", color = NA),
    plot.background     = element_rect(fill = "white", color = NA),
    panel.border        = element_rect(color = "grey60", fill = NA,
                                       linewidth = 0.3)
  )


# ============================================================================
# 7. (c) Top-25 FUSE stability between k=5 (reference) and other k -------- ---------
# ============================================================================

stability_df <- bind_rows(lapply(PARAMS$dim_order, function(d) {
  res <- results[[d]]
  ref_col <- paste0("FUSE_k", PARAMS$k_reference)
  ref_top <- res %>% arrange(desc(.data[[ref_col]])) %>%
    slice_head(n = PARAMS$top_n) %>% pull(species)
  
  out <- tibble(k = PARAMS$k_values, dimension = PARAMS$dim_labels[d],
                shared_with_ref = NA_integer_, jaccard = NA_real_)
  
  for (i in seq_along(PARAMS$k_values)) {
    k_i <- PARAMS$k_values[i]
    col_i <- paste0("FUSE_k", k_i)
    top_i <- res %>% arrange(desc(.data[[col_i]])) %>%
      slice_head(n = PARAMS$top_n) %>% pull(species)
    out$shared_with_ref[i] <- length(intersect(top_i, ref_top))
    out$jaccard[i] <- length(intersect(top_i, ref_top)) /
      length(union(top_i, ref_top))
  }
  out
})) %>%
  mutate(dimension = factor(dimension,
                            levels = unname(PARAMS$dim_labels[PARAMS$dim_order])),
         k = factor(k, levels = PARAMS$k_values))

panel_c <- ggplot(stability_df,
                  aes(x = k, y = shared_with_ref, fill = dimension)) +
  geom_col(width = 0.6, color = "grey20", linewidth = 0.2,
           position = position_dodge(width = 0.7)) +
  geom_text(aes(label = shared_with_ref),
            position = position_dodge(width = 0.7),
            vjust = -0.5,
            size = PARAMS$base_font_size / 3, color = "grey15") +
  scale_fill_manual(values = setNames(unname(PARAMS$dim_colors[PARAMS$dim_order]),
                                      unname(PARAMS$dim_labels[PARAMS$dim_order])),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     limits = c(0, PARAMS$top_n)) +
  geom_hline(yintercept = PARAMS$top_n, linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  annotate("text",
           x = length(PARAMS$k_values) + 0.4,
           y = PARAMS$top_n + 0.5,
           label = sprintf("max = %d", PARAMS$top_n),
           hjust = 1, vjust = 0,
           size = PARAMS$base_font_size / 3, color = "grey50",
           fontface = "italic") +
  labs(x = "k (number of nearest neighbours)",
       y = sprintf("Top-%d species shared with k=%d (reference)",
                   PARAMS$top_n, PARAMS$k_reference),
       title = sprintf("(c) Stability of top-%d FUSE list across k",
                       PARAMS$top_n)) +
  theme_minimal(base_size = PARAMS$base_font_size) +
  theme(
    plot.title          = element_text(size = PARAMS$base_font_size + 1,
                                       face = "bold", hjust = 0,
                                       margin = margin(b = 3)),
    plot.title.position = "plot",
    panel.grid.major.x  = element_blank(),
    panel.grid.minor    = element_blank(),
    panel.grid.major.y  = element_line(color = "grey92", linewidth = 0.25),
    axis.text           = element_text(size = PARAMS$base_font_size - 1.5,
                                       color = "grey20"),
    axis.title          = element_text(size = PARAMS$base_font_size - 0.5),
    axis.ticks          = element_line(color = "grey60", linewidth = 0.25),
    axis.ticks.length   = unit(1, "mm"),
    panel.background    = element_rect(fill = "white", color = NA),
    plot.background     = element_rect(fill = "white", color = NA),
    panel.border        = element_rect(color = "grey60", fill = NA,
                                       linewidth = 0.3),
    legend.position     = "top",
    legend.text         = element_text(size = PARAMS$base_font_size - 1),
    legend.key.size     = unit(3, "mm")
  )


# ============================================================================
# 8. Assemblage final --------------------------------------------------------
# ============================================================================
final_fig <- (panel_a / panel_b / panel_c) +
  patchwork::plot_layout(heights = c(1.3, 1.3, 1.5)) 
  # patchwork::plot_annotation(
  #   caption = sprintf(
  #     "Sensitivity analysis of functional uniqueness (FUn) and FUSE to the choice of k (number of nearest neighbours). Reference k = %d (Pimiento et al. 2020).",
  #     PARAMS$k_reference
  #   ),
  #   theme = theme(
  #     plot.caption    = element_text(size = PARAMS$base_font_size - 1.5,
  #                                    color = "grey25", hjust = 0,
  #                                    margin = margin(t = 4)),
  #     plot.background = element_rect(fill = "white", color = NA)
  #   )
  # )


# ============================================================================
# 9. Export -------------------------------------------------------------------
# ============================================================================
out_pdf <- file.path(PARAMS$out_dir, paste0(PARAMS$fig_basename, ".pdf"))

ggsave(out_pdf, final_fig,
       width = PARAMS$width_mm, height = PARAMS$height_mm,
       units = "mm", device = cairo_pdf)
message("Figure S2 written to:\n  - ", out_pdf,
        "\n  - ", out_png, "\n  - ", out_svg)


# ============================================================================
# 10. Summary table for values to report in the manuscript -------- -----
# ============================================================================
message("\n=== Stabilité du top-", PARAMS$top_n, " FUSE entre k = ",
        PARAMS$k_reference, " et autres k ===")
print(
  stability_df %>%
    arrange(dimension, k) %>%
    select(dimension, k, shared_with_ref, jaccard)
)

message("\n=== Corrélations Spearman médianes entre FUn pour différents k ===")
print(
  cor_long %>%
    filter(as.character(k1) != as.character(k2)) %>%
    group_by(dimension) %>%
    summarise(median_rho = median(rho),
              min_rho    = min(rho),
              max_rho    = max(rho),
              .groups = "drop")
)


# ============================================================================
# 11. Reproducibility -------- --------------------------------------------------------
# ============================================================================
sessionInfo()