#-------------------------------------------------------------------------------
# D. Trait selectivity GLM — which trait modalities are preferentially
#    associated with extinction risk??
#
# For each of the three functional dimensions (locomotion, diet,
# reproduction), fits two complementary models:
#   * Binomial GLMM: P(threatened | traits) with order as random effect
#                     (binary IUCN AT response: 0 = LC/NT, 1 = VU/EN/CR)
#   * Beta GLMM:     P(extinct within 100 yr | traits) (continuous
#                     IUCN 100 on ]0, 1[) — via glmmTMB
#
# Output: RDS with standardised coefficients and 95% CI + figure
#          forest plot (3 dimensions x 2 models, or 3 stacked dimensions
#          coloured by dimension).
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c(
  "dplyr", "tibble", "tidyr", "purrr", "scales",
  "lme4", "glmmTMB", "broom.mixed",
  "ggplot2", "patchwork", "forcats"
)
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# Polyfill for %||% (introduced in base R 4.4.0).
# Some patchwork / ggplot2 functions may call it without namespace
# qualifier; this ensures availability on older R versions.
if (!exists("%||%", envir = baseenv())) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  assign("%||%", `%||%`, envir = .GlobalEnv)
}


# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(
  
  # Categories considered threatened (IUCN AT)
  threatened_cats   = c("VU", "EN", "CR"),
  
  # 100-year extinction probabilities (Mooers/Davis)
  p_extinction_100yr = c(
    LC = 0.0009,
    NT = 0.0071,
    VU = 0.10,
    EN = 0.667,
    CR = 0.999
  ),
  
  # Taxonomic random effect: column name in species_df
  random_effect_col = "order",
  
  # Traits per dimension (names aligned with AVONET / EltonTraits / AMNIOTE
  # as observed in the diagnostic output).
  # IMPORTANT: any trait listed here but absent from PARAMS$no_log_transform
  # will be log10-transformed (must be strictly positive).
  traits = list(
    locomotion = c("Tarsus.Length","Wing.Length","Kipps.Distance","Secondary1",
                   "Hand.Wing.Index","Tail.Length","Mass","adult_svl_cm"),
    reproduction = c("litter_or_clutch_size_n", "egg_mass_g",
                     "incubation_d", "longevity_y", "fledging_age_d", "litters_or_clutches_per_y"),
    diet = c("Diet.Inv", "Diet.Vend","Diet.Vect","Diet.Vfish",
             "Diet.Vunk","Diet.Scav","Diet.Fruit","Diet.Nect",
             "Diet.Seed")
  ),
  
  # Columns that must NOT be log-transformed:
  # - Diet.* sont des proportions [0, 1] avec beaucoup de 0
  no_log_transform = c("Diet.Inv", "Diet.Vend", "Diet.Vect", "Diet.Vfish",
                       "Diet.Vunk", "Diet.Scav", "Diet.Fruit",
                       "Diet.Nect", "Diet.Seed", "Diet.PlantO"),
  
  # Reproducibility
  seed              = 20251028,
  
  # Output
  out_dir           = "data/processed",
  out_models_file   = "trait_selectivity_models.rds",
  out_coefs_file    = "trait_selectivity_coefficients.csv",
  
  # Figure
  fig_dir           = "figures",
  fig_basename      = "FigS1_trait_selectivity_forest",
  fig_width_mm      = 280,
  fig_height_mm     = 130,
  fig_dpi           = 600,
  base_font_size    = 8,
  dim_colors        = c(locomotion   = "#2E7D32",
                        diet         = "#C62828",
                        reproduction = "#1565C0")
)
dir.create(PARAMS$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(PARAMS$fig_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(PARAMS$seed)


# ============================================================================
# 2. USER INPUTS ------------------------------------------------------
# ============================================================================
species_df <- read.csv("data/processed/phenoBirdsImputedREADY.csv")
colnames(species_df)[1] = 'species'
colnames(species_df)[3] = 'order'
colnames(species_df)[67] = 'iucn'
species_df = species_df[!species_df$iucn %in% c("DD","EW","EX","RE"),]
species_df = species_df[!is.na(species_df$iucn),]
stopifnot(
  exists("species_df"),
  all(c("species", "iucn", PARAMS$random_effect_col) %in% names(species_df)),
  all(species_df$category %in% names(PARAMS$p_extinction_100yr))
)

# Check that all traits are present
all_traits <- unique(unlist(PARAMS$traits))
missing_traits <- setdiff(all_traits, names(species_df))
if (length(missing_traits) > 0) {
  stop("Traits absents de species_df : ",
       paste(missing_traits, collapse = ", "),
       "\nAdaptez PARAMS$traits aux noms exacts de vos colonnes.")
}

# ============================================================================
# 3. Prepare response variables -------- ---------------------------------------
# ============================================================================
species_df <- species_df %>%
  mutate(
    # binaire IUCN AT : 0 = LC/NT, 1 = VU/EN/CR
    threatened = as.integer(category %in% PARAMS$threatened_cats),
    # continuous IUCN 100: 100-year extinction probability
    p_ext_100  = unname(PARAMS$p_extinction_100yr[as.character(category)])
  )

# For beta regression, squeeze 0 and 1 into the open interval
# (transformation classique de Smithson & Verkuilen 2006)
n_obs <- nrow(species_df)
species_df$p_ext_100_open <- (species_df$p_ext_100 * (n_obs - 1) + 0.5) / n_obs

# Random effect as factor
species_df[[PARAMS$random_effect_col]] <- as.factor(
  species_df[[PARAMS$random_effect_col]]
)


# ============================================================================
# 4. Fonctions de transformation et standardisation --------------------------
# ============================================================================

#' Standardise a numeric vector (z-score), with optional log-transformation
#' optionnelle des valeurs strictement positives.
#' Retourne aussi un attribut "n_zeros" indiquant combien de valeurs
#' lost at log-transformation (useful for diagnostics).
standardise <- function(x, log_transform = TRUE, trait_name = "") {
  n_input_nonNA <- sum(!is.na(x))
  if (log_transform) {
    n_zeros <- sum(x <= 0, na.rm = TRUE)
    if (n_zeros > 0) {
      message(sprintf(
        "    %s : log10 -> %d valeurs <= 0 mises en NA (%d/%d non-NA conservées).",
        trait_name, n_zeros, n_input_nonNA - n_zeros, n_input_nonNA
      ))
      x[x <= 0] <- NA
    }
    x <- log10(x)
  }
  if (sum(!is.na(x)) == 0) {
    message(sprintf(
      "    %s : 0 valeur valide après transformation. Considérer ajouter à PARAMS$no_log_transform.",
      trait_name
    ))
    return(rep(NA_real_, length(x)))
  }
  if (sd(x, na.rm = TRUE) == 0) {
    message(sprintf("    %s : variance nulle après log10.", trait_name))
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

#' Build the standardised data.frame for a given dimension.
build_dim_data <- function(dim_name, species_df) {
  trait_cols <- PARAMS$traits[[dim_name]]
  out <- species_df %>%
    select(species, iucn, threatened, p_ext_100, p_ext_100_open,
           !!sym(PARAMS$random_effect_col),
           all_of(trait_cols))
  for (tc in trait_cols) {
    do_log <- !(tc %in% PARAMS$no_log_transform)
    out[[tc]] <- standardise(out[[tc]], log_transform = do_log,
                             trait_name = tc)
  }
  # remove rows with NA in predictors or response
  out <- out %>%
    drop_na(all_of(trait_cols), threatened, p_ext_100_open,
            !!sym(PARAMS$random_effect_col))
  out
}


# ============================================================================
# 5. Model fitting -------- ---------------------------------------------------
# ============================================================================

fit_models <- function(dim_name, species_df) {
  
  message("\n--- ", dim_name, " ---")
  
  # ---------- DIAGNOSTICS PRÉ-FILTRAGE ----------------------------------
  message(sprintf("  species_df : %d lignes, %d colonnes",
                  nrow(species_df), ncol(species_df)))
  if (nrow(species_df) == 0) {
    message("  ÉCHEC : species_df est vide. Vérifiez le chargement amont.")
    return(list(data = NULL, m_bin = NULL, m_beta = NULL))
  }
  
  trait_cols <- PARAMS$traits[[dim_name]]
  rand <- PARAMS$random_effect_col
  
  # NA par colonne sur le sous-ensemble dimensionnel
  cols_used <- c(trait_cols, "threatened", "p_ext_100_open", rand)
  cols_present <- intersect(cols_used, names(species_df))
  cols_missing <- setdiff(cols_used, cols_present)
  if (length(cols_missing) > 0) {
    message("  ATTENTION : colonnes manquantes dans species_df : ",
            paste(cols_missing, collapse = ", "))
  }
  na_before <- species_df %>%
    select(all_of(cols_present)) %>%
    summarise(across(everything(), ~ sum(is.na(.x)))) %>%
    unlist()
  message("  NA par colonne (sur ", nrow(species_df), " lignes) :")
  print(na_before)
  
  # Number of complete rows for THIS set of columns
  n_complete <- sum(complete.cases(species_df[, cols_present, drop = FALSE]))
  message(sprintf("  Lignes complètes (sans NA) : %d / %d",
                  n_complete, nrow(species_df)))
  
  if (n_complete == 0) {
    message("  ÉCHEC : aucune ligne complète. ",
            "Causes possibles : (1) jointure cassée en amont (toutes valeurs NA), ",
            "(2) noms de colonnes mal mappés, ",
            "(3) imputation incomplète sur ces traits.")
    return(list(data = NULL, m_bin = NULL, m_beta = NULL))
  }
  
  # ---------- CONSTRUCTION DU TABLEAU STANDARDISÉ -----------------------
  d <- build_dim_data(dim_name, species_df)
  message(sprintf("  n (après drop_na et standardisation) = %d", nrow(d)))
  
  if (nrow(d) == 0) {
    message("  ÉCHEC : 0 ligne après standardisation. ",
            "Probable : un trait a été entièrement converti en NA par la log-transformation ",
            "(valeurs <= 0). Vérifier PARAMS$no_log_transform.")
    return(list(data = NULL, m_bin = NULL, m_beta = NULL))
  }
  
  # Zero variance or SD = 0 (sign of broken standardisation)
  zero_var <- sapply(trait_cols, function(tc) {
    vals <- d[[tc]]
    if (all(is.na(vals))) return(TRUE)
    sd(vals, na.rm = TRUE) == 0
  })
  if (any(zero_var)) {
    message("  ATTENTION : traits à variance nulle après standardisation : ",
            paste(trait_cols[zero_var], collapse = ", "))
  }
  
  # Correlations between predictors (collinearity signal)
  if (length(trait_cols) >= 2) {
    cmat <- abs(cor(d[, trait_cols], use = "pairwise.complete.obs"))
    diag(cmat) <- 0
    if (all(is.na(cmat))) {
      message("  Corrélations toutes NA — données dégénérées.")
    } else {
      max_cor <- max(cmat, na.rm = TRUE)
      if (!is.na(max_cor) && max_cor > 0.85) {
        hi_pairs <- which(cmat > 0.85, arr.ind = TRUE)
        hi_pairs <- hi_pairs[hi_pairs[, 1] < hi_pairs[, 2], , drop = FALSE]
        message("  Forte colinéarité (|r| > 0.85) :")
        for (i in seq_len(nrow(hi_pairs))) {
          message(sprintf("    %s vs %s : r = %.3f",
                          rownames(cmat)[hi_pairs[i, 1]],
                          colnames(cmat)[hi_pairs[i, 2]],
                          cmat[hi_pairs[i, 1], hi_pairs[i, 2]]))
        }
      }
    }
  }
  
  # Effectifs threatened
  n_thr <- sum(d$threatened)
  message(sprintf("  threatened: %d / %d (%.2f%%)",
                  n_thr, nrow(d), 100 * n_thr / nrow(d)))
  if (n_thr < 10) {
    message("  ATTENTION : trop peu de threatened pour un GLMM binomial.")
  }
  # ----------------------------------------------------------------------
  
  # construction de la formule fixe
  fixed_part <- paste(trait_cols, collapse = " + ")
  random_part <- sprintf("(1 | %s)", rand)
  
  # ---- (a) GLMM binomial sur threatened (IUCN AT) ----
  fmla_bin <- as.formula(
    paste("threatened ~", fixed_part, "+", random_part)
  )
  message("  Fitting binomial GLMM (IUCN AT)...")
  m_bin <- tryCatch(
    lme4::glmer(fmla_bin, data = d, family = binomial,
                control = lme4::glmerControl(optimizer = "bobyqa",
                                             optCtrl = list(maxfun = 2e5))),
    error = function(e) {
      message("  ÉCHEC GLMM binomial (", dim_name, ") : ", e$message)
      NULL
    }
  )
  
  # ---- (b) Beta GLMM sur p_ext_100_open (IUCN 100) ----
  fmla_beta <- as.formula(
    paste("p_ext_100_open ~", fixed_part, "+", random_part)
  )
  message("  Fitting beta GLMM (IUCN 100)...")
  m_beta <- tryCatch(
    glmmTMB::glmmTMB(fmla_beta, data = d, family = glmmTMB::beta_family()),
    error = function(e) {
      message("  ÉCHEC beta GLMM (", dim_name, ") : ", e$message)
      NULL
    }
  )
  
  list(data = d, m_bin = m_bin, m_beta = m_beta)
}

results <- lapply(setNames(names(PARAMS$traits), names(PARAMS$traits)),
                  fit_models, species_df = species_df)


# ============================================================================
# 6. Harmonised coefficient extraction -------- ---------------------------------
# ============================================================================

#' Extract coefficients (estimate + 95% CI) from a fitted model.
#' Renvoie un tibble : term, estimate, conf.low, conf.high, p_value, model_type.
extract_coefs <- function(model, model_type) {
  if (is.null(model)) return(NULL)
  
  coefs <- tryCatch(
    broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE,
                      conf.method = "Wald"),
    error = function(e) NULL
  )
  if (is.null(coefs)) return(NULL)
  
  coefs %>%
    filter(term != "(Intercept)") %>%
    transmute(
      term      = term,
      estimate  = estimate,
      conf.low  = conf.low,
      conf.high = conf.high,
      p_value   = if ("p.value" %in% names(.)) p.value else NA_real_,
      model_type = model_type
    )
}

# Global coefficient table (all dimensions, both models)
coefs_all <- bind_rows(lapply(names(results), function(d) {
  bind_rows(
    extract_coefs(results[[d]]$m_bin,  "Binomial (IUCN AT)") %>%
      mutate(dimension = d),
    extract_coefs(results[[d]]$m_beta, "Beta (IUCN 100)") %>%
      mutate(dimension = d)
  )
})) %>%
  mutate(
    dimension = factor(dimension,
                       levels = names(PARAMS$traits),
                       labels = tools::toTitleCase(names(PARAMS$traits))),
    model_type = factor(model_type,
                        levels = c("Binomial (IUCN AT)", "Beta (IUCN 100)")),
    # significance: 95% CI does not cross zero
    significant = (conf.low > 0 & conf.high > 0) |
      (conf.low < 0 & conf.high < 0)
  )

# Re-order terms for display: within each dimension,
# sort by effect magnitude of the binomial IUCN AT model.
# Approach: keep 'term' as character and store a numeric rank
# par (dimension, term). C'est plus robuste que fct_reorder qui partage
# the factor scale across all dimensions breaks if one model
# fails.
#
# If the binomial model failed for a dimension, fall back to
# of the beta model, otherwise alphabetical order.
ranks_bin <- coefs_all %>%
  filter(model_type == "Binomial (IUCN AT)") %>%
  group_by(dimension) %>%
  arrange(estimate, .by_group = TRUE) %>%
  mutate(rank_in_dim = row_number()) %>%
  ungroup() %>%
  select(dimension, term, rank_in_dim)

# For (dimension, term) pairs without a binomial model, fall back to beta
ranks_beta <- coefs_all %>%
  filter(model_type == "Beta (IUCN 100)") %>%
  group_by(dimension) %>%
  arrange(estimate, .by_group = TRUE) %>%
  mutate(rank_beta = row_number()) %>%
  ungroup() %>%
  select(dimension, term, rank_beta)

coefs_all <- coefs_all %>%
  left_join(ranks_bin,  by = c("dimension", "term")) %>%
  left_join(ranks_beta, by = c("dimension", "term")) %>%
  mutate(rank_in_dim = ifelse(is.na(rank_in_dim), rank_beta, rank_in_dim)) %>%
  select(-rank_beta)

# Create a unique 'term_in_dim' factor per dimension, ordered by rank_in_dim.
# Used as the Y axis in the forest plot (facet_wrap by dimension
# en scales = "free_y").
coefs_all <- coefs_all %>%
  group_by(dimension) %>%
  mutate(term_in_dim = factor(
    term,
    levels = unique(term[order(rank_in_dim)])
  )) %>%
  ungroup()


# ============================================================================
# 7. Save models and coefficients -------- ----------------------------------
# ============================================================================
saveRDS(results,
        file.path(PARAMS$out_dir, PARAMS$out_models_file))
write.csv(coefs_all,
          file.path(PARAMS$out_dir, PARAMS$out_coefs_file),
          row.names = FALSE)
message("\nModèles et coefficients sauvegardés dans : ", PARAMS$out_dir)


# ============================================================================
# 8. Forest plot --------------------------------------------------------------
# ============================================================================
make_forest <- function(coefs_dim, dim_name, dim_color) {
  ggplot(coefs_dim,
         aes(x = estimate, y = term_in_dim,
             color = model_type,
             shape = significant)) +
    geom_vline(xintercept = 0, color = "grey60",
               linetype = "dashed", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                   height = 0, position = position_dodge(width = 0.6),
                   linewidth = 0.4) +
    geom_point(position = position_dodge(width = 0.6),
               size = 1.8, fill = "white", stroke = 0.5) +
    scale_color_manual(values = c("Binomial (IUCN AT)" = dim_color,
                                  "Beta (IUCN 100)"    = "grey45"),
                       name = NULL) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                       guide = "none") +
    labs(x = "Effect size",
         y = NULL,
         title = dim_name) +
    theme_minimal(base_size = PARAMS$base_font_size) +
    theme(
      plot.title         = element_text(size = PARAMS$base_font_size + 1,
                                        face = "bold", color = dim_color,
                                        hjust = 0),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.background   = element_rect(fill = "white", color = NA),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.border       = element_rect(color = "grey60", fill = NA,
                                        linewidth = 0.3),
      axis.text          = element_text(size = PARAMS$base_font_size - 1.5,
                                        color = "grey20"),
      axis.title.x       = element_text(size = PARAMS$base_font_size - 0.5),
      axis.ticks         = element_line(color = "grey60", linewidth = 0.25),
      axis.ticks.length  = unit(1, "mm"),
      legend.position    = "bottom",
      legend.text        = element_text(size = PARAMS$base_font_size - 1),
      legend.key.size    = unit(2, "mm"),
      plot.margin        = margin(2, 4, 2, 2)
    )
}

panels <- lapply(names(PARAMS$traits), function(d) {
  coefs_dim <- coefs_all %>%
    filter(dimension == tools::toTitleCase(d)) %>%
    mutate(term_in_dim = droplevels(term_in_dim))   # élague les niveaux des autres dim
  make_forest(
    coefs_dim,
    dim_name = tools::toTitleCase(d),
    dim_color = PARAMS$dim_colors[[d]]
  )
})

# Collected legend (single, at the bottom)
final_fig <- patchwork::wrap_plots(panels, nrow = 1) +
  theme(legend.position = "bottom")

# Annotation
# final_fig <- final_fig +
#   patchwork::plot_annotation(
#     caption = "Standardised coefficients (with 95% CI) from binomial GLMM (threatened: VU+EN+CR vs LC+NT) and beta GLMM (P(extinction in 100 yr)). \nRandom effect: order. Filled symbols: 95% CI excludes zero.",
#     theme = theme(
#       plot.caption = element_text(size = PARAMS$base_font_size - 1.5,
#                                   color = "grey25",
#                                   hjust = 0,
#                                   margin = margin(t = 3))
#     )
#   )


# ============================================================================
# 9. Export -------------------------------------------------------------------
# ============================================================================
out_pdf <- file.path(PARAMS$fig_dir, paste0(PARAMS$fig_basename, ".pdf"))
ggsave(out_pdf, final_fig,
       width = PARAMS$fig_width_mm, height = PARAMS$fig_height_mm,
       units = "mm", device = cairo_pdf)

message("Forest plot écrit dans :\n  - ", out_pdf,
        "\n  - ", out_png, "\n  - ", out_svg)


# ============================================================================
# 10. Summary tables -------- -----------------------------------------------------
# ============================================================================
message("\n=== Coefficients significatifs (IC95 hors de zéro) ===")
print(
  coefs_all %>%
    filter(significant) %>%
    arrange(dimension, model_type, desc(abs(estimate))) %>%
    select(dimension, model_type, term, estimate,
           conf.low, conf.high, p_value),
  n = Inf
)


# ============================================================================
# 11. Reproducibility -------- --------------------------------------------------------
# ============================================================================
sessionInfo()
