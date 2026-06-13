#-------------------------------------------------------------------------------
# SUPP — Confirmatory factor analysis (CFA) on functional traits
#
# Fits two CFA models on bird traits:
#   * 3-latent-factor model: Morphology, Life history, Diet
#   * Global model (single factor)
# Compares the two models via likelihood ratio test.
#
# Prerequisites:
#   - phenoBirdsImputedREADY.csv : complete trait table (imputation already done)
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
required_pkgs <- c(
  "tidyverse", "caret", "psych", "GGally",
  "gridGraphics", "png", "grid", "gridExtra", "RColorBrewer",
  "lavaan"
)
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))


# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(

  morpho_traits = c("Tarsus.Length", "Wing.Length", "Kipps.Distance",
                    "Secondary1", "Tail.Length", "Mass", "adult_svl_cm"),

  lht_traits    = c("litter_or_clutch_size_n", "incubation_d", "longevity_y",
                    "fledging_age_d", "litters_or_clutches_per_y"),

  diet_traits   = c("Diet.Inv", "Diet.Vend", "Diet.Vect", "Diet.Vfish",
                    "Diet.Vunk", "Diet.Scav", "Diet.Fruit", "Diet.Nect",
                    "Diet.Seed", "Diet.PlantO"),

  out_dir       = "tables"
)
dir.create(PARAMS$out_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================================
# 2. USER INPUTS — adjust paths --------------------------------
# ============================================================================
phenoBird <- read_csv("data/processed/phenoBirdsImputedREADY.csv")
rownames(phenoBird) <- phenoBird$scientificNameStd

stopifnot(
  all(PARAMS$morpho_traits %in% names(phenoBird)),
  all(PARAMS$lht_traits    %in% names(phenoBird)),
  all(PARAMS$diet_traits   %in% names(phenoBird))
)


# ============================================================================
# 3. Data preparation --------------------------------------------------
# ============================================================================
all_traits <- c(PARAMS$morpho_traits, PARAMS$lht_traits, PARAMS$diet_traits)

trait_data <- phenoBird[, all_traits] %>%
  mutate(across(everything(), ~ log10(. + 1))) %>%
  scale() %>%
  as.data.frame()


# ============================================================================
# 4. CFA model — 3 latent factors -------- -----------------------------------------
# ============================================================================
model_cfa <- paste0(
  "Morphology  =~ ", paste(PARAMS$morpho_traits, collapse = " + "), "\n",
  "LifeHistory =~ ", paste(PARAMS$lht_traits,    collapse = " + "), "\n",
  "Diet        =~ ", paste(PARAMS$diet_traits,   collapse = " + ")
)

fit_3f <- cfa(model_cfa, data = trait_data, std.lv = TRUE)
summary(fit_3f, fit.measures = TRUE, standardized = TRUE)


# ============================================================================
# 5. CFA model — global factor (benchmark) -------- ---------------------------------
# ============================================================================
model_global <- paste0(
  "Global =~ ", paste(all_traits, collapse = " + ")
)

fit_global <- cfa(model_global, data = trait_data)


# ============================================================================
# 6. Model comparison -------- --------------------------------------------------
# ============================================================================
anova(fit_3f, fit_global)   # test du rapport de vraisemblance
