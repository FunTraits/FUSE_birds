#-------------------------------------------------------------------------------
# SUPP — Contribution of each trait to global functional distinctiveness
#
# For each trait, measures the effect of its removal on global
# functional distinctiveness (Di) via a normalised contribution metric:
# cDi = (Di_all - Di_minus_trait) / Di_all
# High cDi indicates a trait highly informative for functional uniqueness.
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
suppressMessages(suppressWarnings(
  source("scripts/00_START_GeneralScript.R")
))

required_pkgs <- c("ggpubr")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))


# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(

  morpho_traits = c("Tarsus.Length", "Wing.Length", "Kipps.Distance",
                    "Secondary1", "Hand.Wing.Index", "Tail.Length", "Mass"),

  lht_traits    = c("litter_or_clutch_size_n", "adult_body_mass_g", "egg_mass_g",
                    "incubation_d", "longevity_y", "fledging_age_d",
                    "litters_or_clutches_per_y", "adult_svl_cm"),

  diet_traits   = c("Diet.Inv", "Diet.Vend", "Diet.Vect", "Diet.Vfish",
                    "Diet.Vunk", "Diet.Scav", "Diet.Fruit", "Diet.Nect",
                    "Diet.Seed", "Diet.PlantO"),

  beak_traits   = c("Beak.Length_Culmen", "Beak.Length_Nares",
                    "Beak.Width", "Beak.Depth"),

  # Colours per trait group (morpho, LHT, beak, diet)
  colors        = c(rep("#99B898FF", 8), rep("#FCC893FF", 7),
                    rep("#EE6A50",   4), rep("#AAC9EDFF", 10)),

  fig_dir       = "figures"
)
dir.create(PARAMS$fig_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================================
# 2. USER INPUTS — adjust paths --------------------------------
# ============================================================================
phenoBird <- read.csv("data/processed/phenoBirdsImputedREADY.csv")
taxo      <- read.csv("data/processed/Taxo_Birds.csv")
PCA       <- readRDS("data/processed/PCA_Birds.rds")
shortNames <- read.csv("data/processed/Shortnames_Birds.csv")

all_traits <- c(PARAMS$lht_traits, PARAMS$morpho_traits,
                PARAMS$beak_traits, PARAMS$diet_traits)

stopifnot(all(all_traits %in% names(phenoBird)))


# ============================================================================
# 3. Data preparation --------------------------------------------------
# ============================================================================

phenoDiet <- na.omit(as.data.frame(
  prep.fuzzy(phenoBird[, PARAMS$diet_traits],
             col.blocks = length(PARAMS$diet_traits),
             label = "diet")
))
phenoBird <- phenoBird[rownames(phenoDiet), ]

datax <- cbind.data.frame(
  phenoBird[, PARAMS$lht_traits],
  phenoBird[, PARAMS$morpho_traits],
  phenoBird[, PARAMS$beak_traits],
  phenoDiet
)


# ============================================================================
# 4. Compute distinctiveness with/without each trait -------- -----------------------
# ============================================================================

# Reference distinctiveness (all traits)
dist_mat <- compute_dist_matrix(datax)
di_all   <- distinctiveness_global(dist_mat, di_name = "global_di")

# Distinctiveness with each trait removed in turn
for (i in seq_len(ncol(datax))) {
  dist_mat_i <- compute_dist_matrix(datax[, -i])
  di_i       <- distinctiveness_global(dist_mat_i, di_name = "global_di")
  di_all     <- cbind.data.frame(di_all, di_i[, 2])
}

# Normalised contribution per trait: cDi = (Di_all - Di_minus_i) / Di_all
cdi_all           <- (di_all[, 2] - di_all[, -c(1, 2)]) / di_all[, 2]
colnames(cdi_all) <- colnames(datax)
rownames(cdi_all) <- rownames(datax)


# ============================================================================
# 5. Figure — trait contribution -------- ----------------------------------------
# ============================================================================

cdi_plot <- boxPlotCdi_all(cdi_all, title = "", values = PARAMS$colors)

ggpubr::ggarrange(cdi_plot, hjust = 0, align = "v", ncol = 1, nrow = 1) %>%
  ggpubr::ggexport(
    filename = file.path(PARAMS$fig_dir, "Contrib.png"),
    width = 1200, height = 1200, res = 200, pointsize = 5
  )
