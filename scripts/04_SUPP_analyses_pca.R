#-------------------------------------------------------------------------------
# SUPP — Visualisation of trait × PCoA axis correlations
#
# For each functional space, displays the correlations between
# original traits and PCoA axes (axes 1 and 2) as a colour-tile
# heatmap (diverging red -> blue palette).
# Also exports a summary CSV table (rounded correlations).
#
# Author  : A. Toussaint
#-------------------------------------------------------------------------------


# ============================================================================
# 0. Packages -----------------------------------------------------------------
# ============================================================================
suppressMessages(suppressWarnings(
  source("scripts/00_START_GeneralScript.R")
))

required_pkgs <- c("paletteer")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_pkgs, library, character.only = TRUE))


# ============================================================================
# 1. Parameters ---------------------------------------------------------------
# ============================================================================
PARAMS <- list(

  # Diverging red -> blue palette (201 colours)
  palette    = rev(paletteer::paletteer_c("ggthemes::Red-Blue Diverging", 201)),

  out_dir    = "tables",
  fig_dir    = "figures"
)
dir.create(PARAMS$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(PARAMS$fig_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================================
# 2. USER INPUTS — adjust paths --------------------------------
# ============================================================================
PCA        <- readRDS("data/processed/PCA_Birds.rds")
shortNames <- read.csv("data/processed/Shortnames_Birds.csv")

stopifnot(
  is.list(PCA),
  all(c("original", "short") %in% names(shortNames))
)


# ============================================================================
# 3. Result matrix setup -------- ----------------------------------
# ============================================================================

traitNames   <- names(PCA)
n_colors     <- length(PARAMS$palette)

combinations <- expand.grid(c("Axis 1", "Axis 2", "Axis 3", "Axis 4"),
                            traitNames)
result_matrix <- matrix(
  NA,
  ncol = length(traitNames) * 4,
  nrow = nrow(shortNames),
  dimnames = list(
    shortNames$original,
    sprintf("%s.%s", combinations[, 2], combinations[, 1])
  )
)


# ============================================================================
# 4. Figure — trait x axis correlation heatmap -------- -------------------------
# ============================================================================

par(mar = c(1, 2, 3, 1))

plot(0, type = "n", axes = FALSE, ann = FALSE,
     xlim = c(0, length(traitNames)),
     ylim = c(0, nrow(PCA[[1]]$PCoACor)))
graphics::box(which = "plot")

y_pos <- seq(0.5, nrow(PCA[[1]]$PCoACor), 1)

axis(2, at = y_pos, las = 1, tcl = -0.3, lwd = 0.8, labels = FALSE)
mtext(2, at = y_pos,
      text = shortNames$short[match(rownames(PCA[[1]]$PCoACor),
                                    shortNames$original)],
      las = 2, line = 0.6, cex = 0.7)

axis(3, at = seq(0.5, length(traitNames) - 0.5),
     tcl = -0.3, lwd = 0.8, labels = FALSE)
mtext(3, at = seq(0.5, length(traitNames) - 0.5),
      text = traitNames, line = 0.7, cex = 1, las = 1)

for (j in seq_along(traitNames)) {
  corr_matrix <- PCA[[j]]$PCoACor
  min1 <- min(corr_matrix[, 1], na.rm = TRUE)
  min2 <- min(corr_matrix[, 2], na.rm = TRUE)

  for (i in seq_len(nrow(corr_matrix))) {
    row_name <- rownames(corr_matrix)[i]
    if (row_name %in% rownames(PCA[[1]]$PCoACor)) {
      y_idx <- which(rownames(PCA[[1]]$PCoACor) == row_name)

      val1 <- pmin(round((corr_matrix[i, 1] - min1) * 100) + 1, n_colors)
      val2 <- pmin(round((corr_matrix[i, 2] - min2) * 100) + 1, n_colors)

      rect(j - 0.95, y_idx - 1, j - 0.5,  y_idx,
           col = PARAMS$palette[val1], border = NA)   # Axis 1
      rect(j - 0.5,  y_idx - 1, j - 0.05, y_idx,
           col = PARAMS$palette[val2], border = NA)   # Axis 2

      result_matrix[row_name,
                    grep(paste0("^", traitNames[j], "\\.Axis"),
                         colnames(result_matrix))] <- corr_matrix[i, 1:4]
    }
  }
}


# ============================================================================
# 5. Export summary table -------- ------------------------------------------
# ============================================================================

result_matrix           <- result_matrix[, c(1:6, 9, 10, 13:16)]
result_matrix           <- round(result_matrix, 3)
rownames(result_matrix) <- shortNames$short[
  match(rownames(result_matrix), shortNames$original)
]

write.csv(result_matrix,
          file = file.path(PARAMS$out_dir, "tablePCA.csv"))
