# Initialisation de l'environnement renv -------------------------------

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

renv::init(bare = TRUE)  # bare = TRUE : ne capture pas l'environnement actuel

# List of packages for functional ecology analysis --------

packages <- c(
  # Manipulation et import
  "tidyverse", "readxl", "janitor",
  
  # Functional data processing
  "FD", "ade4", "vegan", "cluster", "factoextra",
  
  # Mapping and geospatial
  "sf", "rnaturalearth", "terra", "ggspatial",
  
  # Phylogenetic data
  "ape", "phytools", "picante",
  
  # Visualisation
  "ggplot2", "ggpubr", "patchwork", "cowplot", "viridis",
  
  # Workflow reproductible
  "here", "knitr", "rmarkdown", "quarto"
)

# Installation et enregistrement des packages --------------------------

renv::install(packages)
renv::snapshot()  # enregistre les versions dans renv.lock

message("Environnement renv initialisé et renv.lock généré.")
