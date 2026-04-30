# pkg install script

# Essential setup packages
install.packages(c("pak", "devtools", "BiocManager"))

packages <- c(
  "anndata",
  "ARTool",
  "blandr",
  "bslib",
  "car",
  "circlize",
  "cluster",
  "ComplexHeatmap",
  "corrplot",
  "data.table",
  "doMC",
  "effsize",
  "emmeans",
  "foreach",
  "fpc",
  "ggalluvial",
  "ggcorrplot",
  "ggdark",
  "ggfortify",
  "ggplot2",
  "ggpp",
  "ggpubr",
  "ggrepel",
  "ggreveal",
  "ggsci",
  "ggsignif",
  "ggtext",
  "glue",
  "gprofiler2",
  "gridExtra",
  "gtools",
  "harmony",
  "Hmisc",
  "hdf5r",
  "htmlwidgets",
  "igraph",
  "janitor",
  "kableExtra",
  "lme4",
  "MAST",
  "markdown",
  "multcomp",
  "multcompView",
  "nlme",
  "nortest",
  "nparLD",
  "openxlsx",
  "org.Hs.eg_db",
  "org.Mm.eg_db",
  "patchwork",
  "paleteer",
  "plotly",
  "reticulate",
  "Rcpp",
  "rmarkdown",
  "rsthemes",
  "rstatix",
  "rhdf5",
  "rix",
  "Seurat",
  "SeuratObject",
  "sctransform",
  "shiny",
  "shinythemes",
  "SingleR",
  "styler",
  "survival",
  "survminer",
  "svglite",
  "tidyverse",
  "WebGestaltR"
)


githubs <- c(
  "chris-mcginnis-ucsf/DoubletFinder",
  "ethanbass/ggtukey",
  "Mikata-Project/ggthemr",
  "immunogenomics/presto",
  "mhahsler/dbscan",
  "cole-trapnell-lab/monocle3",
  "mojaveazure/seurat-disk",
  "lazappi/clustree"
)

packages <- c(packages, githubs)


# Install all packages
pak::pkg_install(packages)
pak::pkg_install("data_table")
