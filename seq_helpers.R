#' Choose appropriate number of PCA dimensions based on the amount of variability they explain
#'
#' @param srt A Seurat object after RunPCA has been performed.
#' @returns Numeric value. Number of useful PCA dimensions.
#' @export
choose_pca_dim <- function(srt) {
  ### Adapted from https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
  pct <- srt[["pca"]]@stdev / sum(srt[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  pcs <- min(co1, co2)
  return(pcs)
}

#' Run differentially expressed genes for all clusters between two conditions.
#'
#' @param srt A Seurat object with clusters and multiple conditions.
#' @param group_column String. Name of metadata column containing the conditions.
#' @param group_order String vector, length = 2. Which order the two conditions should be in. This will affect the direction of your Log2FC values.
#' @param ident_column String. Name of metadata column where your cell types or clustering assignments are found.
#' @param pvalcutoff Numeric. What p-value you consider to be significant. This defaults to 0.05.
#' @param assay String. What assay should be used? Defaults to RNA.
#' @param bulk Boolean. Determines if pseudobulking should be applied before analysis. Defaults to F.
#' @param sample_column Name of metadata column containing sample assignments for each cell. Only required for pseudobulking.
#' @param test.use Passed to FindMarkers. Defaults to Wilcoxon for SC and DESeq2 for bulk.
#' @param ... Additional arguments to pass to Seurat's FindMarkers
#' @returns Dataframe of differentially expressed genes.
#' @import Seurat
#' @import magrittr
#' @export
deggerator <- function(srt,
                       group_column,
                       group_order,
                       ident_column,
                       pvalcutoff = 0.05,
                       assay = "RNA",
                       bulk = F,
                       sample_column = NULL,
                       test.use = NULL,
                       ...) {
  if (is.null(test.use)) {
    if (bulk) {
      test.use <- "DESeq2"
    }
    if (!bulk) {
      test.use <- "MAST"
    }
  }

  # Pseudobulk if needed
  if (bulk) {
    srt <- AggregateExpression(
      srt,
      return.seurat = TRUE,
      group.by = c(group_column, sample_column, ident_column)
    )
  }

  # Create new cell by condition metadata column and pull out unique values
  srt$cell_by_condition <- paste0(srt[[ident_column]][, 1], srt[[group_column]][, 1])
  Idents(srt) <- "cell_by_condition"
  # cell_type = unique(srt$cell_by_condition) %>% sort()
  cell_type <- unique(srt[[ident_column]][, 1]) %>% sort()

  # Run DEGs for each unique pair of cluster and condition.
  deg <- list()
  for (type in cell_type) {
    message(paste0("Analyzing ", type))
    g.one <- paste0(type, group_order[1])
    g.two <- paste0(type, group_order[2])

    if ((length(srt@active.ident[srt@active.ident == g.one]) < 3) |
      (length(srt@active.ident[srt@active.ident == g.two]) < 3)) {
      temp_deg <- list(data.frame("Not Enough Cells for Comparison"))
      deg <- c(deg, temp_deg)
      next
    }

    # if (bulk){
    temp_deg <- FindMarkers(
      object = srt,
      ident.1 = g.two,
      ident.2 = g.one,
      assay = assay,
      test.use = test.use,
      ...
    )

    # Combine and filter values

    # Error catching
    if (!"p_val_adj" %in% colnames(temp_deg)) {
      cell_type <- cell_type[!cell_type == type]
      next
    }

    temp_deg <- cbind(gene = rownames(temp_deg), temp_deg)
    temp_deg <- subset(temp_deg, p_val_adj < pvalcutoff)

    temp_deg$FC <- ifelse(
      sign(temp_deg$avg_log2FC) == "-1",
      yes = (1 / (2^temp_deg$avg_log2FC)) * -1,
      no = 2^temp_deg$avg_log2FC
    )

    # temp_deg$FC = 2 ^ temp_deg$avg_log2FC
    temp_deg <- list(temp_deg)
    deg <- c(deg, temp_deg)
  }
  names(deg) <- cell_type
  return(deg)
}

#' #' @param deg List of differentially expressed gene dataframes. Usually the output of the deggerator.
#' #' @param genes String vector. Optional list of genes to collect into a genes of interest sheet at beginning of the list?
#' #' @param pull_genes Boolean. Pull current genes of interest list from GitHub?
#' #' @param gene_pattern String (Regex). Regex matching of gene patterns to collect in the genes of interest.
#' #' @param sort Boolean. Sort the names of the deg list items alphabetically?
#' #' @param stats Boolean. Create a DEG stats item at beginning of list?
#' #' @param rm_mito Boolean. Remove mitochondrial genes?
#' #' @param sort_by_FC Boolean. Sort genes by AvgLog2FC?
#' #' @param gene_name Boolean. Search database to get gene gene name?
#' #' @param db_name String. Name of species gene annotation database from BioConductor. Defaults to "org.Mm.eg.db".
#' #' @param order Sort genes by AvgLog2FC?
#' #' @param ... Arguments to be passed to the select function for the gene annotation database.
#' #' @returns A list of diffentially expressed gene dataframes that have been processed.
#' #' @import AnnotationDbi
#' #' @import dplyr
#' #' @export
#' deg_processor = function(deg,
#'                          genes = NULL,
#'                          pull_genes = T,
#'                          gene_pattern = NULL,
#'                          sort_by_FC = T,
#'                          stats = T,
#'                          rm_mito = F,
#'                          order = T,
#'                          gene_name = T,
#'                          db_name = 'org.Mm.eg.db',
#'                          ...) {
#'   library(dplyr)
#'   # Initial cleanup
#'   for (name in names(deg)) {
#'     if (is.null(deg[[name]]) || nrow(deg[[name]]) == 0) {
#'       deg[[name]] = NULL
#'     }
#'   }
#'
#'   if (order) {
#'     deg = deg[order(names(deg))]
#'   }
#'
#'   if (gene_name) {
#'     data_names = names(deg)
#'     # Check if BiocManager is installed, install if not
#'     if (!requireNamespace("BiocManager", quietly = TRUE)) {
#'       install.packages("BiocManager")
#'     }
#'
#'     # Check if the specified database is installed
#'     if (!requireNamespace(db_name, quietly = TRUE)) {
#'       # Install the database if not installed
#'       BiocManager::install(db_name)
#'     }
#'
#'     # Load the database
#'     library(db_name, character.only = TRUE)
#'
#'     for (name in names(deg)) {
#'       annotations = NULL # try to catch some errors
#'       cluster = deg[[name]]
#'       annotations = tryCatch({
#'         suppressMessages(
#'           select(
#'             org.Mm.eg.db,
#'             keys = cluster$gene,
#'             columns = "GENENAME",
#'             keytype = "SYMBOL"
#'           )
#'         )
#'       }, error = function(e) {
#'         return(NULL)
#'       })
#'
#'       if (is.null(annotations) || nrow(annotations) == 0) {
#'         next
#'       }
#'
#'       deg[[name]] = cbind(deg[[name]], annotations$GENENAME)
#'       colnames(deg[[name]])[colnames(deg[[name]]) == 'annotations$GENENAME'] = 'gene_name'
#'     }
#'
#'   }
#'
#'   if (is.null(genes) & pull_genes) {
#'     goi = read.csv(
#'       url(
#'         'https://raw.githubusercontent.com/JordanWean/jbw_utils/refs/heads/main/genesofinterest.csv'
#'       )
#'     )
#'     genes = goi$gene
#'   }
#'
#'
#'   if (!is.null(genes) | !is.null(gene_pattern)) {
#'     for (name in names(deg)) {
#'       cluster = deg[[name]]
#'       if (length(cluster) == 1 || nrow(cluster) == 0) {
#'         next
#'       }
#'
#'       if (!is.null(genes)) {
#'         x = cluster[cluster$gene %in% genes, ]
#'       }
#'       if (!is.null(gene_pattern)) {
#'         y = cluster[grepl(gene_pattern, cluster$gene), ]
#'       }
#'
#'       if (exists("x") & exists("y")) {
#'         z = rbind(x, y)
#'       } else if (exists("x")) {
#'         z = x
#'       } else if (exists("y")) {
#'         z = y
#'       }
#'
#'       if (name == names(deg)[length(names(deg))] &
#'           nrow(z) == 0 & !exists("genes_of_interest")) {
#'         break
#'       }
#'
#'       if (!exists('z') | nrow(z) == 0) {
#'         next
#'       }
#'       z$cluster = name
#'
#'       if (!exists("genes_of_interest")) {
#'         genes_of_interest = data.frame(z)
#'       } else {
#'         # genes_of_interest = rbind(genes_of_interest, z)
#'         genes_of_interest = bind_rows(genes_of_interest, z)
#'       }
#'     }
#'     if (exists('genes_of_interest')) {
#'       if (sort_by_FC) {
#'         genes_of_interest = genes_of_interest[order(abs(genes_of_interest$avg_log2FC), decreasing = T), ]
#'         genes_of_interest = genes_of_interest[order(genes_of_interest$pct.1, decreasing = T), ]
#'       }
#'
#'       genes_of_interest = list(genes_of_interest)
#'       names(genes_of_interest) = 'Genes of Interest'
#'       deg = append(deg, genes_of_interest, after = 0)
#'     }
#'   }
#'
#'   if (rm_mito) {
#'     deg = lapply(
#'       deg,
#'       FUN = function(x) {
#'         x = x[!grepl('^mt-|^MT-', x$gene), ]
#'       }
#'     )
#'   }
#'
#'
#'   if (sort_by_FC) {
#'     deg = lapply(
#'       deg,
#'       FUN = function(x) {
#'         if (!length(x) == 0 & length(colnames(x)) > 1) {
#'           if (is.numeric(x$avg_log2FC)) {
#'             x = x[order(abs(x$avg_log2FC), decreasing = T), ]
#'           }
#'         }
#'       }
#'     )
#'   }
#'
#'   if (stats) {
#'     colnames = c(
#'       'group',
#'       'nDEG',
#'       'nPositive',
#'       'maxlog2FC',
#'       'avgtop5log2FC',
#'       'nNegative',
#'       'minlog2FC',
#'       'avgbottom5log2FC'
#'     )
#'     stats = data.frame(matrix(ncol = length(colnames), nrow = 0))
#'     for (group in names(deg)) {
#'       x = deg[[group]]
#'       if (is.numeric(x$avg_log2FC)) {
#'         x$avg_log2FC = signif(x$avg_log2FC, digits = 5)
#'       }
#'       if (is.null(x)) {
#'         next
#'       }
#'       if (nrow(x) > 0) {
#'         y = sort(x$avg_log2FC, decreasing = T)
#'         minFC = min(x$avg_log2FC)
#'         maxFC = max(x$avg_log2FC)
#'         nPos = sum(x$avg_log2FC > 0)
#'         nNeg = sum(x$avg_log2FC < 0)
#'       }
#'       if (nrow(x) > 5) {
#'         top = mean(y[1:5])
#'         bottom = mean(tail(y, 5))
#'       }
#'       else {
#'         # There surely is a much better way to do this cleanup
#'         if (!exists('nPos')) {
#'           nPos = ''
#'         }
#'         if (!exists('maxFC')) {
#'           maxFC = ''
#'         }
#'         if (!exists('nNeg')) {
#'           nNeg = ''
#'         }
#'         if (!exists('minFC')) {
#'           minFC = ''
#'         }
#'         top = ''
#'         bottom = ''
#'       }
#'       s = c(group,
#'             as.character(nrow(x)),
#'             nPos,
#'             maxFC,
#'             top,
#'             nNeg,
#'             minFC,
#'             bottom)
#'       stats = rbind(stats, s)
#'
#'       objnames = c('group', 'nPos', 'maxFC', 'top', 'nNeg', 'bottom', 'minFC')
#'       for (obj in objnames) {
#'         if (exists(obj)) {
#'           rm(list = obj)
#'         }
#'       }
#'     }
#'     colnames(stats) = colnames
#'     deg = c(GroupStats = list(stats), deg)
#'   }
#'
#'   # Ending cleanup
#'   for (name in names(deg)) {
#'     if (is.null(deg[[name]]) || nrow(deg[[name]]) == 0) {
#'       deg[[name]] = NULL
#'     }
#'   }
#'
#'   return(deg)
#' }

#### Not my function, just added to a  useful one.
RNAseq_data_import <- function(h5,
                               min.cells = 3,
                               min.features = 200,
                               condition.layer) {
  if (length(h5) == 1) {
    mtx <- Read10X_h5(h5)
    RNA_srt <- CreateSeuratObject(mtx, min.cells = min.cells, min.features = min.features)
    RNA_srt$condition <- strsplit(h5[1], "/")[[1]][length(strsplit(h5[1], "/")[[1]]) - condition.layer]
    return(RNA_srt)
  } else {
    # Initialise list to store seurat objects
    seurat_list <- list()
    # For each h5 file, create seurat object and add "condition" metadata then store in list
    for (i in 1:length(h5)) {
      mtx <- Read10X_h5(paste0(h5[i]))
      RNA_srt <- CreateSeuratObject(mtx, min.cells = min.cells, min.features = min.features)
      RNA_srt$condition <- strsplit(h5[i], "/")[[1]][length(strsplit(h5[i], "/")[[1]]) - condition.layer]
      seurat_list[[i]] <- RNA_srt
    }
    # Merge seurat objects
    RNA_merged <- merge(
      x = seurat_list[[1]],
      y = c(seurat_list[-1]),
      add.cell.ids = sapply(h5, function(x) {
        strsplit(x, "/")[[1]][1]
      })
    )

    return(RNA_merged)
  }
}

#' Store sample ID information as a metadata column.
#' @param srt Seurat object.
#' @import Seurat
store_sample_IDs <- function(srt) {
  sampleIDs <- as.data.frame(srt@assays$RNA@cells)
  sampleIDs[sampleIDs == FALSE] <- ""
  sampleIDs[sampleIDs == TRUE] <- colnames(sampleIDs)[which(sampleIDs == TRUE, arr.ind = TRUE)[, "col"]]
  sampleIDs <- unite(
    col = "comb",
    data = sampleIDs,
    colnames(sampleIDs),
    sep = ""
  )
  srt$sampleIDs <- sampleIDs
  return(srt)
}

#' Check for separation of genes that should not cluster together. Typically for glutamatergic/GABAergic neurons.
#'
#' @param srt Seurat object.
#' @param pattern Regex pattern on which column names to analyze.
#' @param features Which genes to test.
#' @import Seurat
#' @import magrittr
ratiochecker <- function(srt, pattern, features) {
  for (col in names(srt@meta.data)[grepl(pattern, names(srt@meta.data))]) {
    Idents(srt) <- col
    print(paste0("Processing ", col))
    x <- AverageExpression(srt, features = features)
    x <- x$RNA %>%
      as.matrix() %>%
      t()
    x <- ifelse(x[, 1] > x[, 2], yes = x[, 1] / x[, 2], no = x[, 2] / x[, 1])
    stats <- c(col, mean(x), median(x), mean(x) / median(x))

    if (!exists("resframe")) {
      resframe <- as.data.frame(t(stats))
    } else {
      resframe <- rbind(resframe, stats)
    }
  }
  names(resframe) <- c("resolution", "mean", "median", "mean_median_ratio")
  resframe$median <- resframe$median %>% as.numeric()
  resframe$mean <- resframe$mean %>% as.numeric()
  resframe$mean_median_ratio <- resframe$mean_median_ratio %>% as.numeric()
  resframe$medianntile <- ifelse(
    resframe$median > quantile(resframe$median, probs = 0.90),
    yes = TRUE,
    no = FALSE
  )
  return(resframe)
}

#' @param deg Either a path to an excel workbook containing sheets with DEGs in them, a list object where each item is a table of degs, or a data.frame full of degs.
#' @param genes String vector. Optional list of genes to collect into a genes of interest sheet at beginning of the list?
#' @param pull_genes Boolean. Pull current genes of interest list from GitHub?
#' @param gene_pattern String (Regex). Regex matching of gene patterns to collect in the genes of interest.
#' @param sort Boolean. Sort the names of the deg list items alphabetically?
#' @param stats Boolean. Create a DEG stats item at beginning of list?
#' @param rm_mito Boolean. Remove mitochondrial genes?
#' @param sort_by_FC Boolean. Sort genes by AvgLog2FC?
#' @param gene_name Boolean. Search database to get gene names?
#' @param db_name String. Name of species gene annotation database from BioConductor. Defaults to "org.Mm.eg.db".
#' @param order Sort genes by AvgLog2FC?
#' @param pval_thresh Maximum acceptable p value.
#' @param pct.1_thresh  Minimum proportion of cells for pct.1.
#' @param pct.2_thresh  Minimum proportion of cells for pct.2.
#' @param pct.or Minumum proportion of cells in at least one condition.
#' @param log2fc_thresh Mimimum acceptable log2 fold change.
#' @param fc_thresh Minimum acceptable fold change.
#' @param filter_only Filter the input only, no gene of interest or group stats.
#' @param ... Arguments to be passed to the select function for the gene annotation database.
#' @returns A list of diffentially expressed gene dataframes that have been processed.
#' @import AnnotationDbi
#' @import dplyr
#' @import BiocManager
#' @export
deg_processor <- function(deg,
                          genes = NULL,
                          pull_genes = T,
                          gene_pattern = NULL,
                          sort_by_FC = T,
                          stats = T,
                          rm_mito = F,
                          order = T,
                          gene_name = T,
                          db_name = "org.Mm.eg.db",
                          pval_thresh = NULL,
                          pct.1_thresh = NULL,
                          pct.2_thresh = NULL,
                          pct.or_thresh = NULL,
                          log2fc_thresh = NULL,
                          fc_thresh = NULL,
                          filter_only = F,
                          ...) {
  # DEG input type checking
  deg_type <- case_when(
    is.data.frame(deg) ~ "df",
    is_list(deg) ~ "list",
    is.character(deg) ~ "path",
    .default = NA
  )

  unlist_deg <- F

  if (is.na(deg_type)) {
    message("Please enter a valid DEG object.")
    stop()
  }

  if (deg_type == "path") {
    ext <- sub(".*\\.", "", deg)

    if (ext == "xlsx") {
      path2workbook <- deg
      sheet_names <- getSheetNames(path2workbook)
      workbook <- lapply(sheet_names, function(sheet) {
        read.xlsx(path2workbook, sheet = sheet)
      })
      names(workbook) <- sheet_names
      deg <- workbook
      deg_type <- "list"
    }

    if (ext == "csv") {
      deg <- read.csv(deg)
      deg_type <- "df"
    }
  }

  if (deg_type == "df") {
    deg <- list(df = deg)
    deg_type <- "list"
    unlist_deg <- T
  }

  # does this need more conditional statements to make sure it doesn't happen at the wrong times?
  if ("GroupStats" %in% names(deg)) {
    deg$GroupStats <- NULL
  }
  if ("Genes of Interest" %in% names(deg)) {
    deg$`Genes of Interest` <- NULL
  }

  library(dplyr)

  # Initial cleanup
  for (name in names(deg)) {
    if (is.null(deg[[name]]) || nrow(deg[[name]]) == 0) {
      deg[[name]] <- NULL
    }
  }

  if (order) {
    deg <- deg[order(names(deg))]
  }

  if (gene_name) {
    data_names <- names(deg)
    # Check if BiocManager is installed, install if not
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }

    # Check if the specified database is installed
    if (!requireNamespace(db_name, quietly = TRUE)) {
      # Install the database if not installed
      BiocManager::install(db_name)
    }

    # Load the database
    library(db_name, character.only = TRUE)

    for (name in names(deg)) {
      annotations <- NULL # try to catch some errors
      cluster <- deg[[name]]
      annotations <- tryCatch(
        {
          suppressMessages(
            select(
              org.Mm.eg.db,
              keys = cluster$gene,
              columns = "GENENAME",
              keytype = "SYMBOL"
            )
          )
        },
        error = function(e) {
          return(NULL)
        }
      )

      if (is.null(annotations) || nrow(annotations) == 0) {
        next
      }

      deg[[name]] <- cbind(deg[[name]], annotations$GENENAME)
      colnames(deg[[name]])[colnames(deg[[name]]) == "annotations$GENENAME"] <- "gene_name"
    }
  }

  # Filtering

  filt <- function(x,
                   pval_thresh,
                   pct.1_thresh,
                   pct.2_thresh,
                   pct.or_thresh,
                   log2fc_thresh,
                   fc_thresh) {
    if (!is.null(pval_thresh)) {
      x <- x[x$p_val_adj < pval_thresh, ]
    }
    if (!is.null(pct.1_thresh)) {
      x <- x[x$pct.1 > pct.1_thresh, ]
    }
    if (!is.null(pct.2_thresh)) {
      x <- x[x$pct.2 > pct.2_thresh, ]
    }
    if (!is.null(pct.or_thresh)) {
      x <- x[x$pct.1 > pct.or_thresh | x$pct.2 > pct.or_thresh, ]
    }
    if (!is.null(log2fc_thresh)) {
      x <- x[abs(x$log2fc) > log2fc_thresh, ]
    }
    if (!is.null(fc_thresh)) {
      x <- x[abs(x$FC) > fc_thresh, ]
    }
    return(x)
  }

  for (name in names(deg)) {
    tempdeg <- deg[[name]]
    tempdeg <- filt(
      x = deg[[name]],
      pval_thresh = pval_thresh,
      pct.1_thresh = pct.1_thresh,
      pct.2_thresh = pct.2_thresh,
      pct.or_thresh = pct.or_thresh,
      log2fc_thresh = log2fc_thresh,
      fc_thresh = fc_thresh
    )
    if (nrow(tempdeg) == 0) {
      message("Filtering removed all genes from ", name)
    }
    deg[[name]] <- tempdeg
  }

  if (filter_only) {
    if (unlist_deg) {
      deg <- deg$df
    }
    return(deg)
  }

  # GOI and group stats.

  if (is.null(genes) & pull_genes) {
    goi <- read.csv(
      url(
        "https://raw.githubusercontent.com/JordanWean/jbw_utils/refs/heads/main/genesofinterest.csv"
      )
    )
    genes <- goi$gene
  }


  if (!is.null(genes) | !is.null(gene_pattern)) {
    for (name in names(deg)) {
      cluster <- deg[[name]]
      if (length(cluster) == 1 || nrow(cluster) == 0) {
        next
      }

      if (!is.null(genes)) {
        x <- cluster[cluster$gene %in% genes, ]
      }
      if (!is.null(gene_pattern)) {
        y <- cluster[grepl(gene_pattern, cluster$gene), ]
      }

      if (exists("x") & exists("y")) {
        z <- rbind(x, y)
      } else if (exists("x")) {
        z <- x
      } else if (exists("y")) {
        z <- y
      }

      if (name == names(deg)[length(names(deg))] &
        nrow(z) == 0 & !exists("genes_of_interest")) {
        break
      }

      if (!exists("z") | nrow(z) == 0) {
        next
      }
      z$cluster <- name

      if (!exists("genes_of_interest")) {
        genes_of_interest <- data.frame(z)
      } else {
        # genes_of_interest = rbind(genes_of_interest, z)
        genes_of_interest <- bind_rows(genes_of_interest, z)
      }
    }
    if (exists("genes_of_interest")) {
      if (sort_by_FC) {
        genes_of_interest <- genes_of_interest[order(abs(genes_of_interest$avg_log2FC), decreasing = T), ]
        genes_of_interest <- genes_of_interest[order(genes_of_interest$pct.1, decreasing = T), ]
      }

      genes_of_interest <- list(genes_of_interest)
      names(genes_of_interest) <- "Genes of Interest"
      deg <- append(deg, genes_of_interest, after = 0)
    }
  }

  if (rm_mito) {
    deg <- lapply(
      deg,
      FUN = function(x) {
        x <- x[!grepl("^mt-|^MT-", x$gene), ]
      }
    )
  }

  if (sort_by_FC) {
    deg <- lapply(
      deg,
      FUN = function(x) {
        if (!length(x) == 0 & length(colnames(x)) > 1) {
          if (is.numeric(x$avg_log2FC)) {
            x <- x[order(abs(x$avg_log2FC), decreasing = T), ]
          }
        }
      }
    )
  }

  if (stats) {
    colnames <- c(
      "group",
      "nDEG",
      "nPositive",
      "maxlog2FC",
      "avgtop5log2FC",
      "nNegative",
      "minlog2FC",
      "avgbottom5log2FC"
    )
    stats <- data.frame(matrix(ncol = length(colnames), nrow = 0))
    for (group in names(deg)) {
      x <- deg[[group]]
      if (is.numeric(x$avg_log2FC)) {
        x$avg_log2FC <- signif(x$avg_log2FC, digits = 5)
      }
      if (is.null(x)) {
        next
      }
      if (nrow(x) > 0) {
        y <- sort(x$avg_log2FC, decreasing = T)
        minFC <- min(x$avg_log2FC)
        maxFC <- max(x$avg_log2FC)
        nPos <- sum(x$avg_log2FC > 0)
        nNeg <- sum(x$avg_log2FC < 0)
      }
      if (nrow(x) > 5) {
        top <- mean(y[1:5])
        bottom <- mean(tail(y, 5))
      } else {
        # There surely is a much better way to do this cleanup
        if (!exists("nPos")) {
          nPos <- ""
        }
        if (!exists("maxFC")) {
          maxFC <- ""
        }
        if (!exists("nNeg")) {
          nNeg <- ""
        }
        if (!exists("minFC")) {
          minFC <- ""
        }
        top <- ""
        bottom <- ""
      }
      s <- c(
        group,
        as.character(nrow(x)),
        nPos,
        maxFC,
        top,
        nNeg,
        minFC,
        bottom
      )
      stats <- rbind(stats, s)

      objnames <- c("group", "nPos", "maxFC", "top", "nNeg", "bottom", "minFC")
      for (obj in objnames) {
        if (exists(obj)) {
          rm(list = obj)
        }
      }
    }
    colnames(stats) <- colnames
    deg <- c(GroupStats = list(stats), deg)
  }

  # Ending cleanup
  for (name in names(deg)) {
    if (is.null(deg[[name]]) || nrow(deg[[name]]) == 0) {
      deg[[name]] <- NULL
    }
  }

  if (unlist_deg) {
    deg <- deg$df
  }

  return(deg)
}


log2fc_xl_highlighter <- function(input_path, xlsx_files = NULL) {
  library(openxlsx)

  style_blue <- createStyle(fgFill = "#0D6EFD")
  style_green <- createStyle(fgFill = "#33b249")
  style_orange <- createStyle(fgFill = "#ffbe30")
  style_red <- createStyle(fgFill = "#ff3838")


  if (is.null(xlsx_files)) {
    xlsx_files <- list.files(input_path, pattern = "\\.xlsx$", recursive = TRUE, full.names = TRUE)
  }

  for (file in xlsx_files) {
    wb <- loadWorkbook(file)

    modified <- FALSE
    for (sheet in sheets(wb)) {
      df <- read.xlsx(file, sheet = sheet)

      if ("avg_log2FC" %in% names(df)) {
        col_num <- which(names(df) == "avg_log2FC")
        for (row in 1:nrow(df)) {
          val <- abs(df[[col_num]][row])

          # Only format if numeric
          if (!is.na(val) && is.numeric(val)) {
            # Excel coordinates: remember +1 for header
            cell_row <- row + 1
            cell_col <- col_num
            if (val >= 1 && val < 2) {
              addStyle(wb, sheet, style_blue, rows = cell_row, cols = cell_col, gridExpand = FALSE, stack = TRUE)
            } else if (val >= 2 && val < 3) {
              addStyle(wb, sheet, style_green, rows = cell_row, cols = cell_col, gridExpand = FALSE, stack = TRUE)
            } else if (val >= 3 && val < 4) {
              addStyle(wb, sheet, style_orange, rows = cell_row, cols = cell_col, gridExpand = FALSE, stack = TRUE)
            } else if (val >= 4) {
              addStyle(wb, sheet, style_red, rows = cell_row, cols = cell_col, gridExpand = FALSE, stack = TRUE)
            }
          }
        }
        modified <- TRUE
        cat(sprintf("Manually formatted column 'avg_log2FC' in sheet '%s' of file '%s'.\n", sheet, file))
      }
    }
    if (modified) {
      saveWorkbook(wb, file, overwrite = TRUE)
      cat(sprintf("File '%s' updated and saved.\n", file))
    }
  }
}


gzip_all_files <- function(directory) {
  files <- list.files(directory, recursive = TRUE, full.names = TRUE)
  for (f in files) {
    if (!grepl("\\.gz$", f)) {
      system(paste("gzip", shQuote(f)))
    }
  }
}


ungzip_all_files <- function(directory) {
  pkg <- "R.utils"
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
  gz_files <- list.files(directory, pattern = "\\.gz$", recursive = TRUE, full.names = TRUE)
  for (fz in gz_files) {
    message("Unzipping: ", fz)
    # Unzip in-place, remove gz file after
    gunzip(fz, overwrite = TRUE, remove = TRUE)
  }
  invisible(NULL)
}



check_3d_umap <- function(srt,
                          labels,
                          embeddings = NULL,
                          dims = NULL,
                          pt.size = 2) {
  if (is.null(embeddings & !is.null(dims))) {
    srt <- Seurat::RunUMAP(srt, n.components = 3, reduction.name = "umap_3d", dims = dims)
    df_3d <- srt@reductions$umap_3d@cell.embeddings %>% as.data.frame()
  } else if (is.null(embeddings)) {
    df_3d <- as.data.frame(embeddings)
  }
  
  colnames(df_3d) <- c("x", "y", "z")
  df_3d$label <- labels
  
  plotly::plot_ly(
    data = df_3d,
    x = ~x,
    y = ~y,
    z = ~z,
    color = ~label,
    colors = c("red", "blue", "green"),
    type = "scatter3d",
    mode = "markers",
    marker = list(size = pt.size),
    text = ~ paste("Label:", label),
    hoverinfo = "text"
  ) %>%
    layout(
      scene = list(
        xaxis = list(title = "X"),
        yaxis = list(title = "Y"),
        zaxis = list(title = "Z")
      ),
      legend = list(title = list(text = "Cluster"))
    )
}