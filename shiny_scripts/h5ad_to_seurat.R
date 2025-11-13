#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(reticulate)
  library(Matrix)
  library(Seurat)
})

py_is_none <- reticulate::py_is_none

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    action = "append",
    help = "Path(s) to input .h5ad file(s). Provide this flag multiple times for multiple files."
  ),
  make_option(
    c("-o", "--output-dir"),
    type = "character",
    default = ".",
    help = "Directory to write Seurat .rds files (default: current directory)."
  ),
  make_option(
    c("--assay"),
    type = "character",
    default = "RNA",
    help = "Name of the assay to use in the Seurat object (default: RNA)."
  ),
  make_option(
    c("--layer"),
    type = "character",
    default = NULL,
    help = "Optional AnnData layer to use for counts. If not provided, raw counts are used when available, otherwise X."
  ),
  make_option(
    c("--umap-key"),
    type = "character",
    default = "X_umap",
    help = "Key inside adata$obsm containing UMAP coordinates (default: X_umap)."
  )
)

parser <- OptionParser(option_list = option_list,
                       description = "Convert AnnData (.h5ad) files to Seurat objects (.rds) by extracting counts, metadata, and embeddings.")
opts <- parse_args(parser)

if (is.null(opts$input)) {
  print_help(parser)
  stop("At least one --input file must be specified.", call. = FALSE)
}

# Configure Python module imports once.
anndata <- import("anndata", convert = FALSE)
np <- import("numpy", convert = FALSE)
scipy_sparse <- import("scipy.sparse", convert = FALSE)

dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

layer_keys <- function(adata) {
  if (py_is_none(adata$layers)) {
    return(character(0))
  }
  py_to_r(adata$layers$keys())
}

to_dgc <- function(mat) {
  if (inherits(mat, "dgCMatrix")) {
    return(mat)
  }
  if (inherits(mat, "dgRMatrix")) {
    return(as(mat, "dgCMatrix"))
  }
  Matrix::Matrix(mat, sparse = TRUE)
}

fetch_layer <- function(adata, name) {
  keys <- layer_keys(adata)
  if (!(name %in% keys)) {
    return(NULL)
  }
  mat_r <- py_to_r(adata$layers[[name]])
  to_dgc(mat_r)
}

align_matrix <- function(mat, cell_ids, gene_ids) {
  if (nrow(mat) == length(gene_ids) && ncol(mat) == length(cell_ids)) {
    rownames(mat) <- gene_ids
    colnames(mat) <- cell_ids
    return(mat)
  }
  if (nrow(mat) == length(cell_ids) && ncol(mat) == length(gene_ids)) {
    mat <- Matrix::t(mat)
    rownames(mat) <- gene_ids
    colnames(mat) <- cell_ids
    return(mat)
  }
  stop("Matrix dimensions do not match AnnData obs/var.")
}

get_embeddings <- function(adata, key, cells) {
  if (py_is_none(adata$obsm)) {
    return(NULL)
  }
  obsm_names <- py_to_r(adata$obsm$keys())
  if (!(key %in% obsm_names)) {
    warning(sprintf("UMAP key '%s' not present in adata$obsm; skipping embeddings.", key))
    return(NULL)
  }
  emb <- py_to_r(adata$obsm[[key]])
  emb <- as.matrix(emb)
  rownames(emb) <- cells
  colnames(emb) <- paste0("UMAP_", seq_len(ncol(emb)))
  emb
}

safe_make_unique <- function(x) {
  make.unique(as.character(x), sep = "_")
}

convert_one <- function(input_path) {
  message(sprintf("Reading %s", input_path))
  adata <- anndata$read_h5ad(input_path)

  cell_ids <- safe_make_unique(py_to_r(adata$obs_names))
  gene_ids <- safe_make_unique(py_to_r(adata$var_names))
  keys <- layer_keys(adata)

  # Counts matrix
  counts_mat <- NULL
  if (!is.null(opts$layer)) {
    counts_mat <- fetch_layer(adata, opts$layer)
    if (is.null(counts_mat)) {
      stop(sprintf("Requested layer '%s' not found in AnnData object.", opts$layer))
    }
  }
  if (is.null(counts_mat) && "counts" %in% keys) {
    counts_mat <- fetch_layer(adata, "counts")
  }
  if (is.null(counts_mat) && !py_is_none(adata$raw) && !py_is_none(adata$raw$X)) {
    raw_adata <- adata$raw$to_adata()
    raw_counts <- to_dgc(py_to_r(raw_adata$X))
    raw_genes <- safe_make_unique(py_to_r(raw_adata$var_names))
    idx <- match(gene_ids, raw_genes)
    if (anyNA(idx)) {
      stop("Unable to align raw counts with current feature set.")
    }
    counts_mat <- raw_counts[, idx, drop = FALSE]
  }
  if (is.null(counts_mat)) {
    counts_mat <- to_dgc(py_to_r(adata$X))
  }
  counts_mat <- align_matrix(counts_mat, cell_ids, gene_ids)

  # Normalised (log) matrix
  norm_mat <- NULL
  if ("log_normalized" %in% keys) {
    norm_mat <- align_matrix(fetch_layer(adata, "log_normalized"), cell_ids, gene_ids)
  }

  # Scaled matrix (dense)
  scale_mat <- NULL
  if ("scaled" %in% keys) {
    scale_mat <- align_matrix(fetch_layer(adata, "scaled"), cell_ids, gene_ids)
    scale_mat <- as.matrix(scale_mat)
  }

  metadata <- py_to_r(adata$obs)
  metadata <- data.frame(metadata, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(metadata) <- cell_ids

  var_df <- py_to_r(adata$var)
  var_df <- data.frame(var_df, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(var_df) <- gene_ids

  seurat_obj <- CreateSeuratObject(counts = counts_mat, meta.data = metadata, assay = opts$assay)

  if ("feature_name" %in% colnames(var_df)) {
    seurat_obj[[opts$assay]]@meta.features$feature_name <- var_df[colnames(seurat_obj), "feature_name"]
  }

  # Set normalised data if available
  if (!is.null(norm_mat)) {
    norm_mat <- norm_mat[rownames(seurat_obj), colnames(seurat_obj), drop = FALSE]
    seurat_obj <- SetAssayData(seurat_obj, slot = "data", new.data = norm_mat)
  } else {
    message("Normalised layer not found; running Seurat NormalizeData().")
    seurat_obj <- NormalizeData(seurat_obj, assay = opts$assay, normalization.method = "LogNormalize",
                                scale.factor = 1e4, verbose = FALSE)
  }

  # Set scaled data if available
  if (!is.null(scale_mat)) {
    scale_mat <- scale_mat[rownames(seurat_obj), colnames(seurat_obj), drop = FALSE]
    seurat_obj <- SetAssayData(seurat_obj, slot = "scale.data", new.data = scale_mat)
  }

  umap <- get_embeddings(adata, opts$umap_key, cell_ids)
  if (!is.null(umap)) {
    umap <- umap[colnames(seurat_obj), , drop = FALSE]
    seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = umap, key = "UMAP_", assay = opts$assay)
  }

  seurat_obj@misc$h5ad_source <- normalizePath(input_path)

  output_name <- file.path(
    opts$output_dir,
    sprintf("%s.rds", tools::file_path_sans_ext(basename(input_path)))
  )
  saveRDS(seurat_obj, file = output_name)
  message(sprintf("Saved Seurat object to %s", output_name))
}

invisible(lapply(opts$input, convert_one))

