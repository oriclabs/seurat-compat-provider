#!/usr/bin/env Rscript

# GPL-3.0-only process adapter. It calls separately installed R packages and
# writes dependency-neutral matrices; no R or package code is linked into the
# BioLang executable.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("missing provider command")
command <- args[[1L]]
provider_version <- Sys.getenv("BL_SEURAT_PROVIDER_VERSION", unset = "unknown")

required_versions <- c(Seurat = "5.5.1", irlba = "2.3.7", RcppAnnoy = "0.0.23")
versions <- vapply(names(required_versions), function(package) {
  if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(package))
}, character(1L))
missing <- names(versions)[is.na(versions)]
if (length(missing)) stop("missing R packages: ", paste(missing, collapse = ", "))
mismatch <- names(versions)[versions != required_versions]
allow_mismatch <- identical(tolower(Sys.getenv("BL_SEURAT_ALLOW_VERSION_MISMATCH")), "true")
if (length(mismatch) && !allow_mismatch) {
  details <- paste0(mismatch, "=", versions[mismatch], " (required ",
                    required_versions[mismatch], ")")
  stop("unvalidated R package versions: ", paste(details, collapse = ", "),
       "; set BL_SEURAT_ALLOW_VERSION_MISMATCH=true only for a disclosed run")
}

if (command == "doctor") {
  cat("provider=bl-seurat-provider\n")
  cat("provider_version=", provider_version, "\n", sep = "")
  cat("protocol=BLSEURAT-1\n")
  cat("license=GPL-3.0-only\n")
  cat("R=", as.character(getRversion()), "\n", sep = "")
  for (package in names(versions)) cat(package, "=", versions[[package]], "\n", sep = "")
  quit(status = 0L)
}

suppressPackageStartupMessages(library(Seurat))

read_blmat <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  magic <- rawToChar(readBin(con, what = "raw", n = 8L))
  if (magic != "BLMATF64") stop("invalid BLMATF64 magic: ", path)
  read_u64 <- function() {
    bytes <- as.double(readBin(con, what = "raw", n = 8L))
    sum(bytes * 256^(0:7))
  }
  rows <- read_u64()
  columns <- read_u64()
  values <- readBin(con, what = "numeric", n = rows * columns,
                    size = 8L, endian = "little")
  if (length(values) != rows * columns) stop("truncated BLMATF64 payload: ", path)
  matrix(values, nrow = rows, ncol = columns, byrow = TRUE)
}

write_u64 <- function(con, value) {
  bytes <- raw(8L)
  remaining <- value
  for (index in seq_len(8L)) {
    bytes[[index]] <- as.raw(remaining %% 256)
    remaining <- floor(remaining / 256)
  }
  writeBin(bytes, con)
}

write_blmat <- function(path, matrix) {
  con <- file(path, "wb")
  on.exit(close(con))
  writeBin(charToRaw("BLMATF64"), con)
  write_u64(con, nrow(matrix))
  write_u64(con, ncol(matrix))
  writeBin(as.double(t(matrix)), con, size = 8L, endian = "little")
}

integer_arg <- function(index, default) {
  if (length(args) < index || !nzchar(args[[index]])) return(default)
  value <- suppressWarnings(as.integer(args[[index]]))
  if (is.na(value) || value < 1L) stop("argument ", index, " must be a positive integer")
  value
}

logical_arg <- function(index, default) {
  if (length(args) < index || !nzchar(args[[index]])) return(default)
  value <- tolower(args[[index]])
  if (!value %in% c("true", "false")) stop("argument ", index, " must be true or false")
  value == "true"
}

if (command == "cca") {
  if (length(args) < 4L) stop("cca requires LEFT.f64 RIGHT.f64 OUTPUT_DIR")
  left_path <- args[[2L]]
  right_path <- args[[3L]]
  output_dir <- args[[4L]]
  dimensions <- integer_arg(5L, 30L)
  seed <- integer_arg(6L, 42L)
  max_features <- integer_arg(7L, 200L)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  left_cells <- read_blmat(left_path)
  right_cells <- read_blmat(right_path)
  if (ncol(left_cells) != ncol(right_cells)) stop("CCA feature counts differ")
  features <- paste0("feature_", seq_len(ncol(left_cells)))
  ctrl <- t(left_cells)
  stim <- t(right_cells)
  rownames(ctrl) <- rownames(stim) <- features
  colnames(ctrl) <- paste0("left_", seq_len(ncol(ctrl)))
  colnames(stim) <- paste0("right_", seq_len(ncol(stim)))

  started <- proc.time()[["elapsed"]]
  cca <- Seurat:::RunCCA.default(
    object1 = ctrl, object2 = stim, standardize = TRUE,
    num.cc = dimensions, seed.use = seed, verbose = FALSE
  )
  raw <- cca$ccv
  norms <- sqrt(rowSums(raw * raw))
  embedding <- raw / pmax(norms, .Machine$double.eps)
  left_embedding <- embedding[colnames(ctrl), , drop = FALSE]
  right_embedding <- embedding[colnames(stim), , drop = FALSE]

  left_raw <- raw[colnames(ctrl), , drop = FALSE]
  right_raw <- raw[colnames(stim), , drop = FALSE]
  loadings <- ctrl %*% left_raw + stim %*% right_raw
  counts <- vapply(seq_len(100L), function(number) {
    length(unique(unlist(lapply(seq_len(dimensions), function(component) {
      unlist(Seurat:::Top(
        data = loadings[, component, drop = FALSE],
        num = number, balanced = TRUE
      ))
    }))))
  }, integer(1L))
  eligible <- counts[counts < max_features]
  per_dimension <- if (length(eligible)) which.max(eligible) else 1L
  top_features <- unique(unlist(lapply(seq_len(dimensions), function(component) {
    unlist(Seurat:::Top(
      data = loadings[, component, drop = FALSE],
      num = per_dimension, balanced = TRUE
    ))
  })))
  filter_indices <- match(top_features, features) - 1L

  # The integration-weight PCA is a separate decomposition over the centered
  # merged residual matrix. Emitting it here avoids a second provider launch
  # and preserves the exact numeric boundary used by strict integration.
  merged <- cbind(ctrl, stim)
  centered <- merged - rowMeans(merged)
  weight_pca <- Seurat:::RunPCA.default(
    object = centered, npcs = dimensions, seed.use = seed,
    approx = TRUE, verbose = FALSE
  )
  weight_scores <- SeuratObject::Embeddings(weight_pca)[colnames(stim), , drop = FALSE]

  write_blmat(file.path(output_dir, "left-embedding.f64"), left_embedding)
  write_blmat(file.path(output_dir, "right-embedding.f64"), right_embedding)
  write_blmat(file.path(output_dir, "weight-reduction.f64"), weight_scores)
  utils::write.csv(
    data.frame(feature_index = filter_indices),
    file.path(output_dir, "filter-features.csv"),
    row.names = FALSE, quote = FALSE
  )
  utils::write.csv(data.frame(
    protocol = "BLSEURAT-1", provider = "bl-seurat-provider",
    provider_version = provider_version,
    seurat = versions[["Seurat"]], irlba = versions[["irlba"]],
    left_cells = nrow(left_cells), right_cells = nrow(right_cells),
    features = ncol(left_cells), dimensions = ncol(left_embedding),
    filter_features = length(filter_indices), seed = seed,
    elapsed_seconds = proc.time()[["elapsed"]] - started
  ), file.path(output_dir, "manifest.csv"), row.names = FALSE, quote = FALSE)
  quit(status = 0L)
}

if (command == "pca") {
  if (length(args) < 3L) stop("pca requires INPUT.f64 OUTPUT_DIR")
  input_path <- args[[2L]]
  output_dir <- args[[3L]]
  components <- integer_arg(4L, 50L)
  seed <- integer_arg(5L, 42L)
  center <- logical_arg(6L, FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cells <- read_blmat(input_path)
  object <- t(cells)
  rownames(object) <- paste0("feature_", seq_len(nrow(object)))
  colnames(object) <- paste0("cell_", seq_len(ncol(object)))
  if (center) object <- object - rowMeans(object)
  started <- proc.time()[["elapsed"]]
  pca <- Seurat:::RunPCA.default(
    object = object, npcs = components, seed.use = seed,
    approx = TRUE, verbose = FALSE
  )
  scores <- SeuratObject::Embeddings(pca)
  loadings <- SeuratObject::Loadings(pca)
  write_blmat(file.path(output_dir, "scores.f64"), scores)
  write_blmat(file.path(output_dir, "loadings.f64"), loadings)
  utils::write.csv(data.frame(
    protocol = "BLSEURAT-1", provider = "bl-seurat-provider",
    provider_version = provider_version,
    seurat = versions[["Seurat"]], irlba = versions[["irlba"]],
    cells = nrow(cells), features = ncol(cells),
    components = ncol(scores), centered = center, seed = seed,
    elapsed_seconds = proc.time()[["elapsed"]] - started
  ), file.path(output_dir, "manifest.csv"), row.names = FALSE, quote = FALSE)
  quit(status = 0L)
}

stop("unknown provider command: ", command)
