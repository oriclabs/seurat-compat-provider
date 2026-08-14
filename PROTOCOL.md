# BLSEURAT-1 protocol

The provider is invoked as a normal child process. A zero exit status means all
declared outputs were written successfully. Diagnostics go to standard error.

## CCA

```text
bl-seurat-provider cca LEFT.f64 RIGHT.f64 OUTPUT_DIR DIMS SEED MAX_FEATURES
```

Both inputs are cells by features and must have identical feature columns.
Outputs are:

- `left-embedding.f64`: left cells by CCA dimensions;
- `right-embedding.f64`: right cells by CCA dimensions;
- `filter-features.csv`: zero-based input-column indices;
- `weight-reduction.f64`: right cells by CCA dimensions;
- `manifest.csv`: package versions, dimensions, seed, and elapsed time.

## PCA

```text
bl-seurat-provider pca INPUT.f64 OUTPUT_DIR COMPONENTS SEED CENTER
```

Input is cells by features. Outputs are `scores.f64`, `loadings.f64`, and
`manifest.csv`.

Both manifests identify the `BLSEURAT-1` protocol, provider name and version,
and the pinned Seurat/irlba versions. Callers should copy these fields into the
analysis provenance before temporary interchange files are removed.

No object serialization, R expression, executable path, or source code crosses
the protocol. All scientific payloads use neutral matrices or integer CSV.
