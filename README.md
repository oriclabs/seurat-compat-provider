# bl-seurat-provider

`bl-seurat-provider` is an optional, separately installed process for strict
numeric reproduction of selected Seurat 5.5.1 single-cell boundaries. It is
not linked into BioLang and is not required for BioLang's native MIT workflow.

The provider accepts generic row-major `BLMATF64` matrices and emits CCA
embeddings, filter-feature indices, the query weighting reduction, or PCA
scores. BioLang performs anchor search/scoring, integration correction, graph
construction, and clustering itself.

## Licence boundary

This repository is GPL-3.0-only. It invokes separately installed R packages;
their licences and source-distribution obligations remain applicable. BioLang
communicates with this executable only through files and a child-process exit
status. Installing this provider does not add its code to `bl.exe`.

## Requirements

- R with `Rscript` on `PATH`, or `BIOLANG_RSCRIPT` set to its absolute path
- Seurat 5.5.1
- irlba 2.3.7
- RcppAnnoy 0.0.23

The versions are deliberately pinned because approximate decompositions and
nearest-neighbour boundaries can change the final community partition.

BioLang setup after installing the executable:

```text
set BIOLANG_SEURAT_PROVIDER=C:\path\to\bl-seurat-provider.exe
bl doctor
```

In BioLang, request it explicitly; the ordinary call remains native:

```biolang
let anchors = sc.find_integration_anchors(
    control, stimulated, compatibility: "external"
)
let integrated = sc.integrate_data(anchors)
```

## Build and verify

Install directly from the public source repository:

```text
cargo install --git https://github.com/oriclabs/seurat-compat-provider --locked
bl-seurat-provider doctor
```

Or build a checkout locally:

```text
cargo build --release
target/release/bl-seurat-provider doctor
target/release/bl-seurat-provider self-test
```

Commands:

```text
bl-seurat-provider cca left.f64 right.f64 output 30 42 200
bl-seurat-provider pca integrated.f64 output 50 42 false
```

Set `BL_SEURAT_ALLOW_VERSION_MISMATCH=true` only for an explicitly disclosed
experiment. Such a run is not covered by the exact-parity validation.

## Protocol

`BLMATF64` consists of the eight bytes `BLMATF64`, little-endian `u64` rows and
columns, then row-major IEEE-754 `f64` values. See [PROTOCOL.md](PROTOCOL.md).
