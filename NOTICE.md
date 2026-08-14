# Notices

`bl-seurat-provider` is original GPL-3.0-only adapter code. It does not contain
or link R, Seurat, irlba, RcppAnnoy, or their dependency source code. At runtime
it starts an independently installed `Rscript` process and asks those packages
to perform numerical operations through their public or namespace APIs.

Validated runtime components:

- Seurat 5.5.1 — MIT + file LICENSE;
- SeuratObject — MIT + file LICENSE;
- irlba 2.3.7 — GPL-3;
- RcppAnnoy 0.0.23 — GPL-2+ wrapper around Apache-2.0 Annoy;
- R — GPL-2 | GPL-3.

Users and redistributors must comply with the licences of the R environment
they install. Those components are not distributed by this repository.

