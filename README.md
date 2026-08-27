# Single-cell analyses

This is a repo of re-analyses of published sc/snRNA-seq datasets by indication.

## `endometriosis/`

scRNA-seq of the Human Endometrial Cell Atlas (HECA), endometriosis cases versus controls.

The headline finding is negative in a useful way: there is no endometriosis-specific cell
type. The disease signal is a rewiring of expression programs inside cell types that are
already present. Differential abundance and the paper's functional GWAS converge on two
compartments — decidualized stromal cells and macrophages — which in turn converge on a
single druggable node, IGF1.

- **Deliverable:** [`endometriosis_analysis_richard_ahn.ipynb`](endometriosis/endometriosis_analysis_richard_ahn.ipynb) — runs top to bottom, config-driven, with the reasoning carried in markdown alongside the code.
- **Data:** ArrayExpress `E-MTAB-14039`, from Mareckova et al., *Nature Genetics* 2024 (`s41588-024-01873-w`).
- **Stack:** Python / scanpy.

## `mash/`

snRNA-seq of human liver across the MASLD spectrum — which cell types are most perturbed in
MASH versus healthy controls, by composition and by within-cell-type expression.

The source study is primarily a single-cell eQTL paper; the cell-type perturbation question
is described there but under-developed, which is the opening this analysis works in. The
cohort is 249,233 nuclei from 48 donors staged across four disease levels, from no-MASLD
controls through advanced MASH.

- **Deliverable:** [`ahn_mash_analysis.qmd`](mash/ahn_mash_analysis.qmd) and its rendered, self-contained [`ahn_mash_analysis.html`](mash/ahn_mash_analysis.html). GitHub will not render the HTML inline — download it and open it in a browser, or view it through a raw-HTML viewer.
- **Data:** GEO `GSE289173` (BioProject `PRJNA1221860`), distributed as a processed Seurat object via Zenodo (`10.5281/zenodo.14586466`), from Hong et al., *Nature Genetics* 2025 (`s41588-025-02237-8`).
- **Stack:** R / Seurat, authored in Quarto.

## Reproducing

Neither folder ships its input data — both datasets are large and publicly available at the
accessions above. Download the data for the analysis you want to run, point the config at
it, and execute the notebook or render the `.qmd`.

The `CLAUDE.md` in each folder documents that project's methods, design decisions, and the
reasoning behind them in more depth than this index does.
