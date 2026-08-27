# CLAUDE.md

Project instructions for Claude Code. Read this fully before proposing changes.

## What this is

A single-nucleus RNA-seq analysis of MASH (metabolic dysfunction-associated
steatohepatitis) human liver, sitting in the fibrosis and immunology space.

The question being answered:
> Find an scRNA-seq or snRNA-seq dataset for MASH and highlight which cell types are most perturbed and how, vs. healthy controls.

Deliverable: a Quarto document (`ahn_mash_analysis.qmd`) rendered to a single
self-contained HTML report, plus a reveal.js deck (`ahn_mash_slides.qmd`) that
reads the report's cached results.

## Deliverable and output format

Author the analysis as `ahn_mash_analysis.qmd` and render to a standalone HTML file.

A note on the extension, because it matters. A `.qmd` renders to `.html`, not `.nb.html`. The `.nb.html` format is RMarkdown's R Notebook output (`output: html_notebook` in an `.Rmd`), not a Quarto artifact. The Quarto-native equivalent of an R Notebook is `format: html` with `embed-resources: true`, which produces one standalone `ahn_mash_analysis.html` that opens in any browser with no external dependencies. Default to that. Only if the literal `.nb.html` artifact is a hard submission requirement, author the file as `ahn_mash_analysis.Rmd` with `output: html_notebook` instead. Do not silently switch formats; if this looks like it needs `.nb.html` specifically, ask.

YAML header to use:
```yaml
---
title: "Cell-type perturbation in MASH vs healthy liver"
author: "Richard Ahn, PhD"
format:
  html:
    embed-resources: true
    toc: true
    toc-location: left
    code-fold: show
    theme: cosmo
    df-print: paged
execute:
  warning: false
  message: false
engine: knitr
---
```

Render with `quarto render ahn_mash_analysis.qmd` (produces `ahn_mash_analysis.html`). Preview live with `quarto preview ahn_mash_analysis.qmd`.

## The dataset is already decided. Do not re-choose it.

Single-cell eQTL MASH study (Hong et al., Nature Genetics 2025, s41588-025-02237-8). Human liver biopsies, MASH plus healthy controls, snRNA-seq. This is the pick because its processed data is downloadable now, not just raw reads.

Cohort structure, confirmed from the paper (this is what the metadata should reflect after QC): 249,233 nuclei from 48 donors, staged into four disease levels: No-MASLD (control, n=23), MASL (isolated steatosis, n=6), early MASH (eMASH, fibrosis F0-2, n=9), and advanced MASH (aMASH, fibrosis F3-4, n=10). The cohort is a single-center Korean cohort (Boramae / Seoul National University), biopsied at age 17-78.

Cell types annotated in the object: Hepatocyte, Cholangiocyte, Stellate cell (HSC), Endothelial cell, Macrophage, Monocyte, DC, CD4 T cell, CD8 T cell, Other T cell, NK, B cell, plus Misc and Unidentified. Nonimmune cells were annotated manually; immune cells with Azimuth.

How the shipped object was built (so Phase 1 reuses it correctly, and Phase 2 can mirror or improve it): Cellranger for alignment/counting, Cellbender for ambient-RNA removal, DoubletFinder for doublets, a Seurat v4 object, Liger for batch correction/integration, Seurat for clustering. Note the integration is Liger-based, not Harmony or scVI. Reference software: Seurat 4.3.0, R 4.1.0, DESeq2 1.34.0.

Important framing: this is primarily an eQTL paper. Its headline is the EFHD1 / FOXO1 / rs13395911 hepatocyte-maladaptation story. The cell-type perturbation question addressed here (composition and within-cell-type expression vs control) is described but under-developed in the paper. That is the opening: doing that perturbation analysis cleanly and rigorously is adding value, not recapitulating the paper's main result.

- Processed data (use this): Zenodo record 14586466, one 6.6 GB zip, CC-BY 4.0.
  `https://zenodo.org/records/14586466/files/Zenodo_250212.zip?download=1`
  Contains a processed snRNA-seq object in Seurat v4 format with cell-level and sample-level metadata, plus per-sample cellranger gene-expression matrices.
- Raw / accession: GEO GSE289173 (BioProject PRJNA1221860). Code: github.com/snu-mchoi-lab/MASLD-sceQTL.

Download and unzip:
```bash
curl -L -o mash_zenodo.zip "https://zenodo.org/records/14586466/files/Zenodo_250212.zip?download=1"
unzip mash_zenodo.zip -d data/
```

Provenance line for the report: snRNA-seq of 48 human liver biopsies across the MASLD spectrum, deposited under GEO GSE289173 (BioProject PRJNA1221860) and distributed as a processed Seurat object via Zenodo (10.5281/zenodo.14586466), from Hong et al., Nature Genetics 2025.

## Two-phase plan

Phase 1 (primary): use the processed Seurat object as-is to answer the prompt. Validate the shipped annotations, then measure cell-type perturbation. This is the fast path and a defensible reading of the task.

Phase 2 (fallback, full reproducibility): only if Phase 1 cannot cleanly identify perturbed cell types vs controls, reprocess from the per-sample cellranger matrices with Seurat v5 + BPCells. Do not jump to Phase 2 by default. Re-clustering a published, already-integrated object risks recapitulating their result rather than adding to it.

## Environment and commands (R)

Use R with Seurat v5. Lock package versions with `renv` for reproducibility (the R analogue of a lockfile). Quarto CLI must be installed.

```r
install.packages("renv"); renv::init()
# core
install.packages(c("Seurat", "tidyverse", "patchwork", "quarto"))
# on-disk matrices for the from-scratch path
install.packages("BPCells", repos = c("https://bnprks.r-universe.dev", getOption("repos")))
# differential abundance + pseudobulk DE (Bioconductor)
install.packages("BiocManager")
BiocManager::install(c("speckle", "DESeq2", "SingleCellExperiment"))
renv::snapshot()
```

Notes:
- Apple Silicon (M1, 16 GB): install the arm64 R build; do not run under Rosetta. BPCells compiles C++, so Xcode Command Line Tools must be present (`xcode-select --install`). Seurat, BPCells, speckle, DESeq2 all build on arm64.
- Commit `ahn_mash_analysis.qmd`, the rendered `ahn_mash_analysis.html`, and `renv.lock`. Gitignore the data: `*.zip`, `*.rds`, `data/`, and any `bpcells_*/` on-disk matrix directories (they are large and regenerable).

## Analysis approach (both phases share this logic)

Answer the prompt along two separate axes, because "most perturbed" has two meanings and they should not be conflated:

1. Differential abundance: which cell types change in proportion between MASH and control. Tool: `speckle::propeller`, with the donor as the unit. Optionally `miloR` for neighbourhood-level resolution as a second pass.
2. Within-cell-type differential expression: how each cell type is transcriptionally perturbed. Pseudobulk, one pseudo-sample per donor per cell type, then DESeq2. Note: the paper itself used cell-level Seurat FindMarkers for its differential expression, which treats each nucleus as an independent replicate and inflates significance. Donor-level pseudobulk is the deliberate, defensible improvement, and it is worth stating that choice plainly in the report.

Contrast: MASH vs healthy control. The metadata carries four disease levels (No-MASLD, MASL, eMASH, aMASH). Map them: Control = No-MASLD (n=23); MASH = eMASH + aMASH (n=19); MASL (n=6) is the intermediate and is not the headline. Start with the two-group MASH-vs-control contrast. Two useful secondary views: report aMASH separately where the signal is strongest (advanced fibrosis), and, since the four levels are ordered, an ordered/trend analysis across No-MASLD to aMASH shows the monotonic progression. Do the primary two-group contrast first.

## Method decisions that must hold

- The replication unit is the donor, never the cell. Pseudobulk aggregates counts per donor per cell type.
- Do NOT put donor in the DE design as a covariate. Each donor is one condition, so donor is nested within and confounded with condition; a donor term would absorb the effect being tested. This is a correctness issue, not a style choice.
- Adjust only for confounders that vary within condition (e.g. sequencing batch), and only if they exist. Before finalizing the design, run `table(batch, condition)` and `table(donor, condition)` and show them. This is a single-center Korean cohort, so batch-condition confounding risk is lower than a multi-dataset atlas, but demonstrating the check is part of the deliverable.
- Annotation is trust-but-verify. Confirm shipped labels express canonical liver markers (dotplot) before building any claim on them. Do not blindly re-cluster in Phase 1.
- Reuse the shipped integration/embedding in Phase 1. Only re-integrate in Phase 2.
- QC in Phase 1 is validation of the processed object, not re-derivation.

## Phase 1 workflow (processed object)

Staged workflow, in order:
1. Load the `.rds`, then `UpdateSeuratObject()` to bring the v4 object to v5.
2. Inspect `obj@meta.data`; print factor/character columns and their value counts. Set config: cell-type column, condition column, donor column, batch column (or NA). Do not assume names.
3. Derive a clean two-level `condition` (MASH / Control) and print the table.
4. Sanity QC (`VlnPlot` of nFeature/nCount/percent.mt by condition) and the confound crosstabs above.
5. Validate annotations with a canonical liver-marker dotplot (see Domain reference).
6. Differential abundance with `propeller` (Control vs MASH), plus `plotCellTypeProps`.
7. Pseudobulk DE: `AggregateExpression(..., slot = "counts", group.by = c(celltype, donor))`, build per-pseudosample metadata, run DESeq2 per cell type with `~ condition` (add `+ batch` only if warranted), FDR filter.
8. Synthesis: one table per cell type with abundance change, DEG count, and top up/down genes, ranked. That table is the answer.

`AggregateExpression` needs a raw `counts` layer. If the processed object shipped only normalized data, pull counts from the included per-sample cellranger matrices instead.

## Phase 2 workflow (from scratch, Seurat v5 + BPCells, M1 16 GB)

Trigger only if Phase 1 fails to resolve perturbed cell types. The point of BPCells is that counts and normalized layers live on disk, so the full atlas never sits in RAM.

Pipeline:
1. For each sample's cellranger matrix, write an on-disk BPCells matrix: `BPCells::write_matrix_dir(mat, dir = "bpcells_mats/<sample>")`, then `open_matrix_dir()`. For 10x h5, `open_matrix_10x_hdf5()`.
2. Build a Seurat v5 object whose counts assay is the (list of) on-disk BPCells matrices, one layer per sample. Add per-cell metadata (donor, condition, batch).
3. QC filter (min genes, max percent.mt); keep it explicit. The reference did ambient-RNA removal with Cellbender and doublet removal with DoubletFinder before clustering; mirror those if reprocessing from raw. The reference integrated with Liger; Harmony (below) is a fine and lighter substitute, but note in the report that the original used Liger so any clustering differences are expected.
4. `NormalizeData` then `FindVariableFeatures` (2000). Run `ScaleData` on variable features only, not all genes; scaling all genes is the memory blow-up to avoid on 16 GB.
5. `RunPCA` on the scaled HVGs.
6. Integrate across sample/donor: split layers (`obj[["RNA"]] <- split(obj[["RNA"]], f = obj$donor)`), then `IntegrateLayers(method = HarmonyIntegration, ...)`. Harmony is the lightest choice and the right default; RPCA/CCA are heavier alternatives.
7. `FindNeighbors` on the integrated reduction, `FindClusters`, `RunUMAP`.
8. Annotate clusters by markers (Domain reference), then run the same differential-abundance and pseudobulk-DE analyses as Phase 1.
9. `JoinLayers()` before pseudobulk aggregation.

Memory rules on 16 GB:
- Keep counts and normalized layers on disk via BPCells; do not coerce them to in-memory matrices.
- Scale only variable features.
- Restart the R session between Phase 1 and Phase 2 so you are not holding two large objects at once.
- Prototype on a 3-4 sample subset before running all donors.

## Domain reference (use these; do not invent markers)

Canonical human liver cell-type markers:
- Hepatocyte: ALB, APOA1, TTR, CYP2E1, HNF4A, ASGR1
- Cholangiocyte: KRT19, KRT7, EPCAM, SOX9, CFTR
- Hepatic stellate cell (HSC): COL1A1, COL3A1, DCN, PDGFRB, RGS5; activated/myofibroblast adds ACTA2, TIMP1
- Liver sinusoidal endothelial cell (LSEC): PECAM1, VWF, CLEC4G, STAB2, LYVE1
- Kupffer cell (homeostatic macrophage): CD68, MARCO, VSIG4, C1QA, CD5L
- Scar / lipid-associated macrophage: TREM2, CD9, GPNMB, SPP1, LYZ
- T cell: CD3D, CD3E, CD8A, IL7R; NK: NKG7, GNLY, KLRD1
- B cell: CD79A, MS4A1; Plasma: MZB1, IGHG1

Expected MASH perturbations, to sanity-check the output. Two sources: general MASH biology, and the specific findings reported in this paper's own figures. If your analysis recovers these, the pipeline is corroborated; if not, that is the trigger to consider Phase 2. Recover them independently, do not assert them from the paper.

From this paper specifically (Fig 1):
- Composition, the headline abundance result: hepatocyte proportion drops roughly by half across the progression (about 48% in No-MASLD to about 21% in aMASH), with hepatocytes replaced by immune cells, and hepatocyte proportion inversely correlated with serum ALT. This is the clearest "which cell type is most perturbed in abundance" signal and your propeller analysis should recover the hepatocyte decline and immune expansion.
- Cell-cell signaling (their CellChat result): IGF signaling from hepatocytes to endothelial cells is downregulated in MASH (an anti-inflammatory/anti-fibrotic signal being lost); Notch signaling from cholangiocytes and CD6 signaling from CD4 T cells to hepatocytes are activated in MASH. Optional to reproduce (CellChat), but useful context for the "how" half of the prompt.
- Hepatocytes show the most transcriptional change with disease severity, and MASLD heritability is enriched in hepatocytes (MAGMA Celltyping). Hepatocyte programs: reduced metabolic function with increased cell adhesion/migration in MASH (their Hep-M6 up), loss of a FOXO1-linked metabolic-adaptation program (Hep-M12 down). A progenitor cholangiocyte state expands in advanced MASH.

From general MASH biology (should also appear):
- Hepatic stellate cells shift toward activated collagen-producing myofibroblasts (COL1A1, COL3A1, ACTA2, TIMP1 up). The central fibrosis axis.
- Expansion of TREM2+/CD9+/GPNMB+/SPP1+ scar- or lipid-associated macrophages replacing homeostatic Kupffer cells (MARCO, VSIG4, CD5L down in the resident pool).
- LSEC capillarization (loss of CLEC4G/STAB2 sinusoidal identity).
- Ductular / cholangiocyte reaction (KRT19, KRT7, SOX9).
- Hepatocyte injury and metabolic reprogramming.

## Writing conventions for the report

Written for a scientific reader, including a reproducibility-minded reviewer.
- No em dashes. No awkward semicolons. Avoid "genuinely", "honestly", "actually".
- Do not sound AI-generated. Plain, direct, technical.
- Never write a number, fold change, p-value, or cell count into prose that has not been computed by a chunk. Use inline R (`` `r ... ` ``) to pull computed values, or write the sentence after the chunk runs. No placeholder statistics.
- Do not paste large passages from the source paper. Paraphrase, cite briefly.
- Use `code-fold: show` so the code is visible but collapsible; comment chunks for a reader following the reasoning.
- Lead the report with the two-axis framing and end with the ranked cell-type table that answers the prompt.

## Guardrails

- Ask before adding heavy dependencies beyond those listed.
- Do not change the method decisions or the DE design (especially donor-not-a-covariate) without flagging it.
- Do not enter Phase 2 without a stated reason that Phase 1 was insufficient.
- Do not commit data artifacts; commit `ahn_mash_analysis.qmd`, `ahn_mash_analysis.html`, `renv.lock`.
- If the analysis contradicts expected MASH biology or yields a weak signal, say so plainly rather than forcing the narrative. An honest negative, with the diagnosis, is a valid and valuable result.
