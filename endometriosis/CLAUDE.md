# CLAUDE.md

Project instructions for Claude Code. Read this fully before proposing changes.

## What this is

The task: analyze the scRNA-seq arm of the Human Endometrial Cell Atlas (HECA), compare endometriosis cases vs controls, select a dysregulated cell type, derive a biological insight, and propose a therapeutic intervention. The deliverable is a single Jupyter notebook with clear markdown reasoning that the team can run to reproduce the findings.

Primary artifact: `endometriosis_analysis_richard_ahn.ipynb` (already scaffolded, runs top to bottom, config-driven). Most work is refining and running it, not rebuilding it.

Reference paper: Mareckova et al., Nature Genetics 2024, `s41588-024-01873-w` (HECA). ArrayExpress `E-MTAB-14039`.

Budget: the challenge suggests ~6 hours. Spend compute where it matters (the disease signal), not on re-integrating the full atlas.

## The scientific thesis is already decided. Do not relitigate it.

Headline insight: endometriosis has no disease-specific cell type. The signal is a rewiring of programs inside existing cells. The two compartments implicated by both differential abundance and functional GWAS in the paper are decidualized stromal cells and macrophages. They converge on one druggable node: IGF1.

The argument the notebook builds:
1. Within-cell-type differential expression (case vs control) recovers IGF1 up in uM2 macrophages and an inflammatory shift (TNFRSF1B, CEBPD) in uM1.
2. The same IGF1 upregulation appears in decidualized stromal cells, so the two GWAS-anchored cell types share a node.
3. IGF1 is targetable and disease-relevant: macrophage-derived IGF1 sensitizes pelvic nerves, peritoneal IGF1 tracks with patient pain scores, and the IGF-1R inhibitor linsitinib reverses pain behavior in preclinical endometriosis models.

Therapeutic hypothesis: IGF-1R inhibition as a pain-directed, genetically anchored intervention with existing clinical-stage chemistry.

If you think a different angle is stronger, raise it as a question first. Do not silently redirect the analysis.

## Data

Not in the repo. Download it. Public object store, no auth.

- Full scRNA-seq atlas with raw counts (use this for the final run):
  `https://cellgeni.cog.sanger.ac.uk/vento/reproductivecellatlas/endometriumAtlasV2_cells_with_counts.h5ad`
- Immune-only subset (smaller, good for iterating the macrophage work):
  `https://cellgeni.cog.sanger.ac.uk/vento/reproductivecellatlas/endometriumAtlasV2_cells_immune.h5ad`
- Downsampled ~90k-cell subset for fast prototyping (raw counts in `.X`, annotations in `.obs['celltype']`):
  `https://zenodo.org/records/15072628/files/HECA-Subset.h5ad`

Ignore the `_nuclei` objects. Single-nuclei data is out of scope per the challenge.

Provenance note for the notebook: this is the HECA scRNA-seq object deposited under `E-MTAB-14039`, distributed via the VentoLab reproductive cell atlas portal. State it that way.

Schema notes (verify with the inspection cell, do not assume):
- Cell type column is `celltype` in the objects seen so far, not `cell_type`. Set `CELLTYPE_COL = "celltype"`.
- In the `_with_counts` file, counts are likely in `.X`, so `RAW_COUNTS_LAYER = None` may be correct. The loader has fallback logic but confirm.
- Condition, donor, dataset, and menstrual-phase column names are unverified. Read them off the inspection cell and set the config before running anything downstream.

## Environment and commands

Use `uv`, not conda. Project mode gives a `pyproject.toml` plus `uv.lock`, which is the reproducibility the assignment asks for. Execution happens through a **Jupyter MCP server**: Claude drives a live kernel (add cells, run them, read outputs, inspect kernel state) rather than batch-running the notebook. The one hard requirement is that the kernel the MCP talks to is this uv environment. If it is not, none of the locked dependencies are importable and everything downstream fails.

```bash
# one-time project setup
uv init --python 3.11 .
uv add scanpy anndata leidenalg scikit-misc decoupler pydeseq2 \
       matplotlib seaborn harmonypy session-info ipykernel jupyterlab

# register the uv environment as a named Jupyter kernel (this is what the MCP must select)
uv run python -m ipykernel install --user --name heca --display-name "HECA (uv)"

# start the Jupyter server the MCP connects to, from inside the uv env
# (set a token and point your Jupyter MCP server config at this URL + token)
uv run jupyter lab --no-browser --port 8888 --IdentityProvider.token=180d74a125efd07aa952ff85a40192ab7e949ede7e27a2c860803e08830b1b48
```

Working through the Jupyter MCP server:
- The notebook `endometriosis_analysis_richard_ahn.ipynb` must run on the **`HECA (uv)`** kernel. Confirm the kernel before executing anything. A cell running `import sys; print(sys.executable)` should point inside `.venv`, not a system Python.
- Do the schema inspection as a live cell, not a shell one-liner. Run the load-and-inspect cell (section 2) first, read the printed `.obs` columns off the kernel output, then set the Config cell. Use the live kernel state to verify column names rather than guessing.
- The kernel is persistent, so a loaded AnnData stays in memory across cell runs. That is good for iterating, but on 16 GB it means memory accumulates. Restart the kernel to free RAM before loading the full atlas, and after the immune-subset work if you then need the full object. Do not hold two large objects in the kernel at once.
- Before treating the notebook as done, do a clean end-to-end pass: restart the kernel and run all cells top to bottom, on the `HECA (uv)` kernel, so the submission is verified to reproduce from a cold start rather than from accumulated interactive state. This is the equivalent of the old headless `nbconvert` check and it still matters for the deliverable.

Notes:
- Everything runs through `uv run` and the `HECA (uv)` kernel, both of which resolve against the locked environment. Do not `pip install` into a stray interpreter, and do not let the MCP fall back to a default system kernel.
- If a package needs to be added mid-session, use `uv add <pkg>`, not `uv pip install`, so it lands in the lockfile. After `uv add`, the running kernel will not see the new package until it is restarted.
- On Apple Silicon (M1, 16 GB), all of the above have native arm64 wheels. Make sure the interpreter is arm64, not running under Rosetta.
- `harmonypy` is only used for the de novo immune re-clustering. `decoupler` and `pydeseq2` drive the differential expression (see below).
- Commit `pyproject.toml` and `uv.lock`. The lockfile is how the team reproduces the environment.

## Memory (M1, 16 GB)

The full atlas will not fit comfortably if loaded naively. Load backed and subset before materializing.

```python
import scanpy as sc, gc
adata = sc.read_h5ad("endometriumAtlasV2_cells_with_counts.h5ad", backed="r")  # X stays on disk
mask = adata.obs["celltype"].str.contains("uM|uNK|NK|T |DC|mast|B |plasma", case=False, na=False)
imm = adata[mask.values].to_memory()   # only the immune compartment in RAM
adata.file.close(); del adata; gc.collect()
```

Prototype the whole notebook against the Zenodo ~90k subset first, then do the final run against the full object in backed mode. Pseudobulk DE is aggregation, so it is memory-light regardless.

## Notebook structure

Config-driven. The only cell that should need editing per environment is the Config cell.

0. Setup and reproducibility (seeds, versions)
1. Config (paths and `.obs` column names)
2. Load and inspect (prints schema so Config can be set)
3. QC and normalization (validation of the shipped object against paper thresholds, plus normalize for expression work; not a re-derivation)
4. Annotation: 4a validate shipped labels with markers, 4b de novo re-annotate the immune compartment
5. Disease signal: 5a composition (exploratory only), 5b within-cell-type pseudobulk DESeq2 (primary), 5c IGF1 convergence check
6. Biological insight and therapeutic hypothesis (markdown)
7. Reproducibility notes

## Method decisions that must hold

These are deliberate and signal the seniority the reviewers are checking for. Do not undo them without flagging it.

- Reuse the shipped scVI embedding. Do not re-run scVI on the 313K-cell atlas. The only from-scratch integration is on the immune subset.
- QC is validation, not re-derivation. The atlas ships filtered. Show the per-cell metrics and confirm the object honors the paper thresholds (>1000 genes, <20% mito), split by donor and dataset. Do not re-filter to rebuild the object. Normalize (normalize_total + log1p) because downstream expression work needs it. Say in markdown that this is a published integrated reference, so QC validates the shipped object and reuses its integration rather than re-deriving either.
- Composition testing on the scRNA-seq atlas is confounded because the seven source datasets differ in cell-type recovery. Treat abundance as exploratory. Anchor the disease claim in within-cell-type differential expression.
- Annotation is marker-anchored. Cluster indices may shift across library versions. Labels should not.
- Retain exogenous-hormone donors and model treatment as a DE covariate; do not exclude them. (Superseded 2026-08-10, at the user's direction: excluding dropped ~53k cells and, because treatment is confounded with case status, preferentially removed cases (14 case vs 2 control donors). The primary contrast is now reported both without and with a `hormone` term, `~ dataset + condition` vs `~ dataset + hormone + condition`. Excluding vs including was checked directly and does not change the IGF1 null.)

## Differential expression recipe (the primary evidence)

Chosen approach: Python pseudobulk aggregation, then DESeq2 via `pydeseq2`. This is a deliberate choice against the authors' limma-voom. Pseudobulk methods (DESeq2, edgeR, limma-voom) are well calibrated and largely concordant, so the top hits (IGF1) should survive the method change. Name the choice in markdown and cite the authors' notebook.

Rules:
- Replication unit is the donor. Build one pseudosample per donor per cell type. Not cells. Aggregate raw counts (decoupler `get_pseudobulk`, or a manual per-donor sum).
- The authors split each donor into three random metacells to raise replicate count. Do not copy that by default. Three splits of one donor are not independent replicates and make the test anticonservative. Use donor-level pseudobulk as primary and note this deviation in one line.
- Design: `~ dataset + condition`, test the condition coefficient (case vs control). Adjust only for confounders that vary within condition (dataset, and optionally menstrual phase).
- Do NOT put donor in the design as a covariate. In this case/control cohort each donor belongs to one condition only, so donor is nested within and confounded with condition. A donor term would absorb the effect being tested. This is the single most important thing to get right.
- Before fixing the formula, run `pd.crosstab(dataset, condition)` and `pd.crosstab(phase, condition)`. If a study is entirely cases or entirely controls, dataset is partially collinear with condition and cannot be fully adjusted. Note that as a limitation.
- Do NOT route DESeq2 through Seurat `FindMarkers(test.use="DESeq2")`. That runs a bare `~group` design with no covariate control and would mean rebuilding a Seurat object around a matrix you already aggregated in Python. Call DESeq2 directly via `pydeseq2`.
- Gene filter before testing: expressed in a reasonable fraction of the compartment (the authors used >10% of cells). Apply something equivalent.

Authors' reference method, for concordance checking: limma-voom pseudobulk (edgeR TMM into voom into lmFit into eBayes), >10% expression gene filter, contrast control vs endometriosis. See `Endometriosis_degs/Endometriosis_degs2_DEGanalysis_cells.ipynb` in `ventolab/Human-Endometrial-Cell-Atlas_HECA`. If exact reproduction is ever wanted, the option is aggregating in Python and fitting limma-voom in R via rpy2 with only `bioconductor-limma` and `bioconductor-edger` installed. Not the chosen path, but keep it in mind if a reviewer asks for a limma-voom cross-check.

## Domain reference (use these, do not invent markers or genes)

QC thresholds (from the paper Methods): drop cells with < 1000 detected genes or > 20% mitochondrial content. HVGs: 2000, Seurat v3 flavor.

Cell-state markers (paper plus canonical):
- Epithelial lineage: EPCAM, KRT8, KRT18, CDH1. SOX9 basalis (CDH2+): SOX9, CDH2, AXIN2, ALDH1A1. Luminal: LGR5. Ciliated: FOXJ1, PIFO, TP73. preGlandular: OPRK1, CBR3, SUFU, HPRT1. preLuminal: SULT1E1. MUC5B (likely cervical): MUC5B, TFF3, SAA1, BPIFB1.
- Stromal lineage: PDGFRB, DCN, LUM, COL6A3. eStromal MMPs: MMP1, MMP10, MMP3, INHBA. dStromal early: PLCL1. dStromal mid: DKK1. dStromal late: LEFTY2, SMAD7.
- Endothelial: PECAM1, VWF, CLDN5. Perivascular: ACTA2, NOTCH3, RGS5, MYH11.
- Immune: PTPRC. Macrophage (uM): CD68, CD14, LYZ, C1QA, C1QB. uM1 (pro-inflammatory): IL1B, EREG, TNF. uM2 (anti-inflammatory / tissue-resident): HMOX1, FOLR2, LYVE1. uNK: NCAM1, KLRD1, NKG7, GNLY. T cell: CD3D, CD3E, CD8A. DC: CLEC9A, CD1C, LILRA4. Mast: TPSAB1, CPA3, MS4A2. B/Plasma: MS4A1, CD79A, MZB1.

Expected case vs control DE (from paper Fig 5, use to sanity-check, not to hardcode):
- Macrophages: IGF1 up in uM2; TNFRSF1B and CEBPD up in uM1.
- Decidualized stromal: IGF1 up, IGF2 down, GREB1 up, DKK1 down.
- IGF-axis genes the authors track: IGF1, IGF2, IGF1R, IGFBP1, IGFBP3, INSR, IRS2, DKK1, GREB1.

## Writing conventions for the notebook markdown

This is a submission that a scientific team will read. The prose should read as written by Richard.

- No em dashes. No awkward semicolons. Do not use "genuinely", "honestly", "actually".
- Do not sound AI-generated. Plain, direct, technical.
- Do not fabricate numbers, figures, or results in markdown. If a value depends on the run, leave a placeholder or write it after the cell executes. Never assert a fold change or p-value that has not been computed.
- Do not paste large passages from the paper. Paraphrase, cite briefly.
- Comment code for a reader following the thought process, per the assignment.

## Guardrails

- Ask before adding heavy dependencies beyond what `uv add` already installed.
- Add packages with `uv add` so they land in the lockfile. Do not `pip install` into a stray interpreter.
- Do not commit the `.h5ad` files. Add `*.h5ad` and `.venv/` to `.gitignore`. Do commit `pyproject.toml` and `uv.lock`.
- Do not change the method decisions or the DE recipe above without flagging it. The donor-as-covariate point in particular is a correctness issue, not a style preference.
- If the run surfaces a result that contradicts the IGF1 thesis, say so plainly rather than forcing the narrative. An honest negative is fine and worth reporting.

## Rules for Claude when working with Jupyter notebooks

### Tool preference
- Use the Jupyter MCP for all `.ipynb` operations — read, edit, insert, delete, execute.
- Do not use your built-in `NotebookEdit` tool; it writes source as a single JSON string, which ruins standard Jupyter formatting.

### Outputs
- Never print secrets, API keys, tokens, or passwords into cell output.
- Large outputs consume tokens and fill up your context window. Prefer summaries (`.head()`, `.shape`) over dumping full DataFrames.

### Execution
- When installing packages, use `%pip install` inside the notebook (not `!pip install`) so packages install into the running kernel.
- Execute cells to verify they work. Do not assume the code is correct.
- If a cell errors, read the actual traceback before attempting a fix. Do not guess.

### State and reproducibility
- Jupyter kernels are stateful. A notebook that runs top-to-bottom after "Restart & Run All" is the only notebook that works — verify this before declaring a task done.

### Data safety
- Do not modify or delete raw data files. Write derived data to a separate path.

