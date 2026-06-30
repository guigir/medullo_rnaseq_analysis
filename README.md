# Medullo RNA-seq analysis

Reproducible R workflow for bulk RNA-seq differential expression and pathway-level interpretation in a small medulloblastoma-related dataset.

The repository currently contains three main analysis scripts:

```text
scripts/
├── new_de_analysis.R
├── build_excel_from_deseq2_outputs_v2.R
└── make_priority_figures.R
```

The analysis compares:

- `HDMB03_GFP_cells`, used as the reference condition;
- `cb_parenchyma`, used as the condition of interest.

The GitHub repository is described as an RNA-seq analysis of 6 medulloblastoma samples, with 3 pre-graft and 3 post-graft samples.

## Overview

The workflow is organized in three steps.

```text
Raw count matrix
      │
      ▼
1. Differential expression
   scripts/new_de_analysis.R
      │
      ▼
2. Excel result table with raw and normalized counts
   scripts/build_excel_from_deseq2_outputs_v2.R
      │
      ▼
3. Priority figures and Reactome GSEA interpretation
   scripts/make_priority_figures.R
```

The final curated figures focus on:

- differential expression summary with an UpSet plot;
- selected axon guidance / G3 genes on a volcano plot;
- Reactome GSEA;
- SLIT/ROBO-related Reactome pathways;
- top positively enriched Reactome pathways.

## Repository structure

```text
medullo_rnaseq_analysis/
├── data/
│   ├── .gitkeep
│   └── README.md
├── figures/
│   └── .gitkeep
├── results/
│   └── .gitkeep
├── scripts/
│   ├── new_de_analysis.R
│   ├── build_excel_from_deseq2_outputs_v2.R
│   └── make_priority_figures.R
├── docs/
│   └── TUTORIAL.md
├── .gitignore
└── README.md
```

## Input data

The required input file is:

```text
data/expr_matrix.txt
```

This file must contain a raw count matrix:

- rows = genes;
- columns = samples;
- row names = Ensembl gene IDs, for example `ENSG00000000003`;
- values = raw, non-normalized RNA-seq counts.

The current scripts infer the biological condition from sample names:

| Sample name pattern | Assigned condition |
|---|---|
| `cb_parenchyma` | `cb_parenchyma` |
| `HDMB03` or `GFP` | `HDMB03_GFP_cells` |

## Main analysis design

The differential expression model is:

```r
design = ~ condition
```

The main contrast is:

```text
cb_parenchyma vs HDMB03_GFP_cells
```

where `HDMB03_GFP_cells` is the reference condition.

A historical batch variable may appear in some script versions for PCA visualization, but the main differential expression model should remain condition-only because the batch information is collinear with condition in this small dataset.

## Differential expression thresholds

Default thresholds:

```r
alpha  <- 0.05
fc_cut <- 1.5
```

A gene is considered UP when:

```text
padj < 0.05 and log2FoldChange >= log2(1.5)
```

A gene is considered DOWN when:

```text
padj < 0.05 and log2FoldChange <= -log2(1.5)
```

The `log2FoldChange` exported by `new_de_analysis.R` is expected to be the shrunken log2 fold-change when the table is built from the `lfcShrink(..., type = "apeglm")` result.

## Installation

### 1. Check the R version used by the terminal

Before running the workflow, make sure the terminal uses the same R installation as the one where packages are installed.

```bash
which R
which Rscript
Rscript -e 'R.version.string'
Rscript -e '.libPaths()'
```

If the terminal accidentally uses conda's R, deactivate conda:

```bash
conda deactivate
```

Then check again.

### 2. Install R packages

CRAN packages:

```r
install.packages(c(
  "tidyverse",
  "readxl",
  "readr",
  "openxlsx",
  "ggplot2",
  "ggrepel",
  "ComplexUpset",
  "patchwork",
  "svglite",
  "BiocManager"
))
```

Bioconductor packages:

```r
BiocManager::install(c(
  "DESeq2",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "apeglm",
  "limma",
  "clusterProfiler",
  "enrichplot",
  "ReactomePA",
  "msigdbr"
))
```

## Quick start

Run all commands from the repository root.

```bash
cd medullo_rnaseq_analysis
```

### Step 1 — Run DESeq2

```bash
Rscript scripts/new_de_analysis.R
```

Expected outputs include:

```text
DESeq2_redo_clean_v2/
├── DESeq2_results_shrunk_cleaned_symbol.csv
├── UP_SYMBOL.txt
├── DOWN_SYMBOL.txt
└── PCA_rlog_condition.svg
```

### Step 2 — Build Excel workbook

```bash
Rscript scripts/build_excel_from_deseq2_outputs_v2.R
```

Expected output:

```text
DESeq2_redo_clean/DESeq2_results_with_counts.xlsx
```

or, depending on the active script version:

```text
DESeq2_redo_clean_v2/DESeq2_results_with_counts.xlsx
```

The Excel workbook should contain at least the sheet:

```text
ALL
```

with columns including:

```text
ensemblID
symbol
padj
log2FoldChange
FoldChange
FoldChange_signed
```

and normalized count columns ending with:

```text
norm
```

### Step 3 — Generate priority figures

```bash
Rscript scripts/make_priority_figures.R
```

The script can also be run with an explicit Excel file:

```bash
Rscript scripts/make_priority_figures.R DESeq2_redo_clean_v2/DESeq2_results_with_counts.xlsx
```

Expected priority figures:

```text
figures/report/
├── ComplexUpSet_l2fc_0_585_no-legend.svg
├── Volcano_Gene_Interests_no-legend.svg
├── GSEA_SLIT_ROBO_01_Signaling_by_ROBO_receptors.svg
├── GSEA_SLIT_ROBO_02_Regulation_of_expression_of_SLITs_and_ROBOs.svg
├── GSEA_TOP30_POS_with_large_legends.svg
└── GSEA_TOP35_POS_with_large_legends.svg
```

PNG versions may also be generated.

## Output directories

### Differential expression

```text
DESeq2_redo_clean_v2/
```

Contains the main DESeq2 result table and UP/DOWN gene lists.

### Excel table

```text
DESeq2_redo_clean/
```

or:

```text
DESeq2_redo_clean_v2/
```

depending on the active script version.

### Priority figures

```text
figures/report/
```

### Intermediate tables

```text
results/report/
```

This directory contains Reactome GSEA tables, ranked genes, and SLIT/ROBO pathway tables.

## Main final figures

| Figure | Description |
|---|---|
| `ComplexUpSet_l2fc_0_585_no-legend.svg` | UpSet summary of significant UP and DOWN genes at `|log2FC| >= 0.585`. |
| `Volcano_Gene_Interests_no-legend.svg` | Volcano plot highlighting selected axon guidance and G3-related genes. |
| `GSEA_SLIT_ROBO_01_Signaling_by_ROBO_receptors.svg` | Reactome GSEA curve for ROBO receptor signaling. |
| `GSEA_SLIT_ROBO_02_Regulation_of_expression_of_SLITs_and_ROBOs.svg` | Reactome GSEA curve for regulation of SLIT/ROBO expression. |
| `GSEA_TOP30_POS_with_large_legends.svg` | Top 30 positively enriched Reactome pathways with category labels. |
| `GSEA_TOP35_POS_with_large_legends.svg` | Top 35 positively enriched Reactome pathways with category labels. |

## Troubleshooting

### `library(...) : aucun package nommé ...`

The terminal is probably not using the R installation where the packages were installed.

Check:

```bash
which Rscript
Rscript -e 'R.version.string'
Rscript -e '.libPaths()'
```

If conda's R is used by mistake:

```bash
conda deactivate
```

or call the correct `Rscript` explicitly.

### Font error with Arial

On Linux, `Arial` may not be available and may cause plotting errors.

Use:

```r
plot_font <- "sans"
```

The curated `make_priority_figures.R` should use a portable font family.

### Blank SLIT/ROBO SVG files

Some versions of `enrichplot::gseaplot2()` can produce blank SVG files when saved directly with `ggsave()`.

The recommended solution is to save those plots by opening the SVG/PNG device explicitly and using:

```r
print(plot)
```

The current priority figure script should already include this robust saving strategy.

### Missing normalized count columns

`make_priority_figures.R` expects normalized count columns ending in `norm` in the Excel workbook.

If these columns are missing, rerun:

```bash
Rscript scripts/build_excel_from_deseq2_outputs_v2.R
```

### Manual inspection before committing

Before committing changes:

```bash
git status
git diff
git diff --cached
```

Generated outputs should usually not be committed unless they are final deliverables.

## Recommended files to version

Track:

```text
README.md
data/README.md
docs/TUTORIAL.md
scripts/new_de_analysis.R
scripts/build_excel_from_deseq2_outputs_v2.R
scripts/make_priority_figures.R
.gitignore
```

Usually do not track:

```text
DESeq2_redo_clean/
DESeq2_redo_clean_v2/
figures/
results/
Rplots.pdf
```

unless final outputs must be archived in the repository.

## Notes

This workflow was designed for a small RNA-seq comparison. Statistical results and enrichment analyses should be interpreted with caution and together with the experimental context.
