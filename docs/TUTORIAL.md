# Tutorial — reproduce the RNA-seq analysis

This tutorial explains how to reproduce the analysis from a clean clone of the repository.

## 1. Clone the repository

```bash
git clone https://github.com/guigir/medullo_rnaseq_analysis.git
cd medullo_rnaseq_analysis
```

## 2. Add the input count matrix

Place the raw count matrix here:

```text
data/expr_matrix.txt
```

The matrix must contain raw counts, not normalized values.

Expected format:

```text
            sample_1    sample_2    sample_3
ENSG000...  123         456         789
ENSG000...  0           12          35
```

The row names should be Ensembl gene IDs.

The column names should contain enough information to assign conditions:

- `cb_parenchyma` → `cb_parenchyma`;
- `HDMB03` or `GFP` → `HDMB03_GFP_cells`.

## 3. Check the R environment

Run:

```bash
which R
which Rscript
Rscript -e 'R.version.string'
Rscript -e '.libPaths()'
```

If the output points to a conda R installation and this is not intended, run:

```bash
conda deactivate
```

Then check again.

## 4. Install dependencies

Open R and install:

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

## 5. Run the DESeq2 analysis

From the repository root:

```bash
Rscript scripts/new_de_analysis.R
```

This performs:

1. count matrix loading;
2. sample condition inference;
3. DESeq2 object construction;
4. low-count gene filtering;
5. DESeq2 model fitting;
6. log2FC shrinkage with `apeglm`;
7. gene symbol annotation;
8. export of UP and DOWN gene lists.

Expected files:

```text
DESeq2_redo_clean_v2/DESeq2_results_shrunk_cleaned_symbol.csv
DESeq2_redo_clean_v2/UP_SYMBOL.txt
DESeq2_redo_clean_v2/DOWN_SYMBOL.txt
```

## 6. Build the Excel workbook

Run:

```bash
Rscript scripts/build_excel_from_deseq2_outputs_v2.R
```

This creates a workbook containing:

- all genes;
- UP genes;
- DOWN genes;
- raw counts;
- normalized counts;
- DESeq2 statistics.

Expected file:

```text
DESeq2_redo_clean/DESeq2_results_with_counts.xlsx
```

or:

```text
DESeq2_redo_clean_v2/DESeq2_results_with_counts.xlsx
```

## 7. Generate priority figures

Run:

```bash
Rscript scripts/make_priority_figures.R
```

or provide the Excel file explicitly:

```bash
Rscript scripts/make_priority_figures.R DESeq2_redo_clean_v2/DESeq2_results_with_counts.xlsx
```

This generates:

```text
figures/report/ComplexUpSet_l2fc_0_585_no-legend.svg
figures/report/Volcano_Gene_Interests_no-legend.svg
figures/report/GSEA_SLIT_ROBO_01_Signaling_by_ROBO_receptors.svg
figures/report/GSEA_SLIT_ROBO_02_Regulation_of_expression_of_SLITs_and_ROBOs.svg
figures/report/GSEA_TOP30_POS_with_large_legends.svg
figures/report/GSEA_TOP35_POS_with_large_legends.svg
```

## 8. Check outputs

```bash
ls -lh DESeq2_redo_clean_v2/
ls -lh DESeq2_redo_clean/
ls -lh figures/report/
ls -lh results/report/
```

Check specifically:

```bash
ls -lh figures/report/*SLIT_ROBO*.svg
ls -lh figures/report/*TOP30*.svg
ls -lh figures/report/*TOP35*.svg
ls -lh figures/report/*Volcano_Gene_Interests*.svg
ls -lh figures/report/*ComplexUpSet_l2fc_0_585*.svg
```

## 9. Optional validation checks

Inspect the Reactome GSEA table:

```bash
head results/report/GSEA_Reactome_full_table.csv
```

Inspect SLIT/ROBO pathways:

```bash
cat results/report/GSEA_Reactome_SLIT_ROBO_pathways.csv
```

Inspect the ranked gene table:

```bash
head results/report/GSEA_ranked_genes_ENTREZ_clean_ENSEMBL_to_ENTREZ.csv
```

## 10. Commit only source files

Recommended:

```bash
git add README.md data/README.md docs/TUTORIAL.md scripts/
git status
git diff --cached
git commit -m "Add reproducible RNA-seq workflow documentation"
git push
```

Generated outputs are usually ignored unless they are intended final deliverables.
