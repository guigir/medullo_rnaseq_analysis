# Data directory

This directory stores the input count matrix required by the RNA-seq analysis.

## Required input

```text
expr_matrix.txt
```

## Expected content

`expr_matrix.txt` should be a raw count matrix:

- rows = genes;
- columns = samples;
- row names = Ensembl gene IDs;
- values = raw RNA-seq counts.

Example row names:

```text
ENSG00000000003
ENSG00000000005
ENSG00000000419
```

## Sample name requirements

The scripts infer biological conditions from sample names.

| Pattern in sample name | Assigned condition |
|---|---|
| `cb_parenchyma` | `cb_parenchyma` |
| `HDMB03` or `GFP` | `HDMB03_GFP_cells` |

## Important

Do not provide normalized counts as input to `DESeq2`.

The input matrix must contain raw counts. Normalized counts are generated later by the analysis scripts.
