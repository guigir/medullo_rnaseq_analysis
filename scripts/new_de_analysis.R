suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(tibble)
  library(dplyr)
  library(readxl)
})

## ========= PARAMS =========
path_expr <- "data/expr_matrix.txt"
out_dir <- "DESeq2_redo_clean"
alpha          <- 0.05                        # FDR
fc_cut         <- 1.5                         # seuil FC
contrast_ref   <- "HDMB03_GFP_cells"          # condition de référence
contrast_treat <- "cb_parenchyma"             # condition d'intérêt
use_symbol_out <- TRUE
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ========= 1) Lecture de la matrice =========
cts_df <- read.csv(path_expr, sep = ";", stringsAsFactors = FALSE)
gene_ids <- rownames(cts_df)

cts_df[] <- lapply(cts_df, function(x) suppressWarnings(as.numeric(x)))
stopifnot(all(!is.na(colnames(cts_df))))
cts <- as.matrix(cts_df)
rownames(cts) <- make.names(gene_ids, unique = TRUE)
samples <- colnames(cts)

## ========= 2) Métadonnées: condition + batch =========
condition <- case_when(
  grepl("(?i)cb[_ ]?parenchyma", samples) ~ "cb_parenchyma",
  grepl("(?i)HDMB03|GFP", samples)        ~ "HDMB03_GFP_cells",
  TRUE ~ NA_character_
)

batch <- rep(1L, length(samples))
batch[which(samples == "X2_Cb_parenchyma_S5_R1_001")] <- 2L
coldata <- tibble(sample = samples,
                  condition = factor(condition),
                  batch     = factor(batch)) %>%
  column_to_rownames("sample")

coldata$condition <- relevel(coldata$condition, ref = contrast_ref)

### colinéarité batch/condition, on n'utilise que condition comme design
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~condition)

keep <- rowSums(counts(dds) >= 10) >= 2
table(keep)
dds <- dds[keep,]

dds_tx <- rlog(dds, blind = TRUE)
tx_label <- "rlog"

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", contrast_treat, contrast_ref), alpha = alpha)

# Trouver le bon coef pour le shrink 
rn <- resultsNames(dds)
coef_name <- rn[grepl(paste0("^condition_", contrast_treat, "_vs_", contrast_ref, "$"), rn)]
if (length(coef_name) != 1) {
  warning("Coef non trouvé de façon univoque dans resultsNames. Utilisation de 'apeglm' via le contraste directement.")
  resLFC <- lfcShrink(dds, contrast = c("condition", contrast_treat, contrast_ref), type = "apeglm")
} else {
  resLFC <- lfcShrink(dds, coef = coef_name, type = "apeglm")
}

## ========= 6) Table résultats =========
tab <- resLFC %>%
  as.data.frame() %>%
  tibble::rownames_to_column("GENE") %>%
  # Ajout du FoldChange signé avant le FoldChange standard
  dplyr::mutate(
    FoldChange_signed = ifelse(
      is.na(log2FoldChange), NA_real_,
      ifelse(log2FoldChange >= 0,  2^log2FoldChange, -2^(-log2FoldChange))
    ),
    FoldChange = 2^log2FoldChange
  ) %>%
  dplyr::select(GENE, log2FoldChange, FoldChange_signed, FoldChange, baseMean, lfcSE, pvalue, padj)

sym <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = tab$GENE,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

tab$SYMBOL <- sym[tab$GENE]


## ========= 6bisbis) Annotation SYMBOL sans duplication =========

old_annot_file <- "~/Bureau/Guillaume/SOUTH_ROCK/medullo_rnaseq/transfer_10437936_files_78dfbf2c/DE_results_cb_parenchyma_vs_HDMB03_GFP_cells_BATCH_CORRECTED_FULL_ANNOTATIONST.xlsx"
old_annot_sheet <- 1  #

old_annot <- readxl::read_excel(old_annot_file, sheet = old_annot_sheet) %>%
  as_tibble() %>%
  dplyr::select(ensemblID, symbol) %>%
  dplyr::mutate(
    ensemblID_clean = sub("\\.\\d+$", "", as.character(ensemblID)),
    symbol_old = as.character(symbol),
    symbol_old = dplyr::na_if(symbol_old, "")
  ) %>%
  dplyr::filter(!is.na(ensemblID_clean), ensemblID_clean != "") %>%
  dplyr::group_by(ensemblID_clean) %>%
  dplyr::arrange(
    dplyr::desc(!is.na(symbol_old)),
    startsWith(symbol_old, "LOC")
  ) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(GENE_clean = ensemblID_clean, SYMBOL_OLD = symbol_old)

# Fallback org.Hs.eg.db, mais sans créer de doublons
symbol_orgdb <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = sub("\\.\\d+$", "", tab$GENE),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

tab <- tab %>%
  dplyr::mutate(
    GENE_clean = sub("\\.\\d+$", "", as.character(GENE)),
    SYMBOL_ORGDB = unname(symbol_orgdb[GENE_clean])
  ) %>%
  dplyr::left_join(old_annot, by = "GENE_clean") %>%
  dplyr::mutate(
    SYMBOL = dplyr::case_when(
      !is.na(SYMBOL_OLD) & SYMBOL_OLD != "" ~ SYMBOL_OLD,
      !is.na(SYMBOL_ORGDB) & SYMBOL_ORGDB != "" ~ SYMBOL_ORGDB,
      TRUE ~ GENE
    )
  ) %>%
  dplyr::select(-GENE_clean, -SYMBOL_OLD, -SYMBOL_ORGDB) %>%
  dplyr::relocate(SYMBOL, .after = GENE)


tab <- tab %>%
  arrange(padj, desc(abs(log2FoldChange)))

write.csv(tab, file.path(out_dir, "DESeq2_results_shrunk.csv"), row.names = FALSE)


## ========= 7) Listes UP / DOWN selon seuils =========
sym_col <- if (use_symbol_out) "SYMBOL" else "GENE"

UP <- tab %>%
  filter(!is.na(padj), padj < alpha, log2FoldChange >= log2(fc_cut)) %>%
  pull(!!sym_col) %>% unique() %>% na.omit()

DOWN <- tab %>%
  filter(!is.na(padj), padj < alpha, log2FoldChange <= -log2(fc_cut)) %>%
  pull(!!sym_col) %>% unique() %>% na.omit()

writeLines(UP,   file.path(out_dir, paste0("UP_", sym_col, ".txt")))
writeLines(DOWN, file.path(out_dir, paste0("DOWN_", sym_col, ".txt")))

cat("---- DE terminée ----\n",
    "Genes testés: ", nrow(tab), "\n",
    "UP: ", length(UP), "  |  DOWN: ", length(DOWN), "\n",
    "Contraste: ", contrast_treat, " vs ", contrast_ref, "\n",
    "Résultats: ", out_dir, "\n", sep = "")

## ========= 8) PCA =========

## OLD (batch-corrected for comparison with new one) 
## this dds object is BEFORE filtering genes

dds_tx <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~batch+condition)

dds_tx <- rlog(dds_tx, blind = FALSE)
mat <- assay(dds_tx)
mm <- model.matrix(~condition, colData(dds_tx))
mat <- limma::removeBatchEffect(mat, batch=dds_tx$batch, design=mm)
assay(dds_tx) <- mat
p <- plotPCA(dds_tx)
p + theme_bw()

## NEW
rld <- rlog(dds, blind = FALSE)

p <- plotPCA(rld, intgroup = "condition") +
  ggplot2::scale_color_manual(values = c(
    "HDMB03_GFP_cells" = "#00BFC4",
    "cb_parenchyma"    = "#F8766D"
  )) +
  ggplot2::theme_bw()

ggplot2::ggsave(file.path(out_dir, "PCA_rlog_condition.svg"), p, width = 5.5, height = 4.2, dpi = 180)
