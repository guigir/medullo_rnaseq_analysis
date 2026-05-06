# Build an Excel report from DESeq2 outputs and raw counts (v2, + FC_signed)
# -------------------------------------------------------------------------
# - Recalcule baseMean depuis la matrice brute (expr_matrix.txt)
# - Ajoute les expressions par échantillon (normalisées & brutes)
# - Ajoute FoldChange_signed = sign(log2FC) * 2^|log2FC| (juste après FoldChange)

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(openxlsx)
  library(readr)
})

# =====================
# 1) Parameters
# =====================
setwd('~/Bureau/Guillaume/SOUTH_ROCK/medullo_rnaseq/medullo_rnaseq_analysis/new_analysis_13_10/')
alpha          <- 0.05
fc_cut         <- 1.5
path_de_csv    <- "../DESeq2_redo_clean_v2/DESeq2_results_shrunk_cleaned_symbol.csv"
path_counts    <- "../data/expr_matrix.txt"
out_xlsx       <- "../DESeq2_redo_clean/DESeq2_results_with_counts.xlsx"

# Ordre souhaité. Les colonnes manquantes seront signalées puis omises.
requested_cols <- c(
  "ensemblID","baseMean","log2FoldChange","lfcSE","pvalue","padj","symbol","FoldChange","FoldChange_signed",
  "1_Cb_parenchyma_S4_R1_001norm","2_Cb_parenchyma_S5_R1_001norm","4_Cb_parenchyma_S6_R1_001norm",
  "A1_HDMB03_S1_R1_001norm","A2_HDMB03_S2_R1_001norm","B1_HDMB03_S3_R1_001norm",
  "1_Cb_parenchyma_S4_R1_001raw","2_Cb_parenchyma_S5_R1_001raw","4_Cb_parenchyma_S6_R1_001raw",
  "A1_HDMB03_S1_R1_001raw","A2_HDMB03_S2_R1_001raw","B1_HDMB03_S3_R1_001raw"
)

# =====================
# 2) Load DE results
# =====================
stopifnot(file.exists(path_de_csv))
DE <- read.csv(path_de_csv, stringsAsFactors = FALSE, check.names = FALSE) %>% as_tibble()

# Harmonisation des noms
DE <- DE %>%
  dplyr::rename_with(~ "ensemblID", .cols = dplyr::any_of("GENE")) %>%
  dplyr::rename_with(~ "symbol", .cols = dplyr::any_of("SYMBOL")) %>%
  dplyr::select(
    dplyr::any_of(c(
      "ensemblID", "symbol", "baseMean", "log2FoldChange",
      "lfcSE", "stat", "pvalue", "padj",
      "FoldChange", "FoldChange_signed"
    ))
  )
# Ajout du FC signé (juste après le FC classique plus tard lors du réordonnancement)
DE <- DE %>%
  mutate(
    FoldChange_signed = dplyr::case_when(
      is.na(log2FoldChange) ~ NA_real_,
      log2FoldChange >= 0   ~  2^log2FoldChange,
      TRUE                  ~ -2^(-log2FoldChange)
    )
  )

# =====================
# 3) Load raw counts & normalized counts (median-of-ratios)
# =====================
stopifnot(file.exists(path_counts))
cts_df <- read.delim(path_counts, sep = ";", stringsAsFactors = FALSE, check.names = FALSE)

# Si la 1ère colonne contient les gènes, l’utiliser comme rownames
if (!is.null(cts_df[[1]]) && !anyDuplicated(cts_df[[1]])) {
  rownames(cts_df) <- make.names(cts_df[[1]], unique = TRUE)
  cts_df[[1]] <- NULL
}

# Passage en numérique
cts_df[] <- lapply(cts_df, function(x) suppressWarnings(as.numeric(x)))
cts <- as.matrix(cts_df)
stopifnot(!any(is.na(cts)))

samples <- colnames(cts)

# colData minimal pour size factors (même logique que dans ton script DE)
condition <- dplyr::case_when(
  grepl("(?i)cb[_ ]?parenchyma", samples) ~ "cb_parenchyma",
  grepl("(?i)HDMB03|GFP", samples)        ~ "HDMB03_GFP_cells",
  TRUE ~ NA_character_
)

batch <- rep(1L, length(samples))
batch[which(samples == "X2_Cb_parenchyma_S5_R1_001")] <- 2L

coldata <- tibble(sample = samples, condition = factor(condition), batch = factor(batch)) %>%
  column_to_rownames("sample")

dds_sf <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
dds_sf <- estimateSizeFactors(dds_sf)

cts_norm <- counts(dds_sf, normalized = TRUE) %>% as.data.frame() %>% rownames_to_column("ensemblID")
cts_raw  <- as.data.frame(cts)                %>% rownames_to_column("ensemblID")

# Nettoyage des noms (enlève un 'X' initial éventuel)
pretty_names <- function(v) sub("^X", "", v)
names(cts_norm)[-1] <- paste0(pretty_names(names(cts_norm)[-1]), "norm")
names(cts_raw)[-1]  <- paste0(pretty_names(names(cts_raw)[-1]),  "raw")

# =====================
# 4) Merge & recompute baseMean from normalized counts
# =====================
base_mean_df <- cts_norm %>%
  mutate(baseMean = apply(dplyr::select(., -ensemblID), 1, mean, na.rm = TRUE)) %>%
  dplyr::select(ensemblID, baseMean)

merged <- DE %>%
  dplyr::select(-any_of("baseMean")) %>%     # on remplace toujours par le recalculé
  left_join(base_mean_df, by = "ensemblID") %>%
  left_join(cts_norm,     by = "ensemblID") %>%
  left_join(cts_raw,      by = "ensemblID")

# Types numériques
num_cols <- intersect(c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","FoldChange","FoldChange_signed"),
                      names(merged))
merged[num_cols] <- lapply(merged[num_cols], function(x) suppressWarnings(as.numeric(x)))

# Réordonnancement : on force l’ordre demandé ; colonnes restantes ensuite
final_cols   <- intersect(requested_cols, names(merged))
missing_cols <- setdiff(requested_cols, names(merged))
if (length(missing_cols)) {
  warning(sprintf("Columns not found and omitted: %s", paste(missing_cols, collapse = ", ")))
}
remaining <- setdiff(names(merged), final_cols)
ALL <- merged %>% dplyr::select(all_of(final_cols), all_of(remaining))

# =====================
# 5) UP / DOWN (selon FC classique, comme avant)
# =====================
UP   <- ALL %>% filter(!is.na(padj), padj < alpha, FoldChange >= fc_cut)   %>% arrange(desc(FoldChange))
DOWN <- ALL %>% filter(!is.na(padj), padj < alpha, FoldChange <= 1/fc_cut) %>% arrange(FoldChange)

# =====================
# 6) Write Excel workbook
# =====================
wb <- createWorkbook()
addWorksheet(wb, "ALL")
addWorksheet(wb, "UP")
addWorksheet(wb, "DOWN")

writeData(wb, sheet = "ALL",  x = ALL)
writeData(wb, sheet = "UP",   x = UP)
writeData(wb, sheet = "DOWN", x = DOWN)

freezePane(wb, "ALL",  firstRow = TRUE)
freezePane(wb, "UP",   firstRow = TRUE)
freezePane(wb, "DOWN", firstRow = TRUE)

setColWidths(wb, "ALL",  cols = 1:10, widths = c(20, 12, 12, 12, 10, 10, 10, 20, 12, 14))
setColWidths(wb, "UP",   cols = 1:10, widths = c(20, 12, 12, 12, 10, 10, 10, 20, 12, 14))
setColWidths(wb, "DOWN", cols = 1:10, widths = c(20, 12, 12, 12, 10, 10, 10, 20, 12, 14))

dir.create(dirname(out_xlsx), recursive = TRUE, showWarnings = FALSE)
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("✅ Excel written: ", normalizePath(out_xlsx))

