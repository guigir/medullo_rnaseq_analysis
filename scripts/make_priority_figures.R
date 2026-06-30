#!/usr/bin/env Rscript

## ============================================================
## Priority RNA-seq figures from DESeq2 output
##
## This script is a script-form extraction of the Quarto report.
## It keeps only the priority figures:
## - GSEA_TOP35_POS_with_large_legends.svg
## - GSEA_TOP30_POS_with_large_legends.svg
## - GSEA_SLIT_ROBO_01_Signaling_by_ROBO_receptors.svg
## - GSEA_SLIT_ROBO_02_Regulation_of_expression_of_SLITs_and_ROBOs.svg
## - Volcano_Gene_Interests_no-legend.svg
## - ComplexUpSet_l2fc_0_585_no-legend.svg
##
## Run from the project root:
##   Rscript scripts/make_priority_figures_from_qmd.R
##
## Optional:
##   Rscript scripts/make_priority_figures_from_qmd.R path/to/DESeq2_results_with_counts.xlsx
## ============================================================

## -------------------------
## Package loading
## -------------------------

required_packages <- c(
  "readxl",
  "dplyr",
  "stringr",
  "tibble",
  "ggplot2",
  "ggrepel",
  "ComplexUpset",
  "clusterProfiler",
  "enrichplot",
  "ReactomePA",
  "org.Hs.eg.db",
  "readr",
  "AnnotationDbi",
  "patchwork",
  "grid"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them first, then rerun this script."
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(ComplexUpset)
  library(clusterProfiler)
  library(enrichplot)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(readr)
  library(AnnotationDbi)
  library(patchwork)
  library(grid)
})

organism_db <- org.Hs.eg.db
plot_font <- "sans"

## -------------------------
## Parameters
## -------------------------

p_thr <- 0.05
primary_l2fc_thr <- 0.585
include_broad_axon_guidance <- FALSE

results_dir <- "results/report"
figures_dir <- "figures/report"

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

## -------------------------
## Input file
## -------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  results_file <- args[[1]]
} else {
  candidate_files <- c(
    "DESeq2_redo_clean/DESeq2_results_with_counts.xlsx",
    "DESeq2_redo_clean_v2/DESeq2_results_with_counts.xlsx"
  )
  existing_candidates <- candidate_files[file.exists(candidate_files)]

  if (length(existing_candidates) == 0) {
    stop(
      "Could not find DESeq2_results_with_counts.xlsx.\n",
      "Expected one of:\n",
      paste(candidate_files, collapse = "\n"),
      "\nOr pass the path explicitly:\n",
      "Rscript scripts/make_priority_figures_from_qmd.R path/to/DESeq2_results_with_counts.xlsx"
    )
  }

  results_file <- existing_candidates[[1]]
}

cat("Using DESeq2 Excel file: ", results_file, "\n", sep = "")

## -------------------------
## Helper functions
## -------------------------

clean_numeric_columns <- function(df, cols) {
  df %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(cols),
        ~ suppressWarnings(as.numeric(stringr::str_replace_all(as.character(.x), ",", ".")))
      )
    )
}

assert_required_columns <- function(df, cols, object_name = deparse(substitute(df))) {
  missing_cols <- setdiff(cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing required columns in ", object_name, ": ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
}

add_de_category <- function(df, p_thr = 0.05, l2fc_thr = 1) {
  df %>%
    dplyr::mutate(
      significant = is.finite(padj) & padj < p_thr,
      de_class = dplyr::case_when(
        significant & log2FoldChange >=  l2fc_thr ~ "Strong up",
        significant & log2FoldChange <= -l2fc_thr ~ "Strong down",
        significant                               ~ "Other sig",
        TRUE                                      ~ "Not sig"
      )
    )
}

make_safe_filename <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("^_+|_+$", "")
}

save_plot_all <- function(plot, name, width, height, units = "cm", dpi = 300, path = figures_dir) {
  stopifnot(inherits(plot, "ggplot") || inherits(plot, "patchwork"))

  dir.create(path, showWarnings = FALSE, recursive = TRUE)

  base_file <- file.path(path, name)
  plot_no_legend <- plot + ggplot2::theme(legend.position = "none")

  ggplot2::ggsave(paste0(base_file, ".svg"), plot, width = width, height = height, units = units)
  ggplot2::ggsave(paste0(base_file, ".png"), plot, width = width, height = height, units = units, dpi = dpi, bg = "white")
  ggplot2::ggsave(paste0(base_file, "_no-legend.svg"), plot_no_legend, width = width, height = height, units = units)
  ggplot2::ggsave(paste0(base_file, "_no-legend.png"), plot_no_legend, width = width, height = height, units = units, dpi = dpi, bg = "white")

  invisible(NULL)
}

save_plot_device <- function(plot, filename_base, width_cm, height_cm, dpi = 300) {
  ## Robust saving for gseaplot2/aplot-like objects.
  ## Direct ggsave() may produce blank SVGs for some enrichplot versions.

  svg_file <- paste0(filename_base, ".svg")
  png_file <- paste0(filename_base, ".png")

  width_in <- width_cm / 2.54
  height_in <- height_cm / 2.54

  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(svg_file, width = width_in, height = height_in)
  } else {
    grDevices::svg(svg_file, width = width_in, height = height_in, family = plot_font)
  }
  print(plot)
  grDevices::dev.off()

  grDevices::png(
    filename = png_file,
    width = width_cm,
    height = height_cm,
    units = "cm",
    res = dpi,
    bg = "white"
  )
  print(plot)
  grDevices::dev.off()

  invisible(NULL)
}

## -------------------------
## Load DESeq2 table
## -------------------------

data_raw <- readxl::read_excel(results_file, sheet = "ALL") %>%
  clean_numeric_columns(c("padj", "log2FoldChange", "FoldChange", "FoldChange_signed"))

assert_required_columns(data_raw, c("ensemblID", "padj", "log2FoldChange"), "data_raw")

if (!"symbol" %in% colnames(data_raw)) {
  data_raw$symbol <- data_raw$ensemblID
}

if ("stat" %in% colnames(data_raw)) {
  data_raw <- clean_numeric_columns(data_raw, "stat")
}

data <- data_raw %>%
  dplyr::mutate(
    ensemblID = as.character(ensemblID),
    ensemblID_clean = sub("\\.\\d+$", "", ensemblID),
    symbol = dplyr::if_else(is.na(symbol) | symbol == "", ensemblID_clean, as.character(symbol)),
    entrez = AnnotationDbi::mapIds(
      organism_db,
      keys = ensemblID_clean,
      keytype = "ENSEMBL",
      column = "ENTREZID",
      multiVals = "first"
    )
  )

data_de <- data %>%
  add_de_category(p_thr = p_thr, l2fc_thr = primary_l2fc_thr)

## ============================================================
## 1) ComplexUpSet_l2fc_0_585_no-legend.svg
## ============================================================

n_total_genes <- nrow(data_de)
n_sig_genes <- sum(data_de$significant, na.rm = TRUE)
sig_pct <- round(100 * n_sig_genes / n_total_genes, 1)

lfc_thr <- primary_l2fc_thr

df_sets_fc <- data_de %>%
  dplyr::mutate(
    StrongUp = significant & log2FoldChange >= lfc_thr,
    StrongDown = significant & log2FoldChange <= -lfc_thr
  ) %>%
  dplyr::filter(significant) %>%
  dplyr::select(StrongUp, StrongDown)

n_up_genes <- sum(df_sets_fc$StrongUp, na.rm = TRUE)
n_down_genes <- sum(df_sets_fc$StrongDown, na.rm = TRUE)

p_upset_fc <- ComplexUpset::upset(
  df_sets_fc,
  intersect = c("StrongUp", "StrongDown"),
  name = "genes",
  width_ratio = 0.15,
  set_sizes = FALSE
) +
  ggplot2::labs(
    title = paste0("Intersection of strong DE classes (|log2FC| >= ", lfc_thr, ")"),
    subtitle = paste0(
      "Total genes tested: ", format(n_total_genes, big.mark = ","),
      " | Significant genes: ", format(n_sig_genes, big.mark = ","),
      " (", sig_pct, "%)",
      " | Upregulated: ", format(n_up_genes, big.mark = ","),
      " | Downregulated: ", format(n_down_genes, big.mark = ",")
    ),
    caption = "Strong classes are defined among significant genes using the selected log2FC threshold."
  ) +
  ggplot2::theme(
    text = ggplot2::element_text(family = plot_font)
  )

save_plot_all(
  p_upset_fc,
  paste0("ComplexUpSet_l2fc_", gsub("\\.", "_", as.character(lfc_thr))),
  width = 18,
  height = 10,
  units = "cm"
)

## ============================================================
## 2) Volcano_Gene_Interests_no-legend.svg
## ============================================================

interest_genes_axon <- c("PLXNA1", "PLXNA4", "CXCR4", "SEMA6B", "SEMA3B", "NRP2", "SEMA4B", "SEMA6C", "ROBO2")
interest_genes_g3 <- c("MYC", "JUN")

interest_df <- data_de %>%
  dplyr::mutate(
    category = dplyr::case_when(
      symbol %in% interest_genes_axon ~ "Axon genes",
      symbol %in% interest_genes_g3 ~ "G3 genes",
      TRUE ~ "Other"
    )
  )

interest_palette <- c(
  "Axon genes" = "coral2",
  "G3 genes" = "gold",
  "Other" = "grey80"
)

interest_yvals <- -log10(pmax(interest_df$padj, 1e-300))
interest_ymax <- max(interest_yvals[is.finite(interest_yvals)], na.rm = TRUE)

p_volcano_interest <- ggplot2::ggplot(
  interest_df,
  ggplot2::aes(x = log2FoldChange, y = -log10(padj), fill = category)
) +
  ggplot2::geom_point(shape = 21, color = "white", size = 1.6, alpha = 0.7, stroke = 0.2) +
  ggrepel::geom_label_repel(
    data = interest_df %>% dplyr::filter(category != "Other"),
    ggplot2::aes(label = symbol),
    size = 3.4,
    max.overlaps = 200,
    box.padding = 0.45,
    point.padding = 0.25,
    segment.size = 0.3,
    min.segment.length = 0,
    label.size = 0.15,
    seed = 42
  ) +
  ggplot2::geom_hline(yintercept = -log10(p_thr), linetype = "dashed") +
  ggplot2::geom_vline(xintercept = c(-primary_l2fc_thr, primary_l2fc_thr), linetype = "dashed") +
  ggplot2::coord_cartesian(ylim = c(-0.6, interest_ymax * 1.02), clip = "off") +
  ggplot2::scale_fill_manual(values = interest_palette, name = NULL) +
  ggplot2::labs(
    x = expression(Log[2] ~ "Fold Change"),
    y = expression(-Log[10] ~ "adjusted p-value"),
    title = "Volcano plot — selected axon and G3 genes"
  ) +
  ggplot2::theme_bw(base_size = 14, base_family = plot_font) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(20, 24, 50, 24)
  )

save_plot_all(p_volcano_interest, "Volcano_Gene_Interests", width = 32, height = 22, units = "cm")

## ============================================================
## 3) Reactome GSEA
## ============================================================

data_de_gsea <- data_de %>%
  dplyr::mutate(
    ensemblID = as.character(ensemblID),
    ensemblID_clean = sub("\\.\\d+$", "", ensemblID)
  ) %>%
  dplyr::distinct(ensemblID_clean, .keep_all = TRUE)

if ("stat" %in% colnames(data_de_gsea) && any(is.finite(data_de_gsea$stat), na.rm = TRUE)) {
  data_de_gsea <- data_de_gsea %>%
    dplyr::mutate(rank_score = stat)
} else {
  data_de_gsea <- data_de_gsea %>%
    dplyr::mutate(
      padj_safe = pmax(padj, 1e-300),
      rank_score = sign(log2FoldChange) * -log10(padj_safe)
    )
}

ensembl_to_entrez <- AnnotationDbi::select(
  organism_db,
  keys = unique(data_de_gsea$ensemblID_clean),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME")
) %>%
  tibble::as_tibble() %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  dplyr::distinct(ENSEMBL, ENTREZID, .keep_all = TRUE)

ensembl_to_entrez_one <- ensembl_to_entrez %>%
  dplyr::group_by(ENSEMBL) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

gsea_entrez_df <- data_de_gsea %>%
  dplyr::left_join(
    ensembl_to_entrez_one,
    by = c("ensemblID_clean" = "ENSEMBL")
  ) %>%
  dplyr::filter(
    !is.na(ENTREZID),
    !is.na(rank_score),
    is.finite(rank_score)
  ) %>%
  dplyr::transmute(
    ensemblID = ensemblID_clean,
    symbol_input = symbol,
    ENTREZID = as.character(ENTREZID),
    SYMBOL_orgdb = SYMBOL,
    GENENAME = GENENAME,
    rank_score = rank_score
  ) %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::slice_max(order_by = abs(rank_score), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

ranked_entrez <- gsea_entrez_df$rank_score
names(ranked_entrez) <- gsea_entrez_df$ENTREZID
ranked_entrez <- sort(ranked_entrez, decreasing = TRUE)

readr::write_csv(
  gsea_entrez_df,
  file.path(results_dir, "GSEA_ranked_genes_ENTREZ_clean_ENSEMBL_to_ENTREZ.csv")
)

write.table(
  ranked_entrez,
  file.path(results_dir, "ranks_entrez_clean_ENSEMBL_to_ENTREZ.txt"),
  quote = FALSE,
  sep = "\t"
)

cat("Number of unique Ensembl in data_de: ", dplyr::n_distinct(data_de_gsea$ensemblID_clean), "\n", sep = "")
cat("Number of Entrez in ranked vector: ", length(ranked_entrez), "\n", sep = "")
cat("Duplicated Entrez in ranked vector: ", anyDuplicated(names(ranked_entrez)), "\n", sep = "")

set.seed(123)

gsea_reactome <- ReactomePA::gsePathway(
  ranked_entrez,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 5,
  maxGSSize = 500,
  by = "fgsea",
  eps = 0
)

gsea_reactome <- clusterProfiler::setReadable(
  gsea_reactome,
  OrgDb = organism_db,
  keyType = "ENTREZID"
)

gsea_reactome_tbl <- tibble::as_tibble(gsea_reactome@result) %>%
  dplyr::arrange(dplyr::desc(NES))

readr::write_csv(
  gsea_reactome_tbl,
  file.path(results_dir, "GSEA_Reactome_full_table.csv")
)

gsea_reactome_sig_tbl <- gsea_reactome_tbl %>%
  dplyr::filter(!is.na(p.adjust), p.adjust < 0.05) %>%
  dplyr::arrange(desc(NES))

gsea_reactome_top_tbl <- gsea_reactome_tbl %>%
  dplyr::filter(!is.na(p.adjust)) %>%
  dplyr::arrange(desc(NES)) %>%
  dplyr::select(
    ID, Description, setSize, enrichmentScore, NES,
    pvalue, p.adjust, qvalue, rank, leading_edge, core_enrichment
  )

readr::write_csv(gsea_reactome_sig_tbl, file.path(results_dir, "GSEA_Reactome_significant_terms.csv"))
readr::write_csv(gsea_reactome_top_tbl, file.path(results_dir, "GSEA_Reactome_compact_table.csv"))

## ============================================================
## 4) GSEA_SLIT_ROBO_01 / 02 robust curves
## ============================================================

slit_robo_core_ids <- c(
  "R-HSA-376176",   # Signaling by ROBO receptors
  "R-HSA-9010553",  # Regulation of expression of SLITs and ROBOs
  "R-HSA-428542"    # Regulation of commissural axon pathfinding by SLIT and ROBO
)

slit_robo_regex <- paste(
  c("SLIT", "ROBO", "Roundabout", "ROBO-SLIT", "SLIT and ROBO"),
  collapse = "|"
)

slit_robo_tbl <- gsea_reactome_tbl %>%
  dplyr::filter(
    ID %in% slit_robo_core_ids |
      grepl(slit_robo_regex, Description, ignore.case = TRUE)
  )

if (include_broad_axon_guidance) {
  slit_robo_tbl <- gsea_reactome_tbl %>%
    dplyr::filter(
      ID %in% slit_robo_core_ids |
        grepl(slit_robo_regex, Description, ignore.case = TRUE) |
        grepl("Axon guidance", Description, ignore.case = TRUE)
    )
}

slit_robo_tbl <- slit_robo_tbl %>%
  dplyr::distinct(ID, .keep_all = TRUE) %>%
  dplyr::arrange(desc(NES), p.adjust, Description)

readr::write_csv(
  slit_robo_tbl,
  file.path(results_dir, "GSEA_Reactome_SLIT_ROBO_pathways.csv")
)

if (nrow(slit_robo_tbl) == 0) {
  warning("No SLIT/ROBO Reactome terms were found in the GSEA result table.")
}

plot_gsea_curve_robust <- function(gsea_obj, result_tbl, pathway_id) {
  pathway_info <- result_tbl %>%
    dplyr::filter(ID == pathway_id) %>%
    dplyr::slice(1)

  if (nrow(pathway_info) != 1) {
    stop("Could not find pathway ID in result table: ", pathway_id)
  }

  pathway_idx <- match(pathway_id, gsea_obj@result$ID)

  if (is.na(pathway_idx)) {
    stop("Could not find pathway ID in gsea object: ", pathway_id)
  }

  ## Index-based geneSetID is the most robust with ReactomePA/enrichplot.
  p <- enrichplot::gseaplot2(
    gsea_obj,
    geneSetID = pathway_idx,
    title = paste0(
      pathway_info$Description,
      "\nNES = ", round(pathway_info$NES, 3),
      " | padj = ", signif(pathway_info$p.adjust, 3),
      " | set size = ", pathway_info$setSize
    )
  )

  p
}

if (nrow(slit_robo_tbl) > 0) {
  for (i in seq_len(nrow(slit_robo_tbl))) {

    pathway_id <- as.character(slit_robo_tbl$ID[i])
    pathway_desc <- as.character(slit_robo_tbl$Description[i])

    p_gsea_slit_robo <- plot_gsea_curve_robust(
      gsea_obj = gsea_reactome,
      result_tbl = gsea_reactome_tbl,
      pathway_id = pathway_id
    )

    file_suffix <- paste0(
      sprintf("%02d", i), "_",
      make_safe_filename(pathway_desc)
    )

    filename_base <- file.path(figures_dir, paste0("GSEA_SLIT_ROBO_", file_suffix))

    save_plot_device(
      plot = p_gsea_slit_robo,
      filename_base = filename_base,
      width_cm = 13,
      height_cm = 10,
      dpi = 300
    )

    cat("Saved SLIT/ROBO GSEA curve: ", filename_base, ".svg\n", sep = "")
  }
}

## ============================================================
## 5) GSEA_TOP30_POS_with_large_legends.svg
##    GSEA_TOP35_POS_with_large_legends.svg
## ============================================================

col_axon     <- rgb(91, 45, 142, maxColorValue = 255)
col_slitrobo <- rgb(204, 15, 0, maxColorValue = 255)
col_other    <- rgb(136, 136, 136, maxColorValue = 255)
col_transl   <- rgb(166, 125, 11, maxColorValue = 255)

cat_palette <- c(
  "Axon guidance" = col_axon,
  "SLIT/ROBO"     = col_slitrobo,
  "Translation"   = col_transl,
  "Other"         = col_other
)

manual_category_map <- c(
  "rRNA processing" = "Translation",
  "Translation" = "Translation",
  "rRNA processing in the nucleus and cytosol" = "Translation",
  "Major pathway of rRNA processing in the nucleolus and cytosol" = "Translation",
  "Metabolism of amino acids and derivatives" = "Translation",
  "Signaling by ROBO receptors" = "SLIT/ROBO",
  "Influenza Viral RNA Transcription and Replication" = "Other",
  "Influenza Infection" = "Other",
  "Cellular response to starvation" = "Other",
  "Regulation of expression of SLITs and ROBOs" = "SLIT/ROBO",
  "Ribosome-associated quality control" = "Translation",
  "Eukaryotic Translation Initiation" = "Translation",
  "Cap-dependent Translation Initiation" = "Translation",
  "SRP-dependent cotranslational protein targeting to membrane" = "Translation",
  "Selenoamino acid metabolism" = "Other",
  "Nonsense-Mediated Decay (NMD)" = "Translation",
  "Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC)" = "Translation",
  "GTP hydrolysis and joining of the 60S ribosomal subunit" = "Translation",
  "L13a-mediated translational silencing of Ceruloplasmin expression" = "Translation",
  "Formation of a pool of free 40S subunits" = "Translation",
  "Response of EIF2AK4 (GCN2) to amino acid deficiency" = "Translation",
  "ZNF598 and the Ribosome-associated Quality Trigger (RQT) complex dissociate a ribosome stalled on a no-go mRNA" = "Translation",
  "PELO:HBS1L and ABCE1 dissociate a ribosome on a non-stop mRNA" = "Translation",
  "Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC)" = "Translation",
  "Peptide chain elongation" = "Translation",
  "Selenocysteine synthesis" = "Other",
  "Viral mRNA Translation" = "Other",
  "Eukaryotic Translation Termination" = "Translation",
  "Eukaryotic Translation Elongation" = "Translation",
  "SARS-CoV-1 Infection" = "Other",
  "Axon guidance" = "Axon guidance",
  "SARS-CoV-1-host interactions" = "Other",
  "Ribosome Quality Control (RQC) complex extracts and degrades nascent peptide" = "Translation",
  "SARS-CoV-2-host interactions" = "Other",
  "SARS-CoV Infections" = "Other"
)

assign_category <- function(description) {
  category <- unname(manual_category_map[description])

  missing_desc <- description[is.na(category)]

  if (length(missing_desc) > 0) {
    warning(
      paste0(
        "Some pathways are missing from manual_category_map and will be labelled as Other:\n",
        paste(unique(missing_desc), collapse = "\n")
      )
    )
    category[is.na(category)] <- "Other"
  }

  category
}

prepare_top_reactome <- function(gsea_tbl, n_top = 30) {
  gsea_tbl %>%
    dplyr::filter(!is.na(NES), !is.na(p.adjust), NES > 0) %>%
    dplyr::arrange(dplyr::desc(NES)) %>%
    dplyr::slice_head(n = n_top) %>%
    dplyr::mutate(
      Description = as.character(Description),
      Category = assign_category(Description),
      Category = factor(
        Category,
        levels = c("Translation", "SLIT/ROBO", "Axon guidance", "Other")
      ),
      Description = factor(Description, levels = rev(Description))
    )
}

make_reactome_top_plot <- function(top_tbl, title = NULL) {

  axis_label_colors <- setNames(
    cat_palette[as.character(top_tbl$Category)],
    as.character(top_tbl$Description)
  )

  axis_label_colors <- axis_label_colors[levels(top_tbl$Description)]

  p_main <- ggplot2::ggplot(
    top_tbl,
    ggplot2::aes(
      x = Description,
      y = NES,
      color = p.adjust
    )
  ) +
    ggplot2::geom_point(size = 3.2) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::scale_color_distiller(
      type = "seq",
      direction = -1,
      values = c(0, 0.05, 0.8, 1),
      palette = "OrRd",
      name = "Adjusted p-value"
    ) +
    ggplot2::theme_classic(base_family = plot_font) +
    ggplot2::theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c("left", "top"),
      legend.background = ggplot2::element_rect(
        fill = "grey90",
        colour = "grey40",
        linewidth = 0.3
      ),
      legend.key.width = grid::unit(1.2, "lines"),
      legend.key.height = grid::unit(0.6, "lines"),
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 10),
      panel.background = ggplot2::element_rect(fill = "grey90"),
      axis.text.y = ggplot2::element_text(
        size = 12,
        face = "bold",
        colour = axis_label_colors,
        margin = ggplot2::margin(r = 10)
      ),
      axis.text.x = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(
        size = 15,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 12)
      ),
      plot.margin = ggplot2::margin(20, 30, 20, 30)
    ) +
    ggplot2::ylab("NES") +
    ggplot2::xlab("")

  if (!is.null(title)) {
    p_main <- p_main + ggplot2::ggtitle(title)
  }

  legend_cat_df <- tibble::tibble(
    Category = factor(
      c("Translation", "SLIT/ROBO", "Axon guidance", "Other"),
      levels = c("Translation", "SLIT/ROBO", "Axon guidance", "Other")
    ),
    x = c(0.03, 0.33, 0.62, 0.88),
    y = c(0.50, 0.50, 0.50, 0.50),
    label = c(
      "Translation-related pathways",
      "SLIT/ROBO-related pathways",
      "Axon guidance pathway",
      "Other pathways"
    )
  )

  p_cat <- ggplot2::ggplot(
    legend_cat_df,
    ggplot2::aes(x = x, y = y)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = Category),
      size = 6
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      hjust = 0,
      nudge_x = 0.025,
      family = plot_font,
      size = 13 / 2.845
    ) +
    ggplot2::scale_color_manual(
      values = cat_palette,
      guide = "none"
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, 1.25),
      ylim = c(0, 1),
      clip = "off"
    ) +
    ggplot2::theme_void(base_family = plot_font) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 80, 20, 80)
    )

  p_main / p_cat +
    patchwork::plot_layout(heights = c(11, 2.3))
}

top30_reactome <- prepare_top_reactome(gsea_reactome_tbl, n_top = 30)
top35_reactome <- prepare_top_reactome(gsea_reactome_tbl, n_top = 35)

readr::write_csv(
  top30_reactome %>%
    dplyr::select(ID, Description, NES, p.adjust, setSize, Category),
  file.path(results_dir, "Reactome_GSEA_TOP30_positive_category_mapping.csv")
)

readr::write_csv(
  top35_reactome %>%
    dplyr::select(ID, Description, NES, p.adjust, setSize, Category),
  file.path(results_dir, "Reactome_GSEA_TOP35_positive_category_mapping.csv")
)

plot_top30 <- make_reactome_top_plot(
  top30_reactome,
  title = "Reactome GSEA — top 30 positive pathways"
)

plot_top35 <- make_reactome_top_plot(
  top35_reactome,
  title = "Reactome GSEA — top 35 positive pathways"
)

ggplot2::ggsave(
  filename = file.path(figures_dir, "GSEA_TOP30_POS_with_large_legends.svg"),
  plot = plot_top30,
  width = 36,
  height = 28,
  units = "cm"
)

ggplot2::ggsave(
  filename = file.path(figures_dir, "GSEA_TOP30_POS_with_large_legends.png"),
  plot = plot_top30,
  width = 36,
  height = 28,
  units = "cm",
  dpi = 300,
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(figures_dir, "GSEA_TOP35_POS_with_large_legends.svg"),
  plot = plot_top35,
  width = 38,
  height = 32,
  units = "cm"
)

ggplot2::ggsave(
  filename = file.path(figures_dir, "GSEA_TOP35_POS_with_large_legends.png"),
  plot = plot_top35,
  width = 38,
  height = 32,
  units = "cm",
  dpi = 300,
  bg = "white"
)

cat("\nPriority figure generation completed.\n")
cat("Figures written to: ", figures_dir, "\n", sep = "")
cat("Results written to: ", results_dir, "\n", sep = "")
