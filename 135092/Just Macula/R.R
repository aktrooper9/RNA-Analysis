# ------------------------
# Setup and Libraries
# ------------------------
library(DESeq2)
library(ggplot2)
library(tidyr)
library(apeglm)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ReactomePA)
library(umap)
library(AnnotationDbi)

# Create results directories
dir.create("results_amd_macula", showWarnings = FALSE)
dir.create("results_amd_macula/plots", showWarnings = FALSE)
dir.create("results_amd_macula/tables", showWarnings = FALSE)
dir.create("results_amd_macula/rds_plots", showWarnings = FALSE)

# ------------------------
# Load Metadata and Filter (Tissue + AMD status)
# ------------------------
pdata <- read.csv("sample_metadata_gse135092.csv")
pdata_clean <- subset(pdata, 
                      tissue %in% c("RPE, Macula", "Retina, Macula") & 
                        amd_status %in% c("Control", "AMD") &
                        !is.na(amd_status))
pdata_clean$amd_status <- factor(pdata_clean$amd_status, levels = c("Control", "AMD"))

# ------------------------
# Load and Align Counts
# ------------------------
counts <- read.delim("GSE135092_merged_raw_counts.tsv")
rownames(counts) <- counts$GENEID
counts$GENEID <- NULL

common_samples <- intersect(pdata_clean$sample_id, colnames(counts))
pdata_matched <- pdata_clean[match(common_samples, pdata_clean$sample_id), ]
counts_matched <- counts[, common_samples]
rownames(pdata_matched) <- pdata_matched$sample_id

# ------------------------
# DESeq2 Analysis: AMD vs Control
# ------------------------
dds <- DESeqDataSetFromMatrix(countData = counts_matched,
                              colData = pdata_matched,
                              design = ~ amd_status)

dds <- DESeq(dds)
res_amd <- results(dds, contrast = c("amd_status", "AMD", "Control"))
norm_counts <- counts(dds, normalized = TRUE)

# ------------------------
# DEG Summary
# ------------------------
lfc_thresholds <- c(log2(1.5), log2(2))
fdr_threshold <- 0.05

filter_deg_summary <- function(lfc_cutoff, fdr_cutoff, res) {
  filtered <- res[!is.na(res$padj) & res$padj < fdr_cutoff & abs(res$log2FoldChange) >= lfc_cutoff, ]
  data.frame(
    LFC_Cutoff = lfc_cutoff,
    FDR_Cutoff = fdr_cutoff,
    Upregulated = sum(filtered$log2FoldChange >= lfc_cutoff),
    Downregulated = sum(filtered$log2FoldChange <= -lfc_cutoff),
    Total_DEGs = nrow(filtered)
  )
}

results_summary <- do.call(rbind, lapply(lfc_thresholds, filter_deg_summary, fdr_cutoff = fdr_threshold, res = res_amd))
write.csv(results_summary, "results_amd_macula/tables/deg_summary_amd_vs_control_macula.csv")

# Barplot
df_long <- pivot_longer(results_summary, cols = c("Upregulated", "Downregulated", "Total_DEGs"),
                        names_to = "Category", values_to = "Count")
df_long$LFC_Cutoff <- factor(df_long$LFC_Cutoff, labels = c("1.5 Fold", "2 Fold"))

png("results_amd_macula/plots/deg_counts_plot_amd_vs_control_macula.png", width = 800, height = 600)
ggplot(df_long, aes(x = LFC_Cutoff, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "DEG Counts: AMD vs Control (RPE & Retina Macula)",
       x = "Fold Change Threshold", y = "Number of DEGs") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()
dev.off()

# ------------------------
# Save Results
# ------------------------
write.csv(as.data.frame(res_amd), "results_amd_macula/tables/deseq2_results_amd_vs_control_macula.csv")

res_sig <- res_amd[!is.na(res_amd$padj) & res_amd$padj < fdr_threshold & abs(res_amd$log2FoldChange) > log2(1.5), ]
write.csv(as.data.frame(res_sig), "results_amd_macula/tables/deseq2_sig_results_amd_vs_control_macula.csv")
write.csv(norm_counts, "results_amd_macula/tables/normalized_counts_amd_vs_control_macula.csv")

# ------------------------
# MA + Volcano
# ------------------------
png("results_amd_macula/plots/ma_plot_amd_vs_control_macula.png", width = 800, height = 600)
plotMA(res_amd, ylim = c(-4, 4), main = "MA Plot: AMD vs Control (Macula)")
dev.off()

res_df <- as.data.frame(res_amd)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > log2(1.5), "Yes", "No")

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot: AMD vs Control (Macula)",
       x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal()

ggsave("results_amd_macula/plots/volcano_plot_amd_vs_control_macula.png", volcano_plot, width = 8, height = 6)
saveRDS(volcano_plot, file = "results_amd_macula/rds_plots/volcano_plot_amd_vs_control_macula.rds")

# ------------------------
# Heatmap
# ------------------------
top_genes_heatmap <- if(nrow(res_sig) >= 30) rownames(head(res_sig[order(res_sig$padj), ], 30)) else rownames(res_sig)

if(length(top_genes_heatmap) > 0) {
  top_counts <- norm_counts[top_genes_heatmap, ]
  top_counts_scaled <- t(scale(t(top_counts)))
  
  annotation_col <- as.data.frame(colData(dds)[, "amd_status", drop = FALSE])
  sample_order <- rownames(annotation_col)[order(annotation_col$amd_status)]
  top_counts_sorted <- top_counts_scaled[, sample_order]
  annotation_col_sorted <- annotation_col[sample_order, , drop = FALSE]
  
  breaks <- seq(-2, 2, length.out = 51)
  
  pheatmap(top_counts_sorted,
           breaks = breaks,
           cluster_rows = TRUE, cluster_cols = FALSE,
           annotation_col = annotation_col_sorted,
           show_rownames = TRUE, show_colnames = TRUE,
           fontsize_row = 8, fontsize_col = 9,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           main = paste("Heatmap of Top", length(top_genes_heatmap), "DEGs: AMD vs Control (Macula)"),
           filename = "results_amd_macula/plots/heatmap_top_DEGs_amd_vs_control_macula.png")
  
  saveRDS(top_counts_sorted, "results_amd_macula/rds_plots/heatmap_matrix_top_DEGs_amd_vs_control_macula.rds")
}

# ------------------------
# GO Enrichment
# ------------------------

# Use rownames directly since they are ENTREZIDs
sig_entrez_ids <- rownames(res_sig)
bg_entrez_ids <- rownames(res_amd)

# Perform GO enrichment
go_enrich <- enrichGO(gene          = sig_entrez_ids,
                      universe      = bg_entrez_ids,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "ALL",
                      keyType       = "ENTREZID",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

# Save GO results
write.csv(as.data.frame(go_enrich), "results_amd_macula/tables/go_enrichment_results.csv")

# GO dotplot
png("results_amd_macula/plots/go_dotplot.png", width = 1000, height = 700)
dotplot(go_enrich, showCategory = 20, split = "ONTOLOGY") + facet_wrap(~ONTOLOGY, scales = "free_y")
dev.off()

# Save enrichment object
saveRDS(go_enrich, file = "results_amd_macula/rds_plots/go_enrichment.rds")


# ------------------------
# KEGG Pathway Enrichment (clusterProfiler)
# ------------------------

# KEGG enrichment for significant DEGs
kegg_enrich <- enrichKEGG(
  gene = sig_entrez_ids,
  universe = bg_entrez_ids,
  organism = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)


# Check if any enriched terms were found
if (nrow(as.data.frame(kegg_enrich)) > 0) {
  kegg_enrich_readable <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  png("results_amd_macula/plots/kegg_dotplot.png", width = 1000, height = 700)
  dotplot(kegg_enrich_readable, showCategory = 20) + ggtitle("KEGG Pathway Enrichment")
  dev.off()
} else {
  message("No significant KEGG pathways found.")
}



# Convert KEGG IDs to gene symbols (optional, for readability)
kegg_enrich_readable <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Save KEGG results
write.csv(as.data.frame(kegg_enrich_readable), "results_amd_macula/tables/kegg_enrichment_results.csv")

# KEGG dotplot
png("results_amd_macula/plots/kegg_dotplot.png", width = 1000, height = 700)
dotplot(kegg_enrich_readable, showCategory = 20) + ggtitle("KEGG Pathway Enrichment")
dev.off()

# Save KEGG object
saveRDS(kegg_enrich_readable, file = "results_amd_macula/rds_plots/kegg_enrichment.rds")


# ------------------------
# Reactome Pathway Enrichment
# ------------------------

reactome_enrich <- enrichPathway(gene = sig_entrez_ids,
                                 universe = bg_entrez_ids,
                                 organism = "human",
                                 pvalueCutoff = 0.05,
                                 readable = TRUE)

# Save table
write.csv(as.data.frame(reactome_enrich), "results_amd_macula/tables/reactome_enrichment_results.csv")

# Barplot
png("results_amd_macula/plots/reactome_barplot.png", width = 1000, height = 700)
barplot(reactome_enrich, showCategory = 20, title = "Reactome Enrichment")
dev.off()

# Save object
saveRDS(reactome_enrich, "results_amd_macula/rds_plots/reactome_enrichment.rds")


# Category-gene network plot (cnetplot)
png("results_amd_macula/plots/cnetplot_go.png", width = 1000, height = 800)
cnetplot(reactome_enrich, showCategory = 10, foldChange = res_sig$log2FoldChange[match(sig_entrez_ids, rownames(res_sig))])
dev.off()

# Save
saveRDS(reactome_enrich, "results_amd_macula/rds_plots/go_enrichment_obj.rds")


# ------------------------
# PCA Plot
# ------------------------

vsd <- vst(dds, blind = TRUE)  # variance stabilizing transformation
pcaData <- plotPCA(vsd, intgroup = "amd_status", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = amd_status)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(title = "PCA: AMD vs Control (Macula Only)") +
  theme_minimal()

ggsave("results_amd_macula/plots/pca_plot.png", pca_plot, width = 8, height = 6)
saveRDS(pca_plot, "results_amd_macula/rds_plots/pca_plot.rds")


# ------------------------
# Session Info
# ------------------------
writeLines(capture.output(sessionInfo()), "results_amd_macula/session_info.txt")

cat("\n✅ Analysis complete: AMD vs Control in Macula tissues only\n")
cat("Samples analyzed:", nrow(pdata_matched), "\n")
cat("Group sizes:\n")
print(table(pdata_matched$amd_status))
