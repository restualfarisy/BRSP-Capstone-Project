# ============================================================
# ANALISIS DIFFERENTIAL GENE EXPRESSION (DEG) - GSE33177
# Organisme : Solanum lycopersicum (Tomat)
# Platform  : GPL4741 | Affymetrix Tomato Genome Array
# Tujuan    : Mengidentifikasi gen yang berubah ekspresinya
#             akibat infeksi Phytophthora infestans (Late Blight)
# ============================================================


# ============================================================
# STEP 1: INSTALL PACKAGE (jalankan SEKALI saja)
# ============================================================
install.packages(c("ggplot2", "ggrepel", "pheatmap", "reshape2", "BiocManager"))

BiocManager::install(c(
  "GEOquery",
  "limma",
  "clusterProfiler",
  "enrichplot",
  "org.At.tair.db",   # Arabidopsis (referensi utama yang berhasil)
  "org.Sl.eg.db"      # Tomat (coba dulu, fallback ke Arabidopsis)
))


# ============================================================
# STEP 2: LOAD LIBRARY
# ============================================================
library(GEOquery)
library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(reshape2)
library(clusterProfiler)
library(enrichplot)
library(org.At.tair.db)


# ============================================================
# STEP 3: LOAD & PREPROCESSING DATA
# ============================================================
options(timeout = 600)
gset <- getGEO("GSE33177", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

# Ambil matriks ekspresi
ex <- exprs(gset)

# Cek apakah data perlu log2 transform
# qx[5] > 100 : nilai ke-99% masih sangat besar → belum di-log
# qx[6]-qx[1] > 50 : rentang data sangat lebar → perlu dikompres
qx     <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
LogC   <- (qx[5] > 100) || (qx[6] - qx[1] > 50)

if (LogC) {
  ex[ex <= 0] <- NA
  ex          <- log2(ex)
  message("✅ Data berhasil di-log2 transform")
} else {
  message("ℹ️  Data tidak perlu di-log2 transform")
}


# ============================================================
# STEP 4: DEFINISI GRUP SAMPEL
# ============================================================
# GSE33177: 3 sampel Healthy + 4 sampel LateBlight
# Urutan: GSM821376, GSM821377, GSM821378 = Healthy
#         GSM821379, GSM821380, GSM821381, GSM821382 = LateBlight
gsms  <- "0001111"
sml   <- strsplit(gsms, split = "")[[1]]
group <- factor(sml, labels = c("Healthy", "LateBlight"))

pData(gset)$group <- group
print(table(pData(gset)$group))


# ============================================================
# STEP 5: ANALISIS DEG DENGAN LIMMA
# ============================================================

# --- 5a. Design Matrix ---
# ~0 berarti TANPA intercept (best practice limma)
design <- model.matrix(~0 + gset$group)
colnames(design) <- levels(gset$group)

# --- 5b. Kontras: LateBlight vs Healthy ---
contrast_formula <- "LateBlight - Healthy"
contrast_matrix  <- makeContrasts(contrasts = contrast_formula, levels = design)
print(paste("Kontras yang dianalisis:", contrast_formula))

# --- 5c. Fit Model ---
fit  <- lmFit(ex, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# --- 5d. Ambil Hasil DEG ---
# topTableResults: gen signifikan saja (FDR < 0.05)
topTableResults <- topTable(
  fit2,
  adjust  = "fdr",
  sort.by = "B",
  number  = Inf,
  p.value = 0.05
)

# allResults: semua gen tanpa filter (untuk heatmap & enrichment)
allResults <- topTable(
  fit2,
  adjust  = "fdr",
  sort.by = "p",
  number  = Inf
)

print(paste("Total DEG signifikan (FDR < 0.05):", nrow(topTableResults)))


# ============================================================
# STEP 6: ANOTASI GEN
# ============================================================
# Karena AnnotGPL = TRUE saat load data, anotasi sudah tersedia
# langsung di fData(gset) tanpa perlu package tambahan

gene_annotation <- fData(gset)[, c("ID", "Gene symbol", "Gene title")]

# Anotasi topTableResults
topTableResults$PROBEID <- rownames(topTableResults)
topTableResults <- merge(
  topTableResults, gene_annotation,
  by.x = "PROBEID", by.y = "ID", all.x = TRUE
)

# Anotasi allResults
allResults$PROBEID <- rownames(allResults)
allResults <- merge(
  allResults, gene_annotation,
  by.x = "PROBEID", by.y = "ID", all.x = TRUE
)

print(paste("Total DEG:", nrow(topTableResults)))
print(paste("Berhasil dianotasi:", sum(!is.na(topTableResults$"Gene symbol"))))


# ============================================================
# STEP 7: VISUALISASI
# ============================================================

# --- 7a. Boxplot Distribusi Ekspresi ---
group_colors <- ifelse(gset$group == "Healthy", "lightgreen", "tomato")
sample_labels <- paste0(gset$group, "_rep",
                        ave(seq_along(gset$group), gset$group, FUN = seq_along))

boxplot(
  ex,
  col     = group_colors,
  las     = 2,
  outline = FALSE,
  names   = sample_labels,
  main    = "Boxplot Distribusi Nilai Ekspresi per Sampel\nTomat Sehat vs Terinfeksi Phytophthora infestans",
  ylab    = "Expression Value (log2)",
  xlab    = ""
)
legend("topright",
       legend = c("Healthy", "LateBlight"),
       fill   = c("lightgreen", "tomato"),
       title  = "Grup Sampel",
       cex    = 0.8)

# --- 7b. Density Plot ---
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group      = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  scale_color_manual(
    values = c("Healthy" = "lightgreen", "LateBlight" = "tomato"),
    labels = c("Healthy" = "Sehat", "LateBlight" = "Terinfeksi P. infestans")
  ) +
  theme_minimal() +
  labs(
    title    = "Distribusi Nilai Ekspresi Gen\nTomat Sehat vs Terinfeksi Phytophthora infestans",
    subtitle = "Dataset GSE33177 | Platform: Affymetrix Tomato Genome Array",
    x        = "Expression Value (log2)",
    y        = "Density",
    color    = "Grup Sampel"
  )

# --- 7c. PCA Plot ---
pca_result <- prcomp(t(ex), scale. = TRUE, center = TRUE)
pca_var    <- round(summary(pca_result)$importance[2, ] * 100, 1)

pca_df <- data.frame(
  PC1    = pca_result$x[, 1],
  PC2    = pca_result$x[, 2],
  Group  = gset$group,
  Sample = rownames(pca_result$x)
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(vjust = -0.8, size = 3) +
  scale_color_manual(
    values = c("Healthy" = "lightgreen", "LateBlight" = "tomato"),
    labels = c("Healthy" = "Sehat", "LateBlight" = "Terinfeksi P. infestans")
  ) +
  theme_minimal() +
  labs(
    title    = "PCA Plot: Tomat Sehat vs Terinfeksi Phytophthora infestans",
    subtitle = "Dataset GSE33177 | 7 Sampel | Affymetrix Tomato Genome Array",
    x        = paste0("PC1 (", pca_var[1], "% variansi)"),
    y        = paste0("PC2 (", pca_var[2], "% variansi)"),
    color    = "Grup Sampel"
  )

# --- 7d. Bar Chart Jumlah DEG ---
summary_deg <- data.frame(
  Comparison = "LateBlight vs Healthy",
  Up         = sum(topTableResults$logFC > 0, na.rm = TRUE),
  Down       = sum(topTableResults$logFC < 0, na.rm = TRUE)
)

summary_long          <- melt(summary_deg[, c("Comparison", "Up", "Down")], id.vars = "Comparison")
colnames(summary_long) <- c("Comparison", "Regulation", "Count")

ggplot(summary_long, aes(x = Comparison, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Up" = "tomato", "Down" = "steelblue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title    = "Jumlah DEG: Tomat Terinfeksi vs Sehat",
    subtitle = "Dataset GSE33177 | FDR < 0.05",
    x        = "",
    y        = "Jumlah Gen"
  )

# --- 7e. Volcano Plot ---
volcano_data <- data.frame(
  logFC     = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene      = topTableResults$"Gene symbol"
)

volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC >  1 & volcano_data$adj.P.Val < 0.05] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.05] <- "DOWN"
print(table(volcano_data$status))

# Label top 10 gen paling signifikan
top_genes <- head(
  volcano_data[volcano_data$status != "NO" & !is.na(volcano_data$Gene), ],
  10
)

ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c("DOWN" = "steelblue", "NO" = "grey70", "UP" = "tomato"),
    labels = c("DOWN" = "Downregulated", "NO" = "Not Significant", "UP" = "Upregulated")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes, aes(label = Gene),
                  size = 3, color = "black", max.overlaps = 20) +
  theme_minimal() +
  labs(
    title    = "Volcano Plot DEG: Tomat Terinfeksi vs Sehat",
    subtitle = "Dataset GSE33177 | FDR < 0.05 | |logFC| > 1",
    x        = "Log2 Fold Change",
    y        = "-log10(Adjusted P-Value)",
    color    = "Status Gen"
  )

# --- 7f. Heatmap Top 50 DEG ---
# Gunakan allResults agar heatmap tidak kosong
top50      <- head(allResults[order(allResults$P.Value), ], 50)
mat_heatmap <- ex[top50$PROBEID, ]

gene_label <- ifelse(
  is.na(top50$"Gene symbol") | top50$"Gene symbol" == "",
  top50$PROBEID,
  top50$"Gene symbol"
)
rownames(mat_heatmap) <- make.unique(gene_label)

# Bersihkan data (wajib agar tidak error saat clustering)
mat_heatmap <- mat_heatmap[rowSums(is.na(mat_heatmap)) == 0, ]
mat_heatmap <- mat_heatmap[apply(mat_heatmap, 1, var) > 0, ]

annotation_col    <- data.frame(Group = gset$group)
rownames(annotation_col) <- colnames(mat_heatmap)
annotation_colors <- list(Group = c("Healthy" = "lightgreen", "LateBlight" = "tomato"))

pheatmap(
  mat_heatmap,
  scale                    = "row",
  annotation_col           = annotation_col,
  annotation_colors        = annotation_colors,
  show_colnames            = TRUE,
  show_rownames            = TRUE,
  fontsize_row             = 7,
  fontsize_col             = 8,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete",
  color                    = colorRampPalette(c("steelblue", "white", "tomato"))(100),
  main                     = "Top 50 DEG: Tomat Terinfeksi P. infestans vs Sehat\nDataset GSE33177"
)


# ============================================================
# STEP 8: GO ENRICHMENT ANALYSIS
# ============================================================

# Ambil gen signifikan (P.Value < 0.05) dan bersihkan
deg_genes  <- allResults[allResults$P.Value < 0.05, ]
genes_all  <- deg_genes$"Gene symbol"
genes_all  <- genes_all[!is.na(genes_all) & genes_all != ""]
genes_up   <- deg_genes[deg_genes$logFC > 0, "Gene symbol"]
genes_up   <- genes_up[!is.na(genes_up) & genes_up != ""]
genes_down <- deg_genes[deg_genes$logFC < 0, "Gene symbol"]
genes_down <- genes_down[!is.na(genes_down) & genes_down != ""]

print(paste("Total gen untuk enrichment:", length(genes_all)))
print(paste("Upregulated:", length(genes_up)))
print(paste("Downregulated:", length(genes_down)))

# GO Enrichment menggunakan referensi Arabidopsis thaliana
# (org.Sl.eg.db tidak tersedia di Bioconductor untuk versi ini)
go_results <- enrichGO(
  gene          = genes_all,
  OrgDb         = org.At.tair.db,
  keyType       = "SYMBOL",
  ont           = "ALL",          # BP + MF + CC sekaligus
  pAdjustMethod = "fdr",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

print(paste("GO terms signifikan:", nrow(as.data.frame(go_results))))

# Visualisasi GO - Dotplot
dotplot(go_results, showCategory = 20, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  labs(
    title    = "GO Enrichment Analysis",
    subtitle = "Dataset GSE33177 | Tomat Terinfeksi P. infestans vs Sehat"
  )

# Visualisasi GO - Barplot manual per ontologi
go_df <- as.data.frame(go_results)
go_top <- rbind(
  head(go_df[go_df$ONTOLOGY == "BP", ], 10),
  head(go_df[go_df$ONTOLOGY == "MF", ], 10),
  head(go_df[go_df$ONTOLOGY == "CC", ], 10)
)
go_top <- go_top[!is.na(go_top$Description), ]

ggplot(go_top, aes(x = Count, y = reorder(Description, Count), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
  scale_fill_manual(values = c("BP" = "steelblue", "MF" = "tomato", "CC" = "lightgreen")) +
  theme_bw() +
  labs(
    title    = "GO Enrichment - Barplot per Ontologi",
    subtitle = "Dataset GSE33177 | Tomat Terinfeksi P. infestans vs Sehat",
    x        = "Gene Count",
    y        = "",
    fill     = "Ontologi"
  )

# Visualisasi GO - Network (emapplot lebih stabil dari cnetplot)
emapplot(
  pairwise_termsim(go_results),
  showCategory = 20
) +
  labs(
    title    = "GO Enrichment Map",
    subtitle = "Dataset GSE33177"
  )


# ============================================================
# STEP 9: KEGG PATHWAY ANALYSIS
# ============================================================

# Konversi Gene Symbol ke ENTREZ ID menggunakan Arabidopsis
# (karena database tomat tidak tersedia)
genes_all_full <- allResults$"Gene symbol"
genes_all_full <- genes_all_full[!is.na(genes_all_full) & genes_all_full != ""]

gene_entrez <- bitr(
  genes_all_full,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.At.tair.db
)
print(paste("Gen berhasil dikonversi ke ENTREZ:", nrow(gene_entrez)))

# KEGG Enrichment (threshold dilonggarkan karena keterbatasan anotasi)
kegg_results <- enrichKEGG(
  gene          = gene_entrez$ENTREZID,
  organism      = "ath",      # Arabidopsis thaliana
  pAdjustMethod = "fdr",
  pvalueCutoff  = 0.25,       # dilonggarkan karena sampel kecil
  qvalueCutoff  = 0.25
)

kegg_df <- as.data.frame(kegg_results)
print(paste("Total pathway KEGG ditemukan:", nrow(kegg_df)))

# Visualisasi KEGG hanya jika ada hasil
if (nrow(kegg_df) > 0) {

  dotplot(kegg_results, showCategory = 20) +
    labs(
      title    = "KEGG Pathway Enrichment Analysis",
      subtitle = "Dataset GSE33177 | Referensi Arabidopsis thaliana"
    )

  ggplot(
    head(kegg_df[order(kegg_df$pvalue), ], 20),
    aes(x = Count, y = reorder(Description, Count), fill = pvalue)
  ) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "tomato", high = "steelblue") +
    theme_bw() +
    labs(
      title    = "KEGG Pathway - Top 20",
      subtitle = "Dataset GSE33177 | Referensi Arabidopsis thaliana",
      x        = "Gene Count",
      y        = "",
      fill     = "p-value"
    )

  emapplot(
    pairwise_termsim(kegg_results),
    showCategory = 20
  ) +
    labs(
      title    = "KEGG Pathway Map",
      subtitle = "Dataset GSE33177 | Referensi Arabidopsis thaliana"
    )

} else {
  message("⚠️  KEGG tidak menghasilkan pathway signifikan.")
  message("   Kemungkinan gen tidak terpetakan ke pathway Arabidopsis.")
  message("   Analisis dilanjutkan hanya dengan GO enrichment.")
}


# ============================================================
# STEP 10: RINGKASAN GEN UP & DOWN REGULATED
# ============================================================

# Gen Upregulated
up_genes <- allResults[
  allResults$logFC > 1 & allResults$P.Value < 0.05,
  c("PROBEID", "Gene symbol", "Gene title", "logFC", "P.Value", "adj.P.Val")
]
up_genes         <- up_genes[order(up_genes$logFC, decreasing = TRUE), ]
up_genes$FoldChange <- round(2^up_genes$logFC, 2)

# Gen Downregulated
down_genes <- allResults[
  allResults$logFC < -1 & allResults$P.Value < 0.05,
  c("PROBEID", "Gene symbol", "Gene title", "logFC", "P.Value", "adj.P.Val")
]
down_genes <- down_genes[order(down_genes$logFC), ]

print(paste("Total gen Upregulated  :", nrow(up_genes)))
print(paste("Total gen Downregulated:", nrow(down_genes)))


# ============================================================
# STEP 11: SIMPAN SEMUA HASIL
# ============================================================

# Buat folder output
dir.create("GSE33177_output", showWarnings = FALSE)

# --- Simpan tabel CSV ---
write.csv(allResults,           "GSE33177_output/GSE33177_allResults.csv",       row.names = FALSE)
write.csv(up_genes,             "GSE33177_output/GSE33177_UP_genes.csv",         row.names = FALSE)
write.csv(down_genes,           "GSE33177_output/GSE33177_DOWN_genes.csv",       row.names = FALSE)
write.csv(as.data.frame(go_results),   "GSE33177_output/GSE33177_GO_enrichment.csv",  row.names = FALSE)
write.csv(kegg_df,              "GSE33177_output/GSE33177_KEGG_enrichment.csv",  row.names = FALSE)

# --- Simpan visualisasi PNG ---
# Boxplot
png("GSE33177_output/01_boxplot.png", width = 1200, height = 800, res = 150)
boxplot(ex, col = group_colors, las = 2, outline = FALSE, names = sample_labels,
        main = "Boxplot Distribusi Nilai Ekspresi per Sampel", ylab = "Expression Value (log2)")
legend("topright", legend = c("Healthy", "LateBlight"),
       fill = c("lightgreen", "tomato"), cex = 0.8)
dev.off()

# Density Plot
png("GSE33177_output/02_density.png", width = 1200, height = 800, res = 150)
print(
  ggplot(expr_long, aes(x = Expression, color = Group)) +
    geom_density(linewidth = 1) +
    scale_color_manual(values = c("Healthy" = "lightgreen", "LateBlight" = "tomato")) +
    theme_minimal() +
    labs(title = "Distribusi Nilai Ekspresi Gen\nTomat Sehat vs Terinfeksi P. infestans")
)
dev.off()

# PCA Plot
png("GSE33177_output/03_pca.png", width = 1200, height = 800, res = 150)
print(
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text(vjust = -0.8, size = 3) +
    scale_color_manual(values = c("Healthy" = "lightgreen", "LateBlight" = "tomato")) +
    theme_minimal() +
    labs(title = "PCA Plot: Tomat Sehat vs Terinfeksi P. infestans",
         x = paste0("PC1 (", pca_var[1], "% variansi)"),
         y = paste0("PC2 (", pca_var[2], "% variansi)"))
)
dev.off()

# Volcano Plot
png("GSE33177_output/04_volcano.png", width = 1200, height = 800, res = 150)
print(
  ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("DOWN" = "steelblue", "NO" = "grey70", "UP" = "tomato")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(data = top_genes, aes(label = Gene), size = 3, color = "black") +
    theme_minimal() +
    labs(title = "Volcano Plot DEG: Tomat Terinfeksi vs Sehat")
)
dev.off()

# Heatmap
png("GSE33177_output/05_heatmap.png", width = 1200, height = 1600, res = 150)
pheatmap(
  mat_heatmap, scale = "row",
  annotation_col = annotation_col, annotation_colors = annotation_colors,
  show_colnames = TRUE, show_rownames = TRUE,
  fontsize_row = 7, clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean", clustering_method = "complete",
  color = colorRampPalette(c("steelblue", "white", "tomato"))(100),
  main = "Top 50 DEG: Tomat Terinfeksi P. infestans vs Sehat"
)
dev.off()

# Bar Chart DEG
png("GSE33177_output/06_barchart_deg.png", width = 1000, height = 800, res = 150)
print(
  ggplot(summary_long, aes(x = Comparison, y = Count, fill = Regulation)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Up" = "tomato", "Down" = "steelblue")) +
    theme_bw() +
    labs(title = "Jumlah DEG: Tomat Terinfeksi vs Sehat")
)
dev.off()

# GO Enrichment Dotplot
png("GSE33177_output/07_go_dotplot.png", width = 1200, height = 1000, res = 150)
print(
  dotplot(go_results, showCategory = 20, split = "ONTOLOGY") +
    facet_grid(ONTOLOGY ~ ., scales = "free") +
    labs(title = "GO Enrichment Analysis - GSE33177")
)
dev.off()

message("✅ Semua hasil berhasil disimpan di folder GSE33177_output/")
message("   - allResults.csv     : Semua gen hasil DEG")
message("   - UP_genes.csv       : Gen upregulated")
message("   - DOWN_genes.csv     : Gen downregulated")
message("   - GO_enrichment.csv  : Hasil GO enrichment")
message("   - KEGG_enrichment.csv: Hasil KEGG pathway")
message("   - 01 s/d 07 .png    : Semua visualisasi")
message("✅ Analisis selesai!")
