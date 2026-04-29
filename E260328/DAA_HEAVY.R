# #### DAA 2 ####
# # Updated Proteomics Analysis - DMSO vs Foretinib vs SJF8240 - HEAVY CHANNEL (ISOFORM-AWARE)
# 
# library(data.table)
# library(ggplot2)
# library(limma)
# 
# cat("🚀 LIMMA Analysis - HEAVY CHANNEL (DMSO / Foretinib / SJF8240)\n")
# cat("=============================================================\n\n")
# 
# #### 📂 PART 1: LIMMA DIFFERENTIAL EXPRESSION ANALYSIS ###
# cat("📊 PART 1: LIMMA DIFFERENTIAL EXPRESSION ANALYSIS - HEAVY CHANNEL\n")
# cat("=================================================================\n")
# 
# #### Load protein quantitations ###
# cat("📂 Loading protein quantitation data...\n")
# DT <- fread("protein_quant.csv.gz")
# 
# # Safety: ensure required columns exist (ADD Genes + Protein_description)
# req_cols <- c("Protein_group", "Sample", "Condition", "Replicate",
#               "LFQ_H", "N_precursors_H", "Genes", "Protein_description")
# missing <- setdiff(req_cols, names(DT))
# if (length(missing) > 0) stop("Missing required columns in protein_quant.csv.gz: ", paste(missing, collapse = ", "))
# 
# # Keep only your 3 conditions (in case other stuff is present)
# DT <- DT[Condition %in% c("DMSO", "Foretinib", "SJF8240")]
# 
# # Guard against empty/NA gene labels
# DT[, Genes := fifelse(is.na(Genes) | !nzchar(Genes), "UNMAPPED_GENE", Genes)]
# 
# # Create a UNIQUE feature ID per isoform/protein_group to avoid Genes collisions (e.g., MAPK14|Q16539 vs MAPK14|Q16539-2)
# DT[, Feature := paste0(Genes, "|", Protein_group)]
# 
# # Filter out NA / non-positive LFQ_H (Heavy Channel)
# DT <- DT[!is.na(LFQ_H) & LFQ_H > 0]
# cat("✅ Filtered to valid LFQ_H values (Heavy Channel)\n")
# 
# # Log2-transform LFQ_H
# DT[, log2_quant_H := log2(LFQ_H)]
# cat("✅ Log2-transformed LFQ_H values (Heavy Channel)\n")
# 
# # Filter: ≥2 precursors (Heavy Channel)
# DT <- DT[!is.na(N_precursors_H) & N_precursors_H >= 2]
# cat("✅ Filtered for ≥2 precursors (Heavy Channel)\n")
# 
# # Filter: ≥2 replicates per condition AND present in ≥2 conditions (AT FEATURE LEVEL)
# keepers <- DT[, .(n_reps = uniqueN(Replicate)), by = .(Feature, Condition)][
#   n_reps >= 2
# ][, .(n_conds = uniqueN(Condition)), by = Feature][
#   n_conds >= 3, Feature
# ]
# DT <- DT[Feature %in% keepers]
# cat("✅ Filtered features: ≥2 replicates in ≥2 conditions (Heavy Channel)\n")
# 
# # SINGLE-LEVEL Normalize across samples (median-centering) - Heavy Channel
# target_median <- DT[, .(med = median(log2_quant_H, na.rm = TRUE)), by = Sample][, median(med)]
# DT[, log2_quant_H_norm := log2_quant_H - median(log2_quant_H, na.rm = TRUE) + target_median, by = Sample]
# cat("✅ Applied sample-wise median centering (Heavy Channel)\n")
# 
# # Reshape to wide format (features x samples)
# wDT <- dcast(DT, Feature ~ Sample, value.var = "log2_quant_H_norm", fill = NA_real_)
# cat("✅ Reshaped data to wide format (Heavy Channel)\n")
# 
# # Sample map (NO Timepoint)
# sample_map <- unique(DT[, .(Sample, Condition, Replicate)])
# sample_map[, Condition := factor(Condition, levels = c("DMSO", "Foretinib", "SJF8240"))]
# 
# cat("\n🔍 Sample map:\n")
# print(sample_map[order(Condition, Replicate, Sample)])
# 
# # Design matrix: ~ 0 + Condition
# design_mat <- model.matrix(~ 0 + Condition, data = sample_map)
# colnames(design_mat) <- gsub("^Condition", "", colnames(design_mat))
# rownames(design_mat) <- sample_map$Sample
# 
# cat("\nDesign matrix:\n")
# print(design_mat)
# 
# # Align samples between design and expression matrix
# common_samples <- intersect(rownames(design_mat), colnames(wDT))
# if (length(common_samples) < 2) stop("❌ Not enough common samples between design matrix and wDT")
# 
# design_mat <- design_mat[common_samples, , drop = FALSE]
# wDT_subset <- wDT[, c("Feature", common_samples), with = FALSE]
# 
# # Convert to matrix for limma
# features <- wDT_subset$Feature
# expr_matrix <- as.matrix(wDT_subset[, -"Feature", with = FALSE])
# rownames(expr_matrix) <- features
# 
# # Remove features with all NA across samples
# expr_matrix <- expr_matrix[rowSums(!is.na(expr_matrix)) > 0, , drop = FALSE]
# 
# # Final alignment check
# if (!identical(rownames(design_mat), colnames(expr_matrix))) {
#   design_mat <- design_mat[colnames(expr_matrix), , drop = FALSE]
# }
# cat("\n✅ All checks passed, proceeding with limma fit\n")
# 
# # Fit model
# fit <- lmFit(expr_matrix, design_mat)
# cat("✅ Fitted linear models with Condition-only design (Heavy Channel)\n")
# 
# # Contrasts for your 3-condition setup
# contrast.matrix <- makeContrasts(
#   Foretinib_vs_DMSO = Foretinib - DMSO,
#   SJF8240_vs_DMSO   = SJF8240   - DMSO,
#   SJF8240_vs_Foretinib = SJF8240 - Foretinib,
#   levels = design_mat
# )
# cat("✅ Contrast matrix created\n")
# print(contrast.matrix)
# 
# fit2 <- contrasts.fit(fit, contrast.matrix)
# fit2 <- eBayes(fit2)
# cat("✅ Applied empirical Bayes moderation (Heavy Channel)\n")
# 
# #### 🧮 SE CALCULATION (robust) ###
# cat("\n🧮 CALCULATING STANDARD ERRORS (SE) - HEAVY CHANNEL\n")
# cat("==================================================\n")
# 
# contrast_names <- colnames(contrast.matrix)
# 
# results_list <- lapply(contrast_names, function(contrast) {
#   tt <- as.data.table(topTable(fit2, coef = contrast, number = Inf, sort.by = "none"))
#   tt[, Contrast := contrast]
#   tt[, Feature := rownames(expr_matrix)]
#   
#   tt[, SE_model := sqrt(fit2$s2.post) * sqrt(fit2$cov.coefficients[contrast, contrast])]
#   tt[, SE_empirical := fifelse(is.finite(t) & t != 0, abs(logFC / t), NA_real_)]
#   tt[, SE := pmax(SE_model, SE_empirical, na.rm = TRUE)]
#   tt[is.na(SE) | SE <= 0 | !is.finite(SE), SE := 0.1]
#   
#   tt[, .(Feature, logFC, SE, SE_model, SE_empirical, t, P.Value, adj.P.Val, Contrast)]
# })
# 
# results_H <- rbindlist(results_list)
# cat("✅ Calculated SEs (Heavy Channel)\n")
# 
# # Add annotations (MERGE BY Feature so isoforms never collide)
# ann <- unique(DT[, .(Feature, Genes, Protein_group, Protein_description)])
# results_H <- merge(results_H, ann, by = "Feature", all.x = TRUE)
# 
# setnames(results_H, "adj.P.Val", "adj_P_Val")
# 
# # Optional: a display logFC (median-centered within contrast) for clustering/heatmaps
# results_H[, logFC_display := logFC - median(logFC, na.rm = TRUE), by = Contrast]
# 
# # Save
# fwrite(results_H, "LIMMA_HEAVY_results_all_contrasts.csv.gz")
# cat("✅ LIMMA results saved: LIMMA_HEAVY_results_all_contrasts.csv.gz\n")
# 
# # Summary stats
# cat("\n📊 Summary stats by contrast:\n")
# print(results_H[, .(
#   N_features = .N,
#   median_logFC = median(logFC, na.rm = TRUE),
#   median_adjP = median(adj_P_Val, na.rm = TRUE),
#   median_SE = median(SE, na.rm = TRUE)
# ), by = Contrast])
# 
# #### 📈 SE distribution plot ###
# cat("\n📈 Creating SE distribution plot (Heavy Channel)...\n")
# tryCatch({
#   se_plot <- ggplot(results_H, aes(x = SE)) +
#     geom_histogram(bins = 50, alpha = 0.7) +
#     facet_wrap(~ Contrast, scales = "free_y") +
#     scale_x_log10() +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     labs(
#       title = "Standard Error (SE) Distribution by Contrast - HEAVY CHANNEL",
#       subtitle = "Contrasts: Foretinib vs DMSO, SJF8240 vs DMSO, SJF8240 vs Foretinib",
#       x = "Standard Error (SE) [log10 scale]",
#       y = "Count"
#     )
#   ggsave("SE_Distribution_HEAVY.png", se_plot, width = 12, height = 8, dpi = 300, bg = "white")
#   cat("✅ SE distribution plot saved: SE_Distribution_HEAVY.png\n")
# }, error = function(e) {
#   cat("⚠️ Could not create SE plot:", e$message, "\n")
# })
# 
# cat("\n🎯 DONE: LIMMA HEAVY CHANNEL (DMSO / Foretinib / SJF8240)\n")

#### DAA 2 ####
# Updated Proteomics Analysis - DMSO vs Foretinib vs SJF8240 - HEAVY CHANNEL (ISOFORM-AWARE)

library(data.table)
library(ggplot2)
library(limma)

cat("🚀 LIMMA Analysis - HEAVY CHANNEL (DMSO / Foretinib / SJF8240)\n")
cat("=============================================================\n\n")

#### 📂 PART 1: LIMMA DIFFERENTIAL EXPRESSION ANALYSIS ###
cat("📊 PART 1: LIMMA DIFFERENTIAL EXPRESSION ANALYSIS - HEAVY CHANNEL\n")
cat("=================================================================\n")

#### Load protein quantitations ###
cat("📂 Loading protein quantitation data...\n")
DT <- fread("protein_quant.csv.gz")

# Safety: ensure required columns exist (ADD Genes + Protein_description)
req_cols <- c("Protein_group", "Sample", "Condition", "Replicate",
              "LFQ_H", "N_precursors_H", "Genes", "Protein_description")
missing <- setdiff(req_cols, names(DT))
if (length(missing) > 0) stop("Missing required columns in protein_quant.csv.gz: ", paste(missing, collapse = ", "))

# Keep only your 3 conditions (in case other stuff is present)
DT <- DT[Condition %in% c("DMSO", "Foretinib", "SJF8240")]

# Guard against empty/NA gene labels
DT[, Genes := fifelse(is.na(Genes) | !nzchar(Genes), "UNMAPPED_GENE", Genes)]

# Create a UNIQUE feature ID per isoform/protein_group to avoid Genes collisions (e.g., MAPK14|Q16539 vs MAPK14|Q16539-2)
DT[, Feature := paste0(Genes, "|", Protein_group)]

# Filter out NA / non-positive LFQ_H (Heavy Channel)
DT <- DT[!is.na(LFQ_H) & LFQ_H > 0]
cat("✅ Filtered to valid LFQ_H values (Heavy Channel)\n")
cat("   Features after LFQ filtering:", uniqueN(DT$Feature), "\n")

# Log2-transform LFQ_H
DT[, log2_quant_H := log2(LFQ_H)]
cat("✅ Log2-transformed LFQ_H values (Heavy Channel)\n")

# Filter: ≥2 precursors (Heavy Channel)
DT <- DT[!is.na(N_precursors_H) & N_precursors_H >= 2]
cat("✅ Filtered for ≥2 precursors (Heavy Channel)\n")
cat("   Features after precursor filtering:", uniqueN(DT$Feature), "\n")

# FIXED: Filter: ≥2 replicates per condition AND present in ≥2 conditions (NOT ≥3)
keepers <- DT[, .(n_reps = uniqueN(Replicate)), by = .(Feature, Condition)][
  n_reps >= 2
][, .(n_conds = uniqueN(Condition)), by = Feature][
  n_conds >= 2, Feature  # ← FIXED: Changed from 3 to 2
]
DT <- DT[Feature %in% keepers]
cat("✅ Filtered features: ≥2 replicates in ≥2 conditions (Heavy Channel)\n")
cat("   Features after replicate/condition filtering:", uniqueN(DT$Feature), "\n")

# SINGLE-LEVEL Normalize across samples (median-centering) - Heavy Channel
target_median <- DT[, .(med = median(log2_quant_H, na.rm = TRUE)), by = Sample][, median(med)]
DT[, log2_quant_H_norm := log2_quant_H - median(log2_quant_H, na.rm = TRUE) + target_median, by = Sample]
cat("✅ Applied sample-wise median centering (Heavy Channel)\n")

# Reshape to wide format (features x samples)
wDT <- dcast(DT, Feature ~ Sample, value.var = "log2_quant_H_norm", fill = NA_real_)
cat("✅ Reshaped data to wide format (Heavy Channel)\n")
cat("   Final feature count in wide format:", nrow(wDT), "\n")

# Sample map (NO Timepoint)
sample_map <- unique(DT[, .(Sample, Condition, Replicate)])
sample_map[, Condition := factor(Condition, levels = c("DMSO", "Foretinib", "SJF8240"))]

cat("\n🔍 Sample map:\n")
print(sample_map[order(Condition, Replicate, Sample)])

# Design matrix: ~ 0 + Condition
design_mat <- model.matrix(~ 0 + Condition, data = sample_map)
colnames(design_mat) <- gsub("^Condition", "", colnames(design_mat))
rownames(design_mat) <- sample_map$Sample

cat("\nDesign matrix:\n")
print(design_mat)

# Align samples between design and expression matrix
common_samples <- intersect(rownames(design_mat), colnames(wDT))
if (length(common_samples) < 2) stop("❌ Not enough common samples between design matrix and wDT")

design_mat <- design_mat[common_samples, , drop = FALSE]
wDT_subset <- wDT[, c("Feature", common_samples), with = FALSE]

# Convert to matrix for limma
features <- wDT_subset$Feature
expr_matrix <- as.matrix(wDT_subset[, -"Feature", with = FALSE])
rownames(expr_matrix) <- features

# Remove features with all NA across samples
expr_matrix <- expr_matrix[rowSums(!is.na(expr_matrix)) > 0, , drop = FALSE]
cat("✅ Expression matrix prepared, rows:", nrow(expr_matrix), "cols:", ncol(expr_matrix), "\n")

# Final alignment check
if (!identical(rownames(design_mat), colnames(expr_matrix))) {
  design_mat <- design_mat[colnames(expr_matrix), , drop = FALSE]
}
cat("\n✅ All checks passed, proceeding with limma fit\n")

# Fit model
fit <- lmFit(expr_matrix, design_mat)
cat("✅ Fitted linear models with Condition-only design (Heavy Channel)\n")

# Contrasts for your 3-condition setup
contrast.matrix <- makeContrasts(
  Foretinib_vs_DMSO = Foretinib - DMSO,
  SJF8240_vs_DMSO   = SJF8240   - DMSO,
  SJF8240_vs_Foretinib = SJF8240 - Foretinib,
  levels = design_mat
)
cat("✅ Contrast matrix created\n")
print(contrast.matrix)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
cat("✅ Applied empirical Bayes moderation (Heavy Channel)\n")

#### 🧮 SE CALCULATION (robust) ###
cat("\n🧮 CALCULATING STANDARD ERRORS (SE) - HEAVY CHANNEL\n")
cat("==================================================\n")

contrast_names <- colnames(contrast.matrix)

results_list <- lapply(contrast_names, function(contrast) {
  tt <- as.data.table(topTable(fit2, coef = contrast, number = Inf, sort.by = "none"))
  tt[, Contrast := contrast]
  tt[, Feature := rownames(expr_matrix)]
  
  tt[, SE_model := sqrt(fit2$s2.post) * sqrt(fit2$cov.coefficients[contrast, contrast])]
  tt[, SE_empirical := fifelse(is.finite(t) & t != 0, abs(logFC / t), NA_real_)]
  tt[, SE := pmax(SE_model, SE_empirical, na.rm = TRUE)]
  tt[is.na(SE) | SE <= 0 | !is.finite(SE), SE := 0.1]
  
  tt[, .(Feature, logFC, SE, SE_model, SE_empirical, t, P.Value, adj.P.Val, Contrast)]
})

results_H <- rbindlist(results_list)
cat("✅ Calculated SEs (Heavy Channel)\n")

# Add annotations (MERGE BY Feature so isoforms never collide)
ann <- unique(DT[, .(Feature, Genes, Protein_group, Protein_description)])
results_H <- merge(results_H, ann, by = "Feature", all.x = TRUE)

setnames(results_H, "adj.P.Val", "adj_P_Val")

# Optional: a display logFC (median-centered within contrast) for clustering/heatmaps
results_H[, logFC_display := logFC - median(logFC, na.rm = TRUE), by = Contrast]

# Diagnostic checks
cat("\n🔍 Diagnostic checks after annotation merge:\n")
cat("Total rows in results_H:", nrow(results_H), "\n")
cat("Unique features:", uniqueN(results_H$Feature), "\n")
cat("Features with NA Genes:", sum(is.na(results_H$Genes)), "\n")
cat("Features with NA Protein_group:", sum(is.na(results_H$Protein_group)), "\n")
cat("Features with empty/unmapped Genes:", sum(results_H$Genes == "" | results_H$Genes == "UNMAPPED_GENE", na.rm = TRUE), "\n")

# Show distribution of significant hits per contrast
cat("\n📊 Significant features by contrast (FDR < 0.05):\n")
sig_summary <- results_H[!is.na(adj_P_Val) & adj_P_Val < 0.05, 
                         .(N_sig = .N), 
                         by = Contrast]
print(sig_summary)

cat("\n📊 Significant features with |logFC| > 0.5 (Heavy channel threshold):\n")
sig_fc_heavy <- results_H[!is.na(adj_P_Val) & adj_P_Val < 0.05 & abs(logFC) > 0.5, 
                          .(N_sig = .N), 
                          by = Contrast]
print(sig_fc_heavy)

# Check if MET exists
met_check <- results_H[Genes == "MET"]
cat("\n🔍 MET entries found:", nrow(met_check), "\n")
if (nrow(met_check) > 0) {
  print(met_check[, .(Feature, Genes, Protein_group, Contrast, logFC, P.Value, adj_P_Val)])
}

# Save
fwrite(results_H, "LIMMA_HEAVY_results_all_contrasts.csv.gz")
cat("\n✅ LIMMA results saved: LIMMA_HEAVY_results_all_contrasts.csv.gz\n")

# Summary stats
cat("\n📊 Summary stats by contrast:\n")
print(results_H[, .(
  N_features = .N,
  median_logFC = median(logFC, na.rm = TRUE),
  median_adjP = median(adj_P_Val, na.rm = TRUE),
  median_SE = median(SE, na.rm = TRUE),
  N_sig_FDR05 = sum(!is.na(adj_P_Val) & adj_P_Val < 0.05)
), by = Contrast])

#### 📈 SE distribution plot ###
cat("\n📈 Creating SE distribution plot (Heavy Channel)...\n")
tryCatch({
  se_plot <- ggplot(results_H, aes(x = SE)) +
    geom_histogram(bins = 50, alpha = 0.7, fill = "darkred") +
    facet_wrap(~ Contrast, scales = "free_y") +
    scale_x_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Standard Error (SE) Distribution by Contrast - HEAVY CHANNEL",
      subtitle = "Contrasts: Foretinib vs DMSO, SJF8240 vs DMSO, SJF8240 vs Foretinib",
      x = "Standard Error (SE) [log10 scale]",
      y = "Count"
    )
  ggsave("SE_Distribution_HEAVY.png", se_plot, width = 12, height = 8, dpi = 300, bg = "white")
  cat("✅ SE distribution plot saved: SE_Distribution_HEAVY.png\n")
}, error = function(e) {
  cat("⚠️ Could not create SE plot:", e$message, "\n")
})

cat("\n🎯 DONE: LIMMA HEAVY CHANNEL (DMSO / Foretinib / SJF8240)\n")
cat("=" %+% strrep("=", 60) %+% "\n")