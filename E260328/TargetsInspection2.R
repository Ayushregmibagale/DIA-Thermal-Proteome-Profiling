# #############################################################################################################
# #
# # CETSA-MS SILAC TARGET INSPECTION SCRIPT  v2.0
# #
# # PURPOSE
# # -------
# # Inspect a predefined kinase/target list in your CETSA-MS SILAC results.
# # Updated to consume output from CETSA_SILAC_analysis_v2.R.
# #
# # KEY CHANGES vs v1
# # -----------------
# #  - Reads Track A (SILAC ratio) files when present:
# #      long_ratio_table.tsv.gz
# #      mean_ratio_curves.tsv.gz
# #      temp_ratio_summary.tsv.gz
# #  - Protein curve plots show a two-panel layout when Track A data exist:
# #      top panel  = absolute melt curves (Track B, Drug vs DMSO fractions)
# #      bottom panel = intra-sample log2(H/L) ratio curve (Track A)
# #  - Hit table column names updated for v2:
# #      max_abs_mean_log2_ratio  ->  max_abs_mean_log2_ratio_trackB
# #      p_value_best             added (combined p-value)
# #      max_abs_mean_log2_HL     added (Track A score)
# #  - Summary scatter updated to use Track A score when available.
# #  - input_dir default updated to "CETSA_SILAC_output_v2".
# #
# # INPUT FILES EXPECTED
# # --------------------
# # From the v2 pipeline output directory:
# #   - long_fraction_table.tsv.gz              (required)
# #   - mean_curves.tsv.gz                      (required)
# #   - CETSA_SILAC_hit_table.tsv.gz            (required)
# #   - temp_effect_summary.tsv.gz              (required)
# #   - fit_predictions_mean_curves.tsv.gz      (optional)
# #   - long_ratio_table.tsv.gz                 (optional, Track A)
# #   - mean_ratio_curves.tsv.gz                (optional, Track A)
# #   - temp_ratio_summary.tsv.gz               (optional, Track A)
# #
# # OUTPUT
# # ------
# #   - matched_targets_all.tsv
# #   - matched_targets_primary.tsv
# #   - unmatched_targets.tsv
# #   - target_summary.tsv
# #   - Target_curves/  individual PNG per target
# #   - All_target_curves.pdf
# #   - Target_heatmap_mean_log2_Drug_vs_DMSO.png
# #   - Target_summary_scatter.png
# #
# #############################################################################################################
# 
# suppressPackageStartupMessages({
#   library(data.table)
#   library(ggplot2)
#   library(patchwork)
#   library(pheatmap)
#   library(matrixStats)
# })
# 
# #############################################################################################################
# # USER SETTINGS
# #############################################################################################################
# 
# # Leave as NULL to auto-detect, or set explicitly e.g. input_dir <- "CETSA_SILAC_output_v2"
# input_dir <- NULL
# 
# # ── Auto-detect input directory ───────────────────────────────────────────────
# # Searches in order: v2 output, v1 output, then any folder containing the
# # required sentinel file.  Stops with a clear message if nothing is found.
# .sentinel <- "CETSA_SILAC_hit_table.tsv.gz"
# .candidates <- c(
#   "CETSA_SILAC_output_v2",
#   "CETSA_SILAC_output"
# )
# 
# if (is.null(input_dir)) {
#   # First try the named candidates
#   .found <- Filter(function(d) file.exists(file.path(d, .sentinel)), .candidates)
#   
#   if (length(.found) == 0) {
#     # Broader scan: any subdirectory of cwd containing the sentinel
#     .subdirs <- list.dirs(".", recursive = FALSE, full.names = FALSE)
#     .found   <- Filter(function(d) file.exists(file.path(d, .sentinel)), .subdirs)
#   }
#   
#   if (length(.found) == 0) {
#     stop(
#       "Cannot find a CETSA pipeline output directory.\n",
#       "Expected to find '", .sentinel, "' inside one of: ",
#       paste(.candidates, collapse = ", "), "\n",
#       "Either run the main analysis script first, or set input_dir manually at the top of this script."
#     )
#   }
#   
#   if (length(.found) > 1) {
#     # Prefer v2 over v1 when both exist
#     .pref <- intersect(.candidates, .found)
#     .found <- if (length(.pref)) .pref else .found
#     message("Multiple CETSA output directories found: ", paste(.found, collapse = ", "),
#             "\nUsing: ", .found[1])
#   }
#   
#   input_dir <- .found[1]
# }
# 
# message("Using input directory: ", input_dir)
# 
# output_dir <- file.path(input_dir, "Target_inspection")
# dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
# 
# plot_width       <- 10
# plot_height      <- 7
# plot_dpi         <- 300
# max_fraction_cap <- 1.5
# 
# targets <- c(
#   "AXL","DDR1","MERTK","HIPK4","LOK","RET(M918T)","RET","TIE1","FLT3(K663Q)","MET(Y1235D)",
#   "FLT3","PDGFRB","EPHA3","MET(M1250T)","EPHA6","EPHB6","FLT3(D835H)","MET","FLT3(N841I)","FLT4",
#   "RET(V804L)","CSF1R","EPHA8","EPHA4","FLT3(D835Y)","TYRO3","KIT(V559D)","KIT(V559D,T670I)","EPHA7","KIT",
#   "MEK5","FLT3(ITD)","KIT(L576P)","RET(V804M)","MST1R","SLK","TIE2","KIT(A829P)","TRKC","PLK4",
#   "DDR2","EPHA2","CDK7","MKNK2","FLT1","FRK","EPHB2","ABL1(Q252H)-nonphosphorylated","MAP4K5","TRKB",
#   "BRK","PDGFRA","TRKA","BLK","YSK4","ABL1-nonphosphorylated","LCK","EPHA1","EPHA5","KIT(V559D,V654A)",
#   "LYN","KIT(D816V)","ABL1(H396P)-nonphosphorylated","FLT3(R834Q)","AURKC","EPHB1","EPHB4","MEK1","VEGFR2",
#   "MEK2","MUSK","PYK2","ROS1","HCK","LZK","ABL1(T315I)-nonphosphorylated","SRMS","ABL1(F317L)-nonphosphorylated","SRC",
#   "RIPK2","TNK1","DLK","ABL2","PIP5K2C","PFTK1","PFTAIRE2","ABL1(M351T)-phosphorylated","KIT(D816H)","ABL1(Y253F)-phosphorylated",
#   "FER","ABL1(F317I)-nonphosphorylated","FGR","AURKB","HIPK3","HPK1","ZAK","YES","INSRR","p38-delta",
#   "ABL1-phosphorylated","CHEK2","ALK","ITK","MKNK1","ERK8","BTK","EPHB3","EGFR(L861Q)","CDK11",
#   "FYN","STK33","ABL1(H396P)-phosphorylated","FES","LTK","MAP4K3","ABL1(Q252H)-phosphorylated","CSK","EGFR(G719C)",
#   "CDK8","TNNI3K","ABL1(E255K)-phosphorylated","EGFR(L747-E749del, A750P)","EGFR(L747-T751del,Sins)","SIK","ABL1(F317L)-phosphorylated","NLK","PCTK2","DMPK2",
#   "FAK","STK36","EGFR(L858R)","MAP4K4","BMX","TNK2","CDC2L5","EGFR(L747-S752del, P753S)","INSR","TNIK",
#   "HIPK2","STK35","CDKL2","ABL1(F317I)-phosphorylated","PCTK3","TESK1","SYK","MELK","ABL1(T315I)-phosphorylated","EGFR(E746-A750del)",
#   "TEC","p38-alpha","MRCKB","MLK3","EGFR(G719S)","IRAK3","p38-gamma","HIPK1","IGF1R","TTK",
#   "AURKA","EGFR","ULK3","EGFR(T790M)","NEK9","CDC2L2","RIPK5","MLK2","TXK","TAK1",
#   "CAMKK1","CLK1","ERBB2","ERBB4","CDC2L1","SIK2","MRCKA","JNK2","LIMK1","FGFR1",
#   "MYO3B","AMPK-alpha1","p38-beta","FGFR2","TBK1","CLK4","RIPK1","PCTK1","RIPK4","FGFR3",
#   "EIF2AK1","MLK1","EGFR(S752-I759del)","MST1","AMPK-alpha2","S6K1","IKK-epsilon","ANKK1","ERBB3","JAK2(JH1domain-catalytic)",
#   "CDK9","FGFR3(G697C)","LIMK2","EGFR(L858R,T790M)","MAK","ROCK2","MAP3K2","MAP3K4","DYRK2","MAP4K2",
#   "BRAF(V600E)","MYLK2","MINK","CDKL5","CSNK2A2","MYO3A","AAK1","ACVR1","ACVR1B","ACVR2A",
#   "ACVR2B","ACVRL1","ADCK3","ADCK4","AKT1","AKT2","AKT3","ARK5","ASK1","ASK2",
#   "BIKE","BMPR1A","BMPR1B","BMPR2","BRSK1","BRSK2","CAMK1G","CAMK2A","CAMK2B","CAMK2D",
#   "CAMK2G","CAMK4","CASK","CDK2","CDK3","CDK4-cyclinD1","CDK4-cyclinD3","CDK5","CDKL1","CHEK1",
#   "CLK2","CLK3","CSNK1A1","CSNK1A1L","CSNK1D","CSNK1G1","CSNK1G2","CSNK1G3","CSNK2A1","CTK",
#   "DAPK1","DAPK2","DCAMKL1","DCAMKL2","DCAMKL3","DMPK","DRAK1","DRAK2","DYRK1A","DYRK1B",
#   "ERK1","ERK2","ERK3","ERK4","ERK5","ERN1","FGFR4","GAK","GCN2(Kin.Dom.2,S808G)","GRK1",
#   "GRK4","GRK7","GSK3A","GSK3B","ICK","IKK-beta","IRAK1","IRAK4","JAK1(JH1domain-catalytic)","JAK1(JH2domain-pseudokinase)",
#   "JAK3(JH1domain-catalytic)","LATS1","LATS2","LKB1","LRRK2","LRRK2(G2019S)","MAP3K1","MAP3K15","MAPKAPK2","MAPKAPK5",
#   "MARK1","MARK2","MARK3","MARK4","MAST1","MEK3","MEK4","MEK6","MKK7","MLCK",
#   "MST2","MST3","MST4","MTOR","MYLK","MYLK4","NDR1","NDR2","NEK1","NEK11",
#   "NEK2","NEK3","NEK4","NEK5","NEK6","NEK7","NIM1","OSR1","PAK1","PAK2",
#   "PAK4","PAK6","PAK7","PDPK1","PHKG1","PHKG2","PIK3C2B","PIK3C2G","PIK3CA","PIK3CB",
#   "PIK3CD","PIK3CG","PIK4CB","PIM1","PIM2","PIM3","PIP5K1A","PIP5K1C","PIP5K2B","PKAC-alpha",
#   "PKAC-beta","PKMYT1","PKN1","PKN2","PLK1","PLK2","PLK3","PRKCD","PRKCE","PRKCH",
#   "PRKCI","PRKCQ","PRKD1","PRKD2","PRKG1","PRKG2","PRKR","PRKX","PRP4","QSK",
#   "RIOK1","RIOK2","RIOK3","ROCK1","RPS6KA4(Kin.Dom.1-N-terminal)","RPS6KA4(Kin.Dom.2-C-terminal)","RPS6KA5(Kin.Dom.1-N-terminal)","RPS6KA5(Kin.Dom.2-C-terminal)","RSK1(Kin.Dom.1-N-terminal)","RSK1(Kin.Dom.2-C-terminal)",
#   "RSK2(Kin.Dom.1-N-terminal)","RSK3(Kin.Dom.2-C-terminal)","RSK4(Kin.Dom.1-N-terminal)","RSK4(Kin.Dom.2-C-terminal)","SBK1","SgK110","SGK3","SNARK","SNRK","SRPK1",
#   "SRPK2","SRPK3","STK16","STK39","TAOK3","TGFBR1","TGFBR2","TLK1","TLK2","TRPM6",
#   "TSSK1B","TYK2(JH2domain-pseudokinase)","ULK1","ULK2","VRK2","WEE1","WEE2","YANK1","YANK2","YANK3",
#   "ZAP70","DAPK3","CIT","BRAF","TAOK1","JNK3","PAK3","CAMKK2","CDKL3","RAF1",
#   "CSNK1E","MAP3K3","TYK2(JH1domain-catalytic)","RSK3(Kin.Dom.1-N-terminal)","CAMK1D","HUNK","JNK1","YSK1","CAMK1","IKK-alpha",
#   "TAOK2","PRKD3"
# )
# 
# #############################################################################################################
# # HELPERS
# #############################################################################################################
# 
# msg <- function(...) cat(paste0(..., "\n"))
# 
# safe_mean <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else mean(x) }
# safe_sd   <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x) }
# safe_sem  <- function(x) { x <- x[is.finite(x)]; n <- length(x); if (n < 2) NA_real_ else sd(x)/sqrt(n) }
# 
# first_non_na <- function(x) {
#   x <- x[!is.na(x) & x != ""]
#   if (!length(x)) NA_character_ else x[1]
# }
# 
# clean_target_name <- function(x) {
#   x <- toupper(trimws(x))
#   gsub("\\s+", "", x)
# }
# 
# extract_base_symbol <- function(x) {
#   x <- clean_target_name(x)
#   x <- sub("\\(.*\\)", "", x)
#   x <- sub("-PHOSPHORYLATED$", "", x)
#   x <- sub("-NONPHOSPHORYLATED$", "", x)
#   x <- sub("-KIN\\.DOM\\..*$", "", x)
#   x <- sub("-CYCLIND[13]$", "", x)
#   x <- sub("-(CATALYTIC|PSEUDOKINASE)$", "", x)
#   x <- sub("-$", "", x)
#   x
# }
# 
# map_target_alias <- function(x) {
#   y <- extract_base_symbol(x)
#   alias_map <- c(
#     "VEGFR2"="KDR","VEGFR1"="FLT1","VEGFR3"="FLT4","TIE2"="TEK",
#     "TRKA"="NTRK1","TRKB"="NTRK2","TRKC"="NTRK3",
#     "PYK2"="PTK2B","FAK"="PTK2","BRK"="PTK6","FRK"="PTK5",
#     "YES"="YES1","LZK"="MAP3K13","DLK"="MAP3K12",
#     "HPK1"="MAP4K1",
#     "MEK1"="MAP2K1","MEK2"="MAP2K2","MEK3"="MAP2K3",
#     "MEK4"="MAP2K4","MEK5"="MAP2K5","MEK6"="MAP2K6",
#     "ERK1"="MAPK3","ERK2"="MAPK1","ERK3"="MAPK6",
#     "ERK4"="MAPK4","ERK5"="MAPK7","ERK8"="MAPK15",
#     "JNK1"="MAPK8","JNK2"="MAPK9","JNK3"="MAPK10",
#     "p38-alpha"="MAPK14","p38-beta"="MAPK11",
#     "p38-gamma"="MAPK12","p38-delta"="MAPK13",
#     "IKK-alpha"="CHUK","IKK-beta"="IKBKB","IKK-epsilon"="IKBKE",
#     "PKAC-alpha"="PRKACA","PKAC-beta"="PRKACB",
#     "PFTK1"="CDK14","PFTAIRE2"="CDK15",
#     "PCTK1"="CDK16","PCTK2"="CDK17","PCTK3"="CDK18",
#     "CDC2L1"="CDK11B","CDC2L2"="CDK11A","CDC2L5"="CDK20",
#     "DMPK2"="CDC42BPG","MRCKA"="CDC42BPA","MRCKB"="CDC42BPB"
#   )
#   if (y %in% names(alias_map)) alias_map[[y]] else y
# }
# 
# split_gene_tokens <- function(x) {
#   if (is.na(x) || x == "") return(character())
#   parts <- unlist(strsplit(toupper(x), ";|,|\\|"))
#   trimws(parts[nzchar(trimws(parts))])
# }
# 
# #############################################################################################################
# # READ INPUTS
# #############################################################################################################
# 
# msg("Reading CETSA v2 result files from: ", input_dir)
# 
# read_if_exists <- function(f) {
#   path <- file.path(input_dir, f)
#   if (file.exists(path)) { msg("  Reading: ", f); fread(path) }
#   else { msg("  Not found (skipping): ", f); data.table() }
# }
# 
# # Required
# for (f in c("long_fraction_table.tsv.gz", "mean_curves.tsv.gz",
#             "CETSA_SILAC_hit_table.tsv.gz", "temp_effect_summary.tsv.gz")) {
#   if (!file.exists(file.path(input_dir, f)))
#     stop("Required file missing: ", file.path(input_dir, f))
# }
# 
# long_dt      <- fread(file.path(input_dir, "long_fraction_table.tsv.gz"))
# mean_dt      <- fread(file.path(input_dir, "mean_curves.tsv.gz"))
# hit_dt       <- fread(file.path(input_dir, "CETSA_SILAC_hit_table.tsv.gz"))
# temp_eff_dt  <- fread(file.path(input_dir, "temp_effect_summary.tsv.gz"))
# 
# # Optional
# fitpred_dt        <- read_if_exists("fit_predictions_mean_curves.tsv.gz")
# ratio_long_dt     <- read_if_exists("long_ratio_table.tsv.gz")
# ratio_mean_dt     <- read_if_exists("mean_ratio_curves.tsv.gz")
# ratio_temp_dt     <- read_if_exists("temp_ratio_summary.tsv.gz")
# 
# has_track_a <- nrow(ratio_long_dt) > 0
# 
# # Detect which score column name the hit table uses (v1 vs v2 compatible)
# score_col_B <- if ("max_abs_mean_log2_ratio_trackB" %in% names(hit_dt)) {
#   "max_abs_mean_log2_ratio_trackB"
# } else if ("max_abs_mean_log2_ratio" %in% names(hit_dt)) {
#   "max_abs_mean_log2_ratio"
# } else NA_character_
# 
# score_col_A <- if ("max_abs_mean_log2_HL" %in% names(hit_dt)) "max_abs_mean_log2_HL" else NA_character_
# has_pvalue  <- "p_value_best" %in% names(hit_dt)
# 
# msg(if (has_track_a) "Track A (SILAC ratio) data found." else "Track A data not found; Track B only.")
# 
# #############################################################################################################
# # BUILD TARGET TABLE & MATCH
# #############################################################################################################
# 
# targets_dt <- data.table(target_input = unique(targets))
# targets_dt[, target_clean := clean_target_name(target_input)]
# targets_dt[, target_base  := vapply(target_input, extract_base_symbol, character(1))]
# targets_dt[, target_gene  := vapply(target_input, map_target_alias, character(1))]
# 
# # Expand hit table by gene token for multi-gene protein groups
# hit_dt[, Genes := as.character(Genes)]
# hit_dt[, Protein_group := as.character(Protein_group)]
# 
# expanded_hits <- hit_dt[, {
#   toks <- split_gene_tokens(Genes)
#   if (!length(toks)) toks <- NA_character_
#   .(gene_token = toks)
# }, by = .(Protein_group, Genes)]
# expanded_hits[, gene_token := toupper(gene_token)]
# 
# # Primary match: gene symbol
# match_gene <- merge(
#   targets_dt,
#   expanded_hits,
#   by.x = "target_gene", by.y = "gene_token",
#   all.x = TRUE, allow.cartesian = TRUE
# )
# match_gene[, match_type := fifelse(!is.na(Protein_group), "gene_match", NA_character_)]
# 
# # Fallback match: protein group string
# pg_upper <- unique(hit_dt[, .(Protein_group, pg_upper = toupper(Protein_group))])
# match_pg  <- merge(targets_dt, pg_upper,
#                    by.x = "target_clean", by.y = "pg_upper",
#                    all.x = TRUE, allow.cartesian = TRUE)
# match_pg  <- merge(match_pg, hit_dt[, .(Protein_group, Genes)],
#                    by = "Protein_group", all.x = TRUE)
# match_pg[, match_type := fifelse(!is.na(Genes), "protein_group_match", NA_character_)]
# 
# matched_all <- unique(rbindlist(list(
#   match_gene[, .(target_input, target_clean, target_base, target_gene, Protein_group, Genes, match_type)],
#   match_pg[,   .(target_input, target_clean, target_base, target_gene, Protein_group, Genes, match_type)]
# ), fill = TRUE)[!is.na(Protein_group)])
# 
# match_priority <- c(gene_match = 1L, protein_group_match = 2L)
# matched_all[, match_priority := match_priority[match_type]]
# setorder(matched_all, target_input, match_priority, Protein_group)
# matched_primary <- matched_all[, .SD[1], by = target_input]
# 
# unmatched_targets <- targets_dt[!target_input %in% matched_primary$target_input]
# 
# #############################################################################################################
# # TARGET SUMMARY TABLE
# #############################################################################################################
# 
# target_summary <- merge(
#   matched_primary,
#   hit_dt,
#   by = c("Protein_group", "Genes"),
#   all.x = TRUE
# )
# 
# # Sort: p-value first (if available), then rank score
# if (has_pvalue) {
#   target_summary[, .p_ord := fifelse(is.finite(p_value_best), p_value_best, Inf)]
#   setorder(target_summary, .p_ord, -rank_score)
#   target_summary[, .p_ord := NULL]
# } else {
#   target_summary[, .rs := fifelse(is.finite(rank_score), rank_score, -Inf)]
#   setorder(target_summary, -.rs)
#   target_summary[, .rs := NULL]
# }
# 
# #############################################################################################################
# # WRITE MATCH TABLES
# #############################################################################################################
# 
# fwrite(matched_all,       file.path(output_dir, "matched_targets_all.tsv"),     sep = "\t")
# fwrite(matched_primary,   file.path(output_dir, "matched_targets_primary.tsv"), sep = "\t")
# fwrite(unmatched_targets, file.path(output_dir, "unmatched_targets.tsv"),        sep = "\t")
# fwrite(target_summary,    file.path(output_dir, "target_summary.tsv"),          sep = "\t")
# 
# #############################################################################################################
# # PLOT FUNCTIONS
# #############################################################################################################
# 
# # ── helper: format a numeric for subtitle ─────────────────────────────────────
# fmt <- function(x, digits = 2) ifelse(is.finite(x), round(x, digits), "NA")
# 
# plot_one_target <- function(protein_group_id, target_label = NULL) {
#   
#   dt_raw   <- long_dt[Protein_group == protein_group_id]
#   dt_mean  <- mean_dt[Protein_group == protein_group_id]
#   dt_fit   <- if (nrow(fitpred_dt) > 0) fitpred_dt[Protein_group == protein_group_id] else data.table()
#   dt_ratio <- if (has_track_a) ratio_long_dt[Protein_group == protein_group_id] else data.table()
#   dt_rmean <- if (has_track_a) ratio_mean_dt[Protein_group == protein_group_id] else data.table()
#   
#   if (nrow(dt_raw) == 0) return(NULL)
#   
#   hit_row <- target_summary[Protein_group == protein_group_id][1]
#   if (is.null(target_label) || is.na(target_label))
#     target_label <- first_non_na(hit_row$target_input)
#   
#   title_gene <- first_non_na(dt_raw$Genes)
#   if (is.na(title_gene)) title_gene <- protein_group_id
#   
#   # Build subtitle with all available metrics
#   dTm   <- if ("deltaTm_mean"        %in% names(hit_row)) hit_row$deltaTm_mean[1]        else NA
#   dTm2  <- if ("deltaTm_mean_curve"  %in% names(hit_row)) hit_row$deltaTm_mean_curve[1]  else NA
#   dAUC  <- if ("deltaAUC_mean"       %in% names(hit_row)) hit_row$deltaAUC_mean[1]       else NA
#   sB    <- if (!is.na(score_col_B) && score_col_B %in% names(hit_row)) hit_row[[score_col_B]][1] else NA
#   sA    <- if (!is.na(score_col_A) && score_col_A %in% names(hit_row)) hit_row[[score_col_A]][1] else NA
#   pval  <- if (has_pvalue) hit_row$p_value_best[1] else NA
#   
#   subtxt <- paste0(
#     "Target: ", target_label,
#     " | gene: ", first_non_na(hit_row$Genes),
#     " [", first_non_na(hit_row$match_type), "]",
#     "\ndeltaTm=", fmt(dTm),
#     " | deltaTm_curve=", fmt(dTm2),
#     " | deltaAUC=", fmt(dAUC, 3),
#     if (!is.na(sB)) paste0(" | maxLog2(B)=", fmt(sB, 3)) else "",
#     if (!is.na(sA)) paste0(" | maxLog2(A)=", fmt(sA, 3)) else "",
#     if (!is.na(pval)) paste0(" | p=", signif(pval, 2)) else ""
#   )
#   
#   # ── Panel 1: absolute melt curves (Track B) ──────────────────────────────
#   p_abs <- ggplot() +
#     geom_point(
#       data = dt_raw[is.finite(frac)],
#       aes(Temp, frac, colour = condition, shape = Replicate),
#       size = 2, alpha = 0.75
#     ) +
#     geom_line(
#       data = dt_raw[is.finite(frac)],
#       aes(Temp, frac, colour = condition,
#           group = interaction(condition, Replicate)),
#       alpha = 0.4
#     ) +
#     geom_errorbar(
#       data = dt_mean[is.finite(mean_frac)],
#       aes(Temp, ymin = mean_frac - sem_frac,
#           ymax = mean_frac + sem_frac, colour = condition),
#       width = 0.25, linewidth = 0.5
#     ) +
#     geom_point(
#       data = dt_mean[is.finite(mean_frac)],
#       aes(Temp, mean_frac, colour = condition),
#       size = 3
#     ) +
#     geom_line(
#       data = dt_mean[is.finite(mean_frac)],
#       aes(Temp, mean_frac, colour = condition),
#       linewidth = 1
#     ) +
#     theme_bw(base_size = 11) +
#     ylim(0, max_fraction_cap) +
#     labs(title = title_gene, subtitle = subtxt,
#          x = "Temperature (°C)", y = "Normalized soluble fraction")
#   
#   if (nrow(dt_fit) > 0)
#     p_abs <- p_abs +
#     geom_line(
#       data = dt_fit[is.finite(frac_pred)],
#       aes(Temp, frac_pred, colour = condition),
#       linewidth = 1.1, linetype = 2
#     )
#   
#   # ── Panel 2: intra-sample SILAC ratio curve (Track A) ────────────────────
#   if (nrow(dt_ratio) > 0 && any(is.finite(dt_ratio$log2_HL_ratio))) {
#     
#     # mean ± SEM from pre-computed table if available, else compute on the fly
#     if (nrow(dt_rmean) > 0 && "mean_log2_HL" %in% names(dt_rmean)) {
#       rmean_plot <- dt_rmean[is.finite(mean_log2_HL)]
#     } else {
#       rmean_plot <- dt_ratio[is.finite(log2_HL_ratio), .(
#         mean_log2_HL = safe_mean(log2_HL_ratio),
#         sem_log2_HL  = safe_sem(log2_HL_ratio)
#       ), by = Temp][is.finite(mean_log2_HL)]
#     }
#     
#     p_ratio <- ggplot() +
#       geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
#       geom_point(
#         data = dt_ratio[is.finite(log2_HL_ratio)],
#         aes(Temp, log2_HL_ratio, shape = Replicate),
#         size = 2, alpha = 0.7, colour = "#E69F00"
#       ) +
#       geom_errorbar(
#         data = rmean_plot,
#         aes(Temp, ymin = mean_log2_HL - sem_log2_HL,
#             ymax = mean_log2_HL + sem_log2_HL),
#         width = 0.25, colour = "#E69F00", linewidth = 0.5
#       ) +
#       geom_line(
#         data = rmean_plot,
#         aes(Temp, mean_log2_HL),
#         colour = "#E69F00", linewidth = 1
#       ) +
#       geom_point(
#         data = rmean_plot,
#         aes(Temp, mean_log2_HL),
#         colour = "#E69F00", size = 3
#       ) +
#       theme_bw(base_size = 11) +
#       labs(x = "Temperature (°C)", y = "log2(Drug/DMSO)",
#            title = "Track A: intra-sample SILAC ratio")
#     
#     return(p_abs / p_ratio + plot_layout(heights = c(2, 1)))
#   }
#   
#   p_abs
# }
# 
# #############################################################################################################
# # INDIVIDUAL TARGET CURVE PLOTS
# #############################################################################################################
# 
# msg("\nMaking target curve plots...")
# 
# curve_dir <- file.path(output_dir, "Target_curves")
# dir.create(curve_dir, showWarnings = FALSE)
# 
# plot_list <- list()
# 
# if (nrow(target_summary) > 0) {
#   for (i in seq_len(nrow(target_summary))) {
#     this_pg    <- target_summary$Protein_group[i]
#     this_label <- target_summary$target_input[i]
#     
#     p <- plot_one_target(this_pg, this_label)
#     if (is.null(p)) next
#     
#     plot_list[[length(plot_list) + 1]] <- p
#     
#     out_name <- paste0(
#       sprintf("%03d", i), "_",
#       gsub("[^A-Za-z0-9_\\-]", "_", this_label), "_",
#       gsub("[^A-Za-z0-9_\\-]", "_", this_pg), ".png"
#     )
#     h_out <- if (has_track_a) plot_height * 1.4 else plot_height
#     ggsave(file.path(curve_dir, out_name), p,
#            width = plot_width, height = h_out, dpi = plot_dpi, bg = "white")
#   }
# }
# 
# #############################################################################################################
# # COMBINED PDF
# #############################################################################################################
# 
# if (length(plot_list) > 0) {
#   msg("Writing combined PDF...")
#   h_out <- if (has_track_a) plot_height * 1.4 else plot_height
#   pdf(file.path(output_dir, "All_target_curves.pdf"), width = plot_width, height = h_out)
#   for (p in plot_list) print(p)
#   dev.off()
# }
# 
# #############################################################################################################
# # TARGET HEATMAP
# # Use Track A ratio summary if available (more reliable); fall back to Track B.
# #############################################################################################################
# 
# msg("Making target heatmap...")
# 
# heat_source_dt <- if (nrow(ratio_temp_dt) > 0 && "mean_log2_ratio" %in% names(ratio_temp_dt)) {
#   msg("  Heatmap using Track A ratio data.")
#   merge(
#     matched_primary[, .(target_input, Protein_group, Genes)],
#     ratio_temp_dt[, .(Protein_group, Genes, Temp, mean_log2_ratio)],
#     by = c("Protein_group", "Genes"), all.x = TRUE
#   )
# } else {
#   msg("  Heatmap using Track B temp_effect_summary.")
#   merge(
#     matched_primary[, .(target_input, Protein_group, Genes)],
#     temp_eff_dt[, .(Protein_group, Genes, Temp, mean_log2_ratio)],
#     by = c("Protein_group", "Genes"), all.x = TRUE
#   )
# }
# 
# if (nrow(heat_source_dt) > 0) {
#   heat_wide <- dcast(
#     heat_source_dt,
#     target_input + Genes + Protein_group ~ Temp,
#     value.var = "mean_log2_ratio"
#   )
#   
#   if (nrow(heat_wide) > 1) {
#     heat_mat <- as.matrix(heat_wide[, -(1:3)])
#     rownames(heat_mat) <- paste0(heat_wide$target_input, " | ", heat_wide$Genes)
#     
#     # Order rows by max absolute shift
#     row_score <- matrixStats::rowMaxs(abs(heat_mat), na.rm = TRUE)
#     row_score[!is.finite(row_score)] <- -Inf
#     heat_mat <- heat_mat[order(row_score, decreasing = TRUE), , drop = FALSE]
#     
#     png(file.path(output_dir, "Target_heatmap_mean_log2_Drug_vs_DMSO.png"),
#         width = 1800, height = max(1200, 28 * nrow(heat_mat)), res = 200)
#     pheatmap(
#       heat_mat,
#       scale = "none", cluster_rows = TRUE, cluster_cols = FALSE,
#       border_color = NA,
#       main = "Matched targets: mean log2(Drug/DMSO) across temperatures",
#       fontsize_row = 7, fontsize_col = 10
#     )
#     dev.off()
#   }
# }
# 
# #############################################################################################################
# # SUMMARY SCATTER
# # x = deltaTm (Track B), y = Track A score if available, else Track B score
# #############################################################################################################
# 
# if (nrow(target_summary) > 0) {
#   
#   # Y axis: prefer Track A SILAC ratio score
#   y_col   <- if (!is.na(score_col_A) && score_col_A %in% names(target_summary)) score_col_A else score_col_B
#   y_label <- if (!is.na(score_col_A) && score_col_A %in% names(target_summary)) {
#     "Max |mean log2(H/L)| — Track A"
#   } else {
#     "Max |mean log2(Drug/DMSO)| — Track B"
#   }
#   
#   if (!is.null(y_col) && !is.na(y_col) && y_col %in% names(target_summary)) {
#     plot_dt <- target_summary[is.finite(deltaTm_mean) | is.finite(get(y_col))]
#     
#     # colour by p-value if available
#     if (has_pvalue && any(is.finite(target_summary$p_value_best))) {
#       plot_dt[, sig := p_value_best < 0.05]
#       p_scatter <- ggplot(plot_dt,
#                           aes(deltaTm_mean, get(y_col),
#                               colour = sig, label = target_input)) +
#         scale_colour_manual(values = c("TRUE" = "#d73027", "FALSE" = "grey50"),
#                             name = "p < 0.05", na.value = "grey80")
#     } else {
#       p_scatter <- ggplot(plot_dt,
#                           aes(deltaTm_mean, get(y_col), label = target_input))
#     }
#     
#     p_scatter <- p_scatter +
#       geom_point(size = 2, alpha = 0.85) +
#       geom_text(size = 2.5, vjust = -0.5, check_overlap = TRUE) +
#       geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
#       theme_bw(base_size = 11) +
#       labs(title = "Requested targets in CETSA-MS SILAC",
#            x = "Mean deltaTm (Drug – DMSO) [°C]",
#            y = y_label)
#     
#     ggsave(file.path(output_dir, "Target_summary_scatter.png"),
#            p_scatter, width = 11, height = 8, dpi = plot_dpi, bg = "white")
#   }
# }
# 
# #############################################################################################################
# # CONSOLE SUMMARY
# #############################################################################################################
# 
# msg("\n══════════════════════════════════════════════════════")
# msg("Done.  Output folder: ", output_dir)
# msg("══════════════════════════════════════════════════════")
# msg("Requested targets  : ", nrow(targets_dt))
# msg("Matched            : ", nrow(matched_primary))
# msg("Unmatched          : ", nrow(unmatched_targets))
# 
# if (nrow(unmatched_targets) > 0) {
#   msg("Unmatched targets  : ",
#       paste(head(unmatched_targets$target_input, 20), collapse = ", "),
#       if (nrow(unmatched_targets) > 20) " ..." else "")
# }
# 
# if (nrow(target_summary) > 0) {
#   show_cols <- intersect(
#     c("target_input","Genes","Protein_group","match_type",
#       "deltaTm_mean","deltaTm_mean_curve","deltaAUC_mean",
#       score_col_A, score_col_B, "p_value_best"),
#     names(target_summary)
#   )
#   msg("\nTop matched targets:")
#   print(target_summary[seq_len(min(30, .N)), ..show_cols])
# }
# 


#############################################################################################################
#
# CETSA-MS SILAC TARGET INSPECTION SCRIPT  v2.1
#
# PURPOSE
# -------
# Inspect a predefined kinase/target list in CETSA-MS SILAC results.
# This version is aligned with CETSA-MS SILAC analysis v2.1, where:
#   - Track A SILAC H/L ratio is primary evidence.
#   - Track B LFQ melt-curve deltaTm is secondary/supportive evidence.
#   - The main analysis output directory defaults to CETSA_SILAC_output_v2_1.
#
# INPUT FILES EXPECTED
# --------------------
# Required:
#   - long_fraction_table.tsv.gz
#   - mean_curves.tsv.gz
#   - CETSA_SILAC_hit_table.tsv.gz
#   - temp_effect_summary.tsv.gz
#
# Optional:
#   - fit_predictions_mean_curves.tsv.gz
#   - long_ratio_table.tsv.gz
#   - mean_ratio_curves.tsv.gz
#   - temp_ratio_summary.tsv.gz
#
# MAIN OUTPUT
# -----------
#   - matched_targets_all.tsv
#   - matched_targets_primary.tsv
#   - unmatched_targets.tsv
#   - target_summary.tsv
#   - Target_curves/ individual PNGs
#   - All_target_curves.pdf
#   - Target_heatmap_mean_log2_Drug_vs_DMSO.png
#   - Target_summary_scatter.png
#
# IMPORTANT INTERPRETATION
# ------------------------
# For SILAC CETSA-MS, prioritize:
#   1. Track A: signed_peak_log2_HL / max_abs_mean_log2_HL
#   2. Track B: deltaTm_mean
#   3. Track B: deltaTm_mean_curve
#
#############################################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(matrixStats)
})

#############################################################################################################
# USER SETTINGS
#############################################################################################################

# Leave as NULL to auto-detect, or set explicitly, e.g.:
# input_dir <- "CETSA_SILAC_output_v2_1"
input_dir <- NULL

.sentinel <- "CETSA_SILAC_hit_table.tsv.gz"
.candidates <- c(
  "CETSA_SILAC_output_v2_1",
  "CETSA_SILAC_output_v2",
  "CETSA_SILAC_output"
)

if (is.null(input_dir)) {
  .found <- Filter(function(d) file.exists(file.path(d, .sentinel)), .candidates)
  
  if (length(.found) == 0) {
    .subdirs <- list.dirs(".", recursive = FALSE, full.names = FALSE)
    .found <- Filter(function(d) file.exists(file.path(d, .sentinel)), .subdirs)
  }
  
  if (length(.found) == 0) {
    stop(
      "Cannot find a CETSA pipeline output directory.\n",
      "Expected to find '", .sentinel, "' inside one of: ",
      paste(.candidates, collapse = ", "), "\n",
      "Run the main analysis script first, or set input_dir manually."
    )
  }
  
  if (length(.found) > 1) {
    .pref <- intersect(.candidates, .found)
    .found <- if (length(.pref)) .pref else .found
    message("Multiple CETSA output directories found: ", paste(.found, collapse = ", "),
            "\nUsing: ", .found[1])
  }
  
  input_dir <- .found[1]
}

message("Using input directory: ", input_dir)

output_dir <- file.path(input_dir, "Target_inspection")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

plot_width       <- 10
plot_height      <- 7
plot_dpi         <- 300
max_fraction_cap <- 1.5

# Replace or edit this vector as needed.
targets <- c(
  "AXL","DDR1","MERTK","HIPK4","LOK","RET(M918T)","RET","TIE1","FLT3(K663Q)","MET(Y1235D)",
  "FLT3","PDGFRB","EPHA3","MET(M1250T)","EPHA6","EPHB6","FLT3(D835H)","MET","FLT3(N841I)","FLT4",
  "RET(V804L)","CSF1R","EPHA8","EPHA4","FLT3(D835Y)","TYRO3","KIT(V559D)","KIT(V559D,T670I)","EPHA7","KIT",
  "MEK5","FLT3(ITD)","KIT(L576P)","RET(V804M)","MST1R","SLK","TIE2","KIT(A829P)","TRKC","PLK4",
  "DDR2","EPHA2","CDK7","MKNK2","FLT1","FRK","EPHB2","ABL1(Q252H)-nonphosphorylated","MAP4K5","TRKB",
  "BRK","PDGFRA","TRKA","BLK","YSK4","ABL1-nonphosphorylated","LCK","EPHA1","EPHA5","KIT(V559D,V654A)",
  "LYN","KIT(D816V)","ABL1(H396P)-nonphosphorylated","FLT3(R834Q)","AURKC","EPHB1","EPHB4","MEK1","VEGFR2",
  "MEK2","MUSK","PYK2","ROS1","HCK","LZK","ABL1(T315I)-nonphosphorylated","SRMS","ABL1(F317L)-nonphosphorylated","SRC",
  "RIPK2","TNK1","DLK","ABL2","PIP5K2C","PFTK1","PFTAIRE2","ABL1(M351T)-phosphorylated","KIT(D816H)","ABL1(Y253F)-phosphorylated",
  "FER","ABL1(F317I)-nonphosphorylated","FGR","AURKB","HIPK3","HPK1","ZAK","YES","INSRR","p38-delta",
  "ABL1-phosphorylated","CHEK2","ALK","ITK","MKNK1","ERK8","BTK","EPHB3","EGFR(L861Q)","CDK11",
  "FYN","STK33","ABL1(H396P)-phosphorylated","FES","LTK","MAP4K3","ABL1(Q252H)-phosphorylated","CSK","EGFR(G719C)",
  "CDK8","TNNI3K","ABL1(E255K)-phosphorylated","EGFR(L747-E749del, A750P)","EGFR(L747-T751del,Sins)","SIK","ABL1(F317L)-phosphorylated","NLK","PCTK2","DMPK2",
  "FAK","STK36","EGFR(L858R)","MAP4K4","BMX","TNK2","CDC2L5","EGFR(L747-S752del, P753S)","INSR","TNIK",
  "HIPK2","STK35","CDKL2","ABL1(F317I)-phosphorylated","PCTK3","TESK1","SYK","MELK","ABL1(T315I)-phosphorylated","EGFR(E746-A750del)",
  "TEC","p38-alpha","MRCKB","MLK3","EGFR(G719S)","IRAK3","p38-gamma","HIPK1","IGF1R","TTK",
  "AURKA","EGFR","ULK3","EGFR(T790M)","NEK9","CDC2L2","RIPK5","MLK2","TXK","TAK1",
  "CAMKK1","CLK1","ERBB2","ERBB4","CDC2L1","SIK2","MRCKA","JNK2","LIMK1","FGFR1",
  "MYO3B","AMPK-alpha1","p38-beta","FGFR2","TBK1","CLK4","RIPK1","PCTK1","RIPK4","FGFR3",
  "EIF2AK1","MLK1","EGFR(S752-I759del)","MST1","AMPK-alpha2","S6K1","IKK-epsilon","ANKK1","ERBB3","JAK2(JH1domain-catalytic)",
  "CDK9","FGFR3(G697C)","LIMK2","EGFR(L858R,T790M)","MAK","ROCK2","MAP3K2","MAP3K4","DYRK2","MAP4K2",
  "BRAF(V600E)","MYLK2","MINK","CDKL5","CSNK2A2","MYO3A","AAK1","ACVR1","ACVR1B","ACVR2A",
  "ACVR2B","ACVRL1","ADCK3","ADCK4","AKT1","AKT2","AKT3","ARK5","ASK1","ASK2",
  "BIKE","BMPR1A","BMPR1B","BMPR2","BRSK1","BRSK2","CAMK1G","CAMK2A","CAMK2B","CAMK2D",
  "CAMK2G","CAMK4","CASK","CDK2","CDK3","CDK4-cyclinD1","CDK4-cyclinD3","CDK5","CDKL1","CHEK1",
  "CLK2","CLK3","CSNK1A1","CSNK1A1L","CSNK1D","CSNK1G1","CSNK1G2","CSNK1G3","CSNK2A1","CTK",
  "DAPK1","DAPK2","DCAMKL1","DCAMKL2","DCAMKL3","DMPK","DRAK1","DRAK2","DYRK1A","DYRK1B",
  "ERK1","ERK2","ERK3","ERK4","ERK5","ERN1","FGFR4","GAK","GCN2(Kin.Dom.2,S808G)","GRK1",
  "GRK4","GRK7","GSK3A","GSK3B","ICK","IKK-beta","IRAK1","IRAK4","JAK1(JH1domain-catalytic)","JAK1(JH2domain-pseudokinase)",
  "JAK3(JH1domain-catalytic)","LATS1","LATS2","LKB1","LRRK2","LRRK2(G2019S)","MAP3K1","MAP3K15","MAPKAPK2","MAPKAPK5",
  "MARK1","MARK2","MARK3","MARK4","MAST1","MEK3","MEK4","MEK6","MKK7","MLCK",
  "MST2","MST3","MST4","MTOR","MYLK","MYLK4","NDR1","NDR2","NEK1","NEK11",
  "NEK2","NEK3","NEK4","NEK5","NEK6","NEK7","NIM1","OSR1","PAK1","PAK2",
  "PAK4","PAK6","PAK7","PDPK1","PHKG1","PHKG2","PIK3C2B","PIK3C2G","PIK3CA","PIK3CB",
  "PIK3CD","PIK3CG","PIK4CB","PIM1","PIM2","PIM3","PIP5K1A","PIP5K1C","PIP5K2B","PKAC-alpha",
  "PKAC-beta","PKMYT1","PKN1","PKN2","PLK1","PLK2","PLK3","PRKCD","PRKCE","PRKCH",
  "PRKCI","PRKCQ","PRKD1","PRKD2","PRKG1","PRKG2","PRKR","PRKX","PRP4","QSK",
  "RIOK1","RIOK2","RIOK3","ROCK1","RPS6KA4(Kin.Dom.1-N-terminal)","RPS6KA4(Kin.Dom.2-C-terminal)","RPS6KA5(Kin.Dom.1-N-terminal)","RPS6KA5(Kin.Dom.2-C-terminal)","RSK1(Kin.Dom.1-N-terminal)","RSK1(Kin.Dom.2-C-terminal)",
  "RSK2(Kin.Dom.1-N-terminal)","RSK3(Kin.Dom.2-C-terminal)","RSK4(Kin.Dom.1-N-terminal)","RSK4(Kin.Dom.2-C-terminal)","SBK1","SgK110","SGK3","SNARK","SNRK","SRPK1",
  "SRPK2","SRPK3","STK16","STK39","TAOK3","TGFBR1","TGFBR2","TLK1","TLK2","TRPM6",
  "TSSK1B","TYK2(JH2domain-pseudokinase)","ULK1","ULK2","VRK2","WEE1","WEE2","YANK1","YANK2","YANK3",
  "ZAP70","DAPK3","CIT","BRAF","TAOK1","JNK3","PAK3","CAMKK2","CDKL3","RAF1",
  "CSNK1E","MAP3K3","TYK2(JH1domain-catalytic)","RSK3(Kin.Dom.1-N-terminal)","CAMK1D","HUNK","JNK1","YSK1","CAMK1","IKK-alpha",
  "TAOK2","PRKD3"
)

#############################################################################################################
# HELPERS
#############################################################################################################

msg <- function(...) cat(paste0(..., "\n"))

safe_mean <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else mean(x) }
safe_sd   <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x) }
safe_sem  <- function(x) { x <- x[is.finite(x)]; n <- length(x); if (n < 2) NA_real_ else sd(x) / sqrt(n) }

first_non_na <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (!length(x)) NA_character_ else as.character(x[1])
}

clean_target_name <- function(x) {
  x <- toupper(trimws(x))
  gsub("\\s+", "", x)
}

extract_base_symbol <- function(x) {
  x <- clean_target_name(x)
  x <- sub("\\(.*\\)", "", x)
  x <- sub("-PHOSPHORYLATED$", "", x)
  x <- sub("-NONPHOSPHORYLATED$", "", x)
  x <- sub("-KIN\\.DOM\\..*$", "", x)
  x <- sub("-CYCLIND[13]$", "", x)
  x <- sub("-(CATALYTIC|PSEUDOKINASE)$", "", x)
  x <- sub("-$", "", x)
  x
}

map_target_alias <- function(x) {
  y <- extract_base_symbol(x)
  alias_map <- c(
    "VEGFR2" = "KDR", "VEGFR1" = "FLT1", "VEGFR3" = "FLT4", "TIE2" = "TEK",
    "TRKA" = "NTRK1", "TRKB" = "NTRK2", "TRKC" = "NTRK3",
    "PYK2" = "PTK2B", "FAK" = "PTK2", "BRK" = "PTK6", "FRK" = "PTK5",
    "YES" = "YES1", "LZK" = "MAP3K13", "DLK" = "MAP3K12",
    "HPK1" = "MAP4K1",
    "MEK1" = "MAP2K1", "MEK2" = "MAP2K2", "MEK3" = "MAP2K3",
    "MEK4" = "MAP2K4", "MEK5" = "MAP2K5", "MEK6" = "MAP2K6",
    "ERK1" = "MAPK3", "ERK2" = "MAPK1", "ERK3" = "MAPK6",
    "ERK4" = "MAPK4", "ERK5" = "MAPK7", "ERK8" = "MAPK15",
    "JNK1" = "MAPK8", "JNK2" = "MAPK9", "JNK3" = "MAPK10",
    "P38-ALPHA" = "MAPK14", "P38-BETA" = "MAPK11",
    "P38-GAMMA" = "MAPK12", "P38-DELTA" = "MAPK13",
    "IKK-ALPHA" = "CHUK", "IKK-BETA" = "IKBKB", "IKK-EPSILON" = "IKBKE",
    "PKAC-ALPHA" = "PRKACA", "PKAC-BETA" = "PRKACB",
    "PFTK1" = "CDK14", "PFTAIRE2" = "CDK15",
    "PCTK1" = "CDK16", "PCTK2" = "CDK17", "PCTK3" = "CDK18",
    "CDC2L1" = "CDK11B", "CDC2L2" = "CDK11A", "CDC2L5" = "CDK20",
    "DMPK2" = "CDC42BPG", "MRCKA" = "CDC42BPA", "MRCKB" = "CDC42BPB"
  )
  if (y %in% names(alias_map)) alias_map[[y]] else y
}

split_gene_tokens <- function(x) {
  if (is.na(x) || x == "") return(character())
  parts <- unlist(strsplit(toupper(x), ";|,|\\|"))
  trimws(parts[nzchar(trimws(parts))])
}

fmt <- function(x, digits = 2) {
  ifelse(is.finite(x), as.character(round(x, digits)), "NA")
}

#############################################################################################################
# READ INPUTS
#############################################################################################################

msg("Reading CETSA result files from: ", input_dir)

read_if_exists <- function(f) {
  path <- file.path(input_dir, f)
  if (file.exists(path)) {
    msg("  Reading: ", f)
    fread(path)
  } else {
    msg("  Not found, skipping: ", f)
    data.table()
  }
}

required_files <- c(
  "long_fraction_table.tsv.gz",
  "mean_curves.tsv.gz",
  "CETSA_SILAC_hit_table.tsv.gz",
  "temp_effect_summary.tsv.gz"
)

for (f in required_files) {
  if (!file.exists(file.path(input_dir, f))) {
    stop("Required file missing: ", file.path(input_dir, f))
  }
}

long_dt     <- fread(file.path(input_dir, "long_fraction_table.tsv.gz"))
mean_dt     <- fread(file.path(input_dir, "mean_curves.tsv.gz"))
hit_dt      <- fread(file.path(input_dir, "CETSA_SILAC_hit_table.tsv.gz"))
temp_eff_dt <- fread(file.path(input_dir, "temp_effect_summary.tsv.gz"))

fitpred_dt    <- read_if_exists("fit_predictions_mean_curves.tsv.gz")
ratio_long_dt <- read_if_exists("long_ratio_table.tsv.gz")
ratio_mean_dt <- read_if_exists("mean_ratio_curves.tsv.gz")
ratio_temp_dt <- read_if_exists("temp_ratio_summary.tsv.gz")

has_track_a <- nrow(ratio_long_dt) > 0

# v2.1-compatible columns
score_col_A <- if ("max_abs_mean_log2_HL" %in% names(hit_dt)) "max_abs_mean_log2_HL" else NA_character_
signed_col_A <- if ("signed_peak_log2_HL" %in% names(hit_dt)) "signed_peak_log2_HL" else NA_character_
score_col_B <- if ("max_abs_mean_log2_ratio_trackB" %in% names(hit_dt)) {
  "max_abs_mean_log2_ratio_trackB"
} else if ("max_abs_mean_log2_ratio" %in% names(hit_dt)) {
  "max_abs_mean_log2_ratio"
} else NA_character_

has_pvalue <- "p_value_best" %in% names(hit_dt)

if (!("rank_score" %in% names(hit_dt))) hit_dt[, rank_score := NA_real_]
if (!("primary_evidence" %in% names(hit_dt)) && "evidence_track" %in% names(hit_dt)) setnames(hit_dt, "evidence_track", "primary_evidence")
if (!("primary_evidence" %in% names(hit_dt))) hit_dt[, primary_evidence := NA_character_]
if (!("stabilization_direction" %in% names(hit_dt))) hit_dt[, stabilization_direction := NA_character_]

msg(if (has_track_a) "Track A SILAC ratio data found." else "Track A data not found; Track B only.")

#############################################################################################################
# BUILD TARGET TABLE AND MATCH
#############################################################################################################

targets_dt <- data.table(target_input = unique(targets))
targets_dt[, target_clean := clean_target_name(target_input)]
targets_dt[, target_base  := vapply(target_input, extract_base_symbol, character(1))]
targets_dt[, target_gene  := vapply(target_input, map_target_alias, character(1))]

hit_dt[, Genes := as.character(Genes)]
hit_dt[, Protein_group := as.character(Protein_group)]

expanded_hits <- hit_dt[, {
  toks <- split_gene_tokens(Genes)
  if (!length(toks)) toks <- NA_character_
  .(gene_token = toks)
}, by = .(Protein_group, Genes)]
expanded_hits[, gene_token := toupper(gene_token)]

# Primary match: target gene symbol/alias to Genes field.
match_gene <- merge(
  targets_dt,
  expanded_hits,
  by.x = "target_gene", by.y = "gene_token",
  all.x = TRUE, allow.cartesian = TRUE
)
match_gene[, match_type := fifelse(!is.na(Protein_group), "gene_match", NA_character_)]

# Fallback match: exact cleaned target against protein group string.
pg_upper <- unique(hit_dt[, .(Protein_group, pg_upper = toupper(Protein_group))])
match_pg <- merge(
  targets_dt, pg_upper,
  by.x = "target_clean", by.y = "pg_upper",
  all.x = TRUE, allow.cartesian = TRUE
)
match_pg <- merge(match_pg, hit_dt[, .(Protein_group, Genes)], by = "Protein_group", all.x = TRUE)
match_pg[, match_type := fifelse(!is.na(Genes), "protein_group_match", NA_character_)]

matched_all <- unique(rbindlist(list(
  match_gene[, .(target_input, target_clean, target_base, target_gene, Protein_group, Genes, match_type)],
  match_pg[, .(target_input, target_clean, target_base, target_gene, Protein_group, Genes, match_type)]
), fill = TRUE)[!is.na(Protein_group)])

match_priority <- c(gene_match = 1L, protein_group_match = 2L)
matched_all[, match_priority := match_priority[match_type]]
setorder(matched_all, target_input, match_priority, Protein_group)
matched_primary <- matched_all[, .SD[1], by = target_input]

unmatched_targets <- targets_dt[!target_input %in% matched_primary$target_input]

#############################################################################################################
# TARGET SUMMARY TABLE
#############################################################################################################

target_summary <- merge(
  matched_primary,
  hit_dt,
  by = c("Protein_group", "Genes"),
  all.x = TRUE
)

# Track A-first sorting, matching the v2.1 main analysis logic.
if (!is.na(score_col_A) && score_col_A %in% names(target_summary)) {
  target_summary[, .has_trackA := is.finite(get(score_col_A))]
  target_summary[, .trackA_score := fifelse(is.finite(get(score_col_A)), get(score_col_A), -Inf)]
} else {
  target_summary[, .has_trackA := FALSE]
  target_summary[, .trackA_score := -Inf]
}

target_summary[, .trackA_p := if ("p_value_ratio_track" %in% names(target_summary)) {
  fifelse(is.finite(p_value_ratio_track), p_value_ratio_track, Inf)
} else Inf]

target_summary[, .deltaTm_rep := if ("deltaTm_mean" %in% names(target_summary)) {
  fifelse(is.finite(deltaTm_mean), abs(deltaTm_mean), -Inf)
} else -Inf]

target_summary[, .deltaTm_mean := if ("deltaTm_mean_curve" %in% names(target_summary)) {
  fifelse(is.finite(deltaTm_mean_curve), abs(deltaTm_mean_curve), -Inf)
} else -Inf]

setorder(target_summary, -.has_trackA, -.trackA_score, .trackA_p, -.deltaTm_rep, -.deltaTm_mean)
target_summary[, c(".has_trackA", ".trackA_score", ".trackA_p", ".deltaTm_rep", ".deltaTm_mean") := NULL]

#############################################################################################################
# WRITE MATCH TABLES
#############################################################################################################

fwrite(matched_all,       file.path(output_dir, "matched_targets_all.tsv"), sep = "\t")
fwrite(matched_primary,   file.path(output_dir, "matched_targets_primary.tsv"), sep = "\t")
fwrite(unmatched_targets, file.path(output_dir, "unmatched_targets.tsv"), sep = "\t")
fwrite(target_summary,    file.path(output_dir, "target_summary.tsv"), sep = "\t")

#############################################################################################################
# PLOT FUNCTIONS
#############################################################################################################

plot_one_target <- function(protein_group_id, target_label = NULL) {
  dt_raw   <- long_dt[Protein_group == protein_group_id]
  dt_mean  <- mean_dt[Protein_group == protein_group_id]
  dt_fit   <- if (nrow(fitpred_dt) > 0) fitpred_dt[Protein_group == protein_group_id] else data.table()
  dt_ratio <- if (has_track_a) ratio_long_dt[Protein_group == protein_group_id] else data.table()
  dt_rmean <- if (has_track_a) ratio_mean_dt[Protein_group == protein_group_id] else data.table()
  
  if (nrow(dt_raw) == 0) return(NULL)
  
  hit_row <- target_summary[Protein_group == protein_group_id][1]
  if (is.null(target_label) || is.na(target_label)) target_label <- first_non_na(hit_row$target_input)
  
  title_gene <- first_non_na(dt_raw$Genes)
  if (is.na(title_gene)) title_gene <- protein_group_id
  
  dTm   <- if ("deltaTm_mean" %in% names(hit_row)) hit_row$deltaTm_mean[1] else NA_real_
  dTm2  <- if ("deltaTm_mean_curve" %in% names(hit_row)) hit_row$deltaTm_mean_curve[1] else NA_real_
  dAUC  <- if ("deltaAUC_mean" %in% names(hit_row)) hit_row$deltaAUC_mean[1] else NA_real_
  sB    <- if (!is.na(score_col_B) && score_col_B %in% names(hit_row)) hit_row[[score_col_B]][1] else NA_real_
  sA    <- if (!is.na(score_col_A) && score_col_A %in% names(hit_row)) hit_row[[score_col_A]][1] else NA_real_
  signedA <- if (!is.na(signed_col_A) && signed_col_A %in% names(hit_row)) hit_row[[signed_col_A]][1] else NA_real_
  pval  <- if (has_pvalue) hit_row$p_value_best[1] else NA_real_
  psource <- if ("p_value_source" %in% names(hit_row)) first_non_na(hit_row$p_value_source) else NA_character_
  evidence <- if ("primary_evidence" %in% names(hit_row)) first_non_na(hit_row$primary_evidence) else NA_character_
  direction <- if ("stabilization_direction" %in% names(hit_row)) first_non_na(hit_row$stabilization_direction) else NA_character_
  
  subtxt <- paste0(
    "Target: ", target_label,
    " | gene: ", first_non_na(hit_row$Genes),
    " [", first_non_na(hit_row$match_type), "]",
    "\nPrimary evidence: ", ifelse(is.na(evidence), "NA", evidence),
    " | direction: ", ifelse(is.na(direction), "NA", direction),
    "\nTrack A signed peak log2HL=", fmt(signedA, 3),
    " | max|log2HL|=", fmt(sA, 3),
    " | p=", ifelse(is.finite(pval), signif(pval, 2), "NA"),
    if (!is.na(psource)) paste0(" [", psource, "]") else "",
    "\nTrack B deltaTm=", fmt(dTm),
    " | deltaTm_curve=", fmt(dTm2),
    " | deltaAUC=", fmt(dAUC, 3),
    if (!is.na(sB)) paste0(" | maxLog2(B)=", fmt(sB, 3)) else ""
  )
  
  p_abs <- ggplot() +
    geom_point(
      data = dt_raw[is.finite(frac)],
      aes(Temp, frac, colour = condition, shape = Replicate),
      size = 2, alpha = 0.75
    ) +
    geom_line(
      data = dt_raw[is.finite(frac)],
      aes(Temp, frac, colour = condition, group = interaction(condition, Replicate)),
      alpha = 0.4
    ) +
    geom_errorbar(
      data = dt_mean[is.finite(mean_frac)],
      aes(Temp, ymin = mean_frac - sem_frac, ymax = mean_frac + sem_frac, colour = condition),
      width = 0.25, linewidth = 0.5
    ) +
    geom_point(
      data = dt_mean[is.finite(mean_frac)],
      aes(Temp, mean_frac, colour = condition),
      size = 3
    ) +
    geom_line(
      data = dt_mean[is.finite(mean_frac)],
      aes(Temp, mean_frac, colour = condition),
      linewidth = 1
    ) +
    theme_bw(base_size = 11) +
    ylim(0, max_fraction_cap) +
    labs(
      title = title_gene,
      subtitle = subtxt,
      x = "Temperature (°C)",
      y = "Fraction vs lowest temperature, Track B"
    )
  
  if (nrow(dt_fit) > 0) {
    p_abs <- p_abs +
      geom_line(
        data = dt_fit[is.finite(frac_pred)],
        aes(Temp, frac_pred, colour = condition),
        linewidth = 1.1, linetype = 2
      )
  }
  
  if (nrow(dt_ratio) > 0 && any(is.finite(dt_ratio$log2_HL_ratio))) {
    if (nrow(dt_rmean) > 0 && "mean_log2_HL" %in% names(dt_rmean)) {
      rmean_plot <- dt_rmean[is.finite(mean_log2_HL)]
    } else {
      rmean_plot <- dt_ratio[is.finite(log2_HL_ratio), .(
        mean_log2_HL = safe_mean(log2_HL_ratio),
        sem_log2_HL  = safe_sem(log2_HL_ratio)
      ), by = Temp][is.finite(mean_log2_HL)]
    }
    
    p_ratio <- ggplot() +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_point(
        data = dt_ratio[is.finite(log2_HL_ratio)],
        aes(Temp, log2_HL_ratio, shape = Replicate),
        size = 2, alpha = 0.7, colour = "#E69F00"
      ) +
      geom_errorbar(
        data = rmean_plot,
        aes(Temp, ymin = mean_log2_HL - sem_log2_HL, ymax = mean_log2_HL + sem_log2_HL),
        width = 0.25, colour = "#E69F00", linewidth = 0.5
      ) +
      geom_line(
        data = rmean_plot,
        aes(Temp, mean_log2_HL),
        colour = "#E69F00", linewidth = 1
      ) +
      geom_point(
        data = rmean_plot,
        aes(Temp, mean_log2_HL),
        colour = "#E69F00", size = 3
      ) +
      theme_bw(base_size = 11) +
      labs(
        x = "Temperature (°C)",
        y = "normalized log2(Drug/DMSO), Track A",
        title = "Track A: intra-sample SILAC ratio"
      )
    
    return(p_abs / p_ratio + plot_layout(heights = c(2, 1)))
  }
  
  p_abs
}

#############################################################################################################
# INDIVIDUAL TARGET CURVE PLOTS
#############################################################################################################

msg("\nMaking target curve plots...")

curve_dir <- file.path(output_dir, "Target_curves")
dir.create(curve_dir, showWarnings = FALSE)

plot_list <- list()

if (nrow(target_summary) > 0) {
  for (i in seq_len(nrow(target_summary))) {
    this_pg <- target_summary$Protein_group[i]
    this_label <- target_summary$target_input[i]
    
    p <- plot_one_target(this_pg, this_label)
    if (is.null(p)) next
    
    plot_list[[length(plot_list) + 1]] <- p
    
    out_name <- paste0(
      sprintf("%03d", i), "_",
      gsub("[^A-Za-z0-9_\\-]", "_", this_label), "_",
      gsub("[^A-Za-z0-9_\\-]", "_", this_pg), ".png"
    )
    h_out <- if (has_track_a) plot_height * 1.4 else plot_height
    ggsave(file.path(curve_dir, out_name), p,
           width = plot_width, height = h_out, dpi = plot_dpi, bg = "white")
  }
}

#############################################################################################################
# COMBINED PDF
#############################################################################################################

if (length(plot_list) > 0) {
  msg("Writing combined PDF...")
  h_out <- if (has_track_a) plot_height * 1.4 else plot_height
  pdf(file.path(output_dir, "All_target_curves.pdf"), width = plot_width, height = h_out)
  for (p in plot_list) print(p)
  dev.off()
}

#############################################################################################################
# TARGET HEATMAP
#############################################################################################################

msg("Making target heatmap...")

heat_source_dt <- if (nrow(ratio_temp_dt) > 0 && "mean_log2_ratio" %in% names(ratio_temp_dt)) {
  msg("  Heatmap using Track A ratio data.")
  merge(
    matched_primary[, .(target_input, Protein_group, Genes)],
    ratio_temp_dt[, .(Protein_group, Genes, Temp, mean_log2_ratio)],
    by = c("Protein_group", "Genes"), all.x = TRUE
  )
} else {
  msg("  Heatmap using Track B temp_effect_summary.")
  merge(
    matched_primary[, .(target_input, Protein_group, Genes)],
    temp_eff_dt[, .(Protein_group, Genes, Temp, mean_log2_ratio)],
    by = c("Protein_group", "Genes"), all.x = TRUE
  )
}

if (nrow(heat_source_dt) > 0) {
  heat_wide <- dcast(
    heat_source_dt,
    target_input + Genes + Protein_group ~ Temp,
    value.var = "mean_log2_ratio"
  )
  
  if (nrow(heat_wide) > 1 && ncol(heat_wide) > 3) {
    heat_mat <- as.matrix(heat_wide[, -(1:3)])
    mode(heat_mat) <- "numeric"
    rownames(heat_mat) <- paste0(heat_wide$target_input, " | ", heat_wide$Genes)
    
    # Keep targets with at least two finite temperature values.
    keep_rows <- rowSums(is.finite(heat_mat)) >= 2
    heat_mat <- heat_mat[keep_rows, , drop = FALSE]
    
    if (nrow(heat_mat) > 1) {
      row_score <- apply(heat_mat, 1, function(x) {
        x <- x[is.finite(x)]
        if (!length(x)) return(-Inf)
        max(abs(x), na.rm = TRUE)
      })
      heat_mat <- heat_mat[order(row_score, decreasing = TRUE), , drop = FALSE]
      
      # pheatmap cannot cluster rows containing NA/NaN/Inf.
      # Use 0 for missing values because 0 means no observed Drug/DMSO shift.
      heat_mat[!is.finite(heat_mat)] <- 0
      
      row_var <- apply(heat_mat, 1, var, na.rm = TRUE)
      heat_mat <- heat_mat[is.finite(row_var) & row_var > 0, , drop = FALSE]
      
      if (nrow(heat_mat) > 1) {
        png(file.path(output_dir, "Target_heatmap_mean_log2_Drug_vs_DMSO.png"),
            width = 1800, height = max(1200, 28 * nrow(heat_mat)), res = 200)
        pheatmap(
          heat_mat,
          scale = "none",
          cluster_rows = TRUE,
          cluster_cols = FALSE,
          border_color = NA,
          main = "Matched targets: mean normalized log2(Drug/DMSO) across temperatures",
          fontsize_row = 7,
          fontsize_col = 10,
          na_col = "grey90"
        )
        dev.off()
      } else {
        msg("  Heatmap skipped: fewer than two variable target rows after cleaning.")
      }
    } else {
      msg("  Heatmap skipped: fewer than two targets with enough finite values.")
    }
  } else {
    msg("  Heatmap skipped: insufficient rows or temperature columns.")
  }
}

#############################################################################################################
# SUMMARY SCATTER
#############################################################################################################

if (nrow(target_summary) > 0) {
  y_col <- if (!is.na(score_col_A) && score_col_A %in% names(target_summary)) score_col_A else score_col_B
  y_label <- if (!is.na(score_col_A) && score_col_A %in% names(target_summary)) {
    "Max |mean normalized log2(H/L)| — Track A"
  } else {
    "Max |mean log2(Drug/DMSO)| — Track B"
  }
  
  x_col <- if (!is.na(signed_col_A) && signed_col_A %in% names(target_summary)) signed_col_A else "deltaTm_mean"
  x_label <- if (!is.na(signed_col_A) && signed_col_A %in% names(target_summary)) {
    "Signed peak normalized log2(H/L) — Track A"
  } else {
    "Mean deltaTm (Drug - DMSO) [°C] — Track B"
  }
  
  if (!is.null(y_col) && !is.na(y_col) && y_col %in% names(target_summary) && x_col %in% names(target_summary)) {
    plot_dt <- target_summary[is.finite(get(x_col)) | is.finite(get(y_col))]
    
    if (nrow(plot_dt) > 0) {
      if (has_pvalue && any(is.finite(plot_dt$p_value_best))) {
        plot_dt[, sig := p_value_best < 0.05]
        p_scatter <- ggplot(plot_dt, aes(x = get(x_col), y = get(y_col), colour = sig, label = target_input)) +
          scale_colour_manual(values = c("TRUE" = "#d73027", "FALSE" = "grey50"),
                              name = "p < 0.05", na.value = "grey80")
      } else {
        p_scatter <- ggplot(plot_dt, aes(x = get(x_col), y = get(y_col), label = target_input))
      }
      
      p_scatter <- p_scatter +
        geom_point(size = 2, alpha = 0.85) +
        geom_text(size = 2.5, vjust = -0.5, check_overlap = TRUE) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
        theme_bw(base_size = 11) +
        labs(
          title = "Requested targets in CETSA-MS SILAC",
          subtitle = "Track A is prioritized; Track B deltaTm is supportive",
          x = x_label,
          y = y_label
        )
      
      ggsave(file.path(output_dir, "Target_summary_scatter.png"),
             p_scatter, width = 11, height = 8, dpi = plot_dpi, bg = "white")
    }
  }
}

#############################################################################################################
# CONSOLE SUMMARY
#############################################################################################################

msg("\n══════════════════════════════════════════════════════")
msg("Done. Output folder: ", output_dir)
msg("══════════════════════════════════════════════════════")
msg("Requested targets  : ", nrow(targets_dt))
msg("Matched            : ", nrow(matched_primary))
msg("Unmatched          : ", nrow(unmatched_targets))

if (nrow(unmatched_targets) > 0) {
  msg("Unmatched targets  : ",
      paste(head(unmatched_targets$target_input, 20), collapse = ", "),
      if (nrow(unmatched_targets) > 20) " ..." else "")
}

if (nrow(target_summary) > 0) {
  show_cols <- intersect(
    c(
      "target_input", "Genes", "Protein_group", "match_type",
      "primary_evidence", "stabilization_direction",
      "signed_peak_log2_HL", "max_abs_mean_log2_HL", "temp_of_max_abs_log2HL",
      "p_value_ratio_track", "deltaTm_mean", "deltaTm_mean_curve", "deltaAUC_mean",
      score_col_B, "p_value_best", "p_value_source", "rank_score"
    ),
    names(target_summary)
  )
  msg("\nTop matched targets, Track A-first ordering:")
  print(target_summary[seq_len(min(30, .N)), ..show_cols])
}

