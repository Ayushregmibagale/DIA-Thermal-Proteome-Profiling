#############################################################################################################
#
# CETSA-MS SILAC TARGET CLASSIFICATION SCRIPT  v2.0
#
# PURPOSE
# -------
# Classify a requested target list into:
#   - strong stabilization
#   - weak stabilization
#   - destabilization
#   - detected but flat/no effect
#   - not detected
#
# KEY CHANGES vs v1
# -----------------
#  - Auto-detects input directory (v2 then v1 then broad scan), same logic
#    as the target inspection script.
#  - Classification now uses BOTH Track A (SILAC ratio) and Track B (MaxLFQ
#    fraction) evidence when available, with Track A taking priority.
#  - p_value_best from v2 hit table is incorporated into classification:
#    strong/weak thresholds require p < p_threshold when a p-value exists.
#  - Column names updated for v2 hit table (max_abs_mean_log2_ratio_trackB,
#    max_abs_mean_log2_HL, p_value_best); falls back gracefully to v1 names.
#  - Heatmap uses Track A temp_ratio_summary when available.
#  - Summary scatter shows Track A score on Y-axis when available, coloured
#    by class, with p < 0.05 targets labelled more prominently.
#  - evidence_score accumulates Track A score in addition to Track B metrics.
#  - matched_targets_primary.tsv looked up in auto-detected dir, not hardcoded.
#
# INPUTS EXPECTED
# ---------------
# From the main CETSA SILAC pipeline output directory:
#   - CETSA_SILAC_hit_table.tsv.gz                         (required)
#   - temp_effect_summary.tsv.gz                           (required)
#   - Target_inspection/matched_targets_primary.tsv        (optional)
#   - temp_ratio_summary.tsv.gz                            (optional, Track A)
#
#############################################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})

#############################################################################################################
# USER SETTINGS
#############################################################################################################

# Leave as NULL to auto-detect, or set explicitly e.g. input_dir <- "CETSA_SILAC_output_v2"
input_dir <- NULL

# Classification thresholds (Track B – deltaTm / deltaAUC)
strong_deltaTm_threshold   <- 2.0
weak_deltaTm_threshold     <- 1.0
strong_deltaAUC_threshold  <- 1.5
weak_deltaAUC_threshold    <- 0.5
destab_deltaTm_threshold   <- -1.0
destab_deltaAUC_threshold  <- -0.5

# Classification thresholds (Track B – MaxLFQ-derived log2 ratio)
strong_log2ratio_B_threshold <- 0.35
weak_log2ratio_B_threshold   <- 0.20
destab_log2ratio_B_threshold <- -0.20

# Classification thresholds (Track A – intra-sample SILAC log2(H/L))
strong_log2ratio_A_threshold <- 0.35   # same scale as Track B ratios
weak_log2ratio_A_threshold   <- 0.20
destab_log2ratio_A_threshold <- -0.20

# p-value gate: when a p_value_best is available, require it to be below
# this threshold for strong/weak classification (set to 1 to disable).
p_threshold <- 0.05

plot_width  <- 10
plot_height <- 7
plot_dpi    <- 300

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
    "VEGFR2"="KDR","VEGFR1"="FLT1","VEGFR3"="FLT4","TIE2"="TEK",
    "TRKA"="NTRK1","TRKB"="NTRK2","TRKC"="NTRK3",
    "PYK2"="PTK2B","FAK"="PTK2","BRK"="PTK6","FRK"="PTK5",
    "YES"="YES1","LZK"="MAP3K13","DLK"="MAP3K12","HPK1"="MAP4K1",
    "MEK1"="MAP2K1","MEK2"="MAP2K2","MEK3"="MAP2K3",
    "MEK4"="MAP2K4","MEK5"="MAP2K5","MEK6"="MAP2K6",
    "ERK1"="MAPK3","ERK2"="MAPK1","ERK3"="MAPK6",
    "ERK4"="MAPK4","ERK5"="MAPK7","ERK8"="MAPK15",
    "JNK1"="MAPK8","JNK2"="MAPK9","JNK3"="MAPK10",
    "p38-alpha"="MAPK14","p38-beta"="MAPK11",
    "p38-gamma"="MAPK12","p38-delta"="MAPK13",
    "IKK-alpha"="CHUK","IKK-beta"="IKBKB","IKK-epsilon"="IKBKE",
    "PKAC-alpha"="PRKACA","PKAC-beta"="PRKACB",
    "PFTK1"="CDK14","PFTAIRE2"="CDK15",
    "PCTK1"="CDK16","PCTK2"="CDK17","PCTK3"="CDK18",
    "CDC2L1"="CDK11B","CDC2L2"="CDK11A","CDC2L5"="CDK20",
    "DMPK2"="CDC42BPG","MRCKA"="CDC42BPA","MRCKB"="CDC42BPB"
  )
  if (y %in% names(alias_map)) alias_map[[y]] else y
}

split_gene_tokens <- function(x) {
  if (is.na(x) || x == "") return(character())
  parts <- unlist(strsplit(toupper(x), ";|,|\\|"))
  trimws(parts[nzchar(trimws(parts))])
}

fmt <- function(x, digits = 2) ifelse(is.finite(x), round(x, digits), "NA")

#############################################################################################################
# AUTO-DETECT INPUT DIRECTORY
#############################################################################################################

.sentinel   <- "CETSA_SILAC_hit_table.tsv.gz"
.candidates <- c("CETSA_SILAC_output_v2", "CETSA_SILAC_output")

if (is.null(input_dir)) {
  .found <- Filter(function(d) file.exists(file.path(d, .sentinel)), .candidates)
  
  if (!length(.found)) {
    .subdirs <- list.dirs(".", recursive = FALSE, full.names = FALSE)
    .found   <- Filter(function(d) file.exists(file.path(d, .sentinel)), .subdirs)
  }
  
  if (!length(.found))
    stop(
      "Cannot find a CETSA pipeline output directory.\n",
      "Expected '", .sentinel, "' inside one of: ",
      paste(.candidates, collapse = ", "), "\n",
      "Either run the main analysis script first, or set input_dir manually."
    )
  
  if (length(.found) > 1) {
    .pref  <- intersect(.candidates, .found)
    .found <- if (length(.pref)) .pref else .found
    message("Multiple CETSA output directories found: ", paste(.found, collapse = ", "),
            "\nUsing: ", .found[1])
  }
  
  input_dir <- .found[1]
}

message("Using input directory: ", input_dir)

output_dir <- file.path(input_dir, "Target_classification")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#############################################################################################################
# READ DATA
#############################################################################################################

msg("Reading input files...")

hit_dt  <- fread(file.path(input_dir, "CETSA_SILAC_hit_table.tsv.gz"))
temp_dt <- fread(file.path(input_dir, "temp_effect_summary.tsv.gz"))

# Track A temperature summary (v2 only; graceful skip for v1)
ratio_temp_path <- file.path(input_dir, "temp_ratio_summary.tsv.gz")
ratio_temp_dt   <- if (file.exists(ratio_temp_path)) {
  msg("  Reading Track A ratio temp summary.")
  fread(ratio_temp_path)
} else {
  msg("  temp_ratio_summary.tsv.gz not found – heatmap will use Track B data.")
  data.table()
}

# Detect v1 vs v2 column names in hit table
has_score_B  <- "max_abs_mean_log2_ratio_trackB" %in% names(hit_dt)
has_score_A  <- "max_abs_mean_log2_HL"            %in% names(hit_dt)
has_pvalue   <- "p_value_best"                     %in% names(hit_dt)

score_col_B <- if (has_score_B) "max_abs_mean_log2_ratio_trackB" else
  if ("max_abs_mean_log2_ratio" %in% names(hit_dt)) "max_abs_mean_log2_ratio" else
    NA_character_
score_col_A <- if (has_score_A) "max_abs_mean_log2_HL" else NA_character_

msg(sprintf(
  "Hit table: Track A score col = %s | Track B score col = %s | p-value col = %s",
  ifelse(is.na(score_col_A), "absent", score_col_A),
  ifelse(is.na(score_col_B), "absent", score_col_B),
  ifelse(has_pvalue, "p_value_best", "absent")
))

#############################################################################################################
# BUILD / LOAD MATCH TABLE
#############################################################################################################

matched_file <- file.path(input_dir, "Target_inspection", "matched_targets_primary.tsv")

if (file.exists(matched_file)) {
  msg("Using existing match table: ", matched_file)
  matched_primary <- fread(matched_file)
} else {
  msg("No previous match table found. Rebuilding matches...")
  
  targets_dt <- data.table(target_input = unique(targets))
  targets_dt[, target_clean := clean_target_name(target_input)]
  targets_dt[, target_base  := vapply(target_input, extract_base_symbol, character(1))]
  targets_dt[, target_gene  := vapply(target_input, map_target_alias, character(1))]
  
  hit_dt[, Genes         := as.character(Genes)]
  hit_dt[, Protein_group := as.character(Protein_group)]
  
  expanded_hits <- hit_dt[, {
    toks <- split_gene_tokens(Genes)
    if (!length(toks)) toks <- NA_character_
    .(gene_token = toks)
  }, by = .(Protein_group, Genes)]
  expanded_hits[, gene_token := toupper(gene_token)]
  
  match_gene <- merge(targets_dt, expanded_hits,
                      by.x = "target_gene", by.y = "gene_token",
                      all.x = TRUE, allow.cartesian = TRUE)
  match_gene[, match_type := fifelse(!is.na(Protein_group), "gene_match", NA_character_)]
  
  pg_upper <- unique(hit_dt[, .(Protein_group, pg_upper = toupper(Protein_group))])
  match_pg <- merge(targets_dt, pg_upper,
                    by.x = "target_clean", by.y = "pg_upper",
                    all.x = TRUE, allow.cartesian = TRUE)
  match_pg <- merge(match_pg, hit_dt[, .(Protein_group, Genes)],
                    by = "Protein_group", all.x = TRUE)
  match_pg[, match_type := fifelse(!is.na(Genes), "protein_group_match", NA_character_)]
  
  matched_all <- unique(rbindlist(list(
    match_gene[, .(target_input, target_clean, target_base, target_gene,
                   Protein_group, Genes, match_type)],
    match_pg[,   .(target_input, target_clean, target_base, target_gene,
                   Protein_group, Genes, match_type)]
  ), fill = TRUE)[!is.na(Protein_group)])
  
  mp <- c(gene_match = 1L, protein_group_match = 2L)
  matched_all[, match_priority := mp[match_type]]
  setorder(matched_all, target_input, match_priority, Protein_group)
  matched_primary <- matched_all[, .SD[1], by = target_input]
}

#############################################################################################################
# MASTER TARGET TABLE
#############################################################################################################

targets_master <- data.table(target_input = unique(targets))
targets_master <- merge(targets_master, matched_primary,
                        by = "target_input", all.x = TRUE)
targets_master <- merge(targets_master, hit_dt,
                        by = c("Protein_group", "Genes"), all.x = TRUE)

#############################################################################################################
# TEMPERATURE-SHAPE SUMMARY
# Pull max/min log2 ratio from Track A if available, else Track B
#############################################################################################################

heat_source_dt <- if (nrow(ratio_temp_dt) > 0) ratio_temp_dt else temp_dt

temp_target <- merge(
  matched_primary[, .(target_input, Protein_group, Genes)],
  heat_source_dt[, .(Protein_group, Genes, Temp, mean_log2_ratio)],
  by = c("Protein_group", "Genes"),
  all.x = TRUE
)

shape_summary <- temp_target[, {
  ok <- is.finite(mean_log2_ratio)
  .(
    n_temp_points  = sum(ok),
    max_ratio      = if (!any(ok)) NA_real_ else max(mean_log2_ratio[ok]),
    min_ratio      = if (!any(ok)) NA_real_ else min(mean_log2_ratio[ok]),
    max_abs_ratio  = if (!any(ok)) NA_real_ else max(abs(mean_log2_ratio[ok])),
    temp_at_max    = if (!any(ok)) NA_real_ else as.numeric(Temp[which.max(mean_log2_ratio)]),
    temp_at_min    = if (!any(ok)) NA_real_ else as.numeric(Temp[which.min(mean_log2_ratio)]),
    # Number of temperatures with a meaningfully negative ratio:
    # used to require consistency before calling destabilization
    n_neg_temps_B  = sum(ok & mean_log2_ratio < destab_log2ratio_B_threshold),
    n_neg_temps_A  = sum(ok & mean_log2_ratio < destab_log2ratio_A_threshold)
  )
}, by = .(target_input, Protein_group, Genes)]

targets_master <- merge(targets_master, shape_summary,
                        by = c("target_input", "Protein_group", "Genes"), all.x = TRUE)

# Convenience accessors for score columns (NA-safe)
.get <- function(dt, col) if (!is.na(col) && col %in% names(dt)) dt[[col]] else rep(NA_real_, nrow(dt))

targets_master[, .score_A := .get(targets_master, score_col_A)]
targets_master[, .score_B := .get(targets_master, score_col_B)]
targets_master[, .pval    := if (has_pvalue) p_value_best else NA_real_]

# p-value gate: TRUE when either no p-value exists or p passes threshold
targets_master[, .p_ok := !is.finite(.pval) | .pval < p_threshold]

#############################################################################################################
# CLASSIFICATION LOGIC
# Track A (SILAC ratio) evidence is primary when available; Track B is secondary.
# A p-value gate is applied when p_value_best is present.
#############################################################################################################

targets_master[, detected := !is.na(Protein_group)]

# ── Positive evidence ────────────────────────────────────────────────────────
targets_master[, strong_positive_evidence := detected & .p_ok & (
  (is.finite(deltaTm_mean)       & deltaTm_mean       >= strong_deltaTm_threshold)  |
    (is.finite(deltaTm_mean_curve) & deltaTm_mean_curve >= strong_deltaTm_threshold)  |
    (is.finite(deltaAUC_mean)      & deltaAUC_mean      >= strong_deltaAUC_threshold) |
    (is.finite(.score_A)           & .score_A            >= strong_log2ratio_A_threshold) |
    (is.finite(.score_B)           & max_ratio           >= strong_log2ratio_B_threshold)
)]

targets_master[, weak_positive_evidence := detected & .p_ok & (
  (is.finite(deltaTm_mean)       & deltaTm_mean       >= weak_deltaTm_threshold)  |
    (is.finite(deltaTm_mean_curve) & deltaTm_mean_curve >= weak_deltaTm_threshold)  |
    (is.finite(deltaAUC_mean)      & deltaAUC_mean      >= weak_deltaAUC_threshold) |
    (is.finite(.score_A)           & .score_A            >= weak_log2ratio_A_threshold) |
    (is.finite(.score_B)           & max_ratio           >= weak_log2ratio_B_threshold)
)]

# ── Negative evidence ────────────────────────────────────────────────────────
# Requires CONSISTENCY: a single dipping temperature is not enough.
# Ratio-based destabilization requires >= 2 temperatures below threshold.
# deltaTm / deltaAUC based criteria are unchanged (already aggregated metrics).
targets_master[, negative_evidence := (detected) & (
  (is.finite(deltaTm_mean)       & deltaTm_mean       <= destab_deltaTm_threshold)  |
    (is.finite(deltaTm_mean_curve) & deltaTm_mean_curve <= destab_deltaTm_threshold)  |
    (is.finite(deltaAUC_mean)      & deltaAUC_mean      <= destab_deltaAUC_threshold) |
    (is.finite(n_neg_temps_A)      & n_neg_temps_A      >= 2L) |
    (is.finite(n_neg_temps_B)      & n_neg_temps_B      >= 2L)
)]

# ── Assign class ─────────────────────────────────────────────────────────────
targets_master[, target_class := fcase(
  !detected,                                         "not detected",
  strong_positive_evidence & !negative_evidence,     "strong stabilization",
  negative_evidence        & !strong_positive_evidence, "destabilization",
  weak_positive_evidence   & !negative_evidence,     "weak stabilization",
  default =                                          "detected but flat/no effect"
)]

targets_master[, target_class := factor(
  target_class,
  levels = c("strong stabilization","weak stabilization",
             "destabilization","detected but flat/no effect","not detected")
)]

# Clean up internal helper columns
targets_master[, c(".score_A",".score_B",".pval",".p_ok") := NULL]

#############################################################################################################
# HUMAN-READABLE REASON
#############################################################################################################

build_reason <- function(prefix, dt) {
  paste0(
    prefix,
    "deltaTm_mean=",        fmt(dt$deltaTm_mean),
    ", deltaTm_curve=",     fmt(dt$deltaTm_mean_curve),
    ", deltaAUC=",          fmt(dt$deltaAUC_mean, 3),
    if (!is.na(score_col_A) && score_col_A %in% names(dt))
      paste0(", scoreA(", score_col_A, ")=", fmt(dt[[score_col_A]], 3)) else "",
    if (!is.na(score_col_B) && score_col_B %in% names(dt))
      paste0(", scoreB(", score_col_B, ")=", fmt(dt[[score_col_B]], 3)) else "",
    if (has_pvalue) paste0(", p=", signif(dt$p_value_best, 2)) else ""
  )
}

targets_master[, classification_reason := ""]
targets_master[detected == FALSE,
               classification_reason := "No CETSA match found in hit table"]
targets_master[target_class == "strong stabilization",
               classification_reason := build_reason("Strong positive CETSA evidence; ", .SD)]
targets_master[target_class == "weak stabilization",
               classification_reason := build_reason("Moderate positive CETSA evidence; ", .SD)]
targets_master[target_class == "destabilization",
               classification_reason := build_reason("Negative CETSA shift; ", .SD)]
targets_master[target_class == "detected but flat/no effect",
               classification_reason := build_reason("Detected but no strong CETSA shift; ", .SD)]

#############################################################################################################
# ORDER & EVIDENCE SCORE
#############################################################################################################

targets_master[, class_priority := as.integer(target_class)]

targets_master[, evidence_score := 0]
targets_master[is.finite(deltaTm_mean),       evidence_score := evidence_score + abs(deltaTm_mean)]
targets_master[is.finite(deltaTm_mean_curve), evidence_score := evidence_score + abs(deltaTm_mean_curve)]
targets_master[is.finite(deltaAUC_mean),      evidence_score := evidence_score + abs(deltaAUC_mean)]
if (!is.na(score_col_A) && score_col_A %in% names(targets_master))
  targets_master[is.finite(get(score_col_A)), evidence_score := evidence_score + abs(get(score_col_A))]
if (!is.na(score_col_B) && score_col_B %in% names(targets_master))
  targets_master[is.finite(get(score_col_B)), evidence_score := evidence_score + abs(get(score_col_B))]

setorder(targets_master, class_priority, -evidence_score, target_input)

#############################################################################################################
# WRITE OUTPUT TABLES
#############################################################################################################

msg("Writing output tables...")

fwrite(targets_master, file.path(output_dir, "target_classification_full.tsv"), sep = "\t")

summary_cols <- intersect(
  c("target_input","target_class","Genes","Protein_group","match_type",
    "deltaTm_mean","deltaTm_mean_curve","deltaAUC_mean",
    score_col_A, score_col_B, "p_value_best", "classification_reason"),
  names(targets_master)
)
fwrite(targets_master[, ..summary_cols],
       file.path(output_dir, "target_classification_summary.tsv"), sep = "\t")

for (cls in levels(targets_master$target_class)) {
  safe_name <- gsub("[^A-Za-z0-9_]", "_", cls)
  fwrite(targets_master[target_class == cls],
         file.path(output_dir, paste0("targets_", safe_name, ".tsv")), sep = "\t")
}

#############################################################################################################
# CLASS COUNT PLOT
#############################################################################################################

class_counts <- targets_master[, .N, by = target_class]
class_colours <- c(
  "strong stabilization"        = "#d73027",
  "weak stabilization"          = "#f46d43",
  "destabilization"             = "#4575b4",
  "detected but flat/no effect" = "#878787",
  "not detected"                = "#cccccc"
)

p_counts <- ggplot(class_counts, aes(x = target_class, y = N, fill = target_class)) +
  geom_col() +
  geom_text(aes(label = N), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = class_colours) +
  theme_bw(base_size = 11) +
  labs(title = "Requested targets by CETSA class", x = NULL, y = "Count") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")

ggsave(file.path(output_dir, "target_class_counts.png"), p_counts,
       width = 9, height = 6, dpi = plot_dpi, bg = "white")

#############################################################################################################
# SUMMARY SCATTER
# x = deltaTm (Track B), y = Track A score if available, else Track B score
#############################################################################################################

y_col   <- if (!is.na(score_col_A) && score_col_A %in% names(targets_master)) score_col_A else score_col_B
y_label <- if (!is.na(score_col_A) && score_col_A %in% names(targets_master))
  "Max |mean log2(H/L)| — Track A (intra-sample SILAC)"  else
    "Max |mean log2(Drug/DMSO)| — Track B (MaxLFQ-derived)"

if (!is.null(y_col) && !is.na(y_col) && y_col %in% names(targets_master)) {
  scatter_dt <- targets_master[
    (detected) & (is.finite(deltaTm_mean) | is.finite(get(y_col)))
  ]
  
  # label only significant or strongly shifted points to avoid overplotting
  scatter_dt[, show_label := (
    (is.finite(deltaTm_mean)   & abs(deltaTm_mean)   >= weak_deltaTm_threshold) |
      (is.finite(get(y_col))     & get(y_col)           >= weak_log2ratio_A_threshold)
  )]
  if (has_pvalue) scatter_dt[, show_label := show_label | (is.finite(p_value_best) & p_value_best < p_threshold)]
  
  p_scatter <- ggplot(scatter_dt,
                      aes(x = deltaTm_mean, y = get(y_col),
                          colour = target_class, label = target_input)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_point(size = 2.2, alpha = 0.85) +
    ggplot2::geom_text(
      data = scatter_dt[show_label == TRUE],
      size = 2.4, vjust = -0.5, check_overlap = TRUE
    ) +
    scale_colour_manual(values = class_colours, name = "CETSA class") +
    theme_bw(base_size = 11) +
    labs(title = "Target classification summary",
         x = "Mean deltaTm (Drug – DMSO) [°C]",
         y = y_label)
  
  ggsave(file.path(output_dir, "target_classification_scatter.png"), p_scatter,
         width = 11, height = 8, dpi = plot_dpi, bg = "white")
}

#############################################################################################################
# HEATMAP OF DETECTED TARGETS
# Use Track A ratio data when available
#############################################################################################################

heat_temp_dt <- if (nrow(ratio_temp_dt) > 0) {
  msg("Heatmap using Track A ratio data.")
  ratio_temp_dt
} else {
  msg("Heatmap using Track B temp_effect_summary.")
  temp_dt
}

detected_temp <- merge(
  targets_master[(detected), .(target_input, target_class, Protein_group, Genes)],
  heat_temp_dt[, .(Protein_group, Genes, Temp, mean_log2_ratio)],
  by = c("Protein_group", "Genes"),
  all.x = TRUE
)

if (nrow(detected_temp) > 0) {
  heat_wide <- dcast(
    detected_temp,
    target_input + target_class + Genes + Protein_group ~ Temp,
    value.var = "mean_log2_ratio"
  )
  
  if (nrow(heat_wide) > 1) {
    heat_mat <- as.matrix(heat_wide[, -(1:4)])
    rownames(heat_mat) <- paste0(heat_wide$target_input,
                                 " | ", heat_wide$target_class)
    
    # Annotation bar: class colour per row
    ann_row <- data.frame(
      Class = as.character(heat_wide$target_class),
      row.names = rownames(heat_mat)
    )
    ann_colours <- list(Class = class_colours[names(class_colours) %in% ann_row$Class])
    
    png(file.path(output_dir, "target_classes_heatmap.png"),
        width = 1800, height = max(1200, 26 * nrow(heat_mat)), res = 200)
    pheatmap(
      heat_mat,
      annotation_row = ann_row,
      annotation_colors = ann_colours,
      scale = "none", cluster_rows = TRUE, cluster_cols = FALSE,
      border_color = NA,
      main = "Detected targets: mean log2(Drug/DMSO) across temperatures",
      fontsize_row = 7, fontsize_col = 10
    )
    dev.off()
  }
}

#############################################################################################################
# CONSOLE SUMMARY
#############################################################################################################

msg("\n══════════════════════════════════════════════════════")
msg("Done.  Output folder: ", output_dir)
msg("══════════════════════════════════════════════════════")
msg("Requested targets : ", uniqueN(targets))
msg("Detected          : ", targets_master[(detected),  .N])
msg("Not detected      : ", targets_master[(!detected), .N])

msg("\nClass counts:")
print(targets_master[, .N, by = target_class][order(target_class)])

show_cols <- intersect(
  c("target_input","Genes","Protein_group","deltaTm_mean",
    "deltaTm_mean_curve","deltaAUC_mean",
    score_col_A, score_col_B, "p_value_best"),
  names(targets_master)
)

msg("\nStrong stabilization targets:")
print(targets_master[target_class == "strong stabilization", ..show_cols][seq_len(min(30L, .N))])

msg("\nDestabilization targets:")
print(targets_master[target_class == "destabilization", ..show_cols][seq_len(min(20L, .N))])