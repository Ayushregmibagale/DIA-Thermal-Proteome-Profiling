#############################################################################################################
#
# CETSA-MS SILAC TARGET CLASSIFICATION SCRIPT v3.0
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
# This version is aligned with the raw-LFQ CETSA-MS SILAC analysis output:
#
#   Main analysis output directory:
#     CETSA_SILAC_output_v3_rawLFQ
#
#   Track A = SILAC H/L ratio, primary evidence
#     - signed_peak_log2_HL > 0  => stabilization
#     - signed_peak_log2_HL < 0  => destabilization
#     - max_abs_mean_log2_HL gives magnitude
#
#   Track B = raw-LFQ melt curves, secondary/supportive evidence
#     - deltaTm_mean
#     - deltaTm_mean_curve
#     - deltaAUC_mean
#     - max_abs_mean_log2_ratio_trackB
#
# IMPORTANT INTERPRETATION
# ------------------------
# Track B now uses raw Precursor.Quantity-derived LFQ_L/LFQ_H from the upstream
# protein-quant script. It is intentionally not globally run-normalized. Therefore
# Track B can look noisier and should not overrule clear Track A SILAC evidence.
#
# P-value handling
# ----------------
# p_value_best is used when present.
# p_value_source is preserved when present.
# If a p-value exists, strong/weak calls require p < p_threshold unless p_threshold = 1.
#
# INPUTS EXPECTED
# ---------------
# Required, from the main CETSA SILAC pipeline output directory:
#   - CETSA_SILAC_hit_table.tsv.gz
#   - temp_effect_summary.tsv.gz
#
# Optional:
#   - Target_inspection/matched_targets_primary.tsv
#   - temp_ratio_summary.tsv.gz
#
# OUTPUT
# ------
#   - target_classification_full.tsv
#   - target_classification_summary.tsv
#   - targets_<class>.tsv
#   - target_class_counts.png
#   - target_classification_scatter.png
#   - target_classes_heatmap.png
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

# Leave NULL to auto-detect, or set explicitly:
# input_dir <- "CETSA_SILAC_output_v3_rawLFQ"
input_dir <- NULL

# Track A thresholds: normalized intra-sample SILAC log2(H/L)
strong_log2ratio_A_threshold <- 0.35
weak_log2ratio_A_threshold   <- 0.20
destab_log2ratio_A_threshold <- -0.20

# Track B thresholds: raw-LFQ deltaTm / deltaAUC
strong_deltaTm_threshold   <- 2.0
weak_deltaTm_threshold     <- 1.0
strong_deltaAUC_threshold  <- 1.5
weak_deltaAUC_threshold    <- 0.5
destab_deltaTm_threshold   <- -1.0
destab_deltaAUC_threshold  <- -0.5

# Track B raw-LFQ fraction-derived log2 ratio thresholds
strong_log2ratio_B_threshold <- 0.35
weak_log2ratio_B_threshold   <- 0.20
destab_log2ratio_B_threshold <- -0.20

# p-value gate. Set to 1 to disable.
p_threshold <- 0.05

plot_width  <- 10
plot_height <- 7
plot_dpi    <- 300

#############################################################################################################
# TARGET LIST
#############################################################################################################

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
    "VEGFR2" = "KDR", "VEGFR1" = "FLT1", "VEGFR3" = "FLT4", "TIE2" = "TEK",
    "TRKA" = "NTRK1", "TRKB" = "NTRK2", "TRKC" = "NTRK3",
    "PYK2" = "PTK2B", "FAK" = "PTK2", "BRK" = "PTK6", "FRK" = "PTK5",
    "YES" = "YES1", "LZK" = "MAP3K13", "DLK" = "MAP3K12", "HPK1" = "MAP4K1",
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

fmt <- function(x, digits = 2) ifelse(is.finite(x), as.character(round(x, digits)), "NA")

get_col_or_na <- function(dt, col, type = "numeric") {
  if (!is.na(col) && col %in% names(dt)) return(dt[[col]])
  if (type == "character") return(rep(NA_character_, nrow(dt)))
  rep(NA_real_, nrow(dt))
}

#############################################################################################################
# AUTO-DETECT INPUT DIRECTORY
#############################################################################################################

.sentinel <- "CETSA_SILAC_hit_table.tsv.gz"
.candidates <- c(
  "CETSA_SILAC_output_v3_rawLFQ",
  "CETSA_SILAC_output_v2_1",
  "CETSA_SILAC_output_v2",
  "CETSA_SILAC_output"
)

if (is.null(input_dir)) {
  .found <- Filter(function(d) file.exists(file.path(d, .sentinel)), .candidates)
  
  if (!length(.found)) {
    .subdirs <- list.dirs(".", recursive = FALSE, full.names = FALSE)
    .found <- Filter(function(d) file.exists(file.path(d, .sentinel)), .subdirs)
  }
  
  if (!length(.found)) {
    stop(
      "Cannot find a CETSA pipeline output directory.\n",
      "Expected '", .sentinel, "' inside one of: ", paste(.candidates, collapse = ", "), "\n",
      "Either run the main analysis script first, or set input_dir manually."
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

output_dir <- file.path(input_dir, "Target_classification")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#############################################################################################################
# READ DATA
#############################################################################################################

msg("Reading input files...")

hit_file <- file.path(input_dir, "CETSA_SILAC_hit_table.tsv.gz")
temp_file <- file.path(input_dir, "temp_effect_summary.tsv.gz")

if (!file.exists(hit_file)) stop("Missing required file: ", hit_file)
if (!file.exists(temp_file)) stop("Missing required file: ", temp_file)

hit_dt  <- fread(hit_file)
temp_dt <- fread(temp_file)

ratio_temp_path <- file.path(input_dir, "temp_ratio_summary.tsv.gz")
ratio_temp_dt <- if (file.exists(ratio_temp_path)) {
  msg("  Reading Track A ratio temp summary.")
  fread(ratio_temp_path)
} else {
  msg("  temp_ratio_summary.tsv.gz not found; heatmap will use Track B raw-LFQ data.")
  data.table()
}

score_col_A <- if ("max_abs_mean_log2_HL" %in% names(hit_dt)) "max_abs_mean_log2_HL" else NA_character_
signed_col_A <- if ("signed_peak_log2_HL" %in% names(hit_dt)) "signed_peak_log2_HL" else NA_character_
score_col_B <- if ("max_abs_mean_log2_ratio_trackB" %in% names(hit_dt)) {
  "max_abs_mean_log2_ratio_trackB"
} else if ("max_abs_mean_log2_ratio" %in% names(hit_dt)) {
  "max_abs_mean_log2_ratio"
} else NA_character_

has_pvalue <- "p_value_best" %in% names(hit_dt)
has_pvalue_source <- "p_value_source" %in% names(hit_dt)

if (!("primary_evidence" %in% names(hit_dt)) && "evidence_track" %in% names(hit_dt)) {
  setnames(hit_dt, "evidence_track", "primary_evidence")
}
if (!("primary_evidence" %in% names(hit_dt))) hit_dt[, primary_evidence := NA_character_]
if (!("stabilization_direction" %in% names(hit_dt))) hit_dt[, stabilization_direction := NA_character_]
if (!("rank_score" %in% names(hit_dt))) hit_dt[, rank_score := NA_real_]
if (!has_pvalue) hit_dt[, p_value_best := NA_real_]
if (!has_pvalue_source) hit_dt[, p_value_source := NA_character_]

msg(sprintf(
  "Hit table columns: Track A score = %s | Track A signed = %s | Track B raw-LFQ score = %s | p-value = %s",
  ifelse(is.na(score_col_A), "absent", score_col_A),
  ifelse(is.na(signed_col_A), "absent", signed_col_A),
  ifelse(is.na(score_col_B), "absent", score_col_B),
  ifelse("p_value_best" %in% names(hit_dt), "p_value_best", "absent")
))

#############################################################################################################
# BUILD / LOAD MATCH TABLE
#############################################################################################################

matched_file <- file.path(input_dir, "Target_inspection", "matched_targets_primary.tsv")

targets_dt <- data.table(target_input = unique(targets))
targets_dt[, target_clean := clean_target_name(target_input)]
targets_dt[, target_base  := vapply(target_input, extract_base_symbol, character(1))]
targets_dt[, target_gene  := vapply(target_input, map_target_alias, character(1))]

if (file.exists(matched_file)) {
  msg("Using existing match table: ", matched_file)
  matched_primary <- fread(matched_file)
} else {
  msg("No previous match table found. Rebuilding matches...")
  
  hit_dt[, Genes := as.character(Genes)]
  hit_dt[, Protein_group := as.character(Protein_group)]
  
  expanded_hits <- hit_dt[, {
    toks <- split_gene_tokens(Genes)
    if (!length(toks)) toks <- NA_character_
    .(gene_token = toks)
  }, by = .(Protein_group, Genes)]
  expanded_hits[, gene_token := toupper(gene_token)]
  
  match_gene <- merge(
    targets_dt, expanded_hits,
    by.x = "target_gene", by.y = "gene_token",
    all.x = TRUE, allow.cartesian = TRUE
  )
  match_gene[, match_type := fifelse(!is.na(Protein_group), "gene_match", NA_character_)]
  
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
  
  mp <- c(gene_match = 1L, protein_group_match = 2L)
  matched_all[, match_priority := mp[match_type]]
  setorder(matched_all, target_input, match_priority, Protein_group)
  matched_primary <- matched_all[, .SD[1], by = target_input]
}

#############################################################################################################
# MASTER TARGET TABLE
#############################################################################################################

targets_master <- data.table(target_input = unique(targets))
targets_master <- merge(targets_master, matched_primary, by = "target_input", all.x = TRUE)
targets_master <- merge(targets_master, hit_dt, by = c("Protein_group", "Genes"), all.x = TRUE)

#############################################################################################################
# TEMPERATURE-SHAPE SUMMARY
#############################################################################################################

heat_source_dt <- if (nrow(ratio_temp_dt) > 0 && "mean_log2_ratio" %in% names(ratio_temp_dt)) ratio_temp_dt else temp_dt
heat_source_label <- if (nrow(ratio_temp_dt) > 0 && "mean_log2_ratio" %in% names(ratio_temp_dt)) "Track A" else "Track B raw-LFQ"

if (nrow(matched_primary) > 0) {
  temp_target <- merge(
    matched_primary[, .(target_input, Protein_group, Genes)],
    heat_source_dt[, .(Protein_group, Genes, Temp, mean_log2_ratio)],
    by = c("Protein_group", "Genes"),
    all.x = TRUE
  )
  
  shape_summary <- temp_target[, {
    ok <- is.finite(mean_log2_ratio)
    if (!any(ok)) {
      .(
        n_temp_points = 0L,
        max_ratio = NA_real_, min_ratio = NA_real_, max_abs_ratio = NA_real_,
        temp_at_max = NA_real_, temp_at_min = NA_real_,
        n_pos_temps_A = 0L, n_neg_temps_A = 0L,
        n_pos_temps_B = 0L, n_neg_temps_B = 0L
      )
    } else {
      vals <- mean_log2_ratio[ok]
      temps <- as.numeric(Temp[ok])
      .(
        n_temp_points = length(vals),
        max_ratio = max(vals, na.rm = TRUE),
        min_ratio = min(vals, na.rm = TRUE),
        max_abs_ratio = max(abs(vals), na.rm = TRUE),
        temp_at_max = temps[which.max(vals)],
        temp_at_min = temps[which.min(vals)],
        n_pos_temps_A = sum(vals > weak_log2ratio_A_threshold, na.rm = TRUE),
        n_neg_temps_A = sum(vals < destab_log2ratio_A_threshold, na.rm = TRUE),
        n_pos_temps_B = sum(vals > weak_log2ratio_B_threshold, na.rm = TRUE),
        n_neg_temps_B = sum(vals < destab_log2ratio_B_threshold, na.rm = TRUE)
      )
    }
  }, by = .(target_input, Protein_group, Genes)]
  
  targets_master <- merge(targets_master, shape_summary,
                          by = c("target_input", "Protein_group", "Genes"), all.x = TRUE)
}

for (cn in c("n_temp_points", "n_pos_temps_A", "n_neg_temps_A", "n_pos_temps_B", "n_neg_temps_B")) {
  if (!(cn %in% names(targets_master))) targets_master[, (cn) := 0L]
  targets_master[is.na(get(cn)), (cn) := 0L]
}
for (cn in c("max_ratio", "min_ratio", "max_abs_ratio", "temp_at_max", "temp_at_min")) {
  if (!(cn %in% names(targets_master))) targets_master[, (cn) := NA_real_]
}

#############################################################################################################
# CONVENIENCE COLUMNS
#############################################################################################################

targets_master[, detected := !is.na(Protein_group)]
targets_master[, score_A_abs := get_col_or_na(targets_master, score_col_A)]
targets_master[, score_A_signed := get_col_or_na(targets_master, signed_col_A)]

targets_master[
  !is.finite(score_A_signed) & is.finite(score_A_abs) & is.finite(max_ratio) & is.finite(min_ratio),
  score_A_signed := fifelse(abs(max_ratio) >= abs(min_ratio), max_ratio, min_ratio)
]

targets_master[, score_B_abs := get_col_or_na(targets_master, score_col_B)]
targets_master[, pval_for_classification := p_value_best]
targets_master[, p_ok := !is.finite(pval_for_classification) | pval_for_classification < p_threshold]

# Track A evidence. Strong/weak positive require p_ok. Negative does not require p_ok,
# because destabilization can also be a real phenotype with noisier p-values in 2-rep designs.
targets_master[, trackA_strong_pos := detected & p_ok & is.finite(score_A_signed) & score_A_signed >= strong_log2ratio_A_threshold]
targets_master[, trackA_weak_pos   := detected & p_ok & is.finite(score_A_signed) & score_A_signed >= weak_log2ratio_A_threshold]
targets_master[, trackA_neg        := detected & is.finite(score_A_signed) & score_A_signed <= destab_log2ratio_A_threshold]
targets_master[, trackA_neg_shape  := detected & !is.finite(score_A_signed) & n_neg_temps_A >= 2L]

# Track B raw-LFQ evidence.
targets_master[, trackB_strong_pos := detected & p_ok & (
  (is.finite(deltaTm_mean)       & deltaTm_mean       >= strong_deltaTm_threshold) |
    (is.finite(deltaTm_mean_curve) & deltaTm_mean_curve >= strong_deltaTm_threshold) |
    (is.finite(deltaAUC_mean)      & deltaAUC_mean      >= strong_deltaAUC_threshold) |
    (is.finite(max_ratio)          & max_ratio          >= strong_log2ratio_B_threshold)
)]

targets_master[, trackB_weak_pos := detected & p_ok & (
  (is.finite(deltaTm_mean)       & deltaTm_mean       >= weak_deltaTm_threshold) |
    (is.finite(deltaTm_mean_curve) & deltaTm_mean_curve >= weak_deltaTm_threshold) |
    (is.finite(deltaAUC_mean)      & deltaAUC_mean      >= weak_deltaAUC_threshold) |
    (is.finite(max_ratio)          & max_ratio          >= weak_log2ratio_B_threshold)
)]

targets_master[, trackB_neg := detected & (
  (is.finite(deltaTm_mean)       & deltaTm_mean       <= destab_deltaTm_threshold) |
    (is.finite(deltaTm_mean_curve) & deltaTm_mean_curve <= destab_deltaTm_threshold) |
    (is.finite(deltaAUC_mean)      & deltaAUC_mean      <= destab_deltaAUC_threshold) |
    (is.finite(min_ratio)          & min_ratio          <= destab_log2ratio_B_threshold & n_neg_temps_B >= 2L)
)]

#############################################################################################################
# CLASSIFICATION LOGIC
#############################################################################################################

# Track A takes priority. Track B is used only when Track A signed evidence is absent,
# or when it supports the same direction.
targets_master[, target_class := fcase(
  !detected,
  "not detected",
  
  trackA_strong_pos,
  "strong stabilization",
  
  trackA_neg | trackA_neg_shape,
  "destabilization",
  
  trackA_weak_pos,
  "weak stabilization",
  
  !is.finite(score_A_signed) & trackB_strong_pos & !trackB_neg,
  "strong stabilization",
  
  !is.finite(score_A_signed) & trackB_neg,
  "destabilization",
  
  !is.finite(score_A_signed) & trackB_weak_pos & !trackB_neg,
  "weak stabilization",
  
  default = "detected but flat/no effect"
)]

targets_master[, evidence_conflict := fifelse(
  detected &
    (((trackA_strong_pos | trackA_weak_pos) & trackB_neg) |
       ((trackA_neg | trackA_neg_shape) & (trackB_strong_pos | trackB_weak_pos))),
  TRUE, FALSE
)]

targets_master[, target_class := factor(
  target_class,
  levels = c("strong stabilization", "weak stabilization", "destabilization",
             "detected but flat/no effect", "not detected")
)]

#############################################################################################################
# HUMAN-READABLE REASON
#############################################################################################################

build_reason_row <- function(row) {
  if (!isTRUE(row$detected)) return("No CETSA match found in hit table")
  
  parts <- c(
    paste0("classification=", as.character(row$target_class)),
    paste0("primary_evidence=", ifelse(is.na(row$primary_evidence), "NA", row$primary_evidence)),
    paste0("direction=", ifelse(is.na(row$stabilization_direction), "NA", row$stabilization_direction)),
    paste0("TrackA_signed_peak_log2HL=", fmt(row$score_A_signed, 3)),
    paste0("TrackA_max_abs_log2HL=", fmt(row$score_A_abs, 3)),
    paste0("TrackB_rawLFQ_deltaTm=", fmt(row$deltaTm_mean, 2)),
    paste0("TrackB_rawLFQ_deltaTm_curve=", fmt(row$deltaTm_mean_curve, 2)),
    paste0("TrackB_rawLFQ_deltaAUC=", fmt(row$deltaAUC_mean, 3)),
    paste0("TrackB_rawLFQ_max_ratio=", fmt(row$max_ratio, 3)),
    paste0("TrackB_rawLFQ_min_ratio=", fmt(row$min_ratio, 3)),
    paste0("p=", ifelse(is.finite(row$p_value_best), signif(row$p_value_best, 2), "NA")),
    paste0("p_source=", ifelse(is.na(row$p_value_source), "NA", row$p_value_source)),
    paste0("heatmap_source=", heat_source_label),
    paste0("conflict=", row$evidence_conflict)
  )
  paste(parts, collapse = "; ")
}

targets_master[, classification_reason := vapply(seq_len(.N), function(i) build_reason_row(.SD[i]), character(1))]

#############################################################################################################
# ORDER AND EVIDENCE SCORE
#############################################################################################################

targets_master[, class_priority := as.integer(target_class)]

targets_master[, evidence_score := fcase(
  is.finite(score_A_abs), score_A_abs,
  is.finite(deltaTm_mean), abs(deltaTm_mean),
  is.finite(deltaTm_mean_curve), abs(deltaTm_mean_curve),
  is.finite(score_B_abs), score_B_abs,
  is.finite(max_abs_ratio), max_abs_ratio,
  default = 0
)]

targets_master[, p_sort := fifelse(is.finite(p_value_best), p_value_best, Inf)]
setorder(targets_master, class_priority, -evidence_score, p_sort, target_input)
targets_master[, p_sort := NULL]

#############################################################################################################
# WRITE OUTPUT TABLES
#############################################################################################################

msg("Writing output tables...")

fwrite(targets_master, file.path(output_dir, "target_classification_full.tsv"), sep = "\t")

summary_cols <- intersect(
  c(
    "target_input", "target_class", "Genes", "Protein_group", "match_type",
    "primary_evidence", "stabilization_direction", "evidence_conflict",
    "signed_peak_log2_HL", "max_abs_mean_log2_HL", "temp_of_max_abs_log2HL",
    "deltaTm_mean", "deltaTm_mean_curve", "deltaAUC_mean",
    score_col_B, "max_ratio", "min_ratio", "n_temp_points",
    "p_value_best", "p_value_source", "classification_reason"
  ),
  names(targets_master)
)
fwrite(targets_master[, ..summary_cols], file.path(output_dir, "target_classification_summary.tsv"), sep = "\t")

for (cls in levels(targets_master$target_class)) {
  safe_name <- gsub("[^A-Za-z0-9_]", "_", cls)
  fwrite(targets_master[target_class == cls], file.path(output_dir, paste0("targets_", safe_name, ".tsv")), sep = "\t")
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
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")

ggsave(file.path(output_dir, "target_class_counts.png"), p_counts,
       width = 9, height = 6, dpi = plot_dpi, bg = "white")

#############################################################################################################
# SUMMARY SCATTER
#############################################################################################################

x_col <- if (!is.na(signed_col_A) && signed_col_A %in% names(targets_master)) signed_col_A else "deltaTm_mean"
y_col <- if (!is.na(score_col_A) && score_col_A %in% names(targets_master)) score_col_A else score_col_B

x_label <- if (!is.na(signed_col_A) && signed_col_A %in% names(targets_master)) {
  "Signed peak normalized log2(H/L) — Track A"
} else {
  "Mean deltaTm (Drug - DMSO) [°C] — Track B raw-LFQ"
}
y_label <- if (!is.na(score_col_A) && score_col_A %in% names(targets_master)) {
  "Max |mean normalized log2(H/L)| — Track A"
} else {
  "Max |mean log2(Drug/DMSO)| — Track B raw-LFQ"
}

if (!is.null(y_col) && !is.na(y_col) && y_col %in% names(targets_master) && x_col %in% names(targets_master)) {
  scatter_dt <- targets_master[detected & (is.finite(get(x_col)) | is.finite(get(y_col)))]
  
  scatter_dt[, show_label := (
    target_class %in% c("strong stabilization", "destabilization") |
      (is.finite(p_value_best) & p_value_best < p_threshold) |
      evidence_conflict == TRUE
  )]
  
  p_scatter <- ggplot(
    scatter_dt,
    aes(x = get(x_col), y = get(y_col), colour = target_class, label = target_input)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_point(size = 2.2, alpha = 0.85) +
    geom_text(data = scatter_dt[show_label == TRUE], size = 2.4, vjust = -0.5, check_overlap = TRUE) +
    scale_colour_manual(values = class_colours, name = "CETSA class") +
    theme_bw(base_size = 11) +
    labs(
      title = "Target classification summary",
      subtitle = "Track A SILAC ratio is prioritized; Track B raw-LFQ deltaTm is supportive/fallback",
      x = x_label,
      y = y_label
    )
  
  ggsave(file.path(output_dir, "target_classification_scatter.png"), p_scatter,
         width = 11, height = 8, dpi = plot_dpi, bg = "white")
}

#############################################################################################################
# HEATMAP OF DETECTED TARGETS
#############################################################################################################

heat_temp_dt <- if (nrow(ratio_temp_dt) > 0 && "mean_log2_ratio" %in% names(ratio_temp_dt)) {
  msg("Heatmap using Track A ratio data.")
  ratio_temp_dt
} else {
  msg("Heatmap using Track B raw-LFQ temp_effect_summary.")
  temp_dt
}

detected_temp <- merge(
  targets_master[detected == TRUE, .(target_input, target_class, Protein_group, Genes)],
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
  
  if (nrow(heat_wide) > 1 && ncol(heat_wide) > 4) {
    heat_mat <- as.matrix(heat_wide[, -(1:4)])
    mode(heat_mat) <- "numeric"
    rownames(heat_mat) <- paste0(heat_wide$target_input, " | ", heat_wide$target_class)
    
    keep_rows <- rowSums(is.finite(heat_mat)) >= 2
    heat_mat <- heat_mat[keep_rows, , drop = FALSE]
    heat_wide_kept <- heat_wide[keep_rows]
    
    if (nrow(heat_mat) > 1) {
      row_score <- apply(heat_mat, 1, function(x) {
        x <- x[is.finite(x)]
        if (!length(x)) return(-Inf)
        max(abs(x), na.rm = TRUE)
      })
      ord <- order(row_score, decreasing = TRUE)
      heat_mat <- heat_mat[ord, , drop = FALSE]
      heat_wide_kept <- heat_wide_kept[ord]
      
      heat_mat[!is.finite(heat_mat)] <- 0
      row_var <- apply(heat_mat, 1, var, na.rm = TRUE)
      keep_var <- is.finite(row_var) & row_var > 0
      heat_mat <- heat_mat[keep_var, , drop = FALSE]
      heat_wide_kept <- heat_wide_kept[keep_var]
      
      if (nrow(heat_mat) > 1) {
        ann_row <- data.frame(
          Class = as.character(heat_wide_kept$target_class),
          row.names = rownames(heat_mat)
        )
        ann_colours <- list(Class = class_colours[names(class_colours) %in% ann_row$Class])
        
        png(file.path(output_dir, "target_classes_heatmap.png"),
            width = 1800, height = max(1200, 26 * nrow(heat_mat)), res = 200)
        pheatmap(
          heat_mat,
          annotation_row = ann_row,
          annotation_colors = ann_colours,
          scale = "none",
          cluster_rows = TRUE,
          cluster_cols = FALSE,
          border_color = NA,
          main = paste0("Detected targets: mean log2(Drug/DMSO) across temperatures — ", heat_source_label),
          fontsize_row = 7,
          fontsize_col = 10,
          na_col = "grey90"
        )
        dev.off()
      } else {
        msg("Heatmap skipped: fewer than two variable detected targets after cleaning.")
      }
    } else {
      msg("Heatmap skipped: fewer than two detected targets with enough finite temperatures.")
    }
  } else {
    msg("Heatmap skipped: insufficient detected target rows or temperature columns.")
  }
}

#############################################################################################################
# CONSOLE SUMMARY
#############################################################################################################

msg("\n══════════════════════════════════════════════════════")
msg("Done. Output folder: ", output_dir)
msg("══════════════════════════════════════════════════════")
msg("Requested targets : ", uniqueN(targets))
msg("Detected          : ", targets_master[detected == TRUE, .N])
msg("Not detected      : ", targets_master[detected == FALSE, .N])

msg("\nClass counts:")
print(targets_master[, .N, by = target_class][order(target_class)])

show_cols <- intersect(
  c(
    "target_input", "Genes", "Protein_group", "primary_evidence", "stabilization_direction",
    "target_class", "evidence_conflict", "signed_peak_log2_HL", "max_abs_mean_log2_HL",
    "deltaTm_mean", "deltaTm_mean_curve", "deltaAUC_mean", score_col_B,
    "p_value_best", "p_value_source"
  ),
  names(targets_master)
)

msg("\nStrong stabilization targets:")
print(targets_master[target_class == "strong stabilization", ..show_cols][seq_len(min(30L, .N))])

msg("\nDestabilization targets:")
print(targets_master[target_class == "destabilization", ..show_cols][seq_len(min(20L, .N))])

msg("\nEvidence conflicts, Track A vs Track B raw-LFQ:")
print(targets_master[evidence_conflict == TRUE, ..show_cols][seq_len(min(20L, .N))])
