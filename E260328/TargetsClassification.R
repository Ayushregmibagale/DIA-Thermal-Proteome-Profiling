#############################################################################################################
#
# CETSA-MS SILAC TARGET CLASSIFICATION SCRIPT
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
# INPUTS EXPECTED
# ---------------
# From the main CETSA SILAC pipeline output directory:
#   - CETSA_SILAC_hit_table.tsv.gz
#   - temp_effect_summary.tsv.gz
#
# OPTIONAL
# --------
#   - matched_targets_primary.tsv from the first target-inspection script
#     If not found, this script will rematch targets itself.
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

input_dir  <- "CETSA_SILAC_output"
output_dir <- file.path(input_dir, "Target_classification")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Classification thresholds
# -------------------------
# Tune these if needed after first pass.
strong_deltaTm_threshold <- 2.0
weak_deltaTm_threshold   <- 1.0

strong_log2ratio_threshold <- 0.35
weak_log2ratio_threshold   <- 0.20

strong_deltaAUC_threshold <- 1.5
weak_deltaAUC_threshold   <- 0.5

# minimum consistency requirement
# if there is clearly negative evidence beyond this, classify as destabilization
destab_deltaTm_threshold     <- -1.0
destab_log2ratio_threshold   <- -0.20
destab_deltaAUC_threshold    <- -0.5

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
  x <- gsub("\\s+", "", x)
  x
}

extract_base_symbol <- function(x) {
  x <- clean_target_name(x)
  x <- sub("\\(.*\\)", "", x)
  x <- sub("-PHOSPHORYLATED$", "", x)
  x <- sub("-NONPHOSPHORYLATED$", "", x)
  x <- sub("-KIN\\.DOM\\..*$", "", x)
  x <- sub("-CYCLIND1$", "", x)
  x <- sub("-CYCLIND3$", "", x)
  x <- sub("-CATALYTIC$", "", x)
  x <- sub("-PSEUDOKINASE$", "", x)
  x <- sub("-$", "", x)
  x
}

map_target_alias <- function(x) {
  y <- extract_base_symbol(x)
  
  alias_map <- c(
    "VEGFR2" = "KDR",
    "VEGFR1" = "FLT1",
    "VEGFR3" = "FLT4",
    "TIE2"   = "TEK",
    "TRKA"   = "NTRK1",
    "TRKB"   = "NTRK2",
    "TRKC"   = "NTRK3",
    "PYK2"   = "PTK2B",
    "FAK"    = "PTK2",
    "BRK"    = "PTK6",
    "FRK"    = "PTK5",
    "YES"    = "YES1",
    "LZK"    = "MAP3K13",
    "DLK"    = "MAP3K12",
    "HPK1"   = "MAP4K1",
    "MEK1"   = "MAP2K1",
    "MEK2"   = "MAP2K2",
    "MEK3"   = "MAP2K3",
    "MEK4"   = "MAP2K4",
    "MEK5"   = "MAP2K5",
    "MEK6"   = "MAP2K6",
    "ERK1"   = "MAPK3",
    "ERK2"   = "MAPK1",
    "ERK3"   = "MAPK6",
    "ERK4"   = "MAPK4",
    "ERK5"   = "MAPK7",
    "ERK8"   = "MAPK15",
    "JNK1"   = "MAPK8",
    "JNK2"   = "MAPK9",
    "JNK3"   = "MAPK10",
    "P38-ALPHA" = "MAPK14",
    "P38-BETA"  = "MAPK11",
    "P38-GAMMA" = "MAPK12",
    "P38-DELTA" = "MAPK13",
    "IKK-ALPHA" = "CHUK",
    "IKK-BETA"  = "IKBKB",
    "IKK-EPSILON" = "IKBKE",
    "PKAC-ALPHA" = "PRKACA",
    "PKAC-BETA"  = "PRKACB",
    "PFTK1" = "CDK14",
    "PFTAIRE2" = "CDK15",
    "PCTK1" = "CDK16",
    "PCTK2" = "CDK17",
    "PCTK3" = "CDK18",
    "CDC2L1" = "CDK11B",
    "CDC2L2" = "CDK11A",
    "CDC2L5" = "CDK20",
    "MRCKA" = "CDC42BPA",
    "MRCKB" = "CDC42BPB",
    "DMPK2" = "CDC42BPG",
    "TRKA" = "NTRK1",
    "TRKB" = "NTRK2",
    "TRKC" = "NTRK3"
  )
  
  if (y %in% names(alias_map)) alias_map[[y]] else y
}

split_gene_tokens <- function(x) {
  if (is.na(x) || x == "") return(character())
  parts <- unlist(strsplit(toupper(x), ";|,|\\|"))
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

first_non_na <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

#############################################################################################################
# READ DATA
#############################################################################################################

hit_file <- file.path(input_dir, "CETSA_SILAC_hit_table.tsv.gz")
temp_file <- file.path(input_dir, "temp_effect_summary.tsv.gz")
matched_file <- file.path(input_dir, "Target_inspection", "matched_targets_primary.tsv")

if (!file.exists(hit_file)) stop("Missing file: ", hit_file)
if (!file.exists(temp_file)) stop("Missing file: ", temp_file)

hit_dt  <- fread(hit_file)
temp_dt <- fread(temp_file)

#############################################################################################################
# BUILD / LOAD MATCH TABLE
#############################################################################################################

if (file.exists(matched_file)) {
  msg("Using existing match table: ", matched_file)
  matched_primary <- fread(matched_file)
} else {
  msg("No previous match table found. Rebuilding matches...")
  
  targets_dt <- data.table(target_input = unique(targets))
  targets_dt[, target_clean := clean_target_name(target_input)]
  targets_dt[, target_base  := vapply(target_input, extract_base_symbol, character(1))]
  targets_dt[, target_gene  := vapply(target_input, map_target_alias, character(1))]
  
  hit_dt[, Genes := as.character(Genes)]
  hit_dt[, Protein_group := as.character(Protein_group)]
  
  expanded_hits <- hit_dt[, {
    toks <- split_gene_tokens(Genes)
    if (length(toks) == 0) toks <- NA_character_
    .(gene_token = toks)
  }, by = .(Protein_group, Genes)]
  
  expanded_hits[, gene_token := toupper(gene_token)]
  
  match_gene <- merge(
    targets_dt,
    expanded_hits,
    by.x = "target_gene",
    by.y = "gene_token",
    all.x = TRUE,
    allow.cartesian = TRUE
  )
  match_gene[, match_type := fifelse(!is.na(Protein_group), "gene_match", NA_character_)]
  
  protein_groups_upper <- unique(hit_dt[, .(Protein_group, Protein_group_upper = toupper(Protein_group))])
  
  match_pg <- merge(
    targets_dt,
    protein_groups_upper,
    by.x = "target_clean",
    by.y = "Protein_group_upper",
    all.x = TRUE,
    allow.cartesian = TRUE
  )
  match_pg <- merge(match_pg, hit_dt, by = "Protein_group", all.x = TRUE)
  match_pg[, match_type := fifelse(!is.na(Genes), "protein_group_match", NA_character_)]
  
  matched_all <- rbindlist(list(
    match_gene[, .(target_input, target_clean, target_base, target_gene, Protein_group, Genes, match_type)],
    match_pg[,   .(target_input, target_clean, target_base, target_gene, Protein_group, Genes, match_type)]
  ), fill = TRUE)
  
  matched_all <- unique(matched_all[!is.na(Protein_group)])
  match_priority <- c(gene_match = 1L, protein_group_match = 2L)
  matched_all[, match_priority := match_priority[match_type]]
  setorder(matched_all, target_input, match_priority, Protein_group)
  matched_primary <- matched_all[, .SD[1], by = target_input]
}

#############################################################################################################
# PREPARE MASTER TARGET TABLE
#############################################################################################################

targets_master <- data.table(target_input = unique(targets))
targets_master <- merge(targets_master, matched_primary, by = "target_input", all.x = TRUE)
targets_master <- merge(targets_master, hit_dt, by = c("Protein_group", "Genes"), all.x = TRUE)

#############################################################################################################
# ADD TEMPERATURE-SHAPE INFORMATION
#############################################################################################################

temp_target <- merge(
  matched_primary[, .(target_input, Protein_group, Genes)],
  temp_dt,
  by = c("Protein_group", "Genes"),
  all.x = TRUE
)

shape_summary <- temp_target[, .(
  n_temp_points = sum(is.finite(mean_log2_ratio)),
  max_ratio = if (all(!is.finite(mean_log2_ratio))) NA_real_ else max(mean_log2_ratio, na.rm = TRUE),
  min_ratio = if (all(!is.finite(mean_log2_ratio))) NA_real_ else min(mean_log2_ratio, na.rm = TRUE),
  max_abs_ratio = if (all(!is.finite(mean_log2_ratio))) NA_real_ else max(abs(mean_log2_ratio), na.rm = TRUE),
  temp_at_max = if (all(!is.finite(mean_log2_ratio))) NA_real_ else Temp[which.max(mean_log2_ratio)],
  temp_at_min = if (all(!is.finite(mean_log2_ratio))) NA_real_ else Temp[which.min(mean_log2_ratio)]
), by = .(target_input, Protein_group, Genes)]

targets_master <- merge(
  targets_master,
  shape_summary,
  by = c("target_input", "Protein_group", "Genes"),
  all.x = TRUE
)

#############################################################################################################
# CLASSIFICATION LOGIC
#############################################################################################################

targets_master[, detected := !is.na(Protein_group)]

# robust positive evidence
targets_master[, strong_positive_evidence :=
                 (
                   (is.finite(deltaTm_mean)       & deltaTm_mean       >= strong_deltaTm_threshold) |
                     (is.finite(deltaTm_mean_curve) & deltaTm_mean_curve >= strong_deltaTm_threshold) |
                     (is.finite(max_abs_mean_log2_ratio) & max_ratio >= strong_log2ratio_threshold) |
                     (is.finite(deltaAUC_mean)      & deltaAUC_mean      >= strong_deltaAUC_threshold)
                 )
]

targets_master[, weak_positive_evidence :=
                 (
                   (is.finite(deltaTm_mean)       & deltaTm_mean       >= weak_deltaTm_threshold) |
                     (is.finite(deltaTm_mean_curve) & deltaTm_mean_curve >= weak_deltaTm_threshold) |
                     (is.finite(max_abs_mean_log2_ratio) & max_ratio >= weak_log2ratio_threshold) |
                     (is.finite(deltaAUC_mean)      & deltaAUC_mean      >= weak_deltaAUC_threshold)
                 )
]

targets_master[, negative_evidence :=
                 (
                   (is.finite(deltaTm_mean)       & deltaTm_mean       <= destab_deltaTm_threshold) |
                     (is.finite(deltaTm_mean_curve) & deltaTm_mean_curve <= destab_deltaTm_threshold) |
                     (is.finite(min_ratio)          & min_ratio          <= destab_log2ratio_threshold) |
                     (is.finite(deltaAUC_mean)      & deltaAUC_mean      <= destab_deltaAUC_threshold)
                 )
]

targets_master[, target_class := fifelse(
  !detected, "not detected",
  fifelse(
    strong_positive_evidence & !negative_evidence, "strong stabilization",
    fifelse(
      negative_evidence & !strong_positive_evidence, "destabilization",
      fifelse(
        weak_positive_evidence & !negative_evidence, "weak stabilization",
        "detected but flat/no effect"
      )
    )
  )
)]

# ordered factor for output and plots
targets_master[, target_class := factor(
  target_class,
  levels = c(
    "strong stabilization",
    "weak stabilization",
    "destabilization",
    "detected but flat/no effect",
    "not detected"
  )
)]

#############################################################################################################
# ADD HUMAN-READABLE REASON
#############################################################################################################

targets_master[, classification_reason := ""]
targets_master[!detected, classification_reason := "No CETSA match found in hit table"]

targets_master[target_class == "strong stabilization", classification_reason :=
                 paste0(
                   "Strong positive CETSA evidence; ",
                   "deltaTm_mean=", ifelse(is.finite(deltaTm_mean), round(deltaTm_mean, 2), "NA"),
                   ", deltaTm_mean_curve=", ifelse(is.finite(deltaTm_mean_curve), round(deltaTm_mean_curve, 2), "NA"),
                   ", max_ratio=", ifelse(is.finite(max_ratio), round(max_ratio, 3), "NA"),
                   ", deltaAUC_mean=", ifelse(is.finite(deltaAUC_mean), round(deltaAUC_mean, 3), "NA")
                 )
]

targets_master[target_class == "weak stabilization", classification_reason :=
                 paste0(
                   "Moderate positive CETSA evidence; ",
                   "deltaTm_mean=", ifelse(is.finite(deltaTm_mean), round(deltaTm_mean, 2), "NA"),
                   ", deltaTm_mean_curve=", ifelse(is.finite(deltaTm_mean_curve), round(deltaTm_mean_curve, 2), "NA"),
                   ", max_ratio=", ifelse(is.finite(max_ratio), round(max_ratio, 3), "NA"),
                   ", deltaAUC_mean=", ifelse(is.finite(deltaAUC_mean), round(deltaAUC_mean, 3), "NA")
                 )
]

targets_master[target_class == "destabilization", classification_reason :=
                 paste0(
                   "Negative CETSA shift; ",
                   "deltaTm_mean=", ifelse(is.finite(deltaTm_mean), round(deltaTm_mean, 2), "NA"),
                   ", deltaTm_mean_curve=", ifelse(is.finite(deltaTm_mean_curve), round(deltaTm_mean_curve, 2), "NA"),
                   ", min_ratio=", ifelse(is.finite(min_ratio), round(min_ratio, 3), "NA"),
                   ", deltaAUC_mean=", ifelse(is.finite(deltaAUC_mean), round(deltaAUC_mean, 3), "NA")
                 )
]

targets_master[target_class == "detected but flat/no effect", classification_reason :=
                 paste0(
                   "Detected but no strong CETSA shift; ",
                   "deltaTm_mean=", ifelse(is.finite(deltaTm_mean), round(deltaTm_mean, 2), "NA"),
                   ", deltaTm_mean_curve=", ifelse(is.finite(deltaTm_mean_curve), round(deltaTm_mean_curve, 2), "NA"),
                   ", max_abs_ratio=", ifelse(is.finite(max_abs_ratio), round(max_abs_ratio, 3), "NA"),
                   ", deltaAUC_mean=", ifelse(is.finite(deltaAUC_mean), round(deltaAUC_mean, 3), "NA")
                 )
]

#############################################################################################################
# ORDER OUTPUT TABLE
#############################################################################################################

targets_master[, class_priority := fifelse(
  target_class == "strong stabilization", 1,
  fifelse(
    target_class == "weak stabilization", 2,
    fifelse(
      target_class == "destabilization", 3,
      fifelse(
        target_class == "detected but flat/no effect", 4, 5
      )
    )
  )
)]

targets_master[, evidence_score := 0]
targets_master[is.finite(deltaTm_mean),       evidence_score := evidence_score + abs(deltaTm_mean)]
targets_master[is.finite(deltaTm_mean_curve), evidence_score := evidence_score + abs(deltaTm_mean_curve)]
targets_master[is.finite(max_abs_mean_log2_ratio), evidence_score := evidence_score + abs(max_abs_mean_log2_ratio)]
targets_master[is.finite(deltaAUC_mean), evidence_score := evidence_score + abs(deltaAUC_mean)]

setorder(targets_master, class_priority, -evidence_score, target_input)

#############################################################################################################
# WRITE OUTPUT TABLES
#############################################################################################################

fwrite(targets_master, file.path(output_dir, "target_classification_full.tsv"), sep = "\t")

fwrite(
  targets_master[, .(
    target_input,
    target_class,
    Genes,
    Protein_group,
    match_type,
    deltaTm_mean,
    deltaTm_mean_curve,
    deltaAUC_mean,
    max_abs_mean_log2_ratio,
    classification_reason
  )],
  file.path(output_dir, "target_classification_summary.tsv"),
  sep = "\t"
)

fwrite(
  targets_master[target_class == "strong stabilization"],
  file.path(output_dir, "targets_strong_stabilization.tsv"),
  sep = "\t"
)

fwrite(
  targets_master[target_class == "weak stabilization"],
  file.path(output_dir, "targets_weak_stabilization.tsv"),
  sep = "\t"
)

fwrite(
  targets_master[target_class == "destabilization"],
  file.path(output_dir, "targets_destabilization.tsv"),
  sep = "\t"
)

fwrite(
  targets_master[target_class == "detected but flat/no effect"],
  file.path(output_dir, "targets_flat_no_effect.tsv"),
  sep = "\t"
)

fwrite(
  targets_master[target_class == "not detected"],
  file.path(output_dir, "targets_not_detected.tsv"),
  sep = "\t"
)

#############################################################################################################
# CLASS COUNT PLOT
#############################################################################################################

class_counts <- targets_master[, .N, by = target_class]

p_counts <- ggplot(class_counts, aes(x = target_class, y = N, fill = target_class)) +
  geom_col() +
  theme_bw(base_size = 11) +
  labs(
    title = "Requested targets by CETSA class",
    x = NULL,
    y = "Count"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")

ggsave(
  file.path(output_dir, "target_class_counts.png"),
  p_counts,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)

#############################################################################################################
# SUMMARY SCATTER
#############################################################################################################

p_scatter <- ggplot(
  targets_master[detected == TRUE],
  aes(x = deltaTm_mean, y = max_abs_mean_log2_ratio, color = target_class, label = target_input)
) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_text(size = 2.4, vjust = -0.4, check_overlap = TRUE) +
  theme_bw(base_size = 11) +
  labs(
    title = "Target classification summary",
    x = "Mean deltaTm (Drug - DMSO)",
    y = "Max absolute mean log2(Drug/DMSO)"
  )

ggsave(
  file.path(output_dir, "target_classification_scatter.png"),
  p_scatter,
  width = 11,
  height = 8,
  dpi = 300,
  bg = "white"
)

#############################################################################################################
# OPTIONAL HEATMAP OF DETECTED TARGETS
#############################################################################################################

detected_temp <- merge(
  targets_master[detected == TRUE, .(target_input, target_class, Protein_group, Genes)],
  temp_dt,
  by = c("Protein_group", "Genes"),
  all.x = TRUE
)

if (nrow(detected_temp) > 0) {
  heat_dt <- dcast(
    detected_temp,
    target_input + target_class + Genes + Protein_group ~ Temp,
    value.var = "mean_log2_ratio"
  )
  
  if (nrow(heat_dt) > 0) {
    heat_mat <- as.matrix(heat_dt[, -(1:4)])
    rownames(heat_mat) <- paste0(heat_dt$target_input, " | ", heat_dt$target_class)
    
    png(
      file.path(output_dir, "target_classes_heatmap.png"),
      width = 1800,
      height = max(1200, 24 * nrow(heat_mat)),
      res = 200
    )
    
    pheatmap(
      heat_mat,
      scale = "none",
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      border_color = NA,
      main = "Detected targets: mean log2(Drug/DMSO) across temperatures",
      fontsize_row = 7,
      fontsize_col = 10
    )
    
    dev.off()
  }
}

#############################################################################################################
# CONSOLE SUMMARY
#############################################################################################################

msg("Done.")
msg("Output folder: ", output_dir)
msg("Requested targets: ", nrow(unique(data.table(target_input = targets))))
msg("Detected targets: ", targets_master[detected == TRUE, .N])
msg("Not detected: ", targets_master[detected == FALSE, .N])

msg("\nClass counts:")
print(targets_master[, .N, by = target_class][order(target_class)])

msg("\nTop strong stabilization targets:")
print(
  targets_master[target_class == "strong stabilization", .(
    target_input,
    Genes,
    Protein_group,
    deltaTm_mean = round(deltaTm_mean, 2),
    deltaTm_mean_curve = round(deltaTm_mean_curve, 2),
    deltaAUC_mean = round(deltaAUC_mean, 3),
    max_abs_mean_log2_ratio = round(max_abs_mean_log2_ratio, 3)
  )][1:min(30, .N)]
)