#############################################################################################################
#
# CETSA-MS SILAC TARGET INSPECTION SCRIPT
#
# PURPOSE
# -------
# Inspect a predefined kinase/target list in your CETSA-MS SILAC results.
#
# INPUT FILES EXPECTED
# --------------------
# From the main CETSA SILAC pipeline output directory:
#   - long_fraction_table.tsv.gz
#   - mean_curves.tsv.gz
#   - fit_predictions_mean_curves.tsv.gz   (optional but recommended)
#   - CETSA_SILAC_hit_table.tsv.gz
#   - temp_effect_summary.tsv.gz
#
# OUTPUT
# ------
#   - matched_targets.tsv
#   - unmatched_targets.tsv
#   - target_summary.tsv
#   - individual target curve PNGs
#   - one combined PDF of target curves
#   - heatmap of mean log2 Drug/DMSO across temperature for matched targets
#
#############################################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(matrixStats)
  library(gridExtra)
})

#############################################################################################################
# USER SETTINGS
#############################################################################################################

input_dir  <- "CETSA_SILAC_output"
output_dir <- file.path(input_dir, "Target_inspection")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

plot_width  <- 10
plot_height <- 7
plot_dpi    <- 300
max_fraction_cap <- 1.5

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

first_non_na <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x)
}

safe_sem <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x) / sqrt(length(x))
}

clean_target_name <- function(x) {
  x <- toupper(trimws(x))
  x <- gsub("\\s+", "", x)
  x
}

extract_base_symbol <- function(x) {
  x <- clean_target_name(x)
  
  # remove anything in parentheses
  x <- sub("\\(.*\\)", "", x)
  
  # remove common suffix annotations
  x <- sub("-PHOSPHORYLATED$", "", x)
  x <- sub("-NONPHOSPHORYLATED$", "", x)
  x <- sub("-KIN\\.DOM\\..*$", "", x)
  x <- sub("-CYCLIND1$", "", x)
  x <- sub("-CYCLIND3$", "", x)
  x <- sub("-CATALYTIC$", "", x)
  x <- sub("-PSEUDOKINASE$", "", x)
  
  # remove trailing punctuation fragments
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
    "MUSK"   = "MUSK",
    "PYK2"   = "PTK2B",
    "FAK"    = "PTK2",
    "BRK"    = "PTK6",
    "FRK"    = "PTK5",
    "BLK"    = "BLK",
    "HCK"    = "HCK",
    "FGR"    = "FGR",
    "YES"    = "YES1",
    "LZK"    = "MAP3K13",
    "DLK"    = "MAP3K12",
    "HPK1"   = "MAP4K1",
    "MST1R"  = "MST1R",
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
    "p38-alpha" = "MAPK14",
    "p38-beta"  = "MAPK11",
    "p38-gamma" = "MAPK12",
    "p38-delta" = "MAPK13",
    "IKK-alpha" = "CHUK",
    "IKK-beta"  = "IKBKB",
    "IKK-epsilon" = "IKBKE",
    "PKAC-alpha" = "PRKACA",
    "PKAC-beta"  = "PRKACB",
    "PFTK1" = "CDK14",
    "PFTAIRE2" = "CDK15",
    "PCTK1" = "CDK16",
    "PCTK2" = "CDK17",
    "PCTK3" = "CDK18",
    "CDC2L1" = "CDK11B",
    "CDC2L2" = "CDK11A",
    "CDC2L5" = "CDK20",
    "DMPK2"  = "CDC42BPG",
    "MRCKA"  = "CDC42BPA",
    "MRCKB"  = "CDC42BPB",
    "MELK"   = "MELK",
    "AURKA"  = "AURKA",
    "AURKB"  = "AURKB",
    "AURKC"  = "AURKC",
    "BMX"    = "BMX",
    "BTK"    = "BTK",
    "ITK"    = "ITK",
    "TXK"    = "TXK",
    "TEC"    = "TEC",
    "MERTK"  = "MERTK",
    "TYRO3"  = "TYRO3",
    "CSF1R"  = "CSF1R",
    "RET"    = "RET",
    "ALK"    = "ALK",
    "ROS1"   = "ROS1",
    "MET"    = "MET",
    "EGFR"   = "EGFR",
    "ERBB2"  = "ERBB2",
    "ERBB3"  = "ERBB3",
    "ERBB4"  = "ERBB4",
    "FGFR1"  = "FGFR1",
    "FGFR2"  = "FGFR2",
    "FGFR3"  = "FGFR3",
    "FGFR4"  = "FGFR4",
    "FLT3"   = "FLT3",
    "FLT1"   = "FLT1",
    "FLT4"   = "FLT4",
    "KIT"    = "KIT",
    "PDGFRA" = "PDGFRA",
    "PDGFRB" = "PDGFRB",
    "TIE1"   = "TIE1",
    "SRC"    = "SRC",
    "FYN"    = "FYN",
    "LYN"    = "LYN",
    "LCK"    = "LCK",
    "ABL1"   = "ABL1",
    "ABL2"   = "ABL2",
    "SYK"    = "SYK",
    "ZAP70"  = "ZAP70",
    "FAK"    = "PTK2",
    "PYK2"   = "PTK2B",
    "TRPM6"  = "TRPM6",
    "MTOR"   = "MTOR",
    "INSR"   = "INSR",
    "INSRR"  = "INSRR",
    "IGF1R"  = "IGF1R"
  )
  
  if (y %in% names(alias_map)) alias_map[[y]] else y
}

split_gene_tokens <- function(x) {
  if (is.na(x) || x == "") return(character())
  parts <- unlist(strsplit(toupper(x), ";|,|\\|"))
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

#############################################################################################################
# READ INPUTS
#############################################################################################################

msg("Reading CETSA result files...")

long_file      <- file.path(input_dir, "long_fraction_table.tsv.gz")
mean_file      <- file.path(input_dir, "mean_curves.tsv.gz")
fitpred_file   <- file.path(input_dir, "fit_predictions_mean_curves.tsv.gz")
hit_file       <- file.path(input_dir, "CETSA_SILAC_hit_table.tsv.gz")
temp_eff_file  <- file.path(input_dir, "temp_effect_summary.tsv.gz")

if (!file.exists(long_file))     stop("Missing file: ", long_file)
if (!file.exists(mean_file))     stop("Missing file: ", mean_file)
if (!file.exists(hit_file))      stop("Missing file: ", hit_file)
if (!file.exists(temp_eff_file)) stop("Missing file: ", temp_eff_file)

long_dt     <- fread(long_file)
mean_dt     <- fread(mean_file)
hit_dt      <- fread(hit_file)
temp_eff_dt <- fread(temp_eff_file)
fitpred_dt  <- if (file.exists(fitpred_file)) fread(fitpred_file) else data.table()

#############################################################################################################
# BUILD TARGET TABLE
#############################################################################################################

targets_dt <- data.table(
  target_input = unique(targets)
)

targets_dt[, target_clean := clean_target_name(target_input)]
targets_dt[, target_base  := vapply(target_input, extract_base_symbol, character(1))]
targets_dt[, target_gene  := vapply(target_input, map_target_alias, character(1))]

#############################################################################################################
# PREPARE CETSA TABLES FOR MATCHING
#############################################################################################################

hit_dt[, Genes := as.character(Genes)]
hit_dt[, Protein_group := as.character(Protein_group)]

# expand hit table by gene tokens so "MAPK1;MAPK3" etc can still match
expanded_hits <- hit_dt[, {
  toks <- split_gene_tokens(Genes)
  if (length(toks) == 0) toks <- NA_character_
  .(gene_token = toks)
}, by = .(Protein_group, Genes)]

expanded_hits[, gene_token := toupper(gene_token)]

#############################################################################################################
# MATCH TARGETS
#############################################################################################################

# primary match by gene symbol
match_gene <- merge(
  targets_dt,
  expanded_hits,
  by.x = "target_gene",
  by.y = "gene_token",
  all.x = TRUE,
  allow.cartesian = TRUE
)

match_gene[, match_type := fifelse(!is.na(Protein_group), "gene_match", NA_character_)]

# secondary fallback match by protein_group text
protein_groups_upper <- unique(hit_dt[, .(Protein_group, Protein_group_upper = toupper(Protein_group))])

match_pg <- merge(
  targets_dt,
  protein_groups_upper,
  by.x = "target_clean",
  by.y = "Protein_group_upper",
  all.x = TRUE,
  allow.cartesian = TRUE
)

match_pg <- merge(
  match_pg,
  hit_dt,
  by = "Protein_group",
  all.x = TRUE
)

match_pg[, match_type := fifelse(!is.na(Genes), "protein_group_match", NA_character_)]

# combine matches
matched_all <- rbindlist(list(
  match_gene[, .(target_input, target_clean, target_base, target_gene, Protein_group, Genes, match_type)],
  match_pg[,   .(target_input, target_clean, target_base, target_gene, Protein_group, Genes, match_type)]
), fill = TRUE)

matched_all <- unique(matched_all[!is.na(Protein_group)])

# choose one preferred match per target if duplicates exist
match_priority <- c(gene_match = 1L, protein_group_match = 2L)
matched_all[, match_priority := match_priority[match_type]]
setorder(matched_all, target_input, match_priority, Protein_group)

matched_primary <- matched_all[, .SD[1], by = target_input]

unmatched_targets <- targets_dt[!target_input %in% matched_primary$target_input]

#############################################################################################################
# MERGE TARGETS WITH HIT SUMMARY
#############################################################################################################

target_summary <- merge(
  matched_primary,
  hit_dt,
  by = c("Protein_group", "Genes"),
  all.x = TRUE
)

# order targets by strongest CETSA evidence
target_summary[, rank_abs := fifelse(is.finite(rank_score), abs(rank_score), -Inf)]
target_summary[, deltaTm_abs := fifelse(is.finite(deltaTm_mean), abs(deltaTm_mean), -Inf)]
target_summary[, ratio_abs := fifelse(is.finite(max_abs_mean_log2_ratio), abs(max_abs_mean_log2_ratio), -Inf)]
setorder(target_summary, -rank_abs, -deltaTm_abs, -ratio_abs)

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
  dt_raw  <- long_dt[Protein_group == protein_group_id]
  dt_mean <- mean_dt[Protein_group == protein_group_id]
  dt_fit  <- fitpred_dt[Protein_group == protein_group_id]
  
  if (nrow(dt_raw) == 0) return(NULL)
  
  hit_row <- target_summary[Protein_group == protein_group_id][1]
  if (is.null(target_label) || is.na(target_label)) {
    target_label <- first_non_na(hit_row$target_input)
  }
  
  title_gene <- first_non_na(dt_raw$Genes)
  if (is.na(title_gene)) title_gene <- protein_group_id
  
  subtxt <- paste0(
    "Requested target: ", target_label,
    " | matched gene: ", first_non_na(hit_row$Genes),
    " | match type: ", first_non_na(hit_row$match_type),
    "\n",
    "deltaTm_mean=",
    ifelse(is.finite(hit_row$deltaTm_mean), round(hit_row$deltaTm_mean, 2), "NA"),
    " | deltaTm_mean_curve=",
    ifelse(is.finite(hit_row$deltaTm_mean_curve), round(hit_row$deltaTm_mean_curve, 2), "NA"),
    " | deltaAUC_mean=",
    ifelse(is.finite(hit_row$deltaAUC_mean), round(hit_row$deltaAUC_mean, 3), "NA"),
    " | max_abs_mean_log2_ratio=",
    ifelse(is.finite(hit_row$max_abs_mean_log2_ratio), round(hit_row$max_abs_mean_log2_ratio, 3), "NA")
  )
  
  p <- ggplot() +
    geom_point(
      data = dt_raw[is.finite(frac)],
      aes(x = Temp, y = frac, color = condition, shape = Replicate),
      size = 2, alpha = 0.75
    ) +
    geom_line(
      data = dt_raw[is.finite(frac)],
      aes(x = Temp, y = frac, color = condition, group = interaction(condition, Replicate)),
      alpha = 0.4
    ) +
    geom_errorbar(
      data = dt_mean[is.finite(mean_frac)],
      aes(x = Temp, ymin = mean_frac - sem_frac, ymax = mean_frac + sem_frac, color = condition),
      width = 0.25, linewidth = 0.5
    ) +
    geom_point(
      data = dt_mean[is.finite(mean_frac)],
      aes(x = Temp, y = mean_frac, color = condition),
      size = 3
    ) +
    geom_line(
      data = dt_mean[is.finite(mean_frac)],
      aes(x = Temp, y = mean_frac, color = condition),
      linewidth = 1
    ) +
    theme_bw(base_size = 11) +
    labs(
      title = title_gene,
      subtitle = subtxt,
      x = "Temperature (°C)",
      y = "Normalized soluble fraction"
    ) +
    ylim(0, max_fraction_cap)
  
  if (nrow(dt_fit) > 0) {
    p <- p +
      geom_line(
        data = dt_fit[is.finite(frac_pred)],
        aes(x = Temp, y = frac_pred, color = condition),
        linewidth = 1.1,
        linetype = 2
      )
  }
  
  p
}

#############################################################################################################
# INDIVIDUAL TARGET CURVE PLOTS
#############################################################################################################

msg("Making target curve plots...")

curve_dir <- file.path(output_dir, "Target_curves")
dir.create(curve_dir, showWarnings = FALSE)

plot_list <- list()

if (nrow(target_summary) > 0) {
  for (i in seq_len(nrow(target_summary))) {
    this_pg    <- target_summary$Protein_group[i]
    this_label <- target_summary$target_input[i]
    
    p <- plot_one_target(this_pg, this_label)
    if (!is.null(p)) {
      plot_list[[length(plot_list) + 1]] <- p
      
      out_name <- paste0(
        sprintf("%03d", i), "_",
        gsub("[^A-Za-z0-9_\\-]", "_", this_label), "_",
        gsub("[^A-Za-z0-9_\\-]", "_", this_pg),
        ".png"
      )
      
      ggsave(
        file.path(curve_dir, out_name),
        p,
        width = plot_width,
        height = plot_height,
        dpi = plot_dpi,
        bg = "white"
      )
    }
  }
}

#############################################################################################################
# COMBINED PDF OF TARGET CURVES
#############################################################################################################

if (length(plot_list) > 0) {
  pdf(file.path(output_dir, "All_target_curves.pdf"), width = 10, height = 7)
  for (p in plot_list) print(p)
  dev.off()
}

#############################################################################################################
# TARGET HEATMAP
#############################################################################################################

msg("Making target heatmap...")

target_temp_effect <- merge(
  matched_primary[, .(target_input, Protein_group, Genes)],
  temp_eff_dt,
  by = c("Protein_group", "Genes"),
  all.x = TRUE
)

if (nrow(target_temp_effect) > 0) {
  heat_dt <- dcast(
    target_temp_effect,
    target_input + Genes + Protein_group ~ Temp,
    value.var = "mean_log2_ratio"
  )
  
  if (nrow(heat_dt) > 0) {
    heat_mat <- as.matrix(heat_dt[, -(1:3)])
    rownames(heat_mat) <- paste0(heat_dt$target_input, " | ", heat_dt$Genes)
    
    png(file.path(output_dir, "Target_heatmap_mean_log2_Drug_vs_DMSO.png"),
        width = 1800, height = max(1200, 25 * nrow(heat_mat)), res = 200)
    
    pheatmap(
      heat_mat,
      scale = "none",
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      border_color = NA,
      main = "Matched targets: mean log2(Drug/DMSO) across temperatures",
      fontsize_row = 7,
      fontsize_col = 10
    )
    
    dev.off()
  }
}

#############################################################################################################
# OPTIONAL SUMMARY SCATTER
#############################################################################################################

if (nrow(target_summary) > 0) {
  p_scatter <- ggplot(
    target_summary,
    aes(x = deltaTm_mean, y = max_abs_mean_log2_ratio, label = target_input)
  ) +
    geom_point(size = 2, alpha = 0.8) +
    geom_text(size = 2.5, vjust = -0.4, check_overlap = TRUE) +
    theme_bw(base_size = 11) +
    labs(
      title = "Requested targets in CETSA-MS SILAC",
      x = "Mean deltaTm (Drug - DMSO)",
      y = "Max absolute mean log2(Drug/DMSO)"
    )
  
  ggsave(
    file.path(output_dir, "Target_summary_scatter.png"),
    p_scatter,
    width = 11,
    height = 8,
    dpi = plot_dpi,
    bg = "white"
  )
}

#############################################################################################################
# CONSOLE SUMMARY
#############################################################################################################

msg("Done.")
msg("Output folder: ", output_dir)
msg("Requested targets: ", nrow(targets_dt))
msg("Matched targets: ", nrow(matched_primary))
msg("Unmatched targets: ", nrow(unmatched_targets))

if (nrow(target_summary) > 0) {
  msg("\nTop matched targets:")
  print(
    target_summary[, .(
      target_input,
      Genes,
      Protein_group,
      match_type,
      deltaTm_mean = round(deltaTm_mean, 2),
      deltaTm_mean_curve = round(deltaTm_mean_curve, 2),
      deltaAUC_mean = round(deltaAUC_mean, 3),
      max_abs_mean_log2_ratio = round(max_abs_mean_log2_ratio, 3)
    )][1:min(30, .N)]
  )
}