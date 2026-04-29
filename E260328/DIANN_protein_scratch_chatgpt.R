#############################################################################################################
#
# SILAC CETSA-MS / TPP-TR analysis from DIA-NN report.tsv
# Corrected "best-of-both" final script
#
# CORRECTIONS IN THIS VERSION (relative to previous version)
# ----------------------------------------------------------
# BUG 1 [Section 1]  - rename_table column alias 'suffix_num' corrected to 'Suffix_num'
# BUG 2 [Section 5]  - H/L ratio now falls back to Precursor.Quantity when
#                       Precursor.Translated is all-NA (common in first-pass workflows)
# BUG 3 [Section 15] - make_slope_bins() now removes row_id_tmp__ before returning,
#                       preventing it from polluting curve_stats and all output files
# BUG 4 [Section 16] - missing_boot diagnostic replaced: the old code used
#                       [is.na(Replicate.y)] on a merge where Replicate was a join key,
#                       so no .y suffix was added and the column did not exist (crash).
#                       Now uses fsetdiff() for a correct anti-join.
# BUG 5 [Section 19] - pairwise_min_prec moved inside the data.table := expression
#                       to eliminate the order-dependent external vector alignment risk.
# BUG 6 [Section 23] - pNAs y-axis upper limit guarded against empty data_completeness
#                       table (max() returns -Inf on empty input, crashing ggplot).
#
#############################################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(iq)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(minpack.lm)
  library(grid)
})

#############################################################################################################
# USER PARAMETERS
#############################################################################################################

input_report  <- "report.tsv"
workflow_type <- "second"          # "second" or "first" pass DIA-NN workflow
dataset_name  <- "intact"          # change to "extract" for cell extract run

expected_n_temps             <- 11
expected_replicates_per_temp <- 2

min_points_per_curve  <- 6
curve_fit_clamp_upper <- 1.5
curve_fit_clamp_lower <- 0.0

r2_cutoff              <- 0.80
plateau_vehicle_cutoff <- 0.30
min_slope_cutoff       <- -0.06

# Bootstrap curves lack cross-protein global normalization, so their fitted plateaus
# are typically higher than the main analysis curves. Use a relaxed threshold only there.
bootstrap_plateau_cutoff <- 0.60

min_bin_size              <- 300
padj_strict               <- 0.05
padj_lenient              <- 0.10
per_replicate_padj_cutoff <- 0.10   # used for pass_statistical_hit

null_log2HL_abs_cutoff    <- 0.15
null_min_proteins_per_bin <- 50

do_precursor_bootstrap   <- TRUE
bootstrap_B              <- 100     # increase to 300-500 for publication-grade runs
bootstrap_min_precursors <- 2
bootstrap_seed           <- 123

require_two_replicates              <- TRUE
require_min_precursors_each_channel <- 2

top_n_hit_plots     <- 25
write_intermediates <- TRUE
plot_dir            <- "TPPTR_plots"

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

diag_log <- file.path(plot_dir, paste0("diagnostics_", dataset_name, ".txt"))
diag <- function(...) {
  msg <- paste0(...)
  message(msg)
  cat(msg, "\n", file = diag_log, append = TRUE)
}
cat("", file = diag_log)
diag("=== SILAC CETSA-MS DIAGNOSTIC LOG ===")
diag("Dataset:     ", dataset_name)
diag("Run started: ", format(Sys.time()))

#############################################################################################################
# LOAD THE REPORT
#############################################################################################################

diag("\n--- LOADING REPORT ---")
report <- fread(
  input_report,
  select = c(
    "Run", "Protein.Group", "Protein.Ids", "First.Protein.Description", "Genes",
    "Stripped.Sequence", "Precursor.Id", "Proteotypic",
    "Precursor.Quantity", "Precursor.Translated", "Channel.Q.Value",
    "Q.Value", "Global.Q.Value", "PG.Q.Value", "Global.PG.Q.Value",
    "Lib.Q.Value", "Lib.PG.Q.Value", "Precursor.Charge", "RT",
    "Mass.Evidence", "RT.Start", "RT.Stop"
  )
)
report <- as.data.table(report)
diag("Rows loaded:  ", nrow(report))
diag("Unique runs:  ", uniqueN(report$Run))

#############################################################################################################
# SECTION 1. RUN NAME CLEANUP AND SAMPLE ANNOTATION
#############################################################################################################

diag("\n--- SECTION 1: RUN ANNOTATION ---")

report[, Run           := gsub("-", "_", Run)]
report[, Had_timestamp := grepl("_[0-9]{14}$", Run)]
report[, Run_clean     := sub("_[0-9]{14}$", "", Run)]
report[, Temp          := as.integer(sub(".*_CET_(\\d+)_\\d+$", "\\1", Run_clean))]
report[, Suffix_num    := suppressWarnings(as.integer(sub(".*_(\\d+)$", "\\1", Run_clean)))]

runs_dt <- unique(report[, .(Run, Run_clean, Temp, Suffix_num, Had_timestamp)])
setorder(runs_dt, Temp, Suffix_num)

diag("Temperatures found: ", paste(sort(unique(runs_dt$Temp)), collapse = ", "))
diag("Runs per temperature:")
print(runs_dt[, .N, by = Temp])

temp_check <- runs_dt[, .N, by = Temp][N != expected_replicates_per_temp]
if (nrow(temp_check) > 0)
  diag("WARNING - unexpected run counts per temp: ", paste(temp_check$Temp, collapse = ", "))

n_temps <- uniqueN(runs_dt$Temp)
if (n_temps != expected_n_temps) {
  diag("WARNING - expected ", expected_n_temps, " temps but found ", n_temps)
} else {
  diag("OK: correct number of temperatures (", n_temps, ")")
}

# Safer replicate mapping: use run suffix if it is consistent across all temperatures
suffix_pattern_ok <- FALSE
if (all(!is.na(runs_dt$Suffix_num))) {
  suffix_by_temp <- runs_dt[,
                            .(suffix_set = paste(sort(unique(Suffix_num)), collapse = ",")),
                            by = Temp
  ]
  if (nrow(unique(suffix_by_temp[, .(suffix_set)])) == 1) {
    ref_set <- sort(unique(runs_dt$Suffix_num))
    if (length(ref_set) == expected_replicates_per_temp) {
      suffix_pattern_ok <- TRUE
      diag("Replicate assignment: using run suffix as replicate identity (safer mode)")
    }
  }
}

if (suffix_pattern_ok) {
  ref_suffix  <- sort(unique(runs_dt$Suffix_num))
  suffix_map  <- data.table(Suffix_num = ref_suffix, Replicate_num = seq_along(ref_suffix))
  runs_dt     <- merge(runs_dt, suffix_map, by = "Suffix_num", all.x = TRUE, sort = FALSE)
} else {
  diag("WARNING - suffix could not be used as replicate identity. Falling back to within-temp ordering.")
  setorder(runs_dt, Temp, Suffix_num, Run_clean)
  runs_dt[, Replicate_num := seq_len(.N), by = Temp]
}

runs_dt[, Replicate_old := Suffix_num]
runs_dt[, Replicate     := paste0("R", Replicate_num)]
runs_dt[, Sample        := paste0("CET_", Temp, "_", Replicate)]
runs_dt[, Run_base      := sub("_(\\d+)$", "", Run_clean)]
runs_dt[, Run_corrected := paste0(Run_base, "_", Replicate_num)]
runs_dt[, Needs_rename  := Had_timestamp | is.na(Replicate_old) | (Replicate_old != Replicate_num)]

report <- runs_dt[report, on = .(Run, Run_clean, Temp, Suffix_num, Had_timestamp)]

rename_table <- runs_dt[Needs_rename == TRUE,
                        .(old_run           = Run,
                          new_run           = Run_corrected,
                          Temp              = Temp,
                          Suffix_num        = Suffix_num,        # BUG 1 FIX: was lowercase 'suffix_num'
                          assigned_replicate = Replicate)]
if (nrow(rename_table) > 0) {
  fwrite(rename_table, "run_rename_table.tsv", sep = "\t")
  diag("Renamed runs written: ", nrow(rename_table), " entries")
}

#############################################################################################################
# SECTION 2. FDR FILTERING
#############################################################################################################

diag("\n--- SECTION 2: FDR FILTERING ---")
diag("Rows before filter: ", nrow(report))

apply_filters <- function(dt, filters) {
  for (f in filters) {
    col <- f$col
    if (col %in% names(dt)) {
      before <- nrow(dt)
      dt     <- dt[get(col) <= f$thr]
      diag("  ", col, " <= ", f$thr, ": removed ", before - nrow(dt), ", remain ", nrow(dt))
    }
  }
  dt
}

if (workflow_type == "second") {
  filters <- list(
    list(col = "Q.Value",         thr = 0.01),
    list(col = "PG.Q.Value",      thr = 0.05),
    list(col = "Lib.Q.Value",     thr = 0.01),
    list(col = "Lib.PG.Q.Value",  thr = 0.01),
    list(col = "Channel.Q.Value", thr = 0.01)
  )
  report <- apply_filters(report, filters)
  diag("Second-pass FDR applied")
} else if (workflow_type == "first") {
  filters <- list(
    list(col = "Q.Value",           thr = 0.01),
    list(col = "PG.Q.Value",        thr = 0.05),
    list(col = "Global.Q.Value",    thr = 0.01),
    list(col = "Global.PG.Q.Value", thr = 0.01),
    list(col = "Channel.Q.Value",   thr = 0.01)
  )
  report <- apply_filters(report, filters)
  diag("First-pass FDR applied")
} else {
  stop("workflow_type must be 'first' or 'second'")
}

diag("Rows after filter: ", nrow(report))

#############################################################################################################
# SECTION 3. QUANTITY CLEANUP AND NORMALIZATION
#############################################################################################################

diag("\n--- SECTION 3: NORMALIZATION ---")

if (!("Precursor.Quantity" %in% names(report))) stop("Missing column: Precursor.Quantity")
report[Precursor.Quantity == 0, Precursor.Quantity := NA_real_]

if ("Precursor.Translated" %in% names(report)) {
  report[Precursor.Translated == 0, Precursor.Translated := NA_real_]
} else {
  report[, Precursor.Translated := NA_real_]
}

report <- report[!is.na(Precursor.Quantity)]
diag("Rows after zero/NA removal: ", nrow(report))

target_sum <- report[, .(sumQ = sum(Precursor.Quantity, na.rm = TRUE)), by = Run][, median(sumQ)]
diag("Normalization target sum: ", round(target_sum, 0))

report[, Precursor.Quantity :=
         Precursor.Quantity / sum(Precursor.Quantity, na.rm = TRUE) * target_sum,
       by = Run]

# Normalize Precursor.Translated only in runs where it has non-missing values
report[,
       Precursor.Translated := {
         if (all(is.na(Precursor.Translated))) {
           rep(NA_real_, .N)
         } else {
           den <- sum(Precursor.Translated, na.rm = TRUE)
           if (!is.finite(den) || den <= 0) rep(NA_real_, .N)
           else Precursor.Translated / den * target_sum
         }
       },
       by = Run
]

#############################################################################################################
# SECTION 4. SILAC LABEL PARSING
#############################################################################################################

diag("\n--- SECTION 4: SILAC LABELS ---")

if (!("Precursor.Id" %in% names(report))) stop("Missing column: Precursor.Id")

report[, label := NA_character_]
report[grepl("SILAC-[A-Z]-L", Precursor.Id), label := "L"]
report[grepl("SILAC-[A-Z]-H", Precursor.Id), label := "H"]

n_L    <- report[label == "L", .N]
n_H    <- report[label == "H", .N]
n_none <- report[is.na(label), .N]

diag("L (vehicle): ", n_L, "  H (drug): ", n_H, "  unlabelled: ", n_none)

if (n_L == 0 || n_H == 0) {
  diag("CRITICAL: One SILAC channel is empty — check Precursor.Id patterns")
  print(head(report$Precursor.Id, 10))
}

report[, Precursor.Id.nolabels := gsub("SILAC-[A-Z]-[LH]", "SILAC", Precursor.Id)]

#############################################################################################################
# SECTION 5. H/L RATIOS AS SECONDARY QC
#
# BUG 2 FIX: Precursor.Translated is sometimes all-NA in first-pass DIA-NN workflows.
# If it is all-NA, fall back to Precursor.Quantity for ratio computation.
# Log which quantity was used.
#############################################################################################################

diag("\n--- SECTION 5: H/L RATIOS ---")

n_trans_nonNA <- report[!is.na(Precursor.Translated), .N]
if (n_trans_nonNA == 0) {
  diag("WARNING: Precursor.Translated is all-NA. Falling back to Precursor.Quantity for H/L ratios.")
  ratio_col <- "Precursor.Quantity"
} else {
  diag("Using Precursor.Translated for H/L ratio computation (", n_trans_nonNA, " non-NA values)")
  ratio_col <- "Precursor.Translated"
}

SILAC <- dcast(
  report[!is.na(label)],
  formula       = Precursor.Id.nolabels + Protein.Group + Run ~ label,
  value.var     = ratio_col,
  fun.aggregate = mean,
  na.rm         = TRUE
)
SILAC <- SILAC[!is.na(H) & !is.na(L) & H > 0 & L > 0]
SILAC[, precursor_log2_HL_ratio := log2(H / L)]
SILAC <- SILAC[, .(
  log2_HL_ratio = median(precursor_log2_HL_ratio, na.rm = TRUE),
  N_ratios      = .N
), by = .(Run, Protein.Group)]

diag("Protein-run H/L pairs: ", nrow(SILAC))
if (nrow(SILAC) == 0)
  diag("WARNING: No H/L pairs computed. Null protein identification in Section 15 will fall back to bin_global.")

#############################################################################################################
# SECTION 6. PROTEIN QUANTIFICATION (MaxLFQ per channel)
#############################################################################################################

diag("\n--- SECTION 6: MAXLFQ ---")

LFQ_fun <- function(dt, LFQ_colname) {
  if (nrow(dt) == 0) {
    return(data.table(Protein.Group = character(), Run = character())[, (LFQ_colname) := numeric()])
  }
  tmp <- dt[, .(
    protein_list = Protein.Group,
    sample_list  = Run,
    id           = Precursor.Id,
    quant        = log2(Precursor.Quantity)
  )]
  tmp <- fast_MaxLFQ(tmp)
  tmp <- data.table(tmp$estimate, Precursor_group = tmp$annotation, keep.rownames = "Protein.Group")
  tmp <- melt(tmp, id.vars = c("Protein.Group", "Precursor_group"),
              variable.name = "Run", value.name = LFQ_colname)
  tmp <- tmp[Precursor_group == "" & !is.na(get(LFQ_colname))][, -"Precursor_group"]
  tmp[, (LFQ_colname) := 2^get(LFQ_colname)]
  tmp
}

LFQ_T <- LFQ_fun(report,               "LFQ")
LFQ_H <- LFQ_fun(report[label == "H"], "LFQ_H")
LFQ_L <- LFQ_fun(report[label == "L"], "LFQ_L")

diag("LFQ rows — total: ", nrow(LFQ_T), "  H: ", nrow(LFQ_H), "  L: ", nrow(LFQ_L))

#############################################################################################################
# SECTION 7. INTENSITIES, COUNTS, ANNOTATIONS
#############################################################################################################

diag("\n--- SECTION 7: COUNTS AND ANNOTATIONS ---")

int_T <- report[,             .(Intensity   = sum(Precursor.Quantity, na.rm = TRUE)), by = .(Protein.Group, Run)]
int_H <- report[label == "H", .(Intensity_H = sum(Precursor.Quantity, na.rm = TRUE)), by = .(Protein.Group, Run)]
int_L <- report[label == "L", .(Intensity_L = sum(Precursor.Quantity, na.rm = TRUE)), by = .(Protein.Group, Run)]

if (!("Proteotypic" %in% names(report))) report[, Proteotypic := 0L]

# Count unique proteotypic precursor IDs (not row sums, which overcount multi-PSM precursors)
counts_H <- report[label == "H", .(
  N_precursors_H             = as.numeric(uniqueN(Precursor.Id)),
  N_precursors_proteotypic_H = as.numeric(uniqueN(Precursor.Id[Proteotypic %in% c(1, TRUE)]))
), by = .(Protein.Group, Run)]

counts_L <- report[label == "L", .(
  N_precursors_L             = as.numeric(uniqueN(Precursor.Id)),
  N_precursors_proteotypic_L = as.numeric(uniqueN(Precursor.Id[Proteotypic %in% c(1, TRUE)]))
), by = .(Protein.Group, Run)]

ann_cols <- intersect(
  c("Protein.Group", "Run", "Sample", "Temp", "Replicate",
    "First.Protein.Description", "Genes"),
  names(report)
)
annotations <- unique(report[, ..ann_cols])

merge2 <- function(x, y) merge(x, y, by = c("Protein.Group", "Run"), all = TRUE)

proteins <- Reduce(merge2, list(
  LFQ_T, LFQ_H, LFQ_L,
  int_T, int_H, int_L,
  counts_H, counts_L,
  SILAC,
  annotations
))
proteins <- proteins[!is.na(Intensity)]

if ("First.Protein.Description" %in% names(proteins)) {
  setnames(proteins,
           c("Protein.Group", "First.Protein.Description"),
           c("Protein_group", "Protein_description"))
} else {
  setnames(proteins, "Protein.Group", "Protein_group")
}

diag("Protein-run rows: ",    nrow(proteins),
     "  Unique proteins: ",   uniqueN(proteins$Protein_group),
     "  Unique replicates: ", paste(sort(unique(proteins$Replicate)), collapse = ", "))

preferred_first <- intersect(
  c("Run", "Sample", "Temp", "Replicate", "Protein_group", "Protein_description", "Genes",
    "LFQ", "LFQ_H", "LFQ_L", "Intensity", "Intensity_H", "Intensity_L",
    "log2_HL_ratio", "N_ratios", "N_precursors_H", "N_precursors_L",
    "N_precursors_proteotypic_H", "N_precursors_proteotypic_L"),
  names(proteins)
)
setcolorder(proteins, c(preferred_first, setdiff(names(proteins), preferred_first)))

fwrite(proteins, "protein_quant.tsv.gz", sep = "\t")
if (write_intermediates) fwrite(report, "report_filtered.tsv.gz", sep = "\t")

#############################################################################################################
# SECTION 8. PRECURSOR TABLE FOR BOOTSTRAP
#############################################################################################################

diag("\n--- SECTION 8: PRECURSOR TABLE ---")

precursor_run <- report[
  !is.na(label) & !is.na(Precursor.Quantity) & Precursor.Quantity > 0,
  .(
    precursor_qty       = median(Precursor.Quantity,   na.rm = TRUE),
    precursor_qty_trans = median(Precursor.Translated, na.rm = TRUE),
    Proteotypic         = max(Proteotypic,             na.rm = TRUE)
  ),
  by = .(Protein.Group, Precursor.Id.nolabels, label, Run)
]

precursor_run <- merge(
  precursor_run,
  unique(report[, .(Run, Sample, Temp, Replicate)]),
  by = "Run", all.x = TRUE
)
setnames(precursor_run, "Protein.Group", "Protein_group")
precursor_run[, Temp_num := as.numeric(Temp)]
precursor_run <- precursor_run[
  is.finite(precursor_qty) & precursor_qty > 0 &
    !is.na(Protein_group) & !is.na(Replicate) & !is.na(Temp_num)
]

diag("Precursor-run rows: ",       nrow(precursor_run))
diag("Unique proteins in table: ", uniqueN(precursor_run$Protein_group))

if (write_intermediates)
  fwrite(precursor_run, "precursor_run_for_bootstrap.tsv.gz", sep = "\t")

#############################################################################################################
# SECTION 9. BUILD CURVE INPUT
#############################################################################################################

diag("\n--- SECTION 9: CURVE INPUT ---")

proteins[, Temp_num := as.numeric(as.character(Temp))]
proteins <- proteins[!is.na(Temp_num)]

curve_input <- rbindlist(list(
  proteins[!is.na(LFQ_L) & LFQ_L > 0,
           .(Protein_group, Protein_description, Genes, Run, Sample, Replicate,
             Temp = Temp_num, Condition = "L", Quantity = LFQ_L,
             N_precursors = N_precursors_L, log2_HL_ratio)],
  proteins[!is.na(LFQ_H) & LFQ_H > 0,
           .(Protein_group, Protein_description, Genes, Run, Sample, Replicate,
             Temp = Temp_num, Condition = "H", Quantity = LFQ_H,
             N_precursors = N_precursors_H, log2_HL_ratio)]
), use.names = TRUE, fill = TRUE)

setorder(curve_input, Protein_group, Replicate, Condition, Temp)

diag("Curve input rows — L: ",           curve_input[Condition == "L", .N],
     "  H: ",                            curve_input[Condition == "H", .N])
diag("Unique proteins in curve input: ", uniqueN(curve_input$Protein_group))

#############################################################################################################
# SECTION 10. PER-CURVE NORMALIZATION TO LOWEST TEMPERATURE
#############################################################################################################

diag("\n--- SECTION 10: PER-CURVE NORMALIZATION ---")

curve_input[, base_quantity := Quantity[which.min(Temp)],
            by = .(Protein_group, Replicate, Condition)]
curve_input[!is.finite(base_quantity) | base_quantity <= 0, base_quantity := NA_real_]
curve_input[, FC_raw := Quantity / base_quantity]

diag("Rows with missing base quantity: ", curve_input[is.na(base_quantity), .N])

#############################################################################################################
# SECTION 11. GLOBAL NORMALIZATION ACROSS PROTEINS
#############################################################################################################

diag("\n--- SECTION 11: GLOBAL NORMALIZATION ---")

curve_input[, global_median_FC := median(FC_raw, na.rm = TRUE),
            by = .(Replicate, Condition, Temp)]
curve_input[!is.finite(global_median_FC) | global_median_FC <= 0, global_median_FC := NA_real_]
curve_input[, FC_norm := FC_raw / global_median_FC]
curve_input[, FC_fit  := pmin(pmax(FC_norm, curve_fit_clamp_lower), curve_fit_clamp_upper)]
curve_input <- curve_input[is.finite(FC_fit)]

diag("Rows after normalization and clamping: ", nrow(curve_input))
diag("Unique proteins after normalization: ",   uniqueN(curve_input$Protein_group))

if (write_intermediates)
  fwrite(curve_input, "protein_curve_input.tsv.gz", sep = "\t")

#############################################################################################################
# SECTION 12. CURVE FITTING HELPERS
#############################################################################################################

melt_fun <- function(T, a, b, plateau) {
  plateau + (1 - plateau) / (1 + exp(b * (T - a)))
}

calc_tm <- function(a, b, plateau) {
  if (!is.finite(a) || !is.finite(b) || !is.finite(plateau)) return(NA_real_)
  if (plateau >= 0.5) return(NA_real_)
  denom <- 0.5 - plateau
  if (denom <= 0) return(NA_real_)
  ratio <- (1 - plateau) / denom - 1
  if (!is.finite(ratio) || ratio <= 0) return(NA_real_)
  a + log(ratio) / b
}

calc_slope_at_inflection <- function(b, plateau) {
  if (!is.finite(b) || !is.finite(plateau)) return(NA_real_)
  -b * (1 - plateau) / 4
}

calc_r2 <- function(obs, pred) {
  keep <- is.finite(obs) & is.finite(pred)
  obs  <- obs[keep]
  pred <- pred[keep]
  if (length(obs) < 3) return(NA_real_)
  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  if (ss_tot <= 0) return(NA_real_)
  1 - ss_res / ss_tot
}

fit_one_curve <- function(df) {
  df <- as.data.table(copy(df))
  df <- df[order(Temp)]
  df <- df[is.finite(Temp) & is.finite(FC_fit)]
  
  out <- list(
    n_points = nrow(df), a = NA_real_, b = NA_real_, plateau = NA_real_,
    Tm = NA_real_, slope = NA_real_, R2 = NA_real_,
    fit_ok = FALSE, fit_msg = NA_character_
  )
  
  if (nrow(df) < min_points_per_curve) {
    out$fit_msg <- sprintf("Too few points (<%d)", min_points_per_curve)
    return(as.data.table(out))
  }
  
  temp_range <- range(df$Temp)
  temp_span  <- diff(temp_range)
  
  start_plateau <- mean(sort(df$FC_fit)[1:min(3, nrow(df))], na.rm = TRUE)
  start_plateau <- min(max(start_plateau, 0.01), 0.79)
  
  mid_temp <- temp_range[1] + temp_span * 0.4
  df_right <- df[Temp >= mid_temp][order(Temp)]
  if (nrow(df_right) >= 2) {
    start_a <- df_right$Temp[which.min(abs(df_right$FC_fit - 0.5))]
  } else {
    start_a <- temp_range[1] + temp_span * 0.7
  }
  start_a <- min(max(start_a, temp_range[1] - 5), temp_range[2] + 5)
  
  start_list <- list(
    list(a = start_a,                         b = 0.3,  plateau = start_plateau),
    list(a = temp_range[1] + temp_span * 0.7, b = 0.3,  plateau = start_plateau),
    list(a = temp_range[1] + temp_span * 0.6, b = 0.5,  plateau = start_plateau),
    list(a = temp_range[1] + temp_span * 0.8, b = 0.2,  plateau = start_plateau),
    list(a = temp_range[1] + temp_span * 0.5, b = 0.4,  plateau = start_plateau * 0.5),
    list(a = temp_range[1] + temp_span * 0.9, b = 0.15, plateau = start_plateau)
  )
  
  fit      <- NULL
  last_msg <- "All starting points failed"
  
  for (st in start_list) {
    attempt <- try(
      nlsLM(
        FC_fit ~ plateau + (1 - plateau) / (1 + exp(b * (Temp - a))),
        data    = df, start = st,
        lower   = c(a = temp_range[1] - 15, b = 0.001, plateau = 0),
        upper   = c(a = temp_range[2] + 15, b = 10,    plateau = 0.8),
        control = nls.lm.control(maxiter = 300, ftol = 1e-6, ptol = 1e-6)
      ),
      silent = TRUE
    )
    if (!inherits(attempt, "try-error")) {
      fit <- attempt
      break
    }
    last_msg <- as.character(attempt)
  }
  
  if (is.null(fit)) {
    out$fit_msg <- last_msg
    return(as.data.table(out))
  }
  
  co          <- coef(fit)
  pred        <- predict(fit, newdata = df)
  out$a       <- unname(co["a"])
  out$b       <- unname(co["b"])
  out$plateau <- unname(co["plateau"])
  out$Tm      <- calc_tm(out$a, out$b, out$plateau)
  out$slope   <- calc_slope_at_inflection(out$b, out$plateau)
  out$R2      <- calc_r2(df$FC_fit, pred)
  out$fit_ok  <- TRUE
  out$fit_msg <- "OK"
  as.data.table(out)
}

#############################################################################################################
# SECTION 13. FIT ALL CURVES
#############################################################################################################

diag("\n--- SECTION 13: CURVE FITTING ---")

curve_fits_long <- curve_input[
  , fit_one_curve(.SD),
  by = .(Protein_group, Protein_description, Genes, Replicate, Condition)
]

fwrite(curve_fits_long, "protein_curve_fits_long.tsv.gz", sep = "\t")

diag("Curves fitted OK:   ", curve_fits_long[fit_ok == TRUE,  .N])
diag("Curves failed:      ", curve_fits_long[fit_ok == FALSE, .N])
diag("Failure reasons:")
print(curve_fits_long[fit_ok == FALSE, .N, by = fit_msg][order(-N)])

curve_fits <- dcast(
  curve_fits_long,
  Protein_group + Protein_description + Genes + Replicate ~ Condition,
  value.var = c("n_points", "a", "b", "plateau", "Tm", "slope", "R2", "fit_ok", "fit_msg")
)

for (nm in names(curve_fits)) {
  if (grepl("_L$", nm)) setnames(curve_fits, nm, sub("_L$", "_vehicle",   nm))
  if (grepl("_H$", nm)) setnames(curve_fits, nm, sub("_H$", "_treatment", nm))
}

curve_fits[, deltaTm          := Tm_treatment - Tm_vehicle]
curve_fits[, min_slope_pair   := pmin(slope_vehicle, slope_treatment, na.rm = TRUE)]
curve_fits[, max_plateau_pair := pmax(plateau_vehicle, plateau_treatment, na.rm = TRUE)]

diag("Protein-replicate pairs after dcast: ", nrow(curve_fits))
diag("Unique proteins with any curve: ",      uniqueN(curve_fits$Protein_group))

#############################################################################################################
# SECTION 14. CURVE QUALITY FILTERS
#############################################################################################################

diag("\n--- SECTION 14: QUALITY FILTERS ---")

curve_fits[, pass_curve_quality :=
             (fit_ok_vehicle   == TRUE) &
             (fit_ok_treatment == TRUE) &
             is.finite(R2_vehicle)      & R2_vehicle      > r2_cutoff &
             is.finite(R2_treatment)    & R2_treatment    > r2_cutoff &
             is.finite(plateau_vehicle) & plateau_vehicle < plateau_vehicle_cutoff &
             is.finite(min_slope_pair)  & min_slope_pair  < min_slope_cutoff]

n_pass <- curve_fits[pass_curve_quality == TRUE, .N]
n_fail <- curve_fits[pass_curve_quality != TRUE | is.na(pass_curve_quality), .N]
diag("Passing quality: ", n_pass)
diag("Failing quality: ", n_fail)

diag("Failure reasons (can overlap):")
diag("  fit_ok_vehicle != TRUE:   ", curve_fits[!(fit_ok_vehicle   == TRUE) | is.na(fit_ok_vehicle),   .N])
diag("  fit_ok_treatment != TRUE: ", curve_fits[!(fit_ok_treatment == TRUE) | is.na(fit_ok_treatment), .N])
diag("  R2 vehicle failed:        ", curve_fits[!is.finite(R2_vehicle)      | R2_vehicle      <= r2_cutoff, .N])
diag("  R2 treatment failed:      ", curve_fits[!is.finite(R2_treatment)    | R2_treatment    <= r2_cutoff, .N])
diag("  Plateau failed:           ", curve_fits[!is.finite(plateau_vehicle) | plateau_vehicle >= plateau_vehicle_cutoff, .N])
diag("  Slope failed:             ", curve_fits[!is.finite(min_slope_pair)  | min_slope_pair  >= min_slope_cutoff, .N])

diag("Unique proteins passing quality (any replicate): ",
     uniqueN(curve_fits[pass_curve_quality == TRUE, Protein_group]))

fwrite(curve_fits, "protein_curve_fits_per_replicate.tsv.gz", sep = "\t")

#############################################################################################################
# SECTION 15. SLOPE-BINNED Z-TESTS
#
# BUG 3 FIX: make_slope_bins() now removes the row_id_tmp__ helper column
# before returning so it does not pollute curve_stats and all downstream output files.
#############################################################################################################

diag("\n--- SECTION 15: SLOPE-BINNED Z-TESTS ---")

make_slope_bins <- function(dt, slope_col = "min_slope_pair", min_bin_size = 300) {
  x <- copy(dt[is.finite(get(slope_col))])
  x <- x[order(get(slope_col))]
  x[, row_id_tmp__ := .I]
  x[, slope_bin    := ceiling(.I / min_bin_size)]
  x[, row_id_tmp__ := NULL]   # BUG 3 FIX: remove helper column before returning
  x[]
}

null_like <- proteins[, .(
  median_abs_log2HL = {
    v <- log2_HL_ratio[is.finite(log2_HL_ratio)]
    if (length(v) > 0) median(abs(v)) else NA_real_
  }
), by = Protein_group]

curve_stats <- merge(curve_fits, null_like, by = "Protein_group", all.x = TRUE)
curve_stats  <- make_slope_bins(curve_stats, min_bin_size = min_bin_size)

diag("Slope bins created: ", uniqueN(curve_stats$slope_bin))
diag("Proteins with near-zero HL (potential null): ",
     null_like[!is.na(median_abs_log2HL) & median_abs_log2HL <= null_log2HL_abs_cutoff, .N])

ztest_within_bin <- function(dt_bin) {
  dt_bin <- copy(dt_bin)
  
  null_dt <- dt_bin[is.finite(deltaTm) &
                      is.finite(median_abs_log2HL) &
                      median_abs_log2HL <= null_log2HL_abs_cutoff]
  null_source <- "near_zero_HL"
  
  if (nrow(null_dt) < null_min_proteins_per_bin) {
    null_dt     <- dt_bin[is.finite(deltaTm)]
    null_source <- "bin_global"
  }
  
  mu0 <- mean(null_dt$deltaTm, na.rm = TRUE)
  sd0 <- sd(null_dt$deltaTm,   na.rm = TRUE)
  
  if (!is.finite(mu0)) mu0 <- 0
  if (!is.finite(sd0) || sd0 <= 0) {
    dt_bin[, `:=`(
      null_mean = mu0, null_sd = NA_real_, null_source = null_source,
      z_score   = NA_real_, p_value = NA_real_
    )]
    return(dt_bin)
  }
  
  dt_bin[, `:=`(
    null_mean   = mu0,
    null_sd     = sd0,
    null_source = null_source,
    z_score     = (deltaTm - mu0) / sd0,
    p_value     = 2 * pnorm(-abs((deltaTm - mu0) / sd0))
  )]
  dt_bin
}

curve_stats <- curve_stats[, ztest_within_bin(.SD), by = .(Replicate, slope_bin)]
curve_stats[, p_adj_BH := p.adjust(p_value, method = "BH"), by = Replicate]

curve_stats[, deltaTm_exceeds_null :=
              is.finite(deltaTm) & is.finite(null_sd) &
              abs(deltaTm - null_mean) > abs(null_sd)]

diag("Null source distribution:")
print(curve_stats[, .N, by = .(Replicate, null_source)])
diag("p_adj_BH < 0.05 per replicate:")
print(curve_stats[is.finite(p_adj_BH) & p_adj_BH < 0.05, .N, by = Replicate])

#############################################################################################################
# SECTION 16. PRECURSOR BOOTSTRAP NULL
#
# BUG 4 FIX: The missing_boot diagnostic previously did:
#   merge(..., all.x=TRUE)[is.na(Replicate.y)]
# When both tables share the same join columns, data.table does not add .y
# suffixes to join key columns, so Replicate.y does not exist -> crash.
# Fixed by using fsetdiff() for a proper anti-join.
#############################################################################################################

diag("\n--- SECTION 16: BOOTSTRAP ---")

bootstrap_protein_curve <- function(dt_one_protein_rep, B = 100) {
  dt_one_protein_rep <- as.data.table(copy(dt_one_protein_rep))
  if (nrow(dt_one_protein_rep) == 0) return(NULL)
  
  nL <- uniqueN(dt_one_protein_rep[label == "L", Precursor.Id.nolabels])
  nH <- uniqueN(dt_one_protein_rep[label == "H", Precursor.Id.nolabels])
  if (nL < bootstrap_min_precursors || nH < bootstrap_min_precursors) return(NULL)
  
  temps <- sort(unique(dt_one_protein_rep$Temp_num))
  if (length(temps) < min_points_per_curve) return(NULL)
  
  make_mat <- function(dd) {
    dcast(dd, Precursor.Id.nolabels ~ Temp_num,
          value.var = "precursor_qty", fun.aggregate = median)
  }
  
  wideL <- make_mat(dt_one_protein_rep[label == "L"])
  wideH <- make_mat(dt_one_protein_rep[label == "H"])
  if (nrow(wideL) < bootstrap_min_precursors || nrow(wideH) < bootstrap_min_precursors) return(NULL)
  
  temp_cols <- intersect(setdiff(names(wideL), "Precursor.Id.nolabels"),
                         setdiff(names(wideH), "Precursor.Id.nolabels"))
  if (length(temp_cols) < min_points_per_curve) return(NULL)
  
  wideL <- wideL[, c("Precursor.Id.nolabels", temp_cols), with = FALSE]
  wideH <- wideH[, c("Precursor.Id.nolabels", temp_cols), with = FALSE]
  
  out_list <- vector("list", B)
  for (b in seq_len(B)) {
    sampL <- wideL[sample(.N, .N, replace = TRUE)]
    sampH <- wideH[sample(.N, .N, replace = TRUE)]
    
    bootL <- melt(sampL[, lapply(.SD, median, na.rm = TRUE), .SDcols = temp_cols],
                  measure.vars = temp_cols, variable.name = "Temp_num", value.name = "Quantity")
    bootL[, Condition := "L"]
    
    bootH <- melt(sampH[, lapply(.SD, median, na.rm = TRUE), .SDcols = temp_cols],
                  measure.vars = temp_cols, variable.name = "Temp_num", value.name = "Quantity")
    bootH[, Condition := "H"]
    
    boot <- rbindlist(list(bootL, bootH), use.names = TRUE)
    boot[, Temp_num := as.numeric(as.character(Temp_num))]
    boot[, boot_id  := b]
    out_list[[b]] <- boot
  }
  
  rbindlist(out_list, use.names = TRUE, fill = TRUE)
}

fit_bootstrap_deltaTm <- function(boot_curve_dt) {
  if (is.null(boot_curve_dt) || nrow(boot_curve_dt) == 0) return(NULL)
  
  boot_curve_dt <- as.data.table(copy(boot_curve_dt))
  boot_curve_dt <- boot_curve_dt[is.finite(Quantity) & Quantity > 0]
  if (nrow(boot_curve_dt) == 0) return(NULL)
  
  setorder(boot_curve_dt, boot_id, Condition, Temp_num)
  
  boot_curve_dt[, base_quantity := Quantity[which.min(Temp_num)], by = .(boot_id, Condition)]
  boot_curve_dt[!is.finite(base_quantity) | base_quantity <= 0, base_quantity := NA_real_]
  boot_curve_dt[, FC_raw := Quantity / base_quantity]
  boot_curve_dt[, FC_fit := pmin(pmax(FC_raw, curve_fit_clamp_lower), curve_fit_clamp_upper)]
  boot_curve_dt <- boot_curve_dt[is.finite(FC_fit)]
  if (nrow(boot_curve_dt) == 0) return(NULL)
  
  boot_fits <- boot_curve_dt[
    , fit_one_curve(.SD[, .(Temp = Temp_num, FC_fit)]),
    by = .(boot_id, Condition)
  ]
  
  boot_wide <- dcast(boot_fits, boot_id ~ Condition,
                     value.var = c("Tm", "R2", "plateau", "slope", "fit_ok"))
  
  for (nm in names(boot_wide)) {
    if (grepl("_L$", nm)) setnames(boot_wide, nm, sub("_L$", "_vehicle",   nm))
    if (grepl("_H$", nm)) setnames(boot_wide, nm, sub("_H$", "_treatment", nm))
  }
  
  boot_wide[, deltaTm_boot := Tm_treatment - Tm_vehicle]
  
  boot_wide[, pass_boot_quality :=
              (fit_ok_vehicle   == TRUE) &
              (fit_ok_treatment == TRUE) &
              is.finite(R2_vehicle)      & R2_vehicle      > r2_cutoff &
              is.finite(R2_treatment)    & R2_treatment    > r2_cutoff &
              is.finite(plateau_vehicle) & plateau_vehicle < bootstrap_plateau_cutoff &
              is.finite(pmin(slope_vehicle, slope_treatment, na.rm = TRUE)) &
              pmin(slope_vehicle, slope_treatment, na.rm = TRUE) < min_slope_cutoff]
  
  boot_wide
}

bootstrap_summary <- NULL
bootstrap_raw     <- NULL

if (do_precursor_bootstrap) {
  set.seed(bootstrap_seed)
  
  boot_candidates <- curve_fits[
    pass_curve_quality == TRUE & is.finite(deltaTm),
    .(Protein_group, Replicate)
  ]
  
  boot_counts <- precursor_run[, .(
    n_prec_L = uniqueN(Precursor.Id.nolabels[label == "L"]),
    n_prec_H = uniqueN(Precursor.Id.nolabels[label == "H"])
  ), by = .(Protein_group, Replicate)]
  
  boot_candidates <- merge(boot_candidates, boot_counts,
                           by = c("Protein_group", "Replicate"), all.x = TRUE)
  boot_candidates <- boot_candidates[
    n_prec_L >= bootstrap_min_precursors & n_prec_H >= bootstrap_min_precursors
  ]
  
  diag("Bootstrap candidates: ", nrow(boot_candidates), " protein-replicate curves")
  
  boot_list <- vector("list", nrow(boot_candidates))
  
  for (i in seq_len(nrow(boot_candidates))) {
    pid  <- boot_candidates$Protein_group[i]
    repi <- boot_candidates$Replicate[i]
    dt_sub <- precursor_run[Protein_group == pid & Replicate == repi]
    
    boot_fit <- try({
      bc <- bootstrap_protein_curve(dt_sub, B = bootstrap_B)
      fit_bootstrap_deltaTm(bc)
    }, silent = TRUE)
    
    if (!inherits(boot_fit, "try-error") && !is.null(boot_fit) && nrow(boot_fit) > 0) {
      boot_fit[, `:=`(Protein_group = pid, Replicate = repi)]
      boot_list[[i]] <- boot_fit
    }
    
    if (i %% 100 == 0)
      diag("Bootstrap progress: ", i, " / ", nrow(boot_candidates))
  }
  
  boot_list_clean <- Filter(Negate(is.null), boot_list)
  diag("Bootstrap raw fit success: ", length(boot_list_clean),
       " / ", nrow(boot_candidates), " candidates")
  
  if (length(boot_list_clean) > 0) {
    bootstrap_raw <- rbindlist(boot_list_clean, use.names = TRUE, fill = TRUE)
    
    diag("Bootstrap raw rows: ",            nrow(bootstrap_raw))
    diag("pass_boot_quality TRUE rows: ",   bootstrap_raw[pass_boot_quality == TRUE, .N])
    
    bootstrap_summary <- bootstrap_raw[
      pass_boot_quality == TRUE & is.finite(deltaTm_boot),
      .(
        n_boot                = .N,
        deltaTm_boot_mean     = mean(deltaTm_boot,                na.rm = TRUE),
        deltaTm_boot_sd       = sd(deltaTm_boot,                  na.rm = TRUE),
        deltaTm_boot_q025     = quantile(deltaTm_boot, 0.025,     na.rm = TRUE),
        deltaTm_boot_q975     = quantile(deltaTm_boot, 0.975,     na.rm = TRUE),
        deltaTm_boot_abs95    = quantile(abs(deltaTm_boot), 0.95, na.rm = TRUE),
        deltaTm_boot_ci_width = diff(quantile(deltaTm_boot, c(0.025, 0.975), na.rm = TRUE))
      ),
      by = .(Protein_group, Replicate)
    ]
    
    diag("Bootstrap summary entries: ", nrow(bootstrap_summary))
    diag("Candidates with usable bootstrap summaries: ",
         uniqueN(bootstrap_summary[, paste(Protein_group, Replicate)]))
    
    # BUG 4 FIX: use fsetdiff() for anti-join instead of merge()[is.na(Replicate.y)]
    # The old code failed because join-key columns never get .x/.y suffixes in data.table.
    boot_candidates_key <- boot_candidates[, .(Protein_group, Replicate)]
    boot_summary_key    <- bootstrap_summary[, .(Protein_group, Replicate)]
    missing_boot        <- fsetdiff(boot_candidates_key, boot_summary_key)
    diag("Candidates attempted but without usable bootstrap summary: ", nrow(missing_boot))
    
    fwrite(bootstrap_raw,     "protein_deltaTm_bootstrap_raw.tsv.gz",     sep = "\t")
    fwrite(bootstrap_summary, "protein_deltaTm_bootstrap_summary.tsv.gz", sep = "\t")
    
    curve_stats <- merge(curve_stats, bootstrap_summary,
                         by = c("Protein_group", "Replicate"), all.x = TRUE)
    
    curve_stats[, deltaTm_exceeds_bootstrap :=
                  is.finite(deltaTm) & is.finite(deltaTm_boot_abs95) &
                  abs(deltaTm) > deltaTm_boot_abs95]
    
    diag("Proteins exceeding bootstrap threshold per replicate:")
    print(curve_stats[deltaTm_exceeds_bootstrap == TRUE, .N, by = Replicate])
    
  } else {
    diag("WARNING: Bootstrap produced no usable raw results — falling back to null criterion")
    curve_stats[, deltaTm_exceeds_bootstrap := NA]
  }
} else {
  curve_stats[, deltaTm_exceeds_bootstrap := NA]
}

#############################################################################################################
# SECTION 17. REPLICATE-LEVEL STATISTICAL HIT CALL
#############################################################################################################

diag("\n--- SECTION 17: PER-REPLICATE HIT CALL ---")

curve_stats[, pass_statistical_hit :=
              pass_curve_quality == TRUE &
              is.finite(p_adj_BH) &
              p_adj_BH < per_replicate_padj_cutoff &
              fifelse(
                !is.na(deltaTm_exceeds_bootstrap),
                deltaTm_exceeds_bootstrap,
                deltaTm_exceeds_null
              )]

diag("Pass statistical hit per replicate:")
print(curve_stats[, .N, by = .(Replicate, pass_statistical_hit)])

fwrite(curve_stats, "protein_tpptr_stats_per_replicate.tsv.gz", sep = "\t")

#############################################################################################################
# SECTION 18. SUMMARIZE ACROSS REPLICATES
#############################################################################################################

diag("\n--- SECTION 18: CROSS-REPLICATE SUMMARY ---")

summary_hits <- curve_stats[
  ,
  {
    d_by_rep       <- setNames(deltaTm,                  Replicate)
    p_by_rep       <- setNames(p_adj_BH,                 Replicate)
    qf_by_rep      <- setNames(pass_curve_quality,       Replicate)
    en_null_by_rep <- setNames(deltaTm_exceeds_null,     Replicate)
    
    if ("deltaTm_exceeds_bootstrap" %in% names(.SD)) {
      en_boot_by_rep <- setNames(deltaTm_exceeds_bootstrap, Replicate)
    } else {
      en_boot_by_rep <- setNames(rep(NA, .N), Replicate)
    }
    
    d       <- unname(d_by_rep[c("R1", "R2")])
    p       <- unname(p_by_rep[c("R1", "R2")])
    qf      <- unname(qf_by_rep[c("R1", "R2")])
    en_null <- unname(en_null_by_rep[c("R1", "R2")])
    en_boot <- unname(en_boot_by_rep[c("R1", "R2")])
    
    en_final <- ifelse(is.na(en_boot), en_null, en_boot)
    
    # zero/zero deltaTm is not treated as directional agreement
    same_sign <- FALSE
    if (sum(is.finite(d)) == 2 && all(d != 0)) {
      same_sign <- sign(d[1]) == sign(d[2])
    }
    
    both_curve_quality      <- length(qf) == 2 && all(qf %in% TRUE)
    both_exceed_variability <- length(en_final) == 2 && all(en_final %in% TRUE)
    
    padj_rule <- FALSE
    if (sum(is.finite(p)) == 2) {
      padj_rule <- (p[1] < padj_strict & p[2] < padj_lenient) |
        (p[2] < padj_strict & p[1] < padj_lenient)
    }
    
    n_reps   <- uniqueN(Replicate)
    r2v_min  <- if (any(is.finite(R2_vehicle)))      min(R2_vehicle,      na.rm = TRUE) else NA_real_
    r2t_min  <- if (any(is.finite(R2_treatment)))    min(R2_treatment,    na.rm = TRUE) else NA_real_
    plat_max <- if (any(is.finite(plateau_vehicle))) max(plateau_vehicle, na.rm = TRUE) else NA_real_
    slp_min  <- if (any(is.finite(min_slope_pair)))  min(min_slope_pair,  na.rm = TRUE) else NA_real_
    p_min    <- if (any(is.finite(p_adj_BH)))        min(p_adj_BH,        na.rm = TRUE) else NA_real_
    
    .(
      n_replicates            = n_reps,
      deltaTm_R1              = d[1],
      deltaTm_R2              = d[2],
      mean_deltaTm            = mean(deltaTm,   na.rm = TRUE),
      median_deltaTm          = median(deltaTm, na.rm = TRUE),
      sd_deltaTm              = sd(deltaTm,     na.rm = TRUE),
      Tm_vehicle_mean         = mean(Tm_vehicle,   na.rm = TRUE),
      Tm_treatment_mean       = mean(Tm_treatment, na.rm = TRUE),
      R2_vehicle_min          = r2v_min,
      R2_treatment_min        = r2t_min,
      plateau_vehicle_max     = plat_max,
      min_slope_overall       = slp_min,
      p_adj_R1                = p[1],
      p_adj_R2                = p[2],
      best_p_adj              = p_min,
      both_curve_quality      = both_curve_quality,
      same_sign               = same_sign,
      both_exceed_variability = both_exceed_variability,
      padj_rule               = padj_rule,
      TPPTR_hit = {
        base_hit <- both_curve_quality & same_sign & both_exceed_variability & padj_rule
        if (require_two_replicates) base_hit & n_reps >= 2 else base_hit
      },
      hit_direction = fifelse(
        median(deltaTm, na.rm = TRUE) > 0, "stabilized",
        fifelse(median(deltaTm, na.rm = TRUE) < 0, "destabilized", "no_shift")
      )
    )
  },
  by = .(Protein_group, Protein_description, Genes)
]

summary_hits <- summary_hits[!is.na(Protein_group)]

diag("Total proteins in summary: ", nrow(summary_hits))
diag("Criterion breakdown:")
diag("  Both curve quality:      ", summary_hits[both_curve_quality      == TRUE, .N])
diag("  Same sign:               ", summary_hits[same_sign               == TRUE, .N])
diag("  Both exceed variability: ", summary_hits[both_exceed_variability == TRUE, .N])
diag("  Padj rule:               ", summary_hits[padj_rule               == TRUE, .N])
diag("  TPPTR_hit (pre-precursor filter): ", summary_hits[TPPTR_hit == TRUE, .N])

#############################################################################################################
# SECTION 19. SUPPORT METRICS, BOOTSTRAP CI, AND RANKING
#
# BUG 5 FIX: pairwise_min_prec is now computed inside the data.table := expression
# rather than as an external positional vector. The external vector approach was
# order-dependent: if rows were reordered between the vector assignment and the
# := use, values would silently misalign. Moving it inside := eliminates the risk.
#############################################################################################################

diag("\n--- SECTION 19: SUPPORT METRICS ---")

support_counts <- proteins[, .(
  max_precursors_L = {
    v <- N_precursors_L[is.finite(N_precursors_L)]
    if (length(v) > 0) as.numeric(max(v)) else NA_real_
  },
  max_precursors_H = {
    v <- N_precursors_H[is.finite(N_precursors_H)]
    if (length(v) > 0) as.numeric(max(v)) else NA_real_
  },
  median_abs_log2HL = {
    v <- log2_HL_ratio[is.finite(log2_HL_ratio)]
    if (length(v) > 0) median(abs(v)) else NA_real_
  }
), by = .(Protein_group, Protein_description, Genes)]

summary_hits <- merge(summary_hits, support_counts,
                      by = c("Protein_group", "Protein_description", "Genes"),
                      all.x = TRUE)

if (!is.null(bootstrap_summary) && nrow(bootstrap_summary) > 0) {
  boot_ci <- bootstrap_summary[, .(
    mean_boot_ci_width = mean(deltaTm_boot_ci_width, na.rm = TRUE),
    mean_boot_sd       = mean(deltaTm_boot_sd,       na.rm = TRUE),
    n_boot_reps        = .N
  ), by = Protein_group]
  
  summary_hits <- merge(summary_hits, boot_ci, by = "Protein_group", all.x = TRUE)
  summary_hits[, ci_quality := fifelse(
    !is.na(mean_boot_ci_width) & mean_boot_ci_width < 4, "narrow_ci",
    fifelse(!is.na(mean_boot_ci_width), "wide_ci", "no_bootstrap")
  )]
} else {
  summary_hits[, `:=`(
    mean_boot_ci_width = NA_real_,
    mean_boot_sd       = NA_real_,
    n_boot_reps        = NA_integer_,
    ci_quality         = "no_bootstrap"
  )]
}

n_before <- summary_hits[TPPTR_hit == TRUE, .N]
summary_hits[
  TPPTR_hit == TRUE &
    (is.na(max_precursors_L) | max_precursors_L < require_min_precursors_each_channel |
       is.na(max_precursors_H) | max_precursors_H < require_min_precursors_each_channel),
  TPPTR_hit := FALSE
]
diag("TPPTR_hit: ", n_before, " -> ", summary_hits[TPPTR_hit == TRUE, .N],
     " after precursor count filter")

summary_hits[, hit_tier := fifelse(
  TPPTR_hit == TRUE &
    is.finite(mean_deltaTm) & abs(mean_deltaTm) >= 4 &
    is.finite(best_p_adj)   & best_p_adj < 0.01 &
    !is.na(max_precursors_L) & max_precursors_L >= require_min_precursors_each_channel &
    !is.na(max_precursors_H) & max_precursors_H >= require_min_precursors_each_channel &
    (is.na(ci_quality) | ci_quality == "narrow_ci"),
  "Tier_1_high_confidence",
  fifelse(
    TPPTR_hit == TRUE,
    "Tier_2_supported",
    fifelse(
      both_curve_quality == TRUE & same_sign == TRUE &
        is.finite(mean_deltaTm) & abs(mean_deltaTm) >= 2,
      "Tier_3_candidate",
      "Not_significant"
    )
  )
)]

diag("Hit tier distribution:")
print(summary_hits[, .N, by = hit_tier][order(-N)])

safe_rescale <- function(x, w) {
  r <- range(x, na.rm = TRUE)
  if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) return(rep(0, length(x)))
  scales::rescale(x, to = c(0, 1), from = r) * w
}

# BUG 5 FIX: pmin computed inline inside := so it is always aligned to the current row order
summary_hits[, rank_score :=
               safe_rescale(abs(mean_deltaTm), 0.45) +
               safe_rescale(-log10(pmax(best_p_adj, 1e-300)), 0.35) +
               safe_rescale(
                 {
                   v <- pmin(max_precursors_L, max_precursors_H, na.rm = TRUE)
                   v[!is.finite(v)] <- NA_real_
                   v
                 },
                 0.20
               )]

summary_hits[, abs_mean_deltaTm := abs(mean_deltaTm)]
setorder(summary_hits, -TPPTR_hit, -rank_score, -abs_mean_deltaTm)
summary_hits[, abs_mean_deltaTm := NULL]
summary_hits[, dataset := dataset_name]

diag("\nTop 20 hits by rank_score:")
print(summary_hits[1:min(20, .N),
                   .(Protein_group, Genes, mean_deltaTm, best_p_adj,
                     ci_quality, hit_tier, rank_score)])

#############################################################################################################
# SECTION 20. DIRECT / INDIRECT SCAFFOLD
#############################################################################################################

direct_indirect_scaffold <- copy(summary_hits)
direct_indirect_scaffold[, extract_dataset_available := FALSE]
direct_indirect_scaffold[, TPPTR_hit_extract         := NA]
direct_indirect_scaffold[, mean_deltaTm_extract      := NA_real_]
direct_indirect_scaffold[, same_direction_extract    := NA]
direct_indirect_scaffold[, classification := fifelse(
  TPPTR_hit == TRUE, "intact_hit_only_pending_extract", "no_intact_hit"
)]

fwrite(summary_hits,
       sprintf("protein_tpptr_hits_summary_%s.tsv.gz", dataset_name), sep = "\t")
fwrite(direct_indirect_scaffold,
       sprintf("direct_indirect_scaffold_%s.tsv.gz", dataset_name), sep = "\t")
fwrite(summary_hits[TPPTR_hit == TRUE],
       sprintf("protein_tpptr_hits_only_%s.tsv.gz", dataset_name), sep = "\t")

diag("\nOutput files written:")
diag("  protein_tpptr_hits_summary_", dataset_name, ".tsv.gz")
diag("  protein_tpptr_hits_only_", dataset_name, ".tsv.gz")

#############################################################################################################
# SECTION 21. FITTED CURVE LINES FOR PLOTTING
#############################################################################################################

diag("\n--- SECTION 21: FITTED LINES ---")

build_fitted_lines <- function(curve_fits_dt, curve_input_dt, n_grid = 300) {
  x <- unique(curve_fits_dt[, .(
    Protein_group, Protein_description, Genes, Replicate,
    a_vehicle, b_vehicle, plateau_vehicle,
    a_treatment, b_treatment, plateau_treatment
  )])
  
  temp_ranges <- curve_input_dt[, .(
    Temp_min = min(Temp),
    Temp_max = max(Temp)
  ), by = .(Protein_group, Replicate)]
  
  x <- merge(x, temp_ranges, by = c("Protein_group", "Replicate"), all.x = TRUE)
  
  fit_list <- vector("list", nrow(x) * 2)
  idx <- 1L
  
  for (i in seq_len(nrow(x))) {
    if (!is.finite(x$Temp_min[i]) || !is.finite(x$Temp_max[i])) next
    tt <- seq(x$Temp_min[i], x$Temp_max[i], length.out = n_grid)
    
    if (is.finite(x$a_vehicle[i]) && is.finite(x$b_vehicle[i]) &&
        is.finite(x$plateau_vehicle[i])) {
      fit_list[[idx]] <- data.table(
        Protein_group       = x$Protein_group[i],
        Protein_description = x$Protein_description[i],
        Genes               = x$Genes[i],
        Replicate           = x$Replicate[i],
        Condition           = "L",
        Temp                = tt,
        Fitted              = melt_fun(tt, x$a_vehicle[i], x$b_vehicle[i], x$plateau_vehicle[i])
      )
      idx <- idx + 1L
    }
    
    if (is.finite(x$a_treatment[i]) && is.finite(x$b_treatment[i]) &&
        is.finite(x$plateau_treatment[i])) {
      fit_list[[idx]] <- data.table(
        Protein_group       = x$Protein_group[i],
        Protein_description = x$Protein_description[i],
        Genes               = x$Genes[i],
        Replicate           = x$Replicate[i],
        Condition           = "H",
        Temp                = tt,
        Fitted              = melt_fun(tt, x$a_treatment[i], x$b_treatment[i], x$plateau_treatment[i])
      )
      idx <- idx + 1L
    }
  }
  
  if (idx == 1L) {
    return(data.table(
      Protein_group = character(), Protein_description = character(),
      Genes = character(), Replicate = character(),
      Condition = character(), Temp = numeric(), Fitted = numeric()
    ))
  }
  
  rbindlist(fit_list[seq_len(idx - 1L)], use.names = TRUE, fill = TRUE)
}

curve_lines <- build_fitted_lines(curve_fits, curve_input, n_grid = 300)
diag("Fitted line rows: ", nrow(curve_lines))
if (write_intermediates)
  fwrite(curve_lines, "protein_curve_fitted_lines.tsv.gz", sep = "\t")

#############################################################################################################
# SECTION 22. PLOTTING HELPERS
#############################################################################################################

plot_tpptr_protein <- function(protein_id,
                               curve_input_dt = curve_input,
                               curve_lines_dt = curve_lines,
                               summary_dt     = summary_hits,
                               save_file      = NULL,
                               width          = 18,
                               height         = 10) {
  
  if (is.na(protein_id) || is.null(protein_id)) {
    message("Skipping NA/NULL protein_id")
    return(list(plot = NULL, saved = FALSE))
  }
  
  pts <- copy(curve_input_dt[Protein_group == protein_id])
  ln  <- copy(curve_lines_dt[Protein_group == protein_id])
  sm  <- copy(summary_dt[Protein_group    == protein_id])
  
  if (nrow(pts) == 0) {
    message(sprintf("No curve input found for: %s", protein_id))
    return(list(plot = NULL, saved = FALSE))
  }
  
  gene_txt <- if (nrow(sm) > 0 && "Genes" %in% names(sm) &&
                  !is.na(sm$Genes[1])) sm$Genes[1] else ""
  desc_txt <- if (nrow(sm) > 0 && "Protein_description" %in% names(sm) &&
                  !is.na(sm$Protein_description[1])) sm$Protein_description[1] else ""
  
  if (nrow(sm) == 0) {
    title_txt    <- protein_id
    subtitle_txt <- "No summary stats available"
  } else {
    ci_txt <- if ("ci_quality" %in% names(sm) && !is.na(sm$ci_quality[1]))
      sprintf(" | CI=%s", sm$ci_quality[1]) else ""
    subtitle_txt <- sprintf(
      "%s | %s | mean \u0394Tm = %.2f\u00b0C | R1 = %.2f, R2 = %.2f | adjP = %.3g%s | %s",
      gene_txt, desc_txt,
      sm$mean_deltaTm[1],
      ifelse(is.finite(sm$deltaTm_R1[1]), sm$deltaTm_R1[1], NA),
      ifelse(is.finite(sm$deltaTm_R2[1]), sm$deltaTm_R2[1], NA),
      sm$best_p_adj[1],
      ci_txt,
      sm$hit_tier[1]
    )
    title_txt <- protein_id
  }
  
  pts[, Condition := factor(Condition, levels = c("L", "H"))]
  if (nrow(ln) > 0) ln[, Condition := factor(Condition, levels = c("L", "H"))]
  
  p <- ggplot() +
    geom_point(data = pts,
               aes(x = Temp, y = FC_fit, colour = Condition, shape = Condition),
               size = 2.2, alpha = 0.95)
  
  if (nrow(ln) > 0) {
    p <- p + geom_line(data = ln,
                       aes(x = Temp, y = Fitted, colour = Condition),
                       linewidth = 0.8)
  }
  
  p <- p +
    facet_wrap(~Replicate, nrow = 1) +
    scale_colour_manual(values = c(L = "#ca00a3ff", H = "#37abc8ff"),
                        labels = c(L = "DMSO (light)", H = "Drug (heavy)")) +
    scale_shape_manual(values = c(L = 16, H = 17),
                       labels = c(L = "DMSO (light)", H = "Drug (heavy)")) +
    coord_cartesian(ylim = c(0, 1.6)) +
    labs(title    = title_txt,
         subtitle = subtitle_txt,
         x        = "Temperature (\u00b0C)",
         y        = "Normalized soluble fraction") +
    theme_bw(base_size = 11) +
    theme(legend.position  = "top",
          strip.background = element_blank(),
          panel.grid.minor = element_blank())
  
  saved_ok <- FALSE
  if (!is.null(save_file)) {
    saved_ok <- tryCatch({
      ggsave(save_file, p, width = width, height = height,
             units = "cm", dpi = 300, bg = "white")
      TRUE
    }, error = function(e) {
      message(sprintf("Failed to save plot for %s: %s", protein_id, e$message))
      FALSE
    })
  }
  
  list(plot = p, saved = saved_ok)
}

#############################################################################################################
# SECTION 23. GLOBAL QC PLOTS
#
# BUG 6 FIX: pNAs y-axis upper limit is now guarded against empty data_completeness.
# max() on an empty numeric vector returns -Inf, which causes ggplot to error.
#############################################################################################################

diag("\n--- SECTION 23: QC PLOTS ---")

desired_temp_order      <- sort(unique(proteins$Temp_num))
desired_replicate_order <- sort(unique(proteins$Replicate))
temp_levels             <- as.character(desired_temp_order)
temp_color_mapping      <- setNames(
  colorRampPalette(c("#2166ac", "#d73027"))(length(temp_levels)),
  temp_levels
)

CV_function <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x) / mean(x) * 100
}

precursor_by_sample <- report[
  !is.na(Precursor.Quantity) & !is.na(label) & !is.na(Temp),
  .(qty = median(Precursor.Quantity, na.rm = TRUE)),
  by = .(Precursor.Id.nolabels, label, Temp, Sample)
]

precursor_CVs <- precursor_by_sample[
  , .(PrecursorCV = CV_function(qty), N_samples = .N),
  by = .(Precursor.Id.nolabels, label, Temp)
][N_samples >= 2 & is.finite(PrecursorCV)]

precursor_CVs[, precursor_rank := frank(PrecursorCV), by = .(Temp, label)]
precursor_CVs[, label     := factor(label, levels = c("L", "H"))]
precursor_CVs[, Temp_char := as.character(Temp)]

protein_long <- melt(
  proteins,
  id.vars      = c("Protein_group", "Temp_num", "Sample"),
  measure.vars = intersect(c("LFQ_L", "LFQ_H"), names(proteins)),
  variable.name = "label", value.name = "LFQ_intensity", na.rm = TRUE
)
protein_long[, label := gsub("LFQ_", "", label)]

protein_CVs <- protein_long[
  , .(ProteinCV = CV_function(LFQ_intensity), N_samples = uniqueN(Sample)),
  by = .(Protein_group, label, Temp_num)
][N_samples >= 2 & is.finite(ProteinCV)]

protein_CVs[, protein_rank := frank(ProteinCV), by = .(Temp_num, label)]
protein_CVs[, label     := factor(label, levels = c("L", "H"))]
protein_CVs[, Temp_char := as.character(Temp_num)]

proteins[, Temp_plot := factor(Temp_num, levels = desired_temp_order)]
proteins[, Replicate := factor(Replicate, levels = desired_replicate_order)]
report[,   Temp_plot := factor(Temp, levels = desired_temp_order)]
report[,   Replicate := factor(Replicate, levels = desired_replicate_order)]

sample_levels <- proteins[order(Temp_num, Replicate), unique(Sample)]
proteins[, Sample := factor(Sample, levels = sample_levels)]
report[,   Sample := factor(Sample, levels = sample_levels)]

base_theme <- theme(
  line                  = element_line(linewidth = 0.25),
  text                  = element_text(size = 7),
  rect                  = element_rect(colour = "black", linewidth = 0.25, fill = NA),
  panel.border          = element_rect(linewidth = 0.25),
  panel.background      = element_blank(),
  panel.grid            = element_blank(),
  axis.text             = element_text(colour = "black", size = 6),
  legend.background     = element_blank(),
  legend.box.background = element_blank(),
  legend.key.width      = unit(4, "mm"),
  legend.key            = element_rect(colour = NA),
  strip.background      = element_blank(),
  strip.text            = element_text(vjust = 1)
)
theme_set(base_theme)

HL_cols         <- c(H = "#37abc8ff", L = "#ca00a3ff",
                     heavy = "#37abc8ff", light = "#ca00a3ff", total = "black")
scale_HL_labels <- c(L = "DMSO (light)", H = "Drug (heavy)")

plot_int_DT <- rbindlist(list(
  report[!is.na(Precursor.Quantity),
         .(Sample, `log2 quantity` = log2(pmax(Precursor.Quantity, 1)),
           Type = "Precursor quantity")],
  proteins[!is.na(LFQ),
           .(Sample, `log2 quantity` = log2(pmax(LFQ, 1)),
             Type = "Protein MaxLFQ")],
  proteins[!is.na(Intensity),
           .(Sample, `log2 quantity` = log2(pmax(Intensity, 1)),
             Type = "Protein Intensity")]
), use.names = TRUE, fill = TRUE)

plot_int_DT[, Type := factor(Type,
                             levels = c("Precursor quantity", "Protein MaxLFQ", "Protein Intensity"))]

pInt <- ggplot(plot_int_DT, aes(x = Sample, y = `log2 quantity`)) +
  facet_wrap(~Type, scales = "free_x") +
  geom_boxplot(outliers = FALSE, linewidth = 0.25) +
  theme(axis.title.y = element_blank()) +
  coord_flip()

cum_cols <- intersect(c("LFQ_H", "LFQ_L", "Intensity_H", "Intensity_L"), names(proteins))
plot_int <- proteins[, lapply(.SD, function(x) sum(x, na.rm = TRUE)),
                     by = Sample, .SDcols = cum_cols]
plot_int <- melt(plot_int, id.vars = "Sample",
                 value.name = "Cumulative quantity", variable.name = "label")
plot_int[, Type  := fifelse(label %like% "LFQ", "Protein MaxLFQ", "Protein intensity")]
plot_int[, label := gsub(".+_", "", label)]
plot_int[, Type  := factor(Type, levels = c("Protein MaxLFQ", "Protein intensity"))]

pCumInt <- ggplot(plot_int, aes(x = Sample, y = `Cumulative quantity`, fill = label)) +
  facet_wrap(~Type, scales = "free_x") +
  geom_bar(stat = "identity") +
  ylab("Cumulative protein intensity") +
  scale_fill_manual(values = HL_cols,
                    guide  = guide_legend(title = "Channel"),
                    labels = scale_HL_labels) +
  theme(axis.title.y = element_blank()) +
  coord_flip()

count_cols <- intersect(c("Precursor.Id", "Stripped.Sequence", "Protein.Group"), names(report))
counts <- report[!is.na(label),
                 lapply(.SD, uniqueN), by = .(Sample, label), .SDcols = count_cols]
counts <- melt(counts, id.vars = c("Sample", "label"), value.name = "IDs")

pCount <- ggplot(counts, aes(x = Sample, y = IDs / 1000, fill = label,
                             label = format(IDs, big.mark = ",", scientific = FALSE))) +
  facet_wrap(~variable, scales = "free_x",
             labeller = as_labeller(c(Precursor.Id      = "Precursors",
                                      Stripped.Sequence = "Peptides",
                                      Protein.Group     = "Protein groups"))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(size = 2, hjust = 1.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = HL_cols,
                    guide  = guide_legend(title = "Channel"),
                    labels = scale_HL_labels) +
  ylab("Number of IDs [x1,000]") +
  theme(axis.title.y = element_blank()) +
  coord_flip()

f_prot_count <- function(dt, var_name) {
  tmp <- dt[, .(N_samples = .N), by = Protein_group][, .(N_proteins = .N), by = N_samples]
  tmp <- tmp[order(-N_samples)]
  tmp[, cumulative_protein_N := cumsum(N_proteins)]
  tmp[, Precursors := var_name]
  tmp
}

data_completeness <- rbindlist(list(
  f_prot_count(proteins, "all"),
  f_prot_count(proteins[!is.na(N_precursors_L) & N_precursors_L >= 2], "min 2 DMSO (L)"),
  f_prot_count(proteins[!is.na(N_precursors_H) & N_precursors_H >= 2], "min 2 Drug (H)"),
  f_prot_count(proteins[!is.na(N_precursors_L) & !is.na(N_precursors_H) &
                          N_precursors_L >= 2 & N_precursors_H >= 2],
               "min 2 DMSO &\nmin 2 Drug")
), use.names = TRUE, fill = TRUE)

# BUG 6 FIX: guard against empty data_completeness or all-NA cumulative_protein_N
pNAs_ylim_max <- if (nrow(data_completeness) > 0) {
  v <- max(data_completeness$cumulative_protein_N, na.rm = TRUE)
  if (is.finite(v) && v > 0) v / 1000 else 1
} else {
  1
}

pNAs <- ggplot(data_completeness,
               aes(x = N_samples, y = cumulative_protein_N / 1000, colour = Precursors)) +
  geom_point(size = 0.8) + geom_line(linewidth = 0.25) +
  xlab("Number of samples") + ylab("Number of proteins [x1,000]") +
  scale_x_continuous(breaks = seq(1, max(data_completeness$N_samples, 1), 1)) +
  scale_y_continuous(limits = c(0, pNAs_ylim_max)) +
  scale_color_brewer(palette = "Dark2") +
  theme(panel.grid.major = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
        panel.grid.minor = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"))

pCV1 <- ggplot(precursor_CVs,
               aes(x = precursor_rank / 1000, y = PrecursorCV, colour = Temp_char)) +
  facet_wrap(~label, labeller = as_labeller(c(L = "DMSO (light)", H = "Drug (heavy)"))) +
  geom_line(linewidth = 0.25) +
  labs(x = "Number of precursors [x1,000]", y = "Precursor CV [%]") +
  coord_cartesian(ylim = c(0, 50)) +
  scale_color_manual(values = temp_color_mapping, guide = guide_legend(title = "Temp (°C)")) +
  theme(panel.grid.major = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
        panel.grid.minor = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
        legend.position  = "right")

pCV2 <- ggplot(protein_CVs,
               aes(x = protein_rank / 1000, y = ProteinCV, colour = Temp_char)) +
  facet_wrap(~label, labeller = as_labeller(c(L = "DMSO (light)", H = "Drug (heavy)"))) +
  geom_line(linewidth = 0.25) +
  labs(x = "Number of proteins [x1,000]", y = "Protein CV [%]") +
  coord_cartesian(ylim = c(0, 50)) +
  scale_color_manual(values = temp_color_mapping, guide = guide_legend(title = "Temp (°C)")) +
  theme(panel.grid.major = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
        panel.grid.minor = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
        legend.position  = "right")

pca_plot <- function(dt, value_var, title_txt) {
  wide <- dcast(dt, Protein_group ~ Sample, value.var = value_var)
  wide <- wide[complete.cases(wide)]
  if (nrow(wide) < 3) {
    return(ggplot() +
             labs(title = title_txt, subtitle = "Too few complete cases for PCA") +
             theme_bw())
  }
  mat <- t(data.frame(wide[, -"Protein_group"], row.names = wide$Protein_group))
  pca <- prcomp(mat, scale. = TRUE)
  ve  <- summary(pca)$importance
  dtp <- data.table(pca$x, keep.rownames = "Sample")
  dtp <- merge(dtp, unique(dt[, .(Sample, Temp_plot, Replicate)]),
               by = "Sample", all.x = TRUE)
  dtp[, Temp_char := as.character(as.integer(as.character(Temp_plot)))]
  ggplot(dtp, aes(x = PC1, y = PC2, colour = Temp_char, shape = Replicate)) +
    ggtitle(title_txt) + geom_point(size = 2.5) +
    xlab(paste0("PC1 [", round(ve["Proportion of Variance", "PC1"] * 100, 0), "%]")) +
    ylab(paste0("PC2 [", round(ve["Proportion of Variance", "PC2"] * 100, 0), "%]")) +
    scale_color_manual(values = temp_color_mapping, guide = guide_legend(title = "Temp (°C)")) +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "right")
}

pPCA_L <- pca_plot(proteins, "LFQ_L", "DMSO control (Light)")
pPCA_H <- pca_plot(proteins, "LFQ_H", "Drug treated (Heavy)")

pHL <- ggplot(
  proteins[!is.na(log2_HL_ratio)],
  aes(x = Temp_plot, y = log2_HL_ratio,
      fill = as.character(as.integer(as.character(Temp_plot))))
) +
  geom_boxplot(outliers = FALSE, linewidth = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey40") +
  scale_fill_manual(values = temp_color_mapping, guide = "none") +
  labs(x     = "Temperature (°C)",
       y     = "log2(Drug/DMSO) ratio",
       title = "H/L ratio distribution per temperature") +
  theme(panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"))

pShift <- ggplot(
  summary_hits[is.finite(mean_deltaTm)],
  aes(x = mean_deltaTm, y = -log10(pmax(best_p_adj, 1e-300)), colour = TPPTR_hit)
) +
  geom_point(alpha = 0.8, size = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Mean \u0394Tm (H \u2212 L)",
       y = "-log10(best adjusted P)",
       title = "TPP-TR hit overview") +
  theme_bw(base_size = 10)

top_hits_bar <- summary_hits[
  TPPTR_hit == TRUE & is.finite(mean_deltaTm) & !is.na(hit_direction)
][order(-abs(mean_deltaTm))][1:min(.N, 30)]

pTop <- if (nrow(top_hits_bar) > 0) {
  ggplot(top_hits_bar,
         aes(x = reorder(Protein_group, mean_deltaTm),
             y = mean_deltaTm, fill = hit_direction)) +
    geom_col() + coord_flip() +
    labs(x = "Protein", y = "Mean \u0394Tm", title = "Top TPP-TR hits") +
    theme_bw(base_size = 10)
} else {
  ggplot() +
    labs(title = "Top TPP-TR hits", subtitle = "No significant hits found") +
    theme_bw(base_size = 10)
}

combined_plot <- (
  pInt / pCumInt / pCount /
    wrap_elements(pNAs + plot_spacer()) /
    wrap_elements(pCV1 + pCV2) /
    wrap_elements(pPCA_L + pPCA_H) /
    pHL / pShift / pTop
) +
  plot_layout(ncol = 1, heights = c(1, 1, 1.8, 1, 1, 1, 0.8, 0.8, 1.2)) +
  plot_annotation()

qc_file <- file.path(plot_dir, sprintf("QC_plot_%s.png", dataset_name))
ggsave(qc_file, combined_plot, width = 35, height = 115, units = "cm", dpi = 300, bg = "white")
diag("QC plot saved: ", qc_file)

#############################################################################################################
# SECTION 24. INDIVIDUAL HIT PLOTS
#############################################################################################################

diag("\n--- SECTION 24: INDIVIDUAL HIT PLOTS ---")

top_plot_hits <- summary_hits[
  TPPTR_hit == TRUE & !is.na(Protein_group) & is.finite(mean_deltaTm)
][order(-abs(mean_deltaTm))][1:min(.N, top_n_hit_plots), Protein_group]

diag("Proteins to plot: ", length(top_plot_hits))

if (length(top_plot_hits) == 0) {
  diag("No hits to plot.")
  diag("Check: Section 14 quality filter counts, Section 16 bootstrap summary entries,")
  diag("       Section 17 pass_statistical_hit, Section 18 both_exceed_variability.")
} else {
  n_saved <- 0L
  for (pid in top_plot_hits) {
    if (is.na(pid)) next
    safe_name <- gsub("[^A-Za-z0-9._-]", "_", pid)
    out_file  <- file.path(plot_dir,
                           paste0(safe_name, "_TPPTR_curve_", dataset_name, ".png"))
    res <- plot_tpptr_protein(protein_id = pid, save_file = out_file)
    if (isTRUE(res$saved)) n_saved <- n_saved + 1L
  }
  diag("Individual curve plots saved: ", n_saved, " / ", length(top_plot_hits))
}

#############################################################################################################
# SECTION 25. INTACT + EXTRACT MERGE HELPER
#############################################################################################################

merge_intact_extract_results <- function(
    intact_file  = "protein_tpptr_hits_summary_intact.tsv.gz",
    extract_file = "protein_tpptr_hits_summary_extract.tsv.gz",
    out_file     = "direct_indirect_classified.tsv.gz") {
  
  intact  <- fread(intact_file)
  extract <- fread(extract_file)
  
  merged <- merge(intact, extract, by = "Protein_group",
                  suffixes = c("_intact", "_extract"), all = TRUE)
  
  merged[, same_direction_extract := fifelse(
    is.finite(mean_deltaTm_intact) & is.finite(mean_deltaTm_extract),
    sign(mean_deltaTm_intact) == sign(mean_deltaTm_extract),
    NA
  )]
  
  merged[, classification := fifelse(
    TPPTR_hit_intact == TRUE & TPPTR_hit_extract == TRUE,
    "direct_like",
    fifelse(
      TPPTR_hit_intact == TRUE &
        (is.na(TPPTR_hit_extract) | TPPTR_hit_extract == FALSE),
      "indirect_like",
      fifelse(
        (is.na(TPPTR_hit_intact) | TPPTR_hit_intact == FALSE) &
          TPPTR_hit_extract == TRUE,
        "extract_only_ambiguous",
        "no_shift"
      )
    )
  )]
  
  fwrite(merged, out_file, sep = "\t")
  merged
}

#############################################################################################################
diag("\n=== PIPELINE COMPLETE ===")
diag("Run finished:    ", format(Sys.time()))
diag("Diagnostic log:  ", diag_log)
diag("Output directory:", plot_dir)
#############################################################################################################