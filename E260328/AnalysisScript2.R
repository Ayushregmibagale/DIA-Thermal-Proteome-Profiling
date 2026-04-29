#############################################################################################################
#
# CETSA-MS SILAC ANALYSIS SCRIPT  v2.1
#
# DESIGN
# ------
# Input : protein_quant.csv.gz from the updated DIA-NN SILAC protein summarisation script
#
# Channel meaning inherited from upstream script:
#   LFQ_L / label L  = light = DMSO / vehicle / control
#   LFQ_H / label H  = heavy = drug / treated
#
# Important upstream normalization already performed:
#   1. Precursor.Quantity was total-signal normalized per Run before MaxLFQ.
#      Therefore LFQ_L and LFQ_H are already loading-corrected.
#
#   2. Precursor.Translated.raw was kept unnormalized for SILAC H/L ratio calculation.
#      Protein-level log2_HL_ratio_raw was then run-median centred.
#
#   3. In the updated protein-quant script:
#          log2_HL_ratio = log2_HL_ratio_norm
#      Therefore this analysis script treats log2_HL_ratio as the already-normalized primary SILAC ratio.
#
# TWO COMPLEMENTARY ANALYSIS TRACKS
# ---------------------------------
# Track A - SILAC ratio track, PRIMARY
#           Uses log2_HL_ratio per temperature per replicate.
#           This is the cleanest CETSA-MS SILAC readout because H/L is measured within the same run.
#           No additional sample normalization is applied here.
#
# Track B - Absolute melt curve track, SECONDARY
#           Uses LFQ_L and LFQ_H.
#           These LFQ values were already run-normalized upstream.
#           This script does NOT apply another global sample-median normalization by default.
#           It only converts each protein/replicate/condition into fraction remaining:
#               fraction(T) = LFQ(T) / LFQ(lowest temperature)
#
# KEY v2.1 CHANGES
# ----------------
# 1. do_global_sample_median_normalization defaults to FALSE to avoid double-normalizing LFQ values.
# 2. flag_outlier_reps defaults to FALSE because the current design has only two replicates.
# 3. Track A is ranked first in the final hit table.
# 4. The script prefers log2_HL_ratio_norm if present; otherwise falls back to log2_HL_ratio.
# 5. Track B deltaTm is kept as secondary supportive evidence.
#
#############################################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(matrixStats)
  library(pheatmap)
})

#############################################################################################################
# USER SETTINGS
#############################################################################################################

input_file  <- "protein_quant.csv.gz"
output_dir  <- "CETSA_SILAC_output_v2_1"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Filtering ─────────────────────────────────────────────────────────────────
min_precursors_per_channel <- 2
min_ratios_for_HL          <- 2
min_temps_per_curve        <- 6

# ── Outlier replicate detection ───────────────────────────────────────────────
# With only 2 replicates, automatic outlier removal is usually too aggressive.
# Keep FALSE by default. If you later have >=3 biological replicates, TRUE is more reasonable.
flag_outlier_reps          <- FALSE
outlier_mzscore_threshold  <- 3.5

# ── Normalization, Track B only ───────────────────────────────────────────────
# IMPORTANT:
# The upstream protein-quant script already total-signal normalizes Precursor.Quantity before MaxLFQ.
# Therefore LFQ_L and LFQ_H are already loading-corrected.
# Leave this FALSE to avoid double-normalization and avoid flattening real CETSA temperature effects.
do_global_sample_median_normalization <- FALSE
use_log2_for_global_normalization     <- FALSE

# Fraction limits after lowest-temperature normalization.
# A cap avoids unstable fits caused by extreme high-temperature outliers.
max_fraction_cap   <- 1.5
min_fraction_floor <- 0.0

# ── Curve fitting ─────────────────────────────────────────────────────────────
fit_on_replicates       <- TRUE
fit_on_mean_curve       <- TRUE
n_fit_starts            <- 3
min_r2_for_fit          <- 0.7
require_monotonic_soft  <- FALSE

# ── Hit scoring ───────────────────────────────────────────────────────────────
# Track A is primary for SILAC CETSA-MS.
# Track B fitted deltaTm is secondary/supportive.
primary_rank_track <- "ratio"  # currently only "ratio" is recommended for this design

# ── Plotting ──────────────────────────────────────────────────────────────────
proteins_of_interest  <- c(
  # "HSP90AA1", "NAMPT"
)
top_n_deltaTm_plots   <- 24
plot_width  <- 10
plot_height <- 7
plot_dpi    <- 300

#############################################################################################################
# HELPERS
#############################################################################################################

msg <- function(...) cat(paste0(..., "\n"))

safe_mean <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else mean(x) }
safe_sd   <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x) }
safe_sem  <- function(x) { x <- x[is.finite(x)]; n <- length(x); if (n < 2) NA_real_ else sd(x) / sqrt(n) }
safe_n    <- function(x) sum(is.finite(x))

safe_auc <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

cv_percent <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2 || mean(x) == 0) return(NA_real_)
  100 * sd(x) / mean(x)
}

modified_zscore <- function(x) {
  med <- median(x, na.rm = TRUE)
  mad <- median(abs(x - med), na.rm = TRUE)
  if (!is.finite(mad) || mad == 0) mad <- mean(abs(x - med), na.rm = TRUE) / 0.6745
  if (!is.finite(mad) || mad == 0) return(rep(0, length(x)))
  0.6745 * (x - med) / mad
}

welch_t_pvalue <- function(x, mu0 = 0) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2) return(NA_real_)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  t_stat <- (mean(x) - mu0) / (s / sqrt(n))
  2 * pt(-abs(t_stat), df = n - 1)
}

#############################################################################################################
# SIGMOID FIT
#############################################################################################################

sigmoid_fun <- function(temp, top, bottom, Tm, slope) {
  bottom + (top - bottom) / (1 + exp((temp - Tm) / slope))
}

.try_nls <- function(x, y, starts) {
  warn_txt <- character()
  fit <- tryCatch(
    withCallingHandlers(
      nls(
        frac ~ bottom + (top - bottom) / (1 + exp((Temp - Tm) / slope)),
        data      = data.frame(Temp = x, frac = y),
        start     = starts,
        algorithm = "port",
        lower     = c(top = 0.5, bottom = -0.2, Tm = min(x), slope = 0.15),
        upper     = c(top = 1.5, bottom = 1.0, Tm = max(x), slope = 25),
        control   = nls.control(maxiter = 500, warnOnly = TRUE)
      ),
      warning = function(w) {
        warn_txt <<- c(warn_txt, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  list(fit = fit, warnings = warn_txt)
}

fit_sigmoid_curve <- function(dt) {
  dt <- dt[is.finite(Temp) & is.finite(frac)][order(Temp)]
  out <- list(
    success = FALSE,
    n_points = nrow(dt),
    top = NA_real_, bottom = NA_real_, Tm = NA_real_, slope = NA_real_,
    rss = NA_real_, r2 = NA_real_, auc = NA_real_,
    fit = NULL,
    warning_msg = NA_character_
  )
  
  if (nrow(dt) < min_temps_per_curve) return(out)
  if (uniqueN(dt$Temp) < min_temps_per_curve) return(out)
  if (sum(is.finite(dt$frac)) < min_temps_per_curve) return(out)
  if (all(abs(dt$frac - dt$frac[1]) < 1e-8, na.rm = TRUE)) return(out)
  
  x <- dt$Temp
  y <- pmin(pmax(dt$frac, min_fraction_floor), max_fraction_cap)
  
  if (require_monotonic_soft) {
    dy <- diff(y[order(x)])
    if (is.finite(mean(dy > 0, na.rm = TRUE)) && mean(dy > 0, na.rm = TRUE) > 0.50) return(out)
  }
  
  top_start    <- min(max(quantile(y, 0.9, na.rm = TRUE), 0.7), 1.3)
  bottom_start <- min(max(quantile(y, 0.1, na.rm = TRUE), 0.0), 0.6)
  slope_start  <- 2
  
  midpoint <- (top_start + bottom_start) / 2
  Tm_data  <- x[which.min(abs(y - midpoint))]
  if (!is.finite(Tm_data)) Tm_data <- median(x, na.rm = TRUE)
  
  all_candidates <- unique(as.numeric(c(
    Tm_data,
    quantile(x, 0.35, na.rm = TRUE),
    quantile(x, 0.65, na.rm = TRUE),
    median(x, na.rm = TRUE)
  )))
  Tm_candidates <- head(all_candidates[is.finite(all_candidates)], n_fit_starts)
  
  best_rss <- Inf
  best_res <- NULL
  all_warnings <- character()
  
  for (Tm_s in Tm_candidates) {
    starts <- list(top = top_start, bottom = bottom_start, Tm = Tm_s, slope = slope_start)
    res <- .try_nls(x, y, starts)
    all_warnings <- c(all_warnings, res$warnings)
    if (inherits(res$fit, "error")) next
    
    coefs <- coef(res$fit)
    if (!all(is.finite(coefs[c("top", "bottom", "Tm", "slope")]))) next
    
    pred <- tryCatch(predict(res$fit), error = function(e) rep(NA_real_, length(y)))
    if (!all(is.finite(pred))) next
    
    rss_i <- sum((y - pred)^2, na.rm = TRUE)
    if (rss_i < best_rss) {
      best_rss <- rss_i
      best_res <- list(fit = res$fit, coefs = coefs, pred = pred)
    }
  }
  
  out$warning_msg <- paste(unique(all_warnings), collapse = " | ")
  if (is.null(best_res)) return(out)
  
  coefs <- best_res$coefs
  pred  <- best_res$pred
  tss   <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  r2    <- if (tss > 0) 1 - best_rss / tss else NA_real_
  
  if (!is.finite(r2)) return(out)
  if (r2 < min_r2_for_fit) return(out)
  if (coefs["bottom"] >= coefs["top"]) return(out)
  
  out$success <- TRUE
  out$top     <- unname(coefs["top"])
  out$bottom  <- unname(coefs["bottom"])
  out$Tm      <- unname(coefs["Tm"])
  out$slope   <- unname(coefs["slope"])
  out$rss     <- best_rss
  out$r2      <- r2
  out$auc     <- safe_auc(x, pred)
  out$fit     <- best_res$fit
  out
}

predict_fit_grid <- function(fit_obj, temp_grid) {
  if (is.null(fit_obj) || !isTRUE(fit_obj$success) || is.null(fit_obj$fit)) return(NULL)
  pred <- tryCatch(
    predict(fit_obj$fit, newdata = data.frame(Temp = temp_grid)),
    error = function(e) NULL
  )
  if (is.null(pred)) return(NULL)
  data.table(Temp = temp_grid, frac_pred = pred)
}

#############################################################################################################
# READ & VALIDATE INPUT
#############################################################################################################

msg("Reading input: ", input_file)
proteins_raw <- fread(input_file)

required_cols <- c(
  "Run", "Sample", "Temp", "Replicate", "Protein_group", "Genes",
  "LFQ_L", "LFQ_H", "N_precursors_L", "N_precursors_H"
)
missing <- setdiff(required_cols, names(proteins_raw))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))

# Prefer the explicit normalized column from the updated protein quant script.
if ("log2_HL_ratio_norm" %in% names(proteins_raw)) {
  proteins_raw[, log2_HL_ratio_for_analysis := log2_HL_ratio_norm]
  msg("Using log2_HL_ratio_norm as Track A input.")
} else if ("log2_HL_ratio" %in% names(proteins_raw)) {
  proteins_raw[, log2_HL_ratio_for_analysis := log2_HL_ratio]
  msg("Using log2_HL_ratio as Track A input. Confirm upstream script already median-centred it.")
} else {
  proteins_raw[, log2_HL_ratio_for_analysis := NA_real_]
  msg("NOTE: no log2_HL_ratio or log2_HL_ratio_norm column found; Track A will be skipped.")
}

has_HL_ratio <- any(is.finite(proteins_raw$log2_HL_ratio_for_analysis))
has_N_ratios <- "N_ratios" %in% names(proteins_raw)

proteins_raw[, Temp := as.numeric(as.character(Temp))]
proteins_raw[, Replicate := as.character(Replicate)]
proteins_raw[, Genes := fifelse(is.na(Genes) | Genes == "", Protein_group, Genes)]
setorder(proteins_raw, Protein_group, Replicate, Temp)

msg("Rows in input       : ", nrow(proteins_raw))
msg("Unique proteins     : ", uniqueN(proteins_raw$Protein_group))
msg("Temperatures        : ", paste(sort(unique(proteins_raw$Temp)), collapse = ", "))
msg("Replicates          : ", paste(sort(unique(proteins_raw$Replicate)), collapse = ", "))
msg("Global LFQ sample median normalization in this script: ", do_global_sample_median_normalization)
msg("Outlier replicate exclusion in this script           : ", flag_outlier_reps)

#############################################################################################################
# PRECURSOR-COUNT FILTER
#############################################################################################################

proteins_filt <- copy(proteins_raw)

proteins_filt[
  (is.na(N_precursors_L) | N_precursors_L < min_precursors_per_channel) |
    (is.na(N_precursors_H) | N_precursors_H < min_precursors_per_channel),
  c("LFQ_L", "LFQ_H") := .(NA_real_, NA_real_)
]

if (has_HL_ratio && has_N_ratios) {
  proteins_filt[
    is.na(N_ratios) | N_ratios < min_ratios_for_HL,
    log2_HL_ratio_for_analysis := NA_real_
  ]
}

#############################################################################################################
# TRACK A: SILAC RATIO TRACK, PRIMARY
#############################################################################################################

if (has_HL_ratio) {
  msg("\n── Track A: SILAC ratio track, primary ───────────────────────────────")
  
  ratio_long <- proteins_filt[
    is.finite(log2_HL_ratio_for_analysis),
    .(
      Protein_group, Genes, Run, Sample, Temp, Replicate,
      log2_HL_ratio = log2_HL_ratio_for_analysis
    )
  ]
  
  ratio_mean <- ratio_long[, .(
    mean_log2_HL = safe_mean(log2_HL_ratio),
    sd_log2_HL   = safe_sd(log2_HL_ratio),
    sem_log2_HL  = safe_sem(log2_HL_ratio),
    n_rep        = safe_n(log2_HL_ratio)
  ), by = .(Protein_group, Genes, Temp)]
  
  ratio_temp_summary <- ratio_long[, .(
    mean_log2_ratio = safe_mean(log2_HL_ratio),
    sd_log2_ratio   = safe_sd(log2_HL_ratio),
    sem_log2_ratio  = safe_sem(log2_HL_ratio),
    n_rep           = safe_n(log2_HL_ratio)
  ), by = .(Protein_group, Genes, Temp)]
  
  ratio_peak_summary <- ratio_temp_summary[, {
    ok <- is.finite(mean_log2_ratio)
    if (!any(ok)) {
      .(
        max_abs_mean_log2_HL = NA_real_,
        direction_ratio_track = NA_character_,
        signed_peak_log2_HL = NA_real_,
        temp_of_max_abs_log2HL = NA_real_,
        auc_log2_HL_curve = NA_real_,
        p_value_ratio_track = NA_real_
      )
    } else {
      local_idx <- which.max(abs(mean_log2_ratio[ok]))
      signed_peak <- mean_log2_ratio[ok][local_idx]
      peak_temp <- Temp[ok][local_idx]
      .(
        max_abs_mean_log2_HL = abs(signed_peak),
        direction_ratio_track = fifelse(signed_peak > 0, "stabilization", "destabilization"),
        signed_peak_log2_HL = signed_peak,
        temp_of_max_abs_log2HL = peak_temp,
        auc_log2_HL_curve = safe_auc(Temp[ok], mean_log2_ratio[ok]),
        p_value_ratio_track = NA_real_
      )
    }
  }, by = .(Protein_group, Genes)]
  
  ratio_ttest <- ratio_long[
    ratio_peak_summary[, .(Protein_group, temp_of_max_abs_log2HL)],
    on = "Protein_group"
  ][
    is.finite(temp_of_max_abs_log2HL) & Temp == temp_of_max_abs_log2HL,
    .(p_val = welch_t_pvalue(log2_HL_ratio, mu0 = 0)),
    by = Protein_group
  ]
  ratio_peak_summary[ratio_ttest, p_value_ratio_track := i.p_val, on = "Protein_group"]
  
  msg("Proteins with >=1 ratio measurement: ", uniqueN(ratio_long$Protein_group))
  
  fwrite(ratio_long,         file.path(output_dir, "long_ratio_table.tsv.gz"), sep = "\t")
  fwrite(ratio_mean,         file.path(output_dir, "mean_ratio_curves.tsv.gz"), sep = "\t")
  fwrite(ratio_temp_summary, file.path(output_dir, "temp_ratio_summary.tsv.gz"), sep = "\t")
} else {
  ratio_long         <- data.table()
  ratio_mean         <- data.table()
  ratio_temp_summary <- data.table()
  ratio_peak_summary <- data.table()
}

#############################################################################################################
# TRACK B: ABSOLUTE MELT CURVE TRACK, SECONDARY
#############################################################################################################

msg("\n── Track B: absolute melt curve track, secondary ─────────────────────")

long <- rbindlist(list(
  proteins_filt[, .(
    Protein_group, Genes, Run, Sample, Temp, Replicate,
    channel = "L", condition = "DMSO",
    abundance = LFQ_L, N_precursors = N_precursors_L
  )],
  proteins_filt[, .(
    Protein_group, Genes, Run, Sample, Temp, Replicate,
    channel = "H", condition = "Drug",
    abundance = LFQ_H, N_precursors = N_precursors_H
  )]
), use.names = TRUE)

long[abundance <= 0, abundance := NA_real_]
long[, abundance_norm := abundance]

if (do_global_sample_median_normalization) {
  msg("Applying optional global sample median normalization, per channel.")
  msg("WARNING: this is usually not recommended with the updated upstream protein-quant script.")
  
  if (use_log2_for_global_normalization) {
    medtab <- long[is.finite(abundance), .(
      sample_median = median(log2(abundance), na.rm = TRUE)
    ), by = .(channel, Sample)]
    medtab[, global_median := median(sample_median, na.rm = TRUE), by = channel]
    medtab[, shift := global_median - sample_median]
    long <- merge(long, medtab[, .(channel, Sample, shift)], by = c("channel", "Sample"), all.x = TRUE)
    long[, abundance_norm := 2^(log2(abundance) + shift)]
    long[, shift := NULL]
  } else {
    medtab <- long[is.finite(abundance), .(
      sample_median = median(abundance, na.rm = TRUE)
    ), by = .(channel, Sample)]
    medtab[, global_median := median(sample_median, na.rm = TRUE), by = channel]
    medtab[, factor_norm := global_median / sample_median]
    long <- merge(long, medtab[, .(channel, Sample, factor_norm)], by = c("channel", "Sample"), all.x = TRUE)
    long[, abundance_norm := abundance * factor_norm]
    long[, factor_norm := NULL]
  }
} else {
  msg("Skipping global sample median normalization in Track B to avoid double-normalization.")
}

# Reference value: lowest-temperature value for each protein x replicate x condition.
# This is the correct CETSA melt-curve normalization.
long[, temp_min := suppressWarnings(min(Temp[is.finite(abundance_norm)], na.rm = TRUE)),
     by = .(Protein_group, Replicate, condition)]
long[!is.finite(temp_min), temp_min := NA_real_]

long[, ref_value := median(abundance_norm[Temp == temp_min], na.rm = TRUE),
     by = .(Protein_group, Replicate, condition)]
long[!is.finite(ref_value) | ref_value <= 0, ref_value := NA_real_]

if (flag_outlier_reps) {
  ref_dt <- unique(long[is.finite(ref_value), .(Protein_group, condition, Replicate, ref_value)])
  ref_dt[, mzscore := modified_zscore(ref_value), by = .(Protein_group, condition)]
  ref_dt[, outlier_rep := abs(mzscore) > outlier_mzscore_threshold]
  
  long <- merge(
    long,
    ref_dt[, .(Protein_group, condition, Replicate, mzscore_ref = mzscore, outlier_rep)],
    by = c("Protein_group", "condition", "Replicate"),
    all.x = TRUE
  )
  long[is.na(outlier_rep), outlier_rep := FALSE]
  
  n_outlier_reps <- sum(long[, .(is_out = any(outlier_rep)), by = .(Protein_group, condition, Replicate)]$is_out)
  msg("Outlier replicate flags raised and excluded: ", n_outlier_reps)
  long[outlier_rep == TRUE, ref_value := NA_real_]
} else {
  long[, mzscore_ref := NA_real_]
  long[, outlier_rep := FALSE]
  msg("Outlier replicate detection disabled; no Track B replicates excluded automatically.")
}

long[, frac := abundance_norm / ref_value]
long[, frac := pmin(pmax(frac, min_fraction_floor), max_fraction_cap)]

curve_stats <- long[, .(n_temps = sum(is.finite(frac))), by = .(Protein_group, Replicate, condition)]
keep_curves <- curve_stats[n_temps >= min_temps_per_curve]
long <- merge(
  long,
  keep_curves[, .(Protein_group, Replicate, condition)],
  by = c("Protein_group", "Replicate", "condition")
)

msg("Protein x condition x replicate curves passing filter: ", nrow(keep_curves))

curve_mean <- long[, .(
  mean_frac = safe_mean(frac),
  sd_frac   = safe_sd(frac),
  sem_frac  = safe_sem(frac),
  n_rep     = safe_n(frac)
), by = .(Protein_group, Genes, condition, Temp)]

rep_wide <- dcast(long, Protein_group + Genes + Replicate + Temp ~ condition, value.var = "frac")
if (!("DMSO" %in% names(rep_wide))) rep_wide[, DMSO := NA_real_]
if (!("Drug" %in% names(rep_wide))) rep_wide[, Drug := NA_real_]

rep_wide[, log2_ratio_Drug_vs_DMSO := fifelse(
  is.finite(Drug) & is.finite(DMSO) & Drug > 0 & DMSO > 0,
  log2(Drug / DMSO), NA_real_
)]
rep_wide[, diff_Drug_minus_DMSO := Drug - DMSO]

mean_wide <- dcast(curve_mean, Protein_group + Genes + Temp ~ condition, value.var = "mean_frac")
if (!("DMSO" %in% names(mean_wide))) mean_wide[, DMSO := NA_real_]
if (!("Drug" %in% names(mean_wide))) mean_wide[, Drug := NA_real_]

mean_wide[, log2_ratio_Drug_vs_DMSO := fifelse(
  is.finite(Drug) & is.finite(DMSO) & Drug > 0 & DMSO > 0,
  log2(Drug / DMSO), NA_real_
)]
mean_wide[, diff_Drug_minus_DMSO := Drug - DMSO]

temp_effect_summary <- rep_wide[, .(
  mean_log2_ratio = safe_mean(log2_ratio_Drug_vs_DMSO),
  sd_log2_ratio   = safe_sd(log2_ratio_Drug_vs_DMSO),
  sem_log2_ratio  = safe_sem(log2_ratio_Drug_vs_DMSO),
  mean_diff       = safe_mean(diff_Drug_minus_DMSO),
  sd_diff         = safe_sd(diff_Drug_minus_DMSO),
  n_rep           = safe_n(log2_ratio_Drug_vs_DMSO)
), by = .(Protein_group, Genes, Temp)]

#############################################################################################################
# SIGMOID CURVE FITTING
#############################################################################################################

temp_grid <- seq(min(long$Temp, na.rm = TRUE), max(long$Temp, na.rm = TRUE), length.out = 200)

run_fits <- function(input_dt, group_cols) {
  msg("  Fitting ", paste(group_cols, collapse = "/"), " curves...")
  groups <- split(input_dt, by = group_cols, keep.by = TRUE)
  
  results <- rbindlist(lapply(groups, function(dtg) {
    fr <- fit_sigmoid_curve(dtg)
    info <- unique(dtg[, ..group_cols])
    cbind(info, data.table(
      fit_success = fr$success,
      n_points = fr$n_points,
      top = fr$top,
      bottom = fr$bottom,
      Tm = fr$Tm,
      slope = fr$slope,
      rss = fr$rss,
      r2 = fr$r2,
      auc_fit = fr$auc,
      fit_warning = fr$warning_msg
    ))
  }), fill = TRUE)
  
  preds <- rbindlist(lapply(groups, function(dtg) {
    fr <- fit_sigmoid_curve(dtg)
    pred <- predict_fit_grid(fr, temp_grid)
    if (is.null(pred)) return(NULL)
    info <- unique(dtg[, ..group_cols])
    cbind(info, pred)
  }), fill = TRUE)
  
  list(results = results, preds = preds)
}

fit_results_rep  <- data.table()
fit_preds_rep    <- data.table()
fit_results_mean <- data.table()
fit_preds_mean   <- data.table()

if (fit_on_replicates && nrow(long) > 0) {
  rep_input <- long[, .(Protein_group, Genes, Replicate, condition, Temp, frac)]
  r <- run_fits(rep_input, c("Protein_group", "Genes", "Replicate", "condition"))
  fit_results_rep <- r$results
  fit_preds_rep   <- r$preds
}

if (fit_on_mean_curve && nrow(curve_mean) > 0) {
  mean_input <- curve_mean[, .(Protein_group, Genes, condition, Temp, frac = mean_frac)]
  r <- run_fits(mean_input, c("Protein_group", "Genes", "condition"))
  fit_results_mean <- r$results
  fit_preds_mean   <- r$preds
}

#############################################################################################################
# deltaTm / deltaAUC SUMMARIES, TRACK B SUPPORTIVE
#############################################################################################################

deltaTm_rep  <- data.table()
deltaTm_mean <- data.table()
replicate_summary <- data.table()

if (nrow(fit_results_rep) > 0 && any(fit_results_rep$fit_success)) {
  deltaTm_rep <- dcast(
    fit_results_rep[fit_success == TRUE],
    Protein_group + Genes + Replicate ~ condition,
    value.var = c("Tm", "auc_fit", "r2", "n_points")
  )
  
  for (cn in c("Tm_Drug", "Tm_DMSO", "auc_fit_Drug", "auc_fit_DMSO", "r2_Drug", "r2_DMSO")) {
    if (!(cn %in% names(deltaTm_rep))) deltaTm_rep[, (cn) := NA_real_]
  }
  
  deltaTm_rep[, deltaTm := Tm_Drug - Tm_DMSO]
  deltaTm_rep[, deltaAUC := auc_fit_Drug - auc_fit_DMSO]
  
  replicate_summary <- deltaTm_rep[, .(
    n_replicates   = .N,
    n_deltaTm      = sum(is.finite(deltaTm)),
    deltaTm_mean   = safe_mean(deltaTm),
    deltaTm_sd     = safe_sd(deltaTm),
    deltaTm_sem    = safe_sem(deltaTm),
    deltaTm_p_value = welch_t_pvalue(deltaTm, mu0 = 0),
    deltaAUC_mean  = safe_mean(deltaAUC),
    deltaAUC_sd    = safe_sd(deltaAUC),
    Tm_DMSO_mean   = safe_mean(Tm_DMSO),
    Tm_Drug_mean   = safe_mean(Tm_Drug),
    r2_DMSO_mean   = safe_mean(r2_DMSO),
    r2_Drug_mean   = safe_mean(r2_Drug)
  ), by = .(Protein_group, Genes)]
}

if (nrow(fit_results_mean) > 0 && any(fit_results_mean$fit_success)) {
  deltaTm_mean <- dcast(
    fit_results_mean[fit_success == TRUE],
    Protein_group + Genes ~ condition,
    value.var = c("Tm", "auc_fit", "r2", "n_points")
  )
  
  for (cn in c("Tm_Drug", "Tm_DMSO", "auc_fit_Drug", "auc_fit_DMSO")) {
    if (!(cn %in% names(deltaTm_mean))) deltaTm_mean[, (cn) := NA_real_]
  }
  
  deltaTm_mean[, deltaTm_mean_curve := Tm_Drug - Tm_DMSO]
  deltaTm_mean[, deltaAUC_mean_curve := auc_fit_Drug - auc_fit_DMSO]
}

# #############################################################################################################
# # FINAL HIT TABLE
# #############################################################################################################
# 
# msg("\nBuilding hit table...")
# 
# hit_table <- unique(proteins_raw[, .(Protein_group, Genes)])
# 
# if (nrow(replicate_summary) > 0) {
#   hit_table <- merge(hit_table, replicate_summary, by = c("Protein_group", "Genes"), all.x = TRUE)
# }
# if (nrow(deltaTm_mean) > 0) {
#   hit_table <- merge(
#     hit_table,
#     deltaTm_mean[, .(Protein_group, Genes, deltaTm_mean_curve, deltaAUC_mean_curve)],
#     by = c("Protein_group", "Genes"), all.x = TRUE
#   )
# }
# if (nrow(ratio_peak_summary) > 0) {
#   hit_table <- merge(hit_table, ratio_peak_summary, by = c("Protein_group", "Genes"), all.x = TRUE)
# }
# 
# curve_completeness <- long[, .(
#   n_points_total       = sum(is.finite(frac)),
#   n_replicates_present = uniqueN(Replicate[is.finite(frac)]),
#   n_conditions_present = uniqueN(condition[is.finite(frac)])
# ), by = .(Protein_group, Genes)]
# hit_table <- merge(hit_table, curve_completeness, by = c("Protein_group", "Genes"), all.x = TRUE)
# 
# peak_from_trackB <- temp_effect_summary[, {
#   ok <- is.finite(mean_log2_ratio)
#   if (!any(ok)) {
#     .(max_abs_mean_log2_ratio_trackB = NA_real_, signed_peak_log2_ratio_trackB = NA_real_)
#   } else {
#     val <- mean_log2_ratio[ok][which.max(abs(mean_log2_ratio[ok]))]
#     .(max_abs_mean_log2_ratio_trackB = abs(val), signed_peak_log2_ratio_trackB = val)
#   }
# }, by = .(Protein_group, Genes)]
# hit_table <- merge(hit_table, peak_from_trackB, by = c("Protein_group", "Genes"), all.x = TRUE)
# 
# for (cn in c(
#   "max_abs_mean_log2_HL", "signed_peak_log2_HL", "p_value_ratio_track",
#   "deltaTm_mean", "deltaTm_mean_curve", "max_abs_mean_log2_ratio_trackB",
#   "deltaTm_p_value"
# )) {
#   if (!(cn %in% names(hit_table))) hit_table[, (cn) := NA_real_]
# }
# if (!("direction_ratio_track" %in% names(hit_table))) hit_table[, direction_ratio_track := NA_character_]
# 
# # Track A primary ranking.
# # rank_score is intentionally the SILAC ratio effect size first.
# # deltaTm values are supportive but not the first sorting criterion for this SILAC design.
# hit_table[, rank_score := fcase(
#   is.finite(max_abs_mean_log2_HL), max_abs_mean_log2_HL,
#   is.finite(deltaTm_mean), abs(deltaTm_mean),
#   is.finite(deltaTm_mean_curve), abs(deltaTm_mean_curve),
#   is.finite(max_abs_mean_log2_ratio_trackB), max_abs_mean_log2_ratio_trackB,
#   default = NA_real_
# )]
# 
# # Best p-value is reported, but do not treat it as a strict FDR correction.
# p_cols <- intersect(c("p_value_ratio_track", "deltaTm_p_value"), names(hit_table))
# if (length(p_cols) > 0) {
#   hit_table[, p_value_best := do.call(pmin, c(.SD, na.rm = TRUE)), .SDcols = p_cols]
#   hit_table[!is.finite(p_value_best), p_value_best := NA_real_]
# } else {
#   hit_table[, p_value_best := NA_real_]
# }
# 
# # Convenience interpretation columns.
# hit_table[, evidence_track := fcase(
#   is.finite(max_abs_mean_log2_HL), "Track A SILAC ratio",
#   is.finite(deltaTm_mean), "Track B replicate deltaTm",
#   is.finite(deltaTm_mean_curve), "Track B mean-curve deltaTm",
#   is.finite(max_abs_mean_log2_ratio_trackB), "Track B fraction ratio",
#   default = "insufficient evidence"
# )]
# 
# hit_table[, stabilization_direction := fcase(
#   is.finite(signed_peak_log2_HL) & signed_peak_log2_HL > 0, "stabilization",
#   is.finite(signed_peak_log2_HL) & signed_peak_log2_HL < 0, "destabilization",
#   is.finite(deltaTm_mean) & deltaTm_mean > 0, "stabilization",
#   is.finite(deltaTm_mean) & deltaTm_mean < 0, "destabilization",
#   is.finite(deltaTm_mean_curve) & deltaTm_mean_curve > 0, "stabilization",
#   is.finite(deltaTm_mean_curve) & deltaTm_mean_curve < 0, "destabilization",
#   default = NA_character_
# )]
# 
# # Sort by Track A first, then p-value, then supportive Track B shift.
# hit_table[, .has_ratio := is.finite(max_abs_mean_log2_HL)]
# hit_table[, .ratio_score := fifelse(is.finite(max_abs_mean_log2_HL), max_abs_mean_log2_HL, -Inf)]
# hit_table[, .p_ord := fifelse(is.finite(p_value_best), p_value_best, Inf)]
# hit_table[, .trackB_score := pmax(
#   fifelse(is.finite(deltaTm_mean), abs(deltaTm_mean), -Inf),
#   fifelse(is.finite(deltaTm_mean_curve), abs(deltaTm_mean_curve), -Inf),
#   fifelse(is.finite(max_abs_mean_log2_ratio_trackB), max_abs_mean_log2_ratio_trackB, -Inf),
#   na.rm = TRUE
# )]
# hit_table[, .nrep := fifelse(is.finite(n_replicates_present), n_replicates_present, -Inf)]
# setorder(hit_table, -.has_ratio, -.ratio_score, .p_ord, -.trackB_score, -.nrep)
# hit_table[, c(".has_ratio", ".ratio_score", ".p_ord", ".trackB_score", ".nrep") := NULL]
# 
# msg("Proteins in hit table                     : ", nrow(hit_table))
# msg("With SILAC ratio score, Track A           : ", hit_table[is.finite(max_abs_mean_log2_HL), .N])
# msg("With replicate deltaTm, Track B           : ", hit_table[is.finite(deltaTm_mean), .N])
# msg("With mean-curve deltaTm, Track B          : ", hit_table[is.finite(deltaTm_mean_curve), .N])

#############################################################################################################
# FINAL HIT TABLE RANKING
# SILAC CETSA-MS CORRECT PRIORITY:
#   1. Track A: max_abs_mean_log2_HL
#   2. Track B: replicate-level deltaTm_mean
#   3. Track B: mean-curve deltaTm_mean_curve
#############################################################################################################

# Ensure all ranking columns exist
for (cn in c(
  "max_abs_mean_log2_HL",
  "signed_peak_log2_HL",
  "p_value_ratio_track",
  "deltaTm_mean",
  "deltaTm_mean_curve",
  "deltaTm_p_value",
  "max_abs_mean_log2_ratio_trackB"
)) {
  if (!(cn %in% names(hit_table))) hit_table[, (cn) := NA_real_]
}

if (!("direction_ratio_track" %in% names(hit_table))) {
  hit_table[, direction_ratio_track := NA_character_]
}

# Primary rank score:
# Track A first, then Track B only as fallback/support.
hit_table[, rank_score := fcase(
  is.finite(max_abs_mean_log2_HL),       max_abs_mean_log2_HL,
  is.finite(deltaTm_mean),               abs(deltaTm_mean),
  is.finite(deltaTm_mean_curve),         abs(deltaTm_mean_curve),
  is.finite(max_abs_mean_log2_ratio_trackB), max_abs_mean_log2_ratio_trackB,
  default = NA_real_
)]

# Explicitly label what evidence generated the rank score
hit_table[, primary_evidence := fcase(
  is.finite(max_abs_mean_log2_HL),       "Track A: SILAC H/L ratio",
  is.finite(deltaTm_mean),               "Track B: replicate deltaTm",
  is.finite(deltaTm_mean_curve),         "Track B: mean-curve deltaTm",
  is.finite(max_abs_mean_log2_ratio_trackB), "Track B: fraction ratio",
  default = "insufficient evidence"
)]

# Direction call.
# Prefer Track A sign because SILAC H/L is the primary evidence.
hit_table[, stabilization_direction := fcase(
  is.finite(signed_peak_log2_HL) & signed_peak_log2_HL > 0, "stabilization",
  is.finite(signed_peak_log2_HL) & signed_peak_log2_HL < 0, "destabilization",
  is.finite(deltaTm_mean) & deltaTm_mean > 0, "stabilization",
  is.finite(deltaTm_mean) & deltaTm_mean < 0, "destabilization",
  is.finite(deltaTm_mean_curve) & deltaTm_mean_curve > 0, "stabilization",
  is.finite(deltaTm_mean_curve) & deltaTm_mean_curve < 0, "destabilization",
  default = NA_character_
)]

# Best available p-value, but keep source visible.
hit_table[, p_value_best := NA_real_]
hit_table[, p_value_source := NA_character_]

hit_table[
  is.finite(p_value_ratio_track),
  `:=`(
    p_value_best = p_value_ratio_track,
    p_value_source = "Track A: SILAC ratio"
  )
]

hit_table[
  !is.finite(p_value_best) & is.finite(deltaTm_p_value),
  `:=`(
    p_value_best = deltaTm_p_value,
    p_value_source = "Track B: deltaTm"
  )
]

# Secondary p-value option:
# If both exist, store the smaller one separately but do not use it as the primary source.
hit_table[, p_value_min_available := pmin(
  fifelse(is.finite(p_value_ratio_track), p_value_ratio_track, Inf),
  fifelse(is.finite(deltaTm_p_value), deltaTm_p_value, Inf),
  na.rm = TRUE
)]
hit_table[!is.finite(p_value_min_available), p_value_min_available := NA_real_]

# Sort:
# 1. Proteins with Track A evidence first
# 2. Larger Track A effect first
# 3. Better Track A p-value
# 4. Larger replicate deltaTm
# 5. Larger mean-curve deltaTm
hit_table[, sort_has_trackA := is.finite(max_abs_mean_log2_HL)]
hit_table[, sort_trackA_effect := fifelse(is.finite(max_abs_mean_log2_HL), max_abs_mean_log2_HL, -Inf)]
hit_table[, sort_trackA_p := fifelse(is.finite(p_value_ratio_track), p_value_ratio_track, Inf)]
hit_table[, sort_deltaTm_rep := fifelse(is.finite(deltaTm_mean), abs(deltaTm_mean), -Inf)]
hit_table[, sort_deltaTm_mean := fifelse(is.finite(deltaTm_mean_curve), abs(deltaTm_mean_curve), -Inf)]

setorder(
  hit_table,
  -sort_has_trackA,
  -sort_trackA_effect,
  sort_trackA_p,
  -sort_deltaTm_rep,
  -sort_deltaTm_mean
)

hit_table[, c(
  "sort_has_trackA",
  "sort_trackA_effect",
  "sort_trackA_p",
  "sort_deltaTm_rep",
  "sort_deltaTm_mean"
) := NULL]


#############################################################################################################
# QC PLOTS
#############################################################################################################

msg("\nGenerating QC plots...")

p_qc_box <- ggplot(
  long[is.finite(abundance_norm)],
  aes(
    x = interaction(Temp, Replicate, lex.order = TRUE),
    y = log2(abundance_norm), fill = condition
  )
) +
  geom_boxplot(outlier.size = 0.2) +
  theme_bw(base_size = 10) +
  labs(
    x = "Temperature_Replicate",
    y = "log2 abundance, upstream-normalized only",
    title = "Protein abundance distributions, no second global normalization"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p_qc_frac <- ggplot(
  long[is.finite(frac)],
  aes(x = factor(Temp), y = frac, fill = condition)
) +
  geom_boxplot(outlier.size = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  theme_bw(base_size = 10) +
  labs(x = "Temperature", y = "Fraction vs lowest temperature", title = "Track B normalized CETSA fractions")

rep_cv <- long[is.finite(frac), .(cv_pct = cv_percent(frac)), by = .(Protein_group, Temp, condition)]
p_qc_cv <- ggplot(
  rep_cv[is.finite(cv_pct)],
  aes(x = factor(Temp), y = cv_pct, fill = condition)
) +
  geom_boxplot(outlier.size = 0.2) +
  theme_bw(base_size = 10) +
  labs(x = "Temperature", y = "Replicate CV (%)", title = "Replicate CV across temperatures")

if (flag_outlier_reps && "mzscore_ref" %in% names(long)) {
  flag_dt <- unique(long[, .(Protein_group, condition, Replicate, mzscore_ref, outlier_rep)])
  p_qc_flags <- ggplot(flag_dt[is.finite(mzscore_ref)], aes(x = mzscore_ref, fill = outlier_rep)) +
    geom_histogram(bins = 80, position = "identity", alpha = 0.7) +
    geom_vline(xintercept = c(-outlier_mzscore_threshold, outlier_mzscore_threshold), linetype = "dashed", colour = "red") +
    facet_wrap(~condition) +
    theme_bw(base_size = 10) +
    labs(x = "Modified Z-score of reference value", y = "Count", title = "Outlier replicate detection")
  ggsave(file.path(output_dir, "QC_outlier_rep_flags.png"), p_qc_flags,
         width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
}

if (nrow(ratio_long) > 0) {
  p_ratio <- ggplot(
    ratio_long[is.finite(log2_HL_ratio)],
    aes(x = factor(Temp), y = log2_HL_ratio, fill = factor(Temp))
  ) +
    geom_boxplot(outlier.size = 0.2, show.legend = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    theme_bw(base_size = 10) +
    labs(
      x = "Temperature (°C)",
      y = "normalized log2(Drug/DMSO) H/L ratio",
      title = "Track A: SILAC ratio distribution per temperature"
    )
  ggsave(file.path(output_dir, "QC_ratio_distributions.png"), p_ratio,
         width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
}

make_pca_plot <- function(dt_long, cond_name) {
  dsub <- dt_long[condition == cond_name]
  mat <- dcast(dsub, Protein_group ~ Sample, value.var = "frac")
  mat <- mat[complete.cases(mat)]
  if (nrow(mat) < 10 || ncol(mat) < 3) return(NULL)
  mat[, Protein_group := NULL]
  X <- t(as.matrix(mat))
  pca <- tryCatch(prcomp(X, scale. = TRUE), error = function(e) NULL)
  if (is.null(pca)) return(NULL)
  pca_dt <- data.table(Sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
  ann <- unique(dsub[, .(Sample, Temp, Replicate)])
  pca_dt <- merge(pca_dt, ann, by = "Sample", all.x = TRUE)
  vexp <- summary(pca)$importance["Proportion of Variance", 1:2]
  ggplot(pca_dt, aes(PC1, PC2, color = Temp, shape = Replicate)) +
    geom_point(size = 3) +
    theme_bw(base_size = 10) +
    labs(
      title = paste0("PCA - ", cond_name),
      x = paste0("PC1 (", round(100 * vexp[1], 1), "%)"),
      y = paste0("PC2 (", round(100 * vexp[2], 1), "%)")
    )
}

p_pca_dmso <- make_pca_plot(long, "DMSO")
p_pca_drug <- make_pca_plot(long, "Drug")

ggsave(file.path(output_dir, "QC_abundance_boxplots.png"), p_qc_box,
       width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
ggsave(file.path(output_dir, "QC_fraction_boxplots.png"), p_qc_frac,
       width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
ggsave(file.path(output_dir, "QC_replicate_CV.png"), p_qc_cv,
       width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")

if (!is.null(p_pca_dmso)) {
  ggsave(file.path(output_dir, "QC_PCA_DMSO.png"), p_pca_dmso,
         width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
}
if (!is.null(p_pca_drug)) {
  ggsave(file.path(output_dir, "QC_PCA_Drug.png"), p_pca_drug,
         width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
}

#############################################################################################################
# HEATMAP
#############################################################################################################

heat_source <- if (nrow(ratio_temp_summary) > 0) {
  dcast(ratio_temp_summary, Protein_group + Genes ~ Temp, value.var = "mean_log2_ratio")
} else {
  dcast(temp_effect_summary, Protein_group + Genes ~ Temp, value.var = "mean_log2_ratio")
}

if (nrow(heat_source) > 2) {
  
  heat_mat <- as.matrix(heat_source[, -c("Protein_group", "Genes")])
  mode(heat_mat) <- "numeric"
  
  rownames(heat_mat) <- paste(heat_source$Genes, heat_source$Protein_group, sep = " | ")
  
  # Remove rows with no finite values
  finite_row <- rowSums(is.finite(heat_mat)) >= 2
  heat_mat <- heat_mat[finite_row, , drop = FALSE]
  
  if (nrow(heat_mat) > 2) {
    
    # Rank rows by maximum absolute finite effect
    row_score <- apply(heat_mat, 1, function(x) {
      x <- x[is.finite(x)]
      if (length(x) == 0) return(-Inf)
      max(abs(x), na.rm = TRUE)
    })
    
    keep_n <- min(200, nrow(heat_mat))
    heat_mat_top <- heat_mat[order(row_score, decreasing = TRUE)[seq_len(keep_n)], , drop = FALSE]
    
    # Replace remaining NA/NaN/Inf values with 0 for clustering/plotting
    heat_mat_top[!is.finite(heat_mat_top)] <- 0
    
    # Remove rows with zero variance after replacement
    row_var <- apply(heat_mat_top, 1, var, na.rm = TRUE)
    heat_mat_top <- heat_mat_top[is.finite(row_var) & row_var > 0, , drop = FALSE]
    
    if (nrow(heat_mat_top) > 2) {
      png(
        file.path(output_dir, "Heatmap_top_mean_log2_Drug_vs_DMSO.png"),
        width = 1800,
        height = 2200,
        res = 200
      )
      
      pheatmap(
        heat_mat_top,
        scale = "none",
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        main = "Top proteins: mean normalized log2(Drug/DMSO) across temperatures",
        fontsize_row = 6,
        fontsize_col = 10,
        border_color = NA,
        na_col = "grey90"
      )
      
      dev.off()
    } else {
      msg("Heatmap skipped: fewer than 3 variable rows after NA filtering.")
    }
    
  } else {
    msg("Heatmap skipped: fewer than 3 proteins with enough finite temperature values.")
  }
  
} else {
  msg("Heatmap skipped: heat_source has too few rows.")
}
#############################################################################################################
# SUMMARY PLOTS
#############################################################################################################

p_rank <- ggplot(
  hit_table[is.finite(deltaTm_mean) | is.finite(deltaAUC_mean)],
  aes(x = deltaTm_mean, y = deltaAUC_mean)
) +
  geom_point(alpha = 0.5) +
  theme_bw(base_size = 11) +
  labs(
    title = "Track B supportive protein-level CETSA shift",
    x = "Mean deltaTm (Drug - DMSO) [°C]",
    y = "Mean deltaAUC (Drug - DMSO)"
  )

ggsave(file.path(output_dir, "Summary_deltaTm_vs_deltaAUC.png"), p_rank,
       width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")

if (any(is.finite(hit_table$p_value_best)) && any(is.finite(hit_table$max_abs_mean_log2_HL))) {
  hit_table[, neg_log10_p := -log10(p_value_best)]
  p_vol <- ggplot(
    hit_table[is.finite(signed_peak_log2_HL) & is.finite(neg_log10_p)],
    aes(x = signed_peak_log2_HL, y = neg_log10_p)
  ) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "red") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    theme_bw(base_size = 11) +
    labs(
      title = "Track A volcano: SILAC ratio effect vs significance",
      x = "Peak normalized log2(Drug/DMSO)",
      y = expression(-log[10](p-value))
    )
  ggsave(file.path(output_dir, "Summary_volcano.png"), p_vol,
         width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
}

#############################################################################################################
# PROTEIN CURVE PLOTS
#############################################################################################################

plot_one_protein <- function(protein_id = NULL, gene_name = NULL) {
  sel <- if (!is.null(protein_id)) quote(Protein_group == protein_id) else quote(Genes == gene_name)
  
  dt_raw   <- long[eval(sel)]
  dt_mean  <- curve_mean[eval(sel)]
  dt_pred  <- if (nrow(fit_preds_mean) > 0) fit_preds_mean[eval(sel)] else data.table()
  dt_ratio <- if (nrow(ratio_long) > 0) ratio_long[eval(sel)] else data.table()
  hit_row  <- hit_table[eval(sel)]
  
  if (nrow(dt_raw) == 0) return(NULL)
  
  ttl    <- unique(dt_raw$Genes)[1]
  subttl <- unique(dt_raw$Protein_group)[1]
  
  extra <- ""
  if (nrow(hit_row) > 0) {
    peak <- hit_row$signed_peak_log2_HL[1]
    dTm  <- hit_row$deltaTm_mean[1]
    dTm2 <- if ("deltaTm_mean_curve" %in% names(hit_row)) hit_row$deltaTm_mean_curve[1] else NA_real_
    dAUC <- if ("deltaAUC_mean" %in% names(hit_row)) hit_row$deltaAUC_mean[1] else NA_real_
    pval <- if ("p_value_best" %in% names(hit_row)) hit_row$p_value_best[1] else NA_real_
    extra <- paste0(
      "peak log2HL=", ifelse(is.finite(peak), round(peak, 3), "NA"),
      " | deltaTm=", ifelse(is.finite(dTm), round(dTm, 2), "NA"),
      " | deltaTm_mean_curve=", ifelse(is.finite(dTm2), round(dTm2, 2), "NA"),
      " | deltaAUC=", ifelse(is.finite(dAUC), round(dAUC, 3), "NA"),
      " | p=", ifelse(is.finite(pval), signif(pval, 2), "NA")
    )
  }
  
  p_abs <- ggplot() +
    geom_point(
      data = dt_raw[is.finite(frac)],
      aes(Temp, frac, colour = condition, shape = Replicate),
      size = 2, alpha = 0.7
    ) +
    geom_line(
      data = dt_raw[is.finite(frac)],
      aes(Temp, frac, colour = condition, group = interaction(condition, Replicate)),
      alpha = 0.4
    ) +
    geom_errorbar(
      data = dt_mean[is.finite(mean_frac)],
      aes(Temp, ymin = mean_frac - sem_frac, ymax = mean_frac + sem_frac, colour = condition),
      width = 0.3, linewidth = 0.5
    ) +
    geom_point(data = dt_mean[is.finite(mean_frac)], aes(Temp, mean_frac, colour = condition), size = 3) +
    geom_line(data = dt_mean[is.finite(mean_frac)], aes(Temp, mean_frac, colour = condition), linewidth = 1) +
    theme_bw(base_size = 11) +
    ylim(0, max_fraction_cap) +
    labs(
      x = "Temperature (°C)",
      y = "Fraction vs lowest temperature",
      title = ttl,
      subtitle = paste0(subttl, "\n", extra)
    )
  
  if (nrow(dt_pred) > 0) {
    p_abs <- p_abs +
      geom_line(
        data = dt_pred[is.finite(frac_pred)],
        aes(Temp, frac_pred, colour = condition),
        linewidth = 1.2, linetype = 2
      )
  }
  
  if (nrow(dt_ratio) > 0 && any(is.finite(dt_ratio$log2_HL_ratio))) {
    dt_ratio_mean <- dt_ratio[is.finite(log2_HL_ratio), .(
      mean_hl = safe_mean(log2_HL_ratio),
      sem_hl  = safe_sem(log2_HL_ratio)
    ), by = Temp]
    
    p_ratio_panel <- ggplot() +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_point(
        data = dt_ratio[is.finite(log2_HL_ratio)],
        aes(Temp, log2_HL_ratio, shape = Replicate),
        size = 2, alpha = 0.7, colour = "#E69F00"
      ) +
      geom_errorbar(
        data = dt_ratio_mean[is.finite(mean_hl)],
        aes(Temp, ymin = mean_hl - sem_hl, ymax = mean_hl + sem_hl),
        width = 0.3, colour = "#E69F00"
      ) +
      geom_line(
        data = dt_ratio_mean[is.finite(mean_hl)],
        aes(Temp, mean_hl), colour = "#E69F00", linewidth = 1
      ) +
      theme_bw(base_size = 11) +
      labs(
        x = "Temperature (°C)",
        y = "normalized log2(Drug/DMSO)",
        title = "Track A: intra-sample SILAC ratio"
      )
    
    return(p_abs / p_ratio_panel + plot_layout(heights = c(2, 1)))
  }
  
  p_abs
}

if (length(proteins_of_interest) > 0) {
  for (g in proteins_of_interest) {
    p <- plot_one_protein(gene_name = g)
    if (!is.null(p)) {
      ggsave(
        file.path(output_dir, paste0("ProteinCurve_", gsub("[^A-Za-z0-9_\\-]", "_", g), ".png")),
        p, width = plot_width, height = plot_height * 1.4, dpi = plot_dpi, bg = "white"
      )
    }
  }
}

top_hits <- hit_table[seq_len(min(top_n_deltaTm_plots, .N))]
if (nrow(top_hits) > 0) {
  dir.create(file.path(output_dir, "Top_hit_curves"), showWarnings = FALSE)
  for (i in seq_len(nrow(top_hits))) {
    pg <- top_hits$Protein_group[i]
    gn <- top_hits$Genes[i]
    p  <- plot_one_protein(protein_id = pg)
    if (!is.null(p)) {
      fname <- paste0(
        sprintf("%02d", i), "_",
        gsub("[^A-Za-z0-9_\\-]", "_", gn), "_",
        gsub("[^A-Za-z0-9_\\-]", "_", pg), ".png"
      )
      ggsave(
        file.path(output_dir, "Top_hit_curves", fname),
        p,
        width = plot_width,
        height = if (nrow(ratio_long) > 0) plot_height * 1.4 else plot_height,
        dpi = plot_dpi,
        bg = "white"
      )
    }
  }
}

#############################################################################################################
# WRITE OUTPUT TABLES
#############################################################################################################

msg("\nWriting output tables...")

fwrite(long,                file.path(output_dir, "long_fraction_table.tsv.gz"), sep = "\t")
fwrite(curve_mean,          file.path(output_dir, "mean_curves.tsv.gz"), sep = "\t")
fwrite(rep_wide,            file.path(output_dir, "replicate_drug_vs_dmso_by_temp.tsv.gz"), sep = "\t")
fwrite(temp_effect_summary, file.path(output_dir, "temp_effect_summary.tsv.gz"), sep = "\t")

if (nrow(ratio_long) > 0) {
  fwrite(ratio_long,         file.path(output_dir, "long_ratio_table.tsv.gz"), sep = "\t")
  fwrite(ratio_mean,         file.path(output_dir, "mean_ratio_curves.tsv.gz"), sep = "\t")
  fwrite(ratio_temp_summary, file.path(output_dir, "temp_ratio_summary.tsv.gz"), sep = "\t")
}

if (nrow(fit_results_rep) > 0)  fwrite(fit_results_rep,  file.path(output_dir, "fit_results_replicates.tsv.gz"), sep = "\t")
if (nrow(fit_preds_rep) > 0)    fwrite(fit_preds_rep,    file.path(output_dir, "fit_predictions_replicates.tsv.gz"), sep = "\t")
if (nrow(deltaTm_rep) > 0)      fwrite(deltaTm_rep,      file.path(output_dir, "deltaTm_replicates.tsv.gz"), sep = "\t")

if (nrow(fit_results_mean) > 0) fwrite(fit_results_mean, file.path(output_dir, "fit_results_mean_curves.tsv.gz"), sep = "\t")
if (nrow(fit_preds_mean) > 0)   fwrite(fit_preds_mean,   file.path(output_dir, "fit_predictions_mean_curves.tsv.gz"), sep = "\t")
if (nrow(deltaTm_mean) > 0)     fwrite(deltaTm_mean,     file.path(output_dir, "deltaTm_mean_curves.tsv.gz"), sep = "\t")

fwrite(hit_table, file.path(output_dir, "CETSA_SILAC_hit_table.tsv.gz"), sep = "\t")

#############################################################################################################
# CONSOLE SUMMARY
#############################################################################################################

msg("\n══════════════════════════════════════════════════════")
msg("Done. Output directory: ", output_dir)
msg("══════════════════════════════════════════════════════")
msg("Proteins in final hit table              : ", nrow(hit_table))
msg("  with SILAC ratio score  (Track A)      : ", hit_table[is.finite(max_abs_mean_log2_HL), .N])
msg("  with replicate deltaTm  (Track B)      : ", hit_table[is.finite(deltaTm_mean), .N])
msg("  with mean-curve deltaTm (Track B)      : ", hit_table[is.finite(deltaTm_mean_curve), .N])
if (any(is.finite(hit_table$p_value_best))) {
  msg("  with p-value < 0.05                    : ", hit_table[is.finite(p_value_best) & p_value_best < 0.05, .N])
}

msg("\nTop proteins by Track A-first rank score:")
show_cols <- intersect(
  c(
    "Genes", "Protein_group", "evidence_track", "stabilization_direction",
    "signed_peak_log2_HL", "max_abs_mean_log2_HL", "temp_of_max_abs_log2HL",
    "p_value_ratio_track", "deltaTm_mean", "deltaTm_mean_curve",
    "deltaAUC_mean", "p_value_best", "n_replicates_present"
  ),
  names(hit_table)
)
print(hit_table[seq_len(min(20, .N)), ..show_cols])
