#############################################################################################################
#
# CETSA-MS SILAC ANALYSIS SCRIPT  v2.0
#
# DESIGN
# ------
# Input : protein_quant.csv.gz from the DIA-NN SILAC protein summarisation script
#
# Channel meaning (inherited from upstream script):
#   LFQ_L  / label L  =  light  =  DMSO / vehicle / control
#   LFQ_H  / label H  =  heavy  =  drug / treated
#   log2_HL_ratio      =  log2(Drug / DMSO), median across matched precursor pairs
#                         computed from Precursor.Translated by the upstream script
#
# TWO COMPLEMENTARY ANALYSIS TRACKS
# ----------------------------------
#  Track A  –  SILAC ratio track (primary)
#              Uses log2_HL_ratio (intra-sample H/L) per temperature per replicate.
#              Fits a "ratio melt curve": how Drug/DMSO diverges with temperature.
#              Baseline = 0 (no shift).  No cross-sample normalization needed.
#              Most robust: cancels protein-abundance variation within each run.
#
#  Track B  –  Absolute melt curve track (secondary, mirrors original script)
#              Normalizes LFQ_L and LFQ_H to fraction remaining vs lowest temp.
#              Fits separate sigmoid melt curves for DMSO and Drug.
#              deltaTm = Tm(Drug) – Tm(DMSO).
#              More familiar output, but noisier (cross-sample normalization).
#
# NORMALIZATION IMPROVEMENTS vs v1
# ----------------------------------
#  1. Reference normalization uses median (not mean) of lowest-temperature
#     replicates, removing sensitivity to single outlier runs.
#  2. Global sample normalization is applied separately per channel and only
#     to Track B.  Track A is immune by construction.
#  3. Per-protein outlier replicate flagging (modified Z-score on ref value)
#     before fraction calculation; outlier reps are down-weighted / excluded.
#  4. Curve fitting uses multi-start optimization (3 starting Tm values) to
#     escape local minima.  Best fit by RSS is kept.
#  5. Fit quality gate: r² ≥ 0.7 AND bottom < top (monotone decreasing)
#     required for a replicate fit to enter the deltaTm calculation.
#  6. Hit scoring uses a single coherent score combining both tracks, with
#     a Welch t-test p-value on per-replicate deltaTm values.
#  7. Ratio track additionally reports AUC of the ratio curve and the
#     temperature of maximum Drug/DMSO divergence.
#
# OUTPUT
# ------
#  Tables (TSV.GZ):
#    long_fraction_table.tsv.gz            – normalized fractions, both channels
#    long_ratio_table.tsv.gz               – per-rep log2(H/L) ratio per temp
#    mean_curves.tsv.gz                    – mean ± SEM fraction curves
#    mean_ratio_curves.tsv.gz              – mean ± SEM ratio curves
#    replicate_drug_vs_dmso_by_temp.tsv.gz – Drug vs DMSO per temp per rep
#    temp_effect_summary.tsv.gz            – mean log2 ratio per temp
#    fit_results_replicates.tsv.gz         – sigmoid fit params per rep
#    fit_results_mean_curves.tsv.gz        – sigmoid fit params for mean curves
#    deltaTm_replicates.tsv.gz             – per-rep deltaTm
#    deltaTm_mean_curves.tsv.gz            – mean-curve deltaTm
#    CETSA_SILAC_hit_table.tsv.gz          – final ranked hit table
#
#  QC plots (PNG):
#    QC_abundance_boxplots.png
#    QC_fraction_boxplots.png
#    QC_replicate_CV.png
#    QC_ratio_distributions.png
#    QC_PCA_DMSO.png / QC_PCA_Drug.png
#    QC_outlier_rep_flags.png
#    Heatmap_top_mean_log2_Drug_vs_DMSO.png
#    Summary_deltaTm_vs_deltaAUC.png
#    Summary_volcano.png
#    Top_hit_curves/  (individual protein plots)
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
output_dir  <- "CETSA_SILAC_output_v2"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Filtering ─────────────────────────────────────────────────────────────────
min_precursors_per_channel <- 2   # minimum precursors required in each channel
min_ratios_for_HL          <- 2   # minimum matched precursor pairs for log2_HL_ratio
min_temps_per_curve        <- 6   # minimum finite temperature points to attempt fit

# ── Outlier replicate detection ───────────────────────────────────────────────
# Modified Z-score on the per-protein reference (lowest-temp) value.
# Replicates with |mZ| > outlier_mzscore_threshold are flagged.
# Flagged replicates are excluded from Track B fraction curves (Track A unaffected).
flag_outlier_reps     <- TRUE
outlier_mzscore_threshold <- 3.5   # standard Iglewicz-Hoaglin threshold

# ── Normalization (Track B only) ──────────────────────────────────────────────
do_global_sample_median_normalization <- TRUE
use_log2_for_global_normalization     <- FALSE  # FALSE = multiplicative factor (recommended)
max_fraction_cap   <- 1.5   # cap applied after fraction calculation
min_fraction_floor <- 0.0

# ── Curve fitting ─────────────────────────────────────────────────────────────
fit_on_replicates  <- TRUE
fit_on_mean_curve  <- TRUE
n_fit_starts       <- 3     # number of Tm starting values for multi-start fitting
min_r2_for_fit     <- 0.7   # minimum r² to accept a replicate fit
require_monotonic_soft <- FALSE  # if TRUE, rejects curves with >50% upward steps

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
safe_sem  <- function(x) { x <- x[is.finite(x)]; n <- length(x); if (n < 2) NA_real_ else sd(x)/sqrt(n) }
safe_n    <- function(x)   sum(is.finite(x))

safe_auc <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(NA_real_)
  ord <- order(x); x <- x[ord]; y <- y[ord]
  sum(diff(x) * (head(y,-1) + tail(y,-1)) / 2)
}

cv_percent <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2 || mean(x) == 0) return(NA_real_)
  100 * sd(x) / mean(x)
}

# Modified Z-score (Iglewicz & Hoaglin 1993)
modified_zscore <- function(x) {
  med <- median(x, na.rm = TRUE)
  mad <- median(abs(x - med), na.rm = TRUE)
  if (mad == 0) mad <- mean(abs(x - med), na.rm = TRUE) / 0.6745
  if (mad == 0) return(rep(0, length(x)))
  0.6745 * (x - med) / mad
}

welch_t_pvalue <- function(x, mu0 = 0) {
  # one-sample Welch t-test: H0 = mean(x) == mu0
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2) return(NA_real_)
  t_stat <- (mean(x) - mu0) / (sd(x) / sqrt(n))
  2 * pt(-abs(t_stat), df = n - 1)
}

#############################################################################################################
# SIGMOID FIT
#############################################################################################################

sigmoid_fun <- function(temp, top, bottom, Tm, slope) {
  bottom + (top - bottom) / (1 + exp((temp - Tm) / slope))
}

# Single sigmoid fit attempt with given starting values
.try_nls <- function(x, y, starts) {
  warn_txt <- character()
  fit <- tryCatch(
    withCallingHandlers(
      nls(
        frac ~ bottom + (top - bottom) / (1 + exp((Temp - Tm) / slope)),
        data      = data.frame(Temp = x, frac = y),
        start     = starts,
        algorithm = "port",
        lower     = c(top = 0.5,  bottom = -0.2, Tm = min(x), slope = 0.15),
        upper     = c(top = 1.5,  bottom =  1.0, Tm = max(x), slope = 25),
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
  dt   <- dt[is.finite(Temp) & is.finite(frac)][order(Temp)]
  out  <- list(success = FALSE, n_points = nrow(dt),
               top = NA_real_, bottom = NA_real_, Tm = NA_real_, slope = NA_real_,
               rss = NA_real_, r2 = NA_real_, auc = NA_real_,
               fit = NULL, warning_msg = NA_character_)
  
  if (nrow(dt) < min_temps_per_curve)              return(out)
  if (uniqueN(dt$Temp) < min_temps_per_curve)      return(out)
  if (sum(is.finite(dt$frac)) < min_temps_per_curve) return(out)
  if (all(abs(dt$frac - dt$frac[1]) < 1e-8, na.rm = TRUE)) return(out)
  
  x <- dt$Temp
  y <- pmin(pmax(dt$frac, min_fraction_floor), max_fraction_cap)
  
  if (require_monotonic_soft) {
    dy <- diff(y[order(x)])
    if (is.finite(mean(dy > 0, na.rm = TRUE)) && mean(dy > 0, na.rm = TRUE) > 0.50)
      return(out)
  }
  
  top_start    <- min(max(quantile(y, 0.9, na.rm = TRUE), 0.7), 1.3)
  bottom_start <- min(max(quantile(y, 0.1, na.rm = TRUE), 0.0), 0.6)
  slope_start  <- 2
  
  # Multi-start: candidate Tm values spread across the temperature range
  midpoint     <- (top_start + bottom_start) / 2
  Tm_data      <- x[which.min(abs(y - midpoint))]
  if (!is.finite(Tm_data)) Tm_data <- median(x, na.rm = TRUE)
  Tm_candidates <- unique(c(
    Tm_data,
    quantile(x, 0.35, na.rm = TRUE),
    quantile(x, 0.65, na.rm = TRUE)
  ))[seq_len(n_fit_starts)]
  
  best_rss <- Inf
  best_res <- NULL
  all_warnings <- character()
  
  for (Tm_s in Tm_candidates) {
    starts <- list(top = top_start, bottom = bottom_start,
                   Tm = Tm_s, slope = slope_start)
    res <- .try_nls(x, y, starts)
    all_warnings <- c(all_warnings, res$warnings)
    if (inherits(res$fit, "error")) next
    
    coefs <- coef(res$fit)
    if (!all(is.finite(coefs[c("top","bottom","Tm","slope")]))) next
    
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
  
  # Quality gates
  if (!is.finite(r2))            return(out)
  if (r2 < min_r2_for_fit)       return(out)
  if (coefs["bottom"] >= coefs["top"]) return(out)   # must be decreasing
  
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

# Optional but important columns from upstream script
has_HL_ratio <- "log2_HL_ratio" %in% names(proteins_raw)
has_N_ratios <- "N_ratios"      %in% names(proteins_raw)

if (!has_HL_ratio)
  msg("NOTE: log2_HL_ratio column not found – Track A (ratio track) will be skipped.")

proteins_raw[, Temp      := as.numeric(Temp)]
proteins_raw[, Replicate := as.character(Replicate)]
proteins_raw[, Genes     := fifelse(is.na(Genes) | Genes == "", Protein_group, Genes)]
setorder(proteins_raw, Protein_group, Replicate, Temp)

msg("Rows in input       : ", nrow(proteins_raw))
msg("Unique proteins     : ", uniqueN(proteins_raw$Protein_group))
msg("Temperatures        : ", paste(sort(unique(proteins_raw$Temp)), collapse = ", "))
msg("Replicates          : ", paste(sort(unique(proteins_raw$Replicate)), collapse = ", "))

#############################################################################################################
# PRECURSOR-COUNT FILTER
# Zero out LFQ values where either channel has too few precursors
#############################################################################################################

proteins_filt <- copy(proteins_raw)

proteins_filt[
  (is.na(N_precursors_L) | N_precursors_L < min_precursors_per_channel) |
    (is.na(N_precursors_H) | N_precursors_H < min_precursors_per_channel),
  c("LFQ_L", "LFQ_H") := .(NA_real_, NA_real_)
]

# Filter log2_HL_ratio by minimum matched pairs
if (has_HL_ratio && has_N_ratios) {
  proteins_filt[
    is.na(N_ratios) | N_ratios < min_ratios_for_HL,
    log2_HL_ratio := NA_real_
  ]
}

#############################################################################################################
# TRACK A : SILAC RATIO TRACK
# log2(Drug/DMSO) per protein per temperature per replicate
# No cross-sample normalization; cancels run-level variation within each SILAC pair
#############################################################################################################

if (has_HL_ratio) {
  msg("\n── Track A: SILAC ratio track ─────────────────────────────────────────")
  
  ratio_long <- proteins_filt[
    !is.na(log2_HL_ratio),
    .(Protein_group, Genes, Run, Sample, Temp, Replicate, log2_HL_ratio)
  ]
  
  # Mean ratio curve per protein x temperature
  ratio_mean <- ratio_long[, .(
    mean_log2_HL  = safe_mean(log2_HL_ratio),
    sd_log2_HL    = safe_sd(log2_HL_ratio),
    sem_log2_HL   = safe_sem(log2_HL_ratio),
    n_rep         = safe_n(log2_HL_ratio)
  ), by = .(Protein_group, Genes, Temp)]
  
  # Per-temperature summary (same structure as original script for compatibility)
  ratio_temp_summary <- ratio_long[, .(
    mean_log2_ratio = safe_mean(log2_HL_ratio),
    sd_log2_ratio   = safe_sd(log2_HL_ratio),
    sem_log2_ratio  = safe_sem(log2_HL_ratio),
    n_rep           = safe_n(log2_HL_ratio)
  ), by = .(Protein_group, Genes, Temp)]
  
  # Peak ratio summary per protein
  ratio_peak_summary <- ratio_temp_summary[, {
    ok <- is.finite(mean_log2_ratio)
    if (!any(ok)) {
      .(
        max_abs_mean_log2_HL   = NA_real_,
        temp_of_max_abs_log2HL = NA_real_,
        auc_log2_HL_curve      = NA_real_,
        p_value_ratio_track    = NA_real_
      )
    } else {
      idx <- which.max(abs(mean_log2_ratio))
      .(
        max_abs_mean_log2_HL   = max(abs(mean_log2_ratio), na.rm = TRUE),
        temp_of_max_abs_log2HL = Temp[idx],
        auc_log2_HL_curve      = safe_auc(Temp[ok], mean_log2_ratio[ok]),
        p_value_ratio_track    = NA_real_   # filled in below
      )
    }
  }, by = .(Protein_group, Genes)]
  
  # Welch t-test per protein: are per-rep log2(H/L) ratios consistently non-zero?
  # Use the temperature with maximum absolute mean ratio as representative
  ratio_ttest <- ratio_long[ratio_peak_summary[, .(Protein_group, temp_of_max_abs_log2HL)],
                            on = "Protein_group"][
                              !is.na(temp_of_max_abs_log2HL) & Temp == temp_of_max_abs_log2HL,
                              .(p_val = welch_t_pvalue(log2_HL_ratio, mu0 = 0)),
                              by = Protein_group
                            ]
  ratio_peak_summary[ratio_ttest, p_value_ratio_track := i.p_val, on = "Protein_group"]
  
  msg("Proteins with ≥1 ratio measurement: ", uniqueN(ratio_long$Protein_group))
  
  fwrite(ratio_long,          file.path(output_dir, "long_ratio_table.tsv.gz"), sep = "\t")
  fwrite(ratio_mean,          file.path(output_dir, "mean_ratio_curves.tsv.gz"), sep = "\t")
  fwrite(ratio_temp_summary,  file.path(output_dir, "temp_ratio_summary.tsv.gz"), sep = "\t")
} else {
  ratio_long         <- data.table()
  ratio_mean         <- data.table()
  ratio_temp_summary <- data.table()
  ratio_peak_summary <- data.table()
}

#############################################################################################################
# TRACK B : ABSOLUTE MELT CURVE TRACK
#############################################################################################################

msg("\n── Track B: Absolute melt curve track ────────────────────────────────")

# ── Build long table ───────────────────────────────────────────────────────────

long <- rbindlist(list(
  proteins_filt[, .(Protein_group, Genes, Run, Sample, Temp, Replicate,
                    channel = "L", condition = "DMSO",
                    abundance = LFQ_L, N_precursors = N_precursors_L)],
  proteins_filt[, .(Protein_group, Genes, Run, Sample, Temp, Replicate,
                    channel = "H", condition = "Drug",
                    abundance = LFQ_H, N_precursors = N_precursors_H)]
), use.names = TRUE)

long[abundance <= 0, abundance := NA_real_]
long[, abundance_norm := abundance]

# ── Global sample median normalization ────────────────────────────────────────

if (do_global_sample_median_normalization) {
  msg("Applying global sample median normalization (per channel)...")
  
  if (use_log2_for_global_normalization) {
    medtab <- long[is.finite(abundance), .(
      sample_median = median(log2(abundance), na.rm = TRUE)
    ), by = .(channel, Sample)]
    medtab[, global_median := median(sample_median, na.rm = TRUE), by = channel]
    medtab[, shift := global_median - sample_median]
    long <- merge(long, medtab[, .(channel, Sample, shift)],
                  by = c("channel","Sample"), all.x = TRUE)
    long[, abundance_norm := 2^(log2(abundance) + shift)]
    long[, shift := NULL]
  } else {
    medtab <- long[is.finite(abundance), .(
      sample_median = median(abundance, na.rm = TRUE)
    ), by = .(channel, Sample)]
    medtab[, global_median := median(sample_median, na.rm = TRUE), by = channel]
    medtab[, factor_norm := global_median / sample_median]
    long <- merge(long, medtab[, .(channel, Sample, factor_norm)],
                  by = c("channel","Sample"), all.x = TRUE)
    long[, abundance_norm := abundance * factor_norm]
    long[, factor_norm := NULL]
  }
}

# ── Reference value: MEDIAN over lowest-temperature replicates ───────────────
# Using median instead of mean protects against a single bad replicate
# distorting the entire protein's curve.

long[, temp_min := suppressWarnings(
  min(Temp[is.finite(abundance_norm)], na.rm = TRUE)
), by = .(Protein_group, Replicate, condition)]
long[!is.finite(temp_min), temp_min := NA_real_]

long[, ref_value := median(abundance_norm[Temp == temp_min], na.rm = TRUE),
     by = .(Protein_group, Replicate, condition)]
long[!is.finite(ref_value) | ref_value <= 0, ref_value := NA_real_]

# ── Outlier replicate detection ───────────────────────────────────────────────
# Flag replicates whose reference value is an extreme outlier across replicates
# for that protein x condition.  Flagged replicates are excluded from fraction
# curves but their ratio-track data (Track A) is unaffected.

if (flag_outlier_reps) {
  ref_dt <- unique(long[is.finite(ref_value),
                        .(Protein_group, condition, Replicate, ref_value)])
  
  ref_dt[, mzscore := modified_zscore(ref_value),
         by = .(Protein_group, condition)]
  ref_dt[, outlier_rep := abs(mzscore) > outlier_mzscore_threshold]
  
  long <- merge(
    long,
    ref_dt[, .(Protein_group, condition, Replicate, mzscore_ref = mzscore, outlier_rep)],
    by = c("Protein_group", "condition", "Replicate"),
    all.x = TRUE
  )
  long[is.na(outlier_rep), outlier_rep := FALSE]
  
  n_outlier_reps <- sum(long[, .(is_out = any(outlier_rep)), by = .(Protein_group, condition, Replicate)]$is_out)
  msg("Outlier replicate flags raised (protein x condition x rep): ", n_outlier_reps)
  
  # Exclude outlier reps from fraction calculation
  long[outlier_rep == TRUE, ref_value := NA_real_]
} else {
  long[, mzscore_ref := NA_real_]
  long[, outlier_rep := FALSE]
}

# ── Fraction remaining ────────────────────────────────────────────────────────

long[, frac := abundance_norm / ref_value]
long[, frac := pmin(pmax(frac, min_fraction_floor), max_fraction_cap)]

# ── Curve completeness filter ─────────────────────────────────────────────────

curve_stats <- long[, .(n_temps = sum(is.finite(frac))), by = .(Protein_group, Replicate, condition)]
keep_curves <- curve_stats[n_temps >= min_temps_per_curve]
long <- merge(long, keep_curves[, .(Protein_group, Replicate, condition)],
              by = c("Protein_group","Replicate","condition"))

msg("Protein x condition x replicate curves passing filter: ", nrow(keep_curves))

# ── Mean curves ───────────────────────────────────────────────────────────────

curve_mean <- long[, .(
  mean_frac = safe_mean(frac),
  sd_frac   = safe_sd(frac),
  sem_frac  = safe_sem(frac),
  n_rep     = safe_n(frac)
), by = .(Protein_group, Genes, condition, Temp)]

# ── Drug / DMSO comparison tables ────────────────────────────────────────────

rep_wide <- dcast(long, Protein_group + Genes + Replicate + Temp ~ condition,
                  value.var = "frac")
if (!("DMSO" %in% names(rep_wide))) rep_wide[, DMSO := NA_real_]
if (!("Drug" %in% names(rep_wide))) rep_wide[, Drug := NA_real_]

rep_wide[, log2_ratio_Drug_vs_DMSO := fifelse(
  is.finite(Drug) & is.finite(DMSO) & Drug > 0 & DMSO > 0,
  log2(Drug / DMSO), NA_real_
)]
rep_wide[, diff_Drug_minus_DMSO := Drug - DMSO]

mean_wide <- dcast(curve_mean, Protein_group + Genes + Temp ~ condition,
                   value.var = "mean_frac")
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

temp_grid <- seq(min(long$Temp, na.rm = TRUE),
                 max(long$Temp, na.rm = TRUE), length.out = 200)

run_fits <- function(input_dt, group_cols) {
  msg("  Fitting ", paste(group_cols, collapse="/"), " curves...")
  groups <- split(input_dt, by = group_cols, keep.by = TRUE)
  
  results <- rbindlist(lapply(groups, function(dtg) {
    fr <- fit_sigmoid_curve(dtg)
    info <- unique(dtg[, ..group_cols])
    cbind(info, data.table(
      fit_success = fr$success, n_points = fr$n_points,
      top = fr$top, bottom = fr$bottom, Tm = fr$Tm, slope = fr$slope,
      rss = fr$rss, r2 = fr$r2, auc_fit = fr$auc, fit_warning = fr$warning_msg
    ))
  }), fill = TRUE)
  
  preds <- rbindlist(lapply(groups, function(dtg) {
    fr   <- fit_sigmoid_curve(dtg)
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

if (fit_on_replicates) {
  rep_input <- long[, .(Protein_group, Genes, Replicate, condition, Temp, frac)]
  r <- run_fits(rep_input, c("Protein_group","Genes","Replicate","condition"))
  fit_results_rep <- r$results
  fit_preds_rep   <- r$preds
}

if (fit_on_mean_curve) {
  mean_input <- curve_mean[, .(Protein_group, Genes, condition, Temp, frac = mean_frac)]
  r <- run_fits(mean_input, c("Protein_group","Genes","condition"))
  fit_results_mean <- r$results
  fit_preds_mean   <- r$preds
}

#############################################################################################################
# deltaTm / deltaAUC SUMMARIES
#############################################################################################################

deltaTm_rep  <- data.table()
deltaTm_mean <- data.table()
replicate_summary <- data.table()

if (nrow(fit_results_rep) > 0 && any(fit_results_rep$fit_success)) {
  deltaTm_rep <- dcast(
    fit_results_rep[fit_success == TRUE],
    Protein_group + Genes + Replicate ~ condition,
    value.var = c("Tm","auc_fit","r2","n_points")
  )
  for (cn in c("Tm_Drug","Tm_DMSO","auc_fit_Drug","auc_fit_DMSO","r2_Drug","r2_DMSO"))
    if (!(cn %in% names(deltaTm_rep))) deltaTm_rep[, (cn) := NA_real_]
  
  deltaTm_rep[, deltaTm  := Tm_Drug  - Tm_DMSO]
  deltaTm_rep[, deltaAUC := auc_fit_Drug - auc_fit_DMSO]
  
  replicate_summary <- deltaTm_rep[, .(
    n_replicates        = .N,
    n_deltaTm           = sum(is.finite(deltaTm)),
    deltaTm_mean        = safe_mean(deltaTm),
    deltaTm_sd          = safe_sd(deltaTm),
    deltaTm_sem         = safe_sem(deltaTm),
    deltaTm_p_value     = welch_t_pvalue(deltaTm, mu0 = 0),
    deltaAUC_mean       = safe_mean(deltaAUC),
    deltaAUC_sd         = safe_sd(deltaAUC),
    Tm_DMSO_mean        = safe_mean(Tm_DMSO),
    Tm_Drug_mean        = safe_mean(Tm_Drug),
    r2_DMSO_mean        = safe_mean(r2_DMSO),
    r2_Drug_mean        = safe_mean(r2_Drug)
  ), by = .(Protein_group, Genes)]
}

if (nrow(fit_results_mean) > 0 && any(fit_results_mean$fit_success)) {
  deltaTm_mean <- dcast(
    fit_results_mean[fit_success == TRUE],
    Protein_group + Genes ~ condition,
    value.var = c("Tm","auc_fit","r2","n_points")
  )
  for (cn in c("Tm_Drug","Tm_DMSO","auc_fit_Drug","auc_fit_DMSO"))
    if (!(cn %in% names(deltaTm_mean))) deltaTm_mean[, (cn) := NA_real_]
  
  deltaTm_mean[, deltaTm_mean_curve  := Tm_Drug  - Tm_DMSO]
  deltaTm_mean[, deltaAUC_mean_curve := auc_fit_Drug - auc_fit_DMSO]
}

#############################################################################################################
# FINAL HIT TABLE
#############################################################################################################

msg("\nBuilding hit table...")

hit_table <- unique(proteins_raw[, .(Protein_group, Genes)])

# Track B metrics
if (nrow(replicate_summary) > 0)
  hit_table <- merge(hit_table, replicate_summary, by = c("Protein_group","Genes"), all.x = TRUE)
if (nrow(deltaTm_mean) > 0)
  hit_table <- merge(hit_table, deltaTm_mean[, .(
    Protein_group, Genes, deltaTm_mean_curve, deltaAUC_mean_curve
  )], by = c("Protein_group","Genes"), all.x = TRUE)

# Track A metrics
if (nrow(ratio_peak_summary) > 0)
  hit_table <- merge(hit_table, ratio_peak_summary, by = c("Protein_group","Genes"), all.x = TRUE)

# Curve completeness
curve_completeness <- long[, .(
  n_points_total       = sum(is.finite(frac)),
  n_replicates_present = uniqueN(Replicate[is.finite(frac)]),
  n_conditions_present = uniqueN(condition[is.finite(frac)])
), by = .(Protein_group, Genes)]
hit_table <- merge(hit_table, curve_completeness, by = c("Protein_group","Genes"), all.x = TRUE)

# ── Unified rank score ────────────────────────────────────────────────────────
# Primary  : |deltaTm_mean| from replicate-level fits  (Track B, most trusted)
# Secondary: |deltaTm_mean_curve| from mean-curve fit   (Track B fallback)
# Tertiary : |max_abs_mean_log2_HL| from SILAC ratios   (Track A)
# Quaternary: |max_abs_mean_log2_ratio| from MaxLFQ-derived Drug/DMSO fractions

peak_from_trackB <- temp_effect_summary[, {
  ok <- is.finite(mean_log2_ratio)
  if (!any(ok)) .(max_abs_mean_log2_ratio_trackB = NA_real_)
  else .(max_abs_mean_log2_ratio_trackB = max(abs(mean_log2_ratio[ok]), na.rm = TRUE))
}, by = .(Protein_group, Genes)]

hit_table <- merge(hit_table, peak_from_trackB, by = c("Protein_group","Genes"), all.x = TRUE)

# Ensure all score columns exist
for (cn in c("deltaTm_mean","deltaTm_mean_curve","max_abs_mean_log2_HL",
             "max_abs_mean_log2_ratio_trackB")) {
  if (!(cn %in% names(hit_table))) hit_table[, (cn) := NA_real_]
}

hit_table[, rank_score := fcase(
  is.finite(deltaTm_mean),                abs(deltaTm_mean),
  is.finite(deltaTm_mean_curve),          abs(deltaTm_mean_curve),
  is.finite(max_abs_mean_log2_HL),        max_abs_mean_log2_HL,
  is.finite(max_abs_mean_log2_ratio_trackB), max_abs_mean_log2_ratio_trackB,
  default = NA_real_
)]

# Combined p-value: take minimum of available p-values (most significant)
p_cols <- intersect(c("deltaTm_p_value","p_value_ratio_track"), names(hit_table))
if (length(p_cols) > 0) {
  hit_table[, p_value_best := do.call(pmin, c(.SD, na.rm = TRUE)), .SDcols = p_cols]
} else {
  hit_table[, p_value_best := NA_real_]
}

# Sort
hit_table[, .rank_abs   := fifelse(is.finite(rank_score),   abs(rank_score),   -Inf)]
hit_table[, .p_ord      := fifelse(is.finite(p_value_best), p_value_best,        Inf)]
hit_table[, .nrep       := fifelse(is.finite(n_replicates_present), n_replicates_present, -Inf)]
setorder(hit_table, .p_ord, -.rank_abs, -.nrep)
hit_table[, c(".rank_abs",".p_ord",".nrep") := NULL]

msg("Proteins in hit table                     : ", nrow(hit_table))
msg("With replicate deltaTm (Track B)          : ", hit_table[is.finite(deltaTm_mean), .N])
msg("With mean-curve deltaTm (Track B)         : ", hit_table[is.finite(deltaTm_mean_curve), .N])
if ("max_abs_mean_log2_HL" %in% names(hit_table))
  msg("With SILAC ratio score (Track A)          : ", hit_table[is.finite(max_abs_mean_log2_HL), .N])

#############################################################################################################
# QC PLOTS
#############################################################################################################

msg("\nGenerating QC plots...")

# ── Abundance distributions ───────────────────────────────────────────────────
p_qc_box <- ggplot(
  long[is.finite(abundance_norm)],
  aes(x = interaction(Temp, Replicate, lex.order = TRUE),
      y = log2(abundance_norm), fill = condition)
) +
  geom_boxplot(outlier.size = 0.2) +
  theme_bw(base_size = 10) +
  labs(x = "Temperature_Replicate", y = "log2 normalized abundance",
       title = "Protein abundance distributions") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# ── Fraction boxplots ─────────────────────────────────────────────────────────
p_qc_frac <- ggplot(
  long[is.finite(frac)],
  aes(x = factor(Temp), y = frac, fill = condition)
) +
  geom_boxplot(outlier.size = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  theme_bw(base_size = 10) +
  labs(x = "Temperature", y = "Normalized fraction",
       title = "Normalized CETSA fractions")

# ── Replicate CV ──────────────────────────────────────────────────────────────
rep_cv <- long[is.finite(frac), .(cv_pct = cv_percent(frac)),
               by = .(Protein_group, Temp, condition)]
p_qc_cv <- ggplot(
  rep_cv[is.finite(cv_pct)],
  aes(x = factor(Temp), y = cv_pct, fill = condition)
) +
  geom_boxplot(outlier.size = 0.2) +
  theme_bw(base_size = 10) +
  labs(x = "Temperature", y = "Replicate CV (%)", title = "Replicate CV across temperatures")

# ── Outlier rep flag summary ──────────────────────────────────────────────────
if (flag_outlier_reps && "mzscore_ref" %in% names(long)) {
  flag_dt <- unique(long[, .(Protein_group, condition, Replicate, mzscore_ref, outlier_rep)])
  p_qc_flags <- ggplot(flag_dt[is.finite(mzscore_ref)],
                       aes(x = mzscore_ref, fill = outlier_rep)) +
    geom_histogram(bins = 80, position = "identity", alpha = 0.7) +
    geom_vline(xintercept = c(-outlier_mzscore_threshold, outlier_mzscore_threshold),
               linetype = "dashed", colour = "red") +
    facet_wrap(~condition) +
    theme_bw(base_size = 10) +
    labs(x = "Modified Z-score of reference value", y = "Count",
         title = "Outlier replicate detection (ref value distribution)")
  ggsave(file.path(output_dir, "QC_outlier_rep_flags.png"), p_qc_flags,
         width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
}

# ── Track A ratio distributions ───────────────────────────────────────────────
if (nrow(ratio_long) > 0) {
  p_ratio <- ggplot(
    ratio_long[is.finite(log2_HL_ratio)],
    aes(x = factor(Temp), y = log2_HL_ratio, fill = factor(Temp))
  ) +
    geom_boxplot(outlier.size = 0.2, show.legend = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    theme_bw(base_size = 10) +
    labs(x = "Temperature (°C)", y = "log2(Drug/DMSO) H/L ratio",
         title = "Track A: intra-sample SILAC ratio distribution per temperature")
  ggsave(file.path(output_dir, "QC_ratio_distributions.png"), p_ratio,
         width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
}

# ── PCA ───────────────────────────────────────────────────────────────────────
make_pca_plot <- function(dt_long, cond_name) {
  dsub <- dt_long[condition == cond_name]
  mat  <- dcast(dsub, Protein_group ~ Sample, value.var = "frac")
  mat  <- mat[complete.cases(mat)]
  if (nrow(mat) < 10 || ncol(mat) < 3) return(NULL)
  mat[, Protein_group := NULL]
  X   <- t(as.matrix(mat))
  pca <- tryCatch(prcomp(X, scale. = TRUE), error = function(e) NULL)
  if (is.null(pca)) return(NULL)
  pca_dt <- data.table(Sample = rownames(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2])
  ann    <- unique(dsub[, .(Sample, Temp, Replicate)])
  pca_dt <- merge(pca_dt, ann, by = "Sample", all.x = TRUE)
  vexp   <- summary(pca)$importance["Proportion of Variance", 1:2]
  ggplot(pca_dt, aes(PC1, PC2, color = Temp, shape = Replicate)) +
    geom_point(size = 3) +
    theme_bw(base_size = 10) +
    labs(title = paste0("PCA – ", cond_name),
         x = paste0("PC1 (", round(100*vexp[1],1), "%)"),
         y = paste0("PC2 (", round(100*vexp[2],1), "%)"))
}

p_pca_dmso <- make_pca_plot(long, "DMSO")
p_pca_drug <- make_pca_plot(long, "Drug")

# ── Save QC plots ─────────────────────────────────────────────────────────────
ggsave(file.path(output_dir, "QC_abundance_boxplots.png"), p_qc_box,
       width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
ggsave(file.path(output_dir, "QC_fraction_boxplots.png"),  p_qc_frac,
       width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
ggsave(file.path(output_dir, "QC_replicate_CV.png"),       p_qc_cv,
       width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")

if (!is.null(p_pca_dmso))
  ggsave(file.path(output_dir, "QC_PCA_DMSO.png"), p_pca_dmso,
         width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
if (!is.null(p_pca_drug))
  ggsave(file.path(output_dir, "QC_PCA_Drug.png"), p_pca_drug,
         width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")

#############################################################################################################
# HEATMAP
#############################################################################################################

heat_source <- if (nrow(ratio_temp_summary) > 0) {
  dcast(ratio_temp_summary, Protein_group + Genes ~ Temp, value.var = "mean_log2_ratio")
} else {
  dcast(temp_effect_summary, Protein_group + Genes ~ Temp, value.var = "mean_log2_ratio")
}

if (nrow(heat_source) > 2) {
  heat_mat <- as.matrix(heat_source[, -c("Protein_group","Genes")])
  rownames(heat_mat) <- paste(heat_source$Genes, heat_source$Protein_group, sep = " | ")
  row_score <- matrixStats::rowMaxs(abs(heat_mat), na.rm = TRUE)
  row_score[!is.finite(row_score)] <- -Inf
  keep_n <- min(200, nrow(heat_mat))
  heat_mat_top <- heat_mat[order(row_score, decreasing = TRUE)[seq_len(keep_n)], , drop = FALSE]
  
  png(file.path(output_dir, "Heatmap_top_mean_log2_Drug_vs_DMSO.png"),
      width = 1800, height = 2200, res = 200)
  pheatmap(
    heat_mat_top, scale = "none", cluster_rows = TRUE, cluster_cols = FALSE,
    main  = "Top proteins: mean log2(Drug/DMSO) across temperatures",
    fontsize_row = 6, fontsize_col = 10, border_color = NA
  )
  dev.off()
}

#############################################################################################################
# SUMMARY PLOTS
#############################################################################################################

# ── deltaTm vs deltaAUC scatter ───────────────────────────────────────────────
p_rank <- ggplot(
  hit_table[is.finite(deltaTm_mean) | is.finite(deltaAUC_mean)],
  aes(x = deltaTm_mean, y = deltaAUC_mean)
) +
  geom_point(alpha = 0.5) +
  theme_bw(base_size = 11) +
  labs(title  = "Track B: protein-level CETSA shift",
       x = "Mean deltaTm (Drug – DMSO) [°C]",
       y = "Mean deltaAUC (Drug – DMSO)")

ggsave(file.path(output_dir, "Summary_deltaTm_vs_deltaAUC.png"), p_rank,
       width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")

# ── Volcano: deltaTm vs –log10(p-value) ──────────────────────────────────────
if (any(is.finite(hit_table$p_value_best)) && any(is.finite(hit_table$deltaTm_mean))) {
  hit_table[, neg_log10_p := -log10(p_value_best)]
  p_vol <- ggplot(
    hit_table[is.finite(deltaTm_mean) & is.finite(neg_log10_p)],
    aes(x = deltaTm_mean, y = neg_log10_p)
  ) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "red") +
    theme_bw(base_size = 11) +
    labs(title = "Volcano: deltaTm vs significance",
         x = "Mean deltaTm (Drug – DMSO) [°C]",
         y = expression(-log[10](p-value)))
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
  dt_ratio <- if (nrow(ratio_long)     > 0) ratio_long[eval(sel)]     else data.table()
  hit_row  <- hit_table[eval(sel)]
  
  if (nrow(dt_raw) == 0) return(NULL)
  
  ttl    <- unique(dt_raw$Genes)[1]
  subttl <- unique(dt_raw$Protein_group)[1]
  
  extra <- ""
  if (nrow(hit_row) > 0) {
    dTm  <- hit_row$deltaTm_mean[1]
    dTm2 <- if ("deltaTm_mean_curve" %in% names(hit_row)) hit_row$deltaTm_mean_curve[1] else NA
    dAUC <- if ("deltaAUC_mean" %in% names(hit_row))      hit_row$deltaAUC_mean[1]      else NA
    pval <- if ("p_value_best" %in% names(hit_row))        hit_row$p_value_best[1]       else NA
    extra <- paste0(
      "deltaTm=", ifelse(is.finite(dTm),  round(dTm, 2),  "NA"),
      " | deltaTm_mean_curve=", ifelse(is.finite(dTm2), round(dTm2, 2), "NA"),
      " | deltaAUC=",  ifelse(is.finite(dAUC), round(dAUC, 3), "NA"),
      " | p=", ifelse(is.finite(pval), signif(pval, 2), "NA")
    )
  }
  
  # ── Panel 1: absolute melt curves ─────────────────────────────────────────
  p_abs <- ggplot() +
    geom_point(data = dt_raw[is.finite(frac)],
               aes(Temp, frac, colour = condition, shape = Replicate),
               size = 2, alpha = 0.7) +
    geom_line(data = dt_raw[is.finite(frac)],
              aes(Temp, frac, colour = condition,
                  group = interaction(condition, Replicate)),
              alpha = 0.4) +
    geom_errorbar(data = dt_mean[is.finite(mean_frac)],
                  aes(Temp, ymin = mean_frac - sem_frac,
                      ymax = mean_frac + sem_frac, colour = condition),
                  width = 0.3, linewidth = 0.5) +
    geom_point(data = dt_mean[is.finite(mean_frac)],
               aes(Temp, mean_frac, colour = condition), size = 3) +
    geom_line(data = dt_mean[is.finite(mean_frac)],
              aes(Temp, mean_frac, colour = condition), linewidth = 1) +
    theme_bw(base_size = 11) +
    ylim(0, max_fraction_cap) +
    labs(x = "Temperature (°C)", y = "Normalized soluble fraction",
         title = ttl, subtitle = paste0(subttl, "\n", extra))
  
  if (nrow(dt_pred) > 0)
    p_abs <- p_abs +
    geom_line(data = dt_pred[is.finite(frac_pred)],
              aes(Temp, frac_pred, colour = condition),
              linewidth = 1.2, linetype = 2)
  
  # ── Panel 2: SILAC ratio curve ─────────────────────────────────────────────
  if (nrow(dt_ratio) > 0 && any(is.finite(dt_ratio$log2_HL_ratio))) {
    dt_ratio_mean <- dt_ratio[is.finite(log2_HL_ratio), .(
      mean_hl  = safe_mean(log2_HL_ratio),
      sem_hl   = safe_sem(log2_HL_ratio)
    ), by = Temp]
    
    p_ratio_panel <- ggplot() +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_point(data = dt_ratio[is.finite(log2_HL_ratio)],
                 aes(Temp, log2_HL_ratio, shape = Replicate),
                 size = 2, alpha = 0.7, colour = "#E69F00") +
      geom_errorbar(data = dt_ratio_mean[is.finite(mean_hl)],
                    aes(Temp, ymin = mean_hl - sem_hl, ymax = mean_hl + sem_hl),
                    width = 0.3, colour = "#E69F00") +
      geom_line(data = dt_ratio_mean[is.finite(mean_hl)],
                aes(Temp, mean_hl), colour = "#E69F00", linewidth = 1) +
      theme_bw(base_size = 11) +
      labs(x = "Temperature (°C)", y = "log2(Drug/DMSO) H/L ratio",
           title = "Track A: intra-sample SILAC ratio")
    
    return(p_abs / p_ratio_panel + plot_layout(heights = c(2, 1)))
  }
  
  p_abs
}

# ── Save individual protein plots ─────────────────────────────────────────────
if (length(proteins_of_interest) > 0) {
  for (g in proteins_of_interest) {
    p <- plot_one_protein(gene_name = g)
    if (!is.null(p))
      ggsave(
        file.path(output_dir, paste0("ProteinCurve_", gsub("[^A-Za-z0-9_\\-]","_",g), ".png")),
        p, width = plot_width, height = plot_height * 1.4, dpi = plot_dpi, bg = "white"
      )
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
      fname <- paste0(sprintf("%02d", i), "_",
                      gsub("[^A-Za-z0-9_\\-]","_", gn), "_",
                      gsub("[^A-Za-z0-9_\\-]","_", pg), ".png")
      ggsave(file.path(output_dir, "Top_hit_curves", fname),
             p, width = plot_width,
             height = if (nrow(ratio_long) > 0) plot_height * 1.4 else plot_height,
             dpi = plot_dpi, bg = "white")
    }
  }
}

#############################################################################################################
# WRITE OUTPUT TABLES
#############################################################################################################

msg("\nWriting output tables...")

fwrite(long,                file.path(output_dir, "long_fraction_table.tsv.gz"),            sep = "\t")
fwrite(curve_mean,          file.path(output_dir, "mean_curves.tsv.gz"),                    sep = "\t")
fwrite(rep_wide,            file.path(output_dir, "replicate_drug_vs_dmso_by_temp.tsv.gz"), sep = "\t")
fwrite(temp_effect_summary, file.path(output_dir, "temp_effect_summary.tsv.gz"),            sep = "\t")

if (nrow(ratio_long) > 0) {
  fwrite(ratio_long,         file.path(output_dir, "long_ratio_table.tsv.gz"),    sep = "\t")
  fwrite(ratio_mean,         file.path(output_dir, "mean_ratio_curves.tsv.gz"),   sep = "\t")
  fwrite(ratio_temp_summary, file.path(output_dir, "temp_ratio_summary.tsv.gz"),  sep = "\t")
}

if (nrow(fit_results_rep) > 0)  fwrite(fit_results_rep,  file.path(output_dir, "fit_results_replicates.tsv.gz"),    sep = "\t")
if (nrow(fit_preds_rep)   > 0)  fwrite(fit_preds_rep,    file.path(output_dir, "fit_predictions_replicates.tsv.gz"),sep = "\t")
if (nrow(deltaTm_rep)     > 0)  fwrite(deltaTm_rep,      file.path(output_dir, "deltaTm_replicates.tsv.gz"),        sep = "\t")

if (nrow(fit_results_mean) > 0) fwrite(fit_results_mean, file.path(output_dir, "fit_results_mean_curves.tsv.gz"),   sep = "\t")
if (nrow(fit_preds_mean)   > 0) fwrite(fit_preds_mean,   file.path(output_dir, "fit_predictions_mean_curves.tsv.gz"),sep = "\t")
if (nrow(deltaTm_mean)     > 0) fwrite(deltaTm_mean,     file.path(output_dir, "deltaTm_mean_curves.tsv.gz"),       sep = "\t")

fwrite(hit_table,           file.path(output_dir, "CETSA_SILAC_hit_table.tsv.gz"),          sep = "\t")

#############################################################################################################
# CONSOLE SUMMARY
#############################################################################################################

msg("\n══════════════════════════════════════════════════════")
msg("Done.  Output directory: ", output_dir)
msg("══════════════════════════════════════════════════════")
msg("Proteins in final hit table              : ", nrow(hit_table))
msg("  with replicate deltaTm  (Track B)      : ", hit_table[is.finite(deltaTm_mean), .N])
msg("  with mean-curve deltaTm (Track B)      : ", hit_table[is.finite(deltaTm_mean_curve), .N])
if ("max_abs_mean_log2_HL" %in% names(hit_table))
  msg("  with SILAC ratio score  (Track A)      : ", hit_table[is.finite(max_abs_mean_log2_HL), .N])
if (any(is.finite(hit_table$p_value_best)))
  msg("  with p-value < 0.05                    : ", hit_table[is.finite(p_value_best) & p_value_best < 0.05, .N])

msg("\nTop proteins by rank score:")
show_cols <- intersect(
  c("Genes","Protein_group","deltaTm_mean","deltaTm_mean_curve",
    "deltaAUC_mean","max_abs_mean_log2_HL","p_value_best","n_replicates_present"),
  names(hit_table)
)
print(hit_table[seq_len(min(20, .N)), ..show_cols])