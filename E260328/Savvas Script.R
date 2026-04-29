#############################################################################################################
#
# CETSA-MS SILAC ANALYSIS SCRIPT
# Fixed robust version
#
# DESIGN
# ------
# - Input: protein_quant.csv.gz from your DIA-NN SILAC protein summarisation script
# - LFQ_L = light = vehicle / DMSO / control
# - LFQ_H = heavy = drug / treated
# - Each run = one temperature point in one biological replicate
#
# OUTPUT
# ------
# - normalized long table
# - mean curves
# - replicate Drug-vs-DMSO per-temperature summaries
# - fitted curve tables
# - deltaTm / deltaAUC summaries
# - final hit table
# - QC plots
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
output_dir  <- "CETSA_SILAC_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Filtering
min_precursors_per_channel <- 2
min_temps_per_curve        <- 6

# Normalization
do_global_sample_median_normalization <- TRUE
use_log2_for_global_normalization     <- FALSE
max_fraction_cap   <- 1.5
min_fraction_floor <- 0

# Curve fitting
fit_on_replicates <- TRUE
fit_on_mean_curve <- TRUE
require_monotonic_soft <- FALSE

# Plotting
proteins_of_interest <- c(
  # "HSP90AA1", "NAMPT"
)
top_n_deltaTm_plots <- 24
plot_width  <- 10
plot_height <- 7
plot_dpi    <- 300

#############################################################################################################
# HELPERS
#############################################################################################################

msg <- function(...) cat(paste0(..., "\n"))

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
  n <- length(x)
  if (n < 2) return(NA_real_)
  sd(x) / sqrt(n)
}

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
  if (length(x) < 2) return(NA_real_)
  100 * sd(x) / mean(x)
}

safe_whichmax_abs <- function(x) {
  if (length(x) == 0) return(NA_integer_)
  if (all(!is.finite(x))) return(NA_integer_)
  which.max(abs(x))
}

#############################################################################################################
# SIGMOID FIT
#############################################################################################################

# decreasing sigmoid
sigmoid_fun <- function(temp, top, bottom, Tm, slope) {
  bottom + (top - bottom) / (1 + exp((temp - Tm) / slope))
}

fit_sigmoid_curve <- function(dt) {
  dt <- copy(dt)
  dt <- dt[is.finite(Temp) & is.finite(frac)]
  dt <- dt[order(Temp)]
  
  out <- list(
    success = FALSE,
    n_points = nrow(dt),
    top = NA_real_,
    bottom = NA_real_,
    Tm = NA_real_,
    slope = NA_real_,
    rss = NA_real_,
    r2 = NA_real_,
    auc = NA_real_,
    fit = NULL,
    warning_msg = NA_character_
  )
  
  if (nrow(dt) < min_temps_per_curve) return(out)
  
  x <- dt$Temp
  y <- dt$frac
  
  if (length(unique(x)) < min_temps_per_curve) return(out)
  if (sum(is.finite(y)) < min_temps_per_curve) return(out)
  if (all(abs(y - y[1]) < 1e-8, na.rm = TRUE)) return(out)
  
  if (require_monotonic_soft) {
    dy <- diff(y[order(x)])
    frac_up <- mean(dy > 0, na.rm = TRUE)
    if (is.finite(frac_up) && frac_up > 0.50) return(out)
  }
  
  y <- pmin(pmax(y, min_fraction_floor), max_fraction_cap)
  
  top_start <- min(max(quantile(y, 0.9, na.rm = TRUE), 0.7), 1.3)
  bottom_start <- min(max(quantile(y, 0.1, na.rm = TRUE), 0.0), 0.6)
  midpoint <- (top_start + bottom_start) / 2
  
  Tm_start <- x[which.min(abs(y - midpoint))]
  if (!is.finite(Tm_start)) Tm_start <- median(x, na.rm = TRUE)
  
  slope_start <- 2
  
  warn_txt <- character()
  
  fit_try <- tryCatch(
    withCallingHandlers(
      nls(
        frac ~ bottom + (top - bottom) / (1 + exp((Temp - Tm) / slope)),
        data = data.frame(Temp = x, frac = y),
        start = list(
          top = top_start,
          bottom = bottom_start,
          Tm = Tm_start,
          slope = slope_start
        ),
        algorithm = "port",
        lower = c(top = 0.5, bottom = -0.2, Tm = min(x), slope = 0.15),
        upper = c(top = 1.5, bottom = 1.0, Tm = max(x), slope = 25),
        control = nls.control(maxiter = 500, warnOnly = TRUE)
      ),
      warning = function(w) {
        warn_txt <<- c(warn_txt, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  
  if (inherits(fit_try, "error")) {
    out$warning_msg <- paste(unique(c(warn_txt, conditionMessage(fit_try))), collapse = " | ")
    return(out)
  }
  
  coef_fit <- coef(fit_try)
  pred <- tryCatch(predict(fit_try), error = function(e) rep(NA_real_, length(y)))
  
  if (!all(is.finite(coef_fit[c("top", "bottom", "Tm", "slope")]))) {
    out$warning_msg <- paste(unique(warn_txt), collapse = " | ")
    return(out)
  }
  
  rss <- sum((y - pred)^2, na.rm = TRUE)
  tss <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  r2  <- ifelse(tss > 0, 1 - rss / tss, NA_real_)
  
  # reject obviously bad fits
  if (!is.finite(r2)) r2 <- NA_real_
  if (is.finite(r2) && r2 < -0.5) {
    out$warning_msg <- paste(unique(c(warn_txt, "fit rejected: very poor r2")), collapse = " | ")
    return(out)
  }
  
  out$success <- TRUE
  out$top <- unname(coef_fit["top"])
  out$bottom <- unname(coef_fit["bottom"])
  out$Tm <- unname(coef_fit["Tm"])
  out$slope <- unname(coef_fit["slope"])
  out$rss <- rss
  out$r2 <- r2
  out$auc <- safe_auc(x, pred)
  out$fit <- fit_try
  out$warning_msg <- paste(unique(warn_txt), collapse = " | ")
  
  out
}

predict_fit_grid <- function(fit_obj, temp_grid) {
  if (is.null(fit_obj) || isFALSE(fit_obj$success) || is.null(fit_obj$fit)) return(NULL)
  pred <- tryCatch(
    predict(fit_obj$fit, newdata = data.frame(Temp = temp_grid)),
    error = function(e) NULL
  )
  if (is.null(pred)) return(NULL)
  data.table(
    Temp = temp_grid,
    frac_pred = pred
  )
}

#############################################################################################################
# READ INPUT
#############################################################################################################

msg("Reading input: ", input_file)
proteins <- fread(input_file)

required_cols <- c(
  "Run", "Sample", "Temp", "Replicate", "Protein_group", "Genes",
  "LFQ_L", "LFQ_H", "N_precursors_L", "N_precursors_H"
)
missing_cols <- setdiff(required_cols, names(proteins))
if (length(missing_cols) > 0) {
  stop("Input file is missing required columns: ", paste(missing_cols, collapse = ", "))
}

proteins[, Temp := as.numeric(Temp)]
proteins[, Replicate := as.character(Replicate)]
proteins[, Genes := fifelse(is.na(Genes) | Genes == "", Protein_group, Genes)]
setorder(proteins, Protein_group, Replicate, Temp)

msg("Rows in input: ", nrow(proteins))
msg("Unique proteins: ", uniqueN(proteins$Protein_group))
msg("Temperatures: ", paste(sort(unique(proteins$Temp)), collapse = ", "))
msg("Replicates: ", paste(sort(unique(proteins$Replicate)), collapse = ", "))

#############################################################################################################
# FILTER
#############################################################################################################

proteins_filt <- copy(proteins)

proteins_filt[
  (is.na(N_precursors_L) | N_precursors_L < min_precursors_per_channel) |
    (is.na(N_precursors_H) | N_precursors_H < min_precursors_per_channel),
  c("LFQ_L", "LFQ_H") := .(NA_real_, NA_real_)
]

#############################################################################################################
# LONG TABLE
#############################################################################################################

long <- rbindlist(list(
  proteins_filt[, .(
    Protein_group, Genes, Run, Sample, Temp, Replicate,
    channel = "L", condition = "DMSO",
    abundance = LFQ_L,
    N_precursors = N_precursors_L
  )],
  proteins_filt[, .(
    Protein_group, Genes, Run, Sample, Temp, Replicate,
    channel = "H", condition = "Drug",
    abundance = LFQ_H,
    N_precursors = N_precursors_H
  )]
), use.names = TRUE)

long[abundance <= 0, abundance := NA_real_]
long[, abundance_norm := abundance]

#############################################################################################################
# GLOBAL SAMPLE MEDIAN NORMALIZATION
#############################################################################################################

if (do_global_sample_median_normalization) {
  msg("Applying global sample median normalization...")
  
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
}

#############################################################################################################
# NORMALIZE EACH PROTEIN CURVE TO LOWEST TEMPERATURE
#############################################################################################################

long[, temp_min := suppressWarnings(min(Temp[is.finite(abundance_norm)], na.rm = TRUE)),
     by = .(Protein_group, Replicate, condition)]
long[!is.finite(temp_min), temp_min := NA_real_]

long[, ref_value := mean(abundance_norm[Temp == temp_min], na.rm = TRUE),
     by = .(Protein_group, Replicate, condition)]
long[!is.finite(ref_value) | ref_value <= 0, ref_value := NA_real_]

long[, frac := abundance_norm / ref_value]
long[, frac := pmin(pmax(frac, min_fraction_floor), max_fraction_cap)]

#############################################################################################################
# CURVE COMPLETENESS
#############################################################################################################

curve_stats <- long[, .(
  n_temps = sum(is.finite(frac)),
  min_temp = suppressWarnings(min(Temp[is.finite(frac)], na.rm = TRUE)),
  max_temp = suppressWarnings(max(Temp[is.finite(frac)], na.rm = TRUE))
), by = .(Protein_group, Replicate, condition)]

keep_curves <- curve_stats[n_temps >= min_temps_per_curve]

long <- merge(
  long,
  keep_curves[, .(Protein_group, Replicate, condition)],
  by = c("Protein_group", "Replicate", "condition")
)

#############################################################################################################
# MEAN CURVES
#############################################################################################################

curve_mean <- long[, .(
  mean_frac = safe_mean(frac),
  sd_frac   = safe_sd(frac),
  sem_frac  = safe_sem(frac),
  n_rep     = sum(is.finite(frac))
), by = .(Protein_group, Genes, condition, Temp)]

#############################################################################################################
# DRUG / DMSO TABLES
#############################################################################################################

rep_wide <- dcast(
  long,
  Protein_group + Genes + Replicate + Temp ~ condition,
  value.var = "frac"
)

if (!("DMSO" %in% names(rep_wide))) rep_wide[, DMSO := NA_real_]
if (!("Drug" %in% names(rep_wide))) rep_wide[, Drug := NA_real_]

rep_wide[, log2_ratio_Drug_vs_DMSO := fifelse(
  is.finite(Drug) & is.finite(DMSO) & Drug > 0 & DMSO > 0,
  log2(Drug / DMSO),
  NA_real_
)]
rep_wide[, diff_Drug_minus_DMSO := Drug - DMSO]

mean_wide <- dcast(
  curve_mean,
  Protein_group + Genes + Temp ~ condition,
  value.var = "mean_frac"
)

if (!("DMSO" %in% names(mean_wide))) mean_wide[, DMSO := NA_real_]
if (!("Drug" %in% names(mean_wide))) mean_wide[, Drug := NA_real_]

mean_wide[, log2_ratio_Drug_vs_DMSO := fifelse(
  is.finite(Drug) & is.finite(DMSO) & Drug > 0 & DMSO > 0,
  log2(Drug / DMSO),
  NA_real_
)]
mean_wide[, diff_Drug_minus_DMSO := Drug - DMSO]

#############################################################################################################
# FIT CURVES
#############################################################################################################

fit_results_rep <- data.table()
fit_preds_rep   <- data.table()

if (fit_on_replicates) {
  msg("Fitting replicate-level melt curves...")
  
  fit_groups <- split(
    long[, .(Protein_group, Genes, Replicate, condition, Temp, frac)],
    by = c("Protein_group", "Genes", "Replicate", "condition"),
    keep.by = TRUE
  )
  
  temp_grid <- seq(min(long$Temp, na.rm = TRUE), max(long$Temp, na.rm = TRUE), length.out = 200)
  
  fit_results_rep <- rbindlist(lapply(fit_groups, function(dtg) {
    fr <- fit_sigmoid_curve(dtg)
    data.table(
      Protein_group = unique(dtg$Protein_group),
      Genes         = unique(dtg$Genes),
      Replicate     = unique(dtg$Replicate),
      condition     = unique(dtg$condition),
      fit_success   = fr$success,
      n_points      = fr$n_points,
      top           = fr$top,
      bottom        = fr$bottom,
      Tm            = fr$Tm,
      slope         = fr$slope,
      rss           = fr$rss,
      r2            = fr$r2,
      auc_fit       = fr$auc,
      fit_warning   = fr$warning_msg
    )
  }), fill = TRUE)
  
  fit_preds_rep <- rbindlist(lapply(fit_groups, function(dtg) {
    fr <- fit_sigmoid_curve(dtg)
    pred <- predict_fit_grid(fr, temp_grid)
    if (is.null(pred)) return(NULL)
    pred[, `:=`(
      Protein_group = unique(dtg$Protein_group),
      Genes         = unique(dtg$Genes),
      Replicate     = unique(dtg$Replicate),
      condition     = unique(dtg$condition)
    )]
    pred
  }), fill = TRUE)
}

fit_results_mean <- data.table()
fit_preds_mean   <- data.table()

if (fit_on_mean_curve) {
  msg("Fitting mean melt curves...")
  
  mean_input <- curve_mean[, .(
    Protein_group, Genes, condition, Temp, frac = mean_frac
  )]
  
  fit_groups_mean <- split(
    mean_input,
    by = c("Protein_group", "Genes", "condition"),
    keep.by = TRUE
  )
  
  temp_grid <- seq(min(long$Temp, na.rm = TRUE), max(long$Temp, na.rm = TRUE), length.out = 200)
  
  fit_results_mean <- rbindlist(lapply(fit_groups_mean, function(dtg) {
    fr <- fit_sigmoid_curve(dtg)
    data.table(
      Protein_group = unique(dtg$Protein_group),
      Genes         = unique(dtg$Genes),
      condition     = unique(dtg$condition),
      fit_success   = fr$success,
      n_points      = fr$n_points,
      top           = fr$top,
      bottom        = fr$bottom,
      Tm            = fr$Tm,
      slope         = fr$slope,
      rss           = fr$rss,
      r2            = fr$r2,
      auc_fit       = fr$auc,
      fit_warning   = fr$warning_msg
    )
  }), fill = TRUE)
  
  fit_preds_mean <- rbindlist(lapply(fit_groups_mean, function(dtg) {
    fr <- fit_sigmoid_curve(dtg)
    pred <- predict_fit_grid(fr, temp_grid)
    if (is.null(pred)) return(NULL)
    pred[, `:=`(
      Protein_group = unique(dtg$Protein_group),
      Genes         = unique(dtg$Genes),
      condition     = unique(dtg$condition)
    )]
    pred
  }), fill = TRUE)
}

#############################################################################################################
# Tm / deltaTm SUMMARIES
#############################################################################################################

deltaTm_rep <- data.table()
deltaTm_mean <- data.table()

if (nrow(fit_results_rep) > 0) {
  deltaTm_rep <- dcast(
    fit_results_rep[fit_success == TRUE],
    Protein_group + Genes + Replicate ~ condition,
    value.var = c("Tm", "auc_fit", "r2", "n_points")
  )
  
  if (!("Tm_Drug" %in% names(deltaTm_rep)))   deltaTm_rep[, Tm_Drug := NA_real_]
  if (!("Tm_DMSO" %in% names(deltaTm_rep)))   deltaTm_rep[, Tm_DMSO := NA_real_]
  if (!("auc_fit_Drug" %in% names(deltaTm_rep))) deltaTm_rep[, auc_fit_Drug := NA_real_]
  if (!("auc_fit_DMSO" %in% names(deltaTm_rep))) deltaTm_rep[, auc_fit_DMSO := NA_real_]
  if (!("r2_Drug" %in% names(deltaTm_rep)))   deltaTm_rep[, r2_Drug := NA_real_]
  if (!("r2_DMSO" %in% names(deltaTm_rep)))   deltaTm_rep[, r2_DMSO := NA_real_]
  
  deltaTm_rep[, deltaTm := Tm_Drug - Tm_DMSO]
  deltaTm_rep[, deltaAUC := auc_fit_Drug - auc_fit_DMSO]
}

if (nrow(fit_results_mean) > 0) {
  deltaTm_mean <- dcast(
    fit_results_mean[fit_success == TRUE],
    Protein_group + Genes ~ condition,
    value.var = c("Tm", "auc_fit", "r2", "n_points")
  )
  
  if (!("Tm_Drug" %in% names(deltaTm_mean)))   deltaTm_mean[, Tm_Drug := NA_real_]
  if (!("Tm_DMSO" %in% names(deltaTm_mean)))   deltaTm_mean[, Tm_DMSO := NA_real_]
  if (!("auc_fit_Drug" %in% names(deltaTm_mean))) deltaTm_mean[, auc_fit_Drug := NA_real_]
  if (!("auc_fit_DMSO" %in% names(deltaTm_mean))) deltaTm_mean[, auc_fit_DMSO := NA_real_]
  
  deltaTm_mean[, deltaTm_mean_curve := Tm_Drug - Tm_DMSO]
  deltaTm_mean[, deltaAUC_mean_curve := auc_fit_Drug - auc_fit_DMSO]
}

replicate_summary <- data.table()
if (nrow(deltaTm_rep) > 0) {
  replicate_summary <- deltaTm_rep[, .(
    n_replicates   = .N,
    n_deltaTm      = sum(is.finite(deltaTm)),
    deltaTm_mean   = safe_mean(deltaTm),
    deltaTm_sd     = safe_sd(deltaTm),
    deltaTm_sem    = safe_sem(deltaTm),
    deltaAUC_mean  = safe_mean(deltaAUC),
    deltaAUC_sd    = safe_sd(deltaAUC),
    Tm_DMSO_mean   = safe_mean(Tm_DMSO),
    Tm_Drug_mean   = safe_mean(Tm_Drug),
    r2_DMSO_mean   = safe_mean(r2_DMSO),
    r2_Drug_mean   = safe_mean(r2_Drug)
  ), by = .(Protein_group, Genes)]
}

#############################################################################################################
# PER-TEMPERATURE EFFECT SUMMARY
#############################################################################################################

temp_effect_summary <- rep_wide[, .(
  mean_log2_ratio = safe_mean(log2_ratio_Drug_vs_DMSO),
  sd_log2_ratio   = safe_sd(log2_ratio_Drug_vs_DMSO),
  sem_log2_ratio  = safe_sem(log2_ratio_Drug_vs_DMSO),
  mean_diff       = safe_mean(diff_Drug_minus_DMSO),
  sd_diff         = safe_sd(diff_Drug_minus_DMSO),
  sem_diff        = safe_sem(diff_Drug_minus_DMSO),
  n_rep           = sum(is.finite(log2_ratio_Drug_vs_DMSO))
), by = .(Protein_group, Genes, Temp)]

peak_effect_summary <- temp_effect_summary[, {
  idx_ratio <- safe_whichmax_abs(mean_log2_ratio)
  idx_diff  <- safe_whichmax_abs(mean_diff)
  
  .(
    max_abs_mean_log2_ratio = if (all(!is.finite(mean_log2_ratio))) NA_real_ else max(abs(mean_log2_ratio), na.rm = TRUE),
    max_abs_mean_diff       = if (all(!is.finite(mean_diff))) NA_real_ else max(abs(mean_diff), na.rm = TRUE),
    temp_of_max_abs_ratio   = if (is.na(idx_ratio)) NA_real_ else Temp[idx_ratio],
    temp_of_max_abs_diff    = if (is.na(idx_diff))  NA_real_ else Temp[idx_diff]
  )
}, by = .(Protein_group, Genes)]

#############################################################################################################
# FINAL HIT TABLE
#############################################################################################################

hit_table <- unique(proteins[, .(Protein_group, Genes)])
hit_table <- merge(hit_table, replicate_summary, by = c("Protein_group", "Genes"), all.x = TRUE)
hit_table <- merge(hit_table, deltaTm_mean, by = c("Protein_group", "Genes"), all.x = TRUE)
hit_table <- merge(hit_table, peak_effect_summary, by = c("Protein_group", "Genes"), all.x = TRUE)

curve_completeness <- long[, .(
  n_points_total        = sum(is.finite(frac)),
  n_replicates_present  = uniqueN(Replicate[is.finite(frac)]),
  n_conditions_present  = uniqueN(condition[is.finite(frac)])
), by = .(Protein_group, Genes)]

hit_table <- merge(hit_table, curve_completeness, by = c("Protein_group", "Genes"), all.x = TRUE)

hit_table[, rank_score := abs(deltaTm_mean)]
hit_table[!is.finite(rank_score), rank_score := abs(deltaTm_mean_curve)]
hit_table[!is.finite(rank_score), rank_score := max_abs_mean_log2_ratio]

# setorder cannot use expressions directly, so create helper columns
hit_table[, rank_abs := fifelse(is.finite(rank_score), abs(rank_score), -Inf)]
hit_table[, deltaAUC_abs := fifelse(is.finite(deltaAUC_mean), abs(deltaAUC_mean), -Inf)]
hit_table[, nrep_ord := fifelse(is.finite(n_replicates_present), n_replicates_present, -Inf)]

setorder(hit_table, -rank_abs, -deltaAUC_abs, -nrep_ord)

#############################################################################################################
# QC PLOTS
#############################################################################################################

msg("Generating QC plots...")

p_qc_box <- ggplot(
  long[is.finite(abundance_norm)],
  aes(x = interaction(Temp, Replicate, lex.order = TRUE),
      y = log2(abundance_norm),
      fill = condition)
) +
  geom_boxplot(outlier.size = 0.2) +
  theme_bw(base_size = 10) +
  labs(x = "Temperature_Replicate", y = "log2 normalized abundance", title = "Protein abundance distributions") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p_qc_frac <- ggplot(
  long[is.finite(frac)],
  aes(x = factor(Temp), y = frac, fill = condition)
) +
  geom_boxplot(outlier.size = 0.2) +
  theme_bw(base_size = 10) +
  labs(x = "Temperature", y = "Normalized fraction", title = "Normalized CETSA fractions")

rep_cv <- long[is.finite(frac), .(
  cv_percent = cv_percent(frac)
), by = .(Protein_group, Genes, Temp, condition)]

p_qc_cv <- ggplot(
  rep_cv[is.finite(cv_percent)],
  aes(x = factor(Temp), y = cv_percent, fill = condition)
) +
  geom_boxplot(outlier.size = 0.2) +
  theme_bw(base_size = 10) +
  labs(x = "Temperature", y = "Replicate CV (%)", title = "Replicate CV across temperatures")

make_pca_plot <- function(dt_long, cond_name) {
  dsub <- dt_long[condition == cond_name]
  mat <- dcast(dsub, Protein_group ~ Sample, value.var = "frac")
  mat <- mat[complete.cases(mat)]
  if (nrow(mat) < 10 || ncol(mat) < 3) return(NULL)
  
  mat_ids <- mat$Protein_group
  mat[, Protein_group := NULL]
  X <- t(as.matrix(mat))
  pca <- tryCatch(prcomp(X, scale. = TRUE), error = function(e) NULL)
  if (is.null(pca)) return(NULL)
  
  pca_dt <- data.table(
    Sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  )
  
  ann <- unique(dsub[, .(Sample, Temp, Replicate)])
  pca_dt <- merge(pca_dt, ann, by = "Sample", all.x = TRUE)
  
  var_exp <- summary(pca)$importance["Proportion of Variance", 1:2]
  
  ggplot(pca_dt, aes(PC1, PC2, color = Temp, shape = Replicate)) +
    geom_point(size = 3) +
    theme_bw(base_size = 10) +
    labs(
      title = paste0("PCA - ", cond_name),
      x = paste0("PC1 (", round(100 * var_exp[1], 1), "%)"),
      y = paste0("PC2 (", round(100 * var_exp[2], 1), "%)")
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

heat_dt <- dcast(
  temp_effect_summary,
  Protein_group + Genes ~ Temp,
  value.var = "mean_log2_ratio"
)

if (nrow(heat_dt) > 2) {
  heat_mat <- as.matrix(heat_dt[, -c("Protein_group", "Genes")])
  rownames(heat_mat) <- paste(heat_dt$Genes, heat_dt$Protein_group, sep = " | ")
  
  row_score <- matrixStats::rowMaxs(abs(heat_mat), na.rm = TRUE)
  row_score[!is.finite(row_score)] <- -Inf
  row_order <- order(row_score, decreasing = TRUE)
  
  keep_n <- min(200, nrow(heat_mat))
  heat_mat_top <- heat_mat[row_order[seq_len(keep_n)], , drop = FALSE]
  
  png(file.path(output_dir, "Heatmap_top_mean_log2_Drug_vs_DMSO.png"),
      width = 1800, height = 2200, res = 200)
  pheatmap(
    heat_mat_top,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    main = "Top proteins: mean log2(Drug/DMSO) across temperatures",
    fontsize_row = 6,
    fontsize_col = 10,
    border_color = NA
  )
  dev.off()
}

#############################################################################################################
# PROTEIN CURVE PLOTS
#############################################################################################################

plot_one_protein <- function(protein_id = NULL, gene_name = NULL) {
  if (!is.null(protein_id)) {
    dt_raw  <- long[Protein_group == protein_id]
    dt_mean <- curve_mean[Protein_group == protein_id]
    dt_pred <- fit_preds_mean[Protein_group == protein_id]
    hit_row <- hit_table[Protein_group == protein_id]
  } else if (!is.null(gene_name)) {
    dt_raw  <- long[Genes == gene_name]
    dt_mean <- curve_mean[Genes == gene_name]
    dt_pred <- fit_preds_mean[Genes == gene_name]
    hit_row <- hit_table[Genes == gene_name]
  } else {
    return(NULL)
  }
  
  if (nrow(dt_raw) == 0) return(NULL)
  
  ttl <- unique(dt_raw$Genes)[1]
  subttl <- unique(dt_raw$Protein_group)[1]
  
  extra <- ""
  if (nrow(hit_row) > 0) {
    dTm  <- hit_row$deltaTm_mean[1]
    dTm2 <- hit_row$deltaTm_mean_curve[1]
    dAUC <- hit_row$deltaAUC_mean[1]
    extra <- paste0(
      "deltaTm_rep_mean=",
      ifelse(is.finite(dTm), round(dTm, 2), "NA"),
      " | deltaTm_mean_curve=",
      ifelse(is.finite(dTm2), round(dTm2, 2), "NA"),
      " | deltaAUC=",
      ifelse(is.finite(dAUC), round(dAUC, 2), "NA")
    )
  }
  
  p <- ggplot() +
    geom_point(
      data = dt_raw[is.finite(frac)],
      aes(x = Temp, y = frac, color = condition, shape = Replicate),
      size = 2, alpha = 0.75
    ) +
    geom_line(
      data = dt_raw[is.finite(frac)],
      aes(x = Temp, y = frac, color = condition, group = interaction(condition, Replicate)),
      alpha = 0.45
    ) +
    geom_errorbar(
      data = dt_mean[is.finite(mean_frac)],
      aes(x = Temp, ymin = mean_frac - sem_frac, ymax = mean_frac + sem_frac, color = condition),
      width = 0.3, linewidth = 0.5
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
      title = ttl,
      subtitle = paste0(subttl, "\n", extra),
      x = "Temperature (°C)",
      y = "Normalized soluble fraction"
    ) +
    ylim(0, max_fraction_cap)
  
  if (nrow(dt_pred) > 0) {
    p <- p +
      geom_line(
        data = dt_pred[is.finite(frac_pred)],
        aes(x = Temp, y = frac_pred, color = condition),
        linewidth = 1.2,
        linetype = 2
      )
  }
  
  p
}

if (length(proteins_of_interest) > 0) {
  for (g in proteins_of_interest) {
    p <- plot_one_protein(gene_name = g)
    if (!is.null(p)) {
      ggsave(
        file.path(output_dir, paste0("ProteinCurve_", gsub("[^A-Za-z0-9_\\-]", "_", g), ".png")),
        p, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white"
      )
    }
  }
}

top_hits_for_plots <- hit_table[seq_len(min(top_n_deltaTm_plots, .N))]
if (nrow(top_hits_for_plots) > 0) {
  dir.create(file.path(output_dir, "Top_hit_curves"), showWarnings = FALSE)
  for (i in seq_len(nrow(top_hits_for_plots))) {
    pg <- top_hits_for_plots$Protein_group[i]
    gn <- top_hits_for_plots$Genes[i]
    p <- plot_one_protein(protein_id = pg)
    if (!is.null(p)) {
      fname <- paste0(sprintf("%02d", i), "_", gsub("[^A-Za-z0-9_\\-]", "_", gn), "_", gsub("[^A-Za-z0-9_\\-]", "_", pg), ".png")
      ggsave(
        file.path(output_dir, "Top_hit_curves", fname),
        p, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white"
      )
    }
  }
}

#############################################################################################################
# SUMMARY PLOT
#############################################################################################################

p_rank <- ggplot(
  hit_table[is.finite(deltaTm_mean) | is.finite(deltaAUC_mean)],
  aes(x = deltaTm_mean, y = deltaAUC_mean)
) +
  geom_point(alpha = 0.5) +
  theme_bw(base_size = 11) +
  labs(
    title = "Protein-level CETSA shift summary",
    x = "Mean deltaTm (Drug - DMSO)",
    y = "Mean deltaAUC (Drug - DMSO)"
  )

ggsave(file.path(output_dir, "Summary_deltaTm_vs_deltaAUC.png"), p_rank,
       width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")

#############################################################################################################
# WRITE OUTPUT
#############################################################################################################

fwrite(long,                file.path(output_dir, "long_fraction_table.tsv.gz"), sep = "\t")
fwrite(curve_mean,          file.path(output_dir, "mean_curves.tsv.gz"), sep = "\t")
fwrite(rep_wide,            file.path(output_dir, "replicate_drug_vs_dmso_by_temp.tsv.gz"), sep = "\t")
fwrite(temp_effect_summary, file.path(output_dir, "temp_effect_summary.tsv.gz"), sep = "\t")

if (nrow(fit_results_rep)  > 0) fwrite(fit_results_rep,  file.path(output_dir, "fit_results_replicates.tsv.gz"), sep = "\t")
if (nrow(fit_preds_rep)    > 0) fwrite(fit_preds_rep,    file.path(output_dir, "fit_predictions_replicates.tsv.gz"), sep = "\t")
if (nrow(deltaTm_rep)      > 0) fwrite(deltaTm_rep,      file.path(output_dir, "deltaTm_replicates.tsv.gz"), sep = "\t")

if (nrow(fit_results_mean) > 0) fwrite(fit_results_mean, file.path(output_dir, "fit_results_mean_curves.tsv.gz"), sep = "\t")
if (nrow(fit_preds_mean)   > 0) fwrite(fit_preds_mean,   file.path(output_dir, "fit_predictions_mean_curves.tsv.gz"), sep = "\t")
if (nrow(deltaTm_mean)     > 0) fwrite(deltaTm_mean,     file.path(output_dir, "deltaTm_mean_curves.tsv.gz"), sep = "\t")

fwrite(hit_table,          file.path(output_dir, "CETSA_SILAC_hit_table.tsv.gz"), sep = "\t")

#############################################################################################################
# CONSOLE SUMMARY
#############################################################################################################

msg("Done.")
msg("Output directory: ", output_dir)
msg("Proteins in final hit table: ", nrow(hit_table))
msg("Proteins with replicate-level deltaTm: ", hit_table[is.finite(deltaTm_mean), .N])
msg("Proteins with mean-curve deltaTm: ", hit_table[is.finite(deltaTm_mean_curve), .N])

msg("\nTop proteins by rank score:")
print(
  hit_table[, .(
    Genes,
    Protein_group,
    deltaTm_mean = round(deltaTm_mean, 2),
    deltaTm_mean_curve = round(deltaTm_mean_curve, 2),
    deltaAUC_mean = round(deltaAUC_mean, 3),
    max_abs_mean_log2_ratio = round(max_abs_mean_log2_ratio, 3),
    n_replicates_present
  )][1:min(20, .N)]
)