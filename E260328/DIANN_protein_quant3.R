#############################################################################################################
#
# Transform SILAC-labelled precursor-level output from DIA-NN into protein quantities
# CETSA-MS SILAC version with CETSA-safe normalization
#
# Updated normalization logic, v3
# -------------------------------
# 1. Precursor.Quantity is NOT run-level total-signal normalized for the main analysis.
#    Raw Precursor.Quantity is used for LFQ/intensity-style CETSA melt curves.
#    This avoids artificially scaling high-temperature samples upward, where true loss of soluble
#    protein is expected in CETSA-MS.
#
# 2. Precursor.Quantity.run_norm is calculated only as an optional QC quantity.
#    It is NOT used for MaxLFQ, Intensity, LFQ_L, LFQ_H, fraction_L, or fraction_H.
#
# 3. Precursor.Translated is NOT total-sum normalized before H/L ratio calculation.
#    For SILAC CETSA-MS, H/L is an intra-run ratio and should preserve the DIA-NN translated
#    channel relationship.
#
# 4. Protein-level log2(H/L) ratios are median-centred per Run.
#    This corrects small SILAC mixing/channel bias under the assumption that most proteins are
#    unchanged by compound treatment at a given temperature.
#
# Output columns of interest
# --------------------------
# log2_HL_ratio_raw          = raw protein-level median log2(H/L), H = drug, L = DMSO
# log2_HL_ratio_norm         = run-median-centred log2(H/L), recommended primary CETSA-MS SILAC readout
# log2_HL_ratio              = same as log2_HL_ratio_norm, for downstream compatibility
# LFQ_H                      = Heavy/drug protein abundance from raw Precursor.Quantity
# LFQ_L                      = Light/DMSO protein abundance from raw Precursor.Quantity
# Intensity_H / Intensity_L  = summed raw precursor quantities
# Precursor.Quantity.run_norm = optional QC-only run-normalized precursor quantity
#
#############################################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(iq)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

#############################################################################################################
# USER SETTINGS
#############################################################################################################

input_file <- "report.tsv"
workflow_type <- "second"  # Either "second" or "first"
expected_n_temps <- 11L
expected_replicates_per_temp <- 2L

#############################################################################################################
# LOAD DIA-NN REPORT
#############################################################################################################

report <- fread(
  input_file,
  select = c(
    "Run", "Protein.Group", "Protein.Ids", "First.Protein.Description", "Genes",
    "Stripped.Sequence", "Precursor.Id", "Proteotypic",
    "Precursor.Quantity",
    "Precursor.Translated", "Channel.Q.Value",
    "Q.Value", "Global.Q.Value", "PG.Q.Value", "Global.PG.Q.Value",
    "Lib.Q.Value", "Lib.PG.Q.Value",
    "Precursor.Charge", "RT", "Mass.Evidence", "RT.Start", "RT.Stop"
  )
)

report <- as.data.table(report)

message("Runs detected in report:")
print(report[, .N, by = Run][order(Run)])

#############################################################################################################
# PARSE RUN NAMES, TEMPERATURES, AND REPLICATES
#############################################################################################################

report[, Run := gsub("-", "_", Run)]
report[, Had_timestamp := grepl("_[0-9]{14}$", Run)]
report[, Run_clean := sub("_[0-9]{14}$", "", Run)]

report[, Temp := as.integer(sub(".*_CET_(\\d+)_\\d+$", "\\1", Run_clean))]
report[, Inj := as.integer(sub("^[^_]+_(\\d+)_.*$", "\\1", Run_clean))]

runs_dt <- unique(report[, .(Run, Run_clean, Temp, Inj, Had_timestamp)])
setorder(runs_dt, Temp, Inj)
runs_dt[, Replicate_num := seq_len(.N), by = Temp]

temp_check <- runs_dt[, .N, by = Temp][N != expected_replicates_per_temp]
if (nrow(temp_check) > 0) {
  warning(sprintf(
    "These temperatures do not have exactly %d replicates: %s",
    expected_replicates_per_temp,
    paste(temp_check$Temp, collapse = ", ")
  ))
}

bad_rep <- runs_dt[Replicate_num > expected_replicates_per_temp, unique(Temp)]
if (length(bad_rep) > 0) {
  warning(sprintf(
    "These temperatures have more than %d entries: %s",
    expected_replicates_per_temp,
    paste(bad_rep, collapse = ", ")
  ))
}

n_temps <- uniqueN(runs_dt$Temp)
if (n_temps != expected_n_temps) {
  warning(sprintf(
    "Expected %d temperatures but found %d: %s",
    expected_n_temps,
    n_temps,
    paste(sort(unique(runs_dt$Temp)), collapse = ", ")
  ))
} else {
  message(sprintf(
    "OK: %d temperatures found: %s",
    expected_n_temps,
    paste(sort(unique(runs_dt$Temp)), collapse = ", ")
  ))
}

runs_dt[, Replicate_old := as.integer(sub(".*_(\\d+)$", "\\1", Run_clean))]
runs_dt[, Replicate     := paste0("R", Replicate_num)]
runs_dt[, Sample        := paste0("CET_", Temp, "_", Replicate)]
runs_dt[, Run_base      := sub("_(\\d+)$", "", Run_clean)]
runs_dt[, Run_corrected := paste0(Run_base, "_", Replicate_num)]
runs_dt[, Needs_rename  := (Replicate_old != Replicate_num) | Had_timestamp]

report <- runs_dt[report, on = .(Run, Run_clean, Temp, Inj, Had_timestamp)]

message("Run parsing preview:")
print(unique(report[, .(
  Run, Run_corrected, Temp, Inj, Replicate_old, Replicate, Sample, Needs_rename, Had_timestamp
)]))

rename_table <- runs_dt[Needs_rename == TRUE, .(
  old_run = Run,
  new_run = Run_corrected,
  Temp,
  Inj,
  Had_timestamp
)]

message("Rename table:")
print(rename_table)

#############################################################################################################
# HELPERS
#############################################################################################################

apply_filters <- function(dt, filters) {
  for (f in filters) {
    col <- f$col
    if (col %in% names(dt)) {
      before_n <- nrow(dt)
      dt <- dt[get(col) <= f$thr]
      after_n <- nrow(dt)
      message(sprintf("Filter %s <= %.4g: %d -> %d rows", col, f$thr, before_n, after_n))
    } else {
      message(sprintf("Filter skipped, column missing: %s", col))
    }
  }
  dt
}

CV_function <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  m <- mean(x)
  if (!is.finite(m) || m == 0) return(NA_real_)
  sd(x) / m * 100
}

#############################################################################################################
# FDR FILTERING
#############################################################################################################

if (workflow_type == "second") {
  filters <- list(
    list(col = "Q.Value",         thr = 0.01),
    list(col = "PG.Q.Value",      thr = 0.05),
    list(col = "Lib.Q.Value",     thr = 0.01),
    list(col = "Lib.PG.Q.Value",  thr = 0.01),
    list(col = "Channel.Q.Value", thr = 0.01)
  )
  report <- apply_filters(report, filters)
  message("Second-pass FDR filtering applied")
  
} else if (workflow_type == "first") {
  filters <- list(
    list(col = "Q.Value",           thr = 0.01),
    list(col = "PG.Q.Value",        thr = 0.05),
    list(col = "Global.Q.Value",    thr = 0.01),
    list(col = "Global.PG.Q.Value", thr = 0.01),
    list(col = "Channel.Q.Value",   thr = 0.01)
  )
  report <- apply_filters(report, filters)
  message("Single-pass FDR filtering applied")
  
} else {
  stop("workflow_type must be 'first' or 'second'")
}

#############################################################################################################
# QUANTITY CLEANUP, NO MAIN RUN-LEVEL NORMALIZATION
#############################################################################################################

if (!"Precursor.Quantity" %in% names(report)) stop("Missing column: Precursor.Quantity")
if (!"Precursor.Translated" %in% names(report)) {
  warning("Missing column: Precursor.Translated; SILAC H/L ratios cannot be calculated.")
  report[, Precursor.Translated := NA_real_]
}

report[Precursor.Quantity == 0, Precursor.Quantity := NA_real_]
report[Precursor.Translated == 0, Precursor.Translated := NA_real_]

report[, Precursor.Quantity.raw := Precursor.Quantity]
report[, Precursor.Translated.raw := Precursor.Translated]

# Main analysis uses raw Precursor.Quantity. Remove rows without raw quantity.
report <- report[!is.na(Precursor.Quantity.raw)]

# Optional run-level normalization for QC only.
run_sums <- report[, .(sumQ = sum(Precursor.Quantity.raw, na.rm = TRUE)), by = Run]
target_sum <- run_sums[, median(sumQ[is.finite(sumQ) & sumQ > 0], na.rm = TRUE)]

report[, Precursor.Quantity.run_norm := NA_real_]
if (is.finite(target_sum) && target_sum > 0) {
  report[, run_sumQ := sum(Precursor.Quantity.raw, na.rm = TRUE), by = Run]
  report[is.finite(run_sumQ) & run_sumQ > 0,
         Precursor.Quantity.run_norm := Precursor.Quantity.raw / run_sumQ * target_sum]
  report[, run_sumQ := NULL]
}

message("Precursor.Quantity.raw is used for main LFQ/intensity CETSA analysis.")
message("Precursor.Quantity.run_norm was calculated only for QC.")
message("Precursor.Translated.raw is used for SILAC H/L ratio calculation.")

#############################################################################################################
# PROCESS SILAC LABELS
#############################################################################################################

if (!"Precursor.Id" %in% names(report)) stop("Missing column: Precursor.Id")

report[, label := NA_character_]
report[grepl("SILAC-[A-Z]-L", Precursor.Id), label := "L"]
report[grepl("SILAC-[A-Z]-H", Precursor.Id), label := "H"]

message(sprintf(
  "Precursors with both L and H label, should be 0: %d",
  report[grepl("SILAC-[A-Z]-L", Precursor.Id) & grepl("SILAC-[A-Z]-H", Precursor.Id), .N]
))
message(sprintf("L-labelled precursors, DMSO: %d", report[label == "L", .N]))
message(sprintf("H-labelled precursors, drug: %d", report[label == "H", .N]))
message(sprintf("Unlabelled precursors: %d", report[is.na(label), .N]))

report[, Precursor.Id.nolabels := gsub("SILAC-[A-Z]-[LH]", "SILAC", Precursor.Id)]

#############################################################################################################
# INTRA-RUN SILAC H/L RATIOS, TRACK A INPUT
#############################################################################################################

SILAC <- dcast(
  report[!is.na(label)],
  formula = Precursor.Id.nolabels + Protein.Group + Run ~ label,
  value.var = "Precursor.Translated.raw",
  fun.aggregate = mean,
  na.rm = TRUE
)

SILAC <- SILAC[!is.na(H) & !is.na(L) & H > 0 & L > 0]
SILAC[, precursor_log2_HL_ratio := log2(H / L)]

SILAC <- SILAC[
  is.finite(precursor_log2_HL_ratio),
  .(
    log2_HL_ratio_raw = median(precursor_log2_HL_ratio, na.rm = TRUE),
    N_ratios = .N
  ),
  by = .(Run, Protein.Group)
]

message(sprintf("Protein-run SILAC H/L ratio rows: %d", nrow(SILAC)))

#############################################################################################################
# MAXLFQ ON RAW PRECURSOR.QUANTITY, TRACK B INPUT
#############################################################################################################

LFQ_fun <- function(dt, LFQ_colname) {
  if (nrow(dt) == 0) {
    warning(sprintf("No rows provided for %s", LFQ_colname))
    return(data.table(Protein.Group = character(), Run = character()))
  }
  
  tmp <- dt[
    is.finite(Precursor.Quantity.raw) & Precursor.Quantity.raw > 0,
    .(
      protein_list = Protein.Group,
      sample_list  = Run,
      id           = Precursor.Id,
      quant        = log2(Precursor.Quantity.raw)
    )
  ]
  
  if (nrow(tmp) == 0) {
    warning(sprintf("No finite raw quantities for %s", LFQ_colname))
    return(data.table(Protein.Group = character(), Run = character()))
  }
  
  tmp <- fast_MaxLFQ(tmp)
  
  tmp <- data.table(
    tmp$estimate,
    Precursor_group = tmp$annotation,
    keep.rownames = "Protein.Group"
  )
  
  tmp <- melt(
    tmp,
    id.vars = c("Protein.Group", "Precursor_group"),
    variable.name = "Run",
    value.name = LFQ_colname
  )
  
  tmp <- tmp[Precursor_group == "" & !is.na(get(LFQ_colname))][, -"Precursor_group"]
  tmp[, (LFQ_colname) := 2^get(LFQ_colname)]
  tmp[]
}

LFQ_T <- LFQ_fun(report,               "LFQ")
LFQ_H <- LFQ_fun(report[label == "H"], "LFQ_H")
LFQ_L <- LFQ_fun(report[label == "L"], "LFQ_L")

#############################################################################################################
# SUMMED RAW INTENSITIES AND COUNTS
#############################################################################################################

int_T <- report[,             .(Intensity   = sum(Precursor.Quantity.raw, na.rm = TRUE)), by = .(Protein.Group, Run)]
int_H <- report[label == "H", .(Intensity_H = sum(Precursor.Quantity.raw, na.rm = TRUE)), by = .(Protein.Group, Run)]
int_L <- report[label == "L", .(Intensity_L = sum(Precursor.Quantity.raw, na.rm = TRUE)), by = .(Protein.Group, Run)]

# QC-only run-normalized intensities
int_T_run_norm <- report[,             .(Intensity_run_norm   = sum(Precursor.Quantity.run_norm, na.rm = TRUE)), by = .(Protein.Group, Run)]
int_H_run_norm <- report[label == "H", .(Intensity_H_run_norm = sum(Precursor.Quantity.run_norm, na.rm = TRUE)), by = .(Protein.Group, Run)]
int_L_run_norm <- report[label == "L", .(Intensity_L_run_norm = sum(Precursor.Quantity.run_norm, na.rm = TRUE)), by = .(Protein.Group, Run)]

if (!"Proteotypic" %in% names(report)) report[, Proteotypic := 0L]

counts_H <- report[label == "H",
                   .(
                     N_precursors_H = uniqueN(Precursor.Id),
                     N_precursors_proteotypic_H = sum(as.numeric(Proteotypic), na.rm = TRUE)
                   ),
                   by = .(Protein.Group, Run)
]

counts_L <- report[label == "L",
                   .(
                     N_precursors_L = uniqueN(Precursor.Id),
                     N_precursors_proteotypic_L = sum(as.numeric(Proteotypic), na.rm = TRUE)
                   ),
                   by = .(Protein.Group, Run)
]

#############################################################################################################
# ANNOTATIONS
#############################################################################################################

ann_cols <- c(
  "Protein.Group", "Run", "Sample", "Temp", "Replicate",
  "First.Protein.Description", "Genes"
)
ann_cols <- ann_cols[ann_cols %in% names(report)]
annotations <- unique(report[, ..ann_cols])

#############################################################################################################
# MERGE PROTEIN OUTPUTS
#############################################################################################################

merge_function <- function(x, y) merge(x, y, by = c("Protein.Group", "Run"), all = TRUE)

proteins <- Reduce(merge_function, list(
  LFQ_T, LFQ_H, LFQ_L,
  int_T, int_H, int_L,
  int_T_run_norm, int_H_run_norm, int_L_run_norm,
  counts_H, counts_L,
  SILAC,
  annotations
))

proteins <- proteins[!is.na(Intensity)]

if ("First.Protein.Description" %in% names(proteins)) {
  setnames(
    proteins,
    old = c("Protein.Group", "First.Protein.Description"),
    new = c("Protein_group", "Protein_description")
  )
} else {
  setnames(proteins, old = "Protein.Group", new = "Protein_group")
}

#############################################################################################################
# RUN-MEDIAN NORMALIZATION OF PROTEIN-LEVEL H/L RATIOS
#############################################################################################################

if ("log2_HL_ratio_raw" %in% names(proteins)) {
  proteins[, log2_HL_median_shift := median(log2_HL_ratio_raw, na.rm = TRUE), by = Run]
  proteins[, log2_HL_ratio_norm := log2_HL_ratio_raw - log2_HL_median_shift]
  proteins[, log2_HL_ratio := log2_HL_ratio_norm]
  
  message("Added log2_HL_ratio_raw, log2_HL_median_shift, and log2_HL_ratio_norm.")
  message("Column log2_HL_ratio is set equal to log2_HL_ratio_norm for downstream compatibility.")
} else {
  warning("No SILAC ratio column found; skipping H/L median normalization.")
}

#############################################################################################################
# OPTIONAL CURVE-FRIENDLY FRACTION NON-DENATURED VALUES
#############################################################################################################

if (all(c("LFQ_L", "LFQ_H", "Temp", "Replicate") %in% names(proteins))) {
  proteins[, Temp_numeric := as.numeric(as.character(Temp))]
  
  proteins[, LFQ_L_baseline := {
    x <- LFQ_L[Temp_numeric == min(Temp_numeric, na.rm = TRUE)]
    if (length(x) == 0 || all(!is.finite(x))) NA_real_ else median(x, na.rm = TRUE)
  }, by = .(Protein_group, Replicate)]
  
  proteins[, LFQ_H_baseline := {
    x <- LFQ_H[Temp_numeric == min(Temp_numeric, na.rm = TRUE)]
    if (length(x) == 0 || all(!is.finite(x))) NA_real_ else median(x, na.rm = TRUE)
  }, by = .(Protein_group, Replicate)]
  
  proteins[, fraction_L := fifelse(
    is.finite(LFQ_L) & is.finite(LFQ_L_baseline) & LFQ_L_baseline > 0,
    LFQ_L / LFQ_L_baseline, NA_real_
  )]
  
  proteins[, fraction_H := fifelse(
    is.finite(LFQ_H) & is.finite(LFQ_H_baseline) & LFQ_H_baseline > 0,
    LFQ_H / LFQ_H_baseline, NA_real_
  )]
  
  proteins[, delta_fraction_H_minus_L := fraction_H - fraction_L]
  
  message("Added fraction_L, fraction_H, and delta_fraction_H_minus_L for curve inspection.")
}

#############################################################################################################
# ORDER COLUMNS AND WRITE OUTPUT
#############################################################################################################

preferred_first <- intersect(
  c(
    "Run", "Sample", "Temp", "Temp_numeric", "Replicate",
    "Protein_group", "Protein_description", "Genes",
    "LFQ", "LFQ_H", "LFQ_L",
    "LFQ_H_baseline", "LFQ_L_baseline", "fraction_H", "fraction_L", "delta_fraction_H_minus_L",
    "Intensity", "Intensity_H", "Intensity_L",
    "Intensity_run_norm", "Intensity_H_run_norm", "Intensity_L_run_norm",
    "log2_HL_ratio", "log2_HL_ratio_norm", "log2_HL_ratio_raw", "log2_HL_median_shift", "N_ratios",
    "N_precursors_H", "N_precursors_L",
    "N_precursors_proteotypic_H", "N_precursors_proteotypic_L"
  ),
  names(proteins)
)
setcolorder(proteins, c(preferred_first, setdiff(names(proteins), preferred_first)))

fwrite(proteins, "protein_quant.csv.gz")
fwrite(report, "reportfile.tsv")

message("Wrote protein_quant.csv.gz")
message("Wrote reportfile.tsv")

#############################################################################################################
# QC PLOTTING SECTION
#############################################################################################################

if (!"Temp_numeric" %in% names(proteins)) {
  proteins[, Temp_numeric := as.numeric(as.character(Temp))]
}
report[, Temp_numeric := as.numeric(as.character(Temp))]

desired_temp_order <- sort(unique(proteins$Temp_numeric[is.finite(proteins$Temp_numeric)]))
desired_replicate_order <- paste0("R", seq_len(expected_replicates_per_temp))

temp_levels <- as.character(desired_temp_order)
temp_color_mapping <- setNames(colorRampPalette(c("#2166ac", "#d73027"))(length(temp_levels)), temp_levels)

#############################################################################################################
# CV CALCULATIONS BEFORE FACTOR CONVERSION
#############################################################################################################

precursor_by_sample <- report[
  !is.na(Precursor.Quantity.raw) & !is.na(label) & !is.na(Temp_numeric),
  .(qty = median(Precursor.Quantity.raw, na.rm = TRUE)),
  by = .(Precursor.Id.nolabels, label, Temp_numeric, Sample)
]

precursor_CVs <- precursor_by_sample[
  ,
  .(PrecursorCV = CV_function(qty), N_samples = .N),
  by = .(Precursor.Id.nolabels, label, Temp_numeric)
][N_samples >= 2 & is.finite(PrecursorCV)]

precursor_CVs[, precursor_rank := frank(PrecursorCV), by = .(Temp_numeric, label)]
precursor_CVs[, label := factor(label, levels = c("L", "H"))]
precursor_CVs[, Temp_char := as.character(Temp_numeric)]

protein_long <- melt(
  proteins,
  id.vars = c("Protein_group", "Temp_numeric", "Sample"),
  measure.vars = intersect(c("LFQ_L", "LFQ_H"), names(proteins)),
  variable.name = "label",
  value.name = "LFQ_intensity",
  na.rm = TRUE
)
protein_long[, label := gsub("LFQ_", "", label)]

protein_CVs <- protein_long[
  ,
  .(ProteinCV = CV_function(LFQ_intensity), N_samples = uniqueN(Sample)),
  by = .(Protein_group, label, Temp_numeric)
][N_samples >= 2 & is.finite(ProteinCV)]

protein_CVs[, protein_rank := frank(ProteinCV), by = .(Temp_numeric, label)]
protein_CVs[, label := factor(label, levels = c("L", "H"))]
protein_CVs[, Temp_char := as.character(Temp_numeric)]

message(sprintf(
  "precursor_CVs rows: %d, label levels: %s",
  nrow(precursor_CVs), paste(levels(precursor_CVs$label), collapse = ", ")
))

#############################################################################################################
# FACTORS FOR PLOTTING
#############################################################################################################

report[, Temp := factor(Temp_numeric, levels = desired_temp_order)]
report[, Replicate := factor(Replicate, levels = desired_replicate_order)]
proteins[, Temp := factor(Temp_numeric, levels = desired_temp_order)]
proteins[, Replicate := factor(Replicate, levels = desired_replicate_order)]

sample_levels <- proteins[order(Temp_numeric, Replicate), unique(Sample)]
report[, Sample := factor(Sample, levels = sample_levels)]
proteins[, Sample := factor(Sample, levels = sample_levels)]

#############################################################################################################
# THEME AND COLOURS
#############################################################################################################

theme_set(
  theme(
    line = element_line(linewidth = 0.25),
    text = element_text(size = 7),
    rect = element_rect(colour = "black", linewidth = 0.25, fill = NA),
    panel.border = element_rect(linewidth = 0.25),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(colour = "black", size = 6),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key.width = unit(4, "mm"),
    legend.key = element_rect(colour = NA),
    strip.background = element_blank(),
    strip.text = element_text(vjust = 1)
  )
)

HL_cols <- c(
  heavy = "#37abc8ff",
  H = "#37abc8ff",
  light = "#ca00a3ff",
  L = "#ca00a3ff",
  total = "black"
)
scale_HL_labels <- c(L = "DMSO (light)", H = "Drug (heavy)")

#############################################################################################################
# QC PLOT 1: INTENSITY DISTRIBUTIONS
#############################################################################################################

plot_int_DT <- rbindlist(list(
  report[!is.na(Precursor.Quantity.raw),
         .(Sample, `log2 quantity` = log2(pmax(Precursor.Quantity.raw, 1)), Type = "Precursor quantity raw")],
  
  report[!is.na(Precursor.Quantity.run_norm),
         .(Sample, `log2 quantity` = log2(pmax(Precursor.Quantity.run_norm, 1)), Type = "Precursor quantity run-normalized QC")],
  
  proteins[!is.na(LFQ),
           .(Sample, `log2 quantity` = log2(pmax(LFQ, 1)), Type = "Protein MaxLFQ raw")],
  
  proteins[!is.na(Intensity),
           .(Sample, `log2 quantity` = log2(pmax(Intensity, 1)), Type = "Protein intensity raw")]
), use.names = TRUE, fill = TRUE)

plot_int_DT[, Type := factor(Type, levels = c(
  "Precursor quantity raw",
  "Precursor quantity run-normalized QC",
  "Protein MaxLFQ raw",
  "Protein intensity raw"
))]

pInt <- ggplot(plot_int_DT, aes(x = Sample, y = `log2 quantity`)) +
  facet_wrap(~Type, scales = "free_x") +
  geom_boxplot(outliers = FALSE, linewidth = 0.25) +
  theme(axis.title.y = element_blank()) +
  coord_flip()

plot(pInt)

#############################################################################################################
# QC PLOT 2: CUMULATIVE RAW INTENSITY DISTRIBUTIONS
#############################################################################################################

cum_cols <- intersect(c("LFQ_H", "LFQ_L", "Intensity_H", "Intensity_L"), names(proteins))

plot_int <- proteins[, lapply(.SD, function(x) sum(x, na.rm = TRUE)), by = Sample, .SDcols = cum_cols]
plot_int <- melt(plot_int, id.vars = "Sample", value.name = "Cumulative quantity", variable.name = "label")

plot_int[, Type := fifelse(label %like% "LFQ", "Protein MaxLFQ raw", "Protein intensity raw")]
plot_int[, label := gsub(".+_", "", label)]
plot_int[, Type := factor(Type, levels = c("Protein MaxLFQ raw", "Protein intensity raw"))]

pCumInt <- ggplot(plot_int, aes(x = Sample, y = `Cumulative quantity`, fill = label)) +
  facet_wrap(~Type, scales = "free_x") +
  geom_bar(stat = "identity") +
  ylab("Cumulative raw protein intensity") +
  scale_fill_manual(values = HL_cols, guide = guide_legend(title = "Channel"), labels = scale_HL_labels) +
  theme(axis.title.y = element_blank()) +
  coord_flip()

plot(pCumInt)

#############################################################################################################
# QC PLOT 3: COUNTS
#############################################################################################################

count_cols <- intersect(c("Precursor.Id", "Stripped.Sequence", "Protein.Group"), names(report))

counts <- report[!is.na(label), lapply(.SD, uniqueN), by = .(Sample, label), .SDcols = count_cols]
counts <- melt(counts, id.vars = c("Sample", "label"), value.name = "IDs")

facet_labels <- as_labeller(c(
  Precursor.Id = "Precursors",
  Stripped.Sequence = "Peptides",
  Protein.Group = "Protein groups"
))

pCount <- ggplot(counts, aes(x = Sample, y = IDs / 1000, fill = label,
                             label = format(IDs, big.mark = ",", scientific = FALSE))) +
  facet_wrap(~variable, scales = "free_x", labeller = facet_labels) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(size = 2, hjust = 1.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = HL_cols, guide = guide_legend(title = "Channel"), labels = scale_HL_labels) +
  ylab("Number of IDs [x1,000]") +
  theme(axis.title.y = element_blank()) +
  coord_flip()

plot(pCount)

#############################################################################################################
# QC PLOT 4: DATA COMPLETENESS
#############################################################################################################

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
  f_prot_count(proteins[!is.na(N_precursors_L) & !is.na(N_precursors_H) & N_precursors_L >= 2 & N_precursors_H >= 2],
               "min 2 DMSO &\nmin 2 Drug")
), use.names = TRUE, fill = TRUE)

pNAs <- ggplot(data_completeness, aes(x = N_samples, y = cumulative_protein_N / 1000, colour = Precursors)) +
  geom_point(size = 0.8) +
  geom_line(linewidth = 0.25) +
  xlab("Number of samples") +
  ylab("Number of proteins [x1,000]") +
  scale_x_continuous(breaks = seq(1, max(data_completeness$N_samples), 1)) +
  scale_y_continuous(limits = c(0, max(data_completeness$cumulative_protein_N, na.rm = TRUE) / 1000)) +
  scale_color_brewer(palette = "Dark2") +
  theme(
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
    panel.grid.minor = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted")
  )

plot(pNAs)

#############################################################################################################
# QC PLOT 5A: PRECURSOR CVs, RAW QUANTITY
#############################################################################################################

pCV1 <- ggplot(precursor_CVs, aes(x = precursor_rank / 1000, y = PrecursorCV, colour = Temp_char)) +
  facet_wrap(~label, labeller = as_labeller(c(L = "DMSO (light)", H = "Drug (heavy)"))) +
  geom_line(linewidth = 0.25) +
  labs(x = "Number of precursors [x1,000]", y = "Precursor CV, raw quantity [%]") +
  coord_cartesian(ylim = c(0, 50)) +
  scale_color_manual(values = temp_color_mapping, guide = guide_legend(title = "Temp (°C)")) +
  theme(
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
    panel.grid.minor = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
    legend.position = "right"
  )

plot(pCV1)

#############################################################################################################
# QC PLOT 5B: PROTEIN CVs, RAW LFQ
#############################################################################################################

pCV2 <- ggplot(protein_CVs, aes(x = protein_rank / 1000, y = ProteinCV, colour = Temp_char)) +
  facet_wrap(~label, labeller = as_labeller(c(L = "DMSO (light)", H = "Drug (heavy)"))) +
  geom_line(linewidth = 0.25) +
  labs(x = "Number of proteins [x1,000]", y = "Protein CV, raw LFQ [%]") +
  coord_cartesian(ylim = c(0, 50)) +
  scale_color_manual(values = temp_color_mapping, guide = guide_legend(title = "Temp (°C)")) +
  theme(
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
    panel.grid.minor = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
    legend.position = "right"
  )

plot(pCV2)

#############################################################################################################
# QC PLOT 6: PCA
#############################################################################################################

run_pca_plot <- function(proteins_dt, value_col, title_text, temp_color_mapping) {
  if (!value_col %in% names(proteins_dt)) {
    warning(sprintf("Skipping PCA: missing %s", value_col))
    return(NULL)
  }
  
  PCA_input <- dcast(proteins_dt, Protein_group ~ Sample, value.var = value_col)
  PCA_input <- PCA_input[complete.cases(PCA_input)]
  
  if (nrow(PCA_input) < 3 || ncol(PCA_input) < 4) {
    warning(sprintf("Skipping PCA for %s: not enough complete data", value_col))
    return(NULL)
  }
  
  mat <- t(data.frame(PCA_input, row.names = "Protein_group"))
  DT_PCA <- prcomp(mat, scale. = TRUE)
  var_explained <- summary(DT_PCA)$importance
  DT_PCA <- data.table(DT_PCA$x, keep.rownames = "Sample")
  DT_PCA <- merge(DT_PCA, unique(proteins_dt[, .(Sample, Temp, Temp_numeric, Replicate)]), by = "Sample", all.x = TRUE)
  DT_PCA[, Temp_char := as.character(Temp_numeric)]
  
  ggplot(DT_PCA, aes(x = PC1, y = PC2, colour = Temp_char, shape = Replicate)) +
    ggtitle(title_text) +
    geom_point(size = 2.5) +
    xlab(paste0("PC1 [", round(var_explained["Proportion of Variance", "PC1"] * 100, 0), "%]")) +
    ylab(paste0("PC2 [", round(var_explained["Proportion of Variance", "PC2"] * 100, 0), "%]")) +
    scale_color_manual(values = temp_color_mapping, guide = guide_legend(title = "Temp (°C)")) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right"
    )
}

pPCA_L <- run_pca_plot(proteins, "LFQ_L", "DMSO control (Light), raw LFQ", temp_color_mapping)
pPCA_H <- run_pca_plot(proteins, "LFQ_H", "Drug treated (Heavy), raw LFQ", temp_color_mapping)

if (!is.null(pPCA_L)) plot(pPCA_L)
if (!is.null(pPCA_H)) plot(pPCA_H)

#############################################################################################################
# QC PLOT 7: LOG2(H/L) RATIO DISTRIBUTIONS
#############################################################################################################

pHL_raw <- ggplot(proteins[!is.na(log2_HL_ratio_raw)],
                  aes(x = Temp, y = log2_HL_ratio_raw, fill = as.character(Temp_numeric))) +
  geom_boxplot(outliers = FALSE, linewidth = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey40") +
  scale_fill_manual(values = temp_color_mapping, guide = "none") +
  labs(x = "Temperature (°C)", y = "raw log2(Drug/DMSO)",
       title = "Raw H/L ratio distribution per temperature") +
  theme(panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"))

pHL_norm <- ggplot(proteins[!is.na(log2_HL_ratio_norm)],
                   aes(x = Temp, y = log2_HL_ratio_norm, fill = as.character(Temp_numeric))) +
  geom_boxplot(outliers = FALSE, linewidth = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey40") +
  scale_fill_manual(values = temp_color_mapping, guide = "none") +
  labs(x = "Temperature (°C)", y = "median-centred log2(Drug/DMSO)",
       title = "Normalized H/L ratio distribution per temperature") +
  theme(panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"))

plot(pHL_raw)
plot(pHL_norm)

#############################################################################################################
# QC PLOT 8: RUN MEDIAN SHIFT
#############################################################################################################

pHL_shift <- ggplot(unique(proteins[!is.na(log2_HL_median_shift), .(Sample, Temp, Temp_numeric, Replicate, log2_HL_median_shift)]),
                    aes(x = Sample, y = log2_HL_median_shift, fill = as.character(Temp_numeric))) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey40") +
  scale_fill_manual(values = temp_color_mapping, guide = guide_legend(title = "Temp (°C)")) +
  labs(x = NULL, y = "run median raw log2(H/L)", title = "SILAC H/L median shift removed during normalization") +
  coord_flip()

plot(pHL_shift)

#############################################################################################################
# QC PLOT 9: FRACTION NON-DENATURED DISTRIBUTIONS, RAW LFQ
#############################################################################################################

fraction_long <- melt(
  proteins,
  id.vars = c("Protein_group", "Sample", "Temp", "Temp_numeric", "Replicate"),
  measure.vars = intersect(c("fraction_L", "fraction_H"), names(proteins)),
  variable.name = "Channel",
  value.name = "Fraction",
  na.rm = TRUE
)
fraction_long[, Channel := fifelse(Channel == "fraction_L", "DMSO (light)", "Drug (heavy)")]

pFrac <- ggplot(fraction_long[is.finite(Fraction) & Fraction > 0],
                aes(x = Temp, y = log2(Fraction), fill = Channel)) +
  geom_boxplot(outliers = FALSE, linewidth = 0.25, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey40") +
  labs(x = "Temperature (°C)", y = "log2 fraction relative to lowest temperature",
       title = "Protein melt-curve normalization check, raw LFQ") +
  scale_fill_manual(values = c("DMSO (light)" = "#ca00a3ff", "Drug (heavy)" = "#37abc8ff")) +
  theme(panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"))

plot(pFrac)

#############################################################################################################
# COMBINE AND SAVE QC PLOTS
#############################################################################################################

pPCA_combined <- if (!is.null(pPCA_L) && !is.null(pPCA_H)) {
  pPCA_L + pPCA_H
} else if (!is.null(pPCA_L)) {
  pPCA_L
} else if (!is.null(pPCA_H)) {
  pPCA_H
} else {
  plot_spacer()
}

combined_plot <- (
  pInt /
    pCumInt /
    pCount /
    patchwork::wrap_elements(pNAs + plot_spacer()) /
    patchwork::wrap_elements(pCV1 + pCV2) /
    patchwork::wrap_elements(pPCA_combined) /
    pHL_raw /
    pHL_norm /
    pHL_shift /
    pFrac
) +
  plot_layout(ncol = 1, heights = c(1, 1, 1.8, 1, 1, 1, 0.8, 0.8, 0.8, 0.8)) +
  plot_annotation(theme = theme(plot.margin = margin(20, 20, 20, 20)))

plot(combined_plot)

ggsave(
  "DIANNprotein_quant_plot_LARGE.png",
  combined_plot,
  width = 35,
  height = 110,
  units = "cm",
  dpi = 300,
  bg = "white"
)

message("Wrote DIANNprotein_quant_plot_LARGE.png")

#############################################################################################################
# FINAL NOTES
#############################################################################################################

message("Done.")
message("Recommended downstream CETSA-MS SILAC signal: log2_HL_ratio_norm")
message("Use raw LFQ_L/LFQ_H or raw-derived fraction_L/fraction_H as secondary melt-curve evidence.")
message("Precursor.Quantity.run_norm and run-normalized intensities are QC-only, not used for downstream CETSA statistics.")
