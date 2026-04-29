############################################################################################################# #
#
# Transform SILAC-labelled precursor-level output from DIA-NN into protein quantities
# Date: Updated Feb 2025
############################################################################################################# #
# Load required libraries
library(data.table); library(iq); library(ggplot2); library(patchwork); library(scales)
#### Load the report and parameterise the script ####
# Load the precursor report (relevant columns only)
report <- fread("report.tsv",  # Default is report.tsv.gz (second-pass/MBR). Alternatively use report-first-pass.tsv.gz
                select = c("Run", "Protein.Group", "Protein.Ids", "First.Protein.Description", "Genes", "Stripped.Sequence", "Precursor.Id", "Proteotypic", 
                           "Precursor.Quantity",    # Note: we use unnormalised quantities for SILAC data
                           "Precursor.Translated", "Channel.Q.Value",  # These are SILAC-specific columns
                           "Q.Value", "Global.Q.Value", "PG.Q.Value", "Global.PG.Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value", "Precursor.Charge", "RT", "Mass.Evidence", "RT.Start", "RT.Stop"))
# Is this a first- or second-pass report? Required for FDR filtering.
workflow_type <- "second"  # Either "second" (default) or "first"
# Show full sample names
report[, .N, Run]
library(data.table)
report <- as.data.table(report)
# clean names: remove dashes, flag and strip trailing timestamp (14 digits)
report[, Run := gsub("-", "_", Run)]
report[, Had_timestamp := grepl("_[0-9]{14}$", Run)]
report[, Run_clean := sub("_[0-9]{14}$", "", Run)]
# extract temperature
report[, Temp := as.integer(sub(".*_CET_(\\d+)_\\d+$", "\\1", Run_clean))]
# extract injection number from the leading numeric field (e.g. E260328_06_... -> 06)
report[, Inj := as.integer(sub("^[^_]+_(\\d+)_.*$", "\\1", Run_clean))]
# ── assign replicates at run level, not row level ──────────────────────────────
runs_dt <- unique(report[, .(Run, Run_clean, Temp, Inj, Had_timestamp)])
setorder(runs_dt, Temp, Inj)
runs_dt[, Replicate_num := seq_len(.N), by = Temp]
# checks
temp_check <- runs_dt[, .N, by = Temp][N != 2]
if (nrow(temp_check) > 0) {
  warning(
    sprintf(
      "These temperatures do not have exactly 2 replicates: %s",
      paste(temp_check$Temp, collapse = ", ")
    )
  )
}
bad_rep <- runs_dt[Replicate_num > 2, unique(Temp)]
if (length(bad_rep) > 0) {
  warning(
    sprintf(
      "These temperatures have more than 2 entries: %s",
      paste(bad_rep, collapse = ", ")
    )
  )
}
# verify exactly 11 temperatures
n_temps <- uniqueN(runs_dt$Temp)
if (n_temps != 11) {
  warning(sprintf("Expected 11 temperatures but found %d: %s", n_temps, paste(sort(unique(runs_dt$Temp)), collapse = ", ")))
} else {
  message(sprintf("OK: 11 temperatures found: %s", paste(sort(unique(runs_dt$Temp)), collapse = ", ")))
}
# old replicate from original file name
runs_dt[, Replicate_old := as.integer(sub(".*_(\\d+)$", "\\1", Run_clean))]
# new names
runs_dt[, Replicate        := paste0("R", Replicate_num)]
runs_dt[, Sample           := paste0("CET_", Temp, "_", Replicate)]
runs_dt[, Run_base         := sub("_(\\d+)$", "", Run_clean)]
runs_dt[, Run_corrected    := paste0(Run_base, "_", Replicate_num)]
# flag any run where original replicate suffix differs from corrected, OR timestamp was stripped
runs_dt[, Needs_rename := (Replicate_old != Replicate_num) | Had_timestamp]
# join derived columns back to the full report
report <- runs_dt[report, on = .(Run, Run_clean, Temp, Inj, Had_timestamp)]
# preview
report[, .(
  Run,
  Run_corrected,
  Temp,
  Inj,
  Replicate_old,
  Replicate,
  Sample,
  Needs_rename,
  Had_timestamp
)]
# rename table only
rename_table <- runs_dt[Needs_rename == TRUE, .(
  old_run = Run,
  new_run = Run_corrected,
  Temp,
  Inj,
  Had_timestamp
)]
rename_table

# NB: Code below here can usually be executed without further adjustment

#### Data processing section ####

# ----------------------------
# Helper: safe filter only if column exists
# ----------------------------
apply_filters <- function(dt, filters) {
  for (f in filters) {
    col <- f$col
    if (col %in% names(dt)) {
      dt <- dt[ get(col) <= f$thr ]
    }
  }
  dt
}

# ----------------------------
# Apply FDR filters dynamically based on workflow type
# ----------------------------
if (workflow_type == "second") {
  filters <- list(
    list(col = "Q.Value",          thr = 0.01),
    list(col = "PG.Q.Value",       thr = 0.05),
    list(col = "Lib.Q.Value",      thr = 0.01),
    list(col = "Lib.PG.Q.Value",   thr = 0.01),
    list(col = "Channel.Q.Value",  thr = 0.01)
  )
  report <- apply_filters(report, filters)
  message("Second-pass FDR filtering applied")
  
} else if (workflow_type == "first") {
  filters <- list(
    list(col = "Q.Value",          thr = 0.01),
    list(col = "PG.Q.Value",       thr = 0.05),
    list(col = "Global.Q.Value",   thr = 0.01),
    list(col = "Global.PG.Q.Value",thr = 0.01),
    list(col = "Channel.Q.Value",  thr = 0.01)
  )
  report <- apply_filters(report, filters)
  message("Single-pass FDR filtering applied")
} else {
  stop("workflow_type must be 'first' or 'second'")
}

# ----------------------------
# Convert zeroes to NAs for precursor quantities
# ----------------------------
if ("Precursor.Quantity" %in% names(report)) {
  report[ Precursor.Quantity == 0, Precursor.Quantity := NA_real_ ]
} else stop("Missing column: Precursor.Quantity")

if ("Precursor.Translated" %in% names(report)) {
  report[ Precursor.Translated == 0, Precursor.Translated := NA_real_ ]
} else {
  report[, Precursor.Translated := NA_real_ ]
}

# Remove precursors that were not quantified at all (NAs are ok for translated values)
report <- report[ !is.na(Precursor.Quantity) ]

# ----------------------------
# Normalise Precursor.Quantity to adjust for unequal sample loading
# ----------------------------
target_sum <- report[, .(sumQ = sum(Precursor.Quantity, na.rm = TRUE)), by = Run][, median(sumQ)]

report[, Precursor.Quantity :=
         Precursor.Quantity / sum(Precursor.Quantity, na.rm = TRUE) * target_sum,
       by = Run]

report[, Precursor.Translated :=
         Precursor.Translated / sum(Precursor.Translated, na.rm = TRUE) * target_sum,
       by = Run]

# ----------------------------
# Process SILAC labels (L = DMSO control, H = Drug treated)
# Precursor.Id format: PEPTIDESEQ(SILAC-R-L)2 or PEPTIDESEQ(SILAC-K-H)2
# Use grepl with regex (not %like% which uses SQL LIKE pattern matching)
# ----------------------------
if (!("Precursor.Id" %in% names(report))) stop("Missing column: Precursor.Id")

report[ grepl("SILAC-[A-Z]-L", Precursor.Id), label := "L" ]   # L = DMSO (control)
report[ grepl("SILAC-[A-Z]-H", Precursor.Id), label := "H" ]   # H = Drug (treated)

# Control: should be zero (no precursors with both labels)
message(sprintf("Precursors with both L and H label (should be 0): %d",
                report[ grepl("SILAC-[A-Z]-L", Precursor.Id) & grepl("SILAC-[A-Z]-H", Precursor.Id), .N ]))

# Check labelling worked
message(sprintf("L-labelled precursors (DMSO): %d", report[label == "L", .N]))
message(sprintf("H-labelled precursors (Drug): %d", report[label == "H", .N]))
message(sprintf("Unlabelled precursors:        %d", report[is.na(label), .N]))

report[, Precursor.Id.nolabels := gsub("SILAC-[A-Z]-[LH]", "SILAC", Precursor.Id)]

# ----------------------------
# Intra-sample H/L ratios: Drug/Control ratio using translated quantities
# This is the key CETSA-MS readout: thermal stability shift upon drug treatment
# ----------------------------
SILAC <- dcast(
  report[ !is.na(label) ],        # exclude unlabelled precursors to avoid NA column name
  formula       = Precursor.Id.nolabels + Protein.Group + Run ~ label,
  value.var     = "Precursor.Translated",
  fun.aggregate = mean,           # average any within-group duplicates
  na.rm         = TRUE
)

# Avoid log2(H/L) when missing or zero
SILAC <- SILAC[ !is.na(H) & !is.na(L) & H > 0 & L > 0 ]
SILAC[, precursor_log2_HL_ratio := log2(H / L)]   # log2(Drug/DMSO)

SILAC <- SILAC[
  ,
  .(log2_HL_ratio = median(precursor_log2_HL_ratio, na.rm = TRUE),
    N_ratios      = .N),
  by = .(Run, Protein.Group)
]

# ----------------------------
# Cross-sample comparisons: MaxLFQ on non-translated quantities
# Run separately for all precursors, H-only (Drug), and L-only (DMSO)
# ----------------------------
LFQ_fun <- function(dt, LFQ_colname) {
  tmp <- dt[, .(
    protein_list = Protein.Group,
    sample_list  = Run,
    id           = Precursor.Id,
    quant        = log2(Precursor.Quantity)
  )]
  
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
  
  tmp <- tmp[ Precursor_group == "" & !is.na(get(LFQ_colname)) ][, -"Precursor_group"]
  
  tmp[, (LFQ_colname) := 2^get(LFQ_colname)]
  
  return(tmp)
}

LFQ_T <- LFQ_fun(report,                "LFQ")    # all precursors
LFQ_H <- LFQ_fun(report[label == "H"], "LFQ_H")  # Drug (Heavy)
LFQ_L <- LFQ_fun(report[label == "L"], "LFQ_L")  # DMSO (Light)

# ----------------------------
# Crude "total signal": sum precursor quantities per protein+run
# ----------------------------
int_T <- report[,             .(Intensity   = sum(Precursor.Quantity, na.rm = TRUE)), by = .(Protein.Group, Run)]
int_H <- report[label == "H", .(Intensity_H = sum(Precursor.Quantity, na.rm = TRUE)), by = .(Protein.Group, Run)]  # Drug
int_L <- report[label == "L", .(Intensity_L = sum(Precursor.Quantity, na.rm = TRUE)), by = .(Protein.Group, Run)]  # DMSO

# ----------------------------
# Precursor counts per protein+run
# ----------------------------
if (!("Proteotypic" %in% names(report))) report[, Proteotypic := 0L]

counts_H <- report[label == "H",
                   .(N_precursors_H = uniqueN(Precursor.Id),
                     N_precursors_proteotypic_H = sum(Proteotypic, na.rm = TRUE)),
                   by = .(Protein.Group, Run)
]

counts_L <- report[label == "L",
                   .(N_precursors_L = uniqueN(Precursor.Id),
                     N_precursors_proteotypic_L = sum(Proteotypic, na.rm = TRUE)),
                   by = .(Protein.Group, Run)
]

# ----------------------------
# Protein annotations (CETSA-MS SILAC: Temp + Replicate per run; L=DMSO, H=Drug)
# ----------------------------
ann_cols <- c("Protein.Group", "Run", "Sample", "Temp", "Replicate",
              "First.Protein.Description", "Genes")
ann_cols <- ann_cols[ann_cols %in% names(report)]

annotations <- unique(report[, ..ann_cols])

# ----------------------------
# Merge everything
# ----------------------------
merge_function <- function(x, y) merge(x, y, by = c("Protein.Group", "Run"), all = TRUE)

proteins <- Reduce(merge_function, list(
  LFQ_T, LFQ_H, LFQ_L,
  int_T, int_H, int_L,
  counts_H, counts_L,
  SILAC,
  annotations
))

# Remove protein-sample pairs without values
proteins <- proteins[ !is.na(Intensity) ]

# Tidy names
if ("First.Protein.Description" %in% names(proteins)) {
  setnames(proteins,
           old = c("Protein.Group", "First.Protein.Description"),
           new = c("Protein_group", "Protein_description"))
} else {
  setnames(proteins, old = "Protein.Group", new = "Protein_group")
}

# Order columns (by name; robust)
preferred_first <- intersect(
  c("Run", "Sample", "Temp", "Replicate", "Protein_group", "Protein_description", "Genes",
    "LFQ", "LFQ_H", "LFQ_L",
    "Intensity", "Intensity_H", "Intensity_L",
    "log2_HL_ratio", "N_ratios",
    "N_precursors_H", "N_precursors_L",
    "N_precursors_proteotypic_H", "N_precursors_proteotypic_L"),
  names(proteins)
)
setcolorder(proteins, c(preferred_first, setdiff(names(proteins), preferred_first)))

# Write out protein result table (optional)
fwrite(proteins, "protein_quant.csv.gz")
fwrite(report, "reportfile.tsv")




#### QC plotting section ####
# CETSA-MS SILAC design: Light (L) = DMSO control, Heavy (H) = Drug treated
# Each run = one temperature point; 11 temperatures x 2 replicates = 22 runs

# ----------------------------
# Set Sample plotting order: Temp (ascending) then Replicate
# ----------------------------
desired_temp_order      <- sort(unique(proteins$Temp))   # integer at this point
desired_replicate_order <- c("R1", "R2")

# Temperature colours: cold (blue) -> hot (red) gradient
# Keyed on character representation of integer Temp values
temp_levels        <- as.character(desired_temp_order)
temp_color_mapping <- setNames(colorRampPalette(c("#2166ac", "#d73027"))(length(temp_levels)), temp_levels)

# ----------------------------
# CV function
# ----------------------------
CV_function <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x) / mean(x) * 100
}

# ----------------------------
# (5a) Precursor CVs - computed BEFORE Temp is converted to factor
# ----------------------------
precursor_by_sample <- report[
  !is.na(Precursor.Quantity) & !is.na(label) & !is.na(Temp),
  .(qty = median(Precursor.Quantity, na.rm = TRUE)),
  by = .(Precursor.Id.nolabels, label, Temp, Sample)
]

precursor_CVs <- precursor_by_sample[
  ,
  .(PrecursorCV = CV_function(qty), N_samples = .N),
  by = .(Precursor.Id.nolabels, label, Temp)
][N_samples >= 2 & is.finite(PrecursorCV)]

precursor_CVs[, precursor_rank := frank(PrecursorCV), by = .(Temp, label)]
precursor_CVs[, label     := factor(label, levels = c("L", "H"))]
precursor_CVs[, Temp_char := as.character(Temp)]   # for colour mapping

message(sprintf("precursor_CVs rows: %d, label levels: %s",
                nrow(precursor_CVs), paste(levels(precursor_CVs$label), collapse = ", ")))

# ----------------------------
# (5b) Protein CVs - computed BEFORE Temp is converted to factor
# ----------------------------
protein_long <- melt(
  proteins,
  id.vars = c("Protein_group", "Temp", "Sample"),
  measure.vars = intersect(c("LFQ_L", "LFQ_H"), names(proteins)),
  variable.name = "label",
  value.name = "LFQ_intensity",
  na.rm = TRUE
)
protein_long[, label := gsub("LFQ_", "", label)]

protein_CVs <- protein_long[
  ,
  .(ProteinCV = CV_function(LFQ_intensity), N_samples = uniqueN(Sample)),
  by = .(Protein_group, label, Temp)
][N_samples >= 2 & is.finite(ProteinCV)]

protein_CVs[, protein_rank := frank(ProteinCV), by = .(Temp, label)]
protein_CVs[, label     := factor(label, levels = c("L", "H"))]
protein_CVs[, Temp_char := as.character(Temp)]   # for colour mapping

# ----------------------------
# Now convert Temp to factor for plotting axes
# ----------------------------
report[,   Temp      := factor(Temp,      levels = desired_temp_order)]
report[,   Replicate := factor(Replicate, levels = desired_replicate_order)]
proteins[, Temp      := factor(Temp,      levels = desired_temp_order)]
proteins[, Replicate := factor(Replicate, levels = desired_replicate_order)]

sample_levels <- proteins[order(Temp, Replicate), unique(Sample)]
report[,   Sample := factor(Sample, levels = sample_levels)]
proteins[, Sample := factor(Sample, levels = sample_levels)]

# Apply sample order
if ("Sample" %in% names(proteins) && "Sample" %in% names(report)) {
  proteins[, Sample := factor(Sample, levels = sample_levels)]
}

# ----------------------------
# Theme
# ----------------------------
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

# ----------------------------
# Colours
# L = DMSO (control) = purple/magenta; H = Drug (treated) = blue; consistent with SILAC convention
# ----------------------------
HL_cols <- c(heavy = "#37abc8ff", H = "#37abc8ff", light = "#ca00a3ff", L = "#ca00a3ff", total = "black")
scale_HL_labels <- c(L = "DMSO (light)", H = "Drug (heavy)")

# ----------------------------
# (1) Intensity distributions
# ----------------------------
plot_int_DT <- rbindlist(list(
  report[!is.na(Precursor.Quantity),
         .(Sample, `log2 quantity` = log2(pmax(Precursor.Quantity, 1)), Type = "Precursor quantity")],
  
  proteins[!is.na(LFQ),
           .(Sample, `log2 quantity` = log2(pmax(LFQ, 1)), Type = "Protein MaxLFQ")],
  
  proteins[!is.na(Intensity),
           .(Sample, `log2 quantity` = log2(pmax(Intensity, 1)), Type = "Protein Intensity")]
), use.names = TRUE, fill = TRUE)

plot_int_DT[, Type := factor(Type, levels = c("Precursor quantity", "Protein MaxLFQ", "Protein Intensity"))]

pInt <- ggplot(plot_int_DT, aes(x = Sample, y = `log2 quantity`)) +
  facet_wrap(~Type, scales = "free_x") +
  geom_boxplot(outliers = FALSE, linewidth = 0.25) +
  theme(axis.title.y = element_blank()) +
  coord_flip()

plot(pInt)

# ----------------------------
# (2) Cumulative intensity distributions (DMSO/Light vs Drug/Heavy)
# ----------------------------
cum_cols <- intersect(c("LFQ_H", "LFQ_L", "Intensity_H", "Intensity_L"), names(proteins))

plot_int <- proteins[, lapply(.SD, function(x) sum(x, na.rm = TRUE)), by = Sample, .SDcols = cum_cols]
plot_int <- melt(plot_int, id.vars = "Sample", value.name = "Cumulative quantity", variable.name = "label")

plot_int[, Type := fifelse(label %like% "LFQ", "Protein MaxLFQ", "Protein intensity")]
plot_int[, label := gsub(".+_", "", label)]
plot_int[, Type := factor(Type, levels = c("Protein MaxLFQ", "Protein intensity"))]

pCumInt <- ggplot(plot_int, aes(x = Sample, y = `Cumulative quantity`, fill = label)) +
  facet_wrap(~Type, scales = "free_x") +
  geom_bar(stat = "identity") +
  ylab("Cumulative protein intensity") +
  scale_fill_manual(values = HL_cols, guide = guide_legend(title = "Channel"), labels = scale_HL_labels) +
  theme(axis.title.y = element_blank()) +
  coord_flip()

plot(pCumInt)

# ----------------------------
# (3) Counts: precursors, peptides, protein groups (split by L=DMSO / H=Drug)
# ----------------------------
count_cols <- intersect(c("Precursor.Id", "Stripped.Sequence", "Protein.Group"), names(report))

if (!("label" %in% names(report))) report[, label := "total"]

counts <- report[ !is.na(label), lapply(.SD, uniqueN), by = .(Sample, label), .SDcols = count_cols]
counts <- melt(counts, id.vars = c("Sample", "label"), value.name = "IDs")

facet_labels <- as_labeller(c(
  Precursor.Id = "Precursors",
  Stripped.Sequence = "Peptides",
  Protein.Group = "Protein groups"
))

pCount <- ggplot(counts, aes(x = Sample, y = IDs/1000, fill = label,
                             label = format(IDs, big.mark = ",", scientific = FALSE))) +
  facet_wrap(~variable, scales = "free_x", labeller = facet_labels) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(size = 2, hjust = 1.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = HL_cols, guide = guide_legend(title = "Channel"), labels = scale_HL_labels) +
  ylab("Number of IDs [x1,000]") +
  theme(axis.title.y = element_blank()) +
  coord_flip()

plot(pCount)

# ----------------------------
# (4) Data completeness vs precursor cut-offs
# ----------------------------
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

pNAs <- ggplot(data_completeness, aes(x = N_samples, y = cumulative_protein_N/1000, colour = Precursors)) +
  geom_point(size = 0.8) +
  geom_line(linewidth = 0.25) +
  xlab("Number of samples") +
  ylab("Number of proteins [x1,000]") +
  scale_x_continuous(breaks = seq(1, max(data_completeness$N_samples), 1)) +
  scale_y_continuous(limits = c(0, max(data_completeness$cumulative_protein_N, na.rm = TRUE)/1000)) +
  scale_color_brewer(palette = "Dark2") +
  theme(panel.grid.major = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
        panel.grid.minor = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"))

plot(pNAs)

# ----------------------------
# (5a) Precursor CV plot (data computed above before factor conversion)
# ----------------------------
pCV1 <- ggplot(precursor_CVs, aes(x = precursor_rank/1000, y = PrecursorCV, colour = Temp_char)) +
  facet_wrap(~label, labeller = as_labeller(c(L = "DMSO (light)", H = "Drug (heavy)"))) +
  geom_line(linewidth = 0.25) +
  labs(x = "Number of precursors [x1,000]", y = "Precursor CV [%]") +
  coord_cartesian(ylim = c(0, 50)) +
  scale_color_manual(values = temp_color_mapping, guide = guide_legend(title = "Temp (°C)")) +
  theme(panel.grid.major = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
        panel.grid.minor = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
        legend.position = "right")

plot(pCV1)

# ----------------------------
# (5b) Protein CV plot (data computed above before factor conversion)
# ----------------------------
pCV2 <- ggplot(protein_CVs, aes(x = protein_rank/1000, y = ProteinCV, colour = Temp_char)) +
  facet_wrap(~label, labeller = as_labeller(c(L = "DMSO (light)", H = "Drug (heavy)"))) +
  geom_line(linewidth = 0.25) +
  labs(x = "Number of proteins [x1,000]", y = "Protein CV [%]") +
  coord_cartesian(ylim = c(0, 50)) +
  scale_color_manual(values = temp_color_mapping, guide = guide_legend(title = "Temp (°C)")) +
  theme(panel.grid.major = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
        panel.grid.minor = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
        legend.position = "right")

plot(pCV2)

# ----------------------------
# PCA (DMSO/Light and Drug/Heavy separately)
# colour = Temp (cold->hot gradient), shape = Replicate
# ----------------------------
PCA_input_L <- dcast(proteins, Protein_group ~ Sample, value.var = "LFQ_L")
PCA_input_L <- PCA_input_L[complete.cases(PCA_input_L)]
PCA_input_L <- t(data.frame(PCA_input_L, row.names = "Protein_group"))
DT_PCA_L <- prcomp(PCA_input_L, scale. = TRUE)
var_explained_L <- summary(DT_PCA_L)$importance
DT_PCA_L <- data.table(DT_PCA_L$x, keep.rownames = "Sample")
DT_PCA_L <- merge(DT_PCA_L, unique(proteins[, .(Sample, Temp, Replicate)]), by = "Sample", all.x = TRUE)
DT_PCA_L[, Temp_char := as.character(as.integer(as.character(Temp)))]

PCA_input_H <- dcast(proteins, Protein_group ~ Sample, value.var = "LFQ_H")
PCA_input_H <- PCA_input_H[complete.cases(PCA_input_H)]
PCA_input_H <- t(data.frame(PCA_input_H, row.names = "Protein_group"))
DT_PCA_H <- prcomp(PCA_input_H, scale. = TRUE)
var_explained_H <- summary(DT_PCA_H)$importance
DT_PCA_H <- data.table(DT_PCA_H$x, keep.rownames = "Sample")
DT_PCA_H <- merge(DT_PCA_H, unique(proteins[, .(Sample, Temp, Replicate)]), by = "Sample", all.x = TRUE)
DT_PCA_H[, Temp_char := as.character(as.integer(as.character(Temp)))]

pPCA_L <- ggplot(DT_PCA_L, aes(x = PC1, y = PC2, colour = Temp_char, shape = Replicate)) +
  ggtitle("DMSO control (Light)") +
  geom_point(size = 2.5) +
  xlab(paste0("PC1 [", round(var_explained_L["Proportion of Variance", "PC1"] * 100, 0), "%]")) +
  ylab(paste0("PC2 [", round(var_explained_L["Proportion of Variance", "PC2"] * 100, 0), "%]")) +
  scale_color_manual(values = temp_color_mapping, guide = guide_legend(title = "Temp (°C)")) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "right")
plot(pPCA_L)

pPCA_H <- ggplot(DT_PCA_H, aes(x = PC1, y = PC2, colour = Temp_char, shape = Replicate)) +
  ggtitle("Drug treated (Heavy)") +
  geom_point(size = 2.5) +
  xlab(paste0("PC1 [", round(var_explained_H["Proportion of Variance", "PC1"] * 100, 0), "%]")) +
  ylab(paste0("PC2 [", round(var_explained_H["Proportion of Variance", "PC2"] * 100, 0), "%]")) +
  scale_color_manual(values = temp_color_mapping, guide = guide_legend(title = "Temp (°C)")) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom")
plot(pPCA_H)

# ----------------------------
# (6) log2(H/L) ratio distributions across temperatures
# Key CETSA-MS QC: Drug/DMSO ratio should shift at stabilised temperatures
# ----------------------------
pHL <- ggplot(proteins[!is.na(log2_HL_ratio)],
              aes(x = Temp, y = log2_HL_ratio, fill = as.character(as.integer(as.character(Temp))))) +
  geom_boxplot(outliers = FALSE, linewidth = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey40") +
  scale_fill_manual(values = temp_color_mapping, guide = "none") +
  labs(x = "Temperature (°C)", y = "log2(Drug/DMSO) ratio",
       title = "H/L ratio distribution per temperature") +
  theme(panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"))

plot(pHL)

# ----------------------------
# Combine plots (patchwork)
# ----------------------------
combined_plot <- (
  pInt /
    pCumInt /
    pCount /
    patchwork::wrap_elements(pNAs + plot_spacer()) /
    patchwork::wrap_elements(pCV1 + pCV2) /
    patchwork::wrap_elements(pPCA_L + pPCA_H) /
    pHL
) +
  plot_layout(ncol = 1, heights = c(1, 1, 1.8, 1, 1, 1, 0.8)) +
  plot_annotation(theme = theme(plot.margin = margin(20, 20, 20, 20)))

plot(combined_plot)

ggsave(
  "DIANNprotein_quant_plot_LARGE.png",
  combined_plot,
  width = 35, height = 90, units = "cm",
  dpi = 300, bg = "white"
)