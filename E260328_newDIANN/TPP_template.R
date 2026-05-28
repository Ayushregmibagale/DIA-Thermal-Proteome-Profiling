# converting diann-reported LFQs to
# temperature-induced protein loss (FoldChange)
file_location = '/mnt/kustatscher/members/savvas/E260328_NEWDIANN/'
library(stringr)
library(data.table)
library(ggplot2)
library(pheatmap)
PG_report_DIANN = fread(
  glue::glue('{file_location}protein_quant.csv.gz')
)
PG_report_DIANN[
  str_detect(Genes, 'AXL|DDR1|MERTK|HIPK4|LOK|RET$|TIE1|FLT3|MET$'),
  Genes
] |>
  unique()

PG_report_DIANN_cond = PG_report_DIANN |>
  melt(
    id.vars = c(
      'Run',
      'Sample',
      'Protein_group',
      'Condition',
      'Replicate'
    ),
    measure = patterns('LFQ_'),
    variable.name = 'Channel',
    value.name = 'Abundance'
  )
nrow(PG_report_DIANN_cond)

PG_report_DIANN_cond[, `:=`(
  Drug = fifelse(Channel == 'LFQ_H', 'DMSO', 'Drug')
)]
PG_report_DIANN_cond = PG_report_DIANN_cond[order(Abundance, decreasing = F)]
PG_report_DIANN_cond[,
  level_of_detection := mean(head(Abundance, 200), na.rm = T),
  by = Run
]
PG_report_DIANN_cond[, mean_abundance := mean(Abundance), by = Protein_group]

PG_report_DIANN_cond[, Temp := Condition]
PG_report_DIANN_cond |>
  ggplot(aes(x = as.factor(Temp), colour = Drug, y = log2(Abundance))) +
  geom_boxplot() +
  theme_bw() +
  geom_point(aes(y = log2(level_of_detection)), colour = 'red') +
  labs(x = 'Temp', y = 'Log2(MaxLFQ)')

PG_report_DIANN_cond |>
  ggplot(aes(
    x = as.factor(Temp),
    colour = Drug,
    y = Abundance / mean_abundance
  )) +
  geom_boxplot(outliers = F) +
  theme_bw() +
  # geom_point(aes(y = log2(level_of_detection)), colour = 'red') +
  labs(x = 'Temp', y = 'Log2(MaxLFQ)')

PG_report_DIANN_cond |>
  ggplot(aes(
    x = as.factor(Temp),
    colour = as.character(Replicate),
    y = log2(Abundance)
  )) +
  geom_boxplot() +
  theme_bw() +
  geom_point(aes(y = log2(level_of_detection)), colour = 'red') +
  labs(x = 'Temp', y = 'Log2(MaxLFQ)') +
  facet_wrap('Drug')

PG_report_DIANN_cond[, log2_Abundance := log2(Abundance)]
PG_report_DIANN_cond[, Sample := paste(Drug, Temp, Replicate, sep = '_')]
report_wide = PG_report_DIANN_cond |>
  dcast(Protein_group ~ Sample, value.var = 'log2_Abundance') |>
  tibble::column_to_rownames('Protein_group')
NA_heatmap = report_wide
NA_heatmap[is.na(NA_heatmap)] <- 0
pheatmap(NA_heatmap, show_rownames = F)

missingness_map = is.na(report_wide) |> apply(2, as.numeric)
rownames(missingness_map) <- rownames(report_wide)
N_prot = missingness_map |> matrixStats::rowSums2()
N_prot |>
  hist(
    main = 'any proteins with less than 21 present values are removed',
    xlab = 'Number of missing values'
  )

PG_report_DIANN_cond[Protein_group == 'P08581'] |>
  ggplot(aes(x = Temp, y = log2_Abundance, colour = Drug)) +
  geom_point()
# change to allow more or less missing values
N_missing_values_allowed = 44
PG_report_DIANN_cond
report_wide |>
  boxplot(
    main = 'clearly median normalisation would remove the biological trend'
  )
# maybe sum or no normalisation if prenormalised for loading
# report_wide = limma::normalizeMedianValues(report_wide)
# report_wide |> boxplot()

annot_col = PG_report_DIANN_cond[, .(Temp, Replicate, Drug, Sample)] |>
  unique()
annot_col = annot_col |> tibble::column_to_rownames('Sample')

missingness_map = missingness_map[N_prot <= N_missing_values_allowed, ]
pheatmap_clust = pheatmap(
  missingness_map,
  show_rownames = F,
  show_colnames = F,
  annotation_col = annot_col,
  main = '1 is missing'
)

protein_clusters <- data.table(
  Protein_group = missingness_map |> rownames(missingness_map),
  cluster = as.factor(cutree(pheatmap_clust$tree_row, k = 5))
)
protein_clusters |>
  ggplot(aes(x = cluster, fill = cluster)) +
  geom_bar() +
  scale_fill_manual(
    values = c(
      `1` = '#ff84c3',
      `2` = '#ff9289',
      `3` = '#00d1ff',
      `4` = '#e0b400',
      `5` = '#fa8eff'
    )
  ) +
  ggtitle('How many proteins per missingness cluster')
pheatmap(
  missingness_map,
  show_rownames = F,
  show_colnames = F,
  annotation_col = annot_col,
  annotation_row = protein_clusters |>
    tibble::column_to_rownames('Protein_group'),
  main = '1 is missing  \n all clusters seem biologically missing'
)
# correction_values = (1-(missingness_map |> matrixStats::colMeans2()))
# correction_values_df = data.table(names = names(correction_values),
#                                   perc_missingness = correction_values)
# correction_values_df[,Temp:= str_extract(names,'TPP_..')]
# correction_values_df |> ggplot(aes(y= names,
#                                    fill = Temp,
#                                    x = perc_missingness))+
#   geom_col()

# BCA_37 = openxlsx::read.xlsx(
#   here::here('in', 'datasets', 'BCAst.curve_E250320_Cristina.xlsx'),
#   startRow = 29
# ) |>
#   as.data.table()
# BCA_37 = BCA_37[, .(`37-1`, Rep_37, `mg/ml_37`)] |> na.omit()
# setnames(BCA_37, c('Sample', 'Rep', 'Volume'))
# BCA_37[, Temp := '37']
# BCA_52 = openxlsx::read.xlsx(
#   here::here('in', 'datasets', 'BCAst.curve_E250320_Cristina.xlsx'),
#   startRow = 29
# ) |>
#   as.data.table()
# BCA_52 = BCA_52[, .(`52-1`, Rep_52, `mg/ml_52`)] |> na.omit()
# setnames(BCA_52, c('Sample', 'Rep', 'Volume'))
# BCA_52[, Temp := '52']
# BCA = rbind(BCA_37, BCA_52)
# concentrations = data.table(
#   Sample = paste0('S', 1:6),
#   Drug = c('0', '1', '100', '500', '1000', '5000')
# )
# BCA = merge(BCA, concentrations, by = 'Sample')
# BCA[, mgml := as.numeric(Volume)]
# BCA_avg = BCA[, .(Avg_mgml = median(mgml, na.rm = T)), by = .(Temp, Drug)]
# max_mgml = BCA_avg$Avg_mgml |> max()

# BCA[, correction := mgml / max_mgml, by = Sample]

# BCA |> ggplot(aes(x = Temp, colour = Drug, y = correction)) + geom_point()
# BCA = BCA[, .(avg_correction = median(correction)), by = .(Temp, Drug)]

PGs_MNAR <- report_wide[rownames(missingness_map), ]
PGs_MNAR |> dim()
PGs_MNAR_norm = PGs_MNAR |>
  as.data.frame() |>
  tibble::rownames_to_column('Protein_group') |>
  as.data.table() |>
  melt(
    id.vars = 'Protein_group',
    value.name = 'Norm_log2_LFQ',
    variable.name = 'Run'
  ) |>
  as.data.table()
# PGs_MNAR_norm = merge(PGs_MNAR_norm, exp_design, by.x = 'Run', by.y = 'label')
PGs_MNAR_norm[,
  level_of_detection := mean(
    head(sort(Norm_log2_LFQ, decreasing = F), 100),
    na.rm = T
  ),
  by = Run
]


PGs_MNAR_norm[, condition := str_remove(Run, '_R.$')]
PGs_MNAR_norm |>
  ggplot(aes(x = condition, y = Norm_log2_LFQ)) +
  geom_boxplot() +
  geom_point(
    data = unique(PGs_MNAR_norm[, .(condition, Run, level_of_detection)]),
    aes(y = level_of_detection),
    colour = 'red'
  )

PGs_MNAR_norm[,
  N_missing := sum(is.na(Norm_log2_LFQ)),
  by = .(Protein_group, condition)
]

# calcualtes the next/previous Temp series
# shift -1 is previous, shift +1 is next
relative_temp = function(current_temp, Temp = Temp_factor, shift = -1) {
  current_temp = unique(current_temp)
  current_position = which(Temp == current_temp)
  next_positition = current_position + shift
  if (between(next_positition, 1, length(Temp_factor))) {
    return(Temp[next_positition])
  } else {
    return(NA_character_)
  }
}
PG_report_DIANN[, Temp := Condition]
Temp_factor = PG_report_DIANN$Temp |> unique() |> sort() |> as.factor()
PGs_MNAR_norm[,
  Temp := str_extract(condition, '[:digit:]*$') |>
    factor(levels = Temp_factor)
] #
PGs_MNAR_norm[, Drug := str_remove(condition, '_[:print:]*$')]
PGs_MNAR_norm[, Replicate := str_remove(Run, '[:print:]*_')]
# calcualting next Temp point
PGs_MNAR_norm[,
  next_temp := relative_temp(Temp, Temp_factor, shift = 1) |> as.character(),
  by = .(Temp)
]
PGs_MNAR_norm[,
  N_Missing := sum(is.na(Norm_log2_LFQ)),
  by = .(Protein_group, condition)
]
# counting values in next Temppoint
next_temp = PGs_MNAR_norm[,
  .(next_temp_samples = mean(N_Missing)),
  by = .(Temp, Protein_group, Drug)
]
setnames(next_temp, 'Temp', 'next_temp')

PGs_MNAR_norm = merge(
  PGs_MNAR_norm,
  next_temp,
  by = c('next_temp', 'Protein_group', 'Drug'),
  all.x = T
)
# Calculating for the same timepoint
# how many missing in each conditions to see if the missingness
# is biologically informative, and thus would need imputation
Diff_missing = unique(PGs_MNAR_norm[, .(
  Protein_group,
  Drug,
  Temp,
  N_Missing
)]) |>
  dcast(Protein_group + Temp ~ Drug, value.var = 'N_Missing')
Diff_missing[, Diff_missing := Drug - DMSO]

PGs_MNAR_norm = merge(
  PGs_MNAR_norm,
  Diff_missing[, .(Protein_group, Temp, Diff_missing)],
  by = c('Protein_group', 'Temp')
)

PGs_MNAR_norm[
  next_temp_samples > 1 & N_Missing > 1 & abs(Diff_missing) > 1,
  impute := T
]
PGs_MNAR_norm[is.na(impute), impute := F]
PGs_MNAR_norm[,
  imputted_Norm_log2_LFQ := fifelse(
    is.na(Norm_log2_LFQ) & impute == T,
    rnorm(.N, mean = level_of_detection, 0.1),
    Norm_log2_LFQ
  )
]
example_imputation = PGs_MNAR_norm[impute == T][10, Protein_group] |> unlist()

PGs_MNAR_norm[Protein_group == example_imputation] |>
  ggplot(aes(
    x = Temp,
    colour = Drug,
    y = imputted_Norm_log2_LFQ,
    alpha = impute
  )) +
  geom_point() +
  scale_alpha_manual(values = c(1, 0.5))


PGs_MNAR_norm |>
  ggplot(aes(
    x = as.factor(Drug),
    colour = as.factor(Replicate),
    y = imputted_Norm_log2_LFQ
  )) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap('Temp')

PGs_MNAR_norm |>
  ggplot(aes(
    x = as.factor(Drug),
    colour = as.factor(Replicate),
    y = Norm_log2_LFQ
  )) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap('Temp')


# PGs_MNAR_norm = merge(PGs_MNAR_norm, BCA, by = c('Temp', 'Drug'))
# # PGs_MNAR_norm[,temp_cor:= mean(correction),by = Temp]
# PGs_MNAR_norm[, corrected_norm := (2^(imputted_Norm_log2_LFQ)) * avg_correction]
# PGs_MNAR_norm |>
#   ggplot(aes(
#     x = as.factor(Drug),
#     colour = as.factor(Replicate),
#     y = log2(corrected_norm)
#   )) +
#   geom_boxplot() +
#   theme_bw() +
#   ggtitle('All proteins have now values') +
#   facet_wrap('Temp')

PGs_MNAR_norm = merge(
  PGs_MNAR_norm,
  unique(PG_report_DIANN[, .(Protein_group, Genes)]),
  by.x = 'Protein_group',
  by.y = 'Protein_group',
  all.x = T
)
conc = openxlsx::read.xlsx('Equal Conc.xlsx') |>
  as.data.table()
conc = conc[, .(Sample, `Adjusted.concentration.(ug/ul)`)] |> na.omit()
conc[, `:=`(
  Temp = str_extract(Sample, '_[:digit:]*') |> str_remove('_'),
  Drug = str_remove(Sample, '_[:print:]*$')
)]
conc[, mean_conc := mean(`Adjusted.concentration.(ug/ul)`)]

conc |>
  ggplot(aes(
    x = Temp,
    colour = Drug,
    y = `Adjusted.concentration.(ug/ul)` / mean_conc
  )) +
  geom_point()
conc = conc[,
  .(mean_adj_con = mean(`Adjusted.concentration.(ug/ul)` / mean_conc)),
  by = Temp
]
PGs_MNAR_norm = merge(PGs_MNAR_norm, conc, by = 'Temp')
PGs_MNAR_norm[, Abundance := 2^imputted_Norm_log2_LFQ]
PGs_MNAR_norm[, Abundance := Abundance * mean_adj_con]
PGs_MNAR_norm[,
  N_missing := sum(is.na(Abundance)),
  by = .(Protein_group)
]
PGs_MNAR_norm[,
  N_missing_cond := sum(is.na(Abundance)),
  by = .(Protein_group, Drug)
]
PGs_MNAR_norm[, N_missing_cond := max(N_missing_cond), by = Protein_group]
PGs_MNAR_norm_test_GPmelt = PGs_MNAR_norm[N_missing < 25 & N_missing_cond < 15]
PGs_MNAR_norm_test_GPmelt[, Abundance := log2(Abundance)]
PGs_MNAR_norm_test_GPmelt[,
  mean_expression := mean(Abundance, na.rm = T),
  by = .(Protein_group, Drug)
]
PGs_MNAR_norm_test_GPmelt[, rep := str_remove(Replicate, 'R')]

PGs_MNAR_norm_test_GPmelt[
  Protein_group == 'P04083'
] |>
  ggplot(aes(
    x = Temp |> as.character() |> as.numeric(),
    y = Abundance,
    colour = Drug
  )) +
  geom_point(alpha = 0.5) +
  theme_bw()

mean_ratios_fd = PGs_MNAR_norm_test_GPmelt[, .(
  protein_1 = Protein_group,
  condition = Drug,
  Level_1 = paste(Genes, Protein_group, sep = '_'),
  Level_2 = Drug,
  Level_3 = paste(Drug, rep, sep = "."),
  rep = as.numeric(rep),
  x = Temp |> as.character() |> as.numeric(),
  y = Abundance / mean_expression
)]
starting_abundance = mean_ratios_fd[x == min(mean_ratios_fd$x, na.rm = T)]
starting_abundance[
  Level_2 == 'DMSO',
  baseline := mean(y, na.rm = T),
  by = Level_1
]
starting_abundance[, baseline := mean(baseline, na.rm = T), by = Level_1]
starting_abundance = starting_abundance[Level_2 != 'DMSO'][,
  .(mean_baseline = mean(baseline, na.rm = T), mean_y = mean(y, .na.rm = T)),
  by = Level_1
][, condition_diff := mean_y - mean_baseline, by = Level_1]

mean_ratios_fd |>
  ggplot(aes(x = as.factor(x), y = y, colour = Level_2)) +
  geom_boxplot(outliers = F) +
  facet_wrap('rep')

write.csv(
  mean_ratios_fd,
  '/home/v1skourt/gpmelt/NextFlow/dummy_data/E260328_all/dataForGPMelt.csv'
)
write.csv(
  data.table(Level_1 = unique(mean_ratios_fd$Level_1), NSamples = 5),
  '/home/v1skourt/gpmelt/NextFlow/dummy_data/E260328_all/NumberSamples_perID.csv'
)
write.csv(
  data.table(
    Level_1 = unique(mean_ratios_fd$Level_1) |> str_subset('Q9Y6Y8|A0A0B4J2F0')
  ),
  '/home/v1skourt/gpmelt/NextFlow/dummy_data/E260328_all/subset_ID.csv'
)
# cd gpmelt/NextFlow/
# conda activate gpmelt_env_with_dask
#  nextflow run GPMelt_workflow.nf -c E260328_all_params.config  -with-report reports/GPMelt_workflow_report.html
loss_values <- fread(
  "/home/v1skourt/gpmelt/NextFlow/output/E260328_all/full_hgpm_fit_loss_values_df.csv"
)
loss_values[, iteration := 1:.N, by = fitted_levels_1]

# all proteins have reached convergence
ggplot(loss_values[
  fitted_levels_1 %in% unique(loss_values$fitted_levels_1)[1:12]
]) +
  geom_point(aes(x = iteration, y = LossValue), size = 0.5) +
  facet_wrap(. ~ fitted_levels_1, scale = "free") +
  ggtitle("Convergence analysis") +
  theme_bw(base_size = 10)

# good to see no missing values
loss_values[is.na(LossValue)][order(fitted_levels_1)]

# checking convergence in the last iterations
N_last_iterations = 100

FinalIterations <- loss_values[,
  tail(.SD, N_last_iterations),
  by = fitted_levels_1
][,
  .(
    MedianLoss = median(LossValue, na.rm = TRUE),
    MADLoss = mad(LossValue, na.rm = TRUE)
  ),
  by = fitted_levels_1
]

head(FinalIterations)

ggplot(FinalIterations, aes(x = 0, y = log(MADLoss))) +
  geom_violin(width = 0.3) +
  geom_boxplot(width = 0.1, color = "grey", alpha = 0.2) +
  geom_point() +
  #geom_jitter(width=0.1, alpha=0.9) + # geom_jitter can be used when more data points are present
  ggtitle(paste0(
    "Variations around the median loss value for the last ",
    N_last_iterations,
    " iterations"
  )) +
  xlab("") +
  ylab("Median Absolute Deviation [log]") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  )

p_values_subset <- fread(
  "/home/v1skourt/gpmelt/NextFlow/output/E260328_all/p_values_dataset_wise.csv"
)
library(dplyr)
p_values_subset$pValue |> hist()
order_levels_1 = p_values_subset[pValue < 0.005, Level_1]
custom_labeller <- function(value) {
  paste0(Hmisc::capitalize(value), " model")
}

# create a list to store the plots of each level
all_plots <- vector(mode = "list", length = length(order_levels_1))
names(all_plots) <- order_levels_1

for (i in seq_along(order_levels_1)) {
  all_plots[[i]] <- list("Level_1" = NA, "Levels_1to2" = NA, "Levels_1to3" = NA)
}

prediction_full_hgpm_Level_1 <- read.csv(
  "/home/v1skourt/gpmelt/NextFlow/output/E260328_all/prediction_full_hgpm_Level_1.csv"
)

prediction_joint_hgpm_Level_1 <- read.csv(
  "/home/v1skourt/gpmelt/NextFlow/output/E260328_all/prediction_joint_hgpm_Level_1.csv"
)

# if some scaling has been computed, we load the scaled data
if (
  file.exists(
    "/home/v1skourt/gpmelt/NextFlow/output/E260328_all/data_with_scalingFactors.csv"
  )
) {
  real_dataset <- read.csv(
    "/home/v1skourt/gpmelt/NextFlow/output/E260328_all/data_with_scalingFactors.csv"
  )
} else {
  # if no scaling has been computed, then we can use the input data
  real_dataset <- read.csv(
    "/home/v1skourt/gpmelt/NextFlow/dummy_data/E260328_all/dataForGPMelt.csv"
  )
}

prediction_hgpm_Level_1 <- rbind(
  prediction_joint_hgpm_Level_1 %>%
    dplyr::select(
      -comparison_id,
      -treatment_condition,
      -control_condition,
      -Level_1_ForNullHyp
    ),
  prediction_full_hgpm_Level_1
)

for (lev1 in order_levels_1) {
  all_plots[[lev1]][["Level_1"]] <- ggplot(
    prediction_hgpm_Level_1 %>%
      dplyr::filter(Level_1 == lev1) %>%
      dplyr::mutate(
        Level_2 = factor(
          "protein-level",
          levels = c("Drug", "DMSO", "joint_condition", "protein-level")
        )
      )
  ) +
    geom_line(
      mapping = aes(x = x, y = y, color = Level_2),
      show.legend = TRUE
    ) +
    geom_ribbon(
      mapping = aes(
        x = x,
        ymin = conf_lower,
        ymax = conf_upper,
        fill = Level_2
      ),
      show.legend = TRUE,
      alpha = 0.7
    ) +
    geom_point(
      data = real_dataset %>% dplyr::filter(Level_1 == lev1),
      mapping = aes(x = x, y = y, shape = Level_3)
    ) +
    facet_grid(. ~ model, labeller = labeller(model = custom_labeller)) +
    scale_shape_manual(
      name = "Level 3",
      values = c("DMSO.1" = 15, "DMSO.2" = 16, "Drug.1" = 17, "Drug.2" = 18),
      labels = c(
        "DMSO.1" = "DMSO rep 1",
        "DMSO.2" = "DMSO rep 2",
        "Drug.1" = "Treatment rep 1",
        "Drug.2" = "Treatment rep 2"
      ),
      breaks = c("DMSO.1", "DMSO.2", "Drug.1", "Drug.2")
    ) +
    scale_color_manual(
      name = "Condition",
      values = c(
        "Drug" = "darkorange2",
        "DMSO" = "aquamarine3",
        "joint_condition" = "black",
        "protein-level" = "#72D8FF"
      ),
      labels = c(
        "Drug" = "Treatment (Foretinib)",
        "DMSO" = "DMSO",
        "joint_condition" = "Joint",
        "protein-level" = "protein-level"
      ),
      breaks = c("Drug", "DMSO", "joint_condition", "protein-level"),
      drop = FALSE
    ) +
    scale_fill_manual(
      name = "Condition",
      values = c(
        "Drug" = "darkorange2",
        "DMSO" = "aquamarine3",
        "joint_condition" = "black",
        "protein-level" = "#72D8FF"
      ),
      labels = c(
        "Drug" = "Treatment (Foretinib)",
        "DMSO" = "DMSO",
        "joint_condition" = "Joint",
        "protein-level" = "protein-level"
      ),
      breaks = c("Drug", "DMSO", "joint_condition", "protein-level"),
      drop = FALSE
    ) +
    guides(
      shape = guide_legend(
        nrow = 1,
        byrow = FALSE,
        override.aes = list(fill = NA, linetype = 0)
      )
    ) +
    xlab("Temperature [°C]") +
    xlab("Temperature [°C]") +
    ylab("Fold Change") +
    #ggtitle("", subtitle= "Prediction for Level 1")+
    theme_bw(base_size = 12) +
    theme(
      legend.box = "vertical",
      legend.margin = margin(),
      legend.text = element_text(size = 10)
    )
  #aspect.ratio = 0.8)
}


prediction_full_hgpm_Levels_1to2 <- read.csv(
  "/home/v1skourt/gpmelt/NextFlow/output/E260328_all/prediction_full_hgpm_Levels_1to2.csv"
)

prediction_joint_hgpm_Levels_1to2 <- read.csv(
  "/home/v1skourt/gpmelt/NextFlow/output/E260328_all/prediction_joint_hgpm_Levels_1to2.csv"
)

prediction_hgpm_Levels_1to2 <- rbind(
  prediction_joint_hgpm_Levels_1to2 %>%
    dplyr::select(
      -comparison_id,
      -treatment_condition,
      -control_condition,
      -Level_1_ForNullHyp
    ),
  prediction_full_hgpm_Levels_1to2
)

for (lev1 in order_levels_1) {
  all_plots[[lev1]][["Levels_1to2"]] <-
    ggplot(
      prediction_hgpm_Levels_1to2 %>%
        filter(Level_1 == lev1) %>%
        mutate(
          Level_2 = factor(
            Level_2,
            levels = c("Drug", "DMSO", "joint_condition", "protein-level")
          )
        )
    ) +
    geom_line(mapping = aes(x = x, y = y, color = Level_2)) +
    geom_ribbon(
      mapping = aes(
        x = x,
        ymin = conf_lower,
        ymax = conf_upper,
        fill = Level_2
      ),
      alpha = 0.7
    ) +
    geom_point(
      data = real_dataset %>% filter(Level_1 == lev1),
      mapping = aes(x = x, y = y, shape = Level_3)
    ) +
    facet_grid(. ~ model, labeller = labeller(model = custom_labeller)) +
    scale_shape_manual(
      name = "Level 3",
      values = c("DMSO.1" = 15, "DMSO.2" = 16, "Drug.1" = 17, "Drug.2" = 18),
      labels = c(
        "DMSO.1" = "DMSO rep 1",
        "DMSO.2" = "DMSO rep 2",
        "Drug.1" = "Treatment rep 1",
        "Drug.2" = "Treatment rep 2"
      ),
      breaks = c("DMSO.1", "DMSO.2", "Drug.1", "Drug.2")
    ) +
    scale_color_manual(
      name = "Condition",
      values = c(
        "Drug" = "darkorange2",
        "DMSO" = "aquamarine3",
        "joint_condition" = "black",
        "protein-level" = "#72D8FF"
      ),
      labels = c(
        "Drug" = "Treatment (Foretinib)",
        "DMSO" = "DMSO",
        "joint_condition" = "Joint",
        "protein-level" = "protein-level"
      ),
      breaks = c("Drug", "DMSO", "joint_condition", "protein-level"),
      drop = FALSE
    ) +
    scale_fill_manual(
      name = "Condition",
      values = c(
        "Drug" = "darkorange2",
        "DMSO" = "aquamarine3",
        "joint_condition" = "black",
        "protein-level" = "#72D8FF"
      ),
      labels = c(
        "Drug" = "Treatment (Foretinib)",
        "DMSO" = "DMSO",
        "joint_condition" = "Joint",
        "protein-level" = "protein-level"
      ),
      breaks = c("Drug", "DMSO", "joint_condition", "protein-level"),
      drop = FALSE
    ) +
    guides(shape = guide_legend(nrow = 1, byrow = FALSE)) +
    xlab("Temperature [°C]") +
    ylab("Fold Change") +
    #ggtitle("", subtitle= "Prediction for Level 2")+
    theme_bw(base_size = 12) +
    theme(
      legend.box = "vertical",
      legend.margin = margin(),
      legend.text = element_text(size = 10)
    )
  #aspect.ratio = 0.8)
}

prediction_full_hgpm_Levels_1to3 <- read.csv(
  "/home/v1skourt/gpmelt/NextFlow/output/E260328_all/prediction_full_hgpm_Levels_1to3.csv"
)

prediction_joint_hgpm_Levels_1to3 <- read.csv(
  "/home/v1skourt/gpmelt/NextFlow/output/E260328_all/prediction_joint_hgpm_Levels_1to3.csv"
)

prediction_hgpm_Levels_1to3 <- rbind(
  prediction_joint_hgpm_Levels_1to3 %>%
    dplyr::select(
      -comparison_id,
      -treatment_condition,
      -control_condition,
      -Level_1_ForNullHyp
    ),
  prediction_full_hgpm_Levels_1to3
)

Level_3_labeller <- c(
  "DMSO.1" = "DMSO rep 1",
  "DMSO.2" = "DMSO rep 2",
  "Drug.1" = "Treatment rep 1",
  "Drug.2" = "Treatment rep 2"
)


for (lev1 in order_levels_1) {
  all_plots[[lev1]][["Levels_1to3"]] <-
    ggplot(
      prediction_hgpm_Levels_1to3 %>%
        filter(Level_1 == lev1) %>%
        mutate(
          Level_2 = factor(
            Level_2,
            levels = c("Drug", "DMSO", "joint_condition", "protein-level")
          )
        )
    ) +
    geom_line(
      mapping = aes(x = x, y = y, color = Level_2),
      show.legend = TRUE
    ) +
    geom_ribbon(
      mapping = aes(
        x = x,
        ymin = conf_lower,
        ymax = conf_upper,
        fill = Level_2
      ),
      alpha = 0.7,
      show.legend = TRUE
    ) +
    geom_point(
      data = real_dataset %>% filter(Level_1 == lev1),
      mapping = aes(x = x, y = y, shape = Level_3)
    ) +
    facet_grid(
      . ~ model + Level_3,
      labeller = labeller(model = custom_labeller, Level_3 = Level_3_labeller)
    ) +
    scale_shape_manual(
      name = "Level 3",
      values = c("DMSO.1" = 15, "DMSO.2" = 16, "Drug.1" = 17, "Drug.2" = 18),
      labels = c(
        "DMSO.1" = "DMSO rep 1",
        "DMSO.2" = "DMSO rep 2",
        "Drug.1" = "Treatment rep 1",
        "Drug.2" = "Treatment rep 2"
      ),
      breaks = c("DMSO.1", "DMSO.2", "Drug.1", "Drug.2")
    ) +
    scale_color_manual(
      name = "Condition",
      values = c(
        "Drug" = "darkorange2",
        "DMSO" = "aquamarine3",
        "joint_condition" = "black",
        "protein-level" = "#72D8FF"
      ),
      labels = c(
        "Drug" = "Treatment (Foretinib)",
        "DMSO" = "DMSO",
        "joint_condition" = "Joint",
        "protein-level" = "protein-level"
      ),
      breaks = c("Drug", "DMSO", "joint_condition", "protein-level"),
      drop = FALSE
    ) +
    scale_fill_manual(
      name = "Condition",
      values = c(
        "Drug" = "darkorange2",
        "DMSO" = "aquamarine3",
        "joint_condition" = "black",
        "protein-level" = "#72D8FF"
      ),
      labels = c(
        "Drug" = "Treatment (Foretinib)",
        "DMSO" = "DMSO",
        "joint_condition" = "Joint",
        "protein-level" = "protein-level"
      ),
      breaks = c("Drug", "DMSO", "joint_condition", "protein-level"),
      drop = FALSE
    ) +
    guides(
      shape = guide_legend(
        nrow = 1,
        byrow = FALSE,
        override.aes = list(fill = NA, linetype = 0)
      )
    ) +
    xlab("Temperature [°C]") +
    ylab("Fold Change") +
    #ggtitle("", subtitle= "Prediction for Level 3")+
    theme_bw(base_size = 12) +
    theme(
      legend.box = "vertical",
      legend.margin = margin(),
      legend.text = element_text(size = 10)
    )
  #aspect.ratio = 1)
}
library(ggpubr)
for (lev1 in order_levels_1[c(1:5)]) {
  full_hierarchy <- ggarrange(
    all_plots[[lev1]][["Level_1"]] + ggtitle("Predictions for Level 1"),
    all_plots[[lev1]][["Levels_1to2"]] + ggtitle("Predictions for Level 2"),
    all_plots[[lev1]][["Levels_1to3"]] + ggtitle("Predictions for Level 3"),
    ncol = 1,
    nrow = 3,
    common.legend = TRUE,
    legend = "bottom"
  )

  print(annotate_figure(
    full_hierarchy,
    top = text_grob(lev1, face = "bold", size = 14)
  ))
}


for (lev1 in str_subset(
  order_levels_1,
  'AXL|DDR1|MERTK|HIPK4|LOK|RET|TIE1|FLT3|MET'
)) {
  full_hierarchy <- ggarrange(
    all_plots[[lev1]][["Level_1"]] + ggtitle("Predictions for Level 1"),
    all_plots[[lev1]][["Levels_1to2"]] + ggtitle("Predictions for Level 2"),
    all_plots[[lev1]][["Levels_1to3"]] + ggtitle("Predictions for Level 3"),
    ncol = 1,
    nrow = 3,
    common.legend = TRUE,
    legend = "bottom"
  )
}
print(annotate_figure(
  full_hierarchy,
  top = text_grob(lev1, face = "bold", size = 14)
))
CoeffUncertainty = 2 # you can increase this if you want to include only extremely reliable predictions

RegionOfAcceptableUncertainty <- prediction_full_hgpm_Levels_1to2 %>%
  mutate(Uncertainty = abs(conf_upper - conf_lower)) %>%
  group_by(Level_1, Level_2) %>%
  mutate(MedianUncertainty = median(Uncertainty)) %>%
  mutate(
    KeptTemperature = (Uncertainty <= CoeffUncertainty * MedianUncertainty)
  ) %>%
  mutate(Diff = c(0, diff(KeptTemperature)))

index2 = c(
  1,
  which(RegionOfAcceptableUncertainty$Diff != 0),
  nrow(RegionOfAcceptableUncertainty) + 1
)
RegionOfAcceptableUncertainty$grp = rep(1:length(diff(index2)), diff(index2))

RegionOfAcceptableUncertainty_N_temp <- RegionOfAcceptableUncertainty %>%
  group_by(Level_1, Level_2) %>%
  distinct() %>%
  summarise(N_temp_kept = sum(KeptTemperature), N_temp_total = n_distinct(x))

RegionOfAcceptableUncertainty_N_temp %>%
  ggplot() +
  geom_bar(
    aes(x = 100 * (N_temp_kept / N_temp_total))
  ) +
  xlab("Number of temperatures kept per ID [%]") +
  theme_bw()
list()
for (lev1 in c("AXL_P30530", "METTL16_Q86W50")) {
  gg <- ggplot(
    prediction_full_hgpm_Levels_1to2 %>%
      filter(Level_1 == lev1)
  ) +
    # Test points definitions
    geom_vline(
      data = RegionOfAcceptableUncertainty %>%
        filter(Level_1 == lev1),
      mapping = aes(xintercept = x, linetype = KeptTemperature),
      color = "black",
      alpha = 0.7
    ) +
    # condition-level fits
    geom_line(mapping = aes(x = x, y = y, color = Level_2)) +
    geom_ribbon(
      mapping = aes(
        x = x,
        ymin = conf_lower,
        ymax = conf_upper,
        fill = Level_2
      ),
      alpha = 0.7
    ) +
    # real data
    geom_point(
      data = real_dataset %>% filter(Level_1 == lev1),
      mapping = aes(x = x, y = y, shape = Level_3)
    ) +
    facet_grid(. ~ model, labeller = labeller(model = custom_labeller)) +
    scale_shape_manual(
      name = "Level 3",
      values = c("DMSO.1" = 15, "DMSO.2" = 16, "Drug.1" = 17, "Drug.2" = 18),
      labels = c(
        "DMSO.1" = "DMSO rep 1",
        "DMSO.2" = "DMSO rep 2",
        "Drug.1" = "Treatment rep 1",
        "Drug.2" = "Treatment rep 2"
      ),
      breaks = c("DMSO.1", "DMSO.2", "Drug.1", "Drug.2")
    ) +
    scale_color_manual(
      name = "Condition",
      values = c(
        "Drug" = "darkorange2",
        "DMSO" = "aquamarine3",
        "joint_condition" = "black",
        "protein-level" = "#72D8FF"
      ),
      labels = c(
        "Drug" = "Treatment (Foretinib)",
        "DMSO" = "DMSO",
        "joint_condition" = "Joint",
        "protein-level" = "protein-level"
      ),
      breaks = c("Drug", "DMSO", "joint_condition", "protein-level"),
      drop = FALSE
    ) +
    scale_fill_manual(
      name = "Condition",
      values = c(
        "Drug" = "darkorange2",
        "DMSO" = "aquamarine3",
        "joint_condition" = "black",
        "protein-level" = "#72D8FF"
      ),
      labels = c(
        "Drug" = "Treatment (Foretinib)",
        "DMSO" = "DMSO",
        "joint_condition" = "Joint",
        "protein-level" = "protein-level"
      ),
      breaks = c("Drug", "DMSO", "joint_condition", "protein-level"),
      drop = FALSE
    ) +
    guides(shape = guide_legend(nrow = 1, byrow = FALSE)) +
    scale_linetype_manual(
      name = "Test points",
      values = c("TRUE" = "dotted", "FALSE" = "solid"),
      labels = c("TRUE" = "Included", "FALSE" = "Excluded")
    ) +
    xlab("Temperature [°C]") +
    ylab("Fold Change") +
    ggtitle(lev1) +
    theme_bw(base_size = 12) +
    theme(
      legend.box = "vertical",
      legend.margin = margin(),
      legend.text = element_text(size = 10),
      ratio = 0.6
    )

  print(gg)
}


data_for_ABC_computation <- inner_join(
  # we need inner_join to only keep the temperatures which are present in both the control and the TTT

  # correspond to the fit of the control condition
  RegionOfAcceptableUncertainty %>%
    inner_join(data.frame(Level_2 = "DMSO")) %>%
    filter(KeptTemperature) %>%
    dplyr::select(Level_1, Level_2, x, y),
  # correspond to the fit of the TTT condition
  RegionOfAcceptableUncertainty %>%
    anti_join(data.frame(Level_2 = "DMSO")) %>%
    filter(KeptTemperature) %>%
    dplyr::select(Level_1, Level_2, x, y),
  by = c("Level_1", "x"),
  suffix = c("_DMSO", "_TTT")
) %>%
  distinct() # remove repeated temperatures like the first and last.

######## final number of temperatures kept for ABC computation ########
N_temp <- data_for_ABC_computation %>%
  group_by(Level_1) %>%
  distinct() %>%
  summarise(N_temp_kept = n_distinct(x)) %>%
  mutate(
    N_temp_total = uniqueN(RegionOfAcceptableUncertainty_N_temp$N_temp_total)
  )

ComputeIntegral <- function(x, y) {
  tibble::tibble(SimpsonInt = sintegral(x$Temperature, x$Abs_diff)$int)
}
# install.packages('Bolstad2')
library(Bolstad2) #sintegral
######## absolute ABC   ########
abs_ABC <- data_for_ABC_computation %>%
  # compute the abs diff between these fits, at temperatures at which both fits are certain
  mutate(Abs_diff = abs(y_DMSO - y_TTT)) %>%
  dplyr::rename("Temperature" = "x") %>%
  group_by(Level_1) %>%
  # approximate the integral from it
  group_modify(~ ComputeIntegral(.x, .y), .keep = TRUE) %>%
  ungroup()

######## signed ABC   ########
signed_ABC <- data_for_ABC_computation %>%
  # compute the diff between these fits, at temperatures at which both fits are certain
  mutate(Abs_diff = y_TTT - y_DMSO) %>% # the name here is just for the use of ComputeIntegral
  dplyr::rename("Temperature" = "x") %>%
  group_by(Level_1) %>%
  # compute the integral from it
  group_modify(~ ComputeIntegral(.x, .y), .keep = TRUE) %>%
  ungroup()

######## combined signed and full   ########
ABC <- full_join(
  abs_ABC %>% dplyr::rename("abs_ABC" = "SimpsonInt"),
  signed_ABC %>% dplyr::rename("signed_ABC" = "SimpsonInt")
) %>%
  full_join(N_temp)

head(ABC)

p_values_subset <- full_join(
  p_values_subset,
  ABC,
  by = "Level_1"
)
# all_plots
head(p_values_subset)

p_values_subset <- p_values_subset %>%
  mutate(p_adj = p.adjust(pValue, method = "BH"))

p_values_df <- left_join(
  # remove the p-values computed from the subset
  p_values_subset %>%
    dplyr::select(
      -pValue,
      -size_null_distribution_approximation,
      -n_values_as_extreme_in_null_distribution_approximation
    ),
  # add p-values computed on the full dataset, using 10 samples per ID and 4772 IDs.
  p_values_subset %>%
    dplyr::select(
      Level_1,
      pValue,
      size_null_distribution_approximation,
      n_values_as_extreme_in_null_distribution_approximation,
      p_adj
    )
)


BaseSize = 12
th1 = 0.01
th2 = 0.05

VolcanoPlot <- ggplot(p_values_df) +
  geom_point(aes(
    x = log(abs_ABC),
    y = -log10(pValue),
    color = ifelse(p_adj <= th1, "th1", ifelse(p_adj <= th2, "th2", "n.s")),
    text = Level_1
  )) +
  scale_color_manual(
    name = "",
    values = c("th1" = "red", "th2" = "orange", "n.s" = "black"),
    labels = c(
      "th1" = paste0("BH adj. pval <= ", th1),
      "th2" = paste0("BH adj. pval <= ", th2),
      "n.s" = "n.s"
    )
  ) +
  ggtitle(
    "Results of GPMelt on the ATP 2019 dataset",
    subtitle = "Using a subset of 16 IDs"
  ) +
  theme_bw(base_size = BaseSize)

print(VolcanoPlot)

BaseSize = 12
library(ggrepel)
library(viridis)
ABC = merge(ABC, p_values_df[, .(p_adj, Level_1)], by = 'Level_1')
library(data.table)
ABC = as.data.table(ABC)
DiffInABC <- ggplot(ABC) +
  geom_abline(slope = 1, linetype = "dotted") +
  geom_abline(slope = -1, linetype = "dotted") +
  geom_point(
    aes(x = signed_ABC, y = abs_ABC, color = -log10(p_adj)),
    size = 3
  ) +
  geom_label_repel(
    data = ABC[between(signed_ABC, -2, 2) & abs_ABC > 2],
    aes(x = signed_ABC, y = abs_ABC, label = Level_1),
    min.segment.length = 10,
    size = BaseSize / 1.5
  ) +
  #  scale_color_viridis(direction=-1, option="magma", na.value ="black",
  #                     name =  latex2exp::TeX("$-\\log_{10}$(adj p-value)"))+
  # xlab(latex2exp::TeX("$ABC_{GPMelt}$"))+
  # ylab(latex2exp::TeX("$|ABC|_{GPMelt}$"))+
  theme_bw(base_size = BaseSize) +
  scico::scale_color_scico(palette = 'batlow')

DiffInABC
prot_interest = 'ZDHHC13_Q8IUH4'
all_plots[[prot_interest]][['Levels_1to2']] +
  ggtitle(
    prot_interest,
    subtitle = glue::glue(
      'padj = {p_values_df[Level_1  == prot_interest,p_adj] }'
    )
  ) +
  theme_bw() +
  theme(legend.position = 'bottom')


p_values_df = merge(
  p_values_df,
  starting_abundance[, .(Level_1, condition_diff)],
  by = 'Level_1'
)
p_values_df[p_adj < 0.01 & between(condition_diff, -0.01, 0.01)]
fwrite(p_values_df, 'GPMelt_dataset-wise_pvalues.csv')
