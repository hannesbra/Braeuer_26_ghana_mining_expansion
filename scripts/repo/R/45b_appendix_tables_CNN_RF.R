# ================================================================
# 45b_appendix_table.R
# LaTeX appendix tables for the CNN-RF onset consistency check.
# Reads the CSV outputs of 45_CNN_RF_onset_consistency.R (no re-run).
# Run after 00_config.R (for PATHS).
#
# Writes (P5_RESULTS/tables):
#   appendix_cnn_rf_onsets.tex    Table A: onset per mine
#   appendix_cnn_rf_summary.tex   Table B: agreement summary
#
# Format: Taylor & Francis `interact` class (\tbl{caption}{...}, \tabnote)
# with booktabs rules. If \toprule is undefined, add \usepackage{booktabs}.
# ================================================================

pacman::p_load(dplyr, readr, tidyr)

in_dir <- file.path(PATHS$out_validation, "cnn_rf_consistency", "tables")

# Must match CFG$main in script 45
MAIN <- list(radius = 5000, spec = 0.99, k = 3L, patience = 2L, rf_ha = 10,
             thr_mode = "per_year", comparison = "CNN T10 vs RF-P2 T10")
TOL      <- 2      # agreement tolerance (years)
EXT_TOL  <- 1      # tolerance vs documented start of operations
Y_MIN    <- 2000
T10_FRAC <- 0.10

cfg_main <- sprintf("r%g_s%g_k%d_p%d_ha%g_%s", MAIN$radius / 1000, MAIN$spec,
                    MAIN$k, MAIN$patience, MAIN$rf_ha, MAIN$thr_mode)

grid  <- read_csv(file.path(in_dir, "onset_grid_long.csv"),  show_col_types = FALSE)
agree <- read_csv(file.path(in_dir, "agreement_grid.csv"),   show_col_types = FALSE)
msum  <- read_csv(file.path(in_dir, "MAIN_consistency_summary.csv"), show_col_types = FALSE)
stopifnot(cfg_main %in% grid$config)

# ---- helpers ------------------------------------------------------
leq   <- function(x) paste0("$\\leq$", x)
fmt_y <- function(o, cens) ifelse(is.na(o), "n.d.", ifelse(cens, leq(o), as.character(o)))
fmt_n <- function(x) ifelse(is.na(x), "--", formatC(round(x), format = "d", big.mark = ","))
wilson <- function(k, n, z = 1.96) {
  if (n == 0) return(c(NA, NA))
  p <- k / n; d <- 1 + z^2 / n
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d,
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d)
}
pct <- function(x) sprintf("%.0f", 100 * x)

# ---- Table A: onset per mine ---------------------------------------
m <- grid |> filter(config == cfg_main) |>
  mutate(cnn_c = cnn_rel_cens %in% TRUE, rf_c = rf_rel_cens_p2 %in% TRUE)

n_valid_cfg <- grid |> filter(thr_mode == MAIN$thr_mode) |> distinct(config) |> nrow()
rob <- grid |> filter(thr_mode == MAIN$thr_mode, role == "positive") |>
  group_by(name) |>
  summarise(robust = mean(abs(cnn_rel - rf_rel_p2) <= TOL, na.rm = TRUE), .groups = "drop")

tabA <- m |> filter(role != "negative_control") |>
  left_join(rob, by = "name") |>
  mutate(
    group = case_when(role == "excluded" ~ 3L, in_window %in% TRUE ~ 1L, TRUE ~ 2L),
    mincom = ifelse(is.na(mc_first_prod), "--",
                    ifelse(mc_start_censored %in% TRUE, leq(mc_first_prod),
                           as.character(mc_first_prod))),
    rf  = fmt_y(rf_rel_p2, rf_c),
    cnn = fmt_y(cnn_rel,   cnn_c),
    delta = case_when(is.na(cnn_rel) | is.na(rf_rel_p2) ~ "--",
                      cnn_c & rf_c                      ~ "--",
                      cnn_c | rf_c ~ paste0(sprintf("$%+d$", as.integer(cnn_rel - rf_rel_p2)),
                                            "\\textsuperscript{c}"),
                      TRUE ~ sprintf("$%+d$", as.integer(cnn_rel - rf_rel_p2))),
    delta = gsub("\\$\\+0\\$", "$0$", delta),
    robust_txt = ifelse(role == "excluded" | is.na(robust), "--", pct(robust))) |>
  arrange(group, established, name)

row_tex <- function(d) sprintf("%s & %d & %s & %s & %s & %s & %s & %s & %s \\\\",
                               d$name, as.integer(d$established), d$mincom,
                               fmt_n(d$ref_ha), fmt_n(d$t10_threshold_ha),
                               d$rf, d$cnn, d$delta, d$robust_txt)
block <- function(g, label) {
  d <- tabA |> filter(group == g)
  if (!nrow(d)) return(character(0))
  c(sprintf("\\multicolumn{9}{@{}l}{\\textit{%s}}\\\\", label), row_tex(d))
}

ctrl_col <- grep("^ctrl_detected_pre2020_CNN T10$", names(msum), value = TRUE)
ctrl_txt <- if (length(ctrl_col) && !is.na(msum[[ctrl_col]][1]))
  ifelse(isTRUE(as.logical(msum[[ctrl_col]][1])),
         "was detected by the CNN before 2020 (false positive)",
         "was not detected by the CNN before 2020") else "could not be evaluated"

texA <- c(
  "\\begin{table}[htbp]",
  "\\tbl{Onset of industrial mining per mine: Random Forest (RF) classification versus an independent convolutional neural network (CNN).}",
  "{\\small\\setlength{\\tabcolsep}{3pt}\\begin{tabular}{@{}lrlrrccrr@{}}",
  "\\toprule",
  " & & & \\multicolumn{2}{c}{2019 extent (ha)} & \\multicolumn{2}{c}{Onset (T10)} & & \\\\",
  "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
  "Mine & Est. & MinCom\\textsuperscript{a} & Mined & Threshold & RF & CNN & $\\Delta$\\textsuperscript{b} & Robust (\\%)\\textsuperscript{d} \\\\",
  "\\midrule",
  block(1L, "Established 2001--2019"),
  "\\midrule",
  block(2L, "Established before 2000"),
  "\\midrule",
  block(3L, "Excluded from CNN training and agreement statistics\\textsuperscript{e}"),
  "\\bottomrule",
  "\\end{tabular}}",
  paste0(
    "\\tabnote{\\textit{Notes}: Est.\\ = documented start of operations. ",
    "Mined = 2019 RF Mining and Water cells inside the mine footprint (mining polygons within ",
    MAIN$radius / 1000, "~km of the mine, split between neighbouring mines); vegetation inside the polygons is not counted. ",
    "T10 = first year in which at least ", pct(T10_FRAC), "\\% of the 2019 mined extent (Threshold; minimum ",
    MAIN$rf_ha, "~ha) is detected in ", MAIN$patience + 1, " consecutive observed years. ",
    "CNN detection thresholds are calibrated per year to a ", sprintf("%.0f", 100 * (1 - MAIN$spec)),
    "\\% background false-positive rate. $\\leq$ = left-censored (already met in the first observed year). n.d.\\ = not detected.}"),
  "\\tabnote{\\textsuperscript{a}First reported Large-Scale production in Minerals Commission statistics (company-level data from 2013; $\\leq$2013 = producing in the first reported month).}",
  "\\tabnote{\\textsuperscript{b}CNN minus RF onset in years. \\textsuperscript{c}One of the two onsets is left-censored, so the difference is a bound.}",
  sprintf("\\tabnote{\\textsuperscript{d}Share of the %d sensitivity configurations (per-year thresholds) in which CNN and RF onsets differ by at most %d years.}",
          n_valid_cfg, TOL),
  paste0("\\tabnote{\\textsuperscript{e}Obuasi and Prestea are mainly underground operations; Awaso is the only bauxite mine. ",
         "The negative control Namdini (established 2021; outside the RF prediction area) ", ctrl_txt, ".}"),
  "\\label{tab:app_cnn_rf_onsets}",
  "\\end{table}")

# ---- Table B: agreement summary --------------------------------------
pos <- m |> filter(role == "positive")
win <- pos |> filter(in_window %in% TRUE)

agree_row <- function(d, label, subset_name) {
  ok  <- !is.na(d$cnn_rel) & !is.na(d$rf_rel_p2)
  k   <- sum(abs(d$cnn_rel[ok] - d$rf_rel_p2[ok]) <= TOL)
  n   <- sum(ok)
  unc <- ok & !d$cnn_c & !d$rf_c
  ci  <- wilson(k, n)
  g   <- agree |> filter(thr_mode == MAIN$thr_mode, comparison == MAIN$comparison,
                         tol == TOL, subset == subset_name)
  sprintf("%s & %d & %d/%d (%s) & %s--%s & %s & %s & %s--%s--%s \\\\",
          label, nrow(d), k, n, pct(k / n), pct(ci[1]), pct(ci[2]),
          ifelse(sum(unc), sprintf("%.1f", mean(abs(d$cnn_rel[unc] - d$rf_rel_p2[unc]))), "--"),
          ifelse(sum(unc), sprintf("%+.1f", mean(d$cnn_rel[unc] - d$rf_rel_p2[unc])), "--"),
          pct(min(g$consistent_rate, na.rm = TRUE)),
          pct(median(g$consistent_rate, na.rm = TRUE)),
          pct(max(g$consistent_rate, na.rm = TRUE)))
}

ext_row <- function(o, cens, label) {
  ok <- !is.na(o) & !cens
  k  <- sum(abs(o[ok] - win$established[ok]) <= EXT_TOL)
  sprintf("%s & \\multicolumn{3}{r}{%d/%d} & %s & %s & \\\\", label, k, nrow(win),
          ifelse(sum(ok), sprintf("%.1f", mean(abs(o[ok] - win$established[ok]))), "--"),
          ifelse(sum(ok), sub("^\\+0\\.0$", "0.0", sprintf("%+.1f", mean(o[ok] - win$established[ok]))), "--"))
}

texB <- c(
  "\\begin{table}[htbp]",
  "\\tbl{Agreement of mining onset (T10) between the RF classification, the CNN and documented start dates.}",
  "{\\small\\setlength{\\tabcolsep}{3pt}\\begin{tabular}{@{}lrrrrrr@{}}",
  "\\toprule",
  "\\multicolumn{7}{@{}l}{\\textit{(a) CNN versus RF}}\\\\",
  sprintf("Mines & $n$ & Within $\\pm$%d yr (\\%%) & 95\\%% CI\\textsuperscript{a} & MAE\\textsuperscript{b} & Bias\\textsuperscript{b} & Grid\\textsuperscript{c} \\\\", TOL),
  "\\midrule",
  agree_row(pos, "All included", "all positive"),
  agree_row(win, "Established 2001--2019", "window 2001-2019"),
  "\\midrule",
  "\\multicolumn{7}{@{}l}{\\textit{(b) Onset versus documented start of operations (mines established 2001--2019)}}\\\\",
  sprintf("Method & \\multicolumn{3}{r}{Within $\\pm$%d yr} & MAE & Bias & \\\\", EXT_TOL),
  "\\midrule",
  ext_row(win$rf_rel_p2, win$rf_c, "RF"),
  ext_row(win$cnn_rel,   win$cnn_c, "CNN"),
  "\\bottomrule",
  "\\end{tabular}}",
  "\\tabnote{\\textit{Notes}: Main specification as in Table~\\ref{tab:app_cnn_rf_onsets}. MAE = mean absolute error and bias = mean signed difference, both in years (CNN $-$ RF in panel a; onset $-$ start of operations in panel b).}",
  "\\tabnote{\\textsuperscript{a}Wilson score interval. \\textsuperscript{b}Pairs in which neither onset is left-censored.}",
  sprintf("\\tabnote{\\textsuperscript{c}Minimum--median--maximum share (\\%%) within $\\pm$%d years across the %d sensitivity configurations with per-year thresholds (clip 3, 5, 10~km or lease; specificity 0.95--0.995; 1--4 detected tiles; persistence 2--4 years; RF minimum 5--25~ha).}",
          TOL, n_valid_cfg),
  "\\label{tab:app_cnn_rf_summary}",
  "\\end{table}")

out_tex <- results_tab_dir()
writeLines(texA, file.path(out_tex, "appendix_cnn_rf_onsets.tex"))
writeLines(texB, file.path(out_tex, "appendix_cnn_rf_summary.tex"))
cat(texA, sep = "\n"); cat("\n"); cat(texB, sep = "\n")
cat("\nWritten:", file.path(out_tex, c("appendix_cnn_rf_onsets.tex", "appendix_cnn_rf_summary.tex")),
    sep = "\n  ")