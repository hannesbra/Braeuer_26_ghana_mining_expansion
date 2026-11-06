# ================================================================
# 03_Feature_importance_validation.R
# ----------------------------------------------------------------
# Feature audit for the NATIONAL model (10). Documents the band choice:
# the final model uses all 27 candidate bands (locked in 00_config.R),
# because they beat the selected top-10 set under spatial CV. The
# top-10 set written here is the M1 benchmark in 41.
#
# National (main text)
#   1. Permutation importance of all candidate bands (repeated RF fits).
#   2. Correlation families: Spearman |rho| > 0.7 (hierarchical clustering).
#   3. Selection: top_k bands by importance; with the "decorrelated"
#      rule at most one band per correlation family.
#   4. Spatial CV of competing national band sets on the shared folds:
#      all candidates, plain top_k, decorrelated top_k, literature union.
#   5. NESTED selection: steps 1-3 repeated inside every training fold,
#      so the reported accuracy of the selection procedure is not
#      inflated by choosing bands on the test data. Also gives how often
#      each band is chosen (selection stability).
#   6. Drop-family importance: CV with one correlation family removed.
#   7. Paired McNemar tests between the band sets.
#
# Zones (appendix)
#   8. Per-zone permutation importance, literature vs data-driven
#      overlap, and drop-family importance for zone-only models.
#
# Candidate bands exclude quality bands and, by default, the nightlight
# bands (sensor change in 2013; see 00_config.R).
#
# Outputs (PATHS$out_validation)
#   tables/national_selected_features.csv  top-10 set (M1 benchmark in 41;
#                                           read by 10 only if features = NULL)
#   tables/importance_full_kitchen_sink.csv   scope x feature
#   tables/feature_clusters_by_zone.csv       scope x feature -> family
#   tables/nested_selection_by_fold.csv       bands chosen per outer fold
#   tables/selection_stability.csv            selection frequency per band
#   tables/feature_set_cv_summary.csv         national band sets, spatial CV
#   tables/feature_set_mcnemar.csv            paired tests
#   tables/drop_cluster_importance.csv        scope x family deltas
#   tables/feature_set_comparison.csv         literature vs data-driven overlap
#   Figures go to P5_RESULTS/figures/03_feature_audit/:
#   importance_national.png
#   plots/selection_stability.png
#   plots/feature_set_comparison.png
#   plots/corr_heatmap_<scope>.png, dendrogram_<scope>.png
#   plots/drop_cluster_importance_mining_f1.png
#   plots/drop_cluster_importance_multimetric.png
#   plots/feature_overlap_heatmap.png
#   models/cv_folds_national.rds              shared folds (via 01_utils)
# ================================================================

source("R/00_config.R")
source("R/01_utils.R")

pacman::p_load(terra, sf, dplyr, tidyr, readr, ranger, blockCV,
               ggplot2, patchwork, corrplot, dendextend, parallel)

terraOptions(memfrac = 0.6, tempdir = tempdir())
set.seed(CV_CONFIG$seed)
start_time <- Sys.time()

# ================================================================
# 0. CONFIG
# ================================================================

AUDIT_CONFIG <- list(
  top_k             = NATIONAL_CONFIG$top_k,
  selection_rule    = NATIONAL_CONFIG$selection_rule,   # final rule for 10
  cluster_h         = NATIONAL_CONFIG$cluster_h,        # |rho| > 1 - h
  corr_method       = "spearman",
  n_imp_repeats     = 5,     # RF fits averaged for full-data importance
  n_imp_repeats_nested = 3,  # per outer fold (nested selection)
  run_zone_appendix = TRUE,  # per-zone importance, overlap, drop-family
  seed              = CV_CONFIG$seed
)
stopifnot(AUDIT_CONFIG$selection_rule %in% c("decorrelated", "top_k"))

out_dir <- PATHS$out_validation
ensure_output_dirs(out_dir, subdirs = c("tables", "models"))
plot_dir <- results_fig_dir("03_feature_audit")   # figures -> P5_RESULTS
N_THREADS <- max(1L, detectCores() - 1L)

SCOPE_LABELS <- c(national = "National", ZONE_LABELS)
scope_fac <- function(s) factor(s, levels = names(SCOPE_LABELS),
                                labels = unname(SCOPE_LABELS))

cat("================================================================\n")
cat("03 - FEATURE AUDIT (national selection + zone appendix)\n")
cat("================================================================\n")

# ================================================================
# 1. DATA AND SHARED FOLDS
# ================================================================

cat("\n[1/8] Training data and shared spatial folds\n")

df <- prepare_training_data()
y  <- as.integer(df$landcover_3class)
n  <- nrow(df)

candidates <- intersect(CANDIDATE_BANDS, names(df))
zero_var   <- candidates[vapply(df[candidates], function(v)
  isTRUE(stats::var(v, na.rm = TRUE) == 0), logical(1))]
candidates <- setdiff(candidates, zero_var)
cat(sprintf("  candidate bands: %d%s\n", length(candidates),
            if (length(zero_var)) paste0(" (constant, dropped: ",
                                         paste(zero_var, collapse = ", "), ")") else ""))
if (EXCLUDE_UNSTABLE_BANDS)
  cat("  nightlight bands excluded (EXCLUDE_UNSTABLE_BANDS = TRUE)\n")

folds    <- get_national_folds(df)
fold_ids <- folds$fold_ids

scopes <- c("national", if (AUDIT_CONFIG$run_zone_appendix) ZONE_ORDER)
scope_rows <- function(s) if (s == "national") seq_len(n) else which(df$zone == s)
cat("  rows per scope: ",
    paste(sprintf("%s = %d", scopes, vapply(scopes, function(s)
      length(scope_rows(s)), integer(1))), collapse = ", "), "\n", sep = "")

# ================================================================
# 2. HELPERS
# ================================================================

# Permutation importance averaged over n_rep weighted RF fits.
perm_importance <- function(rows, feats, n_rep, seed0) {
  if (length(rows) < 30 || length(unique(y[rows])) < 2) return(NULL)
  X   <- df[rows, feats, drop = FALSE]
  imp <- vapply(seq_len(n_rep), function(r) {
    m <- fit_rf_weighted(X, y[rows], seed = seed0 + r, threads = N_THREADS)
    m$variable.importance[feats]
  }, numeric(length(feats)))
  imp <- matrix(imp, nrow = length(feats), dimnames = list(feats, NULL))
  mu  <- rowMeans(imp)
  data.frame(feature       = feats,
             importance    = mu,
             importance_sd = if (n_rep > 1) apply(imp, 1, stats::sd) else NA_real_,
             rank          = rank(-mu, ties.method = "min"),
             stringsAsFactors = FALSE)
}

# Correlation families on the given rows.
cluster_features <- function(rows, feats,
                             method = AUDIT_CONFIG$corr_method,
                             h = AUDIT_CONFIG$cluster_h) {
  X    <- df[rows, feats, drop = FALSE]
  keep <- vapply(X, function(v) isTRUE(stats::var(v, na.rm = TRUE) > 0), logical(1))
  X    <- X[, keep, drop = FALSE]
  cor_m <- stats::cor(X, method = method, use = "pairwise.complete.obs")
  cor_m[is.na(cor_m)] <- 0
  hc <- stats::hclust(stats::as.dist(1 - abs(cor_m)), method = "average")
  list(features = colnames(X), cor = cor_m, hclust = hc,
       clusters = stats::cutree(hc, h = h))
}

# top_k by importance (positive importance only). "decorrelated":
# walk down the ranking and skip bands whose family is already chosen.
select_features <- function(imp, clusters, k, rule) {
  imp <- imp[imp$importance > 0 & imp$feature %in% names(clusters), , drop = FALSE]
  imp <- imp[order(imp$rank, -imp$importance), , drop = FALSE]
  if (rule == "top_k") return(head(imp$feature, k))
  picked <- character(0); used <- integer(0)
  for (f in imp$feature) {
    cl <- clusters[[f]]
    if (cl %in% used) next
    picked <- c(picked, f); used <- c(used, cl)
    if (length(picked) >= k) break
  }
  picked
}

# Spatial CV on the shared folds, restricted to `rows` (training and
# test). National scope uses the same seeds and settings as 10, so the
# selected set reproduces 10's OOF predictions.
cv_scope <- function(feats, rows = seq_len(n), seed_base = CV_CONFIG$seed) {
  feats    <- intersect(feats, names(df))
  prob     <- matrix(NA_real_, n, 3, dimnames = list(NULL, CLASS_LEVELS))
  in_scope <- seq_len(n) %in% rows
  if (length(feats) >= 2) {
    for (k in seq_len(CV_CONFIG$n_folds)) {
      te <- which(fold_ids == k & in_scope)
      tr <- intersect(folds$buffer_masks[[k]], rows)
      if (length(te) < 10 || length(tr) < 30 || length(unique(y[tr])) < 2) next
      m <- fit_rf_weighted(df[tr, feats, drop = FALSE], y[tr],
                           seed = seed_base + k, importance = "none")
      prob[te, ] <- predict_probs(m, df[te, feats, drop = FALSE])
    }
  }
  pred <- rep(NA_integer_, n)
  ok   <- !is.na(prob[, 1])
  pred[ok] <- argmax_class(prob[ok, , drop = FALSE])
  list(prob = prob, pred = pred)
}

cv_summary <- function(res, rows = seq_len(n)) {
  fm <- fold_metrics(y[rows], res$pred[rows], fold_ids[rows],
                     prob = res$prob[rows, , drop = FALSE])
  s  <- summarise_fold_metrics(fm)
  if (is.null(s)) return(NULL)
  s$n_predicted <- sum(!is.na(res$pred[rows]))
  s
}

scope_seed <- function(s) CV_CONFIG$seed +
  if (s == "national") 0 else 1000 * match(s, ZONE_ORDER)

# ================================================================
# 3. IMPORTANCE AND FAMILIES (national + zones)
# ================================================================

cat("\n[2/8] Permutation importance (all candidate bands)\n")

importance_by_scope <- list()
clusters_by_scope   <- list()
for (s in scopes) {
  rows <- scope_rows(s)
  cl   <- cluster_features(rows, candidates)
  imp  <- perm_importance(rows, cl$features, AUDIT_CONFIG$n_imp_repeats,
                          seed0 = scope_seed(s))
  if (is.null(imp)) { cat(sprintf("  %s: too few rows/classes, skipped\n", s)); next }
  imp$zone    <- s
  imp$cluster <- as.integer(cl$clusters[imp$feature])
  importance_by_scope[[s]] <- imp[order(imp$rank), ]
  clusters_by_scope[[s]]   <- cl
  cat(sprintf("  %-18s %d families; top %d: %s\n", s, length(unique(cl$clusters)),
              AUDIT_CONFIG$top_k,
              paste(head(importance_by_scope[[s]]$feature, AUDIT_CONFIG$top_k),
                    collapse = ", ")))
}

kitchen_sink_imp <- bind_rows(importance_by_scope) |>
  select(zone, feature, importance, importance_sd, rank, cluster)
write_csv(kitchen_sink_imp,
          file.path(out_dir, "tables", "importance_full_kitchen_sink.csv"))

cluster_summary <- bind_rows(lapply(names(clusters_by_scope), function(s)
  data.frame(zone = s, feature = names(clusters_by_scope[[s]]$clusters),
             cluster = as.integer(clusters_by_scope[[s]]$clusters)))) |>
  arrange(zone, cluster, feature)
write_csv(cluster_summary,
          file.path(out_dir, "tables", "feature_clusters_by_zone.csv"))

cat(sprintf("\n  National families (|rho| > %.1f):\n", 1 - AUDIT_CONFIG$cluster_h))
ncl <- clusters_by_scope$national$clusters
for (g in sort(unique(ncl)))
  cat(sprintf("    family %2d: %s\n", g, paste(names(ncl)[ncl == g], collapse = ", ")))

# ================================================================
# 4. NATIONAL SELECTION ON ALL DATA
# ================================================================

cat("\n[3/8] National band selection (all training data)\n")

imp_nat <- importance_by_scope$national
sel_sets <- list(
  top_k        = select_features(imp_nat, ncl, AUDIT_CONFIG$top_k, "top_k"),
  decorrelated = select_features(imp_nat, ncl, AUDIT_CONFIG$top_k, "decorrelated")
)
for (r in names(sel_sets))
  cat(sprintf("  %-12s (%d): %s\n", r, length(sel_sets[[r]]),
              paste(sel_sets[[r]], collapse = ", ")))
if (length(sel_sets$decorrelated) < AUDIT_CONFIG$top_k)
  cat(sprintf("  NOTE: only %d families with positive importance, so the decorrelated set has fewer than %d bands.\n",
              length(sel_sets$decorrelated), AUDIT_CONFIG$top_k))

# ================================================================
# 5. NESTED SELECTION (honest estimate + stability)
# ================================================================

cat("\n[4/8] Nested selection inside each outer fold\n")

nested <- lapply(c("top_k", "decorrelated"), function(rule)
  list(prob = matrix(NA_real_, n, 3, dimnames = list(NULL, CLASS_LEVELS)),
       sets = list()))
names(nested) <- c("top_k", "decorrelated")

for (k in seq_len(CV_CONFIG$n_folds)) {
  te <- which(fold_ids == k)
  tr <- folds$buffer_masks[[k]]
  if (length(te) < 10 || length(tr) < 30 || length(unique(y[tr])) < 2) next
  cl_k  <- cluster_features(tr, candidates)
  imp_k <- perm_importance(tr, cl_k$features, AUDIT_CONFIG$n_imp_repeats_nested,
                           seed0 = CV_CONFIG$seed + 100 * k)
  for (rule in names(nested)) {
    sel_k <- select_features(imp_k, cl_k$clusters, AUDIT_CONFIG$top_k, rule)
    m <- fit_rf_weighted(df[tr, sel_k, drop = FALSE], y[tr],
                         seed = CV_CONFIG$seed + k, importance = "none")
    nested[[rule]]$prob[te, ] <- predict_probs(m, df[te, sel_k, drop = FALSE])
    nested[[rule]]$sets[[as.character(k)]] <- sel_k
  }
  cat(sprintf("  fold %d: %s\n", k,
              paste(nested[[AUDIT_CONFIG$selection_rule]]$sets[[as.character(k)]],
                    collapse = ", ")))
}
for (rule in names(nested)) {
  p  <- nested[[rule]]$prob
  ok <- !is.na(p[, 1])
  nested[[rule]]$pred <- rep(NA_integer_, n)
  nested[[rule]]$pred[ok] <- argmax_class(p[ok, , drop = FALSE])
}

nested_by_fold <- bind_rows(lapply(names(nested), function(rule)
  bind_rows(lapply(names(nested[[rule]]$sets), function(k)
    data.frame(rule = rule, fold = as.integer(k),
               order = seq_along(nested[[rule]]$sets[[k]]),
               feature = nested[[rule]]$sets[[k]])))))
write_csv(nested_by_fold, file.path(out_dir, "tables", "nested_selection_by_fold.csv"))

n_outer <- length(nested[[1]]$sets)
stability <- nested_by_fold |>
  count(rule, feature, name = "n_folds_selected") |>
  mutate(selection_freq = n_folds_selected / n_outer) |>
  tidyr::complete(rule, feature = candidates,
                  fill = list(n_folds_selected = 0L, selection_freq = 0)) |>
  left_join(imp_nat |> select(feature, rank, cluster), by = "feature") |>
  mutate(in_full_data_set = mapply(function(r, f) f %in% sel_sets[[r]], rule, feature)) |>
  arrange(rule, desc(selection_freq), rank)
write_csv(stability, file.path(out_dir, "tables", "selection_stability.csv"))

# ================================================================
# 6. NATIONAL BAND SETS UNDER SPATIAL CV
# ================================================================

cat("\n[5/8] Spatial CV of national band sets (shared folds)\n")

literature_union <- intersect(unique(unlist(REGIONAL_FEATURES)), candidates)
set_defs <- list(
  all_candidates   = candidates,
  literature_union = literature_union,
  top_k            = sel_sets$top_k,
  decorrelated     = sel_sets$decorrelated
)

cv_sets <- lapply(names(set_defs), function(nm) {
  cat(sprintf("  %-18s (%d bands)\n", nm, length(set_defs[[nm]])))
  cv_scope(set_defs[[nm]])
})
names(cv_sets) <- names(set_defs)
cv_sets[[paste0("nested_", "top_k")]]        <- nested$top_k
cv_sets[[paste0("nested_", "decorrelated")]] <- nested$decorrelated

set_summary <- bind_rows(lapply(names(cv_sets), function(nm) {
  s <- cv_summary(cv_sets[[nm]])
  if (is.null(s)) return(NULL)
  cbind(config = nm,
        n_bands = if (nm %in% names(set_defs)) length(set_defs[[nm]]) else AUDIT_CONFIG$top_k,
        selection_on = ifelse(grepl("^nested_", nm), "training folds only",
                              ifelse(nm %in% c("top_k", "decorrelated"),
                                     "all data (optimistic)", "fixed a priori")),
        s)
}))
write_csv(set_summary, file.path(out_dir, "tables", "feature_set_cv_summary.csv"))
print(set_summary |> select(config, n_bands, selection_on, mean_accuracy,
                            ci95_lower, ci95_upper, mean_kappa, mean_mining_f1,
                            mean_log_loss),
      digits = 3, row.names = FALSE)

sel_name <- AUDIT_CONFIG$selection_rule
opt <- set_summary$mean_accuracy[set_summary$config == sel_name] -
  set_summary$mean_accuracy[set_summary$config == paste0("nested_", sel_name)]
cat(sprintf("\n  Selection optimism (%s, all-data minus nested accuracy): %+.3f\n",
            sel_name, opt))

# ---- Paired McNemar tests ---------------------------------------
cat("\n[6/8] Paired McNemar tests\n")
others <- setdiff(names(cv_sets), sel_name)
mcnemar_sets <- bind_rows(lapply(others, function(o)
  run_mcnemar(cv_sets[[sel_name]]$pred, cv_sets[[o]]$pred, y, sel_name, o)))
write_csv(mcnemar_sets, file.path(out_dir, "tables", "feature_set_mcnemar.csv"))
print(mcnemar_sets |> select(comparison, n_paired, a_only_correct, b_only_correct,
                             p_value, test_used, better),
      row.names = FALSE)

# ---- Final selected set for 10 ------------------------------------
final_set <- sel_sets[[sel_name]]
selected_tbl <- imp_nat |>
  filter(feature %in% final_set) |>
  mutate(order = match(feature, final_set)) |>
  arrange(order) |>
  left_join(stability |> filter(rule == sel_name) |>
              select(feature, nested_selection_freq = selection_freq),
            by = "feature") |>
  mutate(rule = sel_name, top_k = AUDIT_CONFIG$top_k,
         cluster_h = AUDIT_CONFIG$cluster_h, created = format(Sys.time())) |>
  select(order, feature, cluster, importance, importance_sd, rank,
         nested_selection_freq, rule, top_k, cluster_h, created)
write_csv(selected_tbl, NATIONAL_CONFIG$selected_file)

cat(sprintf("\n  Selected national band set -> %s\n", NATIONAL_CONFIG$selected_file))
cat("  To freeze it, paste into NATIONAL_CONFIG in 00_config.R:\n")
cat(sprintf("    features = c(%s),\n",
            paste(sprintf('"%s"', final_set), collapse = ", ")))

# ================================================================
# 7. DROP-FAMILY IMPORTANCE (national + zones)
# ================================================================
# Permutation importance splits credit between correlated bands, so
# each band of a strong family can look weak. Here a whole family is
# removed and the model re-fitted under spatial CV. Baseline = all
# candidate bands in the same scope.

cat("\n[7/8] Drop-family importance under spatial CV\n")

metric_list <- function(s) {
  if (is.null(s)) return(list(accuracy = NA_real_, mining_f1 = NA_real_,
                              water_f1 = NA_real_, vegetation_f1 = NA_real_,
                              log_loss = NA_real_, accuracy_sd = NA_real_,
                              mining_f1_sd = NA_real_, n_folds = 0L))
  list(accuracy = s$mean_accuracy, mining_f1 = s$mean_mining_f1,
       water_f1 = s$mean_water_f1, vegetation_f1 = s$mean_vegetation_f1,
       log_loss = s$mean_log_loss, accuracy_sd = s$sd_accuracy,
       mining_f1_sd = s$sd_mining_f1, n_folds = s$n_folds)
}

cluster_importance <- bind_rows(lapply(names(clusters_by_scope), function(s) {
  rows  <- scope_rows(s)
  cl    <- clusters_by_scope[[s]]
  feats <- cl$features
  cat(sprintf("  %s: baseline + %d families\n", s, length(unique(cl$clusters))))
  base_res <- if (s == "national") cv_sets$all_candidates else
    cv_scope(feats, rows, seed_base = scope_seed(s))
  base <- metric_list(cv_summary(base_res, rows))
  bind_rows(lapply(sort(unique(cl$clusters)), function(g) {
    members <- names(cl$clusters)[cl$clusters == g]
    drop <- metric_list(cv_summary(cv_scope(setdiff(feats, members), rows,
                                            seed_base = scope_seed(s)), rows))
    data.frame(
      zone = s, cluster = g, n_features = length(members),
      features = paste(members, collapse = ", "),
      base_accuracy = base$accuracy, base_mining_f1 = base$mining_f1,
      base_water_f1 = base$water_f1, base_vegetation_f1 = base$vegetation_f1,
      base_log_loss = base$log_loss,
      drop_accuracy = drop$accuracy, drop_mining_f1 = drop$mining_f1,
      drop_water_f1 = drop$water_f1, drop_vegetation_f1 = drop$vegetation_f1,
      drop_log_loss = drop$log_loss,
      # positive = the family matters (removing it hurts)
      delta_accuracy      = base$accuracy      - drop$accuracy,
      delta_mining_f1     = base$mining_f1     - drop$mining_f1,
      delta_water_f1      = base$water_f1      - drop$water_f1,
      delta_vegetation_f1 = base$vegetation_f1 - drop$vegetation_f1,
      delta_log_loss      = drop$log_loss      - base$log_loss,
      drop_accuracy_sd = drop$accuracy_sd, drop_mining_f1_sd = drop$mining_f1_sd,
      n_folds = drop$n_folds,
      stringsAsFactors = FALSE)
  }))
})) |> arrange(match(zone, names(SCOPE_LABELS)), desc(delta_mining_f1))
write_csv(cluster_importance,
          file.path(out_dir, "tables", "drop_cluster_importance.csv"))

for (s in unique(cluster_importance$zone)) {
  top <- cluster_importance |> filter(zone == s) |>
    arrange(desc(delta_mining_f1)) |> head(3)
  cat(sprintf("\n  %s (base mining F1 %.3f), top families:\n",
              toupper(s), top$base_mining_f1[1]))
  for (i in seq_len(nrow(top)))
    cat(sprintf("    family %2d  dF1 %+.4f  dAcc %+.4f  dLogLoss %+.4f | %s\n",
                top$cluster[i], top$delta_mining_f1[i], top$delta_accuracy[i],
                top$delta_log_loss[i], top$features[i]))
}

# ================================================================
# 8. ZONE APPENDIX: literature vs data-driven overlap
# ================================================================

cat("\n[8/8] Literature vs data-driven overlap\n")

overlap_row <- function(scope, lit, dat) {
  inter <- intersect(lit, dat)
  rho <- NA_real_
  if (length(inter) >= 3) {
    ranks <- kitchen_sink_imp |> filter(zone == scope)
    rho <- stats::cor(match(inter, lit), ranks$rank[match(inter, ranks$feature)],
                      method = "spearman")
  }
  data.frame(zone = scope, n_literature = length(lit), n_data_driven = length(dat),
             n_overlap = length(inter), jaccard = length(inter) / length(union(lit, dat)),
             spearman_rho = rho,
             overlap_features = paste(inter, collapse = ", "),
             only_literature  = paste(setdiff(lit, dat), collapse = ", "),
             only_data_driven = paste(setdiff(dat, lit), collapse = ", "),
             stringsAsFactors = FALSE)
}

zone_topk <- lapply(setdiff(names(importance_by_scope), "national"), function(z)
  head(importance_by_scope[[z]]$feature, AUDIT_CONFIG$top_k))
names(zone_topk) <- setdiff(names(importance_by_scope), "national")

set_comparison <- bind_rows(
  overlap_row("national", literature_union, final_set),
  bind_rows(lapply(names(zone_topk), function(z)
    overlap_row(z, intersect(REGIONAL_FEATURES[[z]], candidates), zone_topk[[z]])))
)
write_csv(set_comparison, file.path(out_dir, "tables", "feature_set_comparison.csv"))
print(set_comparison |> select(zone, n_literature, n_data_driven, n_overlap,
                               jaccard, spearman_rho),
      digits = 3, row.names = FALSE)

# ================================================================
# 9. PLOTS
# ================================================================

cat("\nPlots\n")

theme_pub <- function(base = 10) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
          strip.background = element_rect(fill = "grey96", colour = NA),
          strip.text       = element_text(face = "bold"),
          plot.title       = element_text(face = "bold", size = base + 1),
          plot.subtitle    = element_text(colour = "grey30", size = base - 1),
          plot.background  = element_rect(fill = "white", colour = NA))
}
png_dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
pdf_dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
# Paper figures carry no title/subtitle; captions are written in LaTeX.
PAPER_FIGURES <- TRUE
save_both <- function(fig, path, width, height) {
  if (PAPER_FIGURES) fig <- fig + labs(title = NULL, subtitle = NULL)
  ggsave(path, fig, device = png_dev, width = width, height = height,
         units = "in", dpi = 600, bg = "white")
  ggsave(sub("\\.png$", ".pdf", path), fig, device = pdf_dev,
         width = width, height = height, units = "in", bg = "white")
}

# (a) National importance, selected bands highlighted
imp_plot <- imp_nat |>
  mutate(selected = feature %in% final_set,
         feature  = factor(feature, levels = feature[order(importance)]))
fig_imp <- ggplot(imp_plot, aes(importance, feature, fill = selected)) +
  geom_col(width = 0.75, colour = "grey20", linewidth = 0.2) +
  geom_errorbar(aes(xmin = pmax(importance - importance_sd, 0),
                    xmax = importance + importance_sd),
                width = 0.3, linewidth = 0.3, colour = "grey30",
                orientation = "y", na.rm = TRUE) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#BBBBBB"),
                    labels = c(`TRUE` = "selected", `FALSE` = "not selected"),
                    name = NULL) +
  labs(x = sprintf("Permutation importance (mean \u00b1 SD, %d RF fits)",
                   AUDIT_CONFIG$n_imp_repeats), y = NULL,
       title = "National band importance",
       subtitle = sprintf("Selected: %s rule, k = %d", sel_name, AUDIT_CONFIG$top_k)) +
  theme_pub() + theme(legend.position = "bottom", panel.grid.major.y = element_blank())
save_both(fig_imp, file.path(plot_dir, "importance_national.png"),
          7, max(4, 0.22 * nrow(imp_plot) + 1.5))

# (b) Selection stability across outer folds
stab_plot <- stability |>
  filter(rule == sel_name, selection_freq > 0 | in_full_data_set) |>
  mutate(feature = factor(feature, levels = unique(feature[order(selection_freq, -rank)])))
fig_stab <- ggplot(stab_plot, aes(selection_freq, feature, fill = in_full_data_set)) +
  geom_col(width = 0.7, colour = "grey20", linewidth = 0.2) +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#BBBBBB"),
                    labels = c(`TRUE` = "in final set", `FALSE` = "not in final set"),
                    name = NULL) +
  labs(x = sprintf("Share of %d outer folds in which the band was selected", n_outer),
       y = NULL, title = "Selection stability (nested)") +
  theme_pub() + theme(legend.position = "bottom", panel.grid.major.y = element_blank())
save_both(fig_stab, file.path(plot_dir, "selection_stability.png"),
          6.5, max(3.5, 0.25 * nrow(stab_plot) + 1.5))

# (c) National band sets: accuracy (95% CI) and mining F1 (SD)
cfg_levels <- c("all_candidates", "literature_union", "top_k", "decorrelated",
                "nested_top_k", "nested_decorrelated")
cmp_plot <- bind_rows(
  set_summary |> transmute(config, metric = "Accuracy", mean = mean_accuracy,
                           lo = ci95_lower, hi = ci95_upper, selection_on),
  set_summary |> transmute(config, metric = "Mining F1", mean = mean_mining_f1,
                           lo = pmax(mean_mining_f1 - sd_mining_f1, 0),
                           hi = pmin(mean_mining_f1 + sd_mining_f1, 1), selection_on)
) |> mutate(config = factor(config, levels = cfg_levels))
cfg_labels <- c(all_candidates      = "All 27 bands",
                literature_union    = "Literature bands",
                top_k               = "Top-10",
                decorrelated        = "Decorrelated top-10",
                nested_top_k        = "Top-10\n(nested)",
                nested_decorrelated = "Decorrelated top-10\n(nested)")
fig_cfg <- ggplot(cmp_plot, aes(config, mean, fill = selection_on)) +
  geom_col(width = 0.65, colour = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", mean), y = hi), vjust = -0.4, size = 2.8) +
  facet_wrap(~ metric, nrow = 1) +
  scale_x_discrete(labels = cfg_labels) +
  scale_fill_manual(values = c("fixed a priori" = "#BBBBBB",
                               "all data (optimistic)" = "#CCBB44",
                               "training folds only" = "#4477AA"),
                    name = "Bands selected on") +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = NULL, y = "Spatial CV (accuracy: mean \u00b1 95% CI; F1: mean \u00b1 SD)",
       title = "National band sets under spatial cross-validation",
       subtitle = sprintf("Shared folds, buffer %.1f km, n = %d",
                          folds$buffer_m / 1000, n)) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
save_both(fig_cfg, file.path(plot_dir, "feature_set_comparison.png"), 10, 4.8)

# (d) Correlation heatmaps and dendrograms per scope
for (s in names(clusters_by_scope)) {
  cl <- clusters_by_scope[[s]]
  n_clust <- length(unique(cl$clusters))
  png(file.path(plot_dir, paste0("corr_heatmap_", s, ".png")),
      width = 1800, height = 1800, res = 200, bg = "white")
  corrplot::corrplot(cl$cor, method = "color", type = "upper",
                     order = "hclust", hclust.method = "average",
                     addrect = n_clust, rect.col = "black", rect.lwd = 1.5,
                     tl.col = "black", tl.srt = 45, tl.cex = 0.8,
                     addCoef.col = "grey25", number.cex = 0.5,
                     title = if (PAPER_FIGURES) "" else
                       sprintf("Spearman correlation - %s (families at |rho| > %.1f)",
                               SCOPE_LABELS[[s]], 1 - AUDIT_CONFIG$cluster_h),
                     mar = c(0, 0, if (PAPER_FIGURES) 0.5 else 2.5, 0))
  dev.off()
  
  png(file.path(plot_dir, paste0("dendrogram_", s, ".png")),
      width = 1800, height = 1100, res = 200, bg = "white")
  dnd <- dendextend::color_branches(stats::as.dendrogram(cl$hclust), k = n_clust)
  par(mar = c(8, 4, 4, 2))
  plot(dnd, main = sprintf("Band families - %s (cut at h = %.2f)",
                           SCOPE_LABELS[[s]], AUDIT_CONFIG$cluster_h),
       ylab = "1 - |rho|")
  abline(h = AUDIT_CONFIG$cluster_h, col = "red", lty = 2, lwd = 1.5)
  dev.off()
}

# (e) Drop-family importance
make_family_label <- function(features_str) {
  m <- sort(strsplit(features_str, ",\\s*")[[1]])
  n_m <- length(m)
  if (n_m == 1) {
    pretty1 <- c(IRON_OXIDE = "Iron oxide", WATER_PROB = "Water probability",
                 NDVI_ANOMALY = "NDVI anomaly", NIR_corr = "NIR correlation",
                 NIR_STD_DEV = "NIR std dev", NIR_contrast = "NIR contrast",
                 NIR_var = "NIR variance", ELEVATION = "Elevation",
                 ASPECT = "Aspect", SLOPE = "Slope", TRI = "TRI",
                 CLAY_INDEX = "Clay index", BSI = "BSI", MNDWI = "MNDWI",
                 NDBI = "NDBI", ENTROPY = "Entropy",
                 NIGHTLIGHT = "Nightlight", NIGHTLIGHT_NORM = "Nightlight (norm)",
                 NIGHTLIGHT_DENSITY = "Nightlight density", IS_LIT = "Is-lit")
    return(unname(if (m %in% names(pretty1)) pretty1[m] else m))
  }
  exact_themes <- list(
    "Bare soil (BSI+CLAY)"   = c("BSI", "CLAY_INDEX"),
    "Topography (slope+TRI)" = c("SLOPE", "TRI"),
    "SWIR + built-up"        = c("SWIR1", "SWIR2", "AWEI_NSH", "NDBI"))
  for (lab in names(exact_themes)) if (setequal(m, exact_themes[[lab]])) return(lab)
  subset_themes <- list(
    "Nightlight"  = c("NIGHTLIGHT", "NIGHTLIGHT_NORM", "IS_LIT", "NIGHTLIGHT_DENSITY"),
    "NIR texture" = c("NIR_STD_DEV", "NIR_contrast", "NIR_var", "NIR_corr"))
  for (lab in names(subset_themes)) if (all(m %in% subset_themes[[lab]])) return(lab)
  tag_rules <- list(Visible = c("BLUE", "GREEN", "RED"), NIR = "NIR",
                    SWIR = c("SWIR1", "SWIR2"),
                    greenness = c("NDVI", "NDWI", "NDVI_ANOMALY"),
                    `water idx.` = c("MNDWI", "AWEI_NSH", "AWEI_SH", "WATER_PROB"),
                    color = c("CI", "BI"))
  tags <- names(tag_rules)[vapply(tag_rules, function(t) any(m %in% t), logical(1))]
  if (length(tags)) return(sprintf("%s (n=%d)", paste(tags, collapse = " + "), n_m))
  if (n_m <= 3) return(paste(m, collapse = " + "))
  sprintf("%s + %d more", paste(m[1:2], collapse = " + "), n_m - 2)
}

plot_data_mF1 <- cluster_importance |>
  mutate(family_label = vapply(features, make_family_label, character(1))) |>
  filter(!is.na(delta_mining_f1)) |>
  group_by(zone) |>
  mutate(family_label = make.unique(family_label, sep = " #")) |>
  ungroup() |>
  arrange(match(zone, names(SCOPE_LABELS)), delta_mining_f1) |>
  mutate(key = factor(paste(zone, family_label, sep = " | "),
                      levels = unique(paste(zone, family_label, sep = " | "))))
strip_scope <- function(b) sub(".*\\| ", "", b)
n_scopes <- length(unique(plot_data_mF1$zone))

theme_bars <- function(base = 11) {
  theme_pub(base) +
    theme(panel.grid.major.y = element_blank(),
          strip.text.y = element_text(angle = 0, hjust = 0),
          panel.spacing = unit(0.6, "lines"))
}

lim <- max(abs(plot_data_mF1$delta_mining_f1), 1e-3, na.rm = TRUE) * 1.25
fig_mF1 <- ggplot(plot_data_mF1,
                  aes(y = key, x = delta_mining_f1, fill = delta_mining_f1 > 0)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey20") +
  geom_col(colour = "grey15", linewidth = 0.25, width = 0.72) +
  geom_text(aes(label = sprintf("%+.3f", delta_mining_f1),
                hjust = ifelse(delta_mining_f1 >= 0, -0.15, 1.15)),
            size = 3, colour = "grey15") +
  facet_grid(rows = vars(scope_fac(zone)), scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = strip_scope) +
  scale_x_continuous(limits = c(-lim, lim), expand = expansion(mult = 0.02)) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#BBBBBB"), guide = "none") +
  labs(y = NULL, x = "\u0394 Mining F1 (baseline \u2212 without family)",
       title = "Drop-family importance for mining detection",
       subtitle = sprintf("Families = correlation clusters at |\u03c1| > %.1f. Positive: removing the family hurts.",
                          1 - AUDIT_CONFIG$cluster_h)) +
  theme_bars(11)
save_both(fig_mF1, file.path(plot_dir, "drop_cluster_importance_mining_f1.png"),
          9, max(3.5, 0.3 * nrow(plot_data_mF1) + 0.5 * n_scopes + 1.6))

plot_data_multi <- plot_data_mF1 |>
  select(zone, key, delta_accuracy, delta_mining_f1, delta_log_loss) |>
  tidyr::pivot_longer(starts_with("delta_"), names_to = "metric", values_to = "delta") |>
  mutate(metric = factor(recode(metric, delta_accuracy = "Accuracy",
                                delta_mining_f1 = "Mining F1",
                                delta_log_loss = "Log-loss"),
                         levels = c("Accuracy", "Mining F1", "Log-loss")))
fig_multi <- ggplot(plot_data_multi, aes(y = key, x = delta, fill = delta > 0)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey20") +
  geom_col(colour = "grey15", linewidth = 0.2, width = 0.72) +
  facet_grid(rows = vars(scope_fac(zone)), cols = vars(metric),
             scales = "free", space = "free_y") +
  scale_y_discrete(labels = strip_scope) +
  scale_x_continuous(n.breaks = 4, expand = expansion(mult = 0.08)) +
  scale_fill_manual(values = c(`TRUE` = "#4477AA", `FALSE` = "#BBBBBB"), guide = "none") +
  labs(y = NULL,
       x = "\u0394 metric (baseline \u2212 without family; log-loss: without \u2212 baseline)",
       title = "Drop-family importance across three metrics",
       subtitle = "Mining F1 and log-loss react to probability shifts that accuracy misses.") +
  theme_bars(10)
save_both(fig_multi, file.path(plot_dir, "drop_cluster_importance_multimetric.png"),
          12, max(4.5, 0.28 * nrow(plot_data_mF1) + 0.5 * n_scopes + 1.8))

# (f) Band membership: literature vs data-driven (zones) and national
overlap_df <- bind_rows(
  data.frame(zone = "national", feature = union(literature_union, final_set),
             in_lit = NA, in_data = NA),
  bind_rows(lapply(names(zone_topk), function(z) {
    lit <- intersect(REGIONAL_FEATURES[[z]], candidates)
    data.frame(zone = z, feature = union(lit, zone_topk[[z]]), in_lit = NA, in_data = NA)
  }))
) |>
  rowwise() |>
  mutate(in_lit  = feature %in% (if (zone == "national") literature_union else
    intersect(REGIONAL_FEATURES[[zone]], candidates)),
    in_data = feature %in% (if (zone == "national") final_set else zone_topk[[zone]])) |>
  ungroup() |>
  mutate(category = case_when(in_lit & in_data ~ "both",
                              in_lit ~ "literature only",
                              TRUE ~ "data-driven only"))
fig_overlap <- ggplot(overlap_df, aes(scope_fac(zone), feature, fill = category)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("both" = "#4477AA", "literature only" = "#228833",
                               "data-driven only" = "#CCBB44"), name = NULL) +
  labs(x = NULL, y = NULL,
       title = "Band membership: literature vs data-driven",
       subtitle = "National: literature union vs selected set. Zones: zone list vs zone top-k.") +
  theme_pub() + theme(panel.grid = element_blank(), legend.position = "bottom")
save_both(fig_overlap, file.path(plot_dir, "feature_overlap_heatmap.png"), 7.5, 8)

# ================================================================
# DONE
# ================================================================

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat("\n================================================================\n")
cat(sprintf("Feature audit complete in %.1f min\n", as.numeric(elapsed)))
cat("================================================================\n")
cat(sprintf("Selected set (%s): %s\n", sel_name, paste(final_set, collapse = ", ")))
cat(sprintf("Nested CV accuracy of this procedure: %.3f\n",
            set_summary$mean_accuracy[set_summary$config == paste0("nested_", sel_name)]))
cat("Next: 10_Training_Prediction.R\n")