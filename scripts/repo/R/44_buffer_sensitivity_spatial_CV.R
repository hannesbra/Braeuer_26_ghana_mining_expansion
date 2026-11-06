# ================================================================
# 44_buffer_sensitivity.R
# ----------------------------------------------------------------
# Sensitivity of the spatial-CV estimates to the variogram-derived
# buffer and to the block size. Run after 10 (needs the shared fold
# cache); independent of 20-50.
#
#   [1] Variogram: exact 95 % practical range of the fitted model
#       (from the fitted curve, not the x4 Matérn rule), and how the
#       range changes with lag cutoff and model family.
#   [2] Buffer sweep on the shared fold ids (0 - 50 km).
#   [3] Block-size sweep (new blockCV folds, several fold seeds).
#   [4] Test-to-train nearest-neighbour distances per configuration,
#       compared with prediction-location-to-train distances
#       (Milà et al. 2022): the configuration whose CV distances
#       resemble the prediction distances is the honest one.
#
# Each CV configuration is repeated with N_REP forest seeds so RF noise
# can be separated from the effect of the buffer.
#
# Outputs (PATHS$out_validation):
#   tables/buffer_sensitivity_variogram.csv
#   tables/buffer_sensitivity.csv            per configuration x rep
#   tables/buffer_sensitivity_summary.csv    mean over reps
#   tables/buffer_sensitivity_nn_distance.csv
#   P5_RESULTS/figures/44_buffer_sensitivity/buffer_sensitivity.png
# ================================================================

source("R/00_config.R")
source("R/01_utils.R")

pacman::p_load(sf, dplyr, tidyr, readr, ranger, blockCV, gstat, automap,
               ggplot2, patchwork, parallel)

sf::sf_use_s2(FALSE)
N_THREADS <- max(1L, detectCores() - 1L)
out_dir   <- PATHS$out_validation
ensure_output_dirs(out_dir, subdirs = c("tables"))
plot_dir  <- results_fig_dir("44_buffer_sensitivity")

BUFFERS_KM    <- c(0, 0.5, 1, NA, 5, 10, 20, 50)   # NA = variogram buffer
BLOCKS_KM     <- c(10, NA, 50, 100)                 # NA = block used in 10
FOLD_SEEDS    <- CV_CONFIG$seed + c(0, 101, 202)    # block sweep only
N_REP         <- 3                                  # forest seeds per config
VGM_CUTOFF_KM <- c(5, 10, 20, 40)
VGM_MODELS    <- c("Sph", "Exp", "Mat", "Gau")

cat("================================================================\n")
cat("44 - BUFFER AND BLOCK-SIZE SENSITIVITY\n")
cat("================================================================\n")

# ================================================================
# DATA AND SHARED FOLDS
# ================================================================
df       <- prepare_training_data()
y        <- as.integer(df$landcover_3class)
n        <- nrow(df)
features <- read_selected_features()$features
folds    <- get_national_folds(df)
data_sf  <- st_as_sf(df, coords = c("x_utm", "y_utm"), crs = GHANA_CRS, remove = FALSE)
data_sf$mining_binary <- as.numeric(y == CLASS_IDS[["Mining"]])

buf_used   <- folds$buffer_m
block_used <- folds$block_size
BUFFERS_KM[is.na(BUFFERS_KM)] <- buf_used / 1000
BLOCKS_KM[is.na(BLOCKS_KM)]   <- block_used / 1000
cat(sprintf("  %d points, %d bands; shared folds: buffer %.0f m, block %.0f m\n",
            n, length(features), buf_used, block_used))

# ================================================================
# [1] VARIOGRAM
# ================================================================
cat("\n[1/4] Variogram range\n")

# Distance at which the fitted model reaches nugget + 95 % of the
# structured (partial) sill.
exact_practical_range <- function(var_model, maxdist) {
  vl  <- gstat::variogramLine(var_model, maxdist = maxdist, n = 20000)
  nug <- var_model$psill[1]; ps <- sum(var_model$psill) - nug
  hit <- which(vl$gamma >= nug + 0.95 * ps)
  if (!length(hit)) NA_real_ else vl$dist[hit[1]]
}
rule_range <- function(model, raw) switch(model, Exp = 3 * raw, Gau = sqrt(3) * raw,
                                          Sph = raw, Mat = 4 * raw, raw)

fit_one <- function(cutoff_m, models) {
  f <- tryCatch(automap::autofitVariogram(mining_binary ~ 1, input_data = data_sf,
                                          cutoff = cutoff_m, model = models),
                error = function(e) NULL)
  if (is.null(f)) return(NULL)
  vm <- f$var_model
  data.frame(cutoff_km  = cutoff_m / 1000,
             model      = as.character(vm$model[2]),
             kappa      = vm$kappa[2],
             nugget     = vm$psill[1],
             sill       = sum(vm$psill),
             nugget_ratio = vm$psill[1] / sum(vm$psill),
             raw_range_m  = vm$range[2],
             rule_range_m = rule_range(as.character(vm$model[2]), vm$range[2]),
             exact_range_m = exact_practical_range(vm, maxdist = 50 * vm$range[2] + 1),
             sserr      = f$sserr)
}

vgm_rows <- list()
for (co in VGM_CUTOFF_KM) {
  r <- fit_one(co * 1000, VGM_MODELS)
  if (!is.null(r)) vgm_rows[[length(vgm_rows) + 1]] <- cbind(r, choice = "best")
  for (m in VGM_MODELS) {
    r <- fit_one(co * 1000, m)
    if (!is.null(r)) vgm_rows[[length(vgm_rows) + 1]] <- cbind(r, choice = "forced")
  }
}
vgm_tab <- bind_rows(vgm_rows)

# Exact range of the model actually used for the buffer
vi <- folds$vgm_info
if (!is.null(vi$fit)) {
  vm_used    <- vi$fit$var_model
  exact_used <- exact_practical_range(vm_used, maxdist = 50 * vi$raw_range + 1)
  vgm_tab <- bind_rows(data.frame(
    cutoff_km = NA, model = vi$model, kappa = vm_used$kappa[2],
    nugget = vi$nugget, sill = vi$sill, nugget_ratio = vi$nugget / vi$sill,
    raw_range_m = vi$raw_range, rule_range_m = vi$practical_range,
    exact_range_m = exact_used, sserr = vi$fit$sserr, choice = "used_in_10"), vgm_tab)
  cat(sprintf("  used: %s (kappa %.2f), raw %.0f m; x-rule %.0f m; exact 95%% range %.0f m\n",
              vi$model, vm_used$kappa[2], vi$raw_range, vi$practical_range, exact_used))
  cat(sprintf("  nugget/sill = %.2f -> %.0f %% of the variance is unstructured at the sample spacing\n",
              vi$nugget / vi$sill, 100 * vi$nugget / vi$sill))
} else exact_used <- NA_real_

print(vgm_tab |> mutate(across(where(is.numeric), ~ round(., 3))), row.names = FALSE)
write_csv(vgm_tab, file.path(out_dir, "tables", "buffer_sensitivity_variogram.csv"))

# Add the exact range to the buffer sweep if it differs noticeably
if (is.finite(exact_used) && abs(exact_used - buf_used) > 250)
  BUFFERS_KM <- sort(unique(c(BUFFERS_KM, round(exact_used / 1000, 3))))
# and the range of the best fit at the largest cutoff
best_far <- vgm_tab |> filter(choice == "best", cutoff_km == max(VGM_CUTOFF_KM))
if (nrow(best_far) && is.finite(best_far$exact_range_m))
  BUFFERS_KM <- sort(unique(c(BUFFERS_KM, round(best_far$exact_range_m / 1000, 3))))

# ================================================================
# CV ENGINE
# ================================================================
buffer_masks_for <- function(fold_ids, buffer_m) {
  lapply(seq_len(CV_CONFIG$n_folds), function(k) {
    te <- which(fold_ids == k); tr <- which(fold_ids != k)
    if (buffer_m <= 0 || !length(te)) return(tr)
    near <- lengths(st_is_within_distance(data_sf[tr, ], data_sf[te, ],
                                          dist = buffer_m)) > 0
    tr[!near]
  })
}

# Nearest training point for every test point, per fold
nn_test_train <- function(fold_ids, masks) {
  unlist(lapply(seq_len(CV_CONFIG$n_folds), function(k) {
    te <- which(fold_ids == k); tr <- masks[[k]]
    if (!length(te) || !length(tr)) return(NULL)
    i <- st_nearest_feature(data_sf[te, ], data_sf[tr, ])
    as.numeric(st_distance(data_sf[te, ], data_sf[tr, ][i, ], by_element = TRUE))
  }))
}

# Local ranger prediction, so a predict_probs() left in the session by
# another script (e.g. 45's CNN version) cannot be picked up here.
rf_probs <- function(model, X) {
  p   <- ranger:::predict.ranger(model, X, num.threads = N_THREADS)$predictions
  out <- matrix(0, nrow(p), length(CLASS_LEVELS), dimnames = list(NULL, CLASS_LEVELS))
  cols <- intersect(colnames(p), CLASS_LEVELS)
  out[, cols] <- p[, cols, drop = FALSE]
  out
}

run_cv <- function(fold_ids, masks, rep) {
  prob <- matrix(NA_real_, n, 3, dimnames = list(NULL, CLASS_LEVELS))
  for (k in seq_len(CV_CONFIG$n_folds)) {
    te <- which(fold_ids == k); tr <- masks[[k]]
    if (length(te) < 10 || length(tr) < 30 || length(unique(y[tr])) < 2) next
    m <- fit_rf_weighted(df[tr, features, drop = FALSE], y[tr],
                         seed = CV_CONFIG$seed + k + 10000 * (rep - 1),
                         threads = N_THREADS, importance = "none")
    prob[te, ] <- rf_probs(m, df[te, features, drop = FALSE])
  }
  pred <- rep(NA_integer_, n); ok <- !is.na(prob[, 1])
  pred[ok] <- argmax_class(prob[ok, , drop = FALSE])
  summarise_fold_metrics(fold_metrics(y, pred, fold_ids, prob))
}

evaluate_config <- function(sweep, fold_ids, buffer_m, block_m, fold_seed) {
  masks  <- buffer_masks_for(fold_ids, buffer_m)
  n_tr   <- sapply(masks, length)
  nn     <- nn_test_train(fold_ids, masks)
  out <- lapply(seq_len(N_REP), function(r) {
    s <- run_cv(fold_ids, masks, r)
    cbind(data.frame(sweep = sweep, buffer_km = buffer_m / 1000,
                     block_km = block_m / 1000, fold_seed = fold_seed, rep = r,
                     n_train_mean = mean(n_tr), n_train_min = min(n_tr),
                     pct_train_removed = 100 * (1 - mean(n_tr) /
                                                  mean(sapply(seq_len(CV_CONFIG$n_folds),
                                                              function(k) sum(fold_ids != k)))),
                     nn_median_km = median(nn) / 1000), s)
  })
  res <- bind_rows(out)
  cat(sprintf("  %-6s buffer %6.2f km  block %6.1f km  seed %d  n_train %4.0f  acc %.3f  mining F1 %.3f  log-loss %.3f\n",
              sweep, buffer_m / 1000, block_m / 1000, fold_seed, mean(n_tr),
              mean(res$mean_accuracy), mean(res$mean_mining_f1), mean(res$mean_log_loss)))
  list(metrics = res,
       nn = data.frame(sweep = sweep, buffer_km = buffer_m / 1000,
                       block_km = block_m / 1000, fold_seed = fold_seed, nn_m = nn))
}

# ================================================================
# [2] BUFFER SWEEP (shared folds)
# ================================================================
cat("\n[2/4] Buffer sweep on the shared folds\n")
buf_res <- lapply(BUFFERS_KM, function(b)
  evaluate_config("buffer", folds$fold_ids, b * 1000, block_used, CV_CONFIG$seed))

# ================================================================
# [3] BLOCK-SIZE SWEEP (variogram buffer)
# ================================================================
cat("\n[3/4] Block-size sweep\n")
blk_res <- list()
for (bk in BLOCKS_KM) for (fs in FOLD_SEEDS) {
  fo <- suppressWarnings(create_spatial_cv_folds(data_sf, buffer_distance = 0,
                                                 block_size = bk * 1000, seed = fs))
  if (length(unique(fo$fold_ids)) < CV_CONFIG$n_folds) {
    cat(sprintf("  block %.1f km seed %d: fewer than %d folds; skipped\n",
                bk, fs, CV_CONFIG$n_folds)); next
  }
  blk_res[[length(blk_res) + 1]] <-
    evaluate_config("block", fo$fold_ids, buf_used, bk * 1000, fs)
}

# Random CV reference (no blocks, no buffer)
set.seed(CV_CONFIG$seed)
rnd_ids <- sample(rep(seq_len(CV_CONFIG$n_folds), length.out = n))
rnd_res <- evaluate_config("random", rnd_ids, 0, NA_real_, CV_CONFIG$seed)

all_res <- c(buf_res, blk_res, list(rnd_res))
metrics <- bind_rows(lapply(all_res, `[[`, "metrics"))
nn_all  <- bind_rows(lapply(all_res, `[[`, "nn"))

summary_tab <- metrics |>
  group_by(sweep, buffer_km, block_km) |>
  summarise(n_runs = n(),
            n_train_mean = mean(n_train_mean),
            pct_train_removed = mean(pct_train_removed),
            nn_median_km = mean(nn_median_km),
            accuracy = mean(mean_accuracy), accuracy_sd_runs = sd(mean_accuracy),
            ci95_lower = mean(ci95_lower), ci95_upper = mean(ci95_upper),
            kappa = mean(mean_kappa),
            mining_f1 = mean(mean_mining_f1), mining_f1_sd_runs = sd(mean_mining_f1),
            log_loss = mean(mean_log_loss), .groups = "drop") |>
  arrange(sweep, buffer_km, block_km)

ref <- summary_tab |> filter(sweep == "buffer", abs(buffer_km - buf_used / 1000) < 1e-6)
summary_tab <- summary_tab |>
  mutate(d_accuracy_vs_used  = accuracy  - ref$accuracy,
         d_mining_f1_vs_used = mining_f1 - ref$mining_f1)

write_csv(metrics,     file.path(out_dir, "tables", "buffer_sensitivity.csv"))
write_csv(summary_tab, file.path(out_dir, "tables", "buffer_sensitivity_summary.csv"))
cat("\n")
print(as.data.frame(summary_tab |> select(sweep, buffer_km, block_km, n_train_mean, nn_median_km,
                                          accuracy, accuracy_sd_runs, mining_f1, log_loss,
                                          d_accuracy_vs_used) |>
                      mutate(across(where(is.numeric), ~ round(., 3)))), row.names = FALSE)

# ================================================================
# [4] PREDICTION-TO-TRAINING DISTANCES
# ================================================================
cat("\n[4/4] Prediction-location distances\n")
# Prediction locations: random points in Ghana (or the training hull)
pred_domain <- tryCatch({
  g <- st_read(PATHS$ghana_districts, quiet = TRUE)
  st_union(st_transform(g, GHANA_CRS))
}, error = function(e) st_convex_hull(st_union(data_sf)))
set.seed(CV_CONFIG$seed)
pp <- st_sample(pred_domain, 5000)
pred_nn <- as.numeric(st_distance(pp, data_sf[st_nearest_feature(pp, data_sf), ],
                                  by_element = TRUE))
# ... and inside the mining polygons, where the areas are reported
pred_nn_poly <- tryCatch({
  pol <- st_union(st_make_valid(st_transform(st_read(PATHS$mining_polygons, quiet = TRUE),
                                             GHANA_CRS)))
  pp2 <- st_sample(pol, 5000)
  as.numeric(st_distance(pp2, data_sf[st_nearest_feature(pp2, data_sf), ],
                         by_element = TRUE))
}, error = function(e) NULL)

nn_ref <- bind_rows(
  data.frame(sweep = "prediction", label = "Prediction: Ghana", nn_m = pred_nn),
  if (!is.null(pred_nn_poly))
    data.frame(sweep = "prediction", label = "Prediction: mining polygons",
               nn_m = pred_nn_poly))
cat(sprintf("  median NN distance, Ghana-wide prediction points: %.1f km\n",
            median(pred_nn) / 1000))
if (!is.null(pred_nn_poly))
  cat(sprintf("  median NN distance, mining-polygon prediction points: %.1f km\n",
              median(pred_nn_poly) / 1000))

nn_plot <- bind_rows(
  nn_all |> filter(sweep == "buffer") |>
    mutate(label = sprintf("CV buffer %.1f km", buffer_km)) |> select(sweep, label, nn_m),
  nn_all |> filter(sweep == "random") |>
    mutate(label = "Random CV") |> select(sweep, label, nn_m),
  nn_ref)
write_csv(nn_plot |> group_by(sweep, label) |>
            summarise(n = n(), median_km = median(nn_m) / 1000,
                      q25_km = quantile(nn_m, 0.25) / 1000,
                      q75_km = quantile(nn_m, 0.75) / 1000, .groups = "drop"),
          file.path(out_dir, "tables", "buffer_sensitivity_nn_distance.csv"))

# ================================================================
# FIGURE
# ================================================================
theme_s <- theme_bw(base_size = 13) +
  theme(plot.title = element_text(size = 20, face = "bold"),
        panel.grid.minor = element_blank(), legend.position = "bottom")

bm <- metrics |> filter(sweep == "buffer") |>
  select(buffer_km, rep, Accuracy = mean_accuracy, `Mining F1` = mean_mining_f1,
         Kappa = mean_kappa) |>
  pivot_longer(-c(buffer_km, rep), names_to = "metric")
bsum <- bm |> group_by(buffer_km, metric) |>
  summarise(m = mean(value), lo = min(value), hi = max(value), .groups = "drop")
rnd_acc <- summary_tab |> filter(sweep == "random")

p_a <- ggplot(bsum, aes(buffer_km + 0.1, m, colour = metric)) +
  geom_vline(xintercept = buf_used / 1000 + 0.1, linetype = "dashed", colour = "#D55E00") +
  geom_hline(yintercept = rnd_acc$accuracy, linetype = "dotted", colour = "grey40") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.04) +
  geom_line() + geom_point(size = 2) +
  scale_x_log10(breaks = c(0, 1, 2, 5, 10, 20, 50) + 0.1,
                labels = c("0", "1", "2", "5", "10", "20", "50")) +
  scale_colour_manual(values = c(Accuracy = "#0072B2", `Mining F1` = "#D55E00",
                                 Kappa = "#009E73"), name = NULL) +
  labs(title = "(a)", x = "Buffer distance (km)", y = "Spatial-CV score") + theme_s

p_b <- ggplot(summary_tab |> filter(sweep == "block"),
              aes(block_km, accuracy)) +
  geom_vline(xintercept = block_used / 1000, linetype = "dashed", colour = "#D55E00") +
  geom_hline(yintercept = rnd_acc$accuracy, linetype = "dotted", colour = "grey40") +
  geom_point(data = metrics |> filter(sweep == "block"),
             aes(block_km, mean_accuracy), colour = "grey60", size = 1.5) +
  geom_line(colour = "#0072B2") + geom_point(colour = "#0072B2", size = 2.5) +
  labs(title = "(b)", x = "Block size (km)", y = "Accuracy") + theme_s

# Panel (c) shows a readable subset; the CSV above keeps every buffer.
show_buf <- sort(unique(c(0, buf_used / 1000, if (is.finite(exact_used)) exact_used / 1000, 20)))
nn_show <- bind_rows(
  nn_all |> filter(sweep == "buffer",
                   sapply(buffer_km, function(b) any(abs(b - show_buf) < 1e-3))) |>
    mutate(label = ifelse(abs(buffer_km - buf_used / 1000) < 1e-3,
                          sprintf("CV buffer %.1f km (used)", buffer_km),
                          sprintf("CV buffer %.1f km", buffer_km)),
           order = buffer_km) |> select(sweep, label, nn_m, order),
  nn_all |> filter(sweep == "random") |>
    mutate(label = "Random CV", order = -1) |> select(sweep, label, nn_m, order),
  nn_ref |> mutate(order = ifelse(grepl("polygons", label), 1001, 1000)))
lab_order <- nn_show |> distinct(label, order) |> arrange(order) |> pull(label)
nn_show$label <- factor(nn_show$label, levels = lab_order)
cv_lab <- grep("^CV buffer", lab_order, value = TRUE)
pal_c  <- c(setNames(grDevices::colorRampPalette(c("#9ECAE1", "#08306B"))(length(cv_lab)), cv_lab),
            "Random CV" = "grey45",
            "Prediction: Ghana" = "#D55E00",
            "Prediction: mining polygons" = "#CC79A7")[lab_order]

p_c <- ggplot(nn_show, aes(nn_m / 1000, colour = label, linetype = sweep)) +
  stat_ecdf(linewidth = 0.8) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 2, 5, 10, 20, 50)) +
  scale_linetype_manual(values = c(buffer = "solid", random = "dotted",
                                   prediction = "longdash"), guide = "none") +
  coord_cartesian(xlim = c(0, 60)) +
  labs(title = "(c)", x = "Nearest training point (km)", y = "ECDF", colour = NULL) +
  scale_colour_manual(values = pal_c) +
  theme_s + guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                                         override.aes = list(linewidth = 1.2)))

fig <- (p_a | p_b) / p_c + plot_layout(heights = c(1, 1))
ggsave(file.path(plot_dir, "buffer_sensitivity.png"), fig,
       width = 14, height = 11, dpi = 400, bg = "white")
cat(sprintf("\n  -> %s\n", file.path(plot_dir, "buffer_sensitivity.png")))

# ================================================================
# VERDICT
# ================================================================
b_rng <- summary_tab |> filter(sweep == "buffer")
cat("\n---- Summary ----\n")
cat(sprintf("  accuracy over buffers 0-%g km: %.3f to %.3f (RF-seed SD at used buffer %.3f)\n",
            max(BUFFERS_KM), min(b_rng$accuracy), max(b_rng$accuracy), ref$accuracy_sd_runs))
cat(sprintf("  mining F1 over buffers: %.3f to %.3f\n",
            min(b_rng$mining_f1), max(b_rng$mining_f1)))
cat(sprintf("  random CV - used config: %+.3f accuracy\n", rnd_acc$accuracy - ref$accuracy))
cat("  A flat curve means the estimate does not hinge on the variogram range;\n")
cat("  a steady drop means the buffer matters and the larger-buffer value\n")
cat("  should be reported alongside (or instead of) the variogram one.\n")