################################################################################
# Reproducible validation script for VI-based LAI Bayesian models (brms)
#
# TO BE COMPLETED after paper acceptance
# Supplementary material to paper:
# Bascietto M., et al "Canopy Closure breakpoints and Vegetation Index
# Saturation for Leaf Area Index Estimation across Stand Age and Phenology
# in Poplar Plantations"
#
# Purpose
# -------
# This script validates fitted Bayesian models (brms) for LAI prediction using
# independent LAI observations and corresponding vegetation indices (NDVI, EVI).
#
# Workflow summary
# ----------------
# 0) Load validation data (vld_l) and fitted models + training data (index_l)
# 1) Scale validation predictors using training-set scaling parameters
# 2) Generate population-level posterior predictions (no random effects)
# 3) Unscale predicted LAI back to original units; bind predictions to validation
# 4) Aggregate to weekly means (observed and predicted)
# 5) Compute validation metrics globally and by groups; build plots
#
# Notes on modelling conventions
# ------------------------------
# - Fitted models were trained on standardized predictors (Age, index_value) and
#   on a standardized LAI response; therefore both scaling and unscaling
#   steps are essential for interpretable validation.
# - Predictions use brms::posterior_epred(..., re_formula = NA) to obtain
#   population-level expectations (fixed effects only), excluding group-level
#   random effects. This is appropriate when validating transferability across
#   plots or when random-effect levels are not meant to be "learned" for new data.
################################################################################

# ------------------------------------------------------------------------------
# Packages
# ------------------------------------------------------------------------------
library(data.table)   # fast data manipulation (by-reference semantics)
library(brms)         # Bayesian regression models using Stan
library(ggplot2)      # plotting
library(patchwork)    # plot composition (layouting multiple ggplots)


# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
# Centralize column names and directories to avoid hard-coding
# throughout script.
cfg_l <- list(
  rsfd_data_dir = "data/",
  models_dir    = "models/",

  # Identifier for experimental unit (plot)
  plot_id = "plotID",

  # Columns that were standardized during model fitting.
  # IMPORTANT: LAI_r (observed LAI in validation) is in real units and must NOT
  # be scaled here; predicted LAI will be unscaled later.
  scaled_cols  = c("Age", "index_value"),

  # Observed LAI column in validation data (real, unscaled)
  lai_obs = "LAI_r",

  # Prediction summary columns produced by pred_pop_age()
  # (posterior mean + 95% credible interval)
  lai_cols     = c("LAI_mean", "LAI_low", "LAI_high")
)


# ==============================================================================
# %% Import data & fitted objects
# ==============================================================================

# ------------------------------------------------------------------------------
# Load validation dataset (independent data)
# ------------------------------------------------------------------------------
# Expected structure:
# - vld_l: named list with elements "NDVI" and "EVI"
#   Each element is a data.table including at least:
#   * plotID, week, phenoMnth, Age, index_value, LAI_r
#
# The dataset must be consistent with the predictors used in the fitted models.
vld_l <- readRDS(paste0(cfg_l$rsfd_data_dir, "validation_LAI.rds"))

# ------------------------------------------------------------------------------
# Load fitted models and training dataset
# ------------------------------------------------------------------------------
# Expected objects in "model.rds" (loaded into global environment):
# - fit_evi_l: list(AML = <brmsfit>, PAS = <brmsfit>)
# - fit_ndvi_l: list(AML = <brmsfit>, PAS = <brmsfit>)
# - index_l: list(NDVI = <training dt>, EVI = <training dt>)
# - model_frms_l, model_priors_l: (not needed for validation here)
#
# NOTE: list2env() intentionally places objects in .GlobalEnv to keep
# this script "single-file runnable" for supplement materials.
# In a package/project context,
# prefer explicit assignment to avoid side effects.
paste0(cfg_l$models_dir, "AML_PAS.rds") |>
  readRDS() |>
  list2env(envir = .GlobalEnv)

# Remove objects not required for this validation script
rm(model_frms_l, model_priors_l)


# ==============================================================================
# %% Helper functions
# ==============================================================================

# ------------------------------------------------------------------------------
# unscale_maker()
# ------------------------------------------------------------------------------
#' Create an unscaling function from training mean and sd
#'
#' @param sigma Numeric. Standard deviation used for scaling (training set).
#' @param mu Numeric. Mean used for centering (training set).
#' @return A function f(v) that maps standardized values back to original units:
#'         f(v) = v * sigma + mu
#'
#' @details
#' Scaling convention assumed is: (x - mu) / sigma.
#' Inverse transform is: x = z * sigma + mu.
unscale_maker <- function(sigma, mu) {
  function(v) {
    v * sigma + mu
  }
}

# ------------------------------------------------------------------------------
# unscale_lai
# ------------------------------------------------------------------------------
# Create unscaling function for LAI based on how LAI was stored in
# training data.
# Here we extract the attributes created by scale() in R:
# - "scaled:center" and "scaled:scale"
#
# IMPORTANT:
# This assumes the training dataset object index_l$NDVI$LAI exists and was
# produced via scale(LAI), so it carries the attributes used below.
unscale_lai <- with(
  index_l$NDVI,
  unscale_maker(
    attr(LAI, "scaled:scale"),
    attr(LAI, "scaled:center")
  )
)

# ------------------------------------------------------------------------------
# pred_pop_age()
# ------------------------------------------------------------------------------
#' Posterior expected predictions for new data (population-level)
#'
#' @param fit A brmsfit object.
#' @param idx_label Character. Label for the vegetation index (e.g., "NDVI", "EVI").
#' @param age_grd_dt data.table (or data.frame). Newdata used for prediction.
#' @return data.table with columns:
#'   - Index: idx_label
#'   - Age: newdata$Age (as provided; typically standardized here)
#'   - LAI_mean: posterior mean of E[LAI | predictors] for each row of newdata
#'   - LAI_low, LAI_high: 2.5% and 97.5% posterior quantiles (95% CrI)
#'
#' @details
#' Uses brms::posterior_epred(), which returns draws of the conditional mean
#' (i.e., expected value, not including observation-level noise).
#' Setting re_formula = NA excludes group-level (random) effects.
pred_pop_age <- function(fit, idx_label, age_grd_dt) {

  ep <- brms::posterior_epred(
    fit,
    newdata    = age_grd_dt,
    re_formula = NA
  )

  data.table::data.table(
    Index    = idx_label,
    Age      = age_grd_dt$Age,
    LAI_mean = colMeans(ep),
    LAI_low  = apply(ep, 2, quantile, 0.025),
    LAI_high = apply(ep, 2, quantile, 0.975)
  )
}

# ------------------------------------------------------------------------------
# val_metrics_dt()
# ------------------------------------------------------------------------------
#' Compute validation metrics for predicted LAI vs observed LAI
#'
#' @param dt data.table/data.frame containing observed and predicted LAI columns.
#' @param obs Character. Column name of observed LAI (real units).
#' @param pred Character. Column name of predicted LAI (real units).
#' @param low Character. Column name of lower credible interval (same units).
#' @param high Character. Column name of upper credible interval (same units).
#' @return A list with:
#'   - metrics: data.table with scalar summary metrics
#'   - dt: data.table copy of input, augmented with diagnostic columns
#'
#' @details
#' Metrics:
#' - RMSE: sqrt(mean squared error)
#' - MAE: mean absolute error
#' - Bias: mean signed error (pred - obs)
#' - MedAE: median absolute error (robust to outliers)
#' - R2_val: out-of-sample R^2 using SSE/SST formulation
# - Coverage95: empirical coverage of the nominal 95% credible interval
# - Mean_IW / Median_IW: mean/median interval width (IW = high - low)
#   (a proxy for predictive uncertainty; smaller is sharper, but must be paired
#   with adequate coverage).
val_metrics_dt <- function(
  dt,
  obs  = cfg_l$lai_obs,
  pred = cfg_l$lai_cols[1],
  low  = cfg_l$lai_cols[2],
  high = cfg_l$lai_cols[3]
) {

  dt <- dt |>
    data.table::as.data.table() |>
    data.table::copy()

  # Row-wise diagnostics:
  # err  : signed error
  # aerr : absolute error
  # serr : squared error
  # in_ci: indicator that observed LAI is within predictive 95% CrI
  # iw   : interval width (uncertainty band width)
  dt[, `:=`(
    err   = get(pred) - get(obs),
    aerr  = abs(get(pred) - get(obs)),
    serr  = (get(pred) - get(obs))^2,
    in_ci = (get(low) <= get(obs)) & (get(obs) <= get(high)),
    iw    = get(high) - get(low)
  )]

  # Global summary metrics
  out <- dt[, .(
    n = .N,
    RMSE = sqrt(mean(serr)),
    MAE  = mean(aerr),
    Bias = mean(err),
    MedAE = median(aerr),

    # Validation R^2 (SSE/SST). Interpretable as fraction of variance explained.
    # Can be negative if predictions are worse than using mean(obs).
    R2_val = {
      y    <- get(obs)
      yhat <- get(pred)
      1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
    },

    # Empirical coverage of nominal 95% credible interval
    Coverage95 = mean(in_ci),

    # Predictive interval width summaries
    Mean_IW   = mean(iw),
    Median_IW = median(iw)
  )]

  # Error distribution quantiles (useful for asymmetric bias / tail behavior)
  out[, `:=`(
    q05_err = quantile(dt$err, 0.05),
    q50_err = quantile(dt$err, 0.50),
    q95_err = quantile(dt$err, 0.95)
  )]

  list(metrics = out, dt = dt)
}


# ==============================================================================
# %% Step 1 — Scale validation predictors using training parameters
# ==============================================================================

# The models were trained with scaled predictors. To ensure comparability,
# we must apply the SAME mu and sigma (from training) to the validation set.
#
# We retrieve scaling parameters from training data attributes:
# attr(x, "scaled:center") and attr(x, "scaled:scale").
#
# IMPORTANT:
# This assumes training columns were created using base::scale() and still carry
# these attributes.
sigmas <- sapply(
  cfg_l$scaled_cols,
  \(x) attr(index_l$NDVI[[x]], "scaled:scale")
)

mus <- sapply(
  cfg_l$scaled_cols,
  \(x) attr(index_l$NDVI[[x]], "scaled:center")
)

stopifnot(
  length(sigmas) == length(cfg_l$scaled_cols),
  length(mus)    == length(cfg_l$scaled_cols)
)

# Apply scaling in-place to each validation data.table inside vld_l.
# We use data.table by-reference assignment for efficiency.
#
# For each column in scaled_cols:
#   z = (x - mu) / sigma
lapply(
  vld_l,
  function(vld_dt, cols) {
    vld_dt[
      ,
      (cols) := Map(
        \(col, sigma, mu) (col - mu) / sigma,
        .SD,
        sigma = sigmas,
        mu    = mus
      ),
      .SDcols = cols
    ]
  },
  cols = cfg_l$scaled_cols
) |>
  invisible()


# ==============================================================================
# %% Step 2 — Posterior predictions on validation grid (population-level)
# ==============================================================================

# Generate posterior expected predictions for the PAS models (selected).
# Output is in the same scale as the brms response used at fit time (often scaled).
pred_l <- list(
  ndvi = pred_pop_age(fit_ndvi_l$PAS, "NDVI", vld_l$NDVI),
  evi  = pred_pop_age(fit_evi_l$PAS,  "EVI",  vld_l$EVI)
)


# ==============================================================================
# %% Step 3 — Unscale predicted LAI, bind to validation, aggregate weekly
# ==============================================================================

# Unscale predicted LAI summaries back to real LAI units.
# This is necessary for comparison with LAI_r.
lapply(
  pred_l,
  \(dt, cols) {
    dt[, (cols) := lapply(.SD, unscale_lai), .SDcols = cols]
  },
  cols = cfg_l$lai_cols
) |>
  invisible()

# Bind predicted columns (LAI_mean/low/high) to the validation grid.
# Map() pairs NDVI validation with NDVI predictions, and EVI with EVI.
vld_fit_l <- Map(
  \(x, y, cols) cbind(x, y[, ..cols]),
  vld_l,
  pred_l,
  MoreArgs = list(cols = cfg_l$lai_cols)
)

# Aggregate to weekly means per plot and phenological month.
# This reduces day-level noise and aligns with the intended temporal support
# of the validation (weekly LAI summaries).
vld_fit_l <- lapply(
  vld_fit_l,
  \(dt, cols, id) {
    dt[
      ,
      lapply(.SD, mean),
      by = c(id, "week", "phenoMnth"),
      .SDcols = cols
    ]
  },
  cols = c(cfg_l$lai_obs, cfg_l$lai_cols),
  id   = cfg_l$plot_id
) |>
  invisible()


# ==============================================================================
# %% Step 4 — Compute validation metrics (overall + grouped)
# ==============================================================================

# Overall (global) metrics
res_evi  <- val_metrics_dt(vld_fit_l$EVI)
res_ndvi <- val_metrics_dt(vld_fit_l$NDVI)

# Grouped metrics:
# by phenological month (seasonal/phenological dependence)
metrics_by_mnth_evi <- res_evi$dt[, .(
  n = .N,
  RMSE = sqrt(mean(serr)),
  MAE  = mean(aerr),
  Bias = mean(err),
  Coverage95 = mean(in_ci),
  Mean_IW = mean(iw)
), by = .(phenoMnth)]

metrics_by_mnth_ndvi <- res_ndvi$dt[, .(
  n = .N,
  RMSE = sqrt(mean(serr)),
  MAE  = mean(aerr),
  Bias = mean(err),
  Coverage95 = mean(in_ci),
  Mean_IW = mean(iw)
), by = .(phenoMnth)]

# Combined overall table for reporting (e.g., manuscript Table 5)
metrics_all <- rbindlist(
  list(
    cbind(Index = "EVI",  res_evi$metrics),
    cbind(Index = "NDVI", res_ndvi$metrics)
  ),
  use.names = TRUE,
  fill      = TRUE
)

# Coverage by phenological month (diagnostic of uncertainty calibration)
cov_ndvi_dt <- res_ndvi$dt[, .(Coverage = mean(in_ci)), by = phenoMnth]
cov_evi_dt  <- res_evi$dt[,  .(Coverage = mean(in_ci)), by = phenoMnth]


# ==============================================================================
# %% Final plot — Agreement scatter + coverage by month
# ==============================================================================

# Panel A: observed vs predicted (NDVI)
# - color by phenoMnth (seasonal partition)
# - dashed 1:1 line for visual assessment
# - coord_equal + fixed limits to standardize visual comparability
plots_l <- list(
  scatter_ndvi = ggplot(
    res_ndvi$dt,
    aes(x = LAI_r, y = LAI_mean, color = phenoMnth)
  ) +
    geom_point(size = 3, show.legend = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_y_continuous(limits = c(0, 1.5)) +
    scale_x_continuous(limits = c(0, 1.5)) +
    coord_equal() +
    labs(
      title = "A) NDVI Model–data agreement",
      x = "Observed LAI",
      y = "Predicted LAI"
    ) +
    theme_light()
)

# Panels C, B, D are built by reusing structure:
# - scatter_evi: same geometry but using EVI residual table
# - cv_*: coverage by phenological month with reference line at 0.95
plots_l <- within(
  plots_l, {

    # Panel C: observed vs predicted (EVI)
    scatter_evi <- scatter_ndvi +
      res_evi$dt +
      labs(title = "C) EVI Model–data agreement")

    # Panel B: coverage by month (NDVI)
    cv_ndvi <- ggplot(
      cov_ndvi_dt,
      aes(phenoMnth, Coverage, color = phenoMnth)
    ) +
      geom_point(size = 3, show.legend = FALSE) +
      geom_hline(yintercept = 0.95, linetype = "dotted") +
      scale_y_continuous(limits = c(0.5, 1.0)) +
      labs(
        title = "B) NDVI Coverage (95%)",
        x = "Phenological month",
        y = "Empirical coverage"
      ) +
      theme_light()

    # Panel D: coverage by month (EVI)
    cv_evi <- cv_ndvi +
      cov_evi_dt +
      labs(title = "D) EVI Coverage (95%)")
  }
)

# Assemble final multi-panel figure (2x2)
with(
  plots_l,
  (scatter_ndvi / cv_ndvi) | (scatter_evi / cv_evi)
)



# > sessionInfo()
# R version 4.5.2 (2025-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS Tahoe 26.2

# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

# locale:
# [1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8

# time zone: Europe/Rome
# tzcode source: internal

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] patchwork_1.3.2   ggplot2_4.0.1     brms_2.23.0       Rcpp_1.1.0       
# [5] data.table_1.18.0

# loaded via a namespace (and not attached):
#  [1] Matrix_1.7-4          bayesplot_1.15.0      gtable_0.3.6         
#  [4] jsonlite_2.0.0        dplyr_1.1.4           compiler_4.5.2       
#  [7] tidyselect_1.2.1      stringr_1.6.0         parallel_4.5.2       
# [10] scales_1.4.0          lattice_0.22-7        coda_0.19-4.1        
# [13] R6_2.6.1              Brobdingnag_1.2-9     generics_0.1.4       
# [16] distributional_0.6.0  backports_1.5.0       checkmate_2.3.3      
# [19] tibble_3.3.0          pillar_1.11.1         RColorBrewer_1.1-3   
# [22] posterior_1.6.1       rlang_1.1.6           stringi_1.8.7        
# [25] S7_0.2.1              RcppParallel_5.1.11-1 cli_3.6.5            
# [28] withr_3.0.2           magrittr_2.0.4        rstantools_2.6.0     
# [31] grid_4.5.2            mvtnorm_1.3-3         lifecycle_1.0.4      
# [34] nlme_3.1-168          vctrs_0.6.5           tensorA_0.36.2.1     
# [37] glue_1.8.0            farver_2.1.2          bridgesampling_1.2-1 
# [40] abind_1.4-8           matrixStats_1.5.0     tools_4.5.2          
# [43] loo_2.9.0             pkgconfig_2.0.3      
