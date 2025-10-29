# ==============================================================================
# PURPOSE:
# This script performs the *full-year* analysis to test the hypothesis
# $\gamma \approx \alpha - 1$.
#
# It calculates one $\gamma$ value (from the raw trade signs) and one $\alpha$
# value (from the simulated metaorders) for *each stock*, using all
# data from the entire year.
#
# It runs this analysis twice in parallel:
# 1. Using the "Original MLB" model to estimate $\alpha$.
# 2. Using the "Adaptive MLB" model to estimate $\alpha$.
#
# The script concludes by generating plots that compare the empirical $\gamma$
# values against the ($\alpha - 1$) values from both models.
#
# PREREQUISITES:
# 1. Must be run AFTER `06_LMF_Validation_load_and_prep_data.R` (to create the
#    "JSE_Top40_workspace.RData" file).
# 2. This script MUST `source("05_LMF_Validation_analysis_functions.R")` at the top
#    to load all the required helper functions.
# * Packages: 'data.table', 'parallel', 'poweRlaw', 'ggplot2', 'cowplot'.
#
# OUTPUT:
# * Several intermediate data.tables (e.g., 'gamma_dt_fullyear',
#   'meta_features_orig', 'results_orig', 'results_adapt').
# * Final plots (saved to the plots pane) comparing $\gamma$ vs. $\alpha - 1$
#   for both the Original and Adaptive models on the full-year data.
# ==============================================================================

source("05_LMF_Validation_analysis_functions.R")

# ==============================================================================
# ESTIMATING GAMMAS
# ==============================================================================
library(parallel)
library(pbapply) # Added for progress bar function
library(ggplot2)
library(poweRlaw)
library(cowplot)

# --- Prepare per-stock ordered sign sequences ---
stock_signs_dt <- DT[, .(signs = list(Sign[order(timestamp)])), by = Stock]

# --- parameters ---
min_signs_for_gamma <- 10000 
lag_min <- 5
lag_max_frac <- 0.2
nbins <- 18
low_freq_frac <- 0.002
ncores <- max(1, detectCores() - 4) 

# --- Run parallel estimation ---
gamma_dt_fullyear <- run_parallel_gamma_estimation(
  grouped_signs_dt = stock_signs_dt,
  grouping_cols = "Stock",
  min_signs_for_gamma = min_signs_for_gamma,
  lag_min = lag_min,
  lag_max_frac = lag_max_frac,
  nbins = nbins,
  low_freq_frac = low_freq_frac,
  ncores_to_use = ncores
)

# ==============================================================================
# Generating realistic metaorders: (Original MLB) 
# ==============================================================================

# --- MLB parameters ---
N     <- 100   # number of traders (choose as needed)
theta <- 0.5  # sampler parameter

meta_features_orig <- generate_metaorders_orig(
  DT = DT,
  N = N,
  theta = theta,
  f_min = 1,
  seed = 123
)

# ==============================================================================
# Estimating alpha (on original MLB assignments) on full year
# =============================================================================
results_orig <- run_parallel_alpha_estimation(
  meta_features_dt = meta_features_orig,
  grouping_cols = "Stock",
  min_sample_size = 50,
  ncores_to_use = ncores,
  xmin_filter = 2
)

# Show the estimates
print(results_orig)

# Plot histogram of results
ggplot(results_orig, aes(x = alpha)) +
  geom_histogram(aes(y = ..density..), bins =42, color = "blue", fill="lightblue", boundary = 0) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "red", size = 1) +
  scale_x_continuous(limits = c(0, NA)) + 
  labs( title = "(a) Original MLB",
        x = expression(alpha),
        y = expression(P(alpha))) +
  theme_minimal()

# ==============================================================================
# Generating realistic metaorders: (Adaptive blocking MLB)
# ==============================================================================

meta_features_adapt <- generate_metaorders_adapt(
  DT = DT,
  N = N,
  theta = theta,
  base_block = 100,
  f_min = 1,
  seed = 123
)

# ==============================================================================
# Estimating alpha (on adaptive MLB assignments) on full year
# ==============================================================================

results_adapt <- run_parallel_alpha_estimation(
  meta_features_dt = meta_features_adapt,
  grouping_cols = "Stock",
  min_sample_size = 50,
  ncores_to_use = ncores,
  xmin_filter = 2
)

# Show the estimates
print(results_adapt)

# Plot histogram of results
ggplot(results_adapt, aes(x = alpha)) +
  geom_histogram(aes(y = ..density..), bins =42, color = "blue", fill="lightblue", boundary = 0) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "red", size = 1) +
  scale_x_continuous(limits = c(0, NA)) + 
  labs(title = "(b) Adaptive MLB",
       x = expression(alpha),
       y = expression(P(alpha))) +
  theme_minimal()

# ==============================================================================
# Plotting gamma=alpha-1 (on original MLB assignments) on full year 
# ==============================================================================

plot_gamma_alpha_comparison(
  alpha_dt = results_orig,
  gamma_dt = gamma_dt_fullyear,
  merge_cols = "Stock",
  coord_limits = c(0, 4.5, 0, 3), # c(xmin, xmax, ymin, ymax)
  bin_width = 0.25
)

# ==============================================================================
# Plotting gamma=alpha-1 (on adaptive MLB assignments) on full year 
# ==============================================================================

plot_gamma_alpha_comparison(
  alpha_dt = results_adapt,
  gamma_dt = gamma_dt_fullyear,
  merge_cols = "Stock",
  coord_limits = c(0, 4.5, 0, 3), # c(xmin, xmax, ymin, ymax)
  bin_width = 0.25
)
