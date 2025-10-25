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
# --- Prepare per-stock ordered sign sequences to avoid repeated subsetting in workers ---
stock_signs_dt <- DT[, .(signs = list(Sign[order(timestamp)])), by = Stock]

# --- parameters ---
min_signs_for_gamma <- 200
lag_min <- 5
lag_max_frac <- 0.2
nbins <- 18
low_freq_frac <- 0.12

# --- cluster setup ---
ncores <- max(1, detectCores() - 4) 
cl <- makeCluster(ncores)

# ensure needed packages are loaded on workers
clusterEvalQ(cl, { library(data.table); NULL })

# Export data, functions and parameter variables to workers
clusterExport(cl, varlist = c(
  "stock_signs_dt",
  "estimate_gamma_acf",
  "estimate_gamma_psd",
  "min_signs_for_gamma",
  "lag_min",
  "lag_max_frac",
  "nbins",
  "low_freq_frac"
), envir = environment())

# --- run computations in parallel (load-balanced) ---
idx_vec <- seq_len(nrow(stock_signs_dt))
res_list <- parLapplyLB(cl, idx_vec, function(idx) {
  # worker-side: use only exported objects and functions
  sname <- stock_signs_dt$Stock[idx]
  signs <- stock_signs_dt$signs[[idx]]
  if (length(signs) < min_signs_for_gamma) {
    gamma_acf <- NA_real_; gamma_psd <- NA_real_
  } else {
    gamma_acf <- tryCatch(estimate_gamma_acf(signs, lag_min = lag_min, lag_max_frac = lag_max_frac, nbins = nbins),
                          error = function(e) NA_real_)
    gamma_psd <- tryCatch(estimate_gamma_psd(signs, low_freq_frac = low_freq_frac), error = function(e) NA_real_)
  }
  list(Stock = sname, gamma_acf = gamma_acf, gamma_psd = gamma_psd, n_orders = length(signs))
})

# --- clean up ---
stopCluster(cl)

# Combine results into a data.table with the same structure as before
gamma_dt_fullyear <- rbindlist(lapply(res_list, as.data.table))

# ==============================================================================
# Generating realistic metaorders: (Original MLB) 
# ==============================================================================

# Compute daily stats
daily_stats_orig <- DT[, .(
  VD   = sum(volume, na.rm = TRUE), 
  p0   = price[1L], 
  pmax = max(price, na.rm = TRUE),
  pmin = min(price, na.rm = TRUE)
), by = .(Stock, trade_date)]

# Daily range volatility proxy
daily_stats_orig[, sigmaD := (pmax - pmin) / p0]

# Order and join back [not sure why this order line is here so I commented it out]
#setorder(daily_stats, Stock, trade_date) # had 'Stock' in this sort, but think it should be removed
#setorder(daily_stats, trade_date) # 'preserve true chronological order of trades' [from step3 of Alg 1]
DT_origMLB <- daily_stats_orig[DT, on = .(Stock, trade_date)]

# Filter small/odd days (not sure why this is here so commented it out)
#min_trades_per_day <- 50L
#DT[, trades_in_day := .N, by = .(Stock, trade_date)]
#DT <- DT[trades_in_day >= min_trades_per_day]

daily_stats_orig <- daily_stats_orig[DT_origMLB[, .(Stock, trade_date)] |> unique(), on=.(Stock, trade_date)]

# Save workspace for future sessions
#save(list = c("DT", "daily_stats"), file = workspace_file)
#message("Workspace saved to ", workspace_file)

# Next, define the 'mapping function' [view key functions code]

# First, preserve true chronological order of trades by sorting by tradestamp
setorder(DT_origMLB, timestamp)
# Next, define the 'sampling function' used within the mapping function [view key functions code]

# call the mapping function to randomly assign trades to traders
# parameters
N     <- 100   # number of traders (choose as needed)
theta <- 0.5  # sampler parameter (higher thetas means more fragmentation between traders.)
# Set the seed ONCE before the group assignment for reproducibility
set.seed(123) 

# Call the mapping function FOR EACH STOCK-DAY GROUP
DT_origMLB[, trader_id := mapping_function_orig(
  N         = N,
  F_sampler = sample_powerlaw,
  theta     = theta,
  n_orders  = .N,     # .N is the number of rows (orders) in the current group
  seed      = NULL,   # Use NULL so the RNG continues and is not reset for every group
  f_min     = 1
)$assignments, 
by = .(Stock, trade_date)] # This performs the assignment per stock and per day

# Step 4: Here we need to now sort trades by traders and timestamp
setorder(DT_origMLB, trader_id, timestamp)
# Step 5: Now we define a metaorder as a sequence of trades of same sign from same traders
DT_origMLB[, metaorder_id := rleid(trader_id, Sign)] 
# Step 6: 
meta_features_orig <- DT_origMLB[, .(
  # I have commented out the features we do not care for
  trader_id        = first(trader_id),
  Stock            = first(Stock),
  Sign             = first(Sign),
  start_time       = first(timestamp),
  end_time         = last(timestamp),
  #start_log_price  = log(first(price)),
  #end_log_price    = log(last(price)),
  n_child_orders   = .N,
  total_volume     = sum(volume)
), by = metaorder_id][n_child_orders > 1]

# Diagnostics check
pareto_pdf_scaled <- function(x, alpha, C) {
  x[x < 1] <- 1
  return(C * alpha / (x^(alpha + 1)))
}
# Generate the plot
ggplot(meta_features_orig, aes(x = n_child_orders)) +
  geom_histogram(bins = 50, fill = "lightblue", color = "black") +
  stat_function(
    fun = pareto_pdf_scaled,
    args = list(alpha = 1.15, C = 10000000),
    aes(color = "Pareto Fit (theta=1.15)"), 
    linetype = "dashed",
    linewidth = 1
  ) +
  scale_color_manual(
    name = "Distribution Fit",
    values = c("Pareto Fit (theta=1.15)" = "red")
  ) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(
    title = "(a) Original MLB", 
    x = "Number of Child Orders (Log Scale)",
    y = "Count (Log Scale)"
  ) +
  theme_minimal()+
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top")
  )

# ==============================================================================
# Estimating alpha (on original MLB assignments) on full year
# ==============================================================================
library(poweRlaw)
setDT(meta_features_orig)
stocks <- unique(meta_features_orig$Stock)

results_orig <- data.table(
  Stock = stocks,
  alpha = as.numeric(NA),
  xmin  = as.integer(NA),
  n_tail = as.integer(NA)
)

for (i in seq_along(stocks)) {
  st <- stocks[i]
  # extract integer metaorder lengths for this stock
  xs <- meta_features_orig[Stock == st, n_child_orders]
  # discard NA or non-positive values, keep integers
  xs <- xs[!is.na(xs) & xs >= 2] # here we choose the smallest metaorder size we wish to use, so that the alpha_estimation works nicely (eg: 2)
  # need at least a small sample to fit; otherwise leave NA
  if (length(xs) < 50) {   # conservative lower bound is at 50; adjust as you see fit
    results[i, `:=`(alpha = NA_real_, xmin = NA_integer_, n_tail = length(xs))]
    next
  }
  # fit discrete power law with poweRlaw
  fit_alpha <- tryCatch({
    m <- displ$new(xs)                 # discrete power-law object
    est <- estimate_xmin(m)            # estimate xmin and pars (Clauset method)
    # set estimated xmin and pars to object (estimate_xmin returns a list-like object)
    m$setXmin(est$xmin)
    m$setPars(est$pars)
    # record results
    xmin_val <- est$xmin
    alpha_val <- est$pars
    n_tail_val <- sum(xs >= xmin_val)
    list(alpha = as.numeric(alpha_val), xmin = as.integer(xmin_val), n_tail = as.integer(n_tail_val))
  }, error = function(e) {
    # on failure, return NA
    list(alpha = NA_real_, xmin = NA_integer_, n_tail = as.integer(length(xs)))
  })
  
  results_orig[i, `:=`(alpha = fit_alpha$alpha, xmin = fit_alpha$xmin, n_tail = fit_alpha$n_tail)]
}

# Show the estimates
print(results_orig)

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
# ===== STEP 2 of Algorithm 1 ===== #
# Compute daily stats
daily_stats_adapt <- DT[, .(
  VD   = sum(volume, na.rm = TRUE), 
  p0   = price[1L], 
  pmax = max(price, na.rm = TRUE),
  pmin = min(price, na.rm = TRUE)
), by = .(Stock, trade_date)]

# Daily range volatility proxy
daily_stats_adapt[, sigmaD := (pmax - pmin) / p0]

# Order and join back [not sure why this order line is here so I commented it out]
#setorder(daily_stats, Stock, trade_date) # had 'Stock' in this sort, but think it should be removed
#setorder(daily_stats, trade_date) # 'preserve true chronological order of trades' [from step3 of Alg 1]
DT_adaptMLB <- daily_stats_adapt[DT, on = .(Stock, trade_date)]

# Filter small/odd days (not sure why this is here so commented it out)
#min_trades_per_day <- 50L
#DT[, trades_in_day := .N, by = .(Stock, trade_date)]
#DT <- DT[trades_in_day >= min_trades_per_day]

daily_stats_adapt <- daily_stats_adapt[DT_adaptMLB[, .(Stock, trade_date)] |> unique(), on=.(Stock, trade_date)]

# Save workspace for future sessions
#save(list = c("DT", "daily_stats"), file = workspace_file)
#message("Workspace saved to ", workspace_file)

# call the mapping function to randomly assign trades to traders
# parameters
N     <- 100   # number of traders (choose as needed)
theta <- 0.5  # sampler parameter (higher thetas means more fragmentation between traders.)
set.seed(123) 

# Call the mapping function FOR EACH STOCK-DAY GROUP
DT_adaptMLB[, trader_id := mapping_function_adapt(
  N         = N,
  F_sampler = sample_powerlaw,
  theta     = theta,
  n_orders  = .N,     # .N is the number of rows (orders) in the current group
  seed      = NULL,   # Use NULL so the RNG continues and is not reset for every group
  adaptive = TRUE,    # enable adaptive block sizes
  base_block = 100,   # base block size (tune this: larger -> longer metaorders)
  f_min     = 1
)$assignments, 
by = .(Stock, trade_date)] # This performs the assignment per stock and per day

# Step 4: Here we need to now sort trades by traders and timestamp
setorder(DT_adaptMLB, trader_id, timestamp)
# Step 5: Now we define a metaorder as a sequence of trades of same sign from same traders
DT_adaptMLB[, metaorder_id := rleid(trader_id, Sign)] 
# Step 6: 
meta_features_adapt <- DT_adaptMLB[, .(
  # I have commented out the features we do not care for
  trader_id        = first(trader_id),
  Stock            = first(Stock),
  Sign             = first(Sign),
  start_time       = first(timestamp),
  end_time         = last(timestamp),
  #start_log_price  = log(first(price)),
  #end_log_price    = log(last(price)),
  n_child_orders   = .N,
  total_volume     = sum(volume)
), by = metaorder_id][n_child_orders > 1]

# Diagnostic checks
library(ggplot2)
pareto_pdf_scaled <- function(x, alpha, C) {
  x[x < 1] <- 1
  return(C * alpha / (x^(alpha + 1)))
}
# Generate the plot
ggplot(meta_features_adapt, aes(x = n_child_orders)) +
  geom_histogram(bins = 50, fill = "lightblue", color = "black") +
  stat_function(
    fun = pareto_pdf_scaled,
    args = list(alpha = 1.15, C = 10000000),
    aes(color = "Pareto Fit (theta=1.15)"), 
    linetype = "dashed",
    linewidth = 1
  ) +
  scale_color_manual(
    name = "Distribution Fit",
    values = c("Pareto Fit (theta=1.15)" = "red")
  ) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(
    title = "(b) Adaptive MLB",
    x = "Number of Child Orders (Log Scale)",
    y = "Count (Log Scale)"
  ) +
  theme_minimal()+
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top")
  )
# ==============================================================================
# Estimating alpha (on adaptive MLB assignments) on full year
# ==============================================================================
library(poweRlaw)
setDT(meta_features_adapt)
stocks <- unique(meta_features_adapt$Stock)

results_adapt <- data.table(
  Stock = stocks,
  alpha = as.numeric(NA),
  xmin  = as.integer(NA),
  n_tail = as.integer(NA)
)

for (i in seq_along(stocks)) {
  st <- stocks[i]
  # extract integer metaorder lengths for this stock
  xs <- meta_features_adapt[Stock == st, n_child_orders]
  # discard NA or non-positive values, keep integers
  xs <- xs[!is.na(xs) & xs >= 2] # here we choose the smallest metaorder size we wish to use, so that the alpha_estimation works nicely (eg: 2)
  # need at least a small sample to fit; otherwise leave NA
  if (length(xs) < 50) {   # conservative lower bound is at 50; adjust as you see fit
    results[i, `:=`(alpha = NA_real_, xmin = NA_integer_, n_tail = length(xs))]
    next
  }
  # fit discrete power law with poweRlaw
  fit_alpha <- tryCatch({
    m <- displ$new(xs)                 # discrete power-law object
    est <- estimate_xmin(m)            # estimate xmin and pars (Clauset method)
    # set estimated xmin and pars to object (estimate_xmin returns a list-like object)
    m$setXmin(est$xmin)
    m$setPars(est$pars)
    # record results
    xmin_val <- est$xmin
    alpha_val <- est$pars
    n_tail_val <- sum(xs >= xmin_val)
    list(alpha = as.numeric(alpha_val), xmin = as.integer(xmin_val), n_tail = as.integer(n_tail_val))
  }, error = function(e) {
    # on failure, return NA
    list(alpha = NA_real_, xmin = NA_integer_, n_tail = as.integer(length(xs)))
  })
  
  results_adapt[i, `:=`(alpha = fit_alpha$alpha, xmin = fit_alpha$xmin, n_tail = fit_alpha$n_tail)]
}

# Show the estimates
print(results_adapt)

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
# Join alpha (results) and gamma (gamma_dt) and prepare plotting objects used by your code
library(data.table)
library(ggplot2)

# 1) Merge into results_FINAL by Stock (keep only stocks present in results)
results_FINAL_orig <- merge(results_orig, gamma_dt_fullyear, by = "Stock", all.x = TRUE)

# ensure numeric and handle missing
results_FINAL_orig[, alpha := as.numeric(alpha)]
results_FINAL_orig[, gamma_acf := as.numeric(gamma_acf)]
results_FINAL_orig[, gamma_psd := as.numeric(gamma_psd)]

# compute alpha_minus1
results_FINAL_orig[, alpha_minus1 := alpha - 1]

# define bins of width 0.05 on alpha_minus1 (midpoint used as bin key)
bin_width <- 0.1
# compute bin index and bin midpoint
results_FINAL_orig[, alpha_bin := (floor(alpha_minus1 / bin_width) + 0.5) * bin_width]

# remove rows with NA alpha_minus1 so boxplots ignore them
res_plot <- results_FINAL_orig[!is.na(alpha_minus1)]

# 3) compute medians per bin for plotting the red median line/points
medians <- res_plot[, .(
  median_acf = median(gamma_acf, na.rm = TRUE),
  median_psd = median(gamma_psd, na.rm = TRUE)
), by = .(alpha_bin)]

# set x column expected by your plotting code
medians[, x := alpha_bin]

# 4) Now run your plotting code 
p_acf <- ggplot(res_plot, aes(x = alpha_minus1, y = gamma_acf)) +
  geom_boxplot(aes(group = alpha_bin), width = 0.05, outlier.shape = NA,
               fill = "gray80", color = "black") +
  geom_point(alpha = 0.3, color = "orange", size = 1.4) +
  geom_line(data = medians, aes(x = x, y = median_acf), color = "red", linewidth = 1) +
  geom_point(data = medians, aes(x = x, y = median_acf), shape = 25, fill = "red", size = 2) +
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.5) +
  labs(x = expression(alpha - 1),
       y = expression(gamma^{(a)}[unbiased]),
       title = "(a) ACF method") +
  coord_cartesian(xlim = c(0, 4.5), ylim = c(0, 3)) + 
  theme_minimal(base_size = 12)

p_psd <- ggplot(res_plot, aes(x = alpha_minus1, y = gamma_psd)) +
  geom_boxplot(aes(group = alpha_bin), width = 0.05, outlier.shape = NA,
               fill = "gray80", color = "black") +
  geom_point(alpha = 0.3, color = "orange", size = 1.4) +
  geom_line(data = medians, aes(x = x, y = median_psd), color = "red", linewidth = 1) +
  geom_point(data = medians, aes(x = x, y = median_psd), shape = 25, fill = "red", size = 2) +
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.5) +
  labs(x = expression(alpha - 1),
       y = expression(gamma^{(s)}[unbiased]),
       title = "(b) PSD method") +
  coord_cartesian(xlim = c(0, 4.5), ylim = c(0, 3)) +
  theme_minimal(base_size = 12)

# combine plots (same as your final line)
cowplot::plot_grid(p_acf, p_psd, ncol = 2, align = "hv")

# ==============================================================================
# Plotting gamma=alpha-1 (on adaptive MLB assignments) on full year 
# ==============================================================================
# Join alpha (results) and gamma (gamma_dt) and prepare plotting objects used by your code
library(data.table)
library(ggplot2)

# 1) Merge into results_FINAL by Stock (keep only stocks present in results)
results_FINAL_adaptive <- merge(results_adapt, gamma_dt_fullyear, by = "Stock", all.x = TRUE)

# ensure numeric and handle missing
results_FINAL_adaptive[, alpha := as.numeric(alpha)]
results_FINAL_adaptive[, gamma_acf := as.numeric(gamma_acf)]
results_FINAL_adaptive[, gamma_psd := as.numeric(gamma_psd)]

# compute alpha_minus1
results_FINAL_adaptive[, alpha_minus1 := alpha - 1]

# define bins of width 0.05 on alpha_minus1 (midpoint used as bin key)
bin_width <- 0.1
# compute bin index and bin midpoint
results_FINAL_adaptive[, alpha_bin := (floor(alpha_minus1 / bin_width) + 0.5) * bin_width]

# remove rows with NA alpha_minus1 so boxplots ignore them
res_plot <- results_FINAL_adaptive[!is.na(alpha_minus1)]

# 3) compute medians per bin for plotting the red median line/points
medians <- res_plot[, .(
  median_acf = median(gamma_acf, na.rm = TRUE),
  median_psd = median(gamma_psd, na.rm = TRUE)
), by = .(alpha_bin)]

# set x column expected by your plotting code
medians[, x := alpha_bin]

# 4) Now run your plotting code 
p_acf <- ggplot(res_plot, aes(x = alpha_minus1, y = gamma_acf)) +
  geom_boxplot(aes(group = alpha_bin), width = 0.05, outlier.shape = NA,
               fill = "gray80", color = "black") +
  geom_point(alpha = 0.3, color = "orange", size = 1.4) +
  geom_line(data = medians, aes(x = x, y = median_acf), color = "red", linewidth = 1) +
  geom_point(data = medians, aes(x = x, y = median_acf), shape = 25, fill = "red", size = 2) +
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.5) +
  labs(x = expression(alpha - 1),
       y = expression(gamma^{(a)}[unbiased]),
       title = "(a) ACF method") +
  coord_cartesian(xlim = c(0, 4.5), ylim = c(0, 3)) + 
  theme_minimal(base_size = 12)

p_psd <- ggplot(res_plot, aes(x = alpha_minus1, y = gamma_psd)) +
  geom_boxplot(aes(group = alpha_bin), width = 0.05, outlier.shape = NA,
               fill = "gray80", color = "black") +
  geom_point(alpha = 0.3, color = "orange", size = 1.4) +
  geom_line(data = medians, aes(x = x, y = median_psd), color = "red", linewidth = 1) +
  geom_point(data = medians, aes(x = x, y = median_psd), shape = 25, fill = "red", size = 2) +
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.5) +
  labs(x = expression(alpha - 1),
       y = expression(gamma^{(s)}[unbiased]),
       title = "(b) PSD method") +
  coord_cartesian(xlim = c(0, 4.5), ylim = c(0, 3)) +
  theme_minimal(base_size = 12)

# combine plots (same as your final line)
cowplot::plot_grid(p_acf, p_psd, ncol = 2, align = "hv")