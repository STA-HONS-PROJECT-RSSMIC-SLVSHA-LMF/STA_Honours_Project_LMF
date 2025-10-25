# ==============================================================================
# SCRIPT: 08_LMF_Validation_analysis_monthly
#
# PURPOSE:
# This script performs the *monthly* analysis to test the hypothesis
# $\gamma \approx \alpha - 1$. This approach is more granular than the
# full-year script, generating a separate data point for each stock
# in each month.
#
# It calculates one $\gamma$ value (from the raw trade signs) and one $\alpha$
# value (from the simulated metaorders) for *each stock-month combination*.
#
# It runs this analysis twice in parallel:
# 1. Using the "Original MLB" model to estimate $\alpha$.
# 2. Using the "Adaptive MLB" model to estimate $\alpha$.
#
# The script concludes by generating plots that compare the empirical $\gamma$
# values against the ($\alpha - 1$) values from both models, using all
# available stock-month data points.
#
# PREREQUISITES:
# 1. Must be run AFTER `06_LMF_Validation_load_and_prep_data.R` (to create the
#    "JSE_Top40_workspace.RData" file), as well as '07_LMF_Validation_analysis_full_year.R'.
# 2. This script MUST `source("05_LMF_Validation_analysis_functions.R")` at the top
#    to load all the required helper functions.
# * Packages: 'data.table', 'parallel', 'pbapply', 'poweRlaw', 'ggplot2',
#   'cowplot'.
#
# OUTPUT:
# * Several intermediate data.tables (e.g., 'gamma_dt_monthly',
#   'results_dt_monthly_orig', 'results_dt_monthly_adapt').
# * Final plots (saved to the plots pane) comparing $\gamma$ vs. $\alpha - 1$
#   for both the Original and Adaptive models on the monthly data.

source("05_LMF_Validation_analysis_functions.R")

# ==============================================================================
# ESTIMATING GAMMAS
# ==============================================================================
library(parallel)
library(pbapply) # <--- ADDED: Load the progress bar library

# --- Prepare per-stock AND per-month ordered sign sequences ---
stock_month_signs_dt <- DT[, .(signs = list(Sign[order(timestamp)])), by = .(Stock, month_year)]

# --- parameters (same as your original loop) ---
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
  "stock_month_signs_dt",
  "estimate_gamma_acf",
  "estimate_gamma_psd",
  "min_signs_for_gamma",
  "lag_min",
  "lag_max_frac",
  "nbins",
  "low_freq_frac"
), envir = environment())

# --- run computations in parallel (load-balanced) ---
idx_vec <- seq_len(nrow(stock_month_signs_dt)) 

# <--- MODIFIED: Replaced parLapplyLB with pblapply
# We provide the cluster 'cl', the vector to iterate 'X = idx_vec',
# and the function 'FUN = ...'
res_list <- pblapply(cl = cl, X = idx_vec, FUN = function(idx) {
  # (The function body is identical to your original)
  
  # <--- worker-side: use only exported objects and functions --->
  sname <- stock_month_signs_dt$Stock[idx]
  m_year <- stock_month_signs_dt$month_year[idx]
  signs <- stock_month_signs_dt$signs[[idx]]
  
  if (length(signs) < min_signs_for_gamma) {
    gamma_acf <- NA_real_; gamma_psd <- NA_real_
  } else {
    gamma_acf <- tryCatch(estimate_gamma_acf(signs, lag_min = lag_min, lag_max_frac = lag_max_frac, nbins = nbins),
                          error = function(e) NA_real_)
    gamma_psd <- tryCatch(estimate_gamma_psd(signs, low_freq_frac = low_freq_frac), error = function(e) NA_real_)
  }
  
  list(Stock = sname, month_year = m_year, gamma_acf = gamma_acf, gamma_psd = gamma_psd, n_orders = length(signs))
})

# --- clean up ---
stopCluster(cl)

# Combine results into a data.table
gamma_dt_monthly <- rbindlist(lapply(res_list, as.data.table))

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
# Estimating alpha (on original blocking MLB assignments) on month-month
# ==============================================================================
# --- Load required libraries ---
library(poweRlaw)
library(data.table)
library(parallel)
library(pbapply) # For the progress bar

# --- Ensure meta_features is a data.table ---
setDT(meta_features_orig)

# --- 1. Add month_year column using 'start_time' ---
# (This is the column we added in the previous step)
meta_features_orig[, month_year := format(start_time, "%Y-%m")]

# --- 2. Pre-group data by Stock and month_year (like we did for gamma) ---
# We also pre-filter the data here to match your original logic (>= 2)
stock_month_data <- meta_features_orig[, 
                                       .(xs_data = list(n_child_orders[!is.na(n_child_orders) & n_child_orders >= 2])), 
                                       by = .(Stock, month_year)
]

# --- 4. Set up cluster ---
ncores <- max(1, detectCores() - 4)
cl <- makeCluster(ncores)

# --- 5. Export necessary data/functions to workers ---
# Ensure workers have data.table and poweRlaw
clusterEvalQ(cl, {
  library(data.table)
  library(poweRlaw)
  NULL
})

# Export the pre-grouped data and the estimation function
clusterExport(cl, varlist = c("stock_month_data", "estimate_alpha_from_vec"), envir = environment())

# --- 6. Run in parallel with progress bar ---
# Iterate over the rows of our new grouped data.table
idx_vec <- seq_len(nrow(stock_month_data))

res_list <- pblapply(cl = cl, X = idx_vec, FUN = function(i) {
  # Get the pre-grouped data vector for this task
  xs <- stock_month_data$xs_data[[i]]
  
  # Run the estimation function
  # You can change the minimum sample size here if needed
  fit_res <- estimate_alpha_from_vec(xs, min_sample_size = 50)
  
  # Return a combined list with Stock/month info
  list(
    Stock = stock_month_data$Stock[i],
    month_year = stock_month_data$month_year[i],
    alpha = fit_res$alpha,
    xmin = fit_res$xmin,
    n_tail = fit_res$n_tail,
    n_obs = fit_res$n_obs
  )
})

# --- 7. Clean up cluster ---
stopCluster(cl)

# --- 8. Combine results ---
# This is your final data.table with alpha estimates per stock-month
results_dt_monthly_orig <- rbindlist(lapply(res_list, as.data.table))

# --- 9. Show the estimates and plot ---
print(results_dt_monthly_orig)

# This ggplot code will now plot the histogram of *all*
# the monthly alpha estimates (e.g., all 485 of them).
ggplot(results_dt_monthly_orig, aes(x = alpha)) +
  geom_histogram(aes(y = ..density..), bins = 42, color = "blue", fill="lightblue", boundary = 0) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "red", size = 1) +
  scale_x_continuous(limits = c(0, NA)) + 
  labs(title = "(a) Original MLB", # <--- Updated title
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
# Estimating alpha (on adaptive MLB assignments) on month-month
# ==============================================================================
# --- Load required libraries ---
library(poweRlaw)
library(data.table)
library(parallel)
library(pbapply) # For the progress bar

# --- Ensure meta_features is a data.table ---
setDT(meta_features_adapt)

# --- 1. Add month_year column using 'start_time' ---
# (This is the column we added in the previous step)
meta_features_adapt[, month_year := format(start_time, "%Y-%m")]

# --- 2. Pre-group data by Stock and month_year (like we did for gamma) ---
# We also pre-filter the data here to match your original logic (>= 2)
stock_month_data <- meta_features_adapt[, 
                                        .(xs_data = list(n_child_orders[!is.na(n_child_orders) & n_child_orders >= 2])), 
                                        by = .(Stock, month_year)
]

# --- 4. Set up cluster ---
ncores <- max(1, detectCores() - 4)
cl <- makeCluster(ncores)

# --- 5. Export necessary data/functions to workers ---
# Ensure workers have data.table and poweRlaw
clusterEvalQ(cl, {
  library(data.table)
  library(poweRlaw)
  NULL
})

# Export the pre-grouped data and the estimation function
clusterExport(cl, varlist = c("stock_month_data", "estimate_alpha_from_vec"), envir = environment())

# --- 6. Run in parallel with progress bar ---
# Iterate over the rows of our new grouped data.table
idx_vec <- seq_len(nrow(stock_month_data))

res_list <- pblapply(cl = cl, X = idx_vec, FUN = function(i) {
  # Get the pre-grouped data vector for this task
  xs <- stock_month_data$xs_data[[i]]
  
  # Run the estimation function
  # You can change the minimum sample size here if needed
  fit_res <- estimate_alpha_from_vec(xs, min_sample_size = 50)
  
  # Return a combined list with Stock/month info
  list(
    Stock = stock_month_data$Stock[i],
    month_year = stock_month_data$month_year[i],
    alpha = fit_res$alpha,
    xmin = fit_res$xmin,
    n_tail = fit_res$n_tail,
    n_obs = fit_res$n_obs
  )
})

# --- 7. Clean up cluster ---
stopCluster(cl)

# --- 8. Combine results ---
# This is your final data.table with alpha estimates per stock-month
results_dt_monthly_adapt <- rbindlist(lapply(res_list, as.data.table))

# --- 9. Show the estimates and plot ---
print(results_dt_monthly_adapt)

# This ggplot code will now plot the histogram of *all*
# the monthly alpha estimates (e.g., all 485 of them).
ggplot(results_dt_monthly_adapt, aes(x = alpha)) +
  geom_histogram(aes(y = ..density..), bins = 42, color = "blue", fill="lightblue", boundary = 0) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "red", size = 1) +
  scale_x_continuous(limits = c(0, NA)) + 
  labs(title = "(b) Adaptive MLB", 
       x = expression(alpha),
       y = expression(P(alpha))) +
  theme_minimal()
# ==============================================================================
# Plotting gamma=alpha-1 (on original MLB assignments) on month-month
# ==============================================================================
# Join alpha (results_dt) and gamma (gamma_dt) by stock AND month
library(data.table)
library(ggplot2)

# 1) Merge into results_FINAL by Stock AND month_year
#    - Use 'results_dt' (your new monthly alpha table)
#    - Merge on both 'Stock' and 'month_year'
results_FINAL_monthly_orig <- merge(results_dt_monthly_orig, gamma_dt_monthly, by = c("Stock", "month_year"))

# ensure numeric and handle missing
results_FINAL_monthly_orig[, alpha := as.numeric(alpha)]
results_FINAL_monthly_orig[, gamma_acf := as.numeric(gamma_acf)]
results_FINAL_monthly_orig[, gamma_psd := as.numeric(gamma_psd)]

# compute alpha_minus1
results_FINAL_monthly_orig[, alpha_minus1 := alpha - 1]

# define bins of width 0.1 on alpha_minus1 (midpoint used as bin key)
bin_width <- 0.1
# compute bin index and bin midpoint
results_FINAL_monthly_orig[, alpha_bin := (floor(alpha_minus1 / bin_width) + 0.5) * bin_width]

# remove rows with NA alpha_minus1 so boxplots ignore them
res_plot <- results_FINAL_monthly_orig[!is.na(alpha_minus1)]

# 3) compute medians per bin for plotting the red median line/points
medians <- res_plot[, .(
  median_acf = median(gamma_acf, na.rm = TRUE),
  median_psd = median(gamma_psd, na.rm = TRUE)
), by = .(alpha_bin)]

# set x column expected by your plotting code
medians[, x := alpha_bin]

# 4) Now run your plotting code (this part was already correct)
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
  coord_cartesian(xlim = c(0, 4.5), ylim = c(0, 4.5)) + 
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
  coord_cartesian(xlim = c(0, 4.5), ylim = c(0, 4.5)) +
  theme_minimal(base_size = 12)

# combine plots
cowplot::plot_grid(p_acf, p_psd, ncol = 2, align = "hv")

# ==============================================================================
# Plotting gamma=alpha-1 (on adaptive MLB assignments) on month-month
# ==============================================================================
# Join alpha (results_dt) and gamma (gamma_dt) by stock AND month
library(data.table)
library(ggplot2)

# 1) Merge into results_FINAL by Stock AND month_year
#    - Use 'results_dt' (your new monthly alpha table)
#    - Merge on both 'Stock' and 'month_year'
results_FINAL_monthly_adapt <- merge(results_dt_monthly_adapt, gamma_dt_monthly, by = c("Stock", "month_year"))

# ensure numeric and handle missing
results_FINAL_monthly_adapt[, alpha := as.numeric(alpha)]
results_FINAL_monthly_adapt[, gamma_acf := as.numeric(gamma_acf)]
results_FINAL_monthly_adapt[, gamma_psd := as.numeric(gamma_psd)]

# compute alpha_minus1
results_FINAL_monthly_adapt[, alpha_minus1 := alpha - 1]

# define bins of width 0.1 on alpha_minus1 (midpoint used as bin key)
bin_width <- 0.1
# compute bin index and bin midpoint
results_FINAL_monthly_adapt[, alpha_bin := (floor(alpha_minus1 / bin_width) + 0.5) * bin_width]

# remove rows with NA alpha_minus1 so boxplots ignore them
res_plot <- results_FINAL_monthly_adapt[!is.na(alpha_minus1)]

# 3) compute medians per bin for plotting the red median line/points
medians <- res_plot[, .(
  median_acf = median(gamma_acf, na.rm = TRUE),
  median_psd = median(gamma_psd, na.rm = TRUE)
), by = .(alpha_bin)]

# set x column expected by your plotting code
medians[, x := alpha_bin]

# 4) Now run your plotting code (this part was already correct)
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
  coord_cartesian(xlim = c(0, 4.5), ylim = c(0, 4.5)) + 
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
  coord_cartesian(xlim = c(0, 4.5), ylim = c(0, 4.5)) +
  theme_minimal(base_size = 12)

# combine plots
cowplot::plot_grid(p_acf, p_psd, ncol = 2, align = "hv")