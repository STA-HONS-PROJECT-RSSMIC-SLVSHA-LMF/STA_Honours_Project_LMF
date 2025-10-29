# ==============================================================================
# SCRIPT: 05_LMF_Validation_analysis_functions
#
# PURPOSE:
# This script is a library of all custom functions required for the
# main $\gamma$ (gamma) vs. $\alpha$ (alpha) analysis. It defines
# the estimators for $\gamma$ and $\alpha$, as well as the "metaorder"
# simulation models (MLB).
#
# This script defines the following key functions:
# - estimate_gamma_acf(): Estimates $\gamma$ from the ACF of a sign series.
# - estimate_gamma_psd(): Estimates $\gamma$ from the PSD of a sign series.
# - sample_powerlaw(): Sampler for the metaorder simulation.
# - mapping_function_orig(): The "Original MLB" model.
# - mapping_function_adapt(): The "Adaptive MLB" model.
# - estimate_alpha_from_vec(): A helper to estimate $\alpha$ using the
#   Clauset method (via the 'poweRlaw' package).
#
# PREREQUISITES:
# * None. This script is a self-contained library of definitions. It does
#   not run any analysis itself.
#
# OUTPUT:
# * None. This script simply loads functions into the R environment when
#   `source()`-ed by other scripts.
# ==============================================================================

estimate_gamma_acf <- function(signs, 
                               lag_min = 5, 
                               lag_max_frac = 0.2, 
                               nbins = 20) {
  N <- length(signs)
  if (N < (lag_min + 10)) return(NA_real_)
  acf_res <- acf(signs, plot = FALSE, lag.max = floor(N * lag_max_frac))
  acf_vals <- as.numeric(acf_res$acf)[-1]
  lags <- as.numeric(acf_res$lag)[-1]
  ok <- which(lags >= lag_min & acf_vals > 0)
  if (length(ok) < 8) return(NA_real_)
  lags_ok <- lags[ok]; acf_ok <- acf_vals[ok]
  
  # Binning data
  bin_edges <- unique(floor(exp(seq(log(min(lags_ok)), log(max(lags_ok)), length.out = nbins + 1))))
  bidx <- findInterval(lags_ok, bin_edges, rightmost.closed = TRUE)
  bmeans <- tapply(acf_ok, bidx, mean)
  blags  <- tapply(lags_ok, bidx, mean)
  keep <- !is.na(bmeans) & !is.na(blags) & (bmeans > 0)
  if (sum(keep) < 6) return(NA_real_)
  df <- data.frame(lag = blags[keep], acf = bmeans[keep])
  
  # --- NLLS FIT (Step 1) ---
  # Get starting values from the log-log linear fit
  fit_lm <- NULL
  try(fit_lm <- lm(log(acf) ~ log(lag), data = df), silent = TRUE)
  if (is.null(fit_lm)) return(NA_real_)
  start_slope <- coef(fit_lm)[2]
  start_intercept <- coef(fit_lm)[1]
  
  if (is.na(start_slope) || !is.finite(start_slope) || 
      is.na(start_intercept) || !is.finite(start_intercept)) {
    return(NA_real_)
  }
  
  start_gamma <- -start_slope
  start_C <- exp(start_intercept)
  if (start_gamma <= 0) start_gamma <- 0.5 
  if (start_C <= 0) start_C <- 1.0       
  # Fit the non-linear model to get gamma_NLLS
  fit_nls <- NULL
  tryCatch({
    fit_nls <- nls(
      formula = acf ~ C * lag^(-gamma),
      data = df,
      start = list(C = start_C, gamma = start_gamma),
      control = nls.control(maxiter = 100, warnOnly = TRUE)
    )
  }, error = function(e) {
    # If nls fails, fit_nls remains NULL
  })
  
  if (is.null(fit_nls)) {
    return(NA_real_)
  }
  # --- HEURISTIC BIAS CORRECTION (Step 2) ---
  # This has been added in for the purposes of the researcher to attempt to run
  # a bias correction if they have sufficient data. We did not run this code when
  # implemting our results in the report, but have decided to add it here to help
  # future researchers see how it would be applied (assuming Data specific Beta's
  # have been found.
  # 1. This is the naive NLLS estimate
  gamma_nlls <- coef(fit_nls)["gamma"]
  # 2. Define heuristic correction parameters (below we use from Sato & Kanazawa (2023), Fig. 8
  # as an example)
  beta_1 <- 0.592
  beta_2 <- 0.147
  # 3. Apply Equation (24) to get the unbiased estimate
  gamma_unbiased <- (gamma_nlls - beta_2) / beta_1
  return(as.numeric(gamma_unbiased))
}

estimate_gamma_psd <- function(signs,
                               min_freq_bin = 16,
                               low_freq_frac = 0.002, # Default changed to 0.002
                               smooth_window = 11) {
  library(zoo)
  if (length(signs) < 250) {
    # Not enough data for a meaningful spectrum
    return(NA_real_)
  }
  
  # 1. Calculate Periodogram
  sp <- try(spec.pgram(signs, taper = 0, plot = FALSE, demean = TRUE, fast = TRUE), silent = TRUE)
  if (inherits(sp, "try-error")) {
    return(NA_real_)
  }
  
  freq <- sp$freq
  spec <- sp$spec
  
  if (length(freq) < min_freq_bin) {
    # Spectrum is shorter than our minimum bin cutoff
    return(NA_real_)
  }
  
  # 2. Apply Smoothing
  spec_smooth <- as.numeric(zoo::rollmean(spec, k = smooth_window, fill = NA))
  
  # 3. Define Frequency Range
  # Lower cutoff: Get the frequency of the 16th bin [cite: 1136]
  min_f <- freq[min_freq_bin]
  
  # Upper cutoff: Paper suggests 10^-3 
  # max(freq) is 0.5, so 0.001 / 0.5 = 0.002
  max_f <- max(freq) * low_freq_frac
  
  if (min_f >= max_f) {
    # This can happen if N is very small or low_freq_frac is tiny
    return(NA_real_)
  }
  
  # 4. Filter for valid data points
  # We add (spec_smooth > 0) to prevent log(0) errors
  ok <- which(!is.na(spec_smooth) & spec_smooth > 0 & freq >= min_f & freq <= max_f)
  
  # Need at least a few points to fit a line (e.g., 10)
  if (length(ok) < 10) {
    return(NA_real_)
  }
  
  # 5. Fit model
  fit <- lm(log(spec_smooth[ok]) ~ log(freq[ok]))
  slope <- coef(fit)[2]
  
  # 6. Calculate gamma
  # S(f) ~ f^(gamma-1), so log(S) = (gamma-1) * log(f)
  gamma_hat <- 1 + slope
  
  # Final check: slope might be NA/Inf if fit is singular
  if (is.na(gamma_hat) | !is.finite(gamma_hat)) {
    return(NA_real_)
  }
  
  return(as.numeric(gamma_hat))
}

run_parallel_gamma_estimation <- function(grouped_signs_dt,
                                          grouping_cols,
                                          min_signs_for_gamma,
                                          lag_min,
                                          lag_max_frac,
                                          nbins,
                                          low_freq_frac,
                                          ncores_to_use) {
  
  # --- Load required libraries ---
  require(parallel)
  require(pbapply)
  
  # --- cluster setup ---
  cl <- makeCluster(ncores_to_use)
  
  # ensure needed packages are loaded on workers
  clusterEvalQ(cl, { library(data.table); NULL })
  
  # Export data, functions and parameter variables to workers
  clusterExport(cl, varlist = c(
    "grouped_signs_dt",
    "estimate_gamma_acf",
    "estimate_gamma_psd",
    "min_signs_for_gamma",
    "lag_min",
    "lag_max_frac",
    "nbins",
    "low_freq_frac",
    "grouping_cols"
  ), envir = environment())
  
  # --- run computations in parallel (load-balanced) ---
  idx_vec <- seq_len(nrow(grouped_signs_dt)) 
  
  res_list <- pblapply(cl = cl, X = idx_vec, FUN = function(idx) {
    # <--- worker-side: use only exported objects and functions --->
    
    # Dynamically get grouping data (e.g., Stock, month_year)
    group_data <- as.list(grouped_signs_dt[idx, ..grouping_cols])
    
    signs <- grouped_signs_dt$signs[[idx]]
    
    if (length(signs) < min_signs_for_gamma) {
      gamma_acf <- NA_real_; gamma_psd <- NA_real_
    } else {
      gamma_acf <- tryCatch(estimate_gamma_acf(signs, lag_min = lag_min, lag_max_frac = lag_max_frac, nbins = nbins),
                            error = function(e) NA_real_)
      gamma_psd <- tryCatch(estimate_gamma_psd(signs, low_freq_frac = low_freq_frac), 
                            error = function(e) NA_real_)  
    }
    
    # Combine dynamic group data with results
    c(group_data, 
      list(gamma_acf = gamma_acf, gamma_psd = gamma_psd, n_orders = length(signs))
    )
  })
  
  # --- clean up ---
  stopCluster(cl)
  
  # Combine results into a data.table
  rbindlist(lapply(res_list, as.data.table))
}

generate_metaorders_orig <- function(DT, 
                                     N, 
                                     theta, 
                                     f_min = 1, 
                                     seed = 123) {
  
  # Compute daily stats
  daily_stats_orig <- DT[, .(
    VD   = sum(volume, na.rm = TRUE), 
    p0   = price[1L], 
    pmax = max(price, na.rm = TRUE),
    pmin = min(price, na.rm = TRUE)
  ), by = .(Stock, trade_date)]
  
  # Daily range volatility proxy
  daily_stats_orig[, sigmaD := (pmax - pmin) / p0]
  
  DT_origMLB <- daily_stats_orig[DT, on = .(Stock, trade_date)]
  
  daily_stats_orig <- daily_stats_orig[DT_origMLB[, .(Stock, trade_date)] |> unique(), on=.(Stock, trade_date)]
  
  # First, preserve true chronological order of trades by sorting by tradestamp
  setorder(DT_origMLB, timestamp)
  
  # Set the seed ONCE before the group assignment for reproducibility
  set.seed(seed) 
  
  # Call the mapping function FOR EACH STOCK-DAY GROUP
  DT_origMLB[, trader_id := mapping_function_orig(
    N         = N,
    F_sampler = sample_powerlaw,
    theta     = theta,
    n_orders  = .N,     # .N is the number of rows (orders) in the current group
    seed      = NULL,   # Use NULL so the RNG continues and is not reset for every group
    f_min     = f_min
  )$assignments, 
  by = .(Stock, trade_date)] # This performs the assignment per stock and per day
  
  # Step 4: Here we need to now sort trades by traders and timestamp
  setorder(DT_origMLB, trader_id, timestamp)
  # Step 5: Now we define a metaorder as a sequence of trades of same sign from same traders
  DT_origMLB[, metaorder_id := rleid(trader_id, Sign)] 
  # Step 6: 
  meta_features_orig <- DT_origMLB[, .(
    trader_id        = first(trader_id),
    Stock            = first(Stock),
    Sign             = first(Sign),
    start_time       = first(timestamp),
    end_time         = last(timestamp),
    n_child_orders   = .N,
    total_volume     = sum(volume)
  ), by = metaorder_id][n_child_orders > 1]
  
  pareto_label <- paste0("Pareto Fit (theta=", theta, ")")
  # Generate the plot
  p <- ggplot(meta_features_orig, aes(x = n_child_orders)) +
    geom_histogram(bins = 50, fill = "lightblue", color = "black") +
    stat_function(
      fun = pareto_pdf_scaled,
      args = list(alpha = theta, C = 1000000),
      # Use the variable directly for the 'color' aesthetic
      aes(color = pareto_label), 
      linetype = "dashed",
      linewidth = 1
    ) +
    scale_color_manual(
      name = "Distribution Fit",
      values = setNames("red", pareto_label)
    ) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    labs(
      title = "(a) Original MLB", 
      x = "Number of Child Orders (Log Scale)",
      y = "Count (Log Scale)"
    ) +
    theme_minimal() +
    theme(
      legend.position = c(0.95, 0.95),
      legend.justification = c("right", "top")
    )
  
  print(p)
  
  return(meta_features_orig)
}

generate_metaorders_adapt <- function(DT, 
                                      N, 
                                      theta, 
                                      base_block, 
                                      f_min = 1, 
                                      seed = 123) {
  
  # Compute daily stats
  daily_stats_adapt <- DT[, .(
    VD   = sum(volume, na.rm = TRUE), 
    p0   = price[1L], 
    pmax = max(price, na.rm = TRUE),
    pmin = min(price, na.rm = TRUE)
  ), by = .(Stock, trade_date)]
  
  # Daily range volatility proxy
  daily_stats_adapt[, sigmaD := (pmax - pmin) / p0]
  
  DT_adaptMLB <- daily_stats_adapt[DT, on = .(Stock, trade_date)]
  
  daily_stats_adapt <- daily_stats_adapt[DT_adaptMLB[, .(Stock, trade_date)] |> unique(), on=.(Stock, trade_date)]
  
  # Set the seed ONCE before the group assignment for reproducibility
  set.seed(seed) 
  
  # Call the mapping function FOR EACH STOCK-DAY GROUP
  DT_adaptMLB[, trader_id := mapping_function_adapt(
    N         = N,
    F_sampler = sample_powerlaw,
    theta     = theta,
    n_orders  = .N,     # .N is the number of rows (orders) in the current group
    seed      = NULL,   # Use NULL so the RNG continues and is not reset for every group
    adaptive = TRUE,    # enable adaptive block sizes
    base_block = base_block,   
    f_min     = f_min
  )$assignments, 
  by = .(Stock, trade_date)] # This performs the assignment per stock and per day
  
  # Step 4: Here we need to now sort trades by traders and timestamp
  setorder(DT_adaptMLB, trader_id, timestamp)
  # Step 5: Now we define a metaorder as a sequence of trades of same sign from same traders
  DT_adaptMLB[, metaorder_id := rleid(trader_id, Sign)] 
  # Step 6: 
  meta_features_adapt <- DT_adaptMLB[, .(
    trader_id        = first(trader_id),
    Stock            = first(Stock),
    Sign             = first(Sign),
    start_time       = first(timestamp),
    end_time         = last(timestamp),
    n_child_orders   = .N,
    total_volume     = sum(volume)
  ), by = metaorder_id][n_child_orders > 1]
  
  pareto_label <- paste0("Pareto Fit (theta=", theta, ")")
  # Generate the plot
  p <- ggplot(meta_features_adapt, aes(x = n_child_orders)) +
    geom_histogram(bins = 50, fill = "lightblue", color = "black") +
    stat_function(
      fun = pareto_pdf_scaled,
      args = list(alpha = theta, C = 1000000),
      aes(color = pareto_label), 
      linetype = "dashed",
      linewidth = 1
    ) +
    scale_color_manual(
      name = "Distribution Fit",
      values = setNames("red", pareto_label)
    ) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    labs(
      title = "(b) Adaptive MLB", 
      x = "Number of Child Orders (Log Scale)",
      y = "Count (Log Scale)"
    ) +
    theme_minimal() +
    theme(
      legend.position = c(0.95, 0.95),
      legend.justification = c("right", "top")
    )
  
  print(p)
  
  return(meta_features_adapt)
}

mapping_function_orig <- function(N, 
                                  F_sampler, 
                                  theta, 
                                  n_orders, 
                                  seed=123, 
                                  ...){
  # ----- Step 1 of Algorithm 2 ----- #
  # Initialisation: set N (number of traders) and choose sampler F
  if (!is.null(seed)) set.seed(seed)
  N <- as.integer(N)
  # ----- Step 2 of Algorithm 2 ----- #
  # Generate f_i from probability law F for i = 1,2,...,N
  f <- F_sampler(N, theta, ...)
  # ----- Step 3 of Algorithm 2 ----- #
  # Define p_i = f_i / sum_j f_i
  p <- f / sum(f)
  # ----- Step 4 of Algorithm 2 ----- #
  # Compute cumulative probabilities c_i with c_0 = 0
  c <- cumsum(p)
  # ----- Step 5 of Algorithm 2 ----- #
  # For each order in the market
  assignments <- sample(1:N, size = n_orders, replace = TRUE, prob = p)
  # ----- Step 9 of Algorithm 2 ----- #
  # End loop and prepare summary
  counts <- tabulate(assignments, nbins = N)
  summary <- data.frame( trader = 1:N, orders = counts)
  return(list(assignments = assignments, summary = summary))
}

mapping_function_adapt <- function(N, 
                                   F_sampler, 
                                   heta, 
                                   n_orders, 
                                   seed = 123,
                                   adaptive = TRUE, 
                                   base_block = 100, 
                                   f_min = 1, 
                                   ...) {
  # ----- Step 1: init -----
  if (!is.null(seed)) set.seed(seed)
  N <- as.integer(N)
  # ----- Step 2: sample weights and probs -----
  f <- F_sampler(N, theta, f_min = f_min, ...)
  p <- f / sum(f)
  # ----- Step 3: determine how many trades each trader gets -----
  n_trades <- round(p * n_orders)
  diff_trades <- n_orders - sum(n_trades)
  if (diff_trades > 0) {
    extra_indices <- sample(1:N, diff_trades, replace = TRUE)
    for (i in extra_indices) n_trades[i] <- n_trades[i] + 1
  } else if (diff_trades < 0) {
    remove_indices <- sample(which(n_trades > 0), -diff_trades)
    for (i in remove_indices) n_trades[i] <- n_trades[i] - 1
  }
  # ----- Step 4: build blocks adaptively (or uniformly if adaptive = FALSE) -----
  # We'll create a list where each element is a block (vector of same trader id)
  blocks <- vector("list", 0)
  for (i in seq_len(N)) {
    ni <- n_trades[i]
    if (ni <= 0) next
    
    if (adaptive) {
      # block_size scales with trader probability relative to mean probability
      # ensure block_size at least 1 and not huge
      rel <- p[i] / mean(p)
      # tuned formula: base_block scaled by rel^(0.5) to avoid extreme blocks
      blk_sz <- max(1, round(base_block * (rel ^ 0.5)))
    } else {
      blk_sz <- base_block
    }
    
    # split ni into contiguous blocks of size ~ blk_sz
    remaining <- ni
    while (remaining > 0) {
      take <- min(remaining, blk_sz)
      blocks[[length(blocks) + 1L]] <- rep(i, take)
      remaining <- remaining - take
      # optionally vary blk_sz slightly to add randomness
      if (adaptive) {
        # small jitter so blocks are not identical sizes
        blk_sz <- max(1, round(blk_sz * runif(1, 0.7, 1.3)))
      }
    }
  }
  
  # ----- Step 5: shuffle the blocks (preserves intra-block contiguity) -----
  perm <- sample(length(blocks))
  blocks_shuffled <- blocks[perm]
  
  # ----- Step 6: flatten blocks into final assignment vector -----
  assignments <- unlist(blocks_shuffled, use.names = FALSE)
  
  # Safety check: ensure length == n_orders (trim or pad if needed)
  L <- length(assignments)
  if (L > n_orders) {
    # trim extra at the end
    assignments <- assignments[1:n_orders]
  } else if (L < n_orders) {
    # pad by sampling according to p (fallback)
    pad_needed <- n_orders - L
    pad <- sample(1:N, pad_needed, replace = TRUE, prob = p)
    assignments <- c(assignments, pad)
  }
  
  # ----- Step 7: prepare summary -----
  counts <- tabulate(assignments, nbins = N)
  summary <- data.frame(trader = seq_len(N), orders = counts)
  
  return(list(assignments = assignments, summary = summary))
}

sample_powerlaw <- function(N, 
                            theta, 
                            f_min = 1, 
                            ...) { # NOTE: higher thetas means more fragmentation between traders.
  u <- runif(N)                                     
  f <- f_min*(u)^(-1/theta)
  return(f)
}

estimate_alpha_by_stock_loop <- function(meta_features_dt, 
                                         min_sample_size = 50, 
                                         xmin_filter = 2) {
  
  require(poweRlaw)
  setDT(meta_features_dt)
  stocks <- unique(meta_features_dt$Stock)
  
  results_dt <- data.table(
    Stock = stocks,
    alpha = as.numeric(NA),
    xmin  = as.integer(NA),
    n_tail = as.integer(NA)
  )
  
  for (i in seq_along(stocks)) {
    st <- stocks[i]
    # extract integer metaorder lengths for this stock
    xs <- meta_features_dt[Stock == st, n_child_orders]
    # discard NA or non-positive values, keep integers
    xs <- xs[!is.na(xs) & xs >= xmin_filter] 
    # need at least a small sample to fit; otherwise leave NA
    if (length(xs) < min_sample_size) {
      results_dt[i, `:=`(alpha = NA_real_, xmin = NA_integer_, n_tail = length(xs))]
      next
    }
    
    # fit discrete power law with poweRlaw
    fit_alpha <- tryCatch({
      m <- displ$new(xs)                 # discrete power-law object
      est <- estimate_xmin(m)            # estimate xmin and pars (Clauset method)
      # set estimated xmin and pars to object (estimate_xmin returns a list-like object)
      m$setXmin(est$xmin)
      m$setPars(est$pars)
      xmin_val <- est$xmin
      alpha_val <- est$pars
      n_tail_val <- sum(xs >= xmin_val)
      list(alpha = as.numeric(alpha_val), xmin = as.integer(xmin_val), n_tail = as.integer(n_tail_val))
    }, error = function(e) {
      # on failure, return NA (with n_tail = length(xs) as per original script)
      list(alpha = NA_real_, xmin = NA_integer_, n_tail = as.integer(length(xs)))
    })
    
    results_dt[i, `:=`(alpha = fit_alpha$alpha, xmin = fit_alpha$xmin, n_tail = fit_alpha$n_tail)]
  }
  
  return(results_dt)
}

estimate_alpha_from_vec <- function(xs, 
                                    min_sample_size = 50) {
  n_xs <- length(xs)
  
  # Check for minimum sample size
  if (n_xs < min_sample_size) {
    # Return NAs but include the observation count
    return(list(alpha = NA_real_, xmin = NA_integer_, n_tail = NA_integer_, n_obs = n_xs))
  }
  
  # Fit discrete power law using tryCatch
  tryCatch({
    m <- displ$new(xs)
    est <- estimate_xmin(m) # Clauset method
    m$setXmin(est$xmin)
    m$setPars(est$pars)
    
    list(
      alpha = as.numeric(est$pars), 
      xmin = as.integer(est$xmin), 
      n_tail = as.integer(sum(xs >= est$xmin)), # Count of obs in tail
      n_obs = n_xs # Total pre-fit observations
    )
  }, error = function(e) {
    # On failure, return NAs
    list(alpha = NA_real_, xmin = NA_integer_, n_tail = NA_integer_, n_obs = n_xs)
  })
}

run_parallel_alpha_estimation <- function(meta_features_dt,
                                                  grouping_cols,
                                                  min_sample_size = 50,
                                                  ncores_to_use,
                                                  xmin_filter = 2) {
  
  # --- Load required libraries ---
  require(poweRlaw)
  require(data.table)
  require(parallel)
  require(pbapply)
  
  # --- Ensure meta_features is a data.table ---
  setDT(meta_features_dt)
  
  # --- Pre-group data ---
  stock_month_data <- meta_features_dt[,
                                       .(xs_data = list(n_child_orders[!is.na(n_child_orders) &
                                                                          n_child_orders >= xmin_filter])), 
                                       by = c(grouping_cols)
  ]
  
  # --- Set up cluster ---
  cl <- makeCluster(ncores_to_use)
  
  # --- Export necessary data/functions to workers ---
  clusterEvalQ(cl, {
    library(data.table)
    library(poweRlaw)
    NULL
  })
  
  clusterExport(cl, varlist = c("stock_month_data", "estimate_alpha_from_vec", 
                                "min_sample_size", "grouping_cols"), 
                envir = environment())
  
  # --- Run in parallel with progress bar ---
  idx_vec <- seq_len(nrow(stock_month_data))
  
  res_list <- pblapply(cl = cl, X = idx_vec, FUN = function(i) {
    # Get the pre-grouped data vector for this task
    xs <- stock_month_data$xs_data[[i]]
    
    # Run the estimation function
    fit_res <- estimate_alpha_from_vec(xs, min_sample_size = min_sample_size)
    
    # Get grouping data
    group_data <- as.list(stock_month_data[i, ..grouping_cols])
    
    # Return a combined list
    c(group_data, fit_res)
  })
  
  # --- Clean up cluster ---
  stopCluster(cl)
  
  # --- Combine results ---
  results_dt <- rbindlist(lapply(res_list, as.data.table))
  return(results_dt)
}

pareto_pdf_scaled <- function(x, 
                              alpha, 
                              C) {
  x[x < 1] <- 1
  return(C * alpha / (x^(alpha + 1)))
}

plot_gamma_alpha_comparison <- function(alpha_dt,
                                        gamma_dt,
                                        merge_cols,
                                        coord_limits,
                                        bin_width = 0.1) {
  
  require(data.table)
  require(ggplot2)
  
  # 1) Merge results
  results_FINAL <- merge(alpha_dt, gamma_dt, by = merge_cols)
  
  # ensure numeric and handle missing
  results_FINAL[, alpha := as.numeric(alpha)]
  results_FINAL[, gamma_acf := as.numeric(gamma_acf)]
  results_FINAL[, gamma_psd := as.numeric(gamma_psd)]
  
  # 2) compute alpha_minus1
  results_FINAL[, alpha_minus1 := alpha - 1]
  
  # compute bin index and bin midpoint
  results_FINAL[, alpha_bin := (floor(alpha_minus1 / bin_width) + 0.5) * bin_width]
  
  # remove rows with NA alpha_minus1 so boxplots ignore them
  res_plot <- results_FINAL[!is.na(alpha_minus1)]
  
  # 3) compute medians per bin for plotting the red median line/points
  medians <- res_plot[, .(
    median_acf = median(gamma_acf, na.rm = TRUE),
    median_psd = median(gamma_psd, na.rm = TRUE)
  ), by = .(alpha_bin)]
  
  # set x column expected by plotting code
  medians[, x := alpha_bin]
  
  # 4) Plot ACF
  p_acf <- ggplot(res_plot, aes(x = alpha_minus1, y = gamma_acf)) +
    geom_boxplot(aes(group = alpha_bin), width = 0.05, outlier.shape = NA,
                 fill = "gray80", color = "black") +
    geom_point(alpha = 0.3, color = "orange", size = 1.4) +
    geom_line(data = medians, aes(x = x, y = median_acf), color = "red", linewidth = 0.5) +
    geom_point(data = medians, aes(x = x, y = median_acf), shape = 25, fill = "red", size = 1) +
    geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.5) +
    labs(x = expression(alpha - 1),
         y = expression(gamma^{(ACF)}[biased]),
         title = "(a) ACF method") +
    coord_cartesian(xlim = c(coord_limits[1], coord_limits[2]), 
                    ylim = c(coord_limits[3], coord_limits[4])) + 
    theme_minimal(base_size = 12)
  
  # 5) Plot PSD
  p_psd <- ggplot(res_plot, aes(x = alpha_minus1, y = gamma_psd)) +
    geom_boxplot(aes(group = alpha_bin), width = 0.05, outlier.shape = NA,
                 fill = "gray80", color = "black") +
    geom_point(alpha = 0.3, color = "orange", size = 1.4) +
    geom_line(data = medians, aes(x = x, y = median_psd), color = "red", linewidth = 0.5) +
    geom_point(data = medians, aes(x = x, y = median_psd), shape = 25, fill = "red", size = 1) +
    geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.5) +
    labs(x = expression(alpha - 1),
         y = expression(gamma^{(PSD)}[biased]),
         title = "(b) PSD method") +
    coord_cartesian(xlim = c(coord_limits[1], coord_limits[2]), 
                    ylim = c(coord_limits[3], coord_limits[4])) +
    theme_minimal(base_size = 12)
  
  # 6) combine plots
  print(cowplot::plot_grid(p_acf, p_psd, ncol = 2, align = "hv"))
}
