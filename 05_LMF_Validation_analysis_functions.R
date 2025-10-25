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

estimate_gamma_acf <- function(signs, lag_min = 5, lag_max_frac = 0.2, nbins = 20) {
  N <- length(signs)
  if (N < (lag_min + 10)) return(NA_real_)
  acf_res <- acf(signs, plot = FALSE, lag.max = floor(N * lag_max_frac))
  acf_vals <- as.numeric(acf_res$acf)[-1]
  lags <- as.numeric(acf_res$lag)[-1]
  ok <- which(lags >= lag_min & acf_vals > 0)
  if (length(ok) < 8) return(NA_real_)
  lags_ok <- lags[ok]; acf_ok <- acf_vals[ok]
  bin_edges <- unique(floor(exp(seq(log(min(lags_ok)), log(max(lags_ok)), length.out = nbins + 1))))
  bidx <- findInterval(lags_ok, bin_edges, rightmost.closed = TRUE)
  bmeans <- tapply(acf_ok, bidx, mean)
  blags  <- tapply(lags_ok, bidx, mean)
  keep <- !is.na(bmeans) & !is.na(blags) & (bmeans > 0)
  if (sum(keep) < 6) return(NA_real_)
  df <- data.frame(lag = blags[keep], acf = bmeans[keep])
  fit <- lm(log(acf) ~ log(lag), data = df)
  slope <- coef(fit)[2]
  gamma_hat <- -slope
  return(as.numeric(gamma_hat))
}

estimate_gamma_psd <- function(signs, low_freq_frac = 0.12) {
  sp <- spec.pgram(signs, taper = 0, plot = FALSE, demean = TRUE, fast = TRUE)
  freq <- sp$freq; spec <- sp$spec
  cutoff <- max(freq) * low_freq_frac
  ok <- which(freq > 0 & freq <= cutoff)
  if (length(ok) < 6) return(NA_real_)
  fit <- lm(log(spec[ok]) ~ log(freq[ok]))
  slope <- coef(fit)[2]
  # relation used in prior script: gamma = 1 + slope (since slope is negative when PSD ~ f^{-beta})
  gamma_hat <- 1 + slope
  return(as.numeric(gamma_hat))
}

mapping_function_orig <- function(N, F_sampler, theta, n_orders, seed=123, ...){
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

mapping_function_adapt <- function(N, F_sampler, theta, n_orders, seed = 123,
                                   adaptive = TRUE, base_block = 100, f_min = 1, ...) {
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

sample_powerlaw <- function(N, theta, f_min = 1, ...) { # NOTE: higher thetas means more fragmentation between traders.
  u <- runif(N)                                     
  #f <- f_min * ((1 - u) ^ (-1 / theta))     
  #f_min * (u ^ (-1 / (theta - 1)))
  f <- f_min*(u)^(-1/theta)
  return(f)
}

estimate_alpha_from_vec <- function(xs, min_sample_size = 50) {
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