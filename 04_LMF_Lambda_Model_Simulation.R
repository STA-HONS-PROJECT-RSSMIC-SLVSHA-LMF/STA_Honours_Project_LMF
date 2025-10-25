# ==============================================================================
# SCRIPT: 04_LMF_Lambda_Model_Simulation
#
# PURPOSE:
# This script simulates the "lambda-model", to replicate
# Figure 4 from a paper (e.g., Lillo, Mike, Farmer). The goal is to
# investigate how the autocorrelation of the *number of active hidden orders*
# changes based on the order arrival rate, lambda.
#
# PREREQUISITES:
# * Packages: 'parallel', 'ggplot2', 'scales'.
# * This script is a self-contained simulation and requires no external data.
# * NOTE: The simulation is long ('T_total'). The caching mechanism is
#   essential. The first run will take time and create a large .rds file.
#
# OUTPUT:
# * A cache file (e.g., "sim_workspace...rds") is created/read.
# * A log-log 'ggplot' object 'p' is created and displayed, comparing the
#   ACF for different 'lambda' values.
# * A 'ggsave' command is commented out but can be activated to save the
#   plot as a .png file.
# ==============================================================================

library(parallel)
library(ggplot2)
library(scales)

# --------------------------- parameters ------------------------------------
alpha   <- 1.50                         # Pareto tail exponent for hidden-order sizes
lambdas <- c(0.33, 0.30, 0.25, 0.20)   
lags    <- unique(as.integer(round(10^seq(0, 4, length.out = 25))))

T_total <- 1e6                          # total measured steps across all workers (start modest)
burn_in <- 1e1                          # steps to reach stationarity (not measured)

# ---------------------- caching / force rerun options ---------------------------
force_rerun <- FALSE
save_file <- file.path(getwd(),
                       paste0("sim_workspace",
                              "_Ttotal", format(as.integer(T_total), big.mark=""),
                              ".rds"))

if (file.exists(save_file) && !force_rerun) {
  message("Loading saved workspace from: ", save_file)
  ws <- readRDS(save_file)
  
  # restore variables (so downstream code expects same names)
  sim_df           <- ws$sim_df
  sim_results_list <- ws$sim_results_list
  lambdas          <- ws$lambdas
  lags             <- ws$lags
  alpha            <- ws$alpha
  T_total          <- ws$T_total
  burn_in          <- ws$burn_in
  T_each           <- ws$T_each
  nrep             <- ws$nrep
  nworkers         <- ws$nworkers
  
  message("Workspace loaded. sim_df and sim_results_list are available.")
} else {
  message("No saved workspace found, or force_rerun = TRUE. Running simulation and caching results...")
  
  # -------------------------- parallel setup (safe) -------------------------------
  nc <- parallel::detectCores()
  if (is.na(nc)) nc <- 1
  nworkers <- max(1, nc - 1)            # leave one core free
  nrep     <- nworkers
  T_each   <- as.integer(max(1, floor(T_total / nrep)))
  
  message("Running ", nrep, " replicate(s) of measured length T_each = ", format(T_each, big.mark=","))
  
  # create cluster and ensure it is stopped on exit / error
  cl <- makeCluster(nworkers)
  on.exit({
    try(stopCluster(cl), silent = TRUE)
  }, add = TRUE)
  
  # reproducible RNG across workers
  set.seed(12345)
  parallel::clusterSetRNGStream(cl, 12345)
  
  # -------------------------- model primitives ------------------------------------
  # Pareto-distributed integer (tail ~ k^(-alpha))
  rpareto_int_one <- function(alpha) pmax(1, floor(runif(1)^(-1/alpha)))
  
  # Main simulator: stream the process; compute centered ACF at requested lags
  simulate_lambda_and_acf_stream <- function(lambda, T, alpha, lags, burn_in = 1e6) {
    T    <- as.integer(T)
    lags <- sort(unique(as.integer(lags)))
    maxL <- max(lags)
    
    # ring buffer to fetch x_{t-L}
    buf_len <- maxL + 1
    buf <- integer(buf_len)
    pos <- 0L
    
    # accumulators for mean/var and cross-products
    sum_x <- 0.0
    sum_x2 <- 0.0
    sum_prod <- numeric(length(lags))
    
    # dynamic array of remaining sizes of active hidden orders
    sizes <- integer(20000L)
    n_hidden <- 0
    
    # run burn-in + measured segment
    for (t in seq_len(burn_in + T)) {
      
      # (1) create a new hidden order with prob λ, drawing its initial size
      if (n_hidden == 0 || runif(1) < lambda) {
        n_hidden <- n_hidden + 1L
        if (n_hidden > length(sizes)) length(sizes) <- as.integer(length(sizes) * 1.5 + 1000)
        sizes[n_hidden] <- rpareto_int_one(alpha)
      }
      
      # (2) execute one unit from a random active hidden order
      if (n_hidden > 0) {
        i <- sample.int(n_hidden, 1)
        sizes[i] <- sizes[i] - 1
        if (sizes[i] <= 0) {                    # remove completed order
          if (i != n_hidden) sizes[i] <- sizes[n_hidden]
          n_hidden <- n_hidden - 1
        }
      }
      
      # (3) current state x_t = number of active hidden orders
      x_t <- as.integer(n_hidden)
      
      # update ring buffer
      pos <- pos %% buf_len + 1
      buf[pos] <- x_t
      
      # skip accumulation during burn-in
      if (t <= burn_in) next
      
      # (4) accumulate moments and cross-products
      sum_x  <- sum_x + x_t
      sum_x2 <- sum_x2 + (x_t * x_t)
      
      for (j in seq_along(lags)) {
        L <- lags[j]
        if ((t - burn_in) > L) {                 # enough history for this lag
          idx <- (pos - L - 1) %% buf_len + 1 # x_{t-L}
          sum_prod[j] <- sum_prod[j] + as.numeric(buf[idx]) * x_t
        }
      }
    }
    
    # centered, normalized autocorrelation at each lag
    mean_x <- sum_x / T
    var_x  <- (sum_x2 - T * mean_x^2) / (T - 1)
    if (!is.finite(var_x) || var_x <= 0) return(setNames(rep(NA_real_, length(lags)), as.character(lags)))
    
    acfs <- sapply(seq_along(lags), function(j) {
      L   <- lags[j]
      num <- sum_prod[j] - (T - L) * mean_x^2      # center: ∑(x_{t-L}x_t) − (T−L)μ^2
      den <- (T - L) * var_x
      if (!is.finite(num) || !is.finite(den) || den <= 0) return(NA_real_)
      num / den
    })
    setNames(acfs, as.character(lags))
  }
  
  # export to workers
  clusterExport(cl,
                c("simulate_lambda_and_acf_stream", "rpareto_int_one", "alpha", "lags", "T_each", "burn_in"),
                envir = environment()
  )
  
  # ------------------------------ run simulation ----------------------------------
  sim_results_list <- lapply(lambdas, function(lambda) {
    message("λ = ", lambda)
    clusterExport(cl, "lambda", envir = environment())
    reps <- parLapply(cl, seq_len(nrep), function(i) {
      simulate_lambda_and_acf_stream(lambda = lambda, T = T_each,
                                     alpha = alpha, lags = lags, burn_in = burn_in)
    })
    mat <- do.call(rbind, reps)
    avg <- colMeans(mat, na.rm = TRUE)
    
    data.frame(
      lag = as.integer(names(avg)),
      acf = as.numeric(avg),
      series = paste0("λ = ", format(lambda, nsmall = 3)),
      legend_series = paste0("λ = ", format(lambda, nsmall = 3)),
      stringsAsFactors = FALSE
    )
  })
  
  # combine
  sim_df <- do.call(rbind, sim_results_list)
  
  # stop the cluster explicitly (on.exit will also try)
  try(stopCluster(cl), silent = TRUE)
  
  # save workspace object
  ws <- list(
    sim_df = sim_df,
    sim_results_list = sim_results_list,
    lambdas = lambdas,
    lags = lags,
    alpha = alpha,
    T_total = T_total,
    burn_in = burn_in,
    T_each = T_each,
    nrep = nrep,
    nworkers = nworkers,
    created_at = Sys.time()
  )
  saveRDS(ws, file = save_file)
  message("Simulation complete. Saved workspace to: ", save_file)
}

# ------------------------- slope guide lines (α − 1) -----------------------------
slope <- alpha - 1
xmin  <- 1
xmax  <- 2e4
ymin  <- 0.03

make_guide <- function(C) {
  # draw a straight line with slope (α−1) on log–log scale, clipped at ymin
  lag_hit <- (C / ymin)^(1 / slope)
  eff_max <- min(xmax, lag_hit)
  if (eff_max <= xmin) return(NULL)
  lseq <- unique(as.integer(round(10^seq(log10(xmin), log10(eff_max), length.out = 400))))
  acf_vals <- C * lseq^(-slope)
  if (lag_hit <= xmax && lag_hit > max(lseq)) { lseq <- c(lseq, lag_hit); acf_vals <- c(acf_vals, ymin) }
  data.frame(lag = lseq, acf = acf_vals,
             series = paste0("guide_", formatC(C, format = "g", digits = 6)),
             legend_series = "Slope = α − 1",
             stringsAsFactors = FALSE)
}

# compute guide constants so each guide sits under each plotted curve
lambda_labels <- unique(sim_df$legend_series)   # <- used by curve plotting later as before
under_factor <- 0.6   # fraction of the curve value where the guide will be anchored (0 < f < 1)

guide_list <- lapply(lambda_labels, function(lbl) {
  sim_sub <- subset(sim_df, legend_series == lbl)
  if (nrow(sim_sub) == 0) return(NULL)
  # choose a reference lag representative of that curve (median within visible range)
  valid_lags <- sim_sub$lag[sim_sub$lag >= xmin & sim_sub$lag <= xmax]
  if (length(valid_lags) == 0) valid_lags <- sim_sub$lag
  ref_lag <- median(valid_lags, na.rm = TRUE)
  # get curve value at reference lag (approx if necessary)
  acf_ref <- approx(x = sim_sub$lag, y = sim_sub$acf, xout = ref_lag, rule = 2)$y
  # set C so that guide(ref_lag) = under_factor * acf_ref  =>  C = (under_factor * acf_ref) * ref_lag^slope
  if (is.na(acf_ref) || acf_ref <= 0) return(NULL)
  C <- (under_factor * acf_ref) * (ref_lag^slope)
  make_guide(C)
})

guide_df <- do.call(rbind, Filter(Negate(is.null), guide_list))
plot_df  <- rbind(sim_df, guide_df)

# combine with simulation data (handle case where no guides were created)
if (is.null(guide_df) || nrow(guide_df) == 0) {
  plot_df <- sim_df
} else {
  plot_df <- rbind(sim_df, guide_df)
}
plot_df  <- rbind(sim_df, guide_df)

# legend aesthetics
lambda_labels     <- unique(sim_df$legend_series)
legend_order      <- c(lambda_labels, "Slope = α − 1")
plot_df$legend_series <- factor(plot_df$legend_series, levels = legend_order)

cols <- c("orange", "red", "green4", "blue", "black"); names(cols) <- legend_order
shps <- c(1, 0, 2, 6, NA);                       names(shps) <- legend_order
ltys <- c(rep("solid", length(lambda_labels)), "dashed"); names(ltys) <- legend_order

# ------------------------------- plot -------------------------------------------
p <- ggplot(plot_df, aes(x = lag, y = acf,
                         color = legend_series,
                         linetype = legend_series,
                         group = series)) +
  geom_line(linewidth = 1) +
  geom_point(data = subset(plot_df, legend_series %in% lambda_labels),
             aes(shape = legend_series), size = 2, stroke = 0.9) +
  scale_x_log10(breaks = c(1,10,100,1000,10000),
                labels = c(1, 10, 100,1000, 10000),
                expand = c(0,0)) +
  scale_y_log10(breaks = c(1, 1e-1, 1e-2),
                labels = c(1, 0.1, expression(0.01)),
                expand = c(0,0)) +
  scale_color_manual(values = cols, breaks = legend_order, name = NULL) +
  scale_shape_manual(values = shps, breaks = legend_order, name = NULL) +
  scale_linetype_manual(values = ltys, breaks = legend_order, name = NULL) +
  coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, 1.05), expand = FALSE) +
  labs(x = "lag", y = "acf") +
  theme_classic(base_size = 16) +
  theme(text = element_text(family = "serif"),
        axis.title = element_text(size = 18),
        axis.text  = element_text(size = 14),
        legend.position = c(0.12, 0.12),
        legend.justification = c("left", "bottom"),
        legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(shape = shps, linetype = ltys,
                                                  linewidth = c(rep(1.2, length(lambda_labels)), 1.0))),
         shape = "none", linetype = "none")

print(p)

#ggsave("Figures/FIG4.png", p, width = 9, height = 6, dpi = 300, bg = "white")
# ================================================================================