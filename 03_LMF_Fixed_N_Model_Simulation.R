# ==============================================================================
# SCRIPT: 03_LMF_Fixed_N_Model_Simulation
#
# PURPOSE:
# This script runs a numerical simulation of a multi-agent model to
# demonstrate the "long-range memory" (power-law autocorrelation) of the
# aggregate "sign" (e.g., market order). It is intended to replicate
# the findings of a paper like Lillo, Mike, & Farmer (2005), 
# specifically how the ACF scales with the number of agents (N).
#
# PREREQUISITES:
# * Packages: 'ggplot2', 'scales'.
# * This script is self-contained and does not require external data files.
# * NOTE: The simulation is extremely long ('T = 1e9'). The caching
#   mechanism is essential. The first run will take a very long time.
#
# OUTPUT:
# * A cache file (e.g., "Simulated_workspace_...T1000000000.rds") is
#   created in the working directory to store simulation results.
# * A log-log 'ggplot' object 'p' is created and displayed, comparing
#   simulated and theoretical ACFs.
# * A 'ggsave' command is commented out but can be activated to save the
#   plot as a .png file.
# ==============================================================================

library(ggplot2)
library(scales)

# PARAMETERS
alpha <- 1.5
T     <- 1e9 # truly 10^9 steps 
Ns    <- c(1, 5, 50)
lags  <- round(10^seq(1, 4.2, length.out = 25))

# Pareto sampler
rpareto <- function(n, alpha) floor((1/runif(n))^(1/alpha))

# STREAMING SIMULATION + ACF
simulate_acf <- function(N, alpha, T, lags) {
  maxlag <- max(lags)
  signs   <- sample(c(-1L, 1L), N, TRUE)
  lengths <- rpareto(N, alpha)
  ring    <- integer(maxlag)
  sums    <- numeric(length(lags))
  counts  <- numeric(length(lags))
  
  for (t in seq_len(T)) {
    i <- sample.int(N, 1)
    s <- signs[i]
    lengths[i] <- lengths[i] - 1L
    if (lengths[i] == 0L) {
      signs[i]   <- if (runif(1) < 0.5) -1L else 1L
      lengths[i] <- rpareto(1, alpha)
    }
    
    pos <- ((t - 1L) %% maxlag) + 1L
    ring[pos] <- s
    
    for (j in seq_along(lags)) {
      Lg <- lags[j]
      if (t > Lg) {
        prevpos     <- ((pos - 1L - Lg) %% maxlag) + 1L
        sums[j]    <- sums[j] + s * ring[prevpos]
        counts[j]  <- counts[j] + 1L
      }
    }
  }
  return(sums / counts)
}

# --- CACHING / WORKSPACE SAVE + LOAD ---
# Set to TRUE if you want to force re-running the simulation even if a cache exists
force_rerun <- FALSE

# Construct a filename that reflects key parameters so different runs don't clobber each other
save_file <- file.path(getwd(),
                       paste0("Simulated_workspace_Fig2_T", format(T, scientific = FALSE), ".rds"))

if (file.exists(save_file) && !force_rerun) {
  cat("Loading saved workspace from:", save_file, "\n")
  workspace <- readRDS(save_file)
  # unpack for convenience (optional)
  sim_list    <- workspace$sim_list
  theory_list <- workspace$theory_list
  lags        <- workspace$lags
  Ns          <- workspace$Ns
  alpha       <- workspace$alpha
  T           <- workspace$T
} else {
  cat("No saved workspace found or force_rerun = TRUE. Running simulations...\n")
  
  # RUN FOR EACH N
  sim_list <- lapply(Ns, function(N) {
    cat("Starting N =", N, "...\n")
    acf_vals <- simulate_acf(N, alpha, T, lags)
    data.frame(lag = lags, acf = acf_vals, N = paste0("N=", N))
  })
  
  # THEORETICAL CURVES (Eq.17)
  theory_list <- lapply(Ns, function(N) {
    pref <- N^(alpha - 2) / alpha
    data.frame(lag = lags,
               acf = pref * lags^(-(alpha - 1)),
               N   = paste0("Theory N=", N))
  })
  
  # package workspace and save
  workspace <- list(
    sim_list    = sim_list,
    theory_list = theory_list,
    lags        = lags,
    Ns          = Ns,
    alpha       = alpha,
    T           = T,
    created_at  = Sys.time()
  )
  saveRDS(workspace, file = save_file)
  cat("Saved workspace to:", save_file, "\n")
}

# COMBINE & PLOT
sim_data    <- do.call(rbind, sim_list)
sim_data_trim <- subset(sim_data, lag < max(lags)) # drop the last point to avoid blowup
theory_data <- do.call(rbind, theory_list)

p <- ggplot() +
  geom_line(data = sim_data_trim,    aes(lag, acf, color = N), size = 1) +
  geom_point(data = sim_data_trim,   aes(lag, acf, color = N), shape = 1) +
  geom_line(data = theory_data, aes(lag, acf, linetype = N),
            color = "black", size = 1) +
  scale_x_log10(
    labels = trans_format("log10", math_format(10^.x)),
    breaks = 10^(1:5)
  ) +
  scale_y_log10(
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_color_manual(
    values = c("N=1" = "green", "N=5" = "red", "N=50" = "blue")
  ) +
  labs(
    x = "time lag", y = "acf",
  ) +
  theme_minimal(base_size = 14)+
  theme(
    legend.position = c(0.95, 0.50),       # move to lower-right area
    legend.justification = c(1, 0),       # anchor point inside the legend box
    legend.box.background = element_rect(color = "black", size = 0.3),
    legend.box.margin = margin(0.5, 0.5, 0.5, 0.5)
  )

plot(p)
#ggsave("Figures/FIG2.png", p, width = 9, height = 6, dpi = 300, bg = "white")