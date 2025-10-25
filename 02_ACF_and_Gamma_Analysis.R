# ==============================================================================
# SCRIPT: 02_ACF_and_Gamma_Analysis
#
# PURPOSE:
# This script performs a complete autocorrelation analysis on the trade-sign
# series for a single, specified stock (hardcoded as Naspers, "NPNJn.J").
# The goal is to replicate the methodology from high-frequency finance
# literature (e.g., Lillo, Mike, Farmer) to demonstrate and quantify
# "long-range memory" in market dynamics. The analysis is conducted in three 
#  main stages, all within this script:
# 1.  Calculate the Autocorrelation Function (ACF) of the trade-sign series.
# 2.  Plot the resulting ACF on a log-log scale to visualize the decay.
# 3.  Estimate the power-law exponent 'gamma' from the ACF decay,
#     quantifying the strength of the long-range memory.
#
# PREREQUISITES:
# * Packages: 'ggplot2', 'dplyr'.
# * Input Data: A *specific, single-instrument CSV file* is hardcoded
#   (".../NPNJn.J28-May-2015_TO_31-Dec-2015.csv").
# * Data Structure: Assumes the CSV has a 'Type' column to filter for
#   trades (specifically " Trade" with a leading space) and a 'Trade.Sign'
#   column.
#
# OUTPUT:
# 1.  A data frame 'df_acf' (ACF results) is created in the R environment.
# 2.  A 'ggplot' object 'p' is created and displayed in the RStudio plots pane.
# 3.  A summary of the fit (Estimated gamma, 95% CI, R-squared) is
#     printed to the console.
# ==============================================================================

# Required packages
library(ggplot2)
library(dplyr)

# === Load data ===
# Adjust filename/path if needed. We assume a CSV; change read.* if your file is a different format.
df <- read.csv("Data/JTOPI CSV files/NPNJn.J/NPNJn.J28-May-2015_TO_31-Dec-2015.csv", stringsAsFactors = FALSE)

ts_col <- "Trade.Sign"
if (length(ts_col) == 0) {
  # Fall back to searching for a column that only contains 0/1 values
  possible <- sapply(df, function(col) {
    is.numeric(col) && all(na.omit(unique(col)) %in% c(0,1))
  })
  if (any(possible)) {
    ts_col <- names(df)[which(possible)[1]]
    message("Using column '", ts_col, "' as trade-sign (inferred by values 0/1).")
  } else {
    stop("Couldn't find a 'trade sign' column. Rename it to 'Trade Sign' or 'Trade.Sign', or ensure it contains only 0/1.")
  }
} else {
  ts_col <- ts_col[1]
}

# === Subset to executed trades ===
# I assume 'Type' column exists and trades have Type == "Trade".
df_trades <- df |>  filter(Type == " Trade") # for some reason there is a blank space in the .csv

# === Prepare sign series: convert 0/1 -> -1/+1 and drop NAs ===
signs <- df_trades[[ts_col]]
signs <- signs[!is.na(signs)]
x <- ifelse(signs == 1, 1L, -1L)

Nevents <- length(x)
message("Number of trade events used: ", Nevents)

# ===  Compute autocorrelation (sample ACF) up to max_lag ===
max_lag_requested <- 3000L
max_lag <- min(max_lag_requested, Nevents - 1L)
if (max_lag < 10) stop("Too few trades to compute meaningful autocorrelation.")

# compute acf
acf_res <- acf(x, lag.max = max_lag, plot = FALSE, demean = TRUE)
# extract lag and rho correctly
lags <- acf_res$lag[,1,1]      # or as.vector(acf_res$lag)
rhos <- acf_res$acf[,1,1]     # or as.vector(acf_res$acf)

# build dataframe and drop lag 0
df_acf <- data.frame(lag = as.numeric(lags), rho = as.numeric(rhos)) %>%
  filter(lag > 0)

# drop lag 0
df_acf <- data.frame(lag = lags, rho = rhos) %>% filter(lag > 0)

# === Plot: log-log axes, ticks at powers of 10 similar to paper ===
x_breaks <- c(1, 10, 100, 1000)
x_breaks <- x_breaks[x_breaks <= max(df_acf$lag)]
y_breaks <- c(1e-3, 1e-2, 1e-1)
y_breaks <- y_breaks[y_breaks <= max(df_acf$rho, na.rm = TRUE) & y_breaks >= min(df_acf$rho, na.rm = TRUE)]
p <- ggplot(df_acf, aes(x = lag, y = rho)) +
  geom_line(size = 0.5) +
  scale_x_log10(breaks = x_breaks, labels = parse(text = paste0("10^", log10(x_breaks)))) +
  scale_y_log10(breaks = y_breaks, labels = parse(text = paste0("10^", round(log10(y_breaks))))) +
  labs(x = "lag (event)", y = "autocorrelation") +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    plot.margin = margin(8, 8, 8, 8)
  )

# show plot
print(p)
# === Save to file resembling paper figure size ===
#ggsave("Figures/FIG1_NPN.png", p, width = 5, height = 4, dpi = 300)

# === Least-squares fit on the log-log ACF to estimate gamma ===
# choose fitting window (adjustable)
lmin <- 10L
lmax <- min(1000L, max(df_acf$lag, na.rm = TRUE))
# prepare data for fit: keep only positive rho (log of non-positive undefined)
df_fit <- df_acf %>%
  filter(lag >= lmin, lag <= lmax, !is.na(rho), rho > 0)
# if the chosen window yields too few points, relax lmin
if (nrow(df_fit) < 10) {
  lmin <- 1L
  df_fit <- df_acf %>% filter(lag >= lmin, lag <= lmax, !is.na(rho), rho > 0)
}
if (nrow(df_fit) < 5) stop("Too few positive ACF points in the fitting window. Inspect df_acf and adjust lmin/lmax.")
# linear regression on log10 scale: log10(rho) = intercept + slope * log10(lag)
fit <- lm(log10(rho) ~ log10(lag), data = df_fit)
# extract gamma (the power-law exponent): rho ~ lag^{-gamma} => slope = -gamma
slope <- coef(fit)["log10(lag)"]
gamma_hat <- -as.numeric(slope)
# 95% CI for slope -> translate to gamma CI
slope_ci <- confint(fit, "log10(lag)", level = 0.95)
gamma_ci <- c(-slope_ci[2], -slope_ci[1])  # invert and flip
# summary stats
fit_summary <- summary(fit)
r_squared <- fit_summary$r.squared
# print results
cat(sprintf("Fitting lag range: %d to %d (using %d points)\n", lmin, lmax, nrow(df_fit)))
cat(sprintf("Estimated gamma = %.5f\n", gamma_hat))
cat(sprintf("95%% CI for gamma = [%.5f, %.5f]\n", gamma_ci[1], gamma_ci[2]))
cat(sprintf("Slope (log10 scale) = %.5f, Intercept = %.5f\n", as.numeric(slope), as.numeric(coef(fit)["(Intercept)"])))
cat(sprintf("R-squared = %.4f\n\n", r_squared))
print(fit_summary)