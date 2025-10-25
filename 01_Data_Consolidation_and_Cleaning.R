# ==============================================================================
# SCRIPT: 01_Data_Consolidation_and_Cleaning
#
# PURPOSE:
# This script consolidates raw trade data from numerous individual .csv files,
# which are organized in ticker-specific subfolders, into a single, large,
# and cleaned CSV file. It is designed to be memory-efficient by processing
# each file one by one (chunk by chunk) rather than loading all data at once.
# 
# PREREQUISITES:
# * Package: 'data.table' must be installed.
# * Directory Structure: This script must be run from a root directory
#     that contains one subfolder for each stock ticker (e.g., ./AGL, ./BTI).
# * Input Data: The raw .csv files are expected to have columns named
#     'Price', 'Volume', 'Trade Sign', and a timestamp (either 'DateTime' or 'DateTimeL').
#
# OUTPUT:
# * A single file named "JSE_Top40_new.csv" in the root directory.
# * Note: Any existing file with this name will be deleted before the script runs.
# ==============================================================================

library(data.table)
# Root folder with one subfolder per ticker, each containing CSVs
ticker_dirs <- list.dirs(".", full.names = TRUE, recursive = FALSE)
# Output file
out_file <- "JSE_Top40_new.csv"
if (file.exists(out_file)) file.remove(out_file)
# Initialize row counter
total_rows <- 0
# Process each CSV one by one
for (folder in ticker_dirs) {
  csvs <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
  for (file in csvs) {
    dt <- fread(file, showProgress = FALSE)
    # Build POSIX timestamp immediately
    if ("DateTime" %in% names(dt) && !inherits(dt$DateTime, "POSIXct")) {
      dt[, DateTime := as.POSIXct(DateTime, tz = "Africa/Johannesburg")]
    } else if (!("DateTime" %in% names(dt)) && "DateTimeL" %in% names(dt)) {
      offset_days <- if (max(dt$DateTimeL, na.rm = TRUE) > 700000) 719529 else 25569
      dt[, DateTime := as.POSIXct((DateTimeL - offset_days) * 86400,
                                  origin = "1970-01-01", tz = "Africa/Johannesburg")]
    }
    # Add Ticker column
    dt[, Stock := basename(folder)]
    # Update counter
    total_rows <- total_rows + nrow(dt)
    # Keep only non-zero trades
    if ("Trade Sign" %in% names(dt)) dt <- dt[dt$`Trade Sign` != 0]
    # Select & rename relevant columns immediately
    dt <- dt[, .(
      Stock     = Stock, # or RIC if preferred
      timestamp = DateTime,
      price     = Price,
      volume    = Volume,
      Sign      = `Trade Sign`
    )]
    # Basic cleaning
    dt <- dt[is.finite(price) & price > 0 & is.finite(volume) & volume > 0]
    # Sort just the chunk (optional, can sort later on full CSV)
    setorder(dt, Stock, timestamp)
    # Append to CSV
    fwrite(dt, out_file, append = file.exists(out_file))
    # Free memory
    rm(dt); gc()
  }
}
cat("Processing complete! Output written to", out_file, "\n")
cat("Total number of rows written:", total_rows, "\n")