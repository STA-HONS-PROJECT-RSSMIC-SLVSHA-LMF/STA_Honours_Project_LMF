# ==============================================================================
# SCRIPT: 06_LMF_Validation_load_and_prep_data
#
# PURPOSE:
# This script is responsible for loading the complete, consolidated trade
# dataset ("JSE_Top40_new.csv"). It performs essential type conversions
# (like for timestamps) and minor cleaning (like removing old,
# conflicting columns from previous analysis runs).
#
# It creates the master data.table 'DT' that all subsequent analysis
# scripts will use. It also saves this 'DT' object into a workspace
# file ("JSE_Top40_workspace.RData") for much faster loading in the future.
#
# PREREQUISITES:
# * Must be run AFTER `01_Data_Consolidation_and_Cleaning.R` (which
#   is the script that *creates* "JSE_Top40_new.csv").
# * Packages: 'data.table'.
#
# OUTPUT:
# * 'DT' (data.table): The master data table is loaded into the
#   R environment.
# * "JSE_Top40_workspace.RData": A file is saved to disk containing
#   the 'DT' object to speed up future sessions.
# ==============================================================================

library(data.table)

workspace_file <- "JSE_Top40_workspace.RData" # Define the workspace file
if (file.exists(workspace_file)) {
  # Load saved workspace if it exists
  load(workspace_file)
  message("Workspace loaded from ", workspace_file)
} else {
  # ===== Load and preprocess the CSV ===== #
  library(data.table)
  DT <- read.csv("JSE_Top40_new.csv")
  setDT(DT)
  # Convert types
  DT[, timestamp := as.POSIXct(timestamp, format="%Y-%m-%dT%H:%M:%OS", tz = "Africa/Johannesburg")]
  DT[, Sign := as.integer(Sign)]
  DT <- DT[Sign %in% c(-1L, 1L)]
  DT[, trade_date := as.IDate(timestamp, tz = "Africa/Johannesburg")]
  DT[, month_year := format(timestamp, "%Y-%m")] # Creates a 'YYYY-MM' column
}
# === REMOVE ANY OLD INFORMATION THAT HAS BEEN SAVED ===
# Check if DT exists and is a data.table
if (exists("DT") && is.data.table(DT)) {
  cols_to_remove <- c("VD", "p0", "pmax", "pmin", "sigmaD", 
                      "trades_in_day", "trader_id", "metaorder_id")
  # Find which of those columns are actually in DT
  cols_present <- intersect(cols_to_remove, names(DT))
  # Remove them if they exist
  if (length(cols_present) > 0) {
    DT[, (cols_present) := NULL]
    message("DT reset")
  }
}
save.image(file = workspace_file)