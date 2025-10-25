# STA Honours Project: LMF Model Validation and Long-Range Memory Analysis
This repository stores the R scripts used by Michael Ross and Shane Silverman for their STA_Honours Project on order-flow and long memory in a simulated financial market. The scripts analyse long-range memory in high-frequency trade data from the Johannesburg Stock Exchange (JSE), specifically focusing on validating the Lillo, Mike, and Farmer (LMF) hypothesis. 

The central goal is to test the relationship ($\gamma \approx \alpha - 1$) between the power-law exponent ($\gamma$) of the trade-sign autocorrelation function (ACF), and the power-law exponent ($\alpha$) of the distribution of "metaorder" sizes.

The project is divided into three main parts:
1. Data Preparation: Consolidating and cleaning raw JSE Top40 trade data.
2. Theoretical Model Replication: Simulating the models mentioned in the LMF literature.
3. Hypothesis Validation: Applying a metaorder generation algorithm to the real JSE data to estimate the power-law exponent of the distribution of "metaorder" sizes and compare it against the empirically measured power-law exponent of the trade-ACF.

## Key Concepts and Functions
The core analysis relies on several custom functions defined in 05_LMF_Validation_analysis_functions.R:
1. estimate_gamma_acf(): Estimates $\gamma$ by fitting a power law to the log-binned Autocorrelation Function (ACF) of the trade-sign series.
2. estimate_gamma_psd(): Estimates $\gamma$ by fitting a power law to the low-frequency end of the Power Spectral Density (PSD) of the trade-sign series.
3. estimate_alpha_from_vec(): A wrapper for the poweRlaw package to estimate the $\alpha$ exponent from a vector of (metaorder) sizes using the Clauset method.
4. mapping_function_orig(): The "Original MLB" simulation model, which assigns trades to traders based on a power-law sample, then shuffles all trades.
5. mapping_function_adapt(): The "Adaptive MLB" simulation model, which creates contiguous blocks of trades for each trader before shuffling the blocks, better preserving metaorder structure.

## Prerequisites
Install the required packages by running the following in your R console:
`install.packages(c(
  "data.table", 
  "ggplot2", 
  "dplyr", 
  "scales", 
  "parallel", 
  "poweRlaw", 
  "cowplot",
  "pbapply"
))`

## Recommended Workflow
### Step 1: Data Preparation
First, you must process your raw data into the master workspace file.
1. `01_Data_Consolidation_and_Cleaning.R`
   
   What it does: Scans all ticker subfolders, reads every .csv file, cleans the data, and consolidates everything into a single, large CSV.
   
   Input: Raw .csv files in ticker subfolders.
   
   Output: JSE_Top40_new.csv
   
2. `06_LMF_Validation_load_and_prep_data.R`
   
   What it does: Loads the massive JSE_Top40_new.csv into a data.table (DT), performs final type conversions, and saves the object as an .RData file for much faster loading in the future.

   Input: JSE_Top40_new.csv (from script 01).

   Output: JSE_Top40_workspace.RData
### Step 2: Main Hypothesis Validation
These scripts test the $\gamma \approx \alpha - 1$ hypothesis on the prepared data. They must be run after Step 1.

0. `05_LMF_Validation_analysis_functions.R`

      This is a library. You do not run this script directly. It is automatically loaded (source()-ed) by scripts 07 and 08.

1. `07_LMF_Validation_analysis_full_year.R`
   
      What it does: Performs the main analysis on a per-stock, full-year basis.

      Process:
   
        a. Loads JSE_Top40_workspace.RData.
     
        b. Calculates one $\gamma$ value for each stock.
     
        c. Runs both the "Original MLB" and "Adaptive MLB" simulations on the data.
     
        d. Calculates one $\alpha$ value for each stock (for each model).
     
        e. Generates plots comparing $\gamma$ vs. $\alpha - 1$.

      Input: JSE_Top40_workspace.RData (from script 06) and functions (from script 05).
      
      Output: Plots displayed in the R plots pane.

2. `08_LMF_Validation_analysis_monthly.R`
   
      What it does: Performs a more granular version of the analysis on a per-stock, per-month basis. This generates many more data points.
      
      Process: Same as script 07, but all calculations are grouped by both Stock and month_year.
      
      Input: JSE_Top40_workspace.RData (from script 06) and functions (from script 05).
      
      Output: Plots displayed in the R plots pane.

## Standalone Scripts
These scripts do not depend on the main data preparation workflow (Steps 1 & 2) and can be run independently.
### Theoretical Model Replications
These scripts replicate the theoretical models from the LMF literature. They are self-contained, long-running simulations that use caching (saving results to an .rds file) to avoid re-running.

`03_LMF_Fixed_N_Model_Simulation.R`

      What it does: Simulates the "Fixed-N" model, demonstrating how the ACF scales with the number of agents (N).
      
      Output: Caches results in a .rds file and generates a log-log plot.

`04_LMF_Lambda_Model_Simulation.R`

      What it does: Simulates the "$\lambda$-model", investigating how the ACF of active orders changes with the order arrival rate ($\lambda$).
      
      Output: Caches results in a .rds file and generates a log-log plot.

### Single-Stock ACF
`02_ACF_and_Gamma_Analysis.R`

      What it does: A simple, standalone script to analyze a single, hardcoded CSV file (Naspers, "NPNJn.J"). It loads the file, plots its ACF on a log-log scale, and estimates $\gamma$.
      
      Use: Good for a quick, preliminary check of auto-correlation in a specific stock.
