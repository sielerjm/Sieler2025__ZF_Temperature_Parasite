# Microbiome Data Analysis Processing -------------------------------------
# 
# This script performs the following main tasks:
# 1. Imports and processes microbiome data from a phyloseq object
# 2. Calculates alpha and beta diversity metrics
# 3. Performs statistical analyses and generates plots
# 4. Saves results and environment for downstream analysis

## Error Handling Setup ---------------------------------------------------

# Custom error reporting function
create_error_report <- function(error, step) {
  # Ensure error directory exists
  if (!dir.exists(file.path("Data", "Errors"))) {
    dir.create(file.path("Data", "Errors"), recursive = TRUE)
  }
  
  # Create detailed error message
  error_msg <- paste0(
    "Error in MicrobiomeProcessing_2025-07-01.R\n",
    "Step: ", step, "\n",
    "Time: ", Sys.time(), "\n",
    "Error Message: ", error$message, "\n",
    "Call: ", deparse(error$call), "\n"
  )
  
  # Save error report with timestamp
  error_file <- file.path("Data", "Errors", 
                         paste0("error_report_", 
                                format(Sys.time(), "%Y%m%d_%H%M%S"), 
                                ".txt"))
  writeLines(error_msg, error_file)
  stop(error_msg)
}

## Main Script Execution --------------------------------------------------

tryCatch({
  # Start timing and logging
  start_time <- Sys.time()
  cat("Starting MicrobiomeProcessing_2025-07-01.R at ", 
      format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
  
  ## Initialize Environment ------------------------------------------------
  cat("INITIALIZING ENVIRONMENT \n")
  
  # Create unique analysis ID
  analysis.ID <- paste0(
    "Sieler2024__ZF_Temperature_Parasite__",  
    Sys.Date()  # Date of last processing
  )
  
  # Define diversity calculation methods
  diversity.method <- list(
    alpha = c("shannon", "inverse_simpson", "observed", "phylogenetic"),
    beta = c("bray", "canberra", "gunifrac")
  )
  
  # Initialize storage for results
  alpha.stats <- list()
  alpha.plots <- list()
  beta.stats <- list()
  beta.plots <- list()
  diffAbnd.stats <- list()
  diffAbnd.plots <- list()
  worm.stats <- list()
  worm.plots <- list()
  randomForest.stats <- list()
  randomForest.plots <- list()
  
  ## Data Import and Processing --------------------------------------------
  
  # Import and subset phyloseq object
  tryCatch({
    cat("CREATING PHYLOSEQ OBJECTS \n")
    ps.list <- list()
    ps.list[["All"]] <- readRDS(
      file.path(path.data, "R_objects", 
               paste0("phyloseq__PostDADA2Cleaning_2024-07-08.rds"))
    )
    
    cat("SUBSETTING PHYLOSEQ OBJECTS \n")
    # Create various subsets based on experimental conditions
    ps.list[["Unexposed"]] <- ps.list[["All"]] %>% 
      ps_filter((Treatment == "Control") | (Treatment == "Exposed" & DPE == 0))
    
    ps.list[["Exposed"]] <- ps.list[["All"]] %>% 
      ps_filter((Treatment == "Exposed") & DPE != 0)
    
    ps.list[["PreExposed"]] <- ps.list[["All"]] %>% 
      ps_filter(DPE == 0)
    
    ps.list[["PostExposed"]] <- ps.list[["All"]] %>% 
      ps_filter(DPE != 0)
    
    ps.list[["TimeFinal"]] <- ps.list[["All"]] %>% 
      ps_filter(DPE == 42)
    
    ps.list[["InitialFinal"]] <- ps.list[["All"]] %>% 
      ps_filter(DPE == 0 | DPE == 42)
  }, error = function(e) create_error_report(e, "Data Import and Processing"))
  
  ## Diversity Analysis ---------------------------------------------------
  
  tryCatch({
    cat("CALCULATING DIVERSITY SCORES \n")
    
    # Calculate alpha diversity metrics
    cat(" - Calculating Alpha Diversity\n")
    source(paste0(proj.path,"/Code/Functions/AnalysisScripts/AlphaDiversity.R"))
    
    # Calculate beta diversity metrics
    cat(" - Calculating Beta Diversity\n")
    source(paste0(proj.path,"/Code/Functions/AnalysisScripts/BetaDiversity.R"))
  }, error = function(e) create_error_report(e, "Diversity Analysis"))
  
  ## Statistical Analysis and Plotting ------------------------------------
  
  tryCatch({
    cat("RUNNING STATISTICAL ANALYSES AND PLOTTING \n")
    
    # Run statistical analyses
    cat(" - Running Statistical Scripts\n")
    sourceFolder(paste0(proj.path,"/Code/Analysis"), T)
    
    # Generate plots
    cat(" - Running Plot Scripts\n")
    sourceFolder(paste0(proj.path,"/Code/Plots"), T)
  }, error = function(e) create_error_report(e, "Statistical Analysis and Plotting"))
  
  ## Save Results ---------------------------------------------------------
  
  tryCatch({
    cat("SAVING RESULTS \n")
    
    # Save statistical results
    cat(" - Saving Statistical Objects\n")
    saveRDS(alpha.stats, file.path(path.results, "Stats_Plots", "alpha_stats.rds"))
    saveRDS(beta.stats, file.path(path.results, "Stats_Plots", "beta_stats.rds"))
    saveRDS(worm.stats, file.path(path.results, "Stats_Plots", "worm_stats.rds"))
    saveRDS(diffAbnd.stats, file.path(path.results, "Stats_Plots", "diffAbnd_stats.rds"))
    saveRDS(randomForest.stats, file.path(path.results, "Stats_Plots", "randomForest_stats.rds"))
    
    # Save plot objects
    cat(" - Saving Plot Objects\n")
    saveRDS(alpha.plots, file.path(path.results, "Stats_Plots", "alpha_plots.rds"))
    saveRDS(beta.plots, file.path(path.results, "Stats_Plots", "beta_plots.rds"))
    saveRDS(worm.plots, file.path(path.results, "Stats_Plots", "worm_plots.rds"))
    saveRDS(diffAbnd.plots, file.path(path.results, "Stats_Plots", "diffAbnd_plots.rds"))
    saveRDS(randomForest.plots, file.path(path.results, "Stats_Plots", "randomForest_plots.rds"))
    
    # Save environment
    cat(" - Saving Environment\n")
    sesh.info <- sessionInfo()
    saveRDS(sesh.info, file.path(path.results, "Stats_Plots", "sessionInfo.rds"))
    save_env(path.objects, ID = analysis.ID, extra_info = "microbiomeProcessing")
  }, error = function(e) create_error_report(e, "Saving Results"))
  
  ## Completion Message ---------------------------------------------------
  
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "secs")
  cat("Completed MicrobiomeProcessing_2025-07-01.R at ", 
      format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Total execution time: ", round(duration, 2), " seconds\n")
  
}, error = function(e) {
  create_error_report(e, "Main Script Execution")
  stop("Script failed. See error_report.txt for details.")
})

