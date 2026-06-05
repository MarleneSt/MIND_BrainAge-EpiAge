compute_epigenetic_ages <- function(samplesheet_path,
                                    methylation_rds_path,
                                    output_path,
                                    clocks = c("PCGrimAge", "ShirebyG2020", "BernabeuE2023c"),
                                    edited_methyAge_path = NULL,
                                    array) {
  # Load data
  samplesheet <- read.csv(samplesheet_path)
  MethylationData <- readRDS(methylation_rds_path)
  
  message("Methylation data dimensions: ")
  print(dim(MethylationData))
  
  # Prepare sample info
  info <- samplesheet[, c('Sample_Name', 'age', 'sex')] # 'cidB3067.qlet
  colnames(info) <- c('Sample','Age', 'Sex') # 'cidB3067.qlet'
  
  # Optionally source edited methyAge
  if ("BernabeuE2023c" %in% clocks) {
    if (is.null(edited_methyAge_path)) {
      stop("Provide the path to edited methyAge function for BernabeuE2023c.")
    } else {
      source(edited_methyAge_path)
    }
  }
  
  # Initialize results list
  results_list <- list()
  
  # Loop over clocks
  for (clock in clocks) {
    message("Running clock: ", clock)
    if (clock == "BernabeuE2023c") {
      age_result <- methyAge_edited(MethylationData, clock = clock, age_info = info)
    } else {
      age_result <- methyAge(MethylationData, clock = clock, age_info = info)
    }
    
    # Rename columns
    age_result[[paste0(clock, "_mAge")]] <- age_result$mAge
    
    # Remove original columns before merging
    results_list[[clock]] <- age_result[, c("Sample", paste0(clock, "_mAge"))]
  }
  
  # Merge all results by Sample
  combined_ages <- Reduce(function(x, y) merge(x, y, by = "Sample"), results_list)
  
  return(combined_ages)
}