#########################################################
############# SCRIPT TO IDENTIFY OUTLIERS ###############
############# based on MRIQC group report ###############
#########################################################
## written by Marlene Staginnus 

# Load necessary library
list.of.packages <- c("data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(data.table)

# Load the MRIQC data
ts_file <- "group_T1w.tsv"  # Replace with your actual file path
data <- fread(ts_file)

# Define function to detect outliers
find_outliers <- function(data, metric, direction) {
  Q1 <- quantile(data[[metric]], 0.25, na.rm = TRUE)
  Q3 <- quantile(data[[metric]], 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  
  lower_fence <- Q1 - 1.5 * IQR_val
  upper_fence <- Q3 + 1.5 * IQR_val
  
  if (direction == "high") {
    return(ifelse(data[[metric]] > upper_fence, "outlier_high", ""))
  } else if (direction == "low") {
    return(ifelse(data[[metric]] < lower_fence, "outlier_low", ""))
  } else if (direction == "both") {
    return(ifelse(data[[metric]] > upper_fence, "outlier_high", 
                  ifelse(data[[metric]] < lower_fence, "outlier_low", "")))
  }
}


# Define the metrics and their outlier directions
metrics <- list(
  "cjv" = "high",
  "cnr" = "low",
  "efc" = "high",
  "fber" = "low",
  "fwhm_avg" = "high",
  "fwhm_x" = "high",
  "fwhm_y" = "high",
  "fwhm_z" = "high",
  "qi_1" = "high",
  "qi_2" = "high",
  "rpve_csf" = "high",
  "rpve_gm" = "high",
  "rpve_wm" = "high",
  "snr_csf" = "low",
  "snr_gm" = "low",
  "snr_total" = "low",
  "snr_wm" = "low",
  "snrd_csf" = "low",
  "snrd_gm" = "low",
  "snrd_total" = "low",
  "snrd_wm" = "low",
  "summary_bg_k" = "both",
  "summary_bg_mad" = "both",
  "summary_bg_mean" = "both",
  "summary_bg_median" = "both",
  "summary_bg_n" = "both",
  "summary_bg_p05" = "both",
  "summary_bg_p95" = "both",
  "summary_bg_stdv" = "both",
  "summary_csf_k" = "both",
  "summary_csf_mad" = "both",
  "summary_csf_mean" = "both",
  "summary_csf_median" = "both",
  "summary_csf_n" = "both",
  "summary_csf_p05" = "both",
  "summary_csf_p95" = "both",
  "summary_csf_stdv" = "both",
  "summary_gm_k" = "both",
  "summary_gm_mad" = "both",
  "summary_gm_mean" = "both",
  "summary_gm_median" = "both",
  "summary_gm_n" = "both",
  "summary_gm_p05" = "both",
  "summary_gm_p95" = "both",
  "summary_gm_stdv" = "both",
  "summary_wm_k" = "both",
  "summary_wm_mad" = "both",
  "summary_wm_mean" = "both",
  "summary_wm_median" = "both",
  "summary_wm_n" = "both",
  "summary_wm_p05" = "both",
  "summary_wm_p95" = "both",
  "summary_wm_stdv" = "both",
  "tpm_overlap_csf" = "low",
  "tpm_overlap_gm" = "low",
  "tpm_overlap_wm" = "low"
)

# Initialize output data frame
output_data <- data[, .(bids_name)]

# Apply outlier detection
for (metric in names(metrics)) {
  if (metric %in% colnames(data)) {
    output_data[[metric]] <- find_outliers(data, metric, metrics[[metric]])
  }
}

# Save the results
write.csv(output_data, "outliers_MRIQC.csv", row.names = FALSE)

print("Outlier detection complete. Results saved in outliers_MRIQC.csv")
