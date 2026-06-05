#### Saving ENIGMA output files with subject IDs ----
#### written by Marlene Staginnus 

cat("loading brain age output files\n\n")
female <- read.table("females_raw_out.csv", header = TRUE, sep="\t")
male <- read.csv("males_raw_out.csv", header = TRUE, sep="\t")

cat("reading brain age input files with subject IDs\n\n")
female_IDs <- read.csv("females_raw_withIDs.csv")
male_IDs <- read.csv("males_raw_withIDs.csv")

# checking 
identical(female$ICV, female_IDs$ICV)
identical(male$ICV, male_IDs$ICV)

cat("adding subject IDs into output files\n\n")
# add SUBJID back in 
female$SUBJID <- female_IDs$SUBJID
male$SUBJID <- male_IDs$SUBJID

#reorder & remove x column 
female <- female[, c("SUBJID", setdiff(names(female), "SUBJID"))]
female$X <- NULL

male <- male[, c("SUBJID", setdiff(names(male), "SUBJID"))]
male$X <- NULL

# also make file with just age prediction 
female_reduced <- female[, c("SUBJID", "age_prediction")]
male_reduced <- male[, c("SUBJID", "age_prediction")]

cat("saving output files with subject IDs:\n- female_raw_out_withIDs.csv\n- male_raw_out_withIDs.csv\n\n")

# save files 
write.csv(female, "female_raw_out_withIDs.csv", row.names = FALSE)
write.csv(male, "male_raw_out_withIDs.csv", row.names = FALSE)

cat("saving files with only ID + brain age estimate:\n- female_ENIGMA_output.csv\n- male_ENIGMA_output.csv\n- all_ENIGMA_output.csv\n\n")

write.csv(female_reduced, "female_ENIGMA_output.csv", row.names = FALSE)
write.csv(male_reduced, "male_ENIGMA_output.csv", row.names = FALSE)

both_sex_reduced <- rbind(female_reduced, male_reduced)
write.csv(both_sex_reduced, "all_ENIGMA_output.csv", row.names = FALSE)


