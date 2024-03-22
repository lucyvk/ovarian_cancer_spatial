library(tidyverse)
library(survminer)
source("survival_helper.R")

# main - main text results
# missing na - na values not replaced with 0 (see SI note)
# primary - only primary tumor samples
# main_th100 - sensitivity test with higher spatial network trimming threshold
data_versions <- c('main','missing_na','primary','main_th100')

for (dv in data_versions) {
  
  mibi_table <- as.data.frame(read.csv(sprintf('Data/mibi_table_%s.csv',dv), header = TRUE, sep = ",", dec = ".", check.names=FALSE))

  # Feature subsets
  clinical_features <- select(mibi_table, Age, H3K14Ace, ATF6, DUSP1, CBX2, BRCA_Mutation)
  # Turn the BRCA mutation column into a factor (categorical)
  clinical_features <- mutate(clinical_features,"BRCA_Mutation" = as.factor(clinical_features$BRCA_Mutation))
  composition_features <- select(mibi_table, matches("*_prop"))
  spatial_features <- select(mibi_table, matches("*_dist"))
  network_features <- cbind(select(mibi_table, matches("*_contact")),select(mibi_table, matches("*_assort")),select(mibi_table, matches("*_region")))
  
  # Change for matching to feature category names in plotting
  colnames(clinical_features) <- gsub(" ","_",colnames(clinical_features), fixed=TRUE)
  colnames(composition_features) <- gsub(" ","_",colnames(composition_features), fixed=TRUE)
  colnames(spatial_features) <- gsub(" ","_",colnames(spatial_features), fixed=TRUE)
  colnames(network_features) <- gsub(" ","_",colnames(network_features), fixed=TRUE)
  # temporarily replace + and - and / and spaces for use in the Surv function
  colnames(mibi_table) <- gsub("+","plus",colnames(mibi_table), fixed=TRUE)
  colnames(mibi_table) <- gsub("-","minus",colnames(mibi_table), fixed=TRUE)
  colnames(mibi_table) <- gsub(" ","_",colnames(mibi_table), fixed=TRUE)
  colnames(mibi_table) <- gsub("/","slash",colnames(mibi_table), fixed=TRUE)
  
  covariates <- colnames(mibi_table)
  covariates <- covariates[-which(covariates %in% c('fov_id','TMA_ID','OS','PFS','OS_high','PFS_high','Death','Recurrence','Primary'))]
  
  # Scale first so we can compare hazard ratios between features on different scales
  covariates_df <- mibi_table[,covariates[-which(covariates %in% c('BRCA_Mutation', 'Age'))]]
  covariates_scaled <- as.data.frame(scale(covariates_df))
  BRCA_Mutation <- as.factor(mibi_table$BRCA_Mutation) # Treat as a factor (not scaled)
  Age <- mibi_table$Age
  
  # Remaining NA values handled by excluding from the regressions 
  
  ### OVERALL SURVIVAL 
  
  time <- mibi_table$OS
  status <- mibi_table$Death 
  df_surv <- cbind(time,status,BRCA_Mutation,Age,covariates_scaled)
  
  # Survival Curve
  filename <- sprintf("Figure_5/kaplan_meier_survival_%s.jpeg",dv)
  survival_curve(time, status, df_surv, filename)

  # Run Cox Regression
  survival_res <- run_cox(time, status, df_surv)
  filename <- sprintf("Figure_5/significant_covariates_%s_%s.csv","Survival",dv)
  # Change names back for printing to file
  rownames(survival_res) <- gsub("plus","+",rownames(survival_res), fixed=TRUE)
  rownames(survival_res) <- gsub("minus","-",rownames(survival_res), fixed=TRUE)
  rownames(survival_res) <- gsub("slash","/",rownames(survival_res), fixed=TRUE)
  write_significant(survival_res, filename) # write the significant ones to file
  write.csv(survival_res,sprintf("Figure_5/cox_reg_survival_%s.csv",dv)) # write the full results to file
  
  # Plotting significant covariates
  filename <- sprintf("Figure_5/top_covariates_%s_%s.jpeg","Survival",dv)
  plot_significant(survival_res,"Survival",FALSE,filename)
  dev.off()
  
  ### PROGRESSION-FREE SURVIVAL
  
  time <- mibi_table$PFS
  status <- mibi_table$Recurrence
  df_surv <- cbind(time,status,BRCA_Mutation,Age,covariates_scaled)
  
  # Survival Curve
  filename <- sprintf("Figure_5/kaplan_meier_recurrence_%s.jpeg",dv)
  survival_curve(time, status, df_surv, filename)

  # Run Cox Regression
  recurrence_res <- run_cox(time, status, df_surv)
  filename <- sprintf("Figure_5/significant_covariates_%s_%s.csv","Recurrence",dv)
  
  # Change names back for printing to file
  rownames(recurrence_res) <- gsub("plus","+",rownames(recurrence_res), fixed=TRUE)
  rownames(recurrence_res) <- gsub("minus","-",rownames(recurrence_res), fixed=TRUE)
  rownames(recurrence_res) <- gsub("slash","/",rownames(recurrence_res), fixed=TRUE)
  write_significant(recurrence_res, filename)
  write.csv(recurrence_res,sprintf("Figure_5/cox_reg_recurrence_%s.csv",dv))
  
  # Plotting significant covariates
  filename <- sprintf("Figure_5/top_covariates_%s_%s.jpeg","Recurrence",dv)
  plot_significant(recurrence_res,"Recurrence",FALSE,filename)
  dev.off()

}
