library(randomForest)
library(pROC)
library(caTools)
library(mltools)
library(data.table)
library(tidyverse)
library(gsubfn)

seed <- 19 # For reproducing the same results as in the paper 
set.seed(seed) 

# main - main text results
# missing na - na values not replaced with 0 (see SI note)
# primary - only primary tumor samples
# main_th100 - sensitivity test with higher spatial network trimming threshold
data_versions <- c('missing_na','main')
#data_versions <- c('primary', 'main_th100')
reps <- 500

# Run the repeated prediction evaluation for versions of the data 
for (dv in data_versions) {

  # Load data table 
  mibi_table <- as.data.frame(read.csv(sprintf('Data/mibi_table_%s.csv',dv), header = TRUE, sep = ",", dec = ".", check.names=FALSE))
  
  # Feature subsets
  clinical_features <- select(mibi_table, Age, H3K14Ace, ATF6, DUSP1, CBX2, BRCA_Mutation)
  # Turn the BRCA mutation column into a factor (categorical)
  clinical_features <- mutate(clinical_features,"BRCA_Mutation" = as.factor(clinical_features$BRCA_Mutation))
  composition_features <- select(mibi_table, matches("*_prop"))
  spatial_features <- select(mibi_table, matches("*_dist"))
  network_features <- cbind(select(mibi_table, matches("*_contact")),select(mibi_table, matches("*_assort")),select(mibi_table, matches("*_region")))
  
  m1_features <- clinical_features
  m2_features <- composition_features
  m3_features <- spatial_features
  m4_features <- network_features
  m5_features <- cbind(clinical_features,composition_features)
  m6_features <- cbind(clinical_features,spatial_features)
  m7_features <- cbind(clinical_features,network_features)
  m8_features <- cbind(composition_features,spatial_features)
  m9_features <- cbind(composition_features,network_features)
  m10_features <- cbind(spatial_features,network_features)
  m11_features <- cbind(clinical_features,composition_features,spatial_features)
  m12_features <- cbind(clinical_features,composition_features,network_features)
  m13_features <- cbind(clinical_features,spatial_features,network_features)
  m14_features <- cbind(spatial_features,composition_features,network_features)
  m15_features <- cbind(clinical_features,composition_features,spatial_features, network_features)
  
  feature_set_sizes <- c(ncol(m1_features),ncol(m2_features),ncol(m3_features),ncol(m4_features),ncol(m5_features),ncol(m6_features),ncol(m7_features),ncol(m8_features),ncol(m9_features),ncol(m10_features),ncol(m11_features),ncol(m12_features),ncol(m13_features),ncol(m14_features),ncol(m15_features))
  
  # Random Forest model training and evaluation
  run_model <- function(feature_subset, Y, return_imp) {
    
    # Split into train and test set, balancing the outcome label
    # If a column value is all NA in the test set, redraw (rare)
    redraw <- 1
    while (redraw == 1) {
      msk <- sample.split(Y, 0.7)
      train_x <- feature_subset[msk,]
      test_x <- feature_subset[!msk,]
      train_y <- as.factor(Y[msk])
      test_y <- Y[!msk]
      
      # check test set
      all_NA_test <- sapply(test_x, function(x) all(is.na(x)))
      if (length(which(all_NA_test==TRUE)) > 0){
        print("all NA columns, redrawing")
        for (i in 1:length(all_NA_test)) {
          if (all_NA_test[i]==TRUE){
            print(names(all_NA_test)[i])
          }
        }
      } else{
        redraw <- 0 
      }
  
    }
    
    # replace NA values 
    train_x <- na.roughfix(train_x)
    # Supply outcome as a factor to force classification instead of regression
    classifier <- randomForest(x = train_x, y = train_y)
    # Predict probabilities for the test dataset
    pred_y <- predict(classifier, newdata = na.roughfix(test_x), type = "prob")
    if (return_imp) {
      # randomForest package - type = 2 returns the mean decrease in node impurity (Gini index)
      importance_res <- importance(classifier, type=2)
    }
    
    # ROC / AUC
    # "direction" forces the direction of the curve / says that positive (1) values should have higher predicted values
    # Use just the second column from the prediction (probability of being a 1)
    prob1 <- pred_y[,2]
    curr_roc <- roc(response = test_y, predictor = prob1, auc=TRUE, direction="<", quiet = TRUE)
    curr_auc <- auc(curr_roc)
    
    if (return_imp) {
      return(list(curr_auc,t(importance_res)))
      
    } else {
      return(list(curr_auc))
    }
  
    
  }
  
  # Run training and evaluation a bunch of times to see average performance / feature importance
  run_a_bunch <- function(Y,N,name) {
    m1s <- numeric(N)
    m2s <- numeric(N)
    m3s <- numeric(N)
    m4s <- numeric(N)
    m5s <- numeric(N)
    m6s <- numeric(N)
    m7s <- numeric(N)
    m8s <- numeric(N)
    m9s <- numeric(N)
    m10s <- numeric(N)
    m11s <- numeric(N)
    m12s <- numeric(N)
    m13s <- numeric(N)
    m14s <- numeric(N)
    m15s <- numeric(N)
    # feature importance results from the full model 
    m15_imp <- matrix(ncol = length(colnames(m15_features)), nrow = N)
    for (x in 1:N) {
      print(sprintf("running rep %d",x))
      list[m1s[x]] <- run_model(m1_features, Y, FALSE)
      list[m2s[x]] <- run_model(m2_features, Y, FALSE)
      list[m3s[x]] <- run_model(m3_features, Y, FALSE)
      list[m4s[x]] <- run_model(m4_features, Y, FALSE)
      list[m5s[x]] <- run_model(m5_features, Y, FALSE)
      list[m6s[x]] <- run_model(m6_features, Y, FALSE)
      list[m7s[x]] <- run_model(m7_features, Y, FALSE)
      list[m8s[x]] <- run_model(m8_features, Y, FALSE)
      list[m9s[x]] <- run_model(m9_features, Y, FALSE)
      list[m10s[x]] <- run_model(m10_features, Y, FALSE)
      list[m11s[x]] <- run_model(m11_features, Y, FALSE)
      list[m12s[x]] <- run_model(m12_features, Y, FALSE)
      list[m13s[x]] <- run_model(m13_features, Y, FALSE)
      list[m14s[x]] <- run_model(m14_features, Y, FALSE)
      list[m15s[x], m15_imp[x,]] <- run_model(m15_features, Y, TRUE)
    }
    
    # summary data frame 
    pred_res <- data.frame(model = 1, mean = mean(m1s), min = min(m1s), max = max(m1s), quant25 = quantile(m1s,0.25), median = median(m1s), quant75 = quantile(m1s, 0.75), sd = sd(m1s), var= var(m1s))
    
    add_row <- function(initial_df, model_num, out){
      new_df <- rbind(initial_df, data.frame(model = model_num, mean = mean(out), min = min(out), max = max(out), quant25 = quantile(out,0.25), median = median(out), quant75 = quantile(out, 0.75), sd = sd(out), var= var(out)))
      return(new_df)
    }
    
    pred_res <- add_row(pred_res, 2, m2s)
    pred_res <- add_row(pred_res, 3, m3s)
    pred_res <- add_row(pred_res, 4, m4s)
    pred_res <- add_row(pred_res, 5, m5s)
    pred_res <- add_row(pred_res, 6, m6s)
    pred_res <- add_row(pred_res, 7, m7s)
    pred_res <- add_row(pred_res, 8, m8s)
    pred_res <- add_row(pred_res, 9, m9s)
    pred_res <- add_row(pred_res, 10, m10s)
    pred_res <- add_row(pred_res, 11, m11s)
    pred_res <- add_row(pred_res, 12, m12s)
    pred_res <- add_row(pred_res, 13, m13s)
    pred_res <- add_row(pred_res, 14, m14s)
    pred_res <- add_row(pred_res, 15, m15s)
    
    m15_imp <- data.frame(m15_imp)
    colnames(m15_imp) <- colnames(m15_features)
    return(list(pred_res, m1s, m2s, m3s, m4s, m5s, m6s, m7s, m8s, m9s, m10s, m11s, m12s, m13s, m14s, m15s, m15_imp))
    
  }
  
  # Outcomes based on splitting at median:
  t1 <- proc.time()
  list[OS_pred, OS_m1s, OS_m2s, OS_m3s, OS_m4s, OS_m5s, OS_m6s, OS_m7s, OS_m8s, OS_m9s, OS_m10s, OS_m11s, OS_m12s, OS_m13s, OS_m14s, OS_m15s, OS_m15_imp] <- run_a_bunch(mibi_table$OS_high,reps,"OS_high")
  list[PFS_pred, PFS_m1s, PFS_m2s, PFS_m3s, PFS_m4s, PFS_m5s, PFS_m6s, PFS_m7s, PFS_m8s, PFS_m9s, PFS_m10s, PFS_m11s, PFS_m12s, PFS_m13s, PFS_m14s, PFS_m15s, PFS_m15_imp] <- run_a_bunch(mibi_table$PFS_high,reps,"PFS_high")
  t2 <- proc.time()
  t2 - t1 
  
  # Save predictive performance results - outcome var / data version / seed / reps 
  write.csv(OS_pred,sprintf("Figure_6/OS_prediction_%s_seed%d_reps%d.csv",dv,seed,reps))
  write.csv(PFS_pred,sprintf("Figure_6/PFS_prediction_%s_seed%d_reps%d.csv",dv,seed,reps))
  
  # Combined output table for displaying in the paper
  combined_prediction <- full_join(OS_pred,PFS_pred, by=c("model")) %>% tibble::add_column(feature_set_sizes,.after="model") %>% select("model","feature_set_sizes","mean.x","sd.x","mean.y","sd.y") %>% rename( "mean.Survival" = "mean.x", "sd.Survival" = "sd.x", "mean.Recurrence" = "mean.y", "sd.Recurrence"= "sd.y")
  write.csv(combined_prediction,sprintf("Figure_6/Combined_prediction_%s_seed%d_reps%d.csv",dv,seed,reps))
  
  # Save full model feature importance results
  write.csv(OS_m15_imp,sprintf("Figure_6/OS_fullmodel_imp_%s_seed%d_reps%d.csv",dv,seed,reps))
  write.csv(PFS_m15_imp,sprintf("Figure_6/PFS_fullmodel_imp_%s_seed%d_reps%d.csv",dv,seed,reps))
  
  # Stack overflow - https://stackoverflow.com/questions/32011873/force-summary-to-report-the-number-of-nas-even-if-none
  even_summary <- function(v){
    if(!any(is.na(v))){
      res <- c(summary(v),"NA's"=0)
    } else{
      res <- summary(v)
    }
    return(res)
  }
  
  # Save aggregations of feature importance results
  OS_imp_agg <- as.data.frame(t(do.call(cbind, lapply(OS_m15_imp, even_summary))))
  OS_imp_agg <- OS_imp_agg[order(OS_imp_agg$Median, decreasing=TRUE),]
  write.csv(OS_imp_agg,sprintf("Figure_6/OS_impagg_%s_seed%d_reps%d.csv",dv,seed,reps))
  PFS_imp_agg <- as.data.frame(t(do.call(cbind, lapply(PFS_m15_imp, even_summary))))
  PFS_imp_agg <- PFS_imp_agg[order(PFS_imp_agg$Median, decreasing=TRUE),]
  write.csv(PFS_imp_agg,sprintf("Figure_6/PFS_impagg_%s_seed%d_reps%d.csv",dv,seed,reps))
  
  # Save all the objects
  save.image(file=sprintf('Figure_6/allfeatures_modelcomparison_%s_seed%d_reps%d.RData',dv,seed,reps))

}


### PLOTTING AUC RESULTS (from loaded in data)

for (dv in data_versions) {

  load(sprintf('Figure_6/allfeatures_modelcomparison_%s_seed%d_reps%d.RData',dv,seed,reps))
  
  stacked_OS_res <- data.frame(rbind(cbind(rep(1,reps), OS_m1s),cbind(rep(2,reps), OS_m2s),cbind(rep(3,reps), OS_m3s),cbind(rep(4,reps), OS_m4s),cbind(rep(5,reps), OS_m5s),cbind(rep(6,reps), OS_m6s),cbind(rep(7,reps), OS_m7s),cbind(rep(8,reps), OS_m8s),cbind(rep(9,reps), OS_m9s),cbind(rep(10,reps), OS_m10s),cbind(rep(11,reps), OS_m11s),cbind(rep(12,reps), OS_m12s),cbind(rep(13,reps), OS_m13s),cbind(rep(14,reps), OS_m14s),cbind(rep(15,reps), OS_m15s)))
  colnames(stacked_OS_res) <- c('Model','AUC')
  stacked_OS_res$Model <- as.factor(stacked_OS_res$Model)
  stacked_OS_res$AUC <- as.numeric(stacked_OS_res$AUC)
  
  ggplot(stacked_OS_res, aes(x = Model, y = AUC)) +
    geom_boxplot() +
    geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.002) +
    ylim(0,1) +
    theme_grey(base_size = 22) +
    geom_hline(yintercept = 0.5,linetype=2, col=2)
  
  ggsave(sprintf("Figure_6/OS_AUCs_%s_seed%d_reps%d.png",dv,seed,reps), width=744/90, height=610/90)
  
  stacked_PFS_res <- data.frame(rbind(cbind(rep(1,reps), PFS_m1s),cbind(rep(2,reps), PFS_m2s),cbind(rep(3,reps), PFS_m3s),cbind(rep(4,reps), PFS_m4s),cbind(rep(5,reps), PFS_m5s),cbind(rep(6,reps), PFS_m6s),cbind(rep(7,reps), PFS_m7s),cbind(rep(8,reps), PFS_m8s),cbind(rep(9,reps), PFS_m9s),cbind(rep(10,reps), PFS_m10s),cbind(rep(11,reps), PFS_m11s),cbind(rep(12,reps), PFS_m12s),cbind(rep(13,reps), PFS_m13s),cbind(rep(14,reps), PFS_m14s),cbind(rep(15,reps), PFS_m15s)))
  colnames(stacked_PFS_res) <- c('Model','AUC')
  stacked_PFS_res$Model <- as.factor(stacked_PFS_res$Model)
  stacked_PFS_res$AUC <- as.numeric(stacked_PFS_res$AUC)
  
  ggplot(stacked_PFS_res, aes(x = Model, y = AUC)) +
    geom_boxplot() +
    geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.002) +
    ylim(0,1) +
    theme_grey(base_size = 22) +
    geom_hline(yintercept = 0.5,linetype=2, col=2)
  
  ggsave(sprintf("Figure_6/PFS_AUCs_%s_seed%d_reps%d.png",dv,seed,reps), width=744/90, height=610/90)

}

