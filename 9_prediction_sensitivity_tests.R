library(randomForest)
library(pROC)
library(caTools)
library(mltools)
library(data.table)
library(tidyverse)
library(gsubfn)

# For reproducing the same results as in the paper 
seed <- 19
set.seed(seed) 
dv <- 'main'

# Load data table 
mibi_table <- read.csv(sprintf('Data/mibi_table_%s.csv', dv), header = TRUE, sep = ",", dec = ".", check.names=FALSE)
# Shuffle data randomly
mibi_table <- mibi_table[sample(1:nrow(mibi_table)), ]

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

nodesizes <- 1:20
mtrys <- 2:12
M <- 15

# K-Fold Cross Validation 
start_time <- Sys.time()
k<-5
folds <- cut(seq(1, nrow(mibi_table)), breaks = k, labels = FALSE)
feature_sets <- list(m1_features, m2_features, m3_features, m4_features, m5_features, m6_features, m7_features, m8_features, m9_features, m10_features, m11_features, m12_features, m13_features, m14_features, m15_features)
o_vars <- list(mibi_table$OS_high, mibi_table$PFS_high)
# AUC results -- output variables x test / train scores x (matrix - folds x model) 
ress_rf <- list(list(matrix(0,k,M),matrix(0,k,M)),list(matrix(0,k,M),matrix(0,k,M)))
# Hyperparameter selection -- output variables x (folds x models)
rf_hyperparams <- list(array(0, dim=c(k,M,2)),array(0, dim=c(k,M,2)))

for (o in 1:2) { # output variables 
  Y <- o_vars[[o]]
  for (j in 1:M) { # models 
    feature_subset <- feature_sets[[j]]
    for(i in 1:k) { # folds 
      
      print(sprintf("Outcome: %d, Model: %d, Fold: %d", o, j, i))
      
      # Split the data into training and testing sets for this fold 
      test_indices <- which(folds == i)
      test_x <- feature_subset[test_indices, ]
      train_x <- feature_subset[-test_indices, ]
      test_y <- Y[test_indices]
      train_y <- Y[-test_indices]
      
      # Hyperparameter selection from training set
      cv_folds <- cut(seq(1, nrow(train_x)), breaks = k, labels = FALSE)
      best_rf_auc = 0
      best_params_rf = c(0,0)
      for(cvi in 1:k) {
        print(sprintf("Cross validation fold %d", cvi))
        cv_test_indices <- which(cv_folds == cvi)
        cv_test_x = train_x[cv_test_indices,]
        cv_train_x = train_x[-cv_test_indices,]
        cv_test_y = train_y[cv_test_indices]
        cv_train_y = train_y[-cv_test_indices]
        
        # Replace NA values 
        cv_test_x <- na.roughfix(cv_test_x)
        cv_train_x <- na.roughfix(cv_train_x)
        
        #for random forest
        for (ns in nodesizes) {
          for (mt in mtrys) {
            #print(sprintf("Nodesize: %d, Mtry: %d", ns, mt))
            cv_classifier = randomForest(x = cv_train_x, y = as.factor(cv_train_y), nodesize=ns, mtry=mt)
            cv_pred_y = predict(cv_classifier, newdata=cv_test_x, type="prob")
            cv_prob = cv_pred_y[,2]
            cv_auc = auc(roc(response = cv_test_y, predictor = cv_prob, auc=TRUE, direction="<", quiet = TRUE))
            if (cv_auc > best_rf_auc){
              best_rf_auc <- cv_auc
              best_params_rf <- c(ns,mt)
            }
          }
        }
        
      }
      
      # Replace NA values
      test_x <- na.roughfix(test_x)
      train_x <- na.roughfix(train_x)
      
      # use factor to force classification
      rf = randomForest(x = train_x, y = as.factor(train_y), nodesize = best_params_rf[1], mtry = best_params_rf[2])
      # Predict probabilities for the test dataset
      pred_rf = predict(rf, newdata = test_x, type = "prob")
      
      # Predict probabilities on the training dataset itself
      pred_rf_train = predict(rf, newdata = train_x, type = "prob")
      
      # Random forest performance on test data
      prob_rf = pred_rf[,2]
      rf_auc = auc(roc(response = test_y, predictor = prob_rf, auc=TRUE, direction="<", quiet = TRUE))
      
      # Random forest performance on train data
      prob_rf_train = pred_rf_train[,2]
      rf_auc_train <- auc(roc(response = train_y, predictor = prob_rf_train, auc=TRUE, direction="<", quiet=TRUE))
      
      ress_rf[[o]][[1]][i,j] <- rf_auc
      ress_rf[[o]][[2]][i,j] <- rf_auc_train
      rf_hyperparams[[o]][i,j,] <- best_params_rf
      print(sprintf("Time elapsed: %f", Sys.time() - start_time))
      
    }
  }
}

print(sprintf("Total time elapsed: %f", Sys.time() - start_time))
save(ress_rf, rf_hyperparams, file = sprintf("k_folds_res_%s_seed%d.RData",dv,seed))

#### RANDOM FOREST #########################################

# Overall Survival Results
# Overall Survival test performance 
res <- ress_rf[[1]][[1]]
res <- as.data.frame(res)
res[nrow(res) + 1,] <- colMeans(res)
rownames(res) <- c('fold1_auc','fold2_auc','fold3_auc','fold4_auc','fold5_auc','avg_auc')
# Hyperparameter selection results - nodesize
ns_res <- rf_hyperparams[[1]][,,1]
ns_res <- as.data.frame(ns_res)
rownames(ns_res) <- c('fold1_nodesize','fold2_nodesize','fold3_nodesize','fold4_nodesize','fold5_nodesize')
# Hyperparameter selection results - mtry
mt_res <- rf_hyperparams[[1]][,,2]
mt_res <- as.data.frame(mt_res)
rownames(mt_res) <- c('fold1_mtry','fold2_mtry','fold3_mtry','fold4_mtry','fold5_mtry')
os_hyp_output <- rbind(res,ns_res,mt_res)
colnames(os_hyp_output) <- c('m1','m2','m3','m4','m5', 'm6', 'm7', 'm8', 'm9', 'm10', 'm11', 'm12', 'm13','m14','m15')
write.csv(os_hyp_output,sprintf("Figure_6/OS_hyp_select_%s_seed%d.csv",dv,seed))

# Progression-free survival Results
# Progression-free survival test performance
res <- ress_rf[[2]][[1]]
res <- as.data.frame(res)
res[nrow(res) + 1,] <- colMeans(res)
rownames(res) <- c('fold1_auc','fold2_auc','fold3_auc','fold4_auc','fold5_auc','avg_auc')
# Hyperparameter selection results - nodesize
ns_res <- rf_hyperparams[[2]][,,1]
ns_res <- as.data.frame(ns_res)
rownames(ns_res) <- c('fold1_nodesize','fold2_nodesize','fold3_nodesize','fold4_nodesize','fold5_nodesize')
# Hyperparameter selection results - mtry
mt_res <- rf_hyperparams[[2]][,,2]
mt_res <- as.data.frame(mt_res)
rownames(mt_res) <- c('fold1_mtry','fold2_mtry','fold3_mtry','fold4_mtry','fold5_mtry')
PFS_hyp_output <- rbind(res,ns_res,mt_res)
colnames(PFS_hyp_output) <- c('m1','m2','m3','m4','m5', 'm6', 'm7', 'm8', 'm9', 'm10', 'm11', 'm12', 'm13','m14','m15')
write.csv(PFS_hyp_output,sprintf("Figure_6/PFS_hyp_select_%s_seed%d.csv",dv,seed))

# Look at Overfitting 

# Overall Survival - Train performance
res_os_train <- ress_rf[[1]][[2]]
res_os_train <- as.data.frame(res_os_train)
res_os_train[nrow(res_os_train) + 1,] <- colMeans(res_os_train)
colnames(res_os_train) <- c('m1','m2','m3','m4','m5', 'm6', 'm7', 'm8', 'm9', 'm10', 'm11', 'm12', 'm13','m14','m15')
rownames(res_os_train) <- c('fold1_train','fold2_train','fold3_train','fold4_train','fold5_train','avg_train')

# Progression-free survival - Train performance
res_PFS_train <- ress_rf[[2]][[2]]
res_PFS_train <- as.data.frame(res_PFS_train)
res_PFS_train[nrow(res_PFS_train) + 1,] <- colMeans(res_PFS_train)
colnames(res_PFS_train) <- c('m1','m2','m3','m4','m5', 'm6', 'm7', 'm8', 'm9', 'm10', 'm11', 'm12', 'm13','m14','m15')
rownames(res_PFS_train) <- c('fold1_train','fold2_train','fold3_train','fold4_train','fold5_train','avg_train')

