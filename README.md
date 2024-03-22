This repository hosts the data and analysis for the paper "The spatial structure of the tumor immune microenvironment can explain and predict patient response in high-grade serous carcinoma."

**Analysis Files:**

1_generate_TSNE_plot.ipynb
- Jupyter notebook that creates the figure 2 tsne plot

2_read_cell_data_and_generate_features.ipynb
- Jupyter notebook that reads in cell data (from Data/cell_table.csv) and produces composition, spatial, and network features for each sample
- Visualizations and intermediate text file results are produced and saved in Figure_* folders
- Intermediate data objects are created and stored in an analysis_objects folder (not included in the repository but created by the scripts if they are re-run)

3_read_clinical_data.ipynb
- Jupyter notebook that reads in clinical/immunohistochemical data (from Data/tma_info.csv) and outcomes
- Visualizations of clinical/immunohistochemical data and outcomes are produced and saved in Figure_* folders

4_combine_features.ipynb
- Jupyter notebook that combines spatial / network features and clinical/immunohistochemical data from to produce final data files which are saved in Data/mibi_table_*.csv

5_cox_regressions_univariate.R
- R script that reads in the processed data (from Data/mibi_table_*.csv) and performs univariate Cox regressions and saves results

6_cox_regressions_multivariate.Rmd
- R markdown that performs multivariate Cox regressions using the top univariate features and saves results

7_allfeatures_modelcomparison.R
- R script that reads in the processed data (from Data/mibi_table_*.csv) and performs comparisons of random forest out-of-sample predictive performance for 15 models including different subsets of features for both OS and PFS outcomes and saves results
 
8_aggregate_feature_importance.R
- R script that produces visualizations of the aggregated Gini importance scores across model evaluation iterations for the model containing all features for both OS and PFS outcomes

9_prediction_sensitivity_tests.R
- Performs hyperparameter selection using 5-fold cross validation for all 15 models and saves results

**Helper Files:**

spatial_net_helper.py
- Includes functions used to generate spatial and network features for the samples

viz_helper.py
- Includes functions used to produce visualizations of the samples

survival_helper.R
- Includes functions used to perform univariate survival analysis and visualize results
