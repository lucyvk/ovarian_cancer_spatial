# Author: Lucy Van Kleunen
# Helper functions for survival analysis 

library(survival)
library(survminer)
library(dplyr)
library(tibble)
library(ggbreak)

get_types <- function(feat_col){
  out <- c()
  for (i in 1:length(feat_col)){
    curr <- feat_col[i]
    if (curr %in% names(clinical_features)){
      out = c(out,"clinical features")
    }
    else if (curr %in% names(network_features)){
      out <- c(out,"network features")
    }
    else if (curr %in% names(spatial_features)){
      out <- c(out,"spatial features")
    }  
    else if (curr %in% names(composition_features)){
      out <- c(out,"composition features")
    }
    else {
      print("not found")
      print(curr)
    }
  }
  return(out)
}

survival_curve <- function(time, status, covariates_scaled, filename){
  # Kaplan-Meier Curve Plot
  fit <- survfit(Surv(time,status) ~ 1, conf.type="log-log")
  df_surv <- cbind(time,status,covariates_scaled)
  plot <- ggsurvplot(fit,data=df_surv)
  
  # save survplot with ggsave
  # from https://stackoverflow.com/questions/70751317/ggsave-error-in-usemethodgrid-draw-no-applicable-method-for-grid-draw-ap
  ggsave_workaround <- function(g){survminer:::.build_ggsurvplot(x = g,
                                                                 surv.plot.height = NULL,
                                                                 risk.table.height = NULL,
                                                                 ncensor.plot.height = NULL)}
  
  g_to_save <- ggsave_workaround(plot)
  
  
  ggsave(filename, plot=g_to_save, width=6, height=5, units="in", dpi=300)
}

# Note - expects all features to be numeric 
univariate_cox <- function(df){
  formulas <- sapply(covariates, function(x) as.formula(paste('Surv(time, status) ~', x)))
  models <- lapply(formulas, function(x){coxph(x, data = df)})
  results <- lapply(models,
                    function(x){ 
                      x <- summary(x)
                      example_x <<- x
                      coef<-format(x$coef[1], scientific=FALSE);#coefficient
                      se <- format(x$coef[3], scientific=FALSE)
                      hazard_ratio <- x$coef[2]
                      hr_confint_lower <- x$conf.int[,"lower .95"]
                      hr_confint_upper <- x$conf.int[,"upper .95"]
                      p_value <- x$coef[5]
                      n <- x$n # number of samples (after NA excluded)
                      nevent <- x$nevent #nevent
                      res<- as.numeric(c(n, nevent, coef, se, hazard_ratio, hr_confint_lower,hr_confint_upper, p_value))
                      names(res)<-c("n","nevent","coef","se","hazard_ratio","hr_confint_lower","hr_confint_upper","p_value")
                      return(res)
                    })
  res_df <- as.data.frame(t(as.data.frame(results, check.names = FALSE)))
  return(res_df)
}

# Separate function to run the regression for BRCA status which is a factor 
univariate_cox_brca <- function(df){
  
  res_brca <- coxph(Surv(time, status)~brca, data = df)
  
  # Add brca results in with the other ones (separate out coefficients for two categories relative to baseline)
  x <- summary(res_brca)
  print(x)
  coef1<-format(x$coef[1], scientific=FALSE);#coefficient
  se1 <- format(x$coef[5], scientific=FALSE)
  hazard_ratio1 <- x$coef[3];#exp(beta)
  hr_confint_lower1 <- x$conf.int[,"lower .95"][1]
  hr_confint_upper1 <- x$conf.int[,"upper .95"][1]
  p_value1<- x$coef[9]
  
  coef2<-format(x$coef[2], scientific=FALSE);#coefficient
  se2 <- format(x$coef[6], scientific=FALSE)
  hazard_ratio2 <- x$coef[4];#exp(beta)
  hr_confint_lower2 <- x$conf.int[,"lower .95"][2]
  hr_confint_upper2 <- x$conf.int[,"upper .95"][2]
  p_value2 <- x$coef[10]
  
  n <- x$n # number of samples (after NA excluded)
  nevent <- x$nevent #nevent
  res1<- as.numeric(c(n, nevent, coef1, se1, hazard_ratio1, hr_confint_lower1, hr_confint_upper1, p_value1))
  res2<- as.numeric(c(n, nevent, coef2, se2, hazard_ratio2, hr_confint_lower2, hr_confint_upper2, p_value2))
  return(list(res1,res2))
  
}

run_cox <- function(time, status, df_surv){
  
  # Running all the Univariate Cox regressions for Overall Survival for numeric features 
  cox_res <- univariate_cox(df_surv)
  p_adjusted <- p.adjust(cox_res$p_value, method="BH")
  cox_res <- cbind(cox_res,p_adjusted)
  cox_res <- cox_res[order(cox_res$p_value),]  
  return(cox_res)
}

# Write the significant covariates to file
write_significant <- function(res_df, filename){
  res_df_trimmed <- res_df %>% dplyr::filter(p_value < 0.05)
  res_df_trimmed <- tibble::rownames_to_column(res_df_trimmed, var = "covariate")
  # sort by hazard ratio rather than p value
  res_df_trimmed <- res_df_trimmed[order(res_df_trimmed$hazard_ratio),]
  res_df_trimmed$covariate <- factor(res_df_trimmed$covariate, levels=res_df_trimmed$covariate)
  write.csv(res_df_trimmed,filename)
}

# Plot the significant covariates
plot_significant <- function(res_df, outcome, BREAK, filename){
  
  plot_colors <- setNames(object = c("#CC79A7","#009E73","#0072B2","#E69F00"), nm = c("clinical features","composition features","spatial features","network features"))
  
  res_df_trimmed <- res_df %>% dplyr::filter(p_value < 0.05)
  res_df_trimmed <- tibble::rownames_to_column(res_df_trimmed, var = "covariate")
  # sort by hazard ratio rather than p value
  res_df_trimmed <- res_df_trimmed[order(res_df_trimmed$hazard_ratio),]
  # add type col
  res_df_trimmed <- mutate(res_df_trimmed, type=get_types(res_df_trimmed$covariate))
  res_df_trimmed$covariate <- display_names(res_df_trimmed$covariate)
  res_df_trimmed$covariate <- factor(res_df_trimmed$covariate, levels=res_df_trimmed$covariate)

  if (BREAK){
    # Hardcoded break in the plot for the outlier HR plotting 
    break_seq <-c(seq(0.5,2,by=0.5),seq(6,26,by=10))
  } else{
    break_seq <- seq(0.5,2.5,by=0.5)
  }
  
  plot1 <- ggplot(data=res_df_trimmed, aes(x=covariate, y=hazard_ratio, colour=type)) +
    geom_point(shape = 15) +
    geom_hline(yintercept = 1, lty=2, size =1) +
    geom_linerange(aes(ymin=hr_confint_lower, ymax=hr_confint_upper)) +
    scale_colour_manual(values = plot_colors) +
    coord_flip() + # flip so all x-y is inverted here
    ggtitle(sprintf("%s Univariate Cox Regression Significant Covariates",outcome)) +
    xlab("Covariates") +
    ylab("Hazard Ratio") +
    scale_y_continuous(breaks = break_seq, limits = c(0, 2.75))
  if (BREAK){
    plot2 <- plot1 + scale_y_break(c(2, 6), scales=0.2)
  } else{
    plot2 <- plot1
  }
  ggsave(filename, plot=plot2, width=8, height=5, units="in", dpi=300)

}

# display name adjustment
display_names <- function(init){
  disp <- gsub("_"," ",init, fixed=TRUE)
  disp <- gsub("med dist","med nn dist",disp, fixed=TRUE)
  disp <- gsub("Blood vessel","vascular endothelial cell",disp, fixed=TRUE)
  disp <- gsub("Lymphatic vessel","lymphatic endothelial cell",disp, fixed=TRUE)
  disp <- gsub("B cells","B cell",disp, fixed=TRUE)
  return(disp)
}
