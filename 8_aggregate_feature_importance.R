library(tidyverse)
library(miscTools) # for colMedians

display_names <- function(init){
  disp <- str_replace_all(init, "_"," ")
  disp <- str_replace_all(disp, "med dist","med nn dist")
  return(disp)
}

get_types <- function(feat_col){
  out <- c()
  for (i in 1:length(feat_col)){
    curr <- feat_col[i]
    if (curr %in% names(clinical_features)){
      out = c(out,"clinical features")
    }
    if (curr %in% names(network_features)){
      out <- c(out,"network features")
    }
    if (curr %in% names(spatial_features)){
      out <- c(out,"spatial features")
    }  
    if (curr %in% names(composition_features)){
      out <- c(out,"composition features")
    }
  }
  return(out)
}

get_types_detailed <- function(feat_col, focus){
  out <- c()
  for (i in 1:length(feat_col)){
    curr <- feat_col[i]
    if (curr %in% names(clinical_features)){
      out = c(out,"clinical features")
    }
    if (curr %in% names(network_features)){
      if (grepl("contact", curr)){
        out <- c(out,"contact enrichment features")
      }
      if (grepl("assort", curr)){
        out <- c(out,"assortativity features")
      }
      if (grepl("region", curr)){
        out <- c(out,"region size features")
      }
    }
    if (curr %in% names(spatial_features)){
      out <- c(out,"spatial features")
    }  
    if (curr %in% names(composition_features)){
      out <- c(out,"composition features")
    }
  }
  out_final <- c()
  for (i in 1:length(out)){
    curr <- out[i]
    if (curr == focus){
      out_final <- c(out_final,curr)
    } else{
      out_final <- c(out_final,"other")
    }
  }
  return(out_final)
}

feature_imp_plot <- function(imp_mat, title, plot_colors, dv){
  
  # Sort by column medians before stacking
  a <- colMedians(imp_mat, na.rm = TRUE)
  b <- sort(a,decreasing=TRUE)
  imp_mat <- imp_mat[,names(b)]
  
  top_ten <- names(head(b,10))
  trimmed <- imp_mat[top_ten]
  
  # Stack and plot feature importances across the replicates
  stacked_imp = stack(imp_mat)
  stacked_imp <- stacked_imp %>% rename(imps = values, features = ind)
  stacked_imp <- mutate(stacked_imp, type=get_types(features))
  p1 <- ggplot(data=stacked_imp) +
    geom_boxplot(aes(x=features, y=imps, color=type)) +
    scale_colour_manual(values = plot_colors) +
    ggtitle(title) +
    xlab("All features") +
    ylab("Gini importance")+
    theme_classic(base_size=22)+
    ylim(0,1.75)+
    #theme(axis.text.x = element_text(angle = 90, size = 5))
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),legend.position = "none")
  ggsave(sprintf("Figure_7/%s_%s_full_importance.png",title,dv), width=744/90, height=610/90)
  
  # Also plot just the top ten features
  stacked_trimmed = stack(trimmed)
  stacked_trimmed <- stacked_trimmed %>% rename(imps = values, features = ind)
  stacked_trimmed <- mutate(stacked_trimmed, type=get_types(features))
  stacked_trimmed <- stacked_trimmed %>% mutate(features = display_names(features))
  levels_adjusted <- display_names(top_ten)
  
  stacked_trimmed$features <- factor(stacked_trimmed$features, levels=levels_adjusted)
  p2 <- ggplot(data=stacked_trimmed) +
    geom_boxplot(aes(x=features, y=imps, color=type)) +
    scale_colour_manual(values = plot_colors) +
    #ggtitle(paste(title))+
    xlab("Top ten features") +
    ylab("Gini importance")+
    theme_classic(base_size=22)+
    ylim(0,1.75)+
    scale_x_discrete(limits = rev(levels(stacked_trimmed$features)))+
    theme(axis.text.x = element_text(angle = 90, size = 14),legend.position = "none") + coord_flip()
  ggsave(sprintf("Figure_7/%s_%s_top_ten_importance.png",title,dv), width=850/90, height=610/90)
}

feature_imp_plot_detail <- function(imp_mat, title, plot_colors, focus, dv){
  
  # Sort by column medians before stacking
  a <- colMedians(imp_mat, na.rm = TRUE)
  b <- sort(a,decreasing=TRUE)
  imp_mat <- imp_mat[,names(b)]
  
  top_ten <- names(head(b,10))
  trimmed <- imp_mat[top_ten]
  
  # Stack and plot feature importances across the replicates
  stacked_imp = stack(imp_mat)
  stacked_imp <- stacked_imp %>% rename(imps = values, features = ind)
  stacked_imp <- mutate(stacked_imp, type=get_types_detailed(features, focus))
  p1 <- ggplot(data=stacked_imp) +
    geom_boxplot(aes(x=features, y=imps, color=type, fill=type)) +
    scale_colour_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    ggtitle(title) +
    xlab("All features") +
    ylab("Gini importance")+
    theme_gray(base_size=22)+
    #theme(axis.text.x = element_text(angle = 90, size = 5))
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),legend.position = "none")
  ggsave(sprintf("Figure_7/%s_%s_full_importance_%s.png",title,dv,focus), width=744/90, height=610/90)
  
}

data_versions <- c('missing_na','main','primary', 'main_th100')

for (dv in data_versions) {

  load(sprintf('Figure_6/allfeatures_modelcomparison_%s_seed19_reps500.RData',dv))
  
  # Unique color for each feature type consistent across plots 
  types <- unique(get_types(names(OS_m15_imp)))
  plot_colors <- setNames(object = c("#CC79A7","#009E73","#0072B2","#E69F00"), nm = types)
  
  # Look at feature importance in full models 
  feature_imp_plot(OS_m15_imp, "Overall survival, full model", plot_colors, dv)
  feature_imp_plot(PFS_m15_imp, "Progression-free survival, full model", plot_colors, dv)
  
  # Separate the different feature types in the aggregated plots (only for main text results)
  if (dv=="main"){
  
    # Clinical 
    plot_colors <- setNames(object = c("#CC79A7","#D3D3D3"), nm = c("clinical features","other"))
    feature_imp_plot_detail(OS_m15_imp, "Overall survival, full model", plot_colors, "clinical features", dv)
    feature_imp_plot_detail(PFS_m15_imp, "Progression-free survival, full model", plot_colors, "clinical features", dv)
    
    # Contact enrichment features
    plot_colors <- setNames(object = c("#E69F00","#D3D3D3"), nm = c("contact enrichment features","other"))
    feature_imp_plot_detail(OS_m15_imp, "Overall survival, full model", plot_colors, "contact enrichment features", dv)
    feature_imp_plot_detail(PFS_m15_imp, "Progression-free survival, full model", plot_colors, "contact enrichment features", dv)
    
    # Assortativity features
    plot_colors <- setNames(object = c("#E69F00","#D3D3D3"), nm = c("assortativity features","other"))
    feature_imp_plot_detail(OS_m15_imp, "Overall survival, full model", plot_colors, "assortativity features", dv)
    feature_imp_plot_detail(PFS_m15_imp, "Progression-free survival, full model", plot_colors, "assortativity features", dv)
    
    # Region size features
    plot_colors <- setNames(object = c("#E69F00","#D3D3D3"), nm = c("region size features","other"))
    feature_imp_plot_detail(OS_m15_imp, "Overall survival, full model", plot_colors, "region size features", dv)
    feature_imp_plot_detail(PFS_m15_imp, "Progression-free survival, full model", plot_colors, "region size features", dv)
    
    
    # Spatial features
    plot_colors <- setNames(object = c("#0072B2","#D3D3D3"), nm = c("spatial features","other"))
    feature_imp_plot_detail(OS_m15_imp, "Overall survival, full model", plot_colors, "spatial features", dv)
    feature_imp_plot_detail(PFS_m15_imp, "Progression-free survival, full model", plot_colors, "spatial features", dv)
    
    # Composition features
    plot_colors <- setNames(object = c("#009E73","#D3D3D3"), nm = c("composition features","other"))
    feature_imp_plot_detail(OS_m15_imp, "Overall survival, full model", plot_colors, "composition features", dv)
    feature_imp_plot_detail(PFS_m15_imp, "Progression-free survival, full model", plot_colors, "composition features", dv)
  
  }

}

