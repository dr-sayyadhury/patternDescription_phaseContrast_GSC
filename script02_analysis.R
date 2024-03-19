
### PACKAGES TO LOAD
###------------------------------------------------------------------#### 
library(ggplot2)  
library(colorspace)
library(tidyr)
library(dplyr)
library(ggpubr)
#library(ggrepel)
library(effectsize)
library(ggthemes)
#library(scales)
library(forcats)
library(PCAtools) ### note this version 2.6.0 requires a downgrade in matrixStats - remotes::install_github("HenrikBengtsson/matrixStats", ref="1.1.0")
#library(matrixStats)

### FILE/FOLDER PATHS
###------------------------------------------------------------------#### 
path_to_repo <- '/Users/mystique27m/Documents/research/PostdoctoralResearch_2020/Projects/PatternRecognitionInGSCs_UsingCV/'
script_path <- paste0(path_to_repo,'/scripts/')
out <- paste0(path_to_repo,'out_intermediate/')
dir.create(out)
figures_filepath <- paste0(path_to_repo, 'results/figures/')
dir.create(figures_filepath)

### LOAD DATA
###------------------------------------------------------------------####
source(paste0(script_path,"source_scripts/sourceData_General_variables_and_functions.R"))

GSC.gsva <- read.csv(paste0(path_to_repo,'/datasets/bulk_RNA/GSC.gsva.csv'), row.names=1)
phase <- read.csv(paste0(path_to_repo, '/results/processed_pixel_feature_values_for_QC_and_GLCM_Grn/final_dataset_afterQC.csv'), row.names = 1)

df_list_tp_subsetted <- readRDS(paste0(out, '/script01_1_image_data_by_confluencies/df_list_tp_subsetted.rds'))


#NORMALIZATION AND STANDARDIZATION BY CONFLUENCY GROUPS
phase <- phase[shorter]

pt.size=1.5
features=34
new_df_list <- list()
for (i in 1:9){
  d <- df_list_tp_subsetted[[i]]
  new_df <- data.frame()
  
  n=length(levels(factor(d$Sample)))
  for (l in 1:n){
    s <- levels(factor(d$Sample))[l]
    df <- subset(d, Sample==s)
    meta <- df[1:5]
    df <- df[6:features]
    df <- apply(df, 2, FUN=function(x)if(mean(x)==0)
    {x=x}
    else{x/max(x)}) %>% as.data.frame(.) ### note keep this
    df <- apply(df,1, standardize)
    df <- cbind(meta,t(df)) 
    df <- df %>% .[order(.$TimePt),]
    
    meta <- df[1:5]
    df <- df[6:features]  %>% t(.) %>% as.data.frame(.)
    rownames(meta) <- colnames(df)
    con_grp <- rep(paste0('C',i), nrow(meta))
    meta <- cbind(con_grp, meta)
    df <- cbind(meta, t(df))
    
    colnames(df) <- c('con_grp', colnames(phase))
    new_df <- rbind(new_df, df)
    
  }
  new_df_list[[i]] <- new_df
}


#PLOT PCA FOR CONFLUENCY GROUPS 1-9
#confluency group 1 was excluded due to less than 10 samples available
library(cowplot)

features=35
pca_list <- list()
pca_loadings_list <- list()
pc1_loadings <- pc2_loadings <- pc3_loadings <- shorter[6:34] %>% as.data.frame(.)
g1 <- list()
g2 <- list()
sc <- list()
pt.size=3.5
var_pc <- as.data.frame(matrix(ncol=1, nrow=0))

for (i in 1:9){
  new_df <- new_df_list[[i]]
  
  new_df[is.na(new_df)] <-0
  meta <- new_df[1:6]
  data <- new_df[7:features]
  rownames(data) <- rownames(meta)
  
  p <- PCAtools::pca(t(data), metadata = meta, center=T, scale=F)
  
  pcs <- which(cumsum(p$variance)>90)
  pc1 <- pcs[1]
  pc2 <- which(cumsum(p$variance)>75)[1]
  
  var <- sort(which((p$variance[1:length(pcs)-1]-p$variance[2:(length(pcs))]) > 0.1), decreasing=T)[1]
  minPC <- min(pc1, var)
  
  sc[[i]] <- PCAtools::screeplot(p, colBar='darkblue', axisLabSize = 9, components = 1:(pc1+1), title = '')
  pca_df <- cbind(meta,p$rotated)
  colnames(pca_df) <- c('con_grp','Media','Sample','TimePt','Well','Area', paste0('PC', 1:ncol(p$rotated)))
  var_pc <- rbind(var_pc, (p$variance[1] + p$variance[2]))
  
  pc1_loadings <- cbind(pc1_loadings, p[["loadings"]][["PC1"]])
  pc2_loadings <- cbind(pc2_loadings, p[["loadings"]][["PC2"]])
  pc3_loadings <- cbind(pc3_loadings, p[["loadings"]][["PC3"]])
  
  pca_df <- cbind(pca_df,data)
  pca_list[[i]] <- pca_df
  
}

min(var_pc)
max(var_pc)



features=35
pt.size=3.5

for (i in 1:9){
  pca_df <- pca_list[[i]]
  g1[[i]] <- ggplot(pca_df, aes(PC2, PC1, color=Sample)) +
    geom_jitter(size=pt.size) +
    theme_minimal() +
    scale_color_manual(values=c(cols_vivid, col_bold)) +
    theme(legend.position = 'none', 
          axis.title = element_text(size=12),
          axis.line = element_line(color='black'),
          axis.text.x = element_blank(),
          axis.text.y = element_text(angle=45, size=9)) +
    #geom_text_repel(aes(label=Sample), size=2.5) +
    xlab("PC2") + ylab("PC1")
  
  g2[[i]] <- ggplot(pca_df, aes(PC2, PC3, color=Sample)) +
    geom_jitter(size=pt.size) +
    theme_minimal() +
    scale_color_manual(values=c(cols_vivid, col_bold)) +
    theme(legend.position = 'none', 
          axis.line = element_line(color='black'),
          axis.text.x = element_blank(),
          axis.text.y = element_text(angle=45, size=9)) +
    #geom_text_repel(aes(label=Sample), size=3) +
    xlab("PC2") + ylab("PC3")
  
}


if(dir.exists(paste0(figures_filepath,'Sup_Fig_S2'))==F){
  dir.create(paste0(figures_filepath,'Sup_Fig_S2'))}


###------------------------------------------------------------------####
### FIGURE S2B
###------------------------------------------------------------------####
pdf(file = paste0(figures_filepath,'Sup_Fig_S2/S2_B_var.pdf'), width=42, height=6)
ggarrange(plotlist=sc, nrow=1, ncol=9)
dev.off()
#--------------------------------------------------------------------------------#
#Fig end



g1 <- list()


for (p in 1:9){
  
  pca_df <- pca_list[[p]] %>% as.data.frame()
  pca_df <- arrange(pca_df, PC1)
  g1[[p]] <- ggplot(pca_df, aes(PC1, reorder(Sample, PC1), fill=Sample)) +
    geom_violin() +
    geom_jitter(size=0.3) +
    theme_minimal() +
    scale_color_manual(values=c(cols_vivid, col_bold)) +
    theme(legend.position = 'none', 
          axis.title = element_text(size=12),
          axis.line = element_line(color='black'),
          axis.text.x = element_blank(),
          axis.text.y = element_text(angle=45, size=9)) +
    #geom_text_repel(aes(label=Sample), size=2.5) +
    xlab("PC1") + ylab("Sample")
  
}



###------------------------------------------------------------------####
### FIGURE S2C1
###------------------------------------------------------------------####
pdf(file = paste0(figures_filepath,'Sup_Fig_S2/S2_C1.pdf'), width=42, height=6)
ggarrange(plotlist=g1, nrow=1, ncol=9)
dev.off()
#--------------------------------------------------------------------------------#
#Fig end


g2 <- list()


for (p in 1:9){
  
  pca_df <- pca_list[[p]] %>% as.data.frame()
  pca_df <- arrange(pca_df, PC2)
  g2[[p]] <- ggplot(pca_df, aes(PC2, reorder(Sample, PC2), fill=Sample)) +
    geom_violin() +
    geom_jitter(size=0.3) +
    theme_minimal() +
    scale_color_manual(values=c(cols_vivid, col_bold)) +
    theme(legend.position = 'none', 
          axis.title = element_text(size=12),
          axis.line = element_line(color='black'),
          axis.text.x = element_blank(),
          axis.text.y = element_text(angle=45, size=9)) +
    #geom_text_repel(aes(label=Sample), size=2.5) +
    xlab("PC2") + ylab("Sample")
  
}

###------------------------------------------------------------------####
### FIGURE S2C2
###------------------------------------------------------------------####
pdf(file = paste0(figures_filepath,'Sup_Fig_S2/S2_C2.pdf'), width=42, height=6)
ggarrange(plotlist=g2, nrow=1, ncol=9)
dev.off()
### ------------------------------------------------------------------####
### Figure end

### BEFORE CHANGING PC DIRECTIONS
#CALCULATE MEAN PC SCORES BY SAMPLE FOR ALL CONFLUENCY GROUPS
pca_df_summary <- list()
for (l in 1:9){
  df <- pca_list[[l]]
  meta_df <- df[1:6]
  df <- cbind(meta_df['Sample'], df[7:20])
  df <- df %>% group_by(Sample) %>% summarise_at(vars(paste0('PC', 1:14)), mean)
  pca_df_summary[[l]] <- df
}


#CORRELATION BETWEEN MEAN PC SCORES WITH GSVA SCORES
pca_list_medPC_corGSVA <- list()
pca_list_medPC_corGSVA_pvalue <- list()
GSC.gsva <- GSC.gsva
sig <- rownames(GSC.gsva)

for (pc in 1:14){
  final_correlation_medPC <- as.data.frame(matrix(nrow=113))
  final_correlation_medPC_pvalue <- as.data.frame(matrix(nrow=113))
  
  for (con in 1:9){
    df <- pca_df_summary[[con]] ### select confluency group
    to_plot1 <- df[pc+1] %>% as.data.frame(.) ### select PC to focus on for the rest of the analysis
    gsva_short <- t(GSC.gsva) %>% as.data.frame(.) %>% .[df$Sample,]
    final_col <- data.frame()  
    for (x in 1:113){
      ctype <- sig[x] ### then select the cell types to focus
      #gsva_short <- gsva_short[,ctype]
      g <- gsva_short[,colnames(gsva_short) %in% ctype]
      #final <- data.frame()
      #to_plot2 <- g[col]
      cor <- cor.test(g, to_plot1[,1])
      #cor <- lm(to_plot1[,1]~to_plot2[,1])
      #r.sq <- summary(cor)$r.squared
      #m=cor$coefficients[2]
      p <- cor$estimate %>% round(.,3)
      d <- data.frame(label=ctype,correlation=p , p_value=cor$p.value)
      final_col <- rbind(final_col, d)
    }
    rownames(final_correlation_medPC_pvalue) <- rownames(final_correlation_medPC) <- final_col$label
    final_correlation_medPC <- cbind(final_correlation_medPC, final_col[2])
    final_correlation_medPC_pvalue <- cbind(final_correlation_medPC_pvalue, final_col[3])
  }
  final_correlation_medPC <- final_correlation_medPC[2:10]
  colnames(final_correlation_medPC) <- c(paste0('pixelBio', 1:9))
  
  pca_list_medPC_corGSVA[[pc]] <- final_correlation_medPC
  
  final_correlation_medPC_pvalue <- final_correlation_medPC_pvalue[2:10]
  colnames(final_correlation_medPC_pvalue) <- paste0('pixelBio', 1:9)
  pca_list_medPC_corGSVA_pvalue[[pc]] <- final_correlation_medPC_pvalue 
  
}


if(dir.exists(paste0(figures_filepath,'Sup_Fig_S3'))==F){
  dir.create(paste0(figures_filepath,'Sup_Fig_S3'))}


## See if we change the signs - materials and methods
#TO CHECK THE CONGRUENCY BETWEEN CONFLUENCY GROUPS
hm_col <- sequential_hcl(45, palette='Inferno', rev=T)
hm_col_neg<- sequential_hcl(45, palette='GnBu')

for (hm in 1:9){
  data <- pca_list_medPC_corGSVA[[hm]]
  #colours <- data$col
  #names(colours) <- data
  meta <- paste0('C', 1:9)
  data.cor <- cor(data)
  rownames(data.cor) <- colnames(data.cor) <- meta


  ###------------------------------------------------------------------####
  ### FIGURE S3B
  ###------------------------------------------------------------------####
  pdf(paste0(figures_filepath, 'Sup_fig_S3/S3_B', hm, '.pdf'), width=6.6, height=6)
  pheatmap::pheatmap(data.cor, 
                     #annotation_col = meta, 
                     #annotation_colors = list(cell_type_group=colours), 
                     show_rownames = T, 
                     annotation_legend = F,
                     legend = T,
                     border_color = NA, 
                     cluster_cols =F,
                     cluster_rows = F,
                     fontsize = 21,
                     color = hm_col
  )
  dev.off()
}
###------------------------------------------------------------------####
### FIGURE END





# CHANGE THE direction for easy comparison across confluency groups
### PC1
d <- pca_list[[5]]
d$PC1 <- (d$PC1)*(-1)
pca_list[[5]] <- d

d <- pca_list[[6]]
d$PC1 <- (d$PC1)*(-1)
pca_list[[6]] <- d

d <- pca_list[[7]]
d$PC1 <- (d$PC1)*(-1)
pca_list[[7]] <- d

d <- pca_list[[8]]
d$PC1 <- (d$PC1)*(-1)
pca_list[[8]] <- d

d <- pca_list[[9]]
d$PC1 <- (d$PC1)*(-1)
pca_list[[9]] <- d

### PC2

d <- pca_list[[2]]
d$PC2 <- (d$PC2)*(-1)
pca_list[[2]] <- d

d <- pca_list[[3]]
d$PC2 <- (d$PC2)*(-1)
pca_list[[3]] <- d

d <- pca_list[[4]]
d$PC2 <- (d$PC2)*(-1)
pca_list[[4]] <- d

d <- pca_list[[6]]
d$PC2 <- (d$PC2)*(-1)
pca_list[[6]] <- d

d <- pca_list[[9]]
d$PC2 <- (d$PC2)*(-1)
pca_list[[9]] <- d



### repeat correlations


### to be re-used in next script.
saveRDS(pca_list, paste0(results_path, 'pca_output_of_wholeImages_by_congrps.rds'))

remove(list=c('pca_list_medPC_corGSVA','pca_list_medPC_corGSVA_pvalue', 'pca_loadings_list'))


### BEFORE CHANGING PC DIRECTIONS
#CALCULATE MEAN PC SCORES BY SAMPLE FOR ALL CONFLUENCY GROUPS
pca_df_summary <- list()
for (l in 1:9){
  df <- pca_list[[l]]
  meta_df <- df[1:6]
  df <- cbind(meta_df['Sample'], df[7:20])
  df <- df %>% group_by(Sample) %>% summarise_at(vars(paste0('PC', 1:14)), mean)
  pca_df_summary[[l]] <- df
}

#CORRELATION BETWEEN MEAN PC SCORES WITH GSVA SCORES
pca_list_medPC_corGSVA <- list()
pca_list_medPC_corGSVA_pvalue <- list()
GSC.gsva <- GSC.gsva
sig <- rownames(GSC.gsva)

for (pc in 1:14){
  final_correlation_medPC <- as.data.frame(matrix(nrow=113))
  final_correlation_medPC_pvalue <- as.data.frame(matrix(nrow=113))
  
  for (con in 1:9){
    df <- pca_df_summary[[con]] ### select confluency group
    to_plot1 <- df[pc+1] %>% as.data.frame(.) ### select PC to focus on for the rest of the analysis
    gsva_short <- t(GSC.gsva) %>% as.data.frame(.) %>% .[df$Sample,]
    final_col <- data.frame()  
    for (x in 1:113){
      ctype <- sig[x] ### then select the cell types to focus
      #gsva_short <- gsva_short[,ctype]
      g <- gsva_short[,colnames(gsva_short) %in% ctype]
      #final <- data.frame()
      #to_plot2 <- g[col]
      cor <- cor.test(g, to_plot1[,1])
      #cor <- lm(to_plot1[,1]~to_plot2[,1])
      #r.sq <- summary(cor)$r.squared
      #m=cor$coefficients[2]
      p <- cor$estimate %>% round(.,3)
      d <- data.frame(label=ctype,correlation=p , p_value=cor$p.value)
      final_col <- rbind(final_col, d)
    }
    rownames(final_correlation_medPC_pvalue) <- rownames(final_correlation_medPC) <- final_col$label
    final_correlation_medPC <- cbind(final_correlation_medPC, final_col[2])
    final_correlation_medPC_pvalue <- cbind(final_correlation_medPC_pvalue, final_col[3])
  }
  final_correlation_medPC <- final_correlation_medPC[2:10]
  colnames(final_correlation_medPC) <- c(paste0('pixelBio', 1:9))
  
  pca_list_medPC_corGSVA[[pc]] <- final_correlation_medPC
  
  final_correlation_medPC_pvalue <- final_correlation_medPC_pvalue[2:10]
  colnames(final_correlation_medPC_pvalue) <- paste0('pixelBio', 1:9)
  pca_list_medPC_corGSVA_pvalue[[pc]] <- final_correlation_medPC_pvalue 
  
}

colnames(pc1_loadings) <- colnames(pc2_loadings) <- colnames(pc3_loadings) <- c('Texture',paste0("C", 1:9)) 


pc1_loadings$C5 <- pc1_loadings$C5 * (-1)
pc1_loadings$C6 <- pc1_loadings$C6 * (-1)
pc1_loadings$C7 <- pc1_loadings$C7 * (-1)
pc1_loadings$C8 <- pc1_loadings$C8 * (-1)
pc1_loadings$C9 <- pc1_loadings$C9 * (-1)

pc2_loadings$C2 <- pc2_loadings$C2 * (-1)
pc2_loadings$C3 <- pc2_loadings$C3 * (-1)
pc2_loadings$C4 <- pc2_loadings$C4 * (-1)
pc2_loadings$C6 <- pc2_loadings$C6 * (-1)
pc2_loadings$C9 <- pc2_loadings$C9 * (-1)


dir.create(paste0(tables_path, 'S_Table2'))

dir.create(paste0(tables_path, '/S_Table2'))
new_table_path <- paste0(tables_path, '/S_Table2/correlation_gsvaScoresWholeImages_PC')

for (cor in 1:14){
  write.csv(pca_list_medPC_corGSVA[[cor]], paste0(new_table_path, cor, ".csv"))  
}

dir.create(paste0(tables_path, 'S_Table3'))

### These are compiled into Table 1
write.csv(pc1_loadings, paste0(tables_path, "/S_Table3/Tbl_3_pc1_wholeImageLoadings.csv"))
write.csv(pc2_loadings, paste0(tables_path, "/S_Table3/Tbl_3_pc2_wholeImageLoadings.csv"))
write.csv(pc3_loadings, paste0(tables_path, "/S_Table3/Tbl_3_pc3_wholeImageLoadings.csv"))


#--------------------------------------------------------------------------------#
### Table end

save(pca_list, pca_list_medPC_corGSVA, pca_list_medPC_corGSVA_pvalue, pc1_loadings, pc2_loadings, pc3_loadings, file=paste0(path_to_repo, '/results/feature_correlations.Rdata'))
remove(list=ls())


