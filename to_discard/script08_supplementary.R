### PACKAGES TO LOAD
###------------------------------------------------------------------#### 
library(ggplot2)  
library(colorspace)
library(tidyr)
library(dplyr)
library(ggthemes)
library(ggpubr)
library(ggrepel)
library(effectsize)
library(ggthemes)
library(scales)
library(forcats)
library(cowplot)


### FILE/FOLDER PATHS
###------------------------------------------------------------------#### 
#path_to_repo <- "/data/"
path_to_repo <- '/Users/tpugh_m04/Documents/Research/Projects/PatternRecognitionInGSCs_UsingCV'
#figures_filepath <-'/results/'
figures_filepath <- paste0(path_to_repo,'/results/')
#script_path <- '/code'
script_path <- paste0(path_to_repo,'/scripts')
#results_path <- '/results/S3_outPutDataFiles/'
results_path <- paste0(figures_filepath,'/S3_outPutDataFiles/')
dir.create(results_path)
#tables_path <- '/results/tables'
tables_path <-  paste0(figures_filepath,'/tables/')

dir.create(tables_path)



### IMPORT DATA FROM S1 & S3
###------------------------------------------------------------------#### 
# S1
#phase <- read.csv('/results/S1_outPutDataFiles/final_dataset_afterQC.csv', row.names=1)
phase <- read.csv(paste0(path_to_repo, '/results/S1_outPutDataFiles/final_dataset_afterQC.csv'), row.names = 1)

# S3
GSC.gsva <- read.csv(paste0(path_to_repo,'/results/S3_outPutDataFiles/GSC.gsva.csv'), row.names=1)
#GSC.gsva <- read.csv('/results/S3_outPutDataFiles/GSC.gsva.csv', row.names=1)

# general script
source(paste0(script_path,"/sourceData_General_variables_and_functions.R"))




#--------------------------------------------------------------------------------#
### --------------------------------- CODE START ------------------------------###
#--------------------------------------------------------------------------------#

phase <- phase[,shorter]

df1 <- phase[!(phase$Sample %in% c('G566','G583')),] %>% as.data.frame(.)
### samples imaged at 12hrly intervals so this was processed separately.
df2 <- phase[(phase$Sample %in% c('G566','G583')),] %>% as.data.frame(.)

### Prepare cellprofiler data for normalization 
df_list_tp_strictArea <- list()

for (confluency in 1:9){
  ### subset data from phase dataset and then perform z-normalization
  tp_AreaMax <- data.frame()
  for (c in 1:length(levels(factor(df1$Sample)))){
    S <- levels(factor(df1$Sample))[c]
    sample <- df1[df1$Sample==S,]
    #tp <- sample %>% group_by(TimePt) %>% summarise(Avg=median(Area))
    t <- subset(sample, Area > confluency*0.1 & Area < (confluency+1)*0.1)
    #s <- subset(sample, TimePt %in% t$TimePt)
    tp_AreaMax <- rbind(tp_AreaMax, t)
  }
  
  for (c in 1:length(levels(factor(df2$Sample)))){
    S <- levels(factor(df2$Sample))[c]
    sample <- df2[df2$Sample==S,]
    #tp <- sample %>% group_by(TimePt) %>% summarise(Avg=median(Area))
    t <- subset(sample, Area > confluency*0.1 & Area < (confluency+1)*0.1)
    #s <- subset(sample, TimePt %in% t$TimePt)
    tp_AreaMax <- rbind(tp_AreaMax, t)
  }
  
  df_list_tp_strictArea[[confluency]] <- tp_AreaMax
  
}

#boxplots to show the pixel-area occupancy of each group
graphs <- list()

for (i in 1:length(df_list_tp_strictArea)){
  x <- df_list_tp_strictArea[[i]]
  g <- ggplot(x, aes(x=Sample, y=Area, fill=Sample)) +
    scale_color_discrete_qualitative(palette='Dark 2') +
    geom_boxplot() +
    geom_jitter(size=0.3, alpha=0.3) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.1)) +
    theme_pubr() + 
    labs(title=paste0("C",i)) +
    theme(axis.text.x=element_blank(),
          legend.position='none', 
          plot.title = element_text(size=18, color="darkblue"),
          axis.text = element_text(size=12),
          axis.title = element_blank(),
          axis.ticks.x = element_blank())
  
  graphs[[i]] <- g
}



###------------------------------------------------------------------###
# SUPPLEMENTARY FIGURE 8A #
###------------------------------------------------------------------###
if(dir.exists(paste0(figures_filepath, 'Sup_Fig_S8'))==F){
  dir.create(paste0(figures_filepath, 'Sup_Fig_S8'))
}

pdf(file = paste0(figures_filepath,'Sup_Fig_S8/S8_A.pdf'), width=18, height=3)
suppressWarnings(print(ggpubr::ggarrange(plotlist = graphs, nrow=1, ncol=9)))
dev.off()
###------------------------------------------------------------------###
#Fig End



for (l in 1:9){
  df <- df_list_tp_strictArea[[l]]
  print(dplyr::count(df, Sample))
  
}

df_list_tp_strictArea[[1]] <- subset(df_list_tp_strictArea[[1]], Sample!='G895')
df_list_tp_strictArea[[1]] <- subset(df_list_tp_strictArea[[1]], Sample!='G800')
df_list_tp_strictArea[[2]] <- subset(df_list_tp_strictArea[[2]], Sample!='G797')

for (l in 1:9){
  df <- df_list_tp_strictArea[[l]]
  print(dplyr::count(df, Sample))
  
}

#NORMALIZATION AND STANDARDIZATION BY CONFLUENCY GROUPS
pt.size=1.5
features=34
new_df_list <- list()
for (i in 1:9){
  d <- df_list_tp_strictArea[[i]]
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
    #rownames(meta) <- colnames(df)
    con_grp <- rep(paste0('C',i), nrow(meta))
    meta <- cbind(con_grp, meta)
    df <- cbind(meta, t(df))
    
    colnames(df) <- c('con_grp', colnames(d))
    new_df <- rbind(new_df, df)
    
  }
  new_df_list[[i]] <- new_df
}



features=35
pca_list <- list()
pca_loadings_list <- list()
pc1_loadings <- pc2_loadings <- pc3_loadings <- shorter[6:34] %>% as.data.frame(.)
g1 <- list()
g2 <- list()
sc <- list()
pt.size=1.5

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
  
  sc[[i]] <- PCAtools::screeplot(p, colBar='darkblue', axisLabSize = 9, components = 1:(pc1+1), title = '', vline=c(var,pc1))
  pca_df <- cbind(meta,p$rotated)
  colnames(pca_df) <- c('con_grp','Media','Sample','TimePt','Well','Area', paste0('PC', 1:ncol(p$rotated)))
  
  pc1_loadings <- cbind(pc1_loadings, p[["loadings"]][["PC1"]])
  pc2_loadings <- cbind(pc2_loadings, p[["loadings"]][["PC2"]])
  pc3_loadings <- cbind(pc3_loadings, p[["loadings"]][["PC3"]])
  
  #common <- intersect(DI_grad_bar_withD_IR$SAMPLE_ID, levels(factor(meta$Sample)))
  #col_log <- c(cols_vivid, col_bold)[DI_grad_bar_withD_IR$SAMPLE_ID %in% common]
  #DI <- DI_grad_bar_withD_IR[DI_grad_bar_withD_IR$SAMPLE_ID %in% common,]
  #DI <- DI[rep(1:nrow(DI),dplyr::count(meta, Sample)$n),]
  
  #pca_df <- cbind(pca_df, DI[2:6])
  pca_df <- cbind(pca_df,data)
  pca_list[[i]] <- pca_df
  
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





###------------------------------------------------------------------###
# SUPPLEMENTARY FIGURE 8B-C #
###------------------------------------------------------------------###

pdf(file = paste0(figures_filepath,'Sup_Fig_S8/S8_B.pdf'), width=21, height=3.5)
ggarrange(plotlist=g1, nrow=1, ncol=9)
dev.off()

pdf(file = paste0(figures_filepath,'Sup_Fig_S8/S8_C.pdf'), width=21, height=3.5)
print(ggarrange(plotlist=sc, nrow=1, ncol=9))
dev.off()
###------------------------------------------------------------------###
### Fig ends



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

for (pc in 1:14){
  final_correlation_medPC <- as.data.frame(matrix(nrow=102))
  final_correlation_medPC_pvalue <- as.data.frame(matrix(nrow=102))
  
  for (con in 1:9){
    df <- pca_df_summary[[con]] ### select confluency group
    to_plot1 <- df[pc+1] %>% as.data.frame(.) ### select PC to focus on for the rest of the analysis
    gsva_short <- t(GSC.gsva) %>% as.data.frame(.) %>% .[df$Sample,]
    final_col <- data.frame()  
    for (x in 1:14){
      ctype <- cell_types_up[[x]] ### then select the cell types to focus
      #gsva_short <- gsva_short[,ctype]
      g <- gsva_short[,colnames(gsva_short) %in% ctype]
      final <- data.frame()
      for (col in 1:ncol(g)){
        to_plot2 <- g[col]
        cor <- lm(to_plot1[,1]~to_plot2[,1])
        r.sq <- summary(cor)$r.squared
        m=cor$coefficients[2]
        p <- summary(cor)$coefficients[2,4] %>% round(.,3)
        d <- data.frame(label=names(to_plot2),r.squared=round(r.sq, digits=3)*(m/abs(m)), p_value=p,  col=color_cellCLass_up[x], cell_type_group=cell_type_names_up[x])
        final <- rbind(final,d)
      }
      final_col <- rbind(final_col, final)
    }
    rownames(final_correlation_medPC_pvalue) <- rownames(final_correlation_medPC) <- final_col$label
    final_correlation_medPC <- cbind(final_correlation_medPC, final_col[2])
    final_correlation_medPC_pvalue <- cbind(final_correlation_medPC_pvalue, final_col[3])
  }
  final_correlation_medPC <- cbind(final_correlation_medPC, final_col[c(4,5)])
  final_correlation_medPC <- final_correlation_medPC[2:12]
  colnames(final_correlation_medPC) <- c(paste0('C', 1:9), 'col', 'cell_type_group')
  
  pca_list_medPC_corGSVA[[pc]] <- final_correlation_medPC
  
  final_correlation_medPC_pvalue <- final_correlation_medPC_pvalue[2:10]
  colnames(final_correlation_medPC_pvalue) <- paste0('C', 1:9)
  pca_list_medPC_corGSVA_pvalue[[pc]] <- final_correlation_medPC_pvalue 
  
}



### FIGURES
###------------------------------------------------------------------###
# SUPPLEMENTARY FIG 8D
###------------------------------------------------------------------###

#TO CHECK THE CONGRUENCY BETWEEN CONFLUENCY GROUPS
hm_col <- sequential_hcl(45, palette='Inferno', rev=T)
hm_col_neg<- sequential_hcl(45, palette='GnBu')

for (hm in 1:3){
  data <- pca_list_medPC_corGSVA[[hm]]
  colours <- data$col
  #names(colours) <- data
  meta <- data[10]
  data.cor <- cor(data[1:9])
  
  pdf(paste0(figures_filepath, 'Sup_Fig_S8/S8_D', hm, '.pdf'), width=6.6, height=6)
  pheatmap::pheatmap(data.cor, 
                     #annotation_row = meta, 
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
#--------------------------------------------------------------------------------#
#Fig end



### FIGURES
###------------------------------------------------------------------###
# SUPPLEMENTARY FIG 8E
###------------------------------------------------------------------###
#ONLY PLOT CORRELATIONS ALONG PC2 FOR ALL CONFLUENCY GROUPS
glist <- list()
for (l in 1:3){
  data <- pca_list_medPC_corGSVA[[l]]
  data <- cbind(rownames(data), data)
  colnames(data) <- c('cell-types', paste0('C', 1:9), 'col', 'Class')
  
  dp <- pca_list_medPC_corGSVA_pvalue[[l]]
  dp <- reshape2::melt(dp[2:8])
  colnames(dp) <- c('c', 'p')
  
  data_c1_c9 <- data[c(1,2,10)]
  data <- data[c(1,3:9)]
  
  data <- reshape2::melt(data)
  data_c1_c9 <- reshape2::melt(data_c1_c9)
  
  data <- cbind(data, dp$p)
  data_c1_c9$p <- rep('p', nrow(data_c1_c9)) 
  colnames(data_c1_c9) <- c(colnames(data_c1_c9)[1:3], 'p')
  colnames(data) <- c(colnames(data)[1:3], 'p')
  
  data_pos <- subset(data, p<0.05 & value > 0)
  data_neg <- subset(data, p<0.05 & value < 0)
  
  #  names <- c('S8_E.pdf', 'S8_F.pdf', 'S8_G.pdf')
  #  fig_name <- paste0('Sup_Fig_S8/S8_E', l, '.pdf')
  glist[[l]] <- ggplot(data, aes(reorder(`cell-types`, value), value, group=variable, colour=variable)) + 
    geom_hline(yintercept = 0, color='lightgrey') +
    #geom_line(data=data_c1_c9, color='lightgrey', linewidth=.9) +
    geom_line() +
    scale_color_manual(values=cols_vivid) +
    geom_point(data=data_pos,color='#DE1D82', size=2.1) +
    geom_point(data=data_neg, color='#0092C8', size=2.1) +
    theme_tufte() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_blank())
  
  
}


pdf(paste0(figures_filepath, "Sup_Fig_S8/S8_E.pdf") , width=21, height=4.5)
ggarrange(plotlist=glist, ncol=3)
dev.off()
#--------------------------------------------------------------------------------#
#Fig end




remove(list=ls())

