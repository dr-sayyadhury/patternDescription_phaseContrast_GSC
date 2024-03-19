

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

### FILE/FOLDER PATHS
###------------------------------------------------------------------#### 
path_to_repo <- '/Users/mystique27m/Documents/research/PostdoctoralResearch_2020/Projects/PatternRecognitionInGSCs_UsingCV/'
script_path <- paste0(path_to_repo,'scripts_final/code/')
out <- paste0(path_to_repo,'out_intermediate/')
figures_filepath <- paste0(path_to_repo, 'results/figures/')
tables_path <-  paste0(path_to_repo,'results/tables/')


### LOAD DATA
###------------------------------------------------------------------#### 
#source(paste0(script_path,"/script04b_.R"))
source(paste0(script_path,"source_scripts/sourceData_General_variables_and_functions.R"))
GSC.gsva <- read.csv(paste0(path_to_repo,'/results/GSC.gsva.csv'), row.names=1)
load(paste0(path_to_repo, '/results/feature_correlations.Rdata'))

rownames(pc2_loadings) <- gsub('_00', '', pc2_loadings$Texture)

pca_list_com <- do.call(rbind,pca_list)
pca_list_com_grn <- pca_list_com[c(1:5, 36:51)]
df_list_tp_subsetted <- readRDS(paste0(path_to_repo, '/results/df_list_tp_subsetted.rds'))

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
    #df <- apply(df,1, standardize)
    df <- cbind(meta,df) 
    df <- df %>% .[order(.$TimePt),]
    
    meta <- df[1:5]
    df <- df[6:features]  %>% t(.) %>% as.data.frame(.)
    rownames(meta) <- colnames(df)
    con_grp <- rep(paste0('C',i), nrow(meta))
    meta <- cbind(con_grp, meta)
    df <- cbind(meta, t(df))
    
    colnames(df) <- c('con_grp', shorter)
    new_df <- rbind(new_df, df)
    
  }
  new_df_list[[i]] <- new_df
}

new_df_list_comb <- do.call(rbind, new_df_list)
new_df_list_comb_grn <- new_df_list_comb[1:22]
new_df_list_comb_grn <- reshape2::melt(new_df_list_comb_grn, id.vars=colnames(new_df_list_comb_grn[1:6]))

descriptive <- c('Contrast_00', 'Correlation_00')
descriptive2 <- c('SumVariance_00', 'Variance_00', 'DifferenceVariance_00', 'InverseDifferenceMoment_00')

entropy <- c('Entropy_00', 'DifferenceEntropy_00', 'SumEntropy_00', 'SumAverage_00')
entropy2 <- c('Angular2ndMoment_00')

correlation <- c('InfoMeas1_00','InfoMeas2_00')

###------------------------------------------------------------------####
### SUPPLEMENTARY FIG 4A
pdf(paste0(figures_filepath, 'Sup_Fig_S4/S4_A.pdf'),  width=9.5, height=7.5)
ggline(new_df_list_comb_grn, x='con_grp', y='value', color='variable', add='mean_sd') +
  scale_color_manual(values=sequential_hcl('viridis', n=16)) +
  theme(legend.position = 'right')
dev.off()
###------------------------------------------------------------------####
### FIG END


ggplot(new_df_list_comb_grn, aes())

if(dir.exists(paste0(figures_filepath,'Sup_Fig_S5'))==F){
  dir.create(paste0(figures_filepath,'Sup_Fig_S5'))}

new_df_list_comb_txt <- new_df_list_comb[c(1:6,23:35)]
new_df_list_comb_txt <- reshape2::melt(new_df_list_comb_txt, id.vars=colnames(new_df_list_comb_txt[1:6]))

###------------------------------------------------------------------####
### SUPPLEMENTARY FIG 5A
###------------------------------------------------------------------####
pdf(paste0(figures_filepath, 'Sup_Fig_S5/S5_A.pdf'),  width=9.5, height=7.5)
ggline(new_df_list_comb_txt, x='con_grp', y='value', color='variable', add='mean_sd') +
  scale_color_manual(values=cols_vivid) +
  theme(legend.position = 'right')
dev.off()
###------------------------------------------------------------------####
### FIG END

new_df_list_comb_txt_rm_infomeas1 <- subset(new_df_list_comb_txt, variable!='InfoMeas1_00')


###------------------------------------------------------------------####
### SUPPLEMENTARY FIG 5B
###------------------------------------------------------------------####
pdf(paste0(figures_filepath, 'Sup_Fig_S5/S5_B.pdf'), width=9.5, height=7.5)
ggline(new_df_list_comb_txt_rm_infomeas1, x='con_grp', y='value', color='variable', add='mean_sd') +
  scale_color_manual(values=cols_vivid) +
  theme(legend.position = 'right')
dev.off()
###------------------------------------------------------------------####
### FIG END


pc2_loadings_melt <- reshape2::melt(pc2_loadings, id.vars='Texture')

max(pc2_loadings_melt$value)
steps_pos <- floor(max(pc2_loadings_melt$value)/0.05)

min(pc2_loadings_melt$value)
steps_neg <- floor(min(pc2_loadings_melt$value)/0.05)*(-1)


if(dir.exists(paste0(figures_filepath,'Figure_3'))==F){
  dir.create(paste0(figures_filepath,'Figure_3'))}

cols <- c(sequential_hcl('Magenta', n=steps_pos*3, rev=F), sequential_hcl('Mako', n=steps_neg*3, rev=T))

###------------------------------------------------------------------####
### FIGURE 3A
###------------------------------------------------------------------####
pdf(paste0(figures_filepath, 'Figure_3/Fig3_A.pdf'), width=13.5, height=3.5)
pheatmap::pheatmap(t(pc2_loadings[2:10]),
                   cluster_rows = F,
                   cluster_cols = F, 
                   color = rev(cols),
                   fontsize = 9,
                   angle_col = 45,
                   cellwidth = 30,
                   cellheight = 15,
                   border_color = NA)
dev.off()
###------------------------------------------------------------------####
### FIG END


cols_blues <- sequential_hcl('Mako', n=9, rev=T) ### low PC2
cols_mag <- sequential_hcl('Magenta', n=9, rev=T) ### PC2


glist <- list()
for (c in 1:9){
  cgroups <- levels(factor(pca_list_com$con_grp))
  df <- pca_list_com[pca_list_com$con_grp==cgroups[c],]
  glist[[c]] <- ggplot(df, aes(reorder(Sample,PC2), Angular2ndMoment_00, fill=con_grp)) +
    geom_violin(color='darkred', fill=cols_mag[c]) +
    theme_minimal() +
    theme(legend.position = 'none', 
          axis.title = element_blank())
  
}

###------------------------------------------------------------------####
### FIGURE 3B1
###------------------------------------------------------------------####
pdf(paste0(figures_filepath, 'Figure_3/Fig3_B1_ang.pdf'), width=6, height=12)
ggarrange(plotlist = glist, nrow=9)
dev.off()
###------------------------------------------------------------------####
### FIG END


### InfoMeas (low Pc2)
pca_list_com <- do.call(rbind, pca_list)
cols_blues <- sequential_hcl('Mako', n=9, rev=T)
glist <- list()
for (c in 1:9){
  cgroups <- levels(factor(pca_list_com$con_grp))
  df <- pca_list_com[pca_list_com$con_grp==cgroups[c],]
  glist[[c]] <- ggplot(df, aes(reorder(Sample,PC2), InfoMeas1_00, fill=con_grp)) +
            geom_violin(color='darkblue', fill=cols_blues[c]) +
    theme_minimal() +
    theme(legend.position = 'none', 
          axis.title = element_blank())
  
}

###------------------------------------------------------------------####
### FIGURE 3B2
pdf(paste0(figures_filepath, 'Figure_3/Fig3_B2_infomeas1.pdf'), width=6, height=12)
ggarrange(plotlist = glist, nrow=9)
dev.off()
###------------------------------------------------------------------####
### FIG END

glist <- list()
for (c in 1:9){
  cgroups <- levels(factor(pca_list_com$con_grp))
  df <- pca_list_com[pca_list_com$con_grp==cgroups[c],]
  glist[[c]] <- ggplot(df, aes(reorder(Sample,PC2), Granularity_16, fill=con_grp)) +
    geom_violin(color='darkblue', fill=cols_blues[c]) +
    theme_minimal() +
    theme(legend.position = 'none', 
          axis.title = element_blank())
  
}

###------------------------------------------------------------------####
### FIGURE 3C2
###------------------------------------------------------------------####
pdf(paste0(figures_filepath, 'Figure_3/Fig3_C2_grn16.pdf'), width=6, height=12)
ggarrange(plotlist = glist, nrow=9)
dev.off()
###------------------------------------------------------------------####
### FIG END


glist <- list()
for (c in 1:9){
  cgroups <- levels(factor(pca_list_com$con_grp))
  df <- pca_list_com[pca_list_com$con_grp==cgroups[c],]
  glist[[c]] <- ggplot(df, aes(reorder(Sample,PC2), Correlation_00, fill=con_grp)) +
    geom_violin(color='darkred', fill=cols_mag[c]) +
    theme_minimal() +
    theme(legend.position = 'none', 
          axis.title = element_blank())
  
}

###------------------------------------------------------------------####
### FIGURE 3C1
###------------------------------------------------------------------####
pdf(paste0(figures_filepath, 'Figure_3/Fig3_C1_correlation.pdf'), width=6, height=12)
ggarrange(plotlist = glist, nrow=9)
dev.off()
###------------------------------------------------------------------####
### FIG END

glist <- list()
for (c in 1:9){
  cgroups <- levels(factor(pca_list_com$con_grp))
  df <- pca_list_com[pca_list_com$con_grp==cgroups[c],]
  glist[[c]] <- ggplot(df, aes(reorder(Sample,PC2), DifferenceVariance_00, fill=con_grp)) +
    geom_violin(color='darkred', fill=cols_mag[c]) +
    theme_minimal() +
    theme(legend.position = 'none', 
          axis.title = element_blank())
  
}

###------------------------------------------------------------------####
### FIGURE 3D1
###------------------------------------------------------------------####
pdf(paste0(figures_filepath, 'Sup_Fig_S5/S5_D_diffV.pdf'), width=6, height=12)
ggarrange(plotlist = glist, nrow=9)
dev.off()
###------------------------------------------------------------------####
### FIG END

glist <- list()
for (c in 1:9){
  cgroups <- levels(factor(pca_list_com$con_grp))
  df <- pca_list_com[pca_list_com$con_grp==cgroups[c],]
  glist[[c]] <- ggplot(df, aes(reorder(Sample,PC2), Contrast_00, fill=con_grp)) +
    geom_violin(color='darkred', fill=cols_mag[c]) +
    theme_minimal() +
    theme(legend.position = 'none', 
          axis.title = element_blank())
  
}

###------------------------------------------------------------------####
### FIGURE 3D2
###------------------------------------------------------------------####
pdf(paste0(figures_filepath, 'Sup_Fig_S5/S5_D_contrast.pdf'), width=6, height=12)
ggarrange(plotlist = glist, nrow=9)
dev.off()
###------------------------------------------------------------------####
### FIG END


high_pc2 <-  subset(pca_list_com, Sample %in% c('G876', 'G895', 'G799'))[c(1:6,8,36:64)]
high_pc2_melt <- reshape2::melt(high_pc2, id.vars=colnames(high_pc2)[1:7])
high_pc2_melt$pc2_dir <- rep('cor', nrow(high_pc2_melt))

low_pc2 <-  subset(pca_list_com, Sample %in% c('G861', 'G564', 'G729'))[c(1:6,8,36:64)]
low_pc2_melt <- reshape2::melt(low_pc2, id.vars=colnames(low_pc2)[1:7])
low_pc2_melt$pc2_dir <- rep('anti_cor', nrow(low_pc2_melt))


pc2_topbot_3 <- rbind(high_pc2_melt, low_pc2_melt)
pc2_topbot_3_grn <- subset(pc2_topbot_3, variable %in% shorter[6:21])
pc2_topbot_3_txt <- subset(pc2_topbot_3, variable %in% shorter[22:34])

pc2_topbot_3_grn$con_grp <- factor(pc2_topbot_3_grn$con_grp, levels=paste0('C', 1:9))

pc2_topbot_3_txt$con_grp <- factor(pc2_topbot_3_txt$con_grp, levels=paste0('C', 1:9))
pc2_topbot_3_txt_disorder <- subset(pc2_topbot_3_txt, variable %in% entropy[1:3])
pc2_topbot_3_txt_hom <- subset(pc2_topbot_3_txt, variable %in% c('Ang2ndMoment','InverseDifferenceMoment', 'DifferenceVariance_00'))
pc2_topbot_3_grn_small <- subset(pc2_topbot_3_grn, variable %in% c('Granularity_1','Granularity_2', 'Granularity_3'))
pc2_topbot_3_grn_mid <- subset(pc2_topbot_3_grn, variable %in% c('Granularity_7','Granularity_8', 'Granularity_9'))
pc2_topbot_3_grn_high <- subset(pc2_topbot_3_grn, variable %in% c('Granularity_14','Granularity_15', 'Granularity_16'))


w=15
h=2.4

###------------------------------------------------------------------####
### FIGURE 3B
###------------------------------------------------------------------####
pdf(paste0(figures_filepath, 'Figure_3/Fig3_B_disorder.pdf'), width=w, height=h)
ggline(pc2_topbot_3_txt_disorder, x='con_grp', y='value', color='pc2_dir', add='mean_se', size=0.5) + 
  scale_color_manual(values=c(cols_blues[6], cols_mag[6])) +
  geom_violin(aes(fill=pc2_dir), position=position_dodge(width=0.3)) +
  scale_fill_manual(values=c(cols_blues[6], cols_mag[6])) +
  theme(legend.position = 'none',
        axis.line = element_line(color='grey'),
        axis.title = element_blank())
dev.off()
###------------------------------------------------------------------####
### FIG END

###------------------------------------------------------------------####
### FIGURE 3B
###------------------------------------------------------------------####
pdf(paste0(figures_filepath, 'Figure_3/Fig3_B_hom.pdf'), width=w, height=h)
ggline(pc2_topbot_3_txt_hom, x='con_grp', y='value', color='pc2_dir', add='mean_se', size=0.5) + 
  scale_color_manual(values=c(cols_blues[6], cols_mag[6])) +
  geom_violin(aes(fill=pc2_dir), position=position_dodge(width=0.3)) +
  scale_fill_manual(values=c(cols_blues[6], cols_mag[6])) +
  theme(legend.position = 'none',
        axis.line = element_line(color='grey'),
        axis.title = element_blank())
dev.off()
###------------------------------------------------------------------####
### FIG END

