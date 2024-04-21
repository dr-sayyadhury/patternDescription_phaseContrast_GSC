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


### source data 
###------------------------------------------------------------------#### 
#source(paste0(script_path,"/script04b_.R"))
source(paste0(script_path,"source_scripts/sourceData_General_variables_and_functions.R"))
GSC.gsva <- read.csv(paste0(path_to_repo,'/results/GSC.gsva.csv'), row.names=1)
load(paste0(out, 'feature_correlations.Rdata'))


if(dir.exists(paste0(figures_filepath,'Figure_2'))==F){
  dir.create(paste0(figures_filepath,'Figure_2'))}


#TO CHECK THE CONGRUENCY BETWEEN CONFLUENCY GROUPS
hm_col <- sequential_hcl(45, palette='Inferno', rev=T)
hm_col_neg<- sequential_hcl(45, palette='GnBu')

for (hm in 1:2){
  data <- pca_list_medPC_corGSVA[[hm]]
  meta <- paste0('C', 1:9)
  colnames(data) <- meta
  

###------------------------------------------------------------------###
# FIGURE 2A
###------------------------------------------------------------------###
  pdf(paste0(figures_filepath, 'Figure_2/Fig2_A_PC', hm, '.pdf'), width=4.5, height=18)
  heatmap_ <- pheatmap::pheatmap(data, 
                     #annotation_row = meta, 
                     #annotation_colors = list(cell_type_group=colours), 
                     show_rownames = F, 
                     annotation_legend = F,
                     legend = T,
                     border_color = NA, 
                     cluster_cols =F,
                     cluster_rows = T,
                     fontsize = 12,
                     color = hm_col
  )
  dev.off()
  ###------------------------------------------------------------------###
  # FIG END

  pdf(paste0(figures_filepath, 'Figure_2/Fig2_A_PC_labeled_', hm, '.pdf'), width=4.5, height=18)
  heatmap_ <- pheatmap::pheatmap(data, 
                     #annotation_row = meta, 
                     #annotation_colors = list(cell_type_group=colours), 
                     show_rownames = T, 
                     annotation_legend = F,
                     legend = T,
                     border_color = NA, 
                     cluster_cols =F,
                     cluster_rows = T,
                     fontsize = 12,
                     color = hm_col
  )
  dev.off()
  ###------------------------------------------------------------------###
  # FIG END


  new_df <- pca_list_medPC_corGSVA[[hm]]
  new_df$heatmap_order <- heatmap_$tree_row$order
  write.csv(new_df, paste0(tables_path, 'S_Table2/correlation_gsvaScoresWholeImages_heatmap_order_PC', hm))
  
  }

#THESE 14 PC1-PC14 CORRELATIONS WERE COMPILED INTO TABLE 2 TO SUPPORT THE ABOVE FIGURE 2B REPRESENTING ONLY PC2 CORRELATIONS
### PC1 correlations
pc1_correlation <- pca_list_medPC_corGSVA[[1]][1:9]
pc1_correlation <- apply(pc1_correlation,2, standardize) %>% as.data.frame(.)
pc1_correlation$sig <- rownames(pca_list_medPC_corGSVA[[1]])

pc1_correlation_p <- pca_list_medPC_corGSVA_pvalue[[1]]
pc1_correlation_p$sig <- rownames(pca_list_medPC_corGSVA_pvalue[[1]])

all(pc1_correlation$sig==rownames(pc1_correlation_p))

pc1_cor_melt <- reshape2::melt(pc1_correlation, id_var='sig')
colnames(pc1_cor_melt) <- c('sig','con_grp','correlation')


pc1_cor_p_melt <- reshape2::melt(pc1_correlation_p, id_var='sig')
colnames(pc1_cor_p_melt) <- c('sig','con_grp','p_value')

pc1_cor_melt$p_value <- pc1_cor_p_melt$p_value

pc1_cor_p_melt <- subset(pc1_cor_melt, p_value < 0.01)

pos <- subset(pc1_cor_p_melt, correlation>0)
pos$direction <- rep('positive_correlation', nrow(pos))

neg <- subset(pc1_cor_p_melt, correlation<0)
neg$direction <- rep('negative_correlation', nrow(neg))

pc1_cor_p_melt <- rbind(pos,neg)


###------------------------------------------------------------------###
# SUPPLEMENTARY FIGURE S3C
###------------------------------------------------------------------###
pdf(paste0(figures_filepath, 'Sup_Fig_S3/S3_C.pdf'), width=9, height=9)
ggplot(pc1_cor_melt, aes(x=reorder(sig, correlation), y=correlation, group=con_grp)) + 
  geom_line(size=0.45) +
  geom_point(data=pc1_cor_p_melt, aes(color=direction), size=3.6) +
  scale_color_manual(values=c('#0092C8','#DE1D82')) +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size=21),
        axis.title = element_text(size=24),
        #axis.line.x = element_blank(),
        legend.position = 'none') +
#  geom_hline(yintercept=0) +
  geom_vline(xintercept = 0) +
  xlab('Gene signatures') +
  ylab('Correlation score')
dev.off()
###------------------------------------------------------------------###
# FIG END


### PC2 correlations
pc2_correlation <- pca_list_medPC_corGSVA[[2]][1:9]
pc2_correlation <- apply(pc2_correlation,2, standardize) %>% as.data.frame(.)
pc2_correlation$sig <- rownames(pca_list_medPC_corGSVA[[2]])

pc2_correlation_p <- pca_list_medPC_corGSVA_pvalue[[2]]
pc2_correlation_p$sig <- rownames(pca_list_medPC_corGSVA_pvalue[[2]])

all(pc2_correlation$sig==rownames(pc2_correlation_p))

pc2_cor_melt <- reshape2::melt(pc2_correlation, id_var='sig')
colnames(pc2_cor_melt) <- c('sig','con_grp','correlation')


pc2_cor_p_melt <- reshape2::melt(pc2_correlation_p, id_var='sig')
colnames(pc2_cor_p_melt) <- c('sig','con_grp','p_value')

pc2_cor_melt$p_value <- pc2_cor_p_melt$p_value

pc2_cor_p_melt <- subset(pc2_cor_melt, p_value < 0.01)

pos <- subset(pc2_cor_p_melt, correlation>0)
pos$direction <- rep('positive_correlation', nrow(pos))

neg <- subset(pc2_cor_p_melt, correlation<0)
neg$direction <- rep('negative_correlation', nrow(neg))

pc2_cor_p_melt <- rbind(pos,neg)


###------------------------------------------------------------------###
# FIGURE 2B1
###------------------------------------------------------------------###

pdf(paste0(figures_filepath, 'Figure_2/Fig2_B1.pdf'), width=9, height=9)
ggplot(pc2_cor_melt, aes(x=reorder(sig, correlation), y=correlation, group=con_grp)) + 
  geom_line(size=0.45) +
  geom_point(data=pc2_cor_p_melt, aes(color=direction), size=3.6) +
  scale_color_manual(values=c('#0092C8','#DE1D82')) +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size=21),
        axis.title = element_text(size=24),
        #axis.line.x = element_blank(),
        legend.position = 'none') +
#  geom_hline(yintercept=0) +
  geom_vline(xintercept = 0) +
  xlab('Gene signatures') +
  ylab('Correlation score')
dev.off()
###------------------------------------------------------------------###
# FIG END


#SUBSETTING TOP AND BOTTOM 25% IMAGES ALONG PC2
g1 <- list()
g2 <- list()
g3 <- list()

df_high_combined <- data.frame()
df_low_combined <- data.frame()
df_rem_combined <- data.frame()

for (c in 1:9){
  
  pca_df <- pca_list[[c]]
  df_high <- subset(pca_df, subset=PC2 > 0.75*(max(PC2)-min(PC2))+min(PC2))
  df_low <- subset(pca_df, subset=PC2 < 0.25*(max(PC2)-min(PC2))+min(PC2))
  df_remaining <- subset(pca_df, subset= PC2 >= 0.25*(max(PC2)-min(PC2))+min(PC2))
  df_remaining <- subset(df_remaining, subset= PC2 <= 0.75*(max(PC2)-min(PC2))+min(PC2))
  
  g <- ggplot(pca_df, aes(PC2, PC1)) +
    geom_jitter(size=1.8, alpha=0.3, color='#B3AA99') +
    theme_minimal() +
    geom_jitter(data=df_high, color='#DE1D82') +
    geom_jitter(data=df_low, color='#0092C8') +
    theme(legend.position = 'none', 
          axis.title = element_text(size=12),
          axis.text.x = element_blank(),
          axis.text.y = element_text(angle=45, size=9),
          #axis.line.x = element_line(color='darkgray'),
          #axis.line.y = element_line(color='darkgray')
    ) +
    #geom_text_repel(aes(label=Sample), size=2.5) +
    xlab("PC2") + ylab("PC1")
  
  g2[[c]] <- g
  
  
  df_high_combined <- rbind(df_high_combined, df_high)
  df_low_combined <- rbind(df_low_combined, df_low)
  df_rem_combined <- rbind(df_rem_combined, df_remaining)
  
}

#SHOWING THE SAMPLES REPRESENTING THE TOP AND BOTTOM 20% AND THE COUNTS
df_high_combined$Group <- rep('high', nrow(df_high_combined))
df_low_combined$Group <- rep('low', nrow(df_low_combined))
df_rem_combined$Group <- rep('rem', nrow(df_rem_combined))

df_hl_combined <- rbind(df_high_combined, df_low_combined, df_rem_combined)

c <- dplyr::count(df_hl_combined, Group, Sample)
cc <- dplyr::count(df_hl_combined, Sample)

df_hl_combined_fraction <- data.frame()
df_hl_combined_fraction2 <- list()

for (g in 1:3){
  group <- levels(factor(c$Group))[g]
  df <- subset(c, Group==group)
  cg <- subset(cc, Sample %in% df$Sample)
  df$n <- df$n/cg$n
  df_hl_combined_fraction <- rbind(df_hl_combined_fraction, df)
  df_hl_combined_fraction2[[g]] <- df
}

df_hl_combined_fraction$Group <- factor(df_hl_combined_fraction$Group, levels = c('rem', 'high', 'low'))
new_order1 <- arrange(df_hl_combined_fraction2[[1]], n) %>% .$Sample 
new_order2 <- arrange(df_hl_combined_fraction2[[2]], n) %>% .$Sample %>% rev(.)
new_order <- c("G861", "G564","G729","G837", "G583", "G566","G549", "G797",  "G800","G851",   "G885","G523","G799", "G895", "G876")
df_hl_combined_fraction$Sample <- factor(df_hl_combined_fraction$Sample, levels = new_order)

###------------------------------------------------------------------###
# FIGURE 2C
###------------------------------------------------------------------###
pdf(paste0(figures_filepath,'Figure_2/Fig2_C.pdf'), width=3, height=3)
ggplot(df_hl_combined_fraction, aes(Sample, n, fill=Group)) + 
  geom_bar(stat='identity') + 
  theme_minimal() + 
  theme(legend.position = 'top',
        axis.text = element_text(size=9, angle=60, hjust=1),
        axis.title = element_blank()) +
  scale_fill_manual(values=c('lightgrey','#DE1D82','#0092C8')) 
#scale_fill_manual(values=cols_vivid) 
dev.off()
###---------------------------------------------------------------------#
#Fig end



#ASSIGNING GSVA SCORES TO THE SAMPLES SUBSAMPLED ABOVE
df_high_samples <- dplyr::count(df_high_combined, Sample) #%>% subset(n>10)
df_low_samples <- dplyr::count(df_low_combined, Sample) #%>% subset(n>10)
df_rem_samples <- dplyr::count(df_rem_combined, Sample)

GSC.gsva_up <- GSC.gsva[rownames(GSC.gsva) %in% rownames(pca_list_medPC_corGSVA[[1]]),] %>% as.data.frame(.)
all(rownames(GSC.gsva_up)==rownames(pca_list_medPC_corGSVA[[1]]))
GSC.gsva_up <- GSC.gsva_up[rownames(pca_list_medPC_corGSVA[[1]]),]
all(rownames(GSC.gsva_up)==rownames(pca_list_medPC_corGSVA[[1]]))


df_high_combined$Group <- rep('high', nrow(df_high_combined))
df_low_combined$Group <- rep('low', nrow(df_low_combined))
df_hl_comb <- rbind(df_high_combined, df_low_combined)
write.csv(df_hl_comb, paste0(out, 'df_top_bot_PC2_images_related_to_figure2E.csv'))

dplyr::count(df_hl_comb, Group, Sample) -> numbers

numbers$Sample <- factor(numbers$Sample, levels=new_order)


###------------------------------------------------------------------###
# FIGURE 2D
###------------------------------------------------------------------###
pdf(paste0(figures_filepath,'Figure_2/Fig2_D.pdf'), width=3, height=1.5)
ggplot(numbers, aes(Sample, n, fill=Group)) +
  geom_bar(stat='identity') +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, size=7, hjust=1),
        axis.text.y=element_text(size=6),
        axis.title = element_blank(),
        legend.position = 'none') +
  scale_fill_manual(values=c('#DE1D82','#0092C8')) 
dev.off()
###------------------------------------------------------------------###
# FIG END


remove(list=ls())


