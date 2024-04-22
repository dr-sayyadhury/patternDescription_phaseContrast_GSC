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
load(paste0(out, '/feature_correlations.Rdata'))

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
    new_df$InfoMeas1_00 <- new_df$InfoMeas1_00 * -1
    
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
  theme(legend.position = 'right',
  legend.text = element_text(size=12),
  legend.title = element_blank())
dev.off()
###------------------------------------------------------------------####
### FIG END



new_df_list_comb_txt <- new_df_list_comb[c(1:6,23:35)]
new_df_list_comb_txt <- reshape2::melt(new_df_list_comb_txt, id.vars=colnames(new_df_list_comb_txt[1:6]))

###------------------------------------------------------------------####
### SUPPLEMENTARY FIG 4
###------------------------------------------------------------------####
pdf(paste0(figures_filepath, 'Sup_Fig_S4/S4_B.pdf'),  width=9.5, height=7.5)
ggline(new_df_list_comb_txt, x='con_grp', y='value', color='variable', add='mean_sd') +
  scale_color_manual(values=cols_vivid) +
  theme(legend.position = 'right',
  legend.text = element_text(size=12),
  legend.title = element_blank())
dev.off()

###------------------------------------------------------------------####
### FIG END

new_df_list_comb_txt_rm_infomeas1 <- subset(new_df_list_comb_txt, variable!='InfoMeas1_00')


###------------------------------------------------------------------####
### SUPPLEMENTARY FIG 5B
###------------------------------------------------------------------####
pdf(paste0(figures_filepath, 'Sup_Fig_S4/S4_C.pdf'), width=9.5, height=7.5)
ggline(new_df_list_comb_txt_rm_infomeas1, x='con_grp', y='value', color='variable', add='mean_sd') +
  scale_color_manual(values=cols_vivid) +
  theme(legend.position = 'right',
  legend.text = element_text(size=12),
  legend.title = element_blank())
dev.off()
###------------------------------------------------------------------####
### FIG END


pc2_loadings_melt <- reshape2::melt(pc2_loadings, id.vars='Texture')
max(pc2_loadings_melt$value)
steps_pos <- floor(max(pc2_loadings_melt$value)/0.05)
min(pc2_loadings_melt$value)
steps_neg <- floor(min(pc2_loadings_melt$value)/0.05)*(-1)
cols <- c(sequential_hcl('Magenta', n=steps_pos*3, rev=F), sequential_hcl('Mako', n=steps_neg*3, rev=T))


if(dir.exists(paste0(figures_filepath,'Figure_3'))==F){
  dir.create(paste0(figures_filepath,'Figure_3'))}


###------------------------------------------------------------------####
### FIGURE 3A
###------------------------------------------------------------------####

pdf(paste0(figures_filepath, 'Figure_3/Fig3_A.pdf'), width=4.0, height=12)
pheatmap::pheatmap((pc2_loadings[2:10]),
                   cluster_rows = F,
                   cluster_cols = F, 
                   color = rev(cols),
                   fontsize = 9,
                   angle_row = 45,
                   cellwidth = 15,
                   cellheight = 27,
                   border_color = NA)
dev.off()
###------------------------------------------------------------------####
### FIG END
sample_cols <- c(
'#ebac23', #yellow
'#b80058', #lipstick
#'#008cf9', #azure
#'#ff9287', #coral
'#56641a', #fern frond
'#c0affb', #perfume
'#e6a176', #apricot
'#00678a', #orient
'#984464', #vin rouge
'#5eccab', #downy
'#cdcdcd', #gray
'#274d52', #plantation
'#c7a2a6', #eunry
'#818b70', #battleship
'#604e3c', #kabul
'#8c9fb7', #bali hai
'#796880') #rum



#cols_blues <- sequential_hcl('Mako', n=9, rev=T) ### low PC2
#cols_mag <- sequential_hcl('Magenta', n=9, rev=T) ### PC2
#cols_neutral <- sequential_hcl('Grays', n=9, rev=T) ### high PC2
#sample_cols <- divergingx_hcl('Zissou 1', n=15)
#all_sample_means <- aggregate(PC2 ~ Sample, data = pca_list_com, mean) %>% arrange(PC2)
#all_sample_means$Sample <- factor(all_sample_means$Sample, levels=all_sample_means$Sample)
all_sample_cols <- data.frame(Sample=unique(pca_list_com$Sample), colors=sample_cols)
rownames(all_sample_cols) <- all_sample_cols$Sample
### plot tje colors from all sample means and label with the color name
all_sample_cols$scale <- 3
colors <- all_sample_cols$colors
names(colors) <- all_sample_cols$Sample

pdf(paste0(figures_filepath, 'Figure_3/Fig3_B_sample_colors.pdf'), width=9, height=1)
ggplot(all_sample_cols, aes(x=Sample, y=scale, fill=Sample)) +
  geom_tile() +
  scale_fill_manual(values=df_cols$colors) +
  theme_minimal() +
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text.x = element_text(angle=90, size=18),
        axis.text.y = element_blank())
dev.off()

###------------------------------------------------------------------####
### Supplementary figure 5A - low confluency
###------------------------------------------------------------------####

### Neurodevelopmental samples
pca_list_com <- do.call(rbind, pca_list)

### granularity 3
glist <- list()
for (c in 1:9){
  cgroups <- levels(factor(pca_list_com$con_grp))
  df <- pca_list_com[pca_list_com$con_grp == cgroups[c],]

  sample_means <- aggregate(Granularity_3 ~ Sample, data = df, mean)
  pc2_means <- aggregate(PC2 ~ Sample, data = df, mean)
  if (all(sample_means$Sample==pc2_means$Sample)==T){
    sample_means$PC2 <- pc2_means$PC2
  }

  sample_means <- sample_means[order(sample_means$PC2),]
 
  cols <- all_sample_cols[sample_means$Sample,]
  cols <- cols[sample_means$Sample,]
  cols$PC2 <- sample_means$PC2
  cols <- cols[order(cols$PC2),]

  cor <- cor.test(sample_means$Granularity_3, sample_means$PC2)
  p_value <- cor$p.value %>% round(.,6)
  r <- (cor$estimate) %>% round(., 3)


if (p_value < 0.001) {
  p <- '***'  # Assign '***' for p-values less than 0.001
} else if (p_value < 0.01) {
  p <- '**'   # Assign '**' for p-values less than 0.01
} else {
  p <- ''     # No asterisks for p-values 0.01 or higher
}



  # Create a ggplot
  glist[[c]] <- ggplot() +
    geom_line(data = sample_means, aes(x = reorder(Sample, PC2), y = Granularity_3, group=1), color = 'black', size = 0.9) +  # Line connecting means
    geom_jitter(data = df, aes(x = Sample, y = Granularity_3, color=Sample), size=0.3) +
    geom_line(data = sample_means, aes(x = Sample, y = Granularity_3, group=1), color = 'black', size = 0.9) +  # Line connecting means
    scale_color_manual(values=colors, breaks = sample_means$Sample) +
    annotate("text", x = Inf, y = Inf, 
         label =  p, 
         #label = formatC(p_value, format = "f", digits = 3),
         #label = paste(sprintf("cor = %.2f,", r), "p =", formatC(p_value, format = "f", digits = 3)), 
         hjust = 1.1, vjust = 1.1, size = 6, color='red') +
    theme_minimal() +
    theme(
  legend.position = 'right', 
  legend.text = element_text(size=9),
  legend.title = element_blank(),
  legend.spacing.y = unit(0.1, 'cm'),  # Reduced spacing
  legend.key.height = unit(0.1, 'cm'),  # Smaller keys
  axis.title = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size=12)
) +   
    guides(color = guide_legend(override.aes = list(size = 3))) +
    guides(fill = guide_legend(byrow = TRUE)) +
    ylim((min(df$Granularity_3)), (max(df$Granularity_3)+0.15))

}



### FIGURE 3B ###
pdf(paste0(figures_filepath, 'Figure_3/Fig3_B1_grn3.pdf'), width = 4.5, height = 18)
ggarrange(plotlist = glist, nrow = 9)
dev.off()



###------------------------------------------------------------------####
### Supplementary figure 5B - low confluency
###------------------------------------------------------------------####


### granularity 3
glist <- list()
for (c in 1:9){
  cgroups <- levels(factor(pca_list_com$con_grp))
  df <- pca_list_com[pca_list_com$con_grp == cgroups[c],]

  sample_means <- aggregate(Granularity_8 ~ Sample, data = df, mean)
  pc2_means <- aggregate(PC2 ~ Sample, data = df, mean)
  if (all(sample_means$Sample==pc2_means$Sample)==T){
    sample_means$PC2 <- pc2_means$PC2
  }

  sample_means <- sample_means[order(sample_means$PC2),]
 
  cols <- all_sample_cols[sample_means$Sample,]
  cols <- cols[sample_means$Sample,]
  cols$PC2 <- sample_means$PC2
  cols <- cols[order(cols$PC2),]

  cor <- cor.test(sample_means$Granularity_8, sample_means$PC2)
  p_value <- cor$p.value
  r <- (cor$estimate) %>% round(., 3)


if (p_value < 0.001) {
  p <- '***'  # Assign '***' for p-values less than 0.001
} else if (p_value < 0.01) {
  p <- '**'   # Assign '**' for p-values less than 0.01
} else {
  p <- ''     # No asterisks for p-values 0.01 or higher
}



  # Create a ggplot
  glist[[c]] <- ggplot() +
    geom_line(data = sample_means, aes(x = reorder(Sample, PC2), y = Granularity_8, group=1), color = 'black', size = 0.9) +  # Line connecting means
    geom_jitter(data = df, aes(x = Sample, y = Granularity_8, color=Sample), size=0.3) +
    geom_line(data = sample_means, aes(x = Sample, y = Granularity_8, group=1), color = 'black', size = 0.9) +  # Line connecting means
    scale_color_manual(values=colors, breaks = sample_means$Sample) +
    annotate("text", x = Inf, y = Inf, 
         label =  p, 
        #label = paste(sprintf("cor = %.2f,", r), "p =", formatC(p_value, format = "e", digits = 2)), 
        hjust = 1.1, vjust = 1.1, size = 9, color='red') +
    theme_minimal() +
    theme(
  legend.position = 'right', 
  legend.text = element_text(size=9),
  legend.title = element_blank(),
  legend.spacing.y = unit(0.1, 'cm'),  # Reduced spacing
  legend.key.height = unit(0.1, 'cm'),  # Smaller keys
  axis.title = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size=12)
) +   
    guides(color = guide_legend(override.aes = list(size = 3))) +
    guides(fill = guide_legend(byrow = TRUE)) +


    ylim((min(df$Granularity_8)), (max(df$Granularity_8)+0.15))

}


### FIGURE 3B ###
pdf(paste0(figures_filepath, 'Figure_3/Fig3_B1_grn8.pdf'), width = 4.5, height = 18)
ggarrange(plotlist = glist, nrow = 9)
dev.off()

###------------------------------------------------------------------####
### Supplementary figure 5C - low confluency
###------------------------------------------------------------------####

### Neurodevelopmental samples
### Neurodevelopmental samples

### granularity 3
glist <- list()
for (c in 1:9){
  cgroups <- levels(factor(pca_list_com$con_grp))
  df <- pca_list_com[pca_list_com$con_grp == cgroups[c],]

  sample_means <- aggregate(InfoMeas1_00 ~ Sample, data = df, mean)
  pc2_means <- aggregate(PC2 ~ Sample, data = df, mean)
  if (all(sample_means$Sample==pc2_means$Sample)==T){
    sample_means$PC2 <- pc2_means$PC2
  }

  sample_means <- sample_means[order(sample_means$PC2),]
 
  cols <- all_sample_cols[sample_means$Sample,]
  cols <- cols[sample_means$Sample,]
  cols$PC2 <- sample_means$PC2
  cols <- cols[order(cols$PC2),]

  cor <- cor.test(sample_means$InfoMeas1_00, sample_means$PC2)
  p_value <- cor$p.value %>% round(.,6)
  r <- (cor$estimate) %>% round(., 3)

if (p_value < 0.001) {
  p <- '***'  # Assign '***' for p-values less than 0.001
} else if (p_value < 0.01) {
  p <- '**'   # Assign '**' for p-values less than 0.01
} else {
  p <- ''     # No asterisks for p-values 0.01 or higher
}

  # Create a ggplot
  glist[[c]] <- ggplot() +
    geom_line(data = sample_means, aes(x = reorder(Sample, PC2), y = InfoMeas1_00, group=1), color = 'black', size = 0.9) +  # Line connecting means
    geom_jitter(data = df, aes(x = Sample, y = InfoMeas1_00, color=Sample), size=0.3) +
    geom_line(data = sample_means, aes(x = Sample, y = InfoMeas1_00, group=1), color = 'black', size = 0.9) +  # Line connecting means
    scale_color_manual(values=colors, breaks = sample_means$Sample) +
    annotate("text", x = Inf, y = Inf, 
         label =  p, 
        #label = paste(sprintf("cor = %.2f,", r), "p =", formatC(p_value, format = "e", digits = 2)), 
        hjust = 1.1, vjust = 1.1, size = 9, color='red') +
    theme_minimal() +
    theme(
  legend.position = 'right', 
  legend.text = element_text(size=9),
  legend.title = element_blank(),
  legend.spacing.y = unit(0.1, 'cm'),  # Reduced spacing
  legend.key.height = unit(0.1, 'cm'),  # Smaller keys
  axis.title = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size=12)
) +   
    guides(color = guide_legend(override.aes = list(size = 3))) +
    guides(fill = guide_legend(byrow = TRUE)) +


    ylim((min(df$InfoMeas1_00)), (max(df$InfoMeas1_00)+0.15))

}


### FIGURE 3B ###
pdf(paste0(figures_filepath, 'Figure_3/Fig3_B1_IMC1.pdf'), width = 4.5, height = 18)
ggarrange(plotlist = glist, nrow = 9)
dev.off()


###------------------------------------------------------------------####
### Figure 3 - low confluency
###------------------------------------------------------------------####

glist <- list()
### Neurodevelopmental samplesglist <- list()
for (c in 1:9){
  cgroups <- levels(factor(pca_list_com$con_grp))
  df <- pca_list_com[pca_list_com$con_grp == cgroups[c],]

  sample_means <- aggregate(InfoMeas2_00 ~ Sample, data = df, mean)
  pc2_means <- aggregate(PC2 ~ Sample, data = df, mean)
  if (all(sample_means$Sample==pc2_means$Sample)==T){
    sample_means$PC2 <- pc2_means$PC2
  }

  sample_means <- sample_means[order(sample_means$PC2),]
 
  cols <- all_sample_cols[sample_means$Sample,]
  cols <- cols[sample_means$Sample,]
  cols$PC2 <- sample_means$PC2
  cols <- cols[order(cols$PC2),]

  cor <- cor.test(sample_means$InfoMeas2_00, sample_means$PC2)
  p_value <- cor$p.value %>% round(.,6)
  r <- (cor$estimate) %>% round(., 3)

if (p_value < 0.001) {
  p <- '***'  # Assign '***' for p-values less than 0.001
} else if (p_value < 0.01) {
  p <- '**'   # Assign '**' for p-values less than 0.01
} else {
  p <- ''     # No asterisks for p-values 0.01 or higher
}


  # Create a ggplot
  glist[[c]] <- ggplot() +
    geom_line(data = sample_means, aes(x = reorder(Sample, PC2), y = InfoMeas2_00, group=1), color = 'black', size = 0.9) +  # Line connecting means
    geom_jitter(data = df, aes(x = Sample, y = InfoMeas2_00, color=Sample), size=0.3) +
    geom_line(data = sample_means, aes(x = Sample, y = InfoMeas2_00, group=1), color = 'black', size = 0.9) +  # Line connecting means
    scale_color_manual(values=colors, breaks = sample_means$Sample) +
    annotate("text", x = Inf, y = Inf, 
        label = p, 
        hjust = 1.1, vjust = 1.1, size = 9, color='red') +
    theme_minimal() +
    theme(
  legend.position = 'right', 
  legend.text = element_text(size=9),
  legend.title = element_blank(),
  legend.spacing.y = unit(0.1, 'cm'),  # Reduced spacing
  legend.key.height = unit(0.1, 'cm'),  # Smaller keys
  axis.title = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size=12)
) +   
    guides(color = guide_legend(override.aes = list(size = 3))) +
    guides(fill = guide_legend(byrow = TRUE)) +


    ylim((min(df$InfoMeas2_00)), (max(df$InfoMeas2_00)+0.15))

}


### FIGURE 3B ###
pdf(paste0(figures_filepath, 'Figure_3/Fig3_B1_IMC2.pdf'), width = 4.5, height = 18)
ggarrange(plotlist = glist, nrow = 9)
dev.off()


###------------------------------------------------------------------####
### Figure 3 - low confluency
###------------------------------------------------------------------####

### Neurodevelopmental samples
pca_list_com <- do.call(rbind, pca_list)
### granularity 3### Neurodevelopmental samplesglist <- list()
glist <- list()
for (c in 1:9){
  cgroups <- levels(factor(pca_list_com$con_grp))
  df <- pca_list_com[pca_list_com$con_grp == cgroups[c],]

  sample_means <- aggregate(InverseDifferenceMoment_00 ~ Sample, data = df, mean)
  pc2_means <- aggregate(PC2 ~ Sample, data = df, mean)
  if (all(sample_means$Sample==pc2_means$Sample)==T){
    sample_means$PC2 <- pc2_means$PC2
  }

  sample_means <- sample_means[order(sample_means$PC2),]
 
  cols <- all_sample_cols[sample_means$Sample,]
  cols <- cols[sample_means$Sample,]
  cols$PC2 <- sample_means$PC2
  cols <- cols[order(cols$PC2),]

  cor <- cor.test(sample_means$InverseDifferenceMoment_00, sample_means$PC2)
  p_value <- cor$p.value %>% round(.,6)
  r <- (cor$estimate) %>% round(., 3)

if (p_value < 0.001) {
  p <- '***'  # Assign '***' for p-values less than 0.001
} else if (p_value < 0.01) {
  p <- '**'   # Assign '**' for p-values less than 0.01
} else {
  p <- ''     # No asterisks for p-values 0.01 or higher
}


  # Create a ggplot
  glist[[c]] <- ggplot() +
    geom_line(data = sample_means, aes(x = reorder(Sample, PC2), y = InverseDifferenceMoment_00, group=1), color = 'black', size = 0.9) +  # Line connecting means
    geom_jitter(data = df, aes(x = Sample, y = InverseDifferenceMoment_00, color=Sample), size=0.3) +
    geom_line(data = sample_means, aes(x = Sample, y = InverseDifferenceMoment_00, group=1), color = 'black', size = 0.9) +  # Line connecting means
    scale_color_manual(values=colors, breaks = sample_means$Sample) +
    annotate("text", x = Inf, y = Inf, 
        label = p, 
        hjust = 1.1, vjust = 1.1, size = 9, color='red') +
    theme_minimal() +
    theme(
  legend.position = 'right', 
  legend.text = element_text(size=9),
  legend.title = element_blank(),
  legend.spacing.y = unit(0.1, 'cm'),  # Reduced spacing
  legend.key.height = unit(0.1, 'cm'),  # Smaller keys
  axis.title = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size=12)
) +   
    guides(color = guide_legend(override.aes = list(size = 3))) +
    guides(fill = guide_legend(byrow = TRUE)) +


    ylim((min(df$InverseDifferenceMoment_00)), (max(df$InverseDifferenceMoment_00)+0.15))

}



### FIGURE 3B ###
pdf(paste0(figures_filepath, 'Figure_3/Fig3_B1_InverseDiffMom.pdf'), width = 4.5, height = 18)
ggarrange(plotlist = glist, nrow = 9)
dev.off()


###------------------------------------------------------------------####
### Figure 3 - low confluency
###------------------------------------------------------------------####

### Neurodevelopmental samples
pca_list_com <- do.call(rbind, pca_list)
### granularity 3
glist <- list()
for (c in 1:9){
  cgroups <- levels(factor(pca_list_com$con_grp))
  df <- pca_list_com[pca_list_com$con_grp == cgroups[c],]

  sample_means <- aggregate(Correlation_00 ~ Sample, data = df, mean)
  pc2_means <- aggregate(PC2 ~ Sample, data = df, mean)
  if (all(sample_means$Sample==pc2_means$Sample)==T){
    sample_means$PC2 <- pc2_means$PC2
  }

  sample_means <- sample_means[order(sample_means$PC2),]
 
  cols <- all_sample_cols[sample_means$Sample,]
  cols <- cols[sample_means$Sample,]
  cols$PC2 <- sample_means$PC2
  cols <- cols[order(cols$PC2),]

  cor <- cor.test(sample_means$Correlation_00, sample_means$PC2)
  p_value <- cor$p.value %>% round(.,6)
  r <- (cor$estimate) %>% round(., 3)

if (p_value < 0.001) {
  p <- '***'  # Assign '***' for p-values less than 0.001
} else if (p_value < 0.01) {
  p <- '**'   # Assign '**' for p-values less than 0.01
} else {
  p <- ''     # No asterisks for p-values 0.01 or higher
}


  # Create a ggplot
  glist[[c]] <- ggplot() +
    geom_line(data = sample_means, aes(x = reorder(Sample, PC2), y = Correlation_00, group=1), color = 'black', size = 0.9) +  # Line connecting means
    geom_jitter(data = df, aes(x = Sample, y = Correlation_00, color=Sample), size=0.3) +
    geom_line(data = sample_means, aes(x = Sample, y = Correlation_00, group=1), color = 'black', size = 0.9) +  # Line connecting means
    scale_color_manual(values=colors, breaks = sample_means$Sample) +
    annotate("text", x = Inf, y = Inf, 
        label = p, 
        hjust = 1.1, vjust = 1.1, size = 9, color='red') +
    theme_minimal() +
    theme(
  legend.position = 'right', 
  legend.text = element_text(size=9),
  legend.title = element_blank(),
  legend.spacing.y = unit(0.1, 'cm'),  # Reduced spacing
  legend.key.height = unit(0.1, 'cm'),  # Smaller keys
  axis.title = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size=12)
) +   
    guides(color = guide_legend(override.aes = list(size = 3))) +
    guides(fill = guide_legend(byrow = TRUE)) +


    ylim((min(df$Correlation_00)), (max(df$Correlation_00)+0.15))

}




### FIGURE 3B ###
pdf(paste0(figures_filepath, 'Figure_3/Fig3_B1_correlation.pdf'), width = 4.5, height = 18)
ggarrange(plotlist = glist, nrow = 9)
dev.off()









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






