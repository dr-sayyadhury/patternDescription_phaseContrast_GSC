### PACKAGES TO LOAD
###------------------------------------------------------------------#### 

library(datawizard)
library(effectsize)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(colorspace)
library(factoextra)
library(NbClust)
suppressPackageStartupMessages(library(dendextend))



### FILE/FOLDER PATHS
###------------------------------------------------------------------#### 
#path_to_repo <- "/data/"
path_to_repo <- '/Users/tpugh_m04/Documents/Research/Projects/PatternRecognitionInGSCs_UsingCV'

#figures_filepath <-'/results/'
figures_filepath <- paste0(path_to_repo,'/results/')
#script_path <- '/code'
script_path <- paste0(path_to_repo,'/scripts')
#tables_path <- '/results/tables'
tables_path <-  paste0(figures_filepath,'/tables/')
#dataset_path <- paste0(path_to_repo, 'cp_results_singleCell/cp_OutputAnalysis') ### to run and download again
dataset_path <- paste0(path_to_repo, '/datasets/cp_results_singleCell/cp_OutputAnalysis') ### to run and download again
#results_path <- '/results/S3_outPutDataFiles/'
results_path <- paste0(figures_filepath,'/S6_outPutDataFiles/')
dir.create(results_path)

source(paste0(script_path,"/sourceData_General_variables_and_functions.R"))



### IMPORT DATA FROM S1/S3
###------------------------------------------------------------------#### 

# S1
#phase <- read.csv('/results/S1_outPutDataFiles/final_dataset_afterQC.csv')
phase <- read.csv(paste0(path_to_repo, '/results/S1_outPutDataFiles/final_dataset_afterQC.csv'), row.names = 1)

# S3
GSC.gsva <- read.csv(paste0(path_to_repo,'/results/S3_outPutDataFiles/GSC.gsva.csv'), row.names=1)
#GSC.gsva <- read.csv('/results/S3_outPutDataFiles/GSC.gsva.csv', row.names=1)

#S4
#df_list_tp_subsetted <- readRDS('/results/S3_outPutDataFiles/df_list_tp_subsetted.rds')
df_list_tp_subsetted <- readRDS(paste0(figures_filepath,'S4_outPutDataFiles/df_list_tp_subsetted.rds'))
#pca_wholeImages_congrps <- readRDS('/results/S3_outPutDataFiles/pca_wholeImages_congrps.rds')
pca_wholeImages_congrps <- readRDS(paste0(figures_filepath,'S4_outPutDataFiles/pca_wholeImages_congrps.rds'))

### Original data
###------------------------------------------------------------------#### 
#rlogC <- rlogC <- read.table('/data/bulk_RNA/bulkGSC-rlog-transformed-counts.txt',sep='\t')
rlogC <- rlogC <- read.table(paste0(path_to_repo, '/datasets/bulk_RNA/bulkGSC-rlog-transformed-counts.txt'),sep='\t')
### rlog normalized bulk gene expression data >>> Standardize colnames for samples and transform matrix for CCA analysis
colnames(rlogC) <- gsub('A[0-9]*_Weiss_|A[0-9]*_Dirks_|A[0-9]*_|_line|_Line','', colnames(rlogC))
rlogC_t <- rlogC[levels(factor(phase$Sample))] %>% t(.) #%>% as.data.frame(.)
rlogC_t <- rlogC_t[, sapply(1:ncol(rlogC_t), function(x) sd(rlogC_t[, x]) > 0)] %>% as.data.frame(.)

scGrn_wholeImage <- read.csv(paste0(dataset_path, "/singleCell_Image.csv"))
scObj <- read.csv(paste0(dataset_path,"/singleCell_shrinkObjects.csv"))

source(paste0(script_path, '/sourceData_forScript06.R'))



#--------------------------------------------------------------------------------#
### --------------------------------- CODE START ------------------------------###
#--------------------------------------------------------------------------------#

counts_con_grps <- do.call('c', lapply(df_list_tp_subsetted, nrow))
counts_single_cell <- dplyr::count(scObjF_shortened, Sample, Well, time)
counts_single_cell$Well <- gsub('T', '', counts_single_cell$Well )
counts_single_cell$time <- gsub('hrs', '', counts_single_cell$time)

wholeImage_list_combined_con_grp <- data.frame()

for (l in 1:9){
  confluency_grp <- rep(paste0('C', l), counts_con_grps[l])
  df <- cbind(confluency_grp, df_list_tp_subsetted[[l]])
  wholeImage_list_combined_con_grp <- rbind(wholeImage_list_combined_con_grp, df)
}

### subset whole images T=12hrs and T=24hrs
wholeImage_list_combined_con_grp <- subset(wholeImage_list_combined_con_grp, TimePt== 'T_24' | TimePt== 'T_12')
wholeImage_list_combined_con_grp$TimePt <- gsub('T_', '', wholeImage_list_combined_con_grp$TimePt)

wholeImage_con_grp <- data.frame()
for (r in 1:nrow(counts_single_cell)){
  row <- subset(wholeImage_list_combined_con_grp, Sample==counts_single_cell$Sample[r] &
                  Well==counts_single_cell$Well[r] &
                  TimePt==counts_single_cell$time[r])
  wholeImage_con_grp <- rbind(wholeImage_con_grp, row)
  
}

### TABLES
###------------------------------------------------------------------#### 
# #SUPPLEMENTARY TABLE 4 #
###------------------------------------------------------------------#### 
dir.create(paste0(tables_path, '/S_Table4'))
write.csv(wholeImage_con_grp[1:5], paste0(tables_path, "/S_Table4/wholeImage_con_grp_selection.csv"))
#--------------------------------------------------------------------------------#
# Table end


### Normalization
scObjN <- as.data.frame(matrix(ncol=ncol(scObjF_shortened), nrow=0))
colnames(scObjN) <- colnames(scObjF_shortened)

colnames(scObjN) <- colnames(scObjF_shortened)
for (s in 1:length(levels(factor(scObjF_shortened$Sample)))){
  id <- levels(factor(scObjF_shortened$Sample))[s]
  sample <- subset(scObjF_shortened, Sample==id)
  sample <- sample %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x))
  meta <- sample[1:7]
  df <- sample[8:ncol(sample)]
  df <- apply(df, 2 ,FUN=function(x){x/max(x)}) %>% as.data.frame(.)
  df <- df %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x)) %>% as.data.frame(.)
  df <- apply(df,1, standardize) %>% t(.) %>% as.data.frame(.)
  df <- cbind(meta, df)
  scObjN <- rbind(scObjN, df)
}
scObjN<- scObjN %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x))
colnames(scObjN) <- colnames(scObjF_shortened)

meta <- scObjN[1:8]
#df <- scObjN[9:ncol(scObjN_reduced)] %>% t(.)
df <- scObjN[9:ncol(scObjN)] %>% t(.)
df[is.na(df)] <- 0

p <- PCAtools::pca(df, metadata=meta)
pcs <- which(cumsum(p$variance)>90)
pc1 <- pcs[1]
pc2 <- which(cumsum(p$variance)>75)[1]

var <- sort(which((p$variance[1:length(pcs)-1]-p$variance[2:(length(pcs))]) > 0.06), decreasing=T)[1]
minPC <- min(pc1, var)
PCAtools::screeplot(p, components=PCAtools::getComponents(p,1:(pc1+3)), vline=c(pc1, pc2), colBar='darkblue')

sc_pca <- c(scObjN, as.data.frame(p$rotated)) %>% as.data.frame(.)

u <- uwot::umap(t(df), n_neighbors = 50, min_dist=0.001, n_components = 2, metric='euclidean', pca=minPC)
colnames(u) <- c('umap1','umap2')
pca.df <- cbind(scObjN, p$rotated, u)

meta <- scObjN[1:8]
df <- scObjN[9:ncol(scObjN)] %>% t(.)
df[is.na(df)] <- 0



### FIGURES
###------------------------------------------------------------------####
# SUPPLEMENTARY FIG 7B #
###------------------------------------------------------------------#### 

df_s <- as.data.frame(scale(t(df)))
dist_mat <- dist(df_s, method='euclidean')
hclust_avg <- hclust(dist_mat, method='ward.D')
cut_avg <- cutree(hclust_avg, k = 3)
cut_avg <- as.data.frame(cut_avg)

pca_df2 <- cbind(pca.df, 
                 cut_avg)

if(dir.exists(paste0(figures_filepath, 'Sup_Fig_S7'))==F){
  dir.create(paste0(figures_filepath, 'Sup_Fig_S7'))
}

pdf(paste0(figures_filepath, 'Sup_Fig_S7/S7_A_umap.pdf'), width=9, height=8)
ggplot(pca_df2, aes(umap1, umap2, color=as.factor(cut_avg))) + 
  geom_point(size=1.2) + 
  theme_classic() + 
  scale_color_manual(values=cols_vivid) + 
  theme(axis.line = element_line(size=0.9, color='grey'),
        axis.title = element_text(size=30)) +
  guides(color=guide_legend(override.aes = list(size=6))) +
  guides(color=guide_legend(override.aes = list(size=4), title = 'Clusters'))
dev.off()

pdf(paste0(figures_filepath, 'Sup_Fig_S7/S7_A_barplot.pdf'), width=6, height=6.5)
ggplot(pca_df2, aes(as.factor(cut_avg),fill=Sample)) + geom_bar(position='fill') + scale_fill_manual(values=cols_vivid) + theme_minimal()
dev.off()
#--------------------------------------------------------------------------------#
#Fig end


### REMOVE CLUSTER 2
scObjF_short_cluster <- subset(cbind(scObjF_shortened,cut_avg), cut_avg!=2) %>% .[1:(ncol(.)-1)]

#Normalization
scObjN_c12 <- as.data.frame(matrix(ncol=ncol(scObjF_short_cluster), nrow=0))
colnames(scObjN_c12) <- colnames(scObjF_short_cluster)

colnames(scObjN_c12) <- colnames(scObjF_short_cluster)
for (s in 1:length(levels(factor(scObjF_short_cluster$Sample)))){
  id <- levels(factor(scObjF_short_cluster$Sample))[s]
  sample <- subset(scObjF_short_cluster, Sample==id)
  sample <- sample %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x))
  meta <- sample[1:7]
  df <- sample[8:ncol(sample)]
  df <- apply(df, 2 ,FUN=function(x){x/max(x)}) %>% as.data.frame(.)
  df <- df %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x)) %>% as.data.frame(.)
  df <- apply(df,1, standardize) %>% t(.) %>% as.data.frame(.)
  df <- cbind(meta, df)
  scObjN_c12 <- rbind(scObjN_c12, df)
}
scObjN_c12 <- scObjN_c12 %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x))
colnames(scObjN_c12) <- colnames(scObjF_short_cluster)

meta <- scObjN_c12[1:8]
df <- scObjN_c12[9:ncol(scObjN_c12)] %>% t(.)
df[is.na(df)] <- 0
  
p <- PCAtools::pca(df, metadata=meta)
pcs <- which(cumsum(p$variance)>90)
pc1 <- pcs[1]
pc2 <- which(cumsum(p$variance)>75)[1]
var <- sort(which((p$variance[1:length(pcs)-1]-p$variance[2:(length(pcs))]) > 0.06), decreasing=T)[1]
minPC <- min(pc1, var)
PCAtools::screeplot(p, components=PCAtools::getComponents(p,1:(pc1+3)), vline=c(pc1, pc2), colBar='darkblue')



sc_pca <- c(scObjN_c12, as.data.frame(p$rotated)) %>% as.data.frame(.)
sample_levels <- levels(factor(sc_pca$Sample))
col_by_proportion <- c('#0092C8','darkgrey','#DE1D82','#DE1D82','#DE1D82','#DE1D82','#0092C8','#0092C8')

### FIGURES
###------------------------------------------------------------------####
# FIGURE 4B #
###------------------------------------------------------------------#### 

pdf(paste0(figures_filepath, 'Figure_4/Fig4_B.pdf'), width=9, height=8)
ggplot(sc_pca, aes(PC1,PC2, color=Sample)) +
  geom_point(size=1.8) + theme_pander() + 
  scale_color_manual(values=cols_vivid) + 
  theme(axis.line = element_line(size=0.9, color='grey'),
        axis.title = element_text(size=30),
        legend.text = element_text(size=9)) 
dev.off()

pdf(paste0(figures_filepath, 'Figure_4/Fig4_B2.pdf'), width=9, height=8)
ggplot(sc_pca, aes(PC1,PC2, color=Sample)) +
  geom_point(size=1.2) + theme_pander() + 
  scale_color_manual(values=col_by_proportion) + 
  theme(axis.line = element_line(size=0.9, color='grey'),
        axis.title = element_text(size=30),
        legend.text = element_text(size=9)) 
dev.off()

#--------------------------------------------------------------------------------#
#Fig end



PCAtools::plotloadings(p,components = PCAtools::getComponents(p, seq_len(pc2)), labSize=1.8, rangeRetain=0.01, shapeSizeRange = c(6,6))
u <- uwot::umap(t(df), n_neighbors = 50, min_dist=0.001, n_components = 2, metric='euclidean', pca=minPC)
colnames(u) <- c('umap1','umap2')
umap_pca.df_sc <- cbind(sc_pca, u)


wholeImage_con_grp <- wholeImage_con_grp[1:5]
pca.df_confGrp <- data.frame()

for (r in 1:nrow(wholeImage_con_grp)){
  row <- subset(umap_pca.df_sc, umap_pca.df_sc$Sample==wholeImage_con_grp$Sample[r] & 
                  umap_pca.df_sc$Well==paste0(wholeImage_con_grp$Well[r],'T') & 
                  umap_pca.df_sc$time==paste0(wholeImage_con_grp$TimePt[r], 'hrs'))
  dup_con_grp <- rep(wholeImage_con_grp$confluency_grp[r], nrow(row))
  row <- cbind(dup_con_grp, row)
  pca.df_confGrp <- rbind(pca.df_confGrp, row)
}



### FIGURES
###------------------------------------------------------------------####
# FIGURE 4C #
###------------------------------------------------------------------#### 

df_s <- as.data.frame(scale(pca.df_confGrp[c(350,351)]))
dist_mat <- dist(df_s)
hclust_avg <- hclust(dist_mat)
cut_avg <- cutree(hclust_avg, k = 2)

cut_avg <- as.data.frame(cut_avg)

pca.df_confGrp <- cbind(pca.df_confGrp[c(1,3,350,351)], 
                        #pca_df2, 
                        cut_avg)

pdf(paste0(figures_filepath, 'Figure_4/Fig4_C_umap.pdf'), width=6, height=5.5)
ggplot(pca.df_confGrp, aes(umap1, umap2, color=as.factor(cut_avg))) + 
  geom_point(size=0.3) + theme_pander() + 
  scale_color_manual(values=cols_vivid) + 
  guides(color=guide_legend(override.aes = list(size=4)))
dev.off()

pca.df_confGrp$Sample <- umap_pca.df_sc$Sample

pdf(paste0(figures_filepath, 'Figure_4/Fig4_C_barplot.pdf'), width=6, height=6.5)
ggplot(pca.df_confGrp, aes(as.factor(cut_avg),fill=Sample)) + geom_bar(position='fill') + scale_fill_manual(values=cols_vivid) + theme_minimal()
dev.off()
#--------------------------------------------------------------------------------#
#Fig end


samples_rem <- levels(factor(pca.df_confGrp$Sample))

g <- list()
g1 <- list()
g2 <- list()

for (c in 1:length(levels(factor(pca.df_confGrp$dup_con_grp)))){
  cgrp <- subset(pca.df_confGrp, dup_con_grp==levels(factor(pca.df_confGrp$dup_con_grp))[c])
  
  g[[c]] <- ggplot(pca.df_confGrp, aes(umap1, umap2)) + 
    geom_point(size=0.6, colour='lightgrey') +
    geom_point(data=cgrp, colour=col_bold[c], size=0.7) +
    theme_classic() + 
    theme(axis.line = element_line(size=0.6, color='grey'),
          axis.title = element_text(size=12))
  
  g2[[c]] <- ggplot(pca.df_confGrp, aes(umap1, umap2)) + 
    geom_point(size=0.6, colour='lightgrey') +
    geom_point(data=cgrp, aes(color=Sample), size=0.7) +
    theme_classic() + 
    scale_color_manual(values=cols_vivid[samples_rem %in% levels(factor(cgrp$Sample))]) + 
    theme(axis.line = element_line(size=0.6, color='grey'),
          axis.title = element_text(size=12),
          legend.title = element_blank())
  
}  



### FIGURES
###------------------------------------------------------------------####
# SUPPLEMENTARY FIG 7B #
###------------------------------------------------------------------#### 


pdf(paste0(figures_filepath, 'Sup_Fig_S7/S7_B_top.pdf'), width=15, height=3)
ggarrange(plotlist=g, ncol=5)
dev.off()

pdf(paste0(figures_filepath, 'Sup_Fig_S7/S7_B_bot.pdf'), width=15, height=3)
ggarrange(plotlist=g2, ncol=5)
dev.off()

#--------------------------------------------------------------------------------#
#Fig end


### FIGURES
###------------------------------------------------------------------####
# SUPPLEMENTARY FIG 4E #
###------------------------------------------------------------------#### 

#pdf(paste0(figures_filepath, 'Sup_Fig_S4/S4_E.pdf'), width=18, height=12)
#ggplot(pca.df_confGrp, aes(Sample, fill=dup_con_grp)) + 
#  geom_bar() +
#  theme_classic() + 
#  scale_fill_manual(values=cols_vivid) + 
#  theme(axis.line = element_line(size=0.6, color='grey'),
#        axis.title = element_blank(),
#        axis.text = element_text(size=36))
##guides(color=guide_legend(override.aes = list(size=6))) +
##guides(color=guide_legend(override.aes = list(size=18), title = 'con grp'))
#dev.off()
#--------------------------------------------------------------------------------#
#Fig end




### SUBSET ONLY CLUSTER 1 AND RE-ANALYZE AGAIN
scObjF_short_cluster$cut_avg <- pca.df_confGrp$cut_avg

scObjF_short_clustera <- subset(scObjF_short_cluster, cut_avg==1)
scObjF_short_clusterb <- subset(scObjF_short_cluster, cut_avg==2)

scObjF_short_clustera <- scObjF_short_clustera[1:178]
scObjF_short_clusterb <- scObjF_short_clusterb[1:178]

### Check for single images only
dplyr::count(scObjF_short_clustera, Sample)
dplyr::count(scObjF_short_clusterb, Sample)

# Remove G729=1 and G851=3 and G799=4
# Remove G549=1, G583=1, G837=4, G885=1
scObjF_short_clustera <- subset(scObjF_short_clustera, Sample !='G729' &Sample !='G799' & Sample !='G851')
scObjF_short_clusterb <- subset(scObjF_short_clusterb, Sample !='G549' & Sample!='G583' & Sample !='G837' & Sample!='G885')


#---PROCESS scObjM_c12a FIRST

scObjN_c12a <- as.data.frame(matrix(ncol=ncol(scObjF_short_clustera), nrow=0))
colnames(scObjN_c12a) <- colnames(scObjF_short_clustera)

colnames(scObjN_c12a) <- colnames(scObjF_short_clustera)
for (s in 1:length(levels(factor(scObjF_short_clustera$Sample)))){
  id <- levels(factor(scObjF_short_clustera$Sample))[s]
  sample <- subset(scObjF_short_clustera, Sample==id)
  sample <- sample %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x))
  meta <- sample[1:7]
  df <- sample[8:ncol(sample)]
  df <- apply(df, 2, FUN=function(x)if(mean(x)==0)
  {x=x}
  else{x/max(x)}) %>% as.data.frame(.) ### note keep this  ##df[is.na(df)] <- 0
  df <- df %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x)) %>% as.data.frame(.)
  df <- apply(df,1, standardize) %>% t(.) %>% as.data.frame(.)
  ##df[is.na(df)] <- 0
  ##df <- df %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x)) %>% as.data.frame(.)
  df <- cbind(meta, df)
  scObjN_c12a <- rbind(scObjN_c12a, df)
}
scObjN_c12a <- scObjN_c12a %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x))
colnames(scObjN_c12a) <- colnames(scObjF_short_clustera)

meta <- scObjN_c12a[1:7]
#df <- scObjN[9:ncol(scObjN_reduced)] %>% t(.)
df <- scObjN_c12a[8:ncol(scObjN_c12a)] %>% t(.)
df[is.na(df)] <- 0

p <- PCAtools::pca(df, metadata=meta)
pcs <- which(cumsum(p$variance)>90)
pc1 <- pcs[1]
pc2 <- which(cumsum(p$variance)>75)[1]

var <- sort(which((p$variance[1:length(pcs)-1]-p$variance[2:(length(pcs))]) > 0.06), decreasing=T)[1]
minPC <- min(pc1, var)
PCAtools::screeplot(p, components=PCAtools::getComponents(p,1:(pc1+3)), vline=c(pc1, pc2), colBar='darkblue')

sc_pca <- c(scObjF_short_clustera, as.data.frame(p$rotated)) %>% as.data.frame(.)

sample_levels <- levels(factor(sc_pca$Sample))
col_by_proportion <- c('#0092C8','black','#DE1D82','#DE1D82','#0092C8')



### FIGURES
###------------------------------------------------------------------####
# FIGURE 5A #
###------------------------------------------------------------------#### 

if(dir.exists(paste0(figures_filepath, 'Figure_5'))==F){
  dir.create(paste0(figures_filepath, 'Figure_5'))
}

pdf(paste0(figures_filepath, 'Figure_5/Fig5_A.pdf'), width=9, height=8)
ggplot(sc_pca, aes(PC2,PC1, color=Sample)) +
  geom_point(size=3, alpha=1) + theme_pander() + 
  scale_color_manual(values=c(cols_vivid)) + 
  theme(axis.line = element_line(size=0.9, color='grey'),
        axis.title = element_text(size=30)) +
  guides(color=guide_legend(override.aes = list(size=6)))
dev.off()

pdf(paste0(figures_filepath, 'Figure_5/Fig5_A2.pdf'), width=9, height=8)
ggplot(sc_pca, aes(PC1,PC2, color=Sample)) +
  geom_point(size=2.4) + theme_pander() + 
  scale_color_manual(values=col_by_proportion) + 
  theme(axis.line = element_line(size=0.9, color='grey'),
        axis.title = element_text(size=30),
        legend.text = element_text(size=9)) 
dev.off()


 #--------------------------------------------------------------------------------#
#Fig end


#---PROCESS scObjM_c12B NEXT
scObjN_c12b <- as.data.frame(matrix(ncol=ncol(scObjF_short_clusterb), nrow=0))
colnames(scObjN_c12b) <- colnames(scObjF_short_clusterb)

colnames(scObjN_c12b) <- colnames(scObjF_short_clusterb)
for (s in 1:length(levels(factor(scObjF_short_clusterb$Sample)))){
  id <- levels(factor(scObjF_short_clusterb$Sample))[s]
  sample <- subset(scObjF_short_clusterb, Sample==id)
  sample <- sample %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x))
  meta <- sample[1:7]
  df <- sample[8:ncol(sample)]
  df <- apply(df, 2, FUN=function(x)if(mean(x)==0)
  {x=x}
  else{x/max(x)}) %>% as.data.frame(.) ### note keep this  ##df[is.na(df)] <- 0
  df <- df %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x)) %>% as.data.frame(.)
  df <- apply(df,1, standardize) %>% t(.) %>% as.data.frame(.)
  ##df[is.na(df)] <- 0
  ##df <- df %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x)) %>% as.data.frame(.)
  df <- cbind(meta, df)
  scObjN_c12b <- rbind(scObjN_c12b, df)
}
scObjN_c12b <- scObjN_c12b %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x))
colnames(scObjN_c12b) <- colnames(scObjF_short_clusterb)

meta <- scObjN_c12b[1:7]
#df <- scObjN[9:ncol(scObjN_reduced)] %>% t(.)
df <- scObjN_c12b[8:ncol(scObjN_c12b)] %>% t(.)
df[is.na(df)] <- 0

p <- PCAtools::pca(df, metadata=meta)
pcs <- which(cumsum(p$variance)>90)
pc1 <- pcs[1]
pc2 <- which(cumsum(p$variance)>75)[1]

var <- sort(which((p$variance[1:length(pcs)-1]-p$variance[2:(length(pcs))]) > 0.06), decreasing=T)[1]
minPC <- min(pc1, var)
PCAtools::screeplot(p, components=PCAtools::getComponents(p,1:(pc1+3)), vline=c(pc1, pc2), colBar='darkblue')

sc_pcab <- c(scObjF_short_clusterb, as.data.frame(p$rotated)) %>% as.data.frame(.)

sample_levels <- levels(factor(sc_pcab$Sample))
col_by_proportion <- c('#DE1D82','#0092C8','#0092C8')



### FIGURES
###------------------------------------------------------------------####
# FIGURE 5A #
###------------------------------------------------------------------#### 

if(dir.exists(paste0(figures_filepath, 'Figure_5'))==F){
  dir.create(paste0(figures_filepath, 'Figure_5'))
}

pdf(paste0(figures_filepath, 'Figure_5/Fig5_A3.pdf'), width=9, height=8)
ggplot(sc_pcab, aes(PC1,PC2, color=Sample)) +
  geom_point(size=2.4) + theme_pander() + 
  scale_color_manual(values=col_by_proportion) + 
  theme(axis.line = element_line(size=0.9, color='grey'),
        axis.title = element_text(size=30),
        legend.text = element_text(size=9)) 
dev.off()



pdf(paste0(figures_filepath, 'Figure_5/Fig5_A4.pdf'), width=9, height=8)
ggplot(sc_pcab, aes(PC3,PC2, color=Sample)) +
  geom_point(size=2.4) + theme_pander() + 
  scale_color_manual(values=col_by_proportion) + 
  theme(axis.line = element_line(size=0.9, color='grey'),
        axis.title = element_text(size=30),
        legend.text = element_text(size=9)) 
dev.off()


pdf(paste0(figures_filepath, 'Figure_5/Fig5_A5.pdf'), width=9, height=8)
ggplot(sc_pcab, aes(PC3,PC4, color=Sample)) +
  geom_point(size=2.4) + theme_pander() + 
  scale_color_manual(values=col_by_proportion) + 
  theme(axis.line = element_line(size=0.9, color='grey'),
        axis.title = element_text(size=30),
        legend.text = element_text(size=9)) 
dev.off()


pdf(paste0(figures_filepath, 'Figure_5/Fig5_A6.pdf'), width=9, height=8)
ggplot(sc_pcab, aes(PC5,PC4, color=Sample)) +
  geom_point(size=2.4) + theme_pander() + 
  scale_color_manual(values=col_by_proportion) + 
  theme(axis.line = element_line(size=0.9, color='grey'),
        axis.title = element_text(size=30),
        legend.text = element_text(size=9)) 
dev.off()


#--------------------------------------------------------------------------------#
#Fig end


# SKIP
# SKIP
# SKIP
########
#######
#######


### FIGURES
###------------------------------------------------------------------####
# FIGURE 4 #
###------------------------------------------------------------------#### 

c3 <- pca_wholeImages_congrps[[3]]
wc3 <- subset(wholeImage_con_grp, confluency_grp=='C3')
wc3 <- subset(c3, rownames(c3) %in% rownames(wc3))

pdf(paste0(figures_filepath, 'Figure_4/Fig4_D_C3.pdf'), width=9, height=6)
ggplot(c3, aes(reorder(Sample, PC2), x=PC2,  fill=Sample)) + 
  geom_boxplot(outlier.color = 'grey') + 
  scale_fill_manual(values=cols_vivid) + 
  geom_point(data=wc3, size=6, shape=18, color='red') +
  theme_tufte() +
  theme(axis.text = element_text(size=18),
        axis.title = element_blank(),
        legend.position = 'none') 
dev.off()


c4 <- pca_wholeImages_congrps[[4]]
wc4 <- subset(wholeImage_con_grp, confluency_grp=='C4')
wc4 <- subset(c4, rownames(c4) %in% rownames(wc4))

pdf(paste0(figures_filepath, 'Figure_4/Fig4_C4.pdf'), width=9, height=6)
ggplot(c4, aes(reorder(Sample, PC2), x=PC2,  fill=Sample)) + 
  geom_boxplot(outlier.color = 'grey') + 
  scale_fill_manual(values=cols_vivid) + 
  geom_point(data=wc4, size=6, shape=18, color='red') +
  theme_tufte() +
  theme(axis.text = element_text(size=18),
        axis.title = element_blank(),
        legend.position = 'none') 
dev.off()


c6 <- pca_wholeImages_congrps[[6]]
wc6 <- subset(wholeImage_con_grp, confluency_grp=='C6')
wc6 <- subset(c6, rownames(c6) %in% rownames(wc6))

pdf(paste0(figures_filepath, 'Figure_4/Fig4_C6.pdf'), width=9, height=6)
ggplot(c6, aes(reorder(Sample, PC2), x=PC2,  fill=Sample)) + 
  geom_boxplot(outlier.color = 'grey') + 
  scale_fill_manual(values=cols_vivid) + 
  geom_point(data=wc6, size=6, shape=18, color='red') +
  theme_tufte() +
  theme(axis.text = element_text(size=18),
        axis.title = element_blank(),
        legend.position = 'none') 
dev.off()
#--------------------------------------------------------------------------------#
#Fig end


### FIGURES
###------------------------------------------------------------------####
# FIGURE 5B #
###------------------------------------------------------------------#### 
c1 <- pca_wholeImages_congrps[[1]]
wc1 <- subset(wholeImage_con_grp, confluency_grp=='C1')
wc1 <- subset(c1, rownames(c1) %in% rownames(wc1))


pdf(paste0(figures_filepath, 'Figure_5/Fig5_B_C1.pdf'), width=9, height=6)
ggplot(c1, aes(reorder(Sample, PC2), x=PC2,  fill=Sample)) + 
  geom_boxplot(outlier.color = 'grey') + 
  scale_fill_manual(values=cols_vivid) + 
  geom_point(data=wc1, size=6, shape=18, color='red') +
  theme_tufte() +
  theme(axis.text = element_text(size=18),
        axis.title = element_blank(),
        legend.position = 'none') 
dev.off()
#--------------------------------------------------------------------------------#
#Fig end



c2 <- pca_wholeImages_congrps[[2]]
wc2 <- subset(wholeImage_con_grp, confluency_grp=='C2')
wc2 <- subset(c2, rownames(c2) %in% rownames(wc2))

pdf(paste0(figures_filepath, 'Figure_5/Fig5_B_C2.pdf'), width=9, height=6)
ggplot(c2, aes(reorder(Sample, PC2), x=PC2,  fill=Sample)) + 
  geom_boxplot(outlier.color = 'grey') + 
  scale_fill_manual(values=cols_vivid) + 
  geom_point(data=wc2, size=6, shape=18, color='red') +
  theme_tufte() +
  theme(axis.text = element_text(size=18),
        axis.title = element_blank(),
        legend.position = 'none') 
dev.off()
#--------------------------------------------------------------------------------#
#Fig end



#remove(list=ls())

### reduce samples to smaller set of features
wc1_wc2 <- rbind(wc1, wc2)
meta <- wc1_wc2[1:6]
wc1_wc2 <- wc1_wc2[c(1:6,36:64)]
wc1_wc2 <- subset(wc1_wc2, Sample != 'G729')

## analyse all images together withoyt grouping by confluency groups
features=35
new_df <- data.frame()
n=length(levels(factor(wc1_wc2$Sample)))
for (l in 1:n){
  s <- levels(factor(wc1_wc2$Sample))[l]
  df <- subset(wc1_wc2, Sample==s)
  meta <- df[1:6]
  df <- df[7:features]
  df <- apply(df, 2, FUN=function(x)if(mean(x)==0)
  {x=x}
  else{x/max(x)}) %>% as.data.frame(.) ### note keep this
  df <- apply(df,1, standardize)
  df <- cbind(meta,t(df)) 
  df <- df %>% .[order(.$TimePt),]
  
  meta <- df[1:6]
  df <- df[7:features]  %>% t(.) %>% as.data.frame(.)
  rownames(meta) <- colnames(df)
  #con_grp <- rep(paste0('C',i), nrow(meta))
  #meta <- cbind(con_grp, meta)
  df <- cbind(meta, t(df))
  
  colnames(df) <- colnames(wc1_wc2)
  new_df <- rbind(new_df, df)
}


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



PCAtools::screeplot(p, colBar='darkblue', axisLabSize = 9, components = 1:(pc1+1), title = '', vline=c(var,pc1))
pca_df <- cbind(meta,p$rotated)
colnames(pca_df) <- c('con_grp', 'Media','Sample','TimePt','Well','Area', paste0('PC', 1:ncol(p$rotated)))

pca_df <- cbind(pca_df,data)


#pdf(paste0(figures_filepath,'Figure_1/Fig1_I.pdf'), width=4.5, height=4.5)
ggplot(pca_df, aes(PC3, PC2, color=con_grp)) +
  geom_jitter(size=3) +
  theme_minimal() +
  scale_color_manual(values=c(cols_vivid, col_bold)) +
  theme(legend.position = 'right', 
        axis.title = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle=45, size=9)) +
  #geom_text_repel(aes(label=Sample), size=2.5) +
  xlab("PC2") + ylab("PC1")
#dev.off()

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
    
    colnames(df) <- c('con_grp', colnames(d))
    new_df <- rbind(new_df, df)
    
  }
  new_df_list[[i]] <- new_df
}

library(cowplot)

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
ggarrange(plotlist=g1, nrow=9, ncol=1)

#remove(list=ls())

### reduce samples to smaller set of features
meta <- wc2[1:6]
wc2_ <- wc2[c(1:6,36:64)]
wc2_ <- subset(wc2_, Sample != 'G729')

## analyse all images together withoyt grouping by confluency groups
features=35
new_df <- data.frame()
n=length(levels(factor(wc2_$Sample)))
for (l in 1:n){
  s <- levels(factor(wc2_$Sample))[l]
  df <- subset(wc2_, Sample==s)
  meta <- df[1:6]
  df <- df[7:features]
  df <- apply(df, 2, FUN=function(x)if(mean(x)==0)
  {x=x}
  else{x/max(x)}) %>% as.data.frame(.) ### note keep this
  df <- apply(df,1, standardize)
  df <- cbind(meta,t(df)) 
  df <- df %>% .[order(.$TimePt),]
  
  meta <- df[1:6]
  df <- df[7:features]  %>% t(.) %>% as.data.frame(.)
  rownames(meta) <- colnames(df)
  #con_grp <- rep(paste0('C',i), nrow(meta))
  #meta <- cbind(con_grp, meta)
  df <- cbind(meta, t(df))
  
  colnames(df) <- colnames(wc1_wc2)
  new_df <- rbind(new_df, df)
}


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



PCAtools::screeplot(p, colBar='darkblue', axisLabSize = 9, components = 1:(pc1+1), title = '', vline=c(var,pc1))
pca_df <- cbind(meta,p$rotated)
colnames(pca_df) <- c('con_grp', 'Media','Sample','TimePt','Well','Area', paste0('PC', 1:ncol(p$rotated)))

pca_df <- cbind(pca_df,data)


#pdf(paste0(figures_filepath,'Figure_1/Fig1_I.pdf'), width=4.5, height=4.5)
ggplot(pca_df, aes(PC2, PC3, color=Sample)) +
  geom_jitter(size=3) +
  theme_minimal() +
  scale_color_manual(values=c(cols_vivid, col_bold)) +
  theme(legend.position = 'right', 
        axis.title = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle=45, size=9)) +
  #geom_text_repel(aes(label=Sample), size=2.5) +
  xlab("PC2") + ylab("PC1")
#dev.off()



########
#######
#######
# SKIP
# SKIP
# SKIP


#wc1_wc2 <- rbind(wc1, wc2)
#wc1_wc2 <- subset(wc1_wc2, Sample != 'G729')
wc1_ <- subset(wc1, Sample != 'G729')
sample_levels <- levels(factor(wc1$Sample))

new_df_list <- list()

for (dframe in 1:length(df_list_tp_subsetted)){
  data <- df_list_tp_subsetted[[dframe]]
  data <- subset(data, Sample %in% levels(factor(wc1_$Sample)))

  ### reduce samples to smaller set of features
  meta <- data[1:5]
  ## analyse all images together withoyt grouping by confluency groups
  features=34
  new_df <- data.frame()
  n=length(levels(factor(data$Sample)))
  for (l in 1:n){
    s <- levels(factor(data$Sample))[l]
    df <- subset(data, Sample==s)
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
    #con_grp <- rep(paste0('C',i), nrow(meta))
    #meta <- cbind(con_grp, meta)
    df <- cbind(meta, t(df))
  
    colnames(df) <- colnames(data)
    new_df <- rbind(new_df, df)
  }
  new_df_list[[dframe]] <- new_df
  
}


pca_new_list <- list()
g <- list()

for (dframe in 1:length(df_list_tp_subsetted)){
  new_df <- new_df_list[[dframe]]
  new_df[is.na(new_df)] <-0
  meta <- new_df[1:5]
  data <- new_df[6:features]
  rownames(data) <- rownames(meta)

  p <- PCAtools::pca(t(data), metadata = meta, center=T, scale=F)

  pcs <- which(cumsum(p$variance)>90)
  pc1 <- pcs[1]
  pc2 <- which(cumsum(p$variance)>75)[1]

  var <- sort(which((p$variance[1:length(pcs)-1]-p$variance[2:(length(pcs))]) > 0.1), decreasing=T)[1]
  minPC <- min(pc1, var)

  PCAtools::screeplot(p, colBar='darkblue', axisLabSize = 9, components = 1:(pc1+1), title = '', vline=c(var,pc1))
  pca_df <- cbind(meta,p$rotated)
  colnames(pca_df) <- c('Media','Sample','TimePt','Well','Area', paste0('PC', 1:ncol(p$rotated)))

  pca_df <- cbind(pca_df,data)
  pca_new_list <- pca_df

  #pdf(paste0(figures_filepath,'Figure_1/Fig1_I.pdf'), width=4.5, height=4.5)
  g[[dframe]] <- ggplot(pca_df, aes(PC2, PC1, color=Sample)) +
    geom_jitter(size=1) +
    theme_minimal() +
    scale_color_manual(values=c(cols_vivid, col_bold)) +
    theme(legend.position = 'right', 
        axis.title = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle=45, size=9)) +
    #geom_text_repel(aes(label=Sample), size=2.5) +
    xlab("PC2") + ylab("PC1")
  #dev.off()
}

ggarrange(plotlist=g, nrow=9)



