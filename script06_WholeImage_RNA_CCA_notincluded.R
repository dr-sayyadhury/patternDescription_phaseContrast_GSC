### PACKAGES TO LOAD
###------------------------------------------------------------------#### 
library(dplyr)
library(ggplot2)
library(ggthemes)
library(tidyr)
library(forcats)
library(colorspace)
library(datawizard)
library(PMA)



### FILE/FOLDER PATHS
###------------------------------------------------------------------#### 
path_to_repo <- '/Users/tpugh_m04/Documents/Research/PostdoctoralResearch_2020/Projects/PatternRecognitionInGSCs_UsingCV'
figures_filepath <- paste0(path_to_repo,'/results/figures/')
script_path <- paste0(path_to_repo,'/scripts/')
results_path <- paste0(path_to_repo,'/results/S4_outPutDataFiles/')
dir.create(results_path)
tables_path <-  paste0(path_to_repo,'/results/tables/')




### source data 
###------------------------------------------------------------------#### 
#source(paste0(script_path,"/script04b_.R"))
source(paste0(script_path,"/sourceData_General_variables_and_functions.R"))
GSC.gsva <- read.csv(paste0(path_to_repo,'/results/GSC.gsva.csv'), row.names=1)
phase <- read.csv(paste0(path_to_repo, '/results/processed_pixel_feature_values_for_QC_and_GLCM_Grn/final_dataset_afterQC.csv'), row.names = 1)
ylist_cca <- readRDS(paste0(path_to_repo,'/results/S3_outPutDataFiles/ylist_cca.rds'))

rlogC <- rlogC <- read.table(paste0(path_to_repo, '/datasets/bulk_RNA/bulkGSC-rlog-transformed-counts.txt'),sep='\t')
colnames(rlogC) <- gsub('A[0-9]*_Weiss_|A[0-9]*_Dirks_|A[0-9]*_|_line|_Line','', colnames(rlogC))
rlogC_t <- rlogC[levels(factor(phase$Sample))] %>% t(.) #%>% as.data.frame(.)



#--------------------------------------------------------------------------------#
### --------------------------------- CODE START ------------------------------###
#--------------------------------------------------------------------------------#
### CCA ANALYSIS 
### Comparing average pixel-derived scores with bulk gene-expression scores
bcca_list <- list() 
bulk.list <- list()
bcca.cov <- numeric(0)
img.list <- list()
for (i in 1:9){
  
  cAreaMax <- data.frame()
  common <-intersect(rownames(rlogC_t), rownames(ylist_cca[[i]]))
  x <- subset(rlogC_t, subset=rownames(rlogC_t) %in% common)
  x[is.na(x)] <- 0
  y <- subset(ylist_cca[[i]], subset=rownames(ylist_cca[[i]]) %in% common)
  y[is.na(y)] <- 0
  bCCA <- CCA(x, y, 
              typex="standard",
              typez="standard", 
              K=1, 
              penaltyz=0.6, 
              penaltyx = 1, 
              niter = 20)
  bulk.list[[i]] <- bCCA$u
  rownames(bulk.list[[i]]) <- colnames(rlogC_t)
  colnames(bulk.list[[i]]) <- paste0('gbCCA', 1)
  img.list[[i]] <- bCCA$v
  bcca.cov <- c(bcca.cov, bCCA$cors[1])
  rownames(img.list[[i]]) <- colnames(ylist_cca[[i]])
  colnames(img.list[[i]]) <- paste0('ibCCA', 1)
  bulk <- as.matrix(x)%*%bCCA$u
  img <- as.matrix(y)%*%bCCA$v
  cAreaMax <- cbind(as.data.frame(common),bulk,img)
  cAreaMax <- cbind(cAreaMax, rep(paste0('C', i), nrow(cAreaMax)))
  
  colnames(cAreaMax) <- c(
    "Sample", 
    paste0('gene_expression_', 1),# canonical variate for gene expression
    paste0('TxtFeat', 1),
    'con_grp') #canonical variates for txt feat
  bcca_list[[i]] <- cAreaMax 
}

bcca.cov <- cbind(paste0('C', 1:9), rep('bulk', length(ylist_cca)), as.data.frame(bcca.cov[1:9]))
colnames(bcca.cov) <- c('Confluency', 'Type', 'Cor')
bcca_rlogC <- bcca.cov

bcca_list2 <- list()

for (l in 1:9){
  bc <- bcca_list[[l]]
  g <- GSC.gsva %>% .[bc$Sample] %>% as.data.frame(.)
  bc <- cbind(bc, t(g))
  bc$gene_expression_1 <- standardize(bc$gene_expression_1)
  bc$TxtFeat1 <- standardize(bc$TxtFeat1)
  bcca_list2[[l]] <- bc
  
}

bcca_list_combined <- do.call(rbind, bcca_list2)

if(dir.exists(paste0(figures_filepath,'Sup_Fig_S5'))==F){
  dir.create(paste0(figures_filepath,'Sup_Fig_S5'))}

pdf(file =paste0(figures_filepath, 'Sup_Fig_S5/S5_A.pdf'), width=6, height=6)
ggplot(bcca_list_combined, aes(con_grp, gene_expression_1), gene_expression_1) + 
  geom_boxplot() + 
  geom_count(aes(color=Neftel_NPC1)) + 
  scale_color_continuous_divergingx('Tropic')


ggplot(bcca_list_combined, aes(con_grp, TxtFeat1), TxtFeat1) + 
  geom_boxplot() + 
  geom_count(aes(color=Neftel_NPC1)) + 
  scale_color_continuous_divergingx('Tropic')
dev.off()



### FIGURES
###------------------------------------------------------------------#### 
# SUPPLEMENTARY FIGURE 5A
###------------------------------------------------------------------#### 
if(dir.exists(paste0(figures_filepath,'Sup_Fig_S5'))==F){
  dir.create(paste0(figures_filepath,'Sup_Fig_S5'))}

pdf(file =paste0(figures_filepath, 'Sup_Fig_S5/S5_B_top.pdf'), width=6, height=6)
ggplot(bcca_list_combined, aes(TxtFeat1, gene_expression_1, grouping=con_grp)) + 
  geom_line() + 
  geom_count(aes(color=TxtFeat1+gene_expression_1)) + 
  theme_tufte() +
  theme(legend.position = 'none') +
  scale_color_continuous_divergingx('Tropic')
dev.off()

pdf(file =paste0(figures_filepath, 'Sup_Fig_S5/S5_B_bot.pdf'), width=6, height=6)
ggplot(bcca_list_combined, aes(TxtFeat1, gene_expression_1, grouping=con_grp, color=con_grp)) + 
  geom_line() + 
  #geom_count(aes(color=Neftel_NPC1)) + 
  theme_tufte() +
  theme(legend.position = 'none') +
  scale_color_manual(values=cols_vivid)
dev.off()


pdf(file =paste0(figures_filepath, 'Sup_Fig_S5/S5_C1.pdf'), width=6, height=6)
ggplot(bcca_list_combined, aes(TxtFeat1, gene_expression_1, grouping=con_grp, color=Neftel_NPC1)) + 
  geom_line() + 
  #geom_count(aes(color=Neftel_NPC1)) + 
  theme_tufte() +
  theme(legend.position = 'none') +
  scale_color_manual(values=cols_vivid)
dev.off()



pdf(file =paste0(figures_filepath, 'Sup_Fig_S5/S5_C1.pdf'), width=6, height=6)
ggplot(bcca_list_combined, aes(TxtFeat1, gene_expression_1, grouping=con_grp, color=bcca_list_combined$)) + 
  geom_line() + 
  #geom_count(aes(color=Neftel_NPC1)) + 
  theme_tufte() +
  theme(legend.position = 'none') +
  scale_color_manual(values=cols_vivid)
dev.off()

#--------------------------------------------------------------------------------#
#Fig end


### bulk-incucyte - textural covariates 
ibulkCCA <- data.frame()
ibulkCCA <- colnames(ylist_cca[[2]]) %>% as.data.frame(.)
for (i in 1:length(img.list)){
  t1 <- img.list[[i]] %>% as.data.frame(.)
  ibulkCCA <- cbind(ibulkCCA,t1$ibCCA1)
}

#ibulkCCA <- ibulkCCA[1:9]
colnames(ibulkCCA) <- c('txt', paste0('C', 2:8))
rownames(ibulkCCA) <- colnames(ylist_cca[[2]])
#ibulkCCA <- ibulkCCA %>% arrange(C6)

for (c in 1:3){
  ibulkCCA[c+1] <- ibulkCCA[c+1] * (-1)
}

## bulk-incucyte - gene covariates 
gbulkCCA <- data.frame()
gbulkCCA <- colnames(rlogC_t) %>% as.data.frame(.)
for (i in 1:length(bulk.list)){
  t1 <- bulk.list[[i]] %>% as.data.frame(.)
  gbulkCCA <- cbind(gbulkCCA,t1$gbCCA1)
}

#gbulkCCA <- gbulkCCA[2:10]
colnames(gbulkCCA) <- c('gene', paste0('C', 2:8))
rownames(gbulkCCA) <- colnames(rlogC_t)

for (c in 1:3){
  gbulkCCA[c+1] <- gbulkCCA[c+1] * (-1)
}


#GET SHARED GENES FROM ACROSS ALL CONFLUENCY GROUPS
genes1 <- list()
genes2 <- list()
for (c in 1:7){
  top <- gbulkCCA %>% arrange(.[c+1]) %>% .[1:90,]
  #top <- top[apply(top!=0,1,all),]
  top$Group <- rep('high', nrow(top))
  
  bot <- gbulkCCA %>% arrange(.[c+1]) %>% .[(nrow(.)-89):nrow(.),]
  bot$Group <- rep('low', nrow(bot))
  
  genes1[[c]] <- top$gene
  genes2[[c]] <- bot$gene
  
}

genes_top <- do.call('c', genes1) %>% unique(.)
genes_bot <- do.call('c', genes2) %>% unique(.)


for (g in 1:6){
  genes_f <- intersect(genes1[[g]], genes1[[g+1]])
  top_common <- subset(gbulkCCA, gene %in% genes_f)
  genes_f <- intersect(genes2[[g]], genes2[[g+1]])
  bot_common <- subset(gbulkCCA, gene %in% genes_f)
}

top_bot_common <- rbind(top_common, bot_common)

gCCA_tb <- apply(top_bot_common[2:8],2,standardize) %>% 
  as.data.frame(.) %>% 
  cbind(top_bot_common$gene, .) %>% 
  reshape2::melt(.)

colnames(gCCA_tb) <- c('Gene', 'con_grp', 'cca')
gCCA_tb %>% mutate(Gene=fct_relevel(Gene, c(rownames(top), rownames(bot)))) -> gCCA_tb

order_ft <- c(
  "Granularity_1", 
  "Granularity_2",              
  "Granularity_3",              
  "Granularity_4",              
  "Granularity_5",              
  "Granularity_6",              
  "Granularity_7",              
  "Granularity_8",             
  "Granularity_9",    
  "Granularity_10",             
  "Granularity_11",             
  "Granularity_12",             
  "Granularity_13",             
  "Granularity_14",             
  "Granularity_15",             
  "Granularity_16",        
  "Contrast_00",                
  "Correlation_00",             
  "SumAverage_00",              
  "Variance_00",      
  "SumVariance_00",            
  "DifferenceVariance_00",      
  "InverseDifferenceMoment_00", 
  "Entropy_00",                 
  "InfoMeas1_00",               
  "InfoMeas2_00",               
  "SumEntropy_00",              
  "DifferenceEntropy_00",       
  "Angular2ndMoment_00")        

iCCA_melt <- apply(ibulkCCA[2:8],2,standardize) %>% 
  as.data.frame(.) %>% 
  cbind(ibulkCCA$txt, .) %>% 
  reshape2::melt(.)

colnames(iCCA_melt) <- c('Texture', 'con_grp', 'cca')
iCCA_melt %>% mutate(Texture=fct_relevel(Texture, rev(order_ft))) -> iCCA_melt




### FIGURES
###------------------------------------------------------------------#### 
# SUPPLEMENTARY FIGURE 5B
###------------------------------------------------------------------#### 
pdf(file =paste0(figures_filepath, 'Sup_Fig_S5/S5_B.pdf'), width=3, height=15)
ggplot(gCCA_tb, aes(y=Gene, x=con_grp, color=cca, size=cca)) + 
  geom_point(shape=15) + 
  scale_color_continuous_diverging('Tropic', name='cca pixel') + 
  scale_size_area(name='cca gene') +
  theme_tufte() +
  theme(axis.text.y = element_text(size=9.9),
        axis.text.x = element_text(size=12),
        axis.title = element_blank(),
        legend.position = 'none') 
dev.off()
#--------------------------------------------------------------------------------#
#Fig end




### FIGURES
###------------------------------------------------------------------#### 
# SUPPLEMENTARY FIGURE 5D
###------------------------------------------------------------------#### 
pdf(file =paste0(figures_filepath, 'Sup_Fig_S5/S5_D.pdf'), width=3, height=6.6)
ggplot(iCCA_melt, aes(y=Texture, x=con_grp, color=cca, size=cca)) + 
  geom_point(shape=15) + 
  scale_color_continuous_diverging('Tropic', name='cca pixel') + 
  scale_size_area(name='cca pixel') +
  theme_tufte() +
  theme(axis.text.y = element_text(size=6.9),
        axis.title = element_blank(),
        legend.position = 'none') 
dev.off()
#--------------------------------------------------------------------------------#
#Fig end


remove(list=ls())

