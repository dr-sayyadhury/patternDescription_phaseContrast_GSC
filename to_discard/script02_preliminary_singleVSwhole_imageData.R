### PACKAGES TO LOAD
###------------------------------------------------------------------#### 
library(ggplot2)  
library(colorspace)
library(tidyr)
library(dplyr)
library(ggthemes)
library(ggpubr)


### FILE/FOLDER PATHS
###------------------------------------------------------------------#### 
#path_to_repo <- "/data/"
path_to_repo <- '/Users/tpugh_m04/Documents/Research/Projects/PatternRecognitionInGSCs_UsingCV'
#dataset_path <- paste0(path_to_repo, 'cp_results_singleCell/cp_OutputAnalysis') ### to run and download again
dataset_path <- paste0(path_to_repo, '/datasets/cp_results_singleCell/cp_OutputAnalysis') ### to run and download again


well_filepath <- paste0(path_to_repo,'/DrugLists.csv') 
#dir.create('/results/runOutput')
#figures_filepath <-'/results/'
figures_filepath <- paste0(path_to_repo,'/results/')
#results_path <- '/results/S3_outPutDataFiles/'
results_path <- paste0(figures_filepath,'/S1_outPutDataFiles/')
dir.create(results_path) 

#script_path <- '/code'
script_path <- paste0(path_to_repo,'/scripts')


### ORIGINAL DATA
###------------------------------------------------------------------#### 
source(paste0(script_path,"/sourceData_General_variables_and_functions.R"))

scGrn_wholeImage <- read.csv(paste0(dataset_path, "/singleCell_Image.csv"))
scObj <- read.csv(paste0(dataset_path,"/singleCell_shrinkObjects.csv"))
source(paste0(script_path,"/sourceData_forScript06.R"))


### IMPORT DATA FROM S1
###------------------------------------------------------------------#### 
#phase <- read.csv('/results/S1_outPutDataFiles/final_dataset_afterQC.csv', row.names=1)
phase <- read.csv(paste0(path_to_repo, '/results/S1_outPutDataFiles/final_dataset_afterQC.csv'), row.names = 1)



#--------------------------------------------------------------------------------#
### --------------------------------- CODE START ------------------------------###
#--------------------------------------------------------------------------------#
### subset out common images
phase_d <-  subset(phase, TimePt == 'T_12' | TimePt == 'T_24')
phase_d <- subset(phase_d, Sample %in% levels(factor(scObjF_shortened_txt_grn$Sample)))
phase_d <- cbind(phase_d[1:5], apply(phase_d[6:73], 2, FUN=function(x)(x-min(x))/(max(x)-min(x)))) 

colnames(scObjF_shortened_txt_grn) <- gsub('_orig_3', '', colnames(scObjF_shortened_txt_grn))
colnames(scObjF_shortened_txt_grn) <- gsub('_orig', '', colnames(scObjF_shortened_txt_grn))

colnames(scObjF_shortened_txt_grn) <- gsub('_256', '', colnames(scObjF_shortened_txt_grn))
colnames(scObjF_shortened_txt_grn) <- gsub('Texture_', '', colnames(scObjF_shortened_txt_grn))

colnames(scObjF_shortened_txt_grn) <- gsub('AngularSecondMoment', 'Angular2ndMoment', colnames(scObjF_shortened_txt_grn))

scObjF_shortened_txt_grn<- cbind(scObjF_shortened_txt_grn[1:5], scObjF_shortened_txt_grn[6:73][colnames(phase_d)[6:73]])
scObjF_txt_grn_melt <- reshape2::melt(scObjF_shortened_txt_grn, colnames(scObjF_shortened_txt_grn)[1:5])
scObjF_txt_grn_melt$type <- rep('sc', nrow(scObjF_txt_grn_melt))


phase_melt <- reshape2::melt(phase_d, id=colnames(phase_d[1:5]))
phase_melt$type <- rep('whole', nrow(phase_melt))

whole_sc <- rbind(phase_melt[c(2,6,7,8)], scObjF_txt_grn_melt[c(3,6,7,8)])


### FIGURES
###------------------------------------------------------------------#### 
# FIGURE 1C
###------------------------------------------------------------------####

if(dir.exists(paste0(figures_filepath,'Figure_1'))==F){
  dir.create(paste0(figures_filepath,'Figure_1'))}

pdf(paste0(figures_filepath,'Figure_1/Fig1_G.pdf'), width=9, height=6)
ggpubr::ggline(phase_melt, x='variable', y='value', add='mean_sd', color='Sample', group='Sample') +
  theme_tufte() +
  theme(axis.text.x = element_text(angle=90, size=9, hjust=1),
        axis.title = element_blank(),
        legend.position = 'top') +
  scale_color_manual(values=cols_vivid)
dev.off()

if(dir.exists(paste0(figures_filepath,'Figure_4'))==F){
  dir.create(paste0(figures_filepath,'Figure_4'))}

pdf(paste0(figures_filepath,'Figure_4/Fig4_A_top.pdf'), width=9, height=4.8)
ggpubr::ggline(scObjF_txt_grn_melt, x='variable', y='value', add='mean_sd', color='Sample', group='Sample') +
  theme_tufte() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none') +
  scale_color_manual(values=cols_vivid)
dev.off()


pdf(paste0(figures_filepath,'Figure_4/Fig4_A_bottom.pdf'), width=9, height=4.8)
ggpubr::ggline(whole_sc, x='variable', y='value', add='mean_sd', group='type', color='type')+
  theme_tufte() +
  theme(axis.text.x = element_text(angle=90, size=9, hjust = 1),
        axis.title = element_blank(),
        legend.position = 'none') +
  scale_color_manual(values=c('#FF6E4E','#21B27B'))
dev.off()
###------------------------------------------------------------------####
#Fig Ends

remove(list=ls())
