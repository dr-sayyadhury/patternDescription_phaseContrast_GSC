
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
dir.create(out)
figures_filepath <- paste0(path_to_repo, 'results/figures/')

tables_path <-  paste0(path_to_repo,'/tables/')

### LOAD DATA
###------------------------------------------------------------------####
phase <- read.csv(paste0(path_to_repo, '/datasets/image_cp/final_dataset_afterQC.csv'), row.names = 1)
source(paste0(script_path,"source_scripts/sourceData_General_variables_and_functions.R")) ### load general variables and functions

#SINCE THERE SEEMS TO BE A DEPENDENCE ON THE AREA ON THE FEATURE SIGNALS, DECIDED TO GROUP IMAGES WITH SIMILAR AVAREAGE AREAS. THIS WAS DETERMINED AFTER TRYING OUT DIFFERENT WAYS OF SPLITTING THE IMAGES.
#THERE WERE 9 BINS , EACH WITH AN INCREMENT OF 10. FOR EXAMPLE THE FIRST BIN WAS ANY IMAGE SET BELONGING TO THE SAME TIME-POINT WITH A MEDIAN AREA OCCUPANCY FALLING WITHIN 0-10%.
#FIRST IMAGES FROM THE SAME SAMPLE WAS SUBSETTED. NEXT THE IMAGES WERE GROUPED TOGETHER BY TIME-POINTS. NEXT THE MEDIAN AREA WAS CALCULATED PER TIME-POINT WAS CALCULATED. FROM THIS, THE IMAGES FALLING WITHIN EACH AREA BIN IN EACH TIME-POINT WAS SELECTED AS A WHOLE  GROUPS WITH SIMILAR CONFLUENCY TOGETHER. 

### reduce samples to smaller set of features
phase <- phase[shorter]
meta <- phase[1:5]
df1 <- phase[!(phase$Sample %in% c('G566','G583')),] %>% as.data.frame(.)

### samples imaged at 12hrly intervals so this was processed separately.
df2 <- phase[(phase$Sample %in% c('G566','G583')),] %>% as.data.frame(.)

### Prepare cellprofiler data for normalization 
df_list_tp_subsetted <- list()

for (confluency in 1:9){
  ### subset data from phase dataset and then perform z-normalization
  tp_AreaMax <- data.frame()
  for (c in 1:length(levels(factor(df1$Sample)))){
    S <- levels(factor(df1$Sample))[c]
    sample <- df1[df1$Sample==S,]
    tp <- sample %>% group_by(TimePt) %>% summarise(Avg=median(Area))
    t <- subset(tp, Avg > confluency*0.1 & Avg < (confluency+1)*0.1)
    s <- subset(sample, TimePt %in% t$TimePt)
    tp_AreaMax <- rbind(tp_AreaMax, s)
  }
  
  for (c in 1:length(levels(factor(df2$Sample)))){
    S <- levels(factor(df2$Sample))[c]
    sample <- df2[df2$Sample==S,]
    tp <- sample %>% group_by(TimePt) %>% summarise(Avg=median(Area))
    t <- subset(tp, Avg > confluency*0.1 & Avg < (confluency+1)*0.1)
    s <- subset(sample, TimePt %in% t$TimePt)
    tp_AreaMax <- rbind(tp_AreaMax, s)
  }
  
  df_list_tp_subsetted[[confluency]] <- tp_AreaMax
  
}

df_list_tp_subsetted[[1]] <- subset(df_list_tp_subsetted[[1]], Sample!='G800')

### plot the distribution of the confluencies
glist <- list()

for (i in 1:length(df_list_tp_subsetted)){  
      glist[[i]] <- df_list_tp_subsetted[[i]] %>% 
        ggplot(aes(x=Sample, y=Area, fill=Sample)) + 
        geom_boxplot() + 
        geom_point(size=0.1) +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12),
              axis.text.y = element_text(size=12),
              axis.title.x = element_text(size=18),
              axis.title.y = element_text(size=18), 
              legend.position='none') + 
        labs(x='Sample', y='Area')
}


###------------------------------------------------------------------####
### FIGURE S2A
###------------------------------------------------------------------####
pdf(file = paste0(figures_filepath,'Sup_Fig_S2/S2_A.pdf'), width=42, height=6)
ggarrange(plotlist=glist, ncol=9, nrow=1)
dev.off()

###------------------------------------------------------------------####
### SAVE INTERMEDIATE DATA
###------------------------------------------------------------------####
dir.create(paste0(out, 'script01_1_image_data_by_confluencies/'))
new_path <- paste0(out, 'script01_1_image_data_by_confluencies/')
#new_path <- paste0(new_dir, 'C')
for (l in 1:length(df_list_tp_subsetted)){
  write.csv(df_list_tp_subsetted[[l]], paste0(new_path, l, '.csv'))
}

#### save nromalized data and proceed with actual biological analysis
saveRDS(df_list_tp_subsetted, paste0(out, '/df_list_tp_subsetted.rds'))

remove(list=ls())




