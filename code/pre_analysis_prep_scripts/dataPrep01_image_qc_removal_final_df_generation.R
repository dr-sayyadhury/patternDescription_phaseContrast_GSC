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



### FILE/FOLDER PATHS
###------------------------------------------------------------------#### 
### directory paths
path_to_repo <- '/Users/mystique27m/Documents/research/PostdoctoralResearch_2020/Projects/PatternRecognitionInGSCs_UsingCV/'

dataset_path <- paste0(path_to_repo, 'datasets/image_cp/original_output_see_dataPrep01_in_scripts/') ### to run and download again
script_path <- paste0(path_to_repo,'/scripts/')

### file paths
well_filepath <- paste0(path_to_repo,'datasets/wells_to_choose.csv') 


### source data 
###------------------------------------------------------------------#### 
source(paste0(script_path,"source_scripts/sourceData_General_variables_and_functions.R"))

### Original data
###------------------------------------------------------------------#### 
# read cellprofiler results
G523 <- read.csv(paste0(dataset_path,"/G523_mediaOnly_fullDatasetImage.csv"))
G549 <- read.csv(paste0(dataset_path,"/G549_mediaOnly_fullDatasetImage.csv"))
G564 <- read.csv(paste0(dataset_path,"/G564_mediaOnly_fullDatasetImage.csv"))
G566 <- read.csv(paste0(dataset_path,"/G566_mediaOnly_fullDatasetImage.csv"))
G583 <- read.csv(paste0(dataset_path,"/G583_mediaOnly_fullDatasetImage.csv"))
G729 <- read.csv(paste0(dataset_path,"/G729_mediaOnly_fullDatasetImage.csv"))
G797 <- read.csv(paste0(dataset_path,"/G797_mediaOnly_fullDatasetImage.csv"))
G799 <- read.csv(paste0(dataset_path,"/G799_mediaOnly_fullDatasetImage.csv"))
G800 <- read.csv(paste0(dataset_path,"/G800_mediaOnly_fullDatasetImage.csv"))
G837 <- read.csv(paste0(dataset_path,"/G837_mediaOnly_fullDatasetImage.csv"))
G861 <- read.csv(paste0(dataset_path,"/G861_mediaOnly_fullDatasetImage.csv"))
G876 <- read.csv(paste0(dataset_path,"/G876_mediaOnly_fullDatasetImage.csv"))
G851 <- read.csv(paste0(dataset_path,"/G851_mediaOnly_fullDatasetImage.csv"))
G885 <- read.csv(paste0(dataset_path,"/G885_mediaOnly_fullDatasetImage.csv"))
G895 <- read.csv(paste0(dataset_path,"/G895_mediaOnly_fullDatasetImage.csv"))

# read drug list
wells <- read.csv(well_filepath)


#--------------------------------------------------------------------------------#
### --------------------------------- CODE START ------------------------------###
#--------------------------------------------------------------------------------#

imageData <- rbind(G523, G549, G564, 
                   G566, G583, G729, 
                   G797, G799, G800, 
                   G837, G851, G861, 
                   G876, G885, G895)

# remove objects once combined
rm(G523, G549, G564, 
   G566, G583, G729, 
   G797, G799, G800, 
   G837, G851, G861, 
   G876, G885, G895)

###RENAME FEATURE NAMES TO MAKE IT TIDIER AND NEATER
### Make column names shorter
colN <- colnames(imageData)
colN <- gsub('AreaOccupied_', 'Area', colN)
colN <- gsub('Threshold', 'Thr', colN)
colN <- gsub('ExecutionTime', 'ExecT', colN)
colN <- gsub('ImageQuality', 'ImgQly', colN)
colN <- gsub('LocalFocusScore', 'LFS', colN)
colN <- gsub('Intensity', 'Int', colN)
colN <- gsub('Texture_', '', colN)
colN <- gsub('AngularSecondMoment', 'Angular2ndMoment', colN)
colN <- gsub('MaskCropRaw', 'MskRaw', colN)
colN <- gsub('CropRaw', 'cRaw', colN)
colN <- gsub('Metadata_', '', colN)
colN <-gsub('_3_', '_', colN)
colN <- gsub('_256', '', colN)
colN <- gsub('Granularity_cRaw', 'Granularity_3_cRaw', colN)
colN <- gsub('Granularity_MskRaw', 'Granularity_3_MskRaw', colN)

colN[1:3] <- c('Area', 'Perimeter', 'TotalArea')
names(imageData) <- colN

###PREPARE DATA FOR ANALYSIS
###SEPARATE IMAGES ANALYSED FOR QUALITY CONTROL FEATURES (WITH AND WITHOUT MASKING)
final_df <-prepare_data_for_Analysis(imageData)

### masked data
df_mask <- final_df[["maskRaw"]]

#df_mask <- drug_added_to_data(df_mask, drug.list)
df_mask <- cbind(as.data.frame(rep('media', nrow(df_mask))), df_mask)
colnames(df_mask) <- c('Drug', colnames(final_df$maskRaw))

### raw data
df_raw <- final_df$raw
df_raw <- drug_added_to_data(df_raw, wells)

df_mask <- df_mask %>% group_by(Sample) %>% 
  arrange(Sample, group_by=T) %>% 
  mutate(TimePt=factor(TimePt, levels=timePt_levels)) %>%
  as.data.frame(.)


### QC from masked Raw images
qc_maskRaw <- final_df$QC_maskRaw
qc_maskRaw <- qc_maskRaw %>% group_by(Sample) %>% 
  arrange(Sample, group_by=T) %>% 
  mutate(TimePt=factor(TimePt, levels=timePt_levels)) %>%
  as.data.frame(.)

### we remove images that are of low quality using the power log spectrum as described in materials & methods
### we then plot the quality of each image by patient sample and remove the low quality images with low pplog values 


### -------------------------------------------------------------------------###
### === Supplementary figure 2A
### -------------------------------------------------------------------------###

###REMOVE LOW QUALITY IMAGES USING RAW UNMASKED IMAGES
### QC from Raw images
qc_Raw <- final_df$QC_raw
qc_Raw <- qc_Raw %>% group_by(Sample) %>% 
  arrange(Sample, group_by=T) %>% 
  mutate(TimePt=factor(TimePt, levels=timePt_levels)) %>%
  as.data.frame(.)

n=15
g1 <- 
  ggplot(qc_Raw, aes(x=Sample, y=ImgQly_PowerLogLogSlope_cRaw, fill=Sample)) + 
  scale_fill_manual(values=c(col_bold, cols_vivid)) +
  #geom_boxplot(size=.3) +
  geom_violin() +
  geom_jitter(size=0.1, alpha=0.1) +
  scale_y_continuous(breaks=scales::breaks_width(0.1)) +
  theme_linedraw() + 
  theme(axis.text.x=element_text(size=9, face='bold', angle=30), 
        axis.text.y=element_text(size=6, angle=30),
        strip.text = element_text(size=12, colour='#80BA5A', face='bold'),
        legend.position = "none")


qc_Raw_cleaned <- data.frame()
for (i in 1:15){
  levels <- levels(factor(qc_Raw$Sample))
  sample <- qc_Raw[qc_Raw$Sample==levels[i],]
  cleaned <- subset(sample, subset=sample$ImgQly_PowerLogLogSlope_cRaw > quantile(sample$ImgQly_PowerLogLogSlope_cRaw)[2])
  qc_Raw_cleaned <- rbind(qc_Raw_cleaned, cleaned)
}

g2 <- 
  ggplot(qc_Raw_cleaned, aes(x=Sample, y=ImgQly_PowerLogLogSlope_cRaw, fill=Sample)) + 
  scale_fill_manual(values=c(col_bold, cols_vivid)) +
  geom_violin() +
  geom_jitter(size=0.1, alpha=0.1) +
  ylim(-1.2, 2.15) +  
  theme_linedraw() + 
  theme(axis.text.x=element_text(size=9, face='bold', angle=30), 
        axis.text.y=element_blank(), 
        strip.text = element_text(size=12, colour='#80BA5A', face='bold'),
        legend.position = "none")

### Figure was included in manuscript but described in paper

#if(dir.exists(paste0(figures_filepath,'Sup_Fig_S4'))==F){
#  dir.create(paste0(figures_filepath,'Sup_Fig_S4'))}

#pdf(file = paste0(figures_filepath,'Sup_Fig_S4/S4_A.pdf'), width=9, height=4 )
#ggarrange(g1,g2, nrow=1)
#dev.off()
### figure end



### OUTPUT DATA FILE
###------------------------------------------------------------------###
            # S1_outPutDataFile - SAVING CLEANED DATA #
###------------------------------------------------------------------###
final_qc_Mask_cleaned <- subset(qc_maskRaw, subset=rownames(qc_maskRaw) %in% rownames(qc_Raw_cleaned))
colN <- colnames(final_qc_Mask_cleaned) 
colN <- gsub('_MskRaw', '', colN) ### 
colnames(final_qc_Mask_cleaned) <- colN
write.csv(final_qc_Mask_cleaned, paste0(path_to_repo, "datasets/image_cp/final_qc_Mask_cleaned.csv"))


final_dataset_afterQC <- subset(df_mask, subset=rownames(df_mask) %in% rownames(qc_Raw_cleaned)) ### subset images with ilastik masks using whole images that passed QC from above
colN <- colnames(final_dataset_afterQC) 
colN <- gsub('_MskRaw', '', colN) ### 
colnames(final_dataset_afterQC) <- colN
qc_failed <- subset(df_mask, subset= rownames(df_mask) %in% setdiff(rownames(df_mask), rownames(qc_Raw_cleaned)))[2:4] %>% .[c(1,3,2)]

### Save QC failed image dataset
write.csv(qc_failed, paste0(path_to_repo, "datasets/image_cp/mediaOnly_failed_QC.csv"))
### output datafile save end


A=1066104 ### size of full image in pixel.sq
n=15 # sample number

### load imaging datasets
phase <- final_dataset_afterQC
phase <- phase[phase.c.nameOrder] ### reorder columns

### Save QC passed image dataset
s <- levels(factor(phase$Sample))
phase <- phase %>% mutate(TimePt=factor(TimePt, levels=timePt_levels)) %>% arrange(TimePt)
phase$Area <- phase$Area/A #convert absolute area to fraction of total image area
write.csv(phase, paste0(path_to_repo, "datasets/image_cp/final_dataset_afterQC.csv"))


remove(list=ls())
