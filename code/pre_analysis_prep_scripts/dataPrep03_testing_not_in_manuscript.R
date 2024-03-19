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
script_path <- paste0(path_to_repo,'/scripts/')

phase <- read.csv(paste0(path_to_repo, '/datasets/image_cp/final_dataset_afterQC.csv'), row.names = 1)
GSC.gsva <- read.csv(paste0(path_to_repo,'/datasets/bulk_RNA/GSC.gsva.csv'), row.names=1)

# general variable script
source(paste0(script_path,"source_scripts/sourceData_General_variables_and_functions.R"))


### Evaluation of pixel features across samples and time points ###
### ------------------------------------------------------------###

#BOXPLOT OF FEATURE SIGNALS PER IMAGE. ALL IMAGES WERE RE-ORDERED BY TIME-POINT TO SHOW THE VARIATION IN SIGNALS WITH TIME AND HENCE CELL GROWTH. 
seq_col <- sequential_hcl(palette='plasma',51, rev=T)
phaseN <- data.frame()
n=15
for (l in 1:n){
  s <- levels(factor(phase$Sample))[l]
  df <- subset(phase, Sample==s)
  meta <- df[1:5]
  df <- df[6:73] 
  df[is.na(df)] <- 0
  df <- apply(df,2,FUN=standardize)
  df <- cbind(meta,df) 
  df <- df %>% .[order(.$TimePt),]
  phaseN <- rbind(phaseN, df)
  
}


### -------------------------------------------------------------------------###
### === Supplementary figure xx
### -------------------------------------------------------------------------###

#PLOT GRANULARITY AND HARALICK FEATURE VECTORS SEPARATELY AND SHOW THEIR VARIATIONS ACROSS IMAGE TIMEPOINTS
t <- phaseN[phaseN$TimePt %in% 'T_0',] %>% as.data.frame(.)
max_confluency <- t %>% group_by(Sample) %>% summarize(.,mean(Area))
max_area <- max_confluency[which(max_confluency$`mean(Area)`==max(max_confluency$`mean(Area)`)),]
max_area <- max_area$`mean(Area)`


#TO SHOW THAT FEATURES THAT BELONG TO THE SAME FAMILY SHOW HIGH CORRELATION WITH EACH OTHER
# ALL SMAPLES COMBINED
cols_grad <- diverging_hcl(palette='Blue-Red', 50)
phaseN[is.na(phaseN)] <- 0


shorter_short <- shorter[22:34] %>%  gsub('_00', '', .)


g_list <- list()
cols <- sequential_hcl(palette='Inferno', n=200, rev=T)
cols <- cols[90:100]
for (t in 1:length(shorter_short)){
  texture <- shorter_short[t]
  txt_fam <- phaseN[ startsWith(colnames(phaseN), texture)] %>% cor(.) %>% as.data.frame(.)
  #txt_fam <- phaseN[,grep(texture, colnames(phaseN))] %>% cor(.) %>% as.data.frame(.)
  print(paste0(texture,' - min:',min(txt_fam)))
  
  txt_fam$txt <- rownames(txt_fam)
  melt_txt <- reshape2::melt(txt_fam, id_var=txt)
  g_l <- ggplot(melt_txt, aes(txt, reorder(value, variable), color=variable)) +
   geom_point(size=6) +
    theme(axis.text=element_text(size=9))
  g_list[[t]] <- g_l
  }

ggarrange(plotlist=g_list, ncol=4, nrow=4)

#calculated by sample


### Figure end -----------------------------------------------------------------
#-------------------------------------------------------------------------------


