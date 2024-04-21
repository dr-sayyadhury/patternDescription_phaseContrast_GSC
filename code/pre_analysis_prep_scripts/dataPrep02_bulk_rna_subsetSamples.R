### PACKAGES TO LOAD
###------------------------------------------------------------------#### 
library(SummarizedExperiment)
library(GSVA)
library(dplyr)
library(tidyr)

### FILE/FOLDER PATHS
###------------------------------------------------------------------#### 
### directory paths
path_to_repo <- '/Users/mystique27m/Documents/research/PostdoctoralResearch_2020/Projects/PatternRecognitionInGSCs_UsingCV/'
script_path <- paste0(path_to_repo,'scripts_final/code/')
tables_path <-  paste0(path_to_repo,'/tables/')
dir.create(tables_path)

### read cellprofiler output
###------------------------------------------------------------------#### 
phase <- read.csv(paste0(path_to_repo, '/datasets/image_cp/final_dataset_afterQC.csv'), row.names = 1)


### Import original bulk RNA datasets
###------------------------------------------------------------------#### 
genesets.and.info <- readRDS(paste0(path_to_repo, 'datasets/bulk_RNA/original_dataset_80_samples/genesets_and_info.rds'))
gsc_line <- readRDS(paste0(path_to_repo, 'datasets/bulk_RNA/original_dataset_80_samples/RNA_fpkm_combat_nonSU2C.rds'))

#--------------------------------------------------------------------------------#
### --------------------------------- CODE START ------------------------------###
#--------------------------------------------------------------------------------#
gsc_line <- gsc_line[gsc_line$Type_bis == 'L']
gsc_line <- gsc_line[,(gsc_line$Patient %in% levels(factor(phase$Sample)))]
as.data.frame(SummarizedExperiment::colData(gsc_line))

genesets <- genesets.and.info$gene_set_list
geneset.info <- genesets.and.info$geneset.info

GSC.log2.fpkm <- log2(SummarizedExperiment::assay(gsc_line, 1))

#GSC.ssgsea <- gsva(expr = GSC.log2.fpkm, 
#                   gset.idx.list = genesets,
#                   method = 'ssgsea',
##                   min.sz = 10,
#                   ssgsea.norm = T,
#                   parallel.sz = 0,
#                   verbose = F)

#colnames(GSC.ssgsea) <- gsc_line$Patient

GSC.gsva <- gsva(expr = GSC.log2.fpkm,
                 gset.idx.list = genesets,
                 method = 'gsva',
                 kcdf = 'Gaussian',
                 min.sz = 10,
                 parallel.sz = 0,
                 verbose = F)
colnames(GSC.gsva) <- gsc_line$Patient

#bulk_Dev_I <- read.csv(paste0(path_to_repo,'datasets/bulk_RNA/SourceData3_Bulk_IR_DEV.csv'), header=T) ### THIS SOURCE 
#DATA CAN BE DOWNLOADED FROM THE RICHARDS PAPER FROM NATURE CANCER 2021

#DI_grad_bar <- data.frame(colorbar=c(rep(3,15)), colorlabel= grad_colours)
#bulk_Dev_I <- bulk_Dev_I[c(1,5,6)]
#bulk_Dev_I$SAMPLE_ID <- gsub('_L', '', bulk_Dev_I$SAMPLE_ID) 
#bulk_Dev_I <- subset(bulk_Dev_I, subset=SAMPLE_ID %in% levels(factor(phase$Sample))) ### reduce to matched sample set

#rownames(bulk_Dev_I) <- bulk_Dev_I$SAMPLE_ID
#all(rownames(bulk_Dev_I)==colnames(GSC.gsva))

#bulk_Dev_I <- bulk_Dev_I[colnames(GSC.gsva),]
#all(rownames(bulk_Dev_I)==colnames(GSC.gsva))

#GSC.gsva <- rbind(GSC.gsva, t(bulk_Dev_I[2:3]))

to_remove <- rownames(GSC.gsva)[grep('downreg', rownames(GSC.gsva))]

#GSC.gsva <- GSC.gsva[setdiff(rownames(GSC.gsva), to_remove),]
GSC.gsva <- GSC.gsva[!rownames(GSC.gsva) %in% to_remove,]
nrow(GSC.gsva)
### OUTPUT DATA FILE
###------------------------------------------------------------------###
### S2_outPutDataFile - SAVING CLEANED DATA & Supplementary Table 1
###------------------------------------------------------------------###

### note RNA.GSC c1 = DEvelopmental and RNA.GSC c2 = Injury Response, rename colname

rownames(GSC.gsva)[which(rownames(GSC.gsva)=='RNA.GSC.c1')] <- 'Developmental'
rownames(GSC.gsva)[which(rownames(GSC.gsva)=='RNA.GSC.c2')] <- 'Injury Response'



write.csv(GSC.gsva, paste0(path_to_repo, "datasets/bulk_RNA/GSC.gsva.csv"))
### also presented in paper as Table S1
dir.create(paste0(path_to_repo, 'results/tables/S_Table1'))
write.csv(paste0(path_to_repo, 'results/tables/S_Table1/sT1_GSC_gsva'))


#library(openxlsx)
#detach("package:openxlsx", unload=TRUE)

### prepare to save gene sets used in gsva scoring

geneset_list <- genesets.and.info$gene_set_list
dir.create(paste0(path_to_repo, 'results/tables/S_Table4'))

# Start with an empty data frame to bind everything together

longest_gene_list <- lengths(geneset_list) %>% max()
print(longest_gene_list)

### note RNA.GSC c1 = DEvelopmental and RNA.GSC c2 = Injury Response, rename colname

names(geneset_list)[which(names(geneset_list)=='RNA.GSC.c1')] <- 'Developmental'
names(geneset_list)[which(names(geneset_list)=='RNA.GSC.c2')] <- 'Injury Response'

names(geneset_list)

combined_wide <- data.frame(matrix(nrow = longest_gene_list+1))

for (i in 1:nrow(GSC.gsva)){
    # Get the gene set
    gene_set <- geneset_list[[rownames(GSC.gsva)[i]]]
    # Get the gene set name
    #gene_set <- c(gene_set_name, gene_set)
    #replace NA with ""
    gene_set[is.na(gene_set)] <- ''

    filler <- rep('', longest_gene_list - length(gene_set))
    #gene_set <- c(names(geneset_list)[i], gene_set)
    gene_set <- c(rownames(GSC.gsva)[i], gene_set)
    
    gene_set <- c(gene_set, filler)
    # Add the gene set to the combined data frame
    combined_wide <- cbind(combined_wide, gene_set)
}
ncol(combined_wide)
colnames(combined_wide) <- combined_wide[1,]
combined_wide <- combined_wide[-1,]
combined_wide <- combined_wide[,-1]

ncol(combined_wide)
#colnames(combined_wide) <- names(geneset_list)

### subset to only the genesets used in the analysis
#to_remove <- colnames(combined_wide)[grep('downreg', colnames(combined_wide))]
#combined_paper <- combined_wide[, !colnames(combined_wide) %in% to_remove]

ncol(combined_wide)
nrow(GSC.gsva)

all(rownames(GSC.gsva) == colnames(combined_wide))
### remove columns with NAs
#select top row


# Add an apostrophe to the start of each cell that is a string to force Excel to treat it as text
combined_wide <- combined_wide %>%
  mutate(across(where(is.character), ~paste0(" ", .)))



# Write the combined data frame to a CSV file
write.csv(combined_wide, file = paste0(path_to_repo, 'results/tables/S_Table4/sT4_combined_gene_sets.csv'), quote=F, row.names = F)


remove(list=ls())



