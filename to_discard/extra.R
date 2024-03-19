#EVALUATING QUALITY MEASURES AGAINST PC1 AND PC2
qc <- read.csv('/results/S1_outPutDataFiles/final_qc_Mask_cleaned.csv', row.names=1)
qc_list <- list()
qc <- qc %>% mutate(TimePt=factor(TimePt, levels=timePt_levels)) %>% arrange(TimePt)

for (l in 1:length(df_list_tp_subsetted)){
  df <- df_list_tp_subsetted[[l]]
  q <- subset(qc, rownames(qc) %in% rownames(df))
  qc_list[[l]] <- q
}

#QC NORMALIZATION AND STANDARDIZATION BY CONFLUENCY GROUPS
pt.size=1.5
new_qc_list <- list()

for (i in 1:9){
  d <- qc_list[[i]]
  new_df <- data.frame()
  #d[is.na(d)] <- 0
  #n=length(levels(factor(d$Sample)))
  for (l in 1:length(levels(factor(d$Sample)))){
    s <- levels(factor(d$Sample))[l]
    df <- subset(d, Sample==s)
    meta <- df[1:3]
    df <- df[4:ncol(df)]
    #df[is.na(df)] <- 0
    df <- apply(df, 2, FUN=function(x){
                                    if(sd(x)==0){
                                      x=x
                                       }else{
                                         standardize(x, remove.na=T)}})### note keep this
    #df <- apply(df, 2, FUN=function(x){x[is.na(x)]<-0
#                                      standardise(x)})
    df <- cbind(meta,df) 
    df <- df %>% .[order(.$TimePt),]

    meta <- df[1:3]
    df <- df[4:ncol(df)]  %>% t(.) %>% as.data.frame(.)
    rownames(meta) <- colnames(df)
    con_grp <- rep(paste0('C',i), nrow(meta))
    meta <- cbind(con_grp, meta)
    df <- cbind(meta, t(df))
    
    #colnames(df) <- c('con_grp', colnames(phase))
    new_df <- rbind(new_df, df)
   
  }
  new_qc_list[[i]] <- new_df
}

qc_list_with_area <- list()
for (i in 1:9){
  q <- new_qc_list[[i]]
  d <- df_list_tp_subsetted[[i]]
  q <- q[rownames(d),]
  df <- cbind(q, d$Area)
  colnames(df) <- c(colnames(new_qc_list[[i]]), 'Area')
  qc_list_with_area[[i]] <- df
}


### Model PC1 against QC paramaters and Area
final_model_qc_pc1 <- data.frame()

for (c in 1:9){
  qc.df <- qc_list_with_area[[c]]
  model_qc_pc1 <- as.data.frame(matrix(ncol=0,nrow=38))
  model_qc_pc1$QCF <- colnames(qc.df)[4:ncol(qc.df)]

  df_ext <- as.data.frame(matrix(ncol=1,nrow=0))  
  colnames(df_ext) <- 'adj.r.sq'
  for (col in 4:ncol(qc.df)){
    d <- cbind(pca_list[[c]]$PC1, qc.df[col])
    colnames(d) <- c('PC1', 'QCF')
    model <- lm(d$PC1~d$QCF, data=d)
    model.summary <- summary(model)
    new_df <- model.summary$adj.r.squared %>% as.data.frame(.)
    #          model.summary$coefficients[2,4]) 
    colnames(new_df) <- 'adj.r.sq'
    df_ext <- rbind(df_ext, new_df)
  }
model_qc_pc1 <- cbind(model_qc_pc1, df_ext)
model_qc_pc1 <- cbind(rep(paste0('C',c),37), model_qc_pc1)
colnames(model_qc_pc1) <- c('ConfluencyGrps', 'QCF', 'adj.r.sq')
final_model_qc_pc1 <- rbind(final_model_qc_pc1, model_qc_pc1)

}

### FIGURES
###------------------------------------------------------------------#### 
                      # SUPPLEMENTARY FIGURE 6C
###------------------------------------------------------------------#### 

pdf(paste0(figures_filepath, '/Sup_Fig_S6/S6_C.pdf'), width=12, height=7)
ggplot(final_model_qc_pc1, aes(QCF, adj.r.sq, fill=QCF)) + 
  theme_tufte() +
  scale_fill_discrete_qualitative() +
  ylim(-0.01,0.2) +
  geom_boxplot() +
  geom_point() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') 
dev.off()

# Fig End

### Model PC2 against QC and Area
final_model_qc_pc2 <- data.frame()

for (c in 1:9){
  qc.df <- qc_list_with_area[[c]]
  model_qc_pc2 <- as.data.frame(matrix(ncol=0,nrow=38))
  model_qc_pc2$QCF <- colnames(qc.df)[4:ncol(qc.df)]
  df_ext <- as.data.frame(matrix(ncol=1,nrow=0))  
  colnames(df_ext) <- 'adj.r.sq'
  for (col in 6:42){
    d <- cbind(pca_list[[c]]$PC2, qc.df[col])
    colnames(d) <- c('PC2', 'QCF')
    model <- lm(d$PC2~d$QCF, data=d)
    model.summary <- summary(model)
    new_df <- model.summary$adj.r.squared %>% as.data.frame(.)
    #          model.summary$coefficients[2,4]) 
    colnames(new_df) <- 'adj.r.sq'
    df_ext <- rbind(df_ext, new_df)
  }
model_qc_pc2 <- cbind(model_qc_pc2, df_ext)
model_qc_pc2 <- cbind(rep(paste0('C',c),37), model_qc_pc2)
colnames(model_qc_pc2) <- c('ConfluencyGrps', 'QCF', 'adj.r.sq')
final_model_qc_pc2 <- rbind(final_model_qc_pc2, model_qc_pc2)

}


### FIGURES
###------------------------------------------------------------------#### 
# SUPPLEMENTARY FIGURE 6D
###------------------------------------------------------------------#### 

pdf(paste0(figures_filepath, '/Sup_Fig_S6/S6D.pdf'), width=12, height=9)
ggplot(final_model_qc_pc2, aes(QCF, adj.r.sq, fill=QCF)) + 
  theme_tufte() +
  scale_fill_discrete_qualitative() +
  geom_boxplot() +
  geom_point() +
  theme(axis.text.x = element_text(angle=60, hjust=1, size=18),
        axis.title.x = element_blank(),
        legend.position = 'none') 
dev.off()


### Model Area against QC
final_model_features_area <- data.frame()

for (c in 1:9){
  pca.df <- pca_list[[c]]
  model_area <- as.data.frame(matrix(ncol=0,nrow=29))
  model_area$Features <- colnames(pca.df)[41:69]

  df_ext <- as.data.frame(matrix(ncol=1,nrow=0))  
  colnames(df_ext) <- 'adj.r.sq'
  for (col in 41:69){
    d <- cbind(pca.df$Area, pca.df[col])
    colnames(d) <- c('Area', 'Features')
    model <- lm(d$Area~d$Features, data=d)
    model.summary <- summary(model)
    new_df <- model.summary$adj.r.squared %>% as.data.frame(.)
    #          model.summary$coefficients[2,4]) 
    colnames(new_df) <- 'adj.r.sq'
    df_ext <- rbind(df_ext, new_df)
  }
model_area <- cbind(model_area, df_ext)
model_area <- cbind(rep(paste0('C',c),29), model_area)
colnames(model_area) <- c('ConfluencyGrps', 'Features', 'adj.r.sq')
final_model_features_area <- rbind(final_model_features_area, model_area)

}

### FIGURES
###------------------------------------------------------------------#### 
# SUPPLEMENTARY FIGURE 6E
###------------------------------------------------------------------#### 

pdf(paste0(figures_filepath, '/Sup_Fig_S6/S6_E.pdf'), width=9, height=6)
ggplot(final_model_features_area , aes(reorder(Features, adj.r.sq), adj.r.sq)) + 
  theme_tufte() +
  geom_boxplot() +
  geom_point() +
  theme(axis.text.x = element_text(angle=60, hjust=1),
        axis.title.x = element_blank( )) 
dev.off()

hl_list <- list()
for (l in 2:8){
  df <- pca_list[[l]]
  df_high <- subset(df, subset=PC2 > 0.7*(max(PC2)-min(PC2))+min(PC2))
  HL <- rep('high', nrow(df_high))
  df_high <- cbind(df_high, HL)

  df_low <- subset(df, subset=PC2 < 0.3*(max(PC2)-min(PC2))+min(PC2))
  HL <- rep('low', nrow(df_low))
  df_low <- cbind(df_low, HL)

  df <- rbind(df_high, df_low)

  hl_list[[l-1]] <- df
}

hl_comb <- do.call(rbind, hl_list)
  
