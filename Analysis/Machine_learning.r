library('tidyverse')
library("phyloseq")
library('qiime2R')
library("ggpubr")
library("vegan")
library("MASS")
library("scales")
library("picante")
library("caret")
library("AppliedPredictiveModeling")
library("ranger")
library("e1071")
library("randomForest")
library("alluvial")

##### Setup, Functions, and loading data #########################################################################################################
# Set seed for analysis
set.seed(10031993)
# Theme, color palette, and working directory
theme_set(theme_pubr())
zoe_palette <- c("gray","#1b9e77", "#7570b3",  "#e6ab02")

setwd("C:/Users/Kenizzer/OneDrive - Donald Danforth Plant Science Center/Grad School/Disseration/Chapter 1/Analysis")

# Function to plot confusion matrix using ggtile plot from a confussion matrix object
# By user: Enrique P?rez Herrero 
# on https://stackoverflow.com/questions/46063234/how-to-produce-a-confusion-matrix-and-find-the-misclassification-rate-of-the-na%C3%AF
ggplotConfusionMatrix <- function(m){
  mytitle <- paste("Accuracy", percent_format()(m$overall[1]),
                   "Kappa", percent_format()(m$overall[2]))
  
  d <- as.data.frame.matrix(m$table)
  drn <- colnames(d)
  drr <- rownames(d)
  drs <- rowSums(d)
  d <- d %>% mutate_if(is.numeric, funs(./drs))
  d <- d %>% gather(x, value)
  Y <- cbind(as.data.frame(m$table), Proportion = d$value)

  p <-
    ggplot(data = Y, aes(x = Reference, y = Prediction, fill= Proportion)) +
    geom_tile( colour = "white") +
    scale_fill_gradient(low = "white", high = "#14A02E", na.value = "white", limits=c(0,1)) +
    ggtitle(mytitle) +
    theme(legend.position = "right") +
    guides(fill = guide_colorbar(frame.colour = "black", ticks = FALSE))
  return(p)
}

MachineLearning_RF_ranger <- function(PHYSEQ_OBJ_1, PHYSEQ_OBJ_2 = FALSE, GROUPING) {
  if (missing(PHYSEQ_OBJ_2)){
    # Remove ASV Table and meta data from phyloseq objects
    # 16s
    ASV.df <- as.data.frame(otu_table(PHYSEQ_OBJ_1))
    ASV_metadata.df <- as.data.frame(sample_data(PHYSEQ_OBJ_1))
    # Format ASV table to be used for machine learning applications and make metadata df
    ASV.df <- t(ASV.df)
    ASV_meta.df <- data.frame(Sample = rownames(ASV_metadata.df), Rootstock = ASV_metadata.df$Rootstock, Tissue = ASV_metadata.df$Tissue, Tissue_Rootstock = paste(ASV_metadata.df$Tissue, ASV_metadata.df$Rootstock, sep = "_"))
    ASV_prefiltered.df <- cbind(ASV.df, ASV_meta.df)
    # ~80:20 split train/test datasets while respecting groups (i.e. sampling the same number of samples from each label)
    # 47 berry, leaf, root, and soil samples
    if ((GROUPING) == "Tissue_Rootstock"){
      train_index <- as.data.frame(ASV_prefiltered.df %>% group_by_(GROUPING) %>% sample_n(9))
    } else {
      train_index <- as.data.frame(ASV_prefiltered.df %>% group_by_(GROUPING) %>% sample_n(36))
    }
    rownames(train_index) <- train_index$Sample
    train_index <- match(rownames(train_index), rownames(ASV_prefiltered.df))
    train_x <- as.data.frame(ASV.df[train_index, ])
    test_y <- as.data.frame(ASV.df[-train_index, ])
    # Train set, 144 samples
    train_x$Sample <- rownames(train_x)
    Training_meta.df <- merge(train_x, ASV_meta.df, by = 'Sample')
    train_x <- subset(Training_meta.df, select = -c(Rootstock, Tissue, Tissue_Rootstock))
    rownames(train_x) <- train_x$Sample
    train_x <- subset(train_x, select = -c(Sample))
    Training_meta.df <- subset(Training_meta.df, select = c(Sample, Rootstock, Tissue, Tissue_Rootstock))
    rownames(Training_meta.df) <- Training_meta.df$Sample 
    # Test set, 40 samples
    test_y$Sample <- rownames(test_y)
    Testing_meta.df <- merge(test_y, ASV_meta.df, by = "Sample")
    test_y <- subset(Testing_meta.df, select = -c(Rootstock, Tissue, Tissue_Rootstock))
    rownames(test_y) <- test_y$Sample
    test_y <- subset(test_y, select = -c(Sample))
    Testing_meta.df <- subset(Testing_meta.df, select = c(Sample, Rootstock, Tissue, Tissue_Rootstock))
    rownames(Testing_meta.df) <- Testing_meta.df$Sample 
    # Training model
    Training_grid <- expand.grid(.mtry = seq(10, length(train_x), round(length(train_x)*0.1)), .splitrule= "gini", .min.node.size = c(1, 5, 10))
    RF_CM <- list()
    RF_CM[["RF_model"]] <- train(x = train_x, y = Training_meta.df[[GROUPING]], method = "ranger", importance = "impurity", tuneGrid = Training_grid, num.trees = 151)
    RF_prediction_3 <- predict(RF_CM[["RF_model"]], test_y)
    RF_CM[["CMatrix"]] <- confusionMatrix(RF_prediction_3, as.factor(Testing_meta.df[[GROUPING]]))
    RF_CM[["CMatrixPLOT"]] <- ggplotConfusionMatrix(RF_CM[["CMatrix"]])
    RF_CM[["VarImporance"]] <- varImp(RF_CM[["RF_model"]])
    RF_CM
  } else {
    # Remove ASV Table and meta data from phyloseq objects
    # 16s
    ASV_16s.df <- as.data.frame(otu_table(PHYSEQ_OBJ_1))
    ASV_16s_metadata.df <- as.data.frame(sample_data(PHYSEQ_OBJ_1))
    # its
    ASV_its.df <- as.data.frame(otu_table(PHYSEQ_OBJ_2))
    ASV_its_metadata.df <- as.data.frame(sample_data(PHYSEQ_OBJ_2))
    # Combine 16s and ITS data into a single matrix and remove samples that are not present both 16s and ITS datasets
    setdiff(rownames(t(ASV_16s.df)), rownames(t(ASV_its.df)))
    ASV_16s_4removed.df <- subset(ASV_16s.df, select = -c(`MV11.3309C.2-3.S`, `MV14.3309C.2-3.R`, `MV15.1103P.10-11.R`, `MV16.3309C.14-15.S`))
    # merge data frames
    Merged_16s_its.df <- rbind(ASV_16s_4removed.df, ASV_its.df)
    # Format ASV table to be used for machine learning applications and make metadata df
    ASV_both.df <- t(Merged_16s_its.df)
    ASV_meta.df <- data.frame(Sample = rownames(ASV_its_metadata.df), Rootstock = ASV_its_metadata.df$Rootstock, Tissue = ASV_its_metadata.df$Tissue, Tissue_Rootstock = paste(ASV_its_metadata.df$Tissue, ASV_its_metadata.df$Rootstock, sep = "_"))
    ASV_both_w_meta.df <- cbind(ASV_both.df, ASV_meta.df)
    # ~80:20 split train/test datasets while respecting groups (i.e. sampling the same number of samples from each label)
    # 47 berry, leaf, root, and soil samples
    if ((GROUPING) == "Tissue_Rootstock"){
      train_index <- as.data.frame(ASV_both_w_meta.df %>% group_by_(GROUPING) %>% sample_n(9))
    } else {
      train_index <- as.data.frame(ASV_both_w_meta.df %>% group_by_(GROUPING) %>% sample_n(36))
    }
    rownames(train_index) <- train_index$Sample
    train_index <- match(rownames(train_index), rownames(ASV_both_w_meta.df))
    train_x <- as.data.frame(ASV_both.df[train_index, ])
    test_y <- as.data.frame(ASV_both.df[-train_index, ])
    # Train set, 144 samples
    train_x$Sample <- rownames(train_x)
    Training_meta.df <- merge(train_x, ASV_meta.df, by = 'Sample')
    train_x <- subset(Training_meta.df, select = -c(Rootstock, Tissue, Tissue_Rootstock))
    rownames(train_x) <- train_x$Sample
    train_x <- subset(train_x, select = -c(Sample))
    Training_meta.df <- subset(Training_meta.df, select = c(Sample, Rootstock, Tissue, Tissue_Rootstock))
    rownames(Training_meta.df) <- Training_meta.df$Sample 
    # Test set, 40 samples
    test_y$Sample <- rownames(test_y)
    Testing_meta.df <- merge(test_y, ASV_meta.df, by = "Sample")
    test_y <- subset(Testing_meta.df, select = -c(Rootstock, Tissue, Tissue_Rootstock))
    rownames(test_y) <- test_y$Sample
    test_y <- subset(test_y, select = -c(Sample))
    Testing_meta.df <- subset(Testing_meta.df, select = c(Sample, Rootstock, Tissue, Tissue_Rootstock))
    rownames(Testing_meta.df) <- Testing_meta.df$Sample 
    #Train model predict, and create confusion matrix
    Training_grid <- expand.grid(.mtry = seq(10, length(train_x), round(length(train_x)*0.1)), .splitrule= "gini", .min.node.size = c(1, 5, 10))
    RF_CM <- list()
    RF_CM[["RF_model"]] <- train(x = train_x, y = Training_meta.df[[GROUPING]], method = "ranger", importance = "impurity", tuneGrid = Training_grid, num.trees = 151)
    RF_prediction_3 <- predict(RF_CM[["RF_model"]], test_y)
    RF_CM[["CMatrix"]] <- confusionMatrix(RF_prediction_3, as.factor(Testing_meta.df[[GROUPING]]))
    RF_CM[["CMatrixPLOT"]] <- ggplotConfusionMatrix(RF_CM[["CMatrix"]])
    RF_CM[["VarImporance"]] <- varImp(RF_CM[["RF_model"]])
    RF_CM
  }
}
#### TESTING
#RF_func_model <- MachineLearning_RF_ranger(physeq_16s, physeq_its, GROUPING = "Rootstock")

# Load QIIME2 objects into a phyloseq class object
physeq_16s <- qza_to_phyloseq('16s/ALL/filtered-table-no-mitochondria-no-chloroplast.qza','16s/ALL/rooted-tree.qza','16s/ALL/taxonomy.qza','16s/ALL/16s_noMockorPosorNeg_metadata.tsv', tmp="C:/tmp")
physeq_its <- qza_to_phyloseq('ITS/ALL/filtered-nocontrol-trimmed-table.qza', taxonomy = 'ITS/ALL/taxonomy.qza', metadata = 'ITS/ALL/its_metadata_noControlsorWine.tsv', tmp="C:/tmp")

# Recode factors OWN to Ungrafted
physeq_16s@sam_data$Rootstock <- recode_factor(physeq_16s@sam_data$Rootstock , OWN = "Ungrafted")
physeq_its@sam_data$Rootstock <- recode_factor(physeq_its@sam_data$Rootstock , OWN = "Ungrafted")

##### 1.0 | 16s and ITS data together, 4 samples removed from 16s, and no preprocessing or prefiltering of ASVs. #####
RF_1.1 <- MachineLearning_RF_ranger(physeq_16s, physeq_its, GROUPING = "Rootstock")
RF_1.2 <- MachineLearning_RF_ranger(physeq_16s, physeq_its, GROUPING = "Tissue")
RF_1.3 <- MachineLearning_RF_ranger(physeq_16s, physeq_its, GROUPING = "Tissue_Rootstock")

##### 2.0 | 16s and ITS data together with prefiltering of ASVs, 4 samples removed from 16s, and no preprocessing. #####
# Filter dataset to eliminate ASVs that are mostly zeros
# Remove taxa that are not found greater than 15 times in 10% of the samples
# 16s: 8199 ASVs => 1084  | ITS: 1195 ASVs => 124
physeq_16s_filtered <- filter_taxa(physeq_16s, function(x) sum(x > 15) > (0.10*length(x)), TRUE)
physeq_its_filtered <- filter_taxa(physeq_its, function(x) sum(x > 15) > (0.10*length(x)), TRUE)
RF_2.1 <- MachineLearning_RF_ranger(physeq_16s_filtered, physeq_its_filtered, GROUPING = "Rootstock")
RF_2.2 <- MachineLearning_RF_ranger(physeq_16s_filtered, physeq_its_filtered, GROUPING = "Tissue")
RF_2.3 <- MachineLearning_RF_ranger(physeq_16s_filtered, physeq_its_filtered, GROUPING = "Tissue_Rootstock")

##### 3.0 | 16s data, no prefiltering or preprocessing of ASVs. #####
RF_3.1 <- MachineLearning_RF_ranger(physeq_16s, GROUPING = "Rootstock")
RF_3.2 <- MachineLearning_RF_ranger(physeq_16s, GROUPING = "Tissue")
RF_3.3 <- MachineLearning_RF_ranger(physeq_16s, GROUPING = "Tissue_Rootstock")

##### 3.0 | 16s data, no preprocessing of ASVs. #####
# Filter dataset to eliminate ASVs that are mostly zeros
# Remove taxa that are not found greater than 15 times in 10% of the samples
# 16s: 8199 ASVs => 1084
physeq_16s_filtered <- filter_taxa(physeq_16s, function(x) sum(x > 15) > (0.10*length(x)), TRUE)
RF_4.1 <- MachineLearning_RF_ranger(physeq_16s_filtered, GROUPING = "Rootstock")
RF_4.2 <- MachineLearning_RF_ranger(physeq_16s_filtered, GROUPING = "Tissue")
RF_4.3 <- MachineLearning_RF_ranger(physeq_16s_filtered, GROUPING = "Tissue_Rootstock")

##### 4.0 | its data, no prefiltering or preprocessing of ASVs. #####
RF_4.1 <- MachineLearning_RF_ranger(physeq_its, GROUPING = "Rootstock")
RF_4.2 <- MachineLearning_RF_ranger(physeq_its, GROUPING = "Tissue")
RF_4.3 <- MachineLearning_RF_ranger(physeq_its, GROUPING = "Tissue_Rootstock")

##### 5.0 | its data, no preprocessing of ASVs. #####
# Filter dataset to eliminate ASVs that are mostly zeros
# Remove taxa that are not found greater than 15 times in 10% of the samples
# 16s: 8199 ASVs => 1084  | ITS: 1195 ASVs => 124
physeq_its_filtered <- filter_taxa(physeq_its, function(x) sum(x > 15) > (0.10*length(x)), TRUE)
RF_5.1 <- MachineLearning_RF_ranger(physeq_its_filtered, GROUPING = "Rootstock")
RF_5.2 <- MachineLearning_RF_ranger(physeq_its_filtered, GROUPING = "Tissue")
RF_5.3 <- MachineLearning_RF_ranger(physeq_its_filtered, GROUPING = "Tissue_Rootstock")