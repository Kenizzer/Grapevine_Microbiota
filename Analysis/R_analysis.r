library('tidyverse')
library("ggpubr")
library("car")
library("vegan")
library("phyloseq")
library('qiime2R')
library("DESeq2")
library("corncob")
library("ALDEx2")
library("anomalize")
library("emmeans")

# Theme set and Color Palette
theme_set(theme_pubr())
zoe_palette <- c("gray","#1b9e77", "#7570b3",  "#e6ab02")

# Functions
Alpha_div_metrics <- function(phyloseq_obj, its_or_16s) {
  if (its_or_16s == "16s"){
    # Calculate alpha diversity metrics and put them into a dataframe with the sample metadata
    Alpha_16s <- estimate_richness(phyloseq_obj)
    R <- picante::pd(samp = t(otu_table(phyloseq_obj)), tree = phy_tree(phyloseq_obj), include.root = FALSE)
    Alpha_16s <- cbind(Alpha_16s, "Faithpd" = R$PD, as.data.frame(phyloseq_obj@sam_data))
  } else if(its_or_16s == "its"){
    # Calculate alpha diversity metrics and put them into a dataframe with the sample metadata
    Alpha_its <- estimate_richness(phyloseq_obj)
    Alpha_its <- cbind(Alpha_its, as.data.frame(phyloseq_obj@sam_data))
  }
}

PLOT_BP_ITS_ROOTSTOCK <- function(plot_data) {
  ggplot(plot_data, aes(Tissue, y= value, fill = variable))+
    geom_bar(stat = "identity", color = "black") +
    ylab("Relative abundance") +
    theme(legend.position = "right", legend.key.size = unit(0.75, "cm"), legend.text=element_text(size=12), legend.title=element_text(size=24)) +
    scale_fill_manual(name ="Fungal class", values=taxapallete, labels = c(expression(italic("Dothideomycetes")), expression(italic("Leotiomycetes")), expression(italic("Pezizomycetes")), expression(italic("Saccharomycetes")), expression(italic("Sordariomycetes")), expression(italic("Agaricomycetes")), expression(italic("Cystobasidiomycetes")), expression(italic("Microbotryomycetes")), expression(italic("Tremellomycetes")), expression(italic("Glomeromycetes")), expression(italic("Paraglomeromycetes")), expression(italic("Mortierellomycetes")), "Low abundance taxa")) +
    theme(legend.text.align = 0)
}

PLOT_BP_16S_ROOTSTOCK <- function(plot_data) {
  ggplot(plot_data, aes(Tissue, y= value, fill = variable))+
    geom_bar(stat = "identity", color = "black") +
    ylab("Relative abundance") +
    theme(legend.position = "right", legend.key.size = unit(0.75, "cm"), legend.text=element_text(size=12), legend.title=element_text(size=24)) +
    scale_fill_manual(name ="Bacterial phylum", values=taxapallete, labels = c(expression(italic("Acidobacteria")), expression(italic("Actinobacteria")), expression(italic("Armatimonadetes")), expression(italic("Bacteroidetes")), expression(italic("Chloroflexi")), expression(italic("Deinococcus-Thermus")), expression(italic("Firmicutes")), expression(italic("Gemmatimonadetes")), expression(italic("Latescibacteria")), expression(italic("Nitrospirae")), expression(italic("Planctomycetes")), expression(italic("Proteobacteria")), expression(italic("Rokubacteria")), expression(italic("Verrucomicrobia")), "Low abundance taxa")) +
    theme(legend.text.align = 0)
}

PLOT_PCoA <- function(plot_data, distance_matrix, axis1, axis2) {
  temp <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2))
  plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2)) +
    geom_point(aes(fill=Rootstock, shape=TissueType), size = 3) +
    scale_shape_manual(values=c(24, 22, 23, 21)) +  
    scale_fill_manual(name = "Rootstock", values=zoe_palette) +
    labs(shape= "Tissue", color= "Rootstock") +
    xlab(paste("PCoA", axis1, sub(".*\\ ", "", temp$labels$x))) +
    ylab(paste("PCoA", axis2, sub(".*\\ ", "", temp$labels$y))) +
    guides(fill = guide_legend(override.aes = list(shape = 21)), shape = guide_legend(override.aes = list(fill = "black")))
}

PLOT_OTU_corncob_anomalized <- function(otu_number, taxa_lvl = 3, grp_by = "Rootstock",  not_plot1 = FALSE, not_plot2 = FALSE, not_plot3 = FALSE) {
  # Function that will given an *otu_number* make a plot that is divded by tissue and rootstock (or other factor if specifed in grp_by).
  # Before plotting the data column is run through anomalize to remove outliers that interfer with plot interpretation.
  # Will assign the title of the plot as user specified taxonomic level (default is 3, i.e. Class)
  # 1 = Kingdom | 2 = Phylum | 3 = Class | 4 = Order | 5 = Family | 6 = Genus | 7 = Species
  # Additionally the user can control which tissues to display by providing the "name" of tissues to remove before plotting (default is plot all).
  
  # Create dataframe and coerce it into the correct form
  TEMP <- as.data.frame(physeq_its@otu_table[otu_number,])
  TEMP_long <- t(TEMP)
  TEMP_long_2 <- cbind(TEMP_long, physeq_its@sam_data)
  # Removing outliers that are over 5 standard devations away from the mean
  column_name <- colnames(TEMP_long_2[1])
  x <- anomalize(as_tibble(TEMP_long_2), target=column_name, method='iqr', alpha=0.03, max_anoms=0.1)
  d <- x[x$anomaly != 'Yes',]
  d <- as.data.frame(d)
  # Generate varible to hold ggtitle
  X <- physeq_its@tax_table[otu_number,]
  X <- as.data.frame(X)
  if (taxa_lvl == 1){
    Y <- paste("ASV", otu_number, " -", " Kingdom:", X[1,1], sep = "")
  } else if(taxa_lvl == 2){
    Y <- paste("ASV", otu_number, " -"," Phylum:",X[1,2], sep = "")
  } else if(taxa_lvl == 3){
    Y <- paste("ASV", otu_number, " -"," Class:",X[1,3], sep = "")
  } else if(taxa_lvl == 4){
    Y <- paste("ASV", otu_number, " -"," Order:",X[1,4], sep = "")
  } else if(taxa_lvl == 5){
    Y <- paste("ASV", otu_number, " -"," Family:",X[1,5], sep = "")
  } else if(taxa_lvl == 6){
    Y <- paste("ASV", otu_number, " -"," Genus:",X[1,6], sep = "")
  } else if(taxa_lvl == 7){
    Y <- paste("ASV", otu_number, " -"," Species:",X[1,7], sep = "")
  }
  # Plot
  # Plots all tissues.
  if (not_plot1 == FALSE){
    ggplot(d, aes_string("TissueType", d[,1], fill= grp_by)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette) + scale_y_continuous(name="Abundance") + xlab("Tissue") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(Y) 
  } else if (not_plot1 != FALSE){
    # Plots tissues that are no specified to be removed from the dataframe prior to plotting.
    d <- d[which(d$TissueType != not_plot1 & d$TissueType != not_plot2 & d$TissueType != not_plot3),]
    ggplot(d, aes_string("TissueType", d[,1], fill= grp_by)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette) + scale_y_continuous(name="Abundance") + xlab("Tissue") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(Y) 
  }
}

# Given a phyloseq object, plot the ASV that was detemined by DESeq2 to be differentailly
# Abundant
plot_deseq2_DiffAbunMicob <- function(physeq_obj, ASV_number){
  # Create dataframe and coerce it into the correct form
  TEMP <- as.data.frame(physeq_its@otu_table[ASV_number,])
  TEMP_long <- t(TEMP)
  TEMP_long_2 <- cbind(TEMP_long, physeq_its@sam_data)
  # Removing outliers that are over 5 standard devations away from the mean
  column_name <- colnames(TEMP_long_2[1])
  x <- anomalize(as_tibble(TEMP_long_2), target=column_name, method='iqr', alpha=0.03, max_anoms=0.1)
  d <- x[x$anomaly != 'Yes',]
  d <- as.data.frame(d)
  # Generate varible to hold ggtitle
  ggplot(d, aes(Tissue, d[,1], fill= Rootstock)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette) + scale_y_continuous(name="Abundance") + xlab("Tissue") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22))
}

# Given a phyloseq object, return a dataframe of the sample metadata 
# From: https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# Given a phyloseq object, plot the ASV that was detemined by DESeq2 to be differentailly
# Abundant
plot_deseq2_DiffAbunMicob <- function(physeq_obj, ASV_number){
  
  
  temp_sample_tab <- pssd2veg(physeq_obj)
  otu_matrix <- as(otu_table(physeq_obj), "matrix")
  TEMP <- data.frame(ASV_count = otu_matrix[ASV_number,])
  TEMP2<- cbind(temp_sample_tab, TEMP)
  ggplot(TEMP2, aes(Tissue, ASV_count, fill= Rootstock)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette) + scale_y_continuous(name="Abundance") + xlab("Tissue") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22))
}


# Home PC Windows        setwd("C:/Users/Kenizzer/OneDrive - Donald Danforth Plant Science Center/Grad School/Disseration/Chapter 1/Analysis")
# Laptop                 setwd("C:/Users/Joel/OneDrive - Donald Danforth Plant Science Center/Grad School/Disseration/Chapter 1")
# Office PC linux        setwd("/media/kenizzer/Storage disk/argonne_sequencing/")

##### 1.0) RAREFIED ANALYSES #####
#Bacterial samples (i.e. 16s) rarefied to 1500 reads
#Fungal samples (i.e. ITS) rarefied to 5000 reads

##### 1.1) Alpha Diversity plots #####
# Load QIIME2 objects into a phyloseq class object.
physeq_16s <- qza_to_phyloseq(features = '16s/ALL/filtered-table-no-mitochondria-no-chloroplast.qza', tree = '16s/ALL/rooted-tree.qza', taxonomy = '16s/ALL/taxonomy.qza', metadata = '16s/ALL/16s_noMockorPosorNeg_metadata.tsv', tmp="C:/tmp")
physeq_its <- qza_to_phyloseq(features = 'ITS/ALL/filtered-nocontrol-trimmed-table.qza', taxonomy = 'ITS/ALL/taxonomy.qza', metadata = 'ITS/ALL/its_metadata_noControlsorWine.tsv', tmp="C:/tmp")

# Recode factors OWN to Ungrafted.
physeq_16s@sam_data$Rootstock <- recode_factor(physeq_16s@sam_data$Rootstock, OWN = "Ungrafted")
physeq_its@sam_data$Rootstock <- recode_factor(physeq_its@sam_data$Rootstock, OWN = "Ungrafted")
physeq_16s@sam_data$Irrigation <- recode_factor(physeq_16s@sam_data$Irrigation, none = "None", rdi = "RDI", full = "Full")
physeq_its@sam_data$Irrigation <- recode_factor(physeq_its@sam_data$Irrigation, none = "None", rdi = "RDI", full = "Full")

# Rarefy the samples to 1500 (required as this was not done in previously in QIIME).
rare_physeq_16s <- rarefy_even_depth(physeq_16s, sample.size = 1500, replace = FALSE, rngseed = 10031993) 
### Removing 7 biological samples
# Rarefy the samples to 5000 (required as this was not done in previously in QIIME).
rare_physeq_its <- rarefy_even_depth(physeq_its, sample.size = 5000, replace = FALSE, rngseed = 10031993) 
### Removing 7 biological samples

# Calculate alpha diversity metrics for 16s and ITS.
ALPHA_div_16s_rare <- Alpha_div_metrics(rare_physeq_16s, its_or_16s = "16s")
ALPHA_div_its_rare <- Alpha_div_metrics(rare_physeq_its, its_or_16s = "its")

# This tibble is used to place the significance letters at a consist hieght over each barplot
labels_df <- tibble(Tissue=levels(ALPHA_div_16s_rare$Tissue), Mfphd=max(ALPHA_div_16s_rare$Faithpd) * 1.2, Mobso=max(ALPHA_div_16s_rare$Observed) * 1.2, Mshan=max(ALPHA_div_16s_rare$Shannon) * 1.2, Msimpinv = max(ALPHA_div_16s_rare$InvSimpson) * 1.2, Mobso_i=max(ALPHA_div_its_rare$Observed) * 1.2, Mshan_i=max(ALPHA_div_its_rare$Shannon) * 1.2, Msimpinv_i = max(ALPHA_div_its_rare$InvSimpson) * 1.2)

# Make plots 16S by tissue
fphd_16s <- ggplot(ALPHA_div_16s_rare, aes(Tissue, Faithpd)) + geom_boxplot(outlier.shape = NA) + ggtitle("16S") + theme(plot.title = element_text(hjust = 0.5)) + geom_jitter(width = 0.25) + xlab ("Tissue") + ylab("Faith's phylogenetic diversity") + geom_text(data=labels_df, aes(Tissue, Mfphd, label=c("a","a","b","c")), size = 6)
obso_16s <- ggplot(ALPHA_div_16s_rare, aes(Tissue, Observed)) + geom_boxplot(outlier.shape = NA) + ggtitle("16S") + theme(plot.title = element_text(hjust = 0.5)) + geom_jitter(width = 0.25) + xlab ("Tissue") + ylab("Observed ASVs") + geom_text(data=labels_df, aes(Tissue, Mobso, label=c("a","a","b","c")), size = 6)
simpI_16s <- ggplot(ALPHA_div_16s_rare, aes(Tissue, InvSimpson)) + geom_boxplot(outlier.shape = NA) + ggtitle("16S") + theme(plot.title = element_text(hjust = 0.5)) + geom_jitter(width = 0.25) + xlab ("Tissue") + ylab(expression("Simpson's D"^-1)) + geom_text(data=labels_df, aes(Tissue, Msimpinv, label=c("a","a","b","c")), size = 6)
shan_16s <- ggplot(ALPHA_div_16s_rare, aes(Tissue, Shannon)) + geom_boxplot(outlier.shape = NA) + ggtitle("16S") + theme(plot.title = element_text(hjust = 0.5)) + geom_jitter(width = 0.25) + xlab ("Tissue") + ylab("Shannon's index") + geom_text(data=labels_df, aes(Tissue, Mshan, label=c("a","a","b","c")), size = 6)
# Tukey tests 16s for significance letters means in alpha diversity charts
TukeyHSD(aov(Faithpd ~ Tissue, data = ALPHA_div_16s_rare), conf.level = 0.95)
TukeyHSD(aov(Observed ~ Tissue, data = ALPHA_div_16s_rare), conf.level = 0.95)
TukeyHSD(aov(InvSimpson ~ Tissue, data = ALPHA_div_16s_rare), conf.level = 0.95)
TukeyHSD(aov(Shannon ~ Tissue, data = ALPHA_div_16s_rare), conf.level = 0.95)
# Make plots ITS by Tissue
obso_its <- ggplot(ALPHA_div_its_rare, aes(Tissue, Observed)) + geom_boxplot(outlier.shape = NA) + ggtitle("ITS") + theme(plot.title = element_text(hjust = 0.5)) + geom_jitter(width = 0.25) + xlab ("Tissue") + ylab("Observed ASVs") + geom_text(data=labels_df, aes(Tissue, Mobso_i, label=c("a","b","c","d")), size = 6)
simpI_its <- ggplot(ALPHA_div_its_rare, aes(Tissue, InvSimpson)) + geom_boxplot(outlier.shape = NA) + ggtitle("ITS") + theme(plot.title = element_text(hjust = 0.5)) + geom_jitter(width = 0.25) + xlab ("Tissue") + ylab(expression("Simpson's D"^-1)) + geom_text(data=labels_df, aes(Tissue, Msimpinv_i, label=c("a","b","c","d")), size = 6)
shan_its <- ggplot(ALPHA_div_its_rare, aes(Tissue, Shannon)) + geom_boxplot(outlier.shape = NA) + ggtitle("ITS") + theme(plot.title = element_text(hjust = 0.5)) + geom_jitter(width = 0.25) + xlab ("Tissue") + ylab("Shannon's index") + geom_text(data=labels_df, aes(Tissue, Mshan_i, label=c("a","b","c","d")), size = 6)
# Tukey tests ITS for significance letters means in alpha diversity charts
TukeyHSD(aov(Observed ~ Tissue, data = ALPHA_div_its_rare), conf.level = 0.95)
TukeyHSD(aov(InvSimpson ~ Tissue, data = ALPHA_div_its_rare), conf.level = 0.95)
TukeyHSD(aov(Shannon ~ Tissue, data = ALPHA_div_its_rare), conf.level = 0.95)

# Make plots 16S by rootstock and tissue
fphd_byR_T_16s <- ggplot(ALPHA_div_16s_rare, aes(Tissue, Faithpd, fill = Rootstock)) + geom_boxplot(outlier.shape = NA) + ggtitle("16S") + theme(plot.title = element_text(hjust = 0.5)) + xlab ("Tissue") + ylab("Faith's phylogenetic diversity") + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette)
simpI_byR_T_16s <- ggplot(ALPHA_div_16s_rare, aes(Tissue, InvSimpson, fill = Rootstock)) + geom_boxplot(outlier.shape = NA) + ggtitle("16S") + theme(plot.title = element_text(hjust = 0.5)) + xlab ("Tissue") + ylab(expression("Simpson's D"^-1)) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette)
obso_byR_T_16s <- ggplot(ALPHA_div_16s_rare, aes(Tissue, Observed, fill = Rootstock)) + geom_boxplot(outlier.shape = NA) + ggtitle("16S") + theme(plot.title = element_text(hjust = 0.5)) + xlab ("Tissue") + ylab("Observed ASVs") + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette)
shan_byR_T_16s <- ggplot(ALPHA_div_16s_rare, aes(Tissue, Shannon, fill = Rootstock)) + geom_boxplot(outlier.shape = NA) + ggtitle("16S") + theme(plot.title = element_text(hjust = 0.5)) + xlab ("Tissue") + ylab("Shannon's index") + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette)
# Make plots ITS by rootstock and tissue
simpI_byR_T_its <- ggplot(ALPHA_div_its_rare, aes(Tissue, InvSimpson, fill = Rootstock)) + geom_boxplot(outlier.shape = NA) + ggtitle("ITS") + theme(plot.title = element_text(hjust = 0.5)) + xlab ("Tissue") + ylab(expression("Simpson's D"^-1)) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette)
obso_byR_T_its <- ggplot(ALPHA_div_its_rare, aes(Tissue, Observed, fill = Rootstock)) + geom_boxplot(outlier.shape = NA) + ggtitle("ITS") + theme(plot.title = element_text(hjust = 0.5)) + xlab ("Tissue") + ylab("Observed ASVs") + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette)
shan_byR_T_its <- ggplot(ALPHA_div_its_rare, aes(Tissue, Shannon, fill = Rootstock)) + geom_boxplot(outlier.shape = NA) + ggtitle("ITS") + theme(plot.title = element_text(hjust = 0.5)) + xlab ("Tissue") + ylab("Shannon's index") + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette)

##### 1.2) Linear models #####
# Set contrasts to deviation coding for factors
options(contrasts = c("contr.sum", "contr.sum"))
# Assign levels to be used as reference
ALPHA_div_16s_rare$Rootstock <- relevel(ALPHA_div_16s_rare$Rootstock, ref = "Ungrafted")
ALPHA_div_16s_rare$Irrigation <- relevel(ALPHA_div_16s_rare$Irrigation, ref = "None")
ALPHA_div_16s_rare$Tissue <- relevel(ALPHA_div_16s_rare$Tissue, ref = "Soil")
ALPHA_div_its_rare$Rootstock <- relevel(ALPHA_div_its_rare$Rootstock, ref = "Ungrafted")
ALPHA_div_its_rare$Irrigation <- relevel(ALPHA_div_its_rare$Irrigation, ref = "None")
ALPHA_div_its_rare$Tissue <- relevel(ALPHA_div_its_rare$Tissue, ref = "Soil")

# linear models 16S FULL
fphd_16s_fit_full <- lm(Faithpd ~ Tissue*Rootstock*Irrigation + Block, data = ALPHA_div_16s_rare)
simpI_16s_fit_full <- lm(InvSimpson ~ Tissue*Rootstock*Irrigation + Block, data = ALPHA_div_16s_rare)
obso_16s_fit_full <- lm(Observed ~ Tissue*Rootstock*Irrigation + Block, data = ALPHA_div_16s_rare)
shan_16s_fit_full <- lm(Shannon ~ Tissue*Rootstock*Irrigation + Block, data = ALPHA_div_16s_rare)

as.data.frame(Anova(fphd_16s_fit_full, type = "III"))
as.data.frame(Anova(simpI_16s_fit_full, type = "III"))
as.data.frame(Anova(obso_16s_fit_full, type = "III"))
as.data.frame(Anova(shan_16s_fit_full, type = "III"))

#linear models ITS FULL
simpI_its_fit_full <- lm(InvSimpson ~ Tissue*Rootstock*Irrigation + Block, data = ALPHA_div_its_rare)
obso_its_fit_full <- lm(Observed ~ Tissue*Rootstock*Irrigation + Block, data = ALPHA_div_its_rare)
shan_its_fit_full <- lm(Shannon ~ Tissue*Rootstock*Irrigation + Block, data = ALPHA_div_its_rare)

as.data.frame(Anova(simpI_its_fit_full, type = "III"))
as.data.frame(Anova(obso_its_fit_full, type = "III"))
as.data.frame(Anova(shan_its_fit_full, type = "III"))

# linear models 16S optimized
fphd_16s_fit <- lm(Faithpd ~ Tissue + Rootstock + Irrigation + Block, data = ALPHA_div_16s_rare)
simpI_16s_fit <- lm(InvSimpson ~ Tissue + Rootstock + Irrigation + Block, data = ALPHA_div_16s_rare)
obso_16s_fit <- lm(Observed ~ Tissue + Rootstock + Irrigation + Block, data = ALPHA_div_16s_rare)
shan_16s_fit <- lm(Shannon ~ Tissue + Rootstock + Irrigation + Block, data = ALPHA_div_16s_rare)

as.data.frame(Anova(fphd_16s_fit, type = "II"))
as.data.frame(Anova(simpI_16s_fit, type = "II"))
as.data.frame(Anova(obso_16s_fit, type = "II"))
as.data.frame(Anova(shan_16s_fit, type = "II"))

#linear models ITS optimized
simpI_its_fit <- lm(InvSimpson ~ Tissue + Rootstock + Irrigation + Block, data = ALPHA_div_its_rare)
obso_its_fit <- lm(Observed ~ Tissue + Rootstock + Irrigation + Block, data = ALPHA_div_its_rare)
shan_its_fit <- lm(Shannon ~ Tissue + Rootstock + Irrigation + Block, data = ALPHA_div_its_rare)

as.data.frame(Anova(simpI_its_fit, type = "II"))
as.data.frame(Anova(obso_its_fit, type = "II"))
as.data.frame(Anova(shan_its_fit, type = "II"))

#Reset contrasts
options(contrasts = c("contr.treatment", "contr.poly"))

##### 1.3) BETA DIVERSITY ANALYSIS #####
# Calculate weighted and unweighted unifrac distances between pairs of samples
out.wunifrac <- ordinate(rare_physeq_16s, method = "MDS", distance = "wunifrac")
out.uunifrac <- ordinate(rare_physeq_16s, method = "MDS", distance = "uunifrac")
out.jaccard  <- ordinate(rare_physeq_16s, method = "MDS", distance = "jaccard")

#Calculate distance matrixes for use in Vegan
out.dist.wunifrac <- phyloseq::distance(rare_physeq_16s, "wunifrac")
out.dist.uunifrac <- phyloseq::distance(rare_physeq_16s, "uunifrac")
out.dist.jaccard  <- phyloseq::distance(rare_physeq_16s, "jaccard", binary = TRUE)

# get variance explained by axes 1-3
sum(out.wunifrac$values$Eigenvalues[1:3])/sum(out.wunifrac$values$Eigenvalues) #72%
sum(out.uunifrac$values$Eigenvalues[1:3])/sum(out.uunifrac$values$Eigenvalues) #41%
sum(out.jaccard$values$Eigenvalues[1:3])/sum(out.jaccard$values$Eigenvalues) #38%

# Plot PCAs 1x2 and 1x3 for weighted and unweighted unifrac 
PCoA_1_2_wunif_16s <- PLOT_PCoA(physeq_16s, out.wunifrac, 1, 2)
PCoA_1_3_wunif_16s <- PLOT_PCoA(physeq_16s, out.wunifrac, 1, 3)
PCoA_1_2_uunif_16s <- PLOT_PCoA(physeq_16s, out.uunifrac, 1, 2)
PCoA_1_3_uunif_16s <- PLOT_PCoA(physeq_16s, out.uunifrac, 1, 3)
PCoA_1_2_jaccd_16s <- PLOT_PCoA(physeq_16s, out.jaccard, 1, 2)
PCoA_1_3_jaccd_16s <- PLOT_PCoA(physeq_16s, out.jaccard, 1, 3)

# Create sample data dataframe from phyloseq object to use in vegan/adonis
physeq_metadata_16s <- pssd2veg(rare_physeq_16s)

# Adonis testing on weighted and unweighted unifrac distance matrix generated using phyloseq::distance
# Full models
adonis2(out.dist.wunifrac ~ Tissue*Rootstock*Irrigation + Block, data = physeq_metadata_16s, permutations = 10000)
adonis2(out.dist.uunifrac ~ Tissue*Rootstock*Irrigation + Block, data = physeq_metadata_16s, permutations = 10000)
## Optimized model
adonis2(out.dist.wunifrac ~ Rootstock + Tissue + Irrigation + Block, data = physeq_metadata_16s, permutations = 10000)
adonis2(out.dist.uunifrac ~ Rootstock + Tissue + Irrigation + Block, data = physeq_metadata_16s, permutations = 10000)

# ITS beta diversity calculations
# Calculate distances between pairs of samples
out.bray_its <- ordinate(rare_physeq_its, method = "MDS", distance = "bray")
out.jaccard_its <- ordinate(rare_physeq_its, method = "MDS", distance = "jaccard")

#Calculate distance matrixes for use in Vegan
out.dist.bray_its <- phyloseq::distance(rare_physeq_its, "bray")
out.dist.jaccard_its <- phyloseq::distance(rare_physeq_its, "jaccard", binary = TRUE)

# get variance explained by axes 1-3
sum(out.bray_its$values$Eigenvalues[1:3])/sum(out.bray_its$values$Eigenvalues) #59%
sum(out.jaccard_its$values$Eigenvalues[1:3])/sum(out.jaccard_its$values$Eigenvalues) #45%

#PCoA plots
PCoA_1_2_bray_its <- PLOT_PCoA(physeq_its, out.bray_its, 1, 2)
PCoA_1_3_bray_its <- PLOT_PCoA(physeq_its, out.bray_its, 1, 3)
PCoA_1_2_jacc_its <- PLOT_PCoA(physeq_its, out.jaccard_its, 1, 2)
PCoA_1_3_jacc_its <- PLOT_PCoA(physeq_its, out.jaccard_its, 1, 3)

# Create sample data dataframe from phyloseq object to use in vegan/adonis
physeq_metadata_its <- pssd2veg(rare_physeq_its)

# Adonis testing on weighted and unweighted unifrac distance matrix generated using
# Full models
adonis2(out.dist.bray_its ~ Tissue*Rootstock*Irrigation + Block, data = physeq_metadata_its, permutations = 10000)
adonis2(out.dist.jaccard_its ~ Tissue*Rootstock*Irrigation + Block, data = physeq_metadata_its, permutations = 10000)
## Optomized model (stepwise)
adonis2(out.dist.bray_its ~ Rootstock + Tissue + Irrigation + Block + Rootstock:Tissue, data = physeq_metadata_its, permutations = 10000)
adonis2(out.dist.jaccard_its ~ Rootstock + Tissue + Irrigation + Block, data = physeq_metadata_its, permutations = 10000)

##### 2.0) Taxonomic Barplots #####
##### 2.1) 16S #####
# Color pallette for ggplot2 taxonomic barplots
taxapallete <- c('#8c8fae', '#584563', '#3e2137', '#9a6348', '#d79b7d', '#f5edba', '#c0c741', '#647d34', '#e4943a', '#9d303b', '#d26471', '#70377f', '#7ec4c1', '#34859d', '#17434b', '#1f0e1c')
# load data
taxa_BP_16s<- read.csv("16s/ALL/taxa_barplot_level2.csv", header = TRUE, sep = ",")
# Get number of reads per taxon
sort(colSums(taxa_BP_16s[2:37]), decreasing = TRUE)
# Subset dataframe to sample x taxon reads dataframe
taxa_sumarized_16s <- taxa_BP_16s[2:37]
# Sum the read count per sample of the taxa that are in low abundance
taxa_other <- (taxa_BP_16s$D_0__Bacteria.D_1__Cyanobacteria +
               taxa_BP_16s$D_0__Bacteria.D_1__Tenericutes +
               taxa_BP_16s$D_0__Bacteria.D_1__GAL15 +
               taxa_BP_16s$D_0__Bacteria.D_1__Zixibacteria +
               taxa_BP_16s$D_0__Bacteria.D_1__FCPU426 +
               taxa_BP_16s$D_0__Bacteria.D_1__Hydrogenedentes +
               taxa_BP_16s$D_0__Bacteria.D_1__WS2 +
               taxa_BP_16s$D_0__Bacteria.D_1__FBP +
               taxa_BP_16s$D_0__Bacteria.D_1__Chlamydiae +
               taxa_BP_16s$D_0__Bacteria.D_1__Entotheonellaeota +
               taxa_BP_16s$D_0__Bacteria.D_1__Omnitrophicaeota +
               taxa_BP_16s$D_0__Bacteria.D_1__Spirochaetes +
               taxa_BP_16s$D_0__Bacteria.D_1__WPS.2 +
               taxa_BP_16s$D_0__Bacteria.D_1__BRC1 +
               taxa_BP_16s$D_0__Bacteria.D_1__Dependentiae +
               taxa_BP_16s$D_0__Bacteria.D_1__Fibrobacteres +
               taxa_BP_16s$D_0__Bacteria.D_1__Elusimicrobia +
               taxa_BP_16s$D_0__Bacteria.D_1__Patescibacteria +
               taxa_BP_16s$D_0__Bacteria.D_1__Dadabacteria)
taxa_other[taxa_other<0] <- 0
# Taxa to remove from dataframe (either summed or were Archaea)
taxa_other_list <- c('D_0__Bacteria.D_1__Cyanobacteria',
                     'D_0__Bacteria.D_1__Tenericutes',
                     'D_0__Bacteria.D_1__GAL15',
                     'D_0__Bacteria.D_1__Zixibacteria',
                     'D_0__Bacteria.D_1__FCPU426',
                     'D_0__Archaea.D_1__Nanoarchaeaeota',
                     'D_0__Bacteria.D_1__Hydrogenedentes',
                     'D_0__Bacteria.D_1__WS2',
                     'D_0__Bacteria.D_1__FBP',
                     'D_0__Archaea.D_1__Euryarchaeota',
                     'D_0__Bacteria.D_1__Chlamydiae',
                     'D_0__Bacteria.D_1__Entotheonellaeota',
                     'D_0__Bacteria.D_1__Omnitrophicaeota',
                     'D_0__Bacteria.D_1__Spirochaetes',
                     'D_0__Bacteria.D_1__WPS.2',
                     'D_0__Bacteria.D_1__BRC1',
                     'D_0__Bacteria.D_1__Dependentiae',
                     'D_0__Bacteria.D_1__Fibrobacteres',
                     'D_0__Bacteria.D_1__Elusimicrobia',
                     'D_0__Bacteria.D_1__Patescibacteria',
                     'D_0__Bacteria.D_1__Dadabacteria',
                     'D_0__Archaea.D_1__Thaumarchaeota')
# Data handling to tidy data and make it into a format useful in ggplot2.
A<- cbind(Tissue=(taxa_BP_16s$Tissue), taxa_sumarized_16s, taxa_other)
B<- A[, -which(names(A) %in% taxa_other_list)]
C<- B %>% group_by(Tissue) %>% summarise_each(funs(sum))
D<- cbind(id = C[, 1], C[, -1]/rowSums(C[, -1]))
BP_16s<- D %>% gather(variable, value, -Tissue)
# For obtaining a vector of taxa names to make labels for the ggplot below.
X<- (C[,2:16])
X<- X[order(colSums(X), decreasing = TRUE)]
rev(colnames(X))

# Barplot 16s by tissue
taxa_barplot16s <- ggplot(BP_16s, aes(Tissue, y=value, fill = reorder(variable, value))) +
  geom_bar(stat = "identity", color = "black") +
  ylab("Relative abundance") +
  xlab("Tissue") +
  theme(legend.position = "right", legend.key.size = unit(0.75, "cm"), legend.text=element_text(size=12), legend.title=element_text(size=24)) +
  scale_fill_manual(name ="Bacterial phylum", values=taxapallete, labels = c(expression(italic("Deinococcus-Thermus")), expression(italic("Firmicutes")), expression(italic("Latescibacteria")), expression(italic("Armatimonadetes")), expression(italic("Nitrospirae")), expression(italic("Rokubacteria")), expression(italic("Gemmatimonadetes")), "Low abundance taxa", expression(italic("Chloroflexi")), expression(italic("Actinobacteria")), expression(italic("Planctomycetes")), expression(italic("Verrucomicrobia")), expression(italic("Bacteroidetes")), expression(italic("Acidobacteria")), expression(italic("Proteobacteria")))) +
  theme(legend.text.align = 0)

##### 2.2) 16S by tissue and rootstock ####
A<- cbind(Rootstock=(taxa_BP_16s$Rootstock), Tissue=(taxa_BP_16s$Tissue), taxa_sumarized_16s, taxa_other)
B<- A[, -which(names(A) %in% taxa_other_list)]
C<- B %>% group_by(Tissue, Rootstock) %>% summarise_each(funs(sum))

#1103P#
R_1103P <- C[c(1,5,9,13),]
R_1103P <- R_1103P[,-2]
R_1103P <- ungroup(R_1103P)
D <- cbind(id = R_1103P[, 1], R_1103P[, -1]/rowSums(R_1103P[, -1]))
R_1103P_BP <- D %>% gather(variable, value, -Tissue)
#PLOT
taxa_barplot16s_1103P <- PLOT_BP_16S_ROOTSTOCK(R_1103P_BP)

#3309C#
R_3309C <- C[c(2,6,10,14),]
R_3309C <- R_3309C[,-2]
R_3309C <- ungroup(R_3309C)
D <- cbind(id = R_3309C[, 1], R_3309C[, -1]/rowSums(R_3309C[, -1]))
R_3309C_BP <- D %>% gather(variable, value, -Tissue)
#PLOT
taxa_barplot16s_3309C <- PLOT_BP_16S_ROOTSTOCK(R_3309C_BP)

#S04#
R_S04 <- C[c(4,8,12,16),]
R_S04 <- R_S04[,-2]
R_S04 <- ungroup(R_S04)
D <- cbind(id = R_S04[, 1], R_S04[, -1]/rowSums(R_S04[, -1]))
R_S04_BP <- D %>% gather(variable, value, -Tissue)
#PLOT
taxa_barplot16s_S04 <- PLOT_BP_16S_ROOTSTOCK(R_S04_BP)

#Ungrafted#
R_Ungrafted <- C[c(3,7,11,15),]
R_Ungrafted <- R_Ungrafted[,-2]
R_Ungrafted <- ungroup(R_Ungrafted)
D <- cbind(id = R_Ungrafted[, 1], R_Ungrafted[, -1]/rowSums(R_Ungrafted[, -1]))
R_Ungrafted_BP <- D %>% gather(variable, value, -Tissue)
#PLOT
taxa_barplot16s_Ungrafted <- PLOT_BP_16S_ROOTSTOCK(R_Ungrafted_BP)


##### 2.3) ITS #####
# Load data
taxa_BP_its<- read.csv("ITS/ALL/taxa_barplot_level3.csv", header = TRUE, sep = ",")
# Get number of reads per taxon
sort(colSums(taxa_BP_its[2:31]), decreasing = TRUE)
# Subset dataframe to sample x taxon reads dataframe
taxa_sumarized_its <- taxa_BP_its[2:31]
# Sum the read count per sample of the taxa that are in low abundance
taxa_other_its <- (taxa_BP_its$k__Fungi.p__Mortierellomycota.__ +
                   taxa_BP_its$k__Fungi.p__Ascomycota.__ +
                   taxa_BP_its$k__Fungi.p__Basidiomycota.c__Pucciniomycetes +
                   taxa_BP_its$k__Fungi.p__unidentified.c__unidentified +
                   taxa_BP_its$k__Fungi.p__Ascomycota.c__Eurotiomycetes +
                   taxa_BP_its$k__Fungi.p__Glomeromycota.__ +
                   taxa_BP_its$k__Fungi.p__Basidiomycota.__ +
                   taxa_BP_its$k__Fungi.p__Rozellomycota.c__unidentified +
                   taxa_BP_its$k__Fungi.p__Olpidiomycota.c__Olpidiomycetes +
                   taxa_BP_its$k__Fungi.p__Rozellomycota.c__Rozellomycotina_cls_Incertae_sedis +
                   taxa_BP_its$k__Fungi.p__Ascomycota.c__Archaeorhizomycetes +
                   taxa_BP_its$k__Fungi.p__Mucoromycota.c__Mucoromycetes +
                   taxa_BP_its$k__Fungi.p__Ascomycota.c__Orbiliomycetes +
                   taxa_BP_its$k__Fungi.p__Entorrhizomycota.c__Entorrhizomycetes +
                   taxa_BP_its$k__Fungi.p__Basidiomycota.c__unidentified +
                   taxa_BP_its$k__Fungi.p__Basidiobolomycota.c__Basidiobolomycetes +
                   taxa_BP_its$k__Fungi.p__Calcarisporiellomycota.c__Calcarisporiellomycetes +
                   taxa_BP_its$k__Fungi.p__Basidiomycota.c__Exobasidiomycetes)
taxa_other_its[taxa_other_its<0] <- 0
# Taxa to remove from dataframe (either summed or were Archaea)
taxa_other_list_its <- c('k__Fungi.p__Mortierellomycota.__',
                         'k__Fungi.p__Ascomycota.__',
                         'k__Fungi.p__Basidiomycota.c__Pucciniomycetes',
                         'k__Fungi.p__unidentified.c__unidentified',
                         'k__Fungi.p__Ascomycota.c__Eurotiomycetes',
                         'k__Fungi.p__Glomeromycota.__',
                         'k__Fungi.p__Basidiomycota.__',
                         'k__Fungi.p__Rozellomycota.c__unidentified',
                         'k__Fungi.p__Olpidiomycota.c__Olpidiomycetes',
                         'k__Fungi.p__Rozellomycota.c__Rozellomycotina_cls_Incertae_sedis',
                         'k__Fungi.p__Ascomycota.c__Archaeorhizomycetes',
                         'k__Fungi.p__Mucoromycota.c__Mucoromycetes',
                         'k__Fungi.p__Ascomycota.c__Orbiliomycetes',
                         'k__Fungi.p__Entorrhizomycota.c__Entorrhizomycetes',
                         'k__Fungi.p__Basidiomycota.c__unidentified',
                         'k__Fungi.p__Basidiobolomycota.c__Basidiobolomycetes',
                         'k__Fungi.p__Calcarisporiellomycota.c__Calcarisporiellomycetes',
                         'k__Fungi.p__Basidiomycota.c__Exobasidiomycetes')

# Data handling to tidy data and make it into a format useful in ggplot2.
A <- cbind(Tissue=(taxa_BP_its$Tissue), taxa_sumarized_its, taxa_other_its)
B <- A[, -which(names(A) %in% taxa_other_list_its)]
C <- B %>% group_by(Tissue) %>% summarise_each(funs(sum))
D <- cbind(id = C[, 1], C[, -1]/rowSums(C[, -1]))
BP_its <- D %>% gather(variable, value, -Tissue)
# For obtaining a vector of taxa names to make labels for the ggplot below.
X<- (C[,2:14])
X<- X[order(colSums(X), decreasing = TRUE)]

### ITS barplot ###
taxa_barplotits <- ggplot(BP_its, aes(Tissue, y=value, fill = reorder(variable, value))) +
  geom_bar(stat = "identity", color = "black") +
  ylab("Relative abundance") +
  xlab("Tissue") +
  theme(legend.position = "right", legend.key.size = unit(0.75, "cm"), legend.text=element_text(size=12), legend.title=element_text(size=24)) +
  scale_fill_manual(name ="Fungal class", values=taxapallete, labels = c(expression(italic("Paraglomeromycetes")), expression(italic("Pezizomycetes")), expression(italic("Glomeromycetes")), "Low abundance taxa", expression(italic("Cystobasidiomycetes")), expression(italic("Leotiomycetes")), expression(italic("Saccharomycetes")), expression(italic("Mortierellomycetes")), expression(italic("Microbotryomycetes")), expression(italic("Agaricomycetes")), expression(italic("Tremellomycetes")), expression(italic("Sordariomycetes")), expression(italic("Dothideomycetes")))) +
  theme(legend.text.align = 0)
                                                                         
##### 2.4) ITS by tissue and rootstock #####
A<- cbind(Rootstock=(taxa_BP_its$Rootstock), Tissue=(taxa_BP_its$Tissue), taxa_sumarized_its, taxa_other_its)
B<- A[, -which(names(A) %in% taxa_other_list_its)]
C<- B %>% group_by(Tissue, Rootstock) %>% summarise_each(funs(sum))

#1103P#
R_1103P_its <- C[c(1,5,9,13),]
R_1103P_its <- R_1103P_its[,-2]
R_1103P_its <- ungroup(R_1103P_its)
D <- cbind(id = R_1103P_its[, 1], R_1103P_its[, -1]/rowSums(R_1103P_its[, -1]))
R_1103P_its_BP <- D %>% gather(variable, value, -Tissue)
#Plot
taxa_barplotits_1103P <- PLOT_BP_ITS_ROOTSTOCK(R_1103P_its_BP)

#3309C#
R_3309C_its <- C[c(2,6,10,14),]
R_3309C_its <- R_3309C_its[,-2]
R_3309C_its <- ungroup(R_3309C_its)
D <- cbind(id = R_3309C_its[, 1], R_3309C_its[, -1]/rowSums(R_3309C_its[, -1]))
R_3309C_its_BP <- D %>% gather(variable, value, -Tissue)
#PLOT
taxa_barplotits_3309C <- PLOT_BP_ITS_ROOTSTOCK(R_3309C_its_BP)

#S04#
R_S04_its <- C[c(4,8,12,16),]
R_S04_its <- R_S04_its[,-2]
R_S04_its <- ungroup(R_S04_its)
D <- cbind(id = R_S04_its[, 1], R_S04_its[, -1]/rowSums(R_S04_its[, -1]))
R_S04_its_BP <- D %>% gather(variable, value, -Tissue)
#PLOT
taxa_barplotits_S04 <- PLOT_BP_ITS_ROOTSTOCK(R_S04_its_BP)

#Ungrafted#
R_Ungrafted_its <- C[c(3,7,11,15),]
R_Ungrafted_its <- R_Ungrafted_its[,-2]
R_Ungrafted_its <- ungroup(R_Ungrafted_its)
D <- cbind(id = R_Ungrafted_its[, 1], R_Ungrafted_its[, -1]/rowSums(R_Ungrafted_its[, -1]))
R_Ungrafted_its_BP <- D %>% gather(variable, value, -Tissue)
#PLOT
taxa_barplotits_Ungrafted <- PLOT_BP_ITS_ROOTSTOCK(R_Ungrafted_its_BP)

# Remove unneeded objects
rm("A" , "B" , "C" , "D" , "X")

##### 3.0) Differential abundance analysis  #####
##### 3.1) 16s #####
#DESEQ2
physeq_16s <- qza_to_phyloseq('16s/ALL/filtered-table-no-mitochondria-no-chloroplast.qza','16s/ALL/rooted-tree.qza','16s/ALL/taxonomy.qza','16s/ALL/16s_noMockorPosorNeg_metadata.tsv', tmp="C:/tmp")
# recode factors OWN to Ungrafted
physeq_16s@sam_data$Rootstock <- recode_factor(physeq_16s@sam_data$Rootstock , OWN = "Ungrafted")
# Fix issue with Silva taxonomy strings before conducting DESEQ2 object.
physeq_16s@tax_table[,1]  <- gsub('D_.__', '', physeq_16s@tax_table[,1])
tax_table(physeq_16s) <- as.matrix(separate(as.data.frame(physeq_16s@tax_table[,1]), col = Kingdom, sep = ";", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")))
# Remove taxa that are not found greater than 25 times in 15% of the samples
physeq_16s <- filter_taxa(physeq_16s, function(x) sum(x > 25) > (0.10*length(x)), TRUE)
# In order for DESEQ2 to run there can be no zeros in the OTU matrix, so I add 1 to each count
otu_table(physeq_16s) <- otu_table(physeq_16s) + 1
# Create Deseq2 object with study desgin formula XXX
# DS2_16s <- phyloseq_to_deseq2(physeq_16s, ~ Rootstock + Tissue + Irrigation + Block + Rootstock:Tissue)
DS2_16s <- phyloseq_to_deseq2(physeq_16s, ~ Tissue*Rootstock + Irrigation + Irrigation:Tissue + Irrigation:Rootstock + Block)
DS2_16s <- DESeq(DS2_16s, test="Wald", fitType = "local") # Local fitype captures dispersion trend better than parametric.

# Results table contrasts
resultsNames(DS2_16s)
# Contrast list
Contrast_list <- c("Tissue_Leaf_vs_Berry", "Tissue_Root_vs_Berry", "Tissue_Soil_vs_Berry", "Rootstock_1103P_vs_Ungrafted", "Rootstock_3309C_vs_Ungrafted", "Rootstock_SO4_vs_Ungrafted", "Irrigation_none_vs_full", "Irrigation_rdi_vs_full", "Block_B_vs_A", "Block_C_vs_A", "TissueLeaf.Rootstock1103P", "TissueRoot.Rootstock1103P", "TissueSoil.Rootstock1103P", "TissueLeaf.Rootstock3309C", "TissueRoot.Rootstock3309C", "TissueSoil.Rootstock3309C", "TissueLeaf.RootstockSO4", "TissueRoot.RootstockSO4", "TissueSoil.RootstockSO4", "TissueLeaf.Irrigationnone", "TissueRoot.Irrigationnone", "TissueSoil.Irrigationnone", "TissueLeaf.Irrigationrdi", "TissueRoot.Irrigationrdi", "TissueSoil.Irrigationrdi", "Rootstock1103P.Irrigationnone", "Rootstock3309C.Irrigationnone", "RootstockSO4.Irrigationnone", "Rootstock1103P.Irrigationrdi", "Rootstock3309C.Irrigationrdi", "RootstockSO4.Irrigationrdi")
# Rough plotting of the number of significantly different abundant ASVs per factor
p <- c()
for (i in Contrast_list){
   p[[i]] <- as.data.frame(results(DS2_16s, cooksCutoff = FALSE, name = i))
   p[[i]] <- p[[i]][which(p[[i]]$padj < 0.05), ]
}

# Main effects Tissue, Rootstock, Irrigation, Block
T_num <- length(unique(c(rownames(p[["Tissue_Leaf_vs_Berry"]]), rownames(p[["Tissue_Root_vs_Berry"]]), rownames(p[["Tissue_Soil_vs_Berry"]]))))
R_num <- length(unique(c(rownames(p[["Rootstock_1103P_vs_Ungrafted"]]), rownames(p[["Rootstock_3309C_vs_Ungrafted"]]), rownames(p[["Rootstock_SO4_vs_Ungrafted"]]))))
I_num <- length(unique(c(rownames(p[["Irrigation_none_vs_full"]]), rownames(p[["Irrigation_rdi_vs_full"]]))))
B_num <- length(unique(c(rownames(p[["Block_B_vs_A"]]), rownames(p[["Block_C_vs_A"]]))))
# Interactions
RxT_num <- length(unique(c(rownames(p[["TissueLeaf.Rootstock1103P"]]), rownames(p[["TissueRoot.Rootstock1103P"]]), rownames(p[["TissueSoil.Rootstock1103P"]]), rownames(p[["TissueLeaf.Rootstock3309C"]]), rownames(p[["TissueRoot.Rootstock3309C"]]), rownames(p[["TissueSoil.Rootstock3309C"]]), rownames(p[["TissueLeaf.RootstockSO4"]]), rownames(p[["TissueRoot.RootstockSO4"]]), rownames(p[["TissueSoil.RootstockSO4"]]))))
TxI_num <- length(unique(c(rownames(p[["TissueLeaf.Irrigationnone"]]), rownames(p[["TissueRoot.Irrigationnone"]]), rownames(p[["TissueSoil.Irrigationnone"]]), rownames(p[["TissueLeaf.Irrigationrdi"]]), rownames(p[["TissueRoot.Irrigationrdi"]]), rownames(p[["TissueSoil.Irrigationrdi"]]))))
RxI_num <- length(unique(c(rownames(p[["Rootstock1103P.Irrigationnone"]]), rownames(p[["Rootstock3309C.Irrigationnone"]]), rownames(p[["RootstockSO4.Irrigationnone"]]), rownames(p[["Rootstock1103P.Irrigationrdi"]]), rownames(p[["Rootstock3309C.Irrigationrdi"]]), rownames(p[["RootstockSO4.Irrigationrdi"]]))))

Number_of_DiffAbund_ASVs_16s <- as.data.frame(c(T_num, R_num, I_num, B_num, RxT_num, TxI_num, RxI_num))
rownames(Number_of_DiffAbund_ASVs_16s) <- c( "Tissue", "Rootstock", "Irrigation", "Block", "Rootstock x Tissue", "Tissue x Irrigation", "Rootstock x Irrigation")
colnames(Number_of_DiffAbund_ASVs_16s) <- "num_ASVs"
Number_of_DiffAbund_ASVs_16s <- cbind(Number_of_DiffAbund_ASVs_16s, Factor = c("Tissue", "Rootstock", "Irrigation", "Block", "Rootstock x Tissue", "Tissue x Irrigation", "Rootstock x Irrigation")) #XXX
Number_of_DiffAbund_ASVs_16s$relabun <- (Number_of_DiffAbund_ASVs_16s$num_ASVs / 757)
TEMP_name <- ggplot(Number_of_DiffAbund_ASVs_16s,aes(x=Factor,y=relabun)) + geom_bar(stat="identity",position="dodge",width=0.9) + ylab("Proportion of ASVs responding to each factor") + ylim(0,1) + xlab("Factor") + ggtitle("16s data")
ggsave(filename = "DESeq_16s_prelim.pdf", plot = TEMP_name, units = "in", height = 5, width = 8, path = "Figures")

# Rough plotting of log fold changes by factor
p <- c()
for (i in Contrast_list){
  p[[i]] <- as.data.frame(results(DS2_16s, cooksCutoff = FALSE, name = i))
  p[[i]] <- p[[i]][which(p[[i]]$padj < 0.05), ]
}
p[[1]]

# Main effects
T_lf2c <- unique(rbind(p[[1]], p[[2]], p[[3]]))
R_lf2c <- unique(rbind(p[[4]], p[[5]], p[[6]]))
I_lf2c <- unique(rbind(p[[7]], p[[8]]))
B_lf2c <- unique(rbind(p[[9]], p[[10]]))
# Interactions
TxR_lf2c <- unique(rbind(p[[11]], p[[12]], p[[13]], p[[14]], p[[15]], p[[16]], p[[17]], p[[18]], p[[19]]))
TxI_lf2c <- unique(rbind(p[[20]], p[[21]], p[[22]], p[[23]], p[[24]], p[[25]]))
RxI_lf2c <- unique(rbind(p[[26]], p[[27]], p[[28]], p[[29]], p[[30]], p[[31]]))
# Extract Log2FoldChange column, Main Effects
R_lf2c <- as.data.frame(abs(R_lf2c$log2FoldChange))
T_lf2c <- as.data.frame(abs(T_lf2c$log2FoldChange))
B_lf2c <- as.data.frame(abs(B_lf2c$log2FoldChange))
I_lf2c <- as.data.frame(abs(I_lf2c$log2FoldChange))
# Extract Log2FoldChange column, Interactions
TxR_lf2c <- as.data.frame(abs(TxR_lf2c$log2FoldChange))
TxI_lf2c <- as.data.frame(abs(TxI_lf2c$log2FoldChange))
RxI_lf2c <- as.data.frame(abs(RxI_lf2c$log2FoldChange))

# Rename 1st column to Log2FoldChange
colnames(R_lf2c)[1] <- "Log2FoldChange"
colnames(T_lf2c)[1] <- "Log2FoldChange"
colnames(B_lf2c)[1] <- "Log2FoldChange"
colnames(I_lf2c)[1] <- "Log2FoldChange"
colnames(TxR_lf2c)[1] <- "Log2FoldChange"
colnames(TxI_lf2c)[1] <- "Log2FoldChange"
colnames(RxI_lf2c)[1] <- "Log2FoldChange"
# Add factor column
R_lf2c["Factor"] = "Rootstock"
T_lf2c["Factor"] = "Tissue"
B_lf2c["Factor"] = "Block"
I_lf2c["Factor"] = "Irrigation"
TxR_lf2c["Factor"] = "Rootstock x Tissue"
TxI_lf2c["Factor"] = "Tissue x Irrigation"
RxI_lf2c["Factor"] = "Rootstock x Irrigation"
# Bind by row all dfs
Log2fold_by_factor_16s <- rbind(R_lf2c, T_lf2c, B_lf2c, I_lf2c, TxR_lf2c, TxI_lf2c, RxI_lf2c)
# Box plot
ggplot(Log2fold_by_factor_16s, aes(x=Factor, y=Log2FoldChange)) + geom_boxplot(outlier.color = "red")
# Violin plot
ggplot(Log2fold_by_factor_16s, aes(x=Factor, y=Log2FoldChange)) + geom_violin(fill= "green", trim = FALSE, scale = "width") + stat_summary(fun.data=mean_sdl, geom="pointrange") + aes(x = fct_inorder(Factor)) + xlab("Factor")

# Creating dataframe for plotting of differentially abundant taxa by rootstock and tissue
otu_table(physeq_16s) <- otu_table(physeq_16s) - 1
otu_matrix <- as(otu_table(physeq_16s), "matrix")
if(taxa_are_rows(physeq_16s)){otu_matrix<- t(otu_matrix)}
otu_df <- as.data.frame(otu_matrix)


which(colnames(otu_df) == "7726eecd3a974455ae6f8c080caf8f5d")

plot_deseq2_DiffAbunMicob(physeq_16s, 214)

##### 3.2) ITS #####
#DESEQ2
physeq_its <- qza_to_phyloseq('ITS/ALL/filtered-nocontrol-trimmed-table.qza', taxonomy = 'ITS/ALL/taxonomy.qza', metadata = 'ITS/ALL/its_metadata_noControlsorWine.tsv', tmp="C:/tmp")
# recode factors OWN to Ungrafted
physeq_its@sam_data$Rootstock <- recode_factor(physeq_its@sam_data$Rootstock , OWN = "Ungrafted")
# Fix issue with UNITE taxonomy strings before conducting DESEQ2 on slices of the phyloseq object
physeq_its@tax_table[,1] <- gsub('.__', '', physeq_its@tax_table[,1])
tax_table(physeq_its) <- as.matrix(separate(as.data.frame(physeq_its@tax_table[,1]), col = Kingdom, sep = ";", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")))
# Remove taxa that are not found greater than 25 times in 10% of the samples
physeq_its <- filter_taxa(physeq_its, function(x) sum(x > 25) > (0.10*length(x)), TRUE)
physeq_its
# In order for DESEQ2 to run there can be no zeros in the OTU matrix, so I add 1 to each count
otu_table(physeq_its) <- otu_table(physeq_its) + 1

# Create Deseq2 object with study desgin formula
DS2_its <- phyloseq_to_deseq2(physeq_its, ~ Tissue*Rootstock + Irrigation + Irrigation:Tissue + Irrigation:Rootstock + Block)
DS2_its <- DESeq(DS2_its, test="Wald", fitType = "local")

# Results table
resultsNames(DS2_its)
# Contrast list
Contrast_list <- c("Tissue_Leaf_vs_Berry", "Tissue_Root_vs_Berry", "Tissue_Soil_vs_Berry", "Rootstock_1103P_vs_Ungrafted", "Rootstock_3309C_vs_Ungrafted", "Rootstock_SO4_vs_Ungrafted", "Irrigation_none_vs_full", "Irrigation_rdi_vs_full", "Block_B_vs_A", "Block_C_vs_A", "TissueLeaf.Rootstock1103P", "TissueRoot.Rootstock1103P", "TissueSoil.Rootstock1103P", "TissueLeaf.Rootstock3309C", "TissueRoot.Rootstock3309C", "TissueSoil.Rootstock3309C", "TissueLeaf.RootstockSO4", "TissueRoot.RootstockSO4", "TissueSoil.RootstockSO4", "TissueLeaf.Irrigationnone", "TissueRoot.Irrigationnone", "TissueSoil.Irrigationnone", "TissueLeaf.Irrigationrdi", "TissueRoot.Irrigationrdi", "TissueSoil.Irrigationrdi", "Rootstock1103P.Irrigationnone", "Rootstock3309C.Irrigationnone", "RootstockSO4.Irrigationnone", "Rootstock1103P.Irrigationrdi", "Rootstock3309C.Irrigationrdi", "RootstockSO4.Irrigationrdi")
# Rough plotting of the number of significantly different abundant ASVs per factor
p <- c()
for (i in Contrast_list){
  p[[i]] <- as.data.frame(results(DS2_its, cooksCutoff = FALSE, name = i))
  p[[i]] <- p[[i]][which(p[[i]]$padj < 0.05), ]
}

# Main effects Tissue, Rootstock, Irrigation, Block
T_num <- length(unique(c(rownames(p[["Tissue_Leaf_vs_Berry"]]), rownames(p[["Tissue_Root_vs_Berry"]]), rownames(p[["Tissue_Soil_vs_Berry"]]))))
R_num <- length(unique(c(rownames(p[["Rootstock_1103P_vs_Ungrafted"]]), rownames(p[["Rootstock_3309C_vs_Ungrafted"]]), rownames(p[["Rootstock_SO4_vs_Ungrafted"]]))))
I_num <- length(unique(c(rownames(p[["Irrigation_none_vs_full"]]), rownames(p[["Irrigation_rdi_vs_full"]]))))
B_num <- length(unique(c(rownames(p[["Block_B_vs_A"]]), rownames(p[["Block_C_vs_A"]]))))
# Interactions
RxT_num <- length(unique(c(rownames(p[["TissueLeaf.Rootstock1103P"]]), rownames(p[["TissueRoot.Rootstock1103P"]]), rownames(p[["TissueSoil.Rootstock1103P"]]), rownames(p[["TissueLeaf.Rootstock3309C"]]), rownames(p[["TissueRoot.Rootstock3309C"]]), rownames(p[["TissueSoil.Rootstock3309C"]]), rownames(p[["TissueLeaf.RootstockSO4"]]), rownames(p[["TissueRoot.RootstockSO4"]]), rownames(p[["TissueSoil.RootstockSO4"]]))))
TxI_num <- length(unique(c(rownames(p[["TissueLeaf.Irrigationnone"]]), rownames(p[["TissueRoot.Irrigationnone"]]), rownames(p[["TissueSoil.Irrigationnone"]]), rownames(p[["TissueLeaf.Irrigationrdi"]]), rownames(p[["TissueRoot.Irrigationrdi"]]), rownames(p[["TissueSoil.Irrigationrdi"]]))))
RxI_num <- length(unique(c(rownames(p[["Rootstock1103P.Irrigationnone"]]), rownames(p[["Rootstock3309C.Irrigationnone"]]), rownames(p[["RootstockSO4.Irrigationnone"]]), rownames(p[["Rootstock1103P.Irrigationrdi"]]), rownames(p[["Rootstock3309C.Irrigationrdi"]]), rownames(p[["RootstockSO4.Irrigationrdi"]]))))

Number_of_DiffAbund_ASVs_its <- as.data.frame(c(T_num, R_num, I_num, B_num, RxT_num, TxI_num, RxI_num))
rownames(Number_of_DiffAbund_ASVs_its) <- c( "Tissue", "Rootstock", "Irrigation", "Block", "Rootstock x Tissue", "Tissue x Irrigation", "Rootstock x Irrigation")
colnames(Number_of_DiffAbund_ASVs_its) <- "num_ASVs"
Number_of_DiffAbund_ASVs_its <- cbind(Number_of_DiffAbund_ASVs_its, Factor = c("Tissue", "Rootstock", "Irrigation", "Block", "Rootstock x Tissue", "Tissue x Irrigation", "Rootstock x Irrigation"))
Number_of_DiffAbund_ASVs_its$relabun <- (Number_of_DiffAbund_ASVs_its$num_ASVs / 111)
TEMP_name <- ggplot(Number_of_DiffAbund_ASVs_its,aes(x=Factor,y=relabun)) + geom_bar(stat="identity",position="dodge",width=0.9) + ylab("Proportion of ASVs responding to each factor") + ylim(0,1) + xlab("Factor") + ggtitle("16s data")
ggsave(filename = "DESeq_its_prelim.pdf", plot = TEMP_name, units = "in", height = 5, width = 8, path = "Figures")

# Rough plotting of log fold changes by factor
p <- c()
for (i in Contrast_list){
  p[[i]] <- as.data.frame(results(DS2_its, cooksCutoff = FALSE, name = i))
  p[[i]] <- p[[i]][which(p[[i]]$padj < 0.05), ]
}

# Main effects
T_lf2c <- unique(rbind(p[[1]], p[[2]], p[[3]]))
R_lf2c <- unique(rbind(p[[4]], p[[5]], p[[6]]))
I_lf2c <- unique(rbind(p[[7]], p[[8]]))
B_lf2c <- unique(rbind(p[[9]], p[[10]]))
# Interactions
TxR_lf2c <- unique(rbind(p[[11]], p[[12]], p[[13]], p[[14]], p[[15]], p[[16]], p[[17]], p[[18]], p[[19]]))
TxI_lf2c <- unique(rbind(p[[20]], p[[21]], p[[22]], p[[23]], p[[24]], p[[25]]))
RxI_lf2c <- unique(rbind(p[[26]], p[[27]], p[[28]], p[[29]], p[[30]], p[[31]]))
# Extract Log2FoldChange column, Main Effects
R_lf2c <- as.data.frame(abs(R_lf2c$log2FoldChange))
T_lf2c <- as.data.frame(abs(T_lf2c$log2FoldChange))
B_lf2c <- as.data.frame(abs(B_lf2c$log2FoldChange))
I_lf2c <- as.data.frame(abs(I_lf2c$log2FoldChange))
# Extract Log2FoldChange column, Interactions
TxR_lf2c <- as.data.frame(abs(TxR_lf2c$log2FoldChange))
TxI_lf2c <- as.data.frame(abs(TxI_lf2c$log2FoldChange))
RxI_lf2c <- as.data.frame(abs(RxI_lf2c$log2FoldChange))

# Rename 1st column to Log2FoldChange
colnames(R_lf2c)[1] <- "Log2FoldChange"
colnames(T_lf2c)[1] <- "Log2FoldChange"
colnames(B_lf2c)[1] <- "Log2FoldChange"
colnames(I_lf2c)[1] <- "Log2FoldChange"
colnames(TxR_lf2c)[1] <- "Log2FoldChange"
colnames(TxI_lf2c)[1] <- "Log2FoldChange"
colnames(RxI_lf2c)[1] <- "Log2FoldChange"
# Add factor column
R_lf2c["Factor"] = "Rootstock"
T_lf2c["Factor"] = "Tissue"
B_lf2c["Factor"] = "Block"
I_lf2c["Factor"] = "Irrigation"
TxR_lf2c["Factor"] = "Rootstock x Tissue"
TxI_lf2c["Factor"] = "Tissue x Irrigation"
RxI_lf2c["Factor"] = "Rootstock x Irrigation"
# Bind by row all dfs
Log2fold_by_factor_its <- rbind(R_lf2c, T_lf2c, B_lf2c, I_lf2c, TxR_lf2c, TxI_lf2c, RxI_lf2c)
# Box plot
ggplot(Log2fold_by_factor_its, aes(x=Factor, y=Log2FoldChange)) + geom_boxplot(outlier.color = "red")
# Violin plot
ggplot(Log2fold_by_factor_its, aes(x=Factor, y=Log2FoldChange)) + geom_violin(fill= "green", trim = FALSE, scale = "width") + stat_summary(fun.data=mean_sdl, geom="pointrange") + aes(x = fct_inorder(Factor)) + xlab("Factor")

# Combination bar plot
Number_of_DiffAbund_ASVs_16s["Marker"] = "16S"
Number_of_DiffAbund_ASVs_its["Marker"] = "ITS"
Number_of_DiffAbund_ASVs_both <- rbind(Number_of_DiffAbund_ASVs_16s, Number_of_DiffAbund_ASVs_its)
Number_of_DiffAbund_ASVs_both$Factor <- as.factor(Number_of_DiffAbund_ASVs_both$Factor)
Number_of_DiffAbund_ASVs_both$Factor <- factor(Number_of_DiffAbund_ASVs_both$Factor, levels = c("Rootstock", "Tissue", "Irrigation", "Rootstock x Tissue", "Rootstock x Irrigation", "Tissue x Irrigation", "Block"))
Numb_ASVs_both_plot <- ggplot(Number_of_DiffAbund_ASVs_both, aes(x=Factor,y=relabun, fill = Marker)) + geom_bar(stat="identity",position="dodge",width=0.9, colour="black") + geom_text(aes(label=num_ASVs), position=position_dodge(width=0.9), vjust=-0.25) + ylab("Proportion of ASVs") + ylim(0,1) + scale_fill_manual(name = "Marker", values=c("gray36", "gray80")) + xlab("Source of variation")


# Combination violin plot
Log2fold_by_factor_16s["Marker"] = "16S"
Log2fold_by_factor_its["Marker"] = "ITS"
Log2fold_by_factor_both <- rbind(Log2fold_by_factor_16s, Log2fold_by_factor_its)
Log2fold_by_factor_both$Factor <- as.factor(Log2fold_by_factor_both$Factor)
Log2fold_by_factor_both$Factor <- factor(Log2fold_by_factor_both$Factor, levels = c("Rootstock", "Tissue", "Irrigation", "Rootstock x Tissue", "Rootstock x Irrigation", "Tissue x Irrigation", "Block"))
Log2fold_both_plot <- ggplot(Log2fold_by_factor_both, aes(x=Factor, y=Log2FoldChange, fill=Marker)) + geom_violin(trim = FALSE, scale = "width") + stat_summary(fun.data=mean_sdl, geom="pointrange", position = position_dodge(width = 0.9))+ xlab("Factor") + scale_fill_manual(name = "Marker", values=c("gray36", "gray80"))  + labs(y=expression(paste(Log[2]," fold change")), x=("Source of variation"))

Bar_plot_and_log2foldchange_plot <- ggarrange(Numb_ASVs_both_plot, Log2fold_both_plot, nrow = 2, common.legend = TRUE, legend = "right")
ggsave(filename = "Bar_plot_and_log2foldchange_plot.pdf", plot = Bar_plot_and_log2foldchange_plot, units = "in", height = 12, width = 16, path = "Figures")


Bar_plot_and_log2foldchange_plot

# Creating dataframe for plotting of differentially abundant taxa by rootstock and tissue
otu_table(physeq_16s) <- otu_table(physeq_16s) - 1
which(colnames(otu_df) == "b467bbf7f9d84b55f768d9100667bb5c")

##### 5.0) Combining plots and saving #####
#Alpha diversity plots#
#16s by Tissue and rootsock
alpha_16s <- ggarrange(fphd_byR_T_16s, even_byR_T_16s, obso_byR_T_16s, shan_byR_T_16s, labels = c("A","B","C","D"), ncol = 2 , nrow = 2, common.legend = TRUE, legend = "right")
ggsave(filename = "16s_alphadiv.pdf", plot = alpha_16s, units = "in", height = 12, width = 12, path = "Figures")
#its by Tissue and rootstock
alpha_its <- ggarrange(even_byR_T_its, obso_byR_T_its, shan_byR_T_its, labels = c("A","B","C"), common.legend = TRUE, legend = "right")
ggsave(filename = "its_alphadiv.pdf", plot = alpha_its, units = "in", height = 12, width = 12, path = "Figures")
#Combine by Tissue
alpha_16s_its <- ggarrange(fphd_16s, shan_16s, obso_16s, even_16s, shan_its, obso_its, even_its, labels = c("A","B","C","D","E","F","G","H","I"))
ggsave(filename = "16s-its_alphadiv.pdf", plot = alpha_16s_its, units = "in", height = 12, width = 12, path = "Figures")
#Breakaway diversity
brek_fig <- ggarrange(brek_16s, brek_its, labels = c("A", "B"), ncol = 2)
ggsave(filename = "Breakaway_div_both.pdf", plot = brek_fig, units = "in", height = 5, width = 8, path = "Figures")

#Taxonomic bar plots#
#16s and ITS by tissue
taxabp_16s_its <- ggarrange(taxa_barplot16s, taxa_barplotits, labels = c("A", "B"), ncol = 2)
ggsave(filename = "16s-its_taxonomic_barplot.pdf", plot = taxabp_16s_its, units = "in", height = 6, width = 12, path = "Figures")

#16s by tissue and rootstock
taxabp_16s <- ggarrange(taxa_barplot16s_S04, taxa_barplot16s_3309C, taxa_barplot16s_1103P, taxa_barplot16s_Ungrafted, labels = c("   S04", "3309C", "1103P", "  Ungrafted"), label.x = 0.15, common.legend = TRUE, legend = "right")
ggsave(filename = "16s_taxonomic_barplot.pdf", plot = taxabp_16s, units = "in", height = 12, width = 12, path = "Figures")

#ITS by tissue and rootstock
taxabp_its<- ggarrange(taxa_barplotits_S04, taxa_barplotits_3309C, taxa_barplotits_1103P, taxa_barplotits_Ungrafted, labels = c("   S04", "3309C", "1103P", "  Ungrafted"), label.x = 0.15, common.legend = TRUE, legend = "right")
ggsave(filename = "its_taxonomic_barplot.pdf", plot = taxabp_its, units = "in", height = 12, width = 12, path = "Figures")

#PCoAs#
#16S
#weighted unifrac
PCoA_16s_wunif <- ggarrange(PCoA_1_2_wunif_16s, PCoA_1_3_wunif_16s, labels = c("A", "B"), ncol = 2, common.legend = TRUE, legend = "right")
ggsave(filename = "16s_PCoA_weightedunifrac.pdf", plot = PCoA_16s_wunif, units = "in", height = 6, width = 12, path = "Figures")

#unweighted unifrac
PCoA_16s_uunif <- ggarrange(PCoA_1_2_uunif_16s, PCoA_1_3_uunif_16s, labels = c("A", "B"), ncol = 2, common.legend = TRUE, legend = "right")
ggsave(filename = "16s_PCoA_unweightedunifrac.pdf", plot = PCoA_16s_uunif, units = "in", height = 6, width = 12, path = "Figures")

#jaccard
PCoA_16s_jaccd <- ggarrange(PCoA_1_2_jaccd_16s, PCoA_1_3_jaccd_16s, labels = c("A", "B"), ncol = 2, common.legend = TRUE, legend = "right")
ggsave(filename = "16s_PCoA_jaccard.pdf", plot = PCoA_16s_jaccd, units = "in", height = 6, width = 12, path = "Figures")

#its
#jaccard
PCoA_its_jaccd <- ggarrange(PCoA_1_2_jacc_its, PCoA_1_3_jacc_its, labels = c("A", "B"), ncol = 2, common.legend = TRUE, legend = "right")
ggsave(filename = "its_PCoA_jaccard.pdf", plot = PCoA_its_jaccd, units = "in", height = 6, width = 12, path = "Figures")

#bray
PCoA_its_bray <- ggarrange(PCoA_1_2_bray_its, PCoA_1_3_bray_its, labels = c("A", "B"), ncol = 2, common.legend = TRUE, legend = "right")
ggsave(filename = "its_PCoA_bray.pdf", plot = PCoA_its_bray, units = "in", height = 6, width = 12, path = "Figures")

# Figure 2
PCoA_fig2 <- ggarrange(PCoA_1_2_uunif_16s, PCoA_1_2_bray_its, taxa_barplot16s, taxa_barplotits,  labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2, legend = "right")
ggsave(filename = "Figure2_unweightedUnifrac1and2_bray1and2.png", plot = PCoA_fig2, units = "in", height = 12, width = 12, path = "Figures")

# Code that might be useful in reviews

# Test between QIIME2 and Phyloseq rarification and alpha diversity calcalutions
# This code might be handy if a reviewer is concerned that alpha diversity measures were calcaulated on a QIIME2 rareified ASV table
# and beta diversity measures were calcautled on a Phyloseq ASV Table. The chart below shows that the difference between these is
# negligible, I choose to stick with QIIME2 as it supports a broader range of alpha diversity measures (i.e. faith's phylo)
# TEMP <- estimate_richness(rare_physeq_16s, split = TRUE, measures = NULL)
# TEMP <- cbind(TEMP, rare_physeq_16s@sam_data)
# ggplot(TEMP, aes(TissueType, Observed)) + geom_boxplot(outlier.shape = NA) + ggtitle("16S") + theme(plot.title = element_text(hjust = 0.5)) + geom_jitter(width = 0.25) + xlab ("Tissue") + ylab("Observed ASVs")

# Function to plot corncob result by otu without removing outliers
# PLOT_OTU_corncob<- function(otu_number) {
#  TEMP <- as.data.frame(physeq_its@otu_table[otu_number,])
#  TEMP_long <- t(TEMP)
#  TEMP_long_2 <- cbind(TEMP_long, physeq_its@sam_data)
#  ggplot(TEMP_long_2, aes(TissueType, TEMP_long_2[,1], fill= Rootstock)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette) + scale_y_continuous(name="Abundance") + xlab("Tissue") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22)) + ggtitle(otu_number)
# }

# # Make subset of data by tissuetype, this way we are testing for differences
# # that are the effect of rootstock without tissue type being a confunding factor.
# Berry_16s<- subset_samples(physeq_16s, TissueType=="Berry")
# Leaf_16s<- subset_samples(physeq_16s, TissueType=="Leaf")
# Soil_16s<- subset_samples(physeq_16s, TissueType=="Soil")
# Root_16s<- subset_samples(physeq_16s, TissueType=="Root")
# 
# # DESEQ2 berry 16s
# DS2_16s<- phyloseq_to_deseq2(Berry_16s, ~ Rootstock)
# DS2_16s<- DESeq(DS2_16s, test="Wald", fitType = "parametric")
# res = results(DS2_16s, cooksCutoff = FALSE)
# sigtab = res[which(res$padj < 0.01), ]
# sigtab_berry_16s = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_16s)[rownames(sigtab), ], "matrix"))
# sigtab_berry_16s
# 
# # DESEQ2 leaf 16s
# DS2_16s<- phyloseq_to_deseq2(Leaf_16s, ~ Rootstock)
# DS2_16s<- DESeq(DS2_16s, test="Wald", fitType = "parametric")
# res = results(DS2_16s, cooksCutoff = FALSE)
# sigtab = res[which(res$padj < 0.01), ]
# sigtab_leaf_16s = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_16s)[rownames(sigtab), ], "matrix"))
# sigtab_leaf_16s
# 
# # DESEQ2 soil 16s
# DS2_16s<- phyloseq_to_deseq2(Soil_16s, ~ Rootstock)
# DS2_16s<- DESeq(DS2_16s, test="Wald", fitType = "parametric")
# res = results(DS2_16s, cooksCutoff = FALSE)
# sigtab = res[which(res$padj < 0.01), ]
# sigtab_soil_16s = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_16s)[rownames(sigtab), ], "matrix"))
# # sigtab_soil_16s
# 
# # DESEQ2 root 16s
# DS2_16s<- phyloseq_to_deseq2(Root_16s, ~ Rootstock)
# DS2_16s<- DESeq(DS2_16s, test="Wald", fitType = "parametric")
# res = results(DS2_16s, cooksCutoff = FALSE)
# sigtab = res[which(res$padj < 0.01), ]
# sigtab_root_16s = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_16s)[rownames(sigtab), ], "matrix"))
# sigtab_root_16s
# 
# 
# # Creating dataframe for plotting of differentially abundant taxa by rootstock and tissue
# otu_table(physeq_16s) <- otu_table(physeq_16s) - 1
# otu_matrix <- as(otu_table(physeq_16s), "matrix")
# if(taxa_are_rows(physeq_16s)){otu_matrix<- t(otu_matrix)}
# otu_df <- as.data.frame(otu_matrix)
# 
# which(colnames(otu_df) == "e9e2b7d11dd4cef36a1d5f3677574931")
# 
# specific_otu_df <- data.frame(row.names = rownames(otu_df), "ADD_NAME" = otu_df[,"e9e2b7d11dd4cef36a1d5f3677574931"], "Rootstock" = physeq_16s@sam_data$Rootstock, "TissueType" = physeq_16s@sam_data$TissueType)
# s16_diff_abun_taxa_bp<- ggplot(specific_otu_df[which(specific_otu_df$TissueType != "Leaf" & specific_otu_df$TissueType != "Berry"),], aes(TissueType, ADD_NAME, fill= Rootstock)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette) + scale_y_continuous(name="Abundance", breaks = c(0,80,160,240,320,400), limits = c(0,400)) + xlab("Tissue") + ggtitle("ASV 452 (Order: Rhizobiales)") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22))
# 
# # good examples ROOT e4e8962b9fe62868270bb1c8b0429891   ebf7dc25ef54e67b73d2c9a53244a838   d9be11723b6d824ccd900c865fddbe41
# # good example Berry 7726eecd3a974455ae6f8c080caf8f5d
# # good example Soil e9e2b7d11dd4cef36a1d5f3677574931

# # Make subset of data by tissuetype, this way we are testing for differences
# # that are the effect of rootstock without tissue type being a cofunding factor.
# Berry_its<- subset_samples(physeq_its, TissueType=="Berry")
# Leaf_its<- subset_samples(physeq_its, TissueType=="Leaf")
# Soil_its<- subset_samples(physeq_its, TissueType=="Soil")
# Root_its<- subset_samples(physeq_its, TissueType=="Root")
# 
# # DESEQ2 berry its
# DS2_its<- phyloseq_to_deseq2(Berry_its, ~ Rootstock)
# DS2_its<- DESeq(DS2_its, test="Wald", fitType = "parametric")
# res = results(DS2_its, cooksCutoff = FALSE)
# sigtab = res[which(res$padj < 0.05), ]
# sigtab_berry_its = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_its)[rownames(sigtab), ], "matrix"))
# sigtab_berry_its
# 
# # DESEQ2 leaf its
# DS2_its<- phyloseq_to_deseq2(Leaf_its, ~ Rootstock)
# DS2_its<- DESeq(DS2_its, test="Wald", fitType = "parametric")
# res = results(DS2_its, cooksCutoff = FALSE)
# sigtab = res[which(res$padj < 0.01), ]
# sigtab_leaf_its = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_its)[rownames(sigtab), ], "matrix"))
# sigtab_leaf_its
# 
# # DESEQ2 soil its
# DS2_its<- phyloseq_to_deseq2(Soil_its, ~ Rootstock)
# DS2_its<- DESeq(DS2_its, test="Wald", fitType = "parametric")
# res = results(DS2_its, cooksCutoff = FALSE)
# sigtab = res[which(res$padj < 0.01), ]
# sigtab_soil_its = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_its)[rownames(sigtab), ], "matrix"))
# sigtab_soil_its
# 
# # DESEQ2 root its
# DS2_its<- phyloseq_to_deseq2(Root_its, ~ Rootstock)
# DS2_its<- DESeq(DS2_its, test="Wald", fitType = "parametric")
# res = results(DS2_its, cooksCutoff = FALSE)
# sigtab = res[which(res$padj < 0.01), ]
# sigtab_root_its = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_its)[rownames(sigtab), ], "matrix"))
# sigtab_root_its
# 
# # Creating dataframe for plotting of differentially abundant taxa by rootstock and tissue
# otu_table(physeq_its) <- otu_table(physeq_its) - 1
# otu_matrix <- as(otu_table(physeq_its), "matrix")
# if(taxa_are_rows(physeq_its)){otu_matrix<- t(otu_matrix)}
# otu_df <- as.data.frame(otu_matrix)
# 
# which(colnames(otu_df) == "ab5e682bab43be7d10eb2082ef8ae281")
# 
# specific_otu_df <- data.frame(row.names = rownames(otu_df), "ADD_NAME" = otu_df[,"ab5e682bab43be7d10eb2082ef8ae281"], "Rootstock" = physeq_its@sam_data$Rootstock, "TissueType" = physeq_its@sam_data$TissueType)
# 
# its_diff_abun_taxa_bp <- ggplot(specific_otu_df[which(specific_otu_df$TissueType != "Leaf"),], aes(TissueType, ADD_NAME, fill= Rootstock)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette) + scale_y_continuous(name="Abundance", breaks = c(0,1200,2400,3600,4800,6000)) + xlab("Tissue") + ggtitle("ASV 58 (Order: Saccharomycetales)") + theme(legend.position="right", axis.title = element_text(size = 14), axis.text = element_text(size = 12), plot.title = element_text(size=22))
# 
# # Good examples berry b93271b97699bdbf400e0eb32b000a45    4a3f37e7e5c8d56bb929058f4efe1cfe    108362b3a9fcb289f305252e77a942f5   08cc4ea90d7e0b633f0e4df857c4ae7b
# # Good example soil ab5e682bab43be7d10eb2082ef8ae281
# 
# # Linear model on taxa that will be presented on the poster XXX
# Anova(lm(ADD_NAME ~ TissueType*Rootstock, data = specific_otu_df), type = "III")
# 
# # Save file for use on poster XXX
# ggsave(filename = "/Figures/its_diff_abun_taxa_bp.pdf", plot = its_diff_abun_taxa_bp, units = "in", height = 6, width = 6)

# ##### 4.0) Differential abundance analysis corncob  #####
# 
# #### 4.1) 16s #####
# # Clean import of phyloseq object
# physeq_16s <- qza_to_phyloseq('16s/ALL/filtered-table-no-mitochondria-no-chloroplast.qza','16s/ALL/rooted-tree.qza','16s/ALL/taxonomy.qza','16s/ALL/16s_noMockorPosorNeg_metadata.tsv', tmp="C:/tmp")
# 
# # recode factors OWN to Ungrafted
# physeq_16s@sam_data$Rootstock <- recode_factor(physeq_16s@sam_data$Rootstock , OWN = "Ungrafted")
# 
# # Fix issue with Silva taxonomy strings before conducting DESEQ2 on slices of the 
# # phyloseq 16s object.
# physeq_16s@tax_table[,1]  <- gsub('D_.__', '', physeq_16s@tax_table[,1])
# tax_table(physeq_16s) <- as.matrix(separate(as.data.frame(physeq_16s@tax_table[,1]), col = Kingdom, sep = ";", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")))
# 
# # Cleaning the names of ASVs to avoid issues with varible names (R does not like varibles that start with numbers)
# physeq_16s <- clean_taxa_names(physeq_16s)
# 
# corncob_null <- bbdml(formula = OTU793 ~ 1 , phi.formula = ~ 1, data = physeq_16s)
# corncob <- bbdml(formula = OTU793 ~ TissueType + Rootstock + IrrigationLevel, phi.formula = ~ 1, data = physeq_16s)
# plot(corncob, total = TRUE, color = "Rootstock")
# summary(corncob)
# 
# lrtest(mod_null = corncob_null, mod = corncob)
# 
# physeq_16s@sam_data
# set.seed(1)
# da_analysis <- differentialTest(formula = ~ Rootstock + Block + TissueType,
#                                 phi.formula = ~ 1,
#                                 formula_null = ~ Rootstock + Block + TissueType,
#                                 phi.formula_null = ~ 1,
#                                 test = "Wald", boot = FALSE,
#                                 data = physeq_16s,
#                                 fdr_cutoff = 0.05)
# 
# da_analysis$significant_taxa
# otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = physeq_16s)
# da_analysis$p_fdr
# plot(da_analysis)
# which(is.na(da_analysis$p)) %>% names
# 
# 
# ##### 4.2) ITS #####
# # clean import of phyloseq ITS object
# physeq_its <- qza_to_phyloseq('ITS/ALL/filtered-nocontrol-trimmed-table.qza', taxonomy = 'ITS/ALL/taxonomy.qza', metadata = 'ITS/ALL/its_metadata_noControlsorWine.tsv', tmp="C:/tmp")
# #LINUX physeq_its <- qza_to_phyloseq('ITS/ALL/filtered-nocontrol-trimmed-table.qza', taxonomy = 'ITS/ALL/taxonomy.qza', metadata = 'ITS/ALL/its_metadata_noControlsorWine.tsv')
# 
# # recode factors OWN to Ungrafted
# physeq_its@sam_data$Rootstock <- recode_factor(physeq_its@sam_data$Rootstock , OWN = "Ungrafted")
# # Fix issue with UNITE taxonomy strings before conducting DESEQ2 on slices of the 
# # phyloseq 16s object.
# physeq_its@tax_table[,1] <- gsub('.__', '', physeq_its@tax_table[,1])
# tax_table(physeq_its) <- as.matrix(separate(as.data.frame(physeq_its@tax_table[,1]), col = Kingdom, sep = ";", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")))
# 
# # Cleaning the names of ASVs to avoid issues with varible names (R does not like varibles that start with numbers)
# physeq_its <- clean_taxa_names(physeq_its)
# 
# corncob_null <- bbdml(formula = OTU307 ~ 1 , phi.formula = ~ 1, data = physeq_its)
# corncob <- bbdml(formula = OTU307 ~ Rootstock + Block + TissueType, phi.formula = ~ 1, data = physeq_its)
# summary(corncob)
# plot(corncob, total = TRUE, color = "Rootstock")
# 
# 
# 
# 
# PLOT_OTU_corncob_anomalized(896, grp_by = "IrrigationLevel")
# 
# da_analysis <- differentialTest(formula = ~ Rootstock + Block + TissueType,
#                                 phi.formula = ~ Rootstock + Block + TissueType,
#                                 formula_null = ~ 1,
#                                 phi.formula_null = ~ Rootstock + Block + TissueType,
#                                 test = "Wald", boot = FALSE,
#                                 data = physeq_its,
#                                 fdr_cutoff = 0.05)
# 
# Y <- da_analysis$significant_taxa
# as.data.frame(Y)
# 
# name_list<- otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = physeq_its)
# da_analysis$p_fdr
# plot(da_analysis, level = "Class")
# which(is.na(da_analysis$p)) %>% names
# 
# name_list
# # For posthoc test of a specific ASV found using DESEQ2
# otu_matrix <- as(otu_table(physeq_16s), "matrix")
# if(taxa_are_rows(physeq_16s)){otu_matrix<- t(otu_matrix)}
# otu_df <- as.data.frame(otu_matrix)
# which(colnames(otu_df) == "7726eecd3a974455ae6f8c080caf8f5d")
# specific_otu_df <- data.frame(row.names = rownames(otu_df), "OTU_NAME" = otu_df[,"7726eecd3a974455ae6f8c080caf8f5d"], "Rootstock" = physeq_16s@sam_data$Rootstock, "Tissue" = physeq_16s@sam_data$Tissue, "Irrigation" = physeq_16s@sam_data$Irrigation, "Block" = physeq_16s@sam_data$Block)
# # For posthoc test of a specific ASV found using DESEQ2
# # OTU_abundance_model <- lm(OTU_NAME ~ Tissue*Rootstock*Irrigation + Block, data = specific_otu_df)
# # Anova(lm(OTU_NAME ~ Tissue*Rootstock*Irrigation + Block, data = specific_otu_df), type = "III")
# # R16s.emm <- as.data.frame(pairs(emmeans(OTU_abundance_model, ~Tissue:Rootstock)))
# # R16s.emm <- transform(R16s.emm, comps=reshape2::colsplit(contrast, pattern = "-", names = c('c1', 'c2')))
# # R16s.emm <- transform(R16s.emm, c1=reshape2::colsplit(comps.c1, pattern = ",", names = c('i1', 'p1')))
# # R16s.emm <- transform(R16s.emm, c2=reshape2::colsplit(comps.c2, pattern = ",", names = c('i2', 'p2')))
# # R16s.emm$c2.p2 <- trimws(R16s.emm$c2.p2, which='left')
# # R16s.emm$c1.p1 <- trimws(R16s.emm$c1.p1, which='right')
# # R16s.emm <- R16s.emm[R16s.emm$c1.p1 == R16s.emm$c2.p2,] %>% select(contrast:p.value)
# # R16s.emm$padj <- p.adjust(R16s.emm$p.value, method='fdr')
# # R16s.emm[R16s.emm$padj < 0.05,]