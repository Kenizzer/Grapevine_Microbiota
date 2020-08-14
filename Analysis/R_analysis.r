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
library("patchwork")

# Theme set and Color Palettes
theme_set(theme_pubr(base_size = 16))
zoe_palette <- c("gray","#1b9e77", "#7570b3",  "#e6ab02")
compartment_pallete <- c("#5a1991", "#139d08", "#f5b784", "#5c3c0d") #https://lospec.com/palette-list/famicube

# Functions

# Calculate alpha diversity metrics and put them into a dataframe with the sample metadata
Alpha_div_metrics <- function(phyloseq_obj, its_or_16s) {
  if (its_or_16s == "16s"){
    # For 16s will also calculate faith's phylogenetic diversity metric using Picante package
    Alpha_16s <- estimate_richness(phyloseq_obj)
    R <- picante::pd(samp = as(t(otu_table(rare_physeq_16s)), "matrix"), tree = phy_tree(phyloseq_obj), include.root = FALSE)
    Alpha_16s <- cbind(Alpha_16s, "Faithpd" = R$PD, as.data.frame(phyloseq_obj@sam_data))
  } else if(its_or_16s == "its"){
    # Calculate alpha diversity metrics and put them into a dataframe with the sample metadata
    Alpha_its <- estimate_richness(phyloseq_obj)
    Alpha_its <- cbind(Alpha_its, as.data.frame(phyloseq_obj@sam_data))
  }
}

# Create taxonomic barplot for ITS
PLOT_BP_ITS_ROOTSTOCK <- function(plot_data) {
  ggplot(plot_data, aes(Compartment, y= value, fill = variable))+
    geom_bar(stat = "identity", color = "black") +
    ylab("Relative abundance") +
    theme(legend.position = "right", legend.key.size = unit(0.75, "cm"), legend.text=element_text(size=12), legend.title=element_text(size=24)) +
    scale_fill_manual(name ="Fungal class", values=taxapallete_ITS, labels = c(expression(italic("Agaricomycetes")), expression(italic("Cystobasidiomycetes")), expression(italic("Dothideomycetes")), expression(italic("Glomeromycetes")), expression(italic("Leotiomycetes")), expression(italic("Microbotryomycetes")), expression(italic("Mortierellomycetes")), expression(italic("Paraglomeromycetes")), expression(italic("Pezizomycetes")), expression(italic("Saccharomycetes")), expression(italic("Sordariomycetes")), expression(italic("Tremellomycetes")), "Low abundance taxa")) +
    theme(legend.text.align = 0)
}

# Create taxonomic barplot for 16S
PLOT_BP_16S_ROOTSTOCK <- function(plot_data) {
  ggplot(plot_data, aes(Compartment, y= value, fill = variable))+
    geom_bar(stat = "identity", color = "black") +
    ylab("Relative abundance") +
    theme(legend.position = "right", legend.key.size = unit(0.75, "cm"), legend.text=element_text(size=12), legend.title=element_text(size=24)) +
    scale_fill_manual(name ="Bacterial phylum", values=taxapallete_16S, labels = c(expression(italic("Acidobacteria")), expression(italic("Actinobacteria")), expression(italic("Armatimonadetes")), expression(italic("Bacteroidetes")), expression(italic("Chloroflexi")), expression(italic("Deinococcus-Thermus")), expression(italic("Firmicutes")), expression(italic("Gemmatimonadetes")), expression(italic("Latescibacteria")), expression(italic("Nitrospirae")), expression(italic("Planctomycetes")), expression(italic("Proteobacteria")), expression(italic("Rokubacteria")), expression(italic("Verrucomicrobia")), "Low abundance taxa")) +
    theme(legend.text.align = 0)
}

# Plot PCoA using ggplot from phyloseq ordination, can color points by rootstock if desired
PLOT_PCoA <- function(plot_data, distance_matrix, axis1, axis2, split_by_rootstock = TRUE) {
  temp <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2))
  if (split_by_rootstock == TRUE){
    plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2)) +
      geom_point(aes(fill=Rootstock, shape=Compartment), size = 3) +
      scale_shape_manual(values=c(24, 22, 23, 21)) +  
      scale_fill_manual(name = "Rootstock", values=zoe_palette) +
      labs(shape= "Compartment", color= "Rootstock") +
      xlab(paste("PCoA", axis1, sub(".*\\ ", "", temp$labels$x))) +
      ylab(paste("PCoA", axis2, sub(".*\\ ", "", temp$labels$y))) +
      guides(fill = guide_legend(override.aes = list(shape = 21)), shape = guide_legend(override.aes = list(fill = "black")))
  } else if(split_by_rootstock == FALSE){
    temp2 <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2), justDF = TRUE)
    ggplot(temp2, aes(x=Axis.1, y= Axis.2)) + geom_point(aes(fill = Compartment), size = 10, color="black",pch=21, alpha=0.8) +
      scale_fill_manual(name = "Compartment", values=compartment_pallete) +
      #scale_shape_manual(values=c(15, 1, 17, 5)) +  
      labs(shape= "Compartment") +
      xlab(paste("PCoA", axis1, sub(".*\\ ", "", temp$labels$x))) +
      ylab(paste("PCoA", axis2, sub(".*\\ ", "", temp$labels$y)))
  }
}

# Given a phyloseq object, return a dataframe of the sample metadata 
# From: https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# Set WD (Change as needed)
setwd("/PATH/TO/WorkingDirectory")

##### 1.0) Diversity ANALYSES #####
##### 1.1) Alpha Diversity #####

# Load QIIME2 objects into a phyloseq class object.
physeq_16s <- qza_to_phyloseq(features = '16s/ALL/filtered-table-no-mitochondria-no-chloroplast.qza', tree = '16s/ALL/rooted-tree.qza', taxonomy = '16s/ALL/taxonomy.qza', metadata = '16s/ALL/16s_noMockorPosorNeg_metadata.tsv', tmp="C:/tmp")
physeq_its <- qza_to_phyloseq(features = 'ITS/ALL/filtered-nocontrol-trimmed-table.qza', taxonomy = 'ITS/ALL/taxonomy.qza', metadata = 'ITS/ALL/its_metadata_noControlsorWine.tsv', tmp="C:/tmp")

# Recode factors: OWN to Ungrafted and correct names of irrigation treatments.
physeq_16s@sam_data$Rootstock <- recode_factor(physeq_16s@sam_data$Rootstock, OWN = "Ungrafted")
physeq_its@sam_data$Rootstock <- recode_factor(physeq_its@sam_data$Rootstock, OWN = "Ungrafted")
physeq_16s@sam_data$Irrigation <- recode_factor(physeq_16s@sam_data$Irrigation, none = "None", rdi = "RDI", full = "Full")
physeq_its@sam_data$Irrigation <- recode_factor(physeq_its@sam_data$Irrigation, none = "None", rdi = "RDI", full = "Full")

# Rarify bacterial (16s) samples to 1500.
rare_physeq_16s <- rarefy_even_depth(physeq_16s, sample.size = 1500, replace = FALSE, rngseed = 10031993) 
### Removing 7 biological samples

# Rarify fungal (ITS) samples to 5000.
rare_physeq_its <- rarefy_even_depth(physeq_its, sample.size = 5000, replace = FALSE, rngseed = 10031993) 
### Removing 7 biological samples

# Calculate alpha diversity metrics for 16s and ITS.
ALPHA_div_16s_rare <- Alpha_div_metrics(rare_physeq_16s, its_or_16s = "16s")
ALPHA_div_its_rare <- Alpha_div_metrics(rare_physeq_its, its_or_16s = "its")

# Means by Compartment
as.data.frame(ALPHA_div_16s_rare %>% group_by(Compartment) %>% summarise(Observed = mean(Observed), Shannon = mean(Shannon), invSimpson = mean(InvSimpson)))
as.data.frame(ALPHA_div_its_rare %>% group_by(Compartment) %>% summarise(Observed = mean(Observed), Shannon = mean(Shannon), invSimpson = mean(InvSimpson)))

# This tibble is used to place the significance letters at a consist height over each barplot
labels_df <- tibble(Compartment=levels(ALPHA_div_16s_rare$Compartment), Mfphd=max(ALPHA_div_16s_rare$Faithpd) * 1.2, Mobso=max(ALPHA_div_16s_rare$Observed) * 1.2, Mshan=max(ALPHA_div_16s_rare$Shannon) * 1.2, Msimpinv = max(ALPHA_div_16s_rare$InvSimpson) * 1.2, Mobso_i=max(ALPHA_div_its_rare$Observed) * 1.2, Mshan_i=max(ALPHA_div_its_rare$Shannon) * 1.2, Msimpinv_i = max(ALPHA_div_its_rare$InvSimpson) * 1.2)

# Make plots 16S by Compartment
fphd_16s <- ggplot(ALPHA_div_16s_rare, aes(Compartment, Faithpd)) + geom_boxplot(outlier.shape = NA) + theme(plot.title = element_text(hjust = 0.5)) + geom_jitter(width = 0.25) + xlab ("Compartment") + ylab("Faith's phylogenetic diversity") + geom_text(data=labels_df, aes(Compartment, Mfphd, label=c("a","a","b","c")), size = 6)
obso_16s <- ggplot(ALPHA_div_16s_rare, aes(Compartment, Observed)) + geom_boxplot(outlier.shape = NA) + theme(plot.title = element_text(hjust = 0.5)) + geom_jitter(width = 0.25) + xlab ("Compartment") + ylab("Observed ASVs") + geom_text(data=labels_df, aes(Compartment, Mobso, label=c("a","a","b","c")), size = 6)
simpI_16s <- ggplot(ALPHA_div_16s_rare, aes(Compartment, InvSimpson)) + geom_boxplot(outlier.shape = NA) + theme(plot.title = element_text(hjust = 0.5)) + geom_jitter(width = 0.25) + xlab ("Compartment") + ylab(expression("Simpson's D"^-1)) + geom_text(data=labels_df, aes(Compartment, Msimpinv, label=c("a","a","b","c")), size = 6)
shan_16s <- ggplot(ALPHA_div_16s_rare, aes(Compartment, Shannon)) + geom_boxplot(outlier.shape = NA) + theme(plot.title = element_text(hjust = 0.5)) + geom_jitter(width = 0.25) + xlab ("Compartment") + ylab("Shannon's index") + geom_text(data=labels_df, aes(Compartment, Mshan, label=c("a","a","b","c")), size = 6)

# Tukey tests 16s for significance letters means in alpha diversity charts
TukeyHSD(aov(Faithpd ~ Compartment, data = ALPHA_div_16s_rare), conf.level = 0.95)
TukeyHSD(aov(Observed ~ Compartment, data = ALPHA_div_16s_rare), conf.level = 0.95)
TukeyHSD(aov(InvSimpson ~ Compartment, data = ALPHA_div_16s_rare), conf.level = 0.95)
TukeyHSD(aov(Shannon ~ Compartment, data = ALPHA_div_16s_rare), conf.level = 0.95)

#Make plots ITS by Compartment
obso_its <- ggplot(ALPHA_div_its_rare, aes(Compartment, Observed)) + geom_boxplot(outlier.shape = NA) + theme(plot.title = element_text(hjust = 0.5)) + geom_jitter(width = 0.25) + xlab ("Compartment") + ylab("Observed ASVs") + geom_text(data=labels_df, aes(Compartment, Mobso_i, label=c("a","b","c","d")), size = 6)
simpI_its <- ggplot(ALPHA_div_its_rare, aes(Compartment, InvSimpson)) + geom_boxplot(outlier.shape = NA) + theme(plot.title = element_text(hjust = 0.5)) + geom_jitter(width = 0.25) + xlab ("Compartment") + ylab(expression("Simpson's D"^-1)) + geom_text(data=labels_df, aes(Compartment, Msimpinv_i, label=c("a","b","a","c")), size = 6)
shan_its <- ggplot(ALPHA_div_its_rare, aes(Compartment, Shannon)) + geom_boxplot(outlier.shape = NA) + theme(plot.title = element_text(hjust = 0.5)) + geom_jitter(width = 0.25) + xlab ("Compartment") + ylab("Shannon's index") + geom_text(data=labels_df, aes(Compartment, Mshan_i, label=c("a","b","c","d")), size = 6)

# Tukey tests ITS for significance letters means in alpha diversity charts
TukeyHSD(aov(Observed ~ Compartment, data = ALPHA_div_its_rare), conf.level = 0.95)
TukeyHSD(aov(InvSimpson ~ Compartment, data = ALPHA_div_its_rare), conf.level = 0.95)
TukeyHSD(aov(Shannon ~ Compartment, data = ALPHA_div_its_rare), conf.level = 0.95)

# Make plots 16S by Rootstock and Compartment
fphd_byR_T_16s <- ggplot(ALPHA_div_16s_rare, aes(Compartment, Faithpd, fill = Rootstock)) + geom_boxplot(outlier.shape = NA) + theme(plot.title = element_text(hjust = 0.5)) + xlab ("Compartment") + ylab("Faith's phylogenetic diversity") + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette)
simpI_byR_T_16s <- ggplot(ALPHA_div_16s_rare, aes(Compartment, InvSimpson, fill = Rootstock)) + geom_boxplot(outlier.shape = NA) + theme(plot.title = element_text(hjust = 0.5)) + xlab ("Compartment") + ylab(expression("Simpson's D"^-1)) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette)
obso_byR_T_16s <- ggplot(ALPHA_div_16s_rare, aes(Compartment, Observed, fill = Rootstock)) + geom_boxplot(outlier.shape = NA) + theme(plot.title = element_text(hjust = 0.5)) + xlab ("Compartment") + ylab("Observed ASVs") + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette)
shan_byR_T_16s <- ggplot(ALPHA_div_16s_rare, aes(Compartment, Shannon, fill = Rootstock)) + geom_boxplot(outlier.shape = NA) + theme(plot.title = element_text(hjust = 0.5)) + xlab ("Compartment") + ylab("Shannon's index") + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette)

# Make plots ITS by Rootstock and Compartment
simpI_byR_T_its <- ggplot(ALPHA_div_its_rare, aes(Compartment, InvSimpson, fill = Rootstock)) + geom_boxplot(outlier.shape = NA) + theme(plot.title = element_text(hjust = 0.5)) + xlab ("Compartment") + ylab(expression("Simpson's D"^-1)) + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette)
obso_byR_T_its <- ggplot(ALPHA_div_its_rare, aes(Compartment, Observed, fill = Rootstock)) + geom_boxplot(outlier.shape = NA) + theme(plot.title = element_text(hjust = 0.5)) + xlab ("Compartment") + ylab("Observed ASVs") + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette)
shan_byR_T_its <- ggplot(ALPHA_div_its_rare, aes(Compartment, Shannon, fill = Rootstock)) + geom_boxplot(outlier.shape = NA) + theme(plot.title = element_text(hjust = 0.5)) + xlab ("Compartment") + ylab("Shannon's index") + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette)

##### 1.2) Linear models #####

# Set contrasts to deviation coding for factors
options(contrasts = c("contr.sum", "contr.sum"))

# Assign levels to be used as reference
# For Rootstock = Ungrafted
# For Compartment = Soil
# For Irrigation = None
ALPHA_div_16s_rare$Rootstock <- relevel(ALPHA_div_16s_rare$Rootstock, ref = "Ungrafted")
ALPHA_div_16s_rare$Irrigation <- relevel(ALPHA_div_16s_rare$Irrigation, ref = "None")
ALPHA_div_16s_rare$Compartment <- relevel(ALPHA_div_16s_rare$Compartment, ref = "Soil")
ALPHA_div_its_rare$Rootstock <- relevel(ALPHA_div_its_rare$Rootstock, ref = "Ungrafted")
ALPHA_div_its_rare$Irrigation <- relevel(ALPHA_div_its_rare$Irrigation, ref = "None")
ALPHA_div_its_rare$Compartment <- relevel(ALPHA_div_its_rare$Compartment, ref = "Soil")

# Linear models 16S FULL
fphd_16s_fit_full <- lm(Faithpd ~ Compartment*Rootstock*Irrigation + Block, data = ALPHA_div_16s_rare)
simpI_16s_fit_full <- lm(InvSimpson ~ Compartment*Rootstock*Irrigation + Block, data = ALPHA_div_16s_rare)
obso_16s_fit_full <- lm(Observed ~ Compartment*Rootstock*Irrigation + Block, data = ALPHA_div_16s_rare)
shan_16s_fit_full <- lm(Shannon ~ Compartment*Rootstock*Irrigation + Block, data = ALPHA_div_16s_rare)

# Print results, copied to Excel
as.data.frame(Anova(fphd_16s_fit_full, type = "III"))
as.data.frame(Anova(simpI_16s_fit_full, type = "III"))
as.data.frame(Anova(obso_16s_fit_full, type = "III"))
as.data.frame(Anova(shan_16s_fit_full, type = "III"))

# Posthoc testing
pairs(emmeans(fphd_16s_fit_full, ~ Compartment))
pairs(emmeans(fphd_16s_fit_full, ~ Block))
pairs(emmeans(simpI_16s_fit_full, ~ Compartment))
pairs(emmeans(simpI_16s_fit_full, ~ Block))
pairs(emmeans(obso_16s_fit_full, ~ Compartment))
pairs(emmeans(shan_16s_fit_full, ~ Compartment))

# Linear models ITS FULL
simpI_its_fit_full <- lm(InvSimpson ~ Compartment*Rootstock*Irrigation + Block, data = ALPHA_div_its_rare)
obso_its_fit_full <- lm(Observed ~ Compartment*Rootstock*Irrigation + Block, data = ALPHA_div_its_rare)
shan_its_fit_full <- lm(Shannon ~ Compartment*Rootstock*Irrigation + Block, data = ALPHA_div_its_rare)

# Print results, copied to Excel
as.data.frame(Anova(simpI_its_fit_full, type = "III"))
as.data.frame(Anova(obso_its_fit_full, type = "III"))
as.data.frame(Anova(shan_its_fit_full, type = "III"))

# Posthoc testing
pairs(emmeans(simpI_its_fit_full, ~ Compartment))
pairs(emmeans(simpI_its_fit_full, ~ Rootstock))
pairs(emmeans(simpI_its_fit_full, ~ Rootstock|Compartment))
pairs(emmeans(obso_its_fit_full, ~ Compartment))
pairs(emmeans(shan_its_fit_full, ~ Compartment))

# Linear models 16S optimized
fphd_16s_fit <- lm(Faithpd ~ Compartment + Rootstock + Irrigation + Block, data = ALPHA_div_16s_rare)
simpI_16s_fit <- lm(InvSimpson ~ Compartment + Rootstock + Irrigation + Block, data = ALPHA_div_16s_rare)
obso_16s_fit <- lm(Observed ~ Compartment + Rootstock + Irrigation + Block, data = ALPHA_div_16s_rare)
shan_16s_fit <- lm(Shannon ~ Compartment + Rootstock + Irrigation + Block, data = ALPHA_div_16s_rare)

# Print results, copied to Excel
as.data.frame(Anova(fphd_16s_fit, type = "II"))
as.data.frame(Anova(simpI_16s_fit, type = "II"))
as.data.frame(Anova(obso_16s_fit, type = "II"))
as.data.frame(Anova(shan_16s_fit, type = "II"))

# Posthoc testing
pairs(emmeans(fphd_16s_fit, ~ Compartment))
pairs(emmeans(fphd_16s_fit, ~ Block))
pairs(emmeans(simpI_16s_fit, ~ Compartment))
pairs(emmeans(simpI_16s_fit, ~ Block))
pairs(emmeans(obso_16s_fit, ~ Compartment))
pairs(emmeans(shan_16s_fit, ~ Compartment))

# Linear models ITS optimized
simpI_its_fit <- lm(InvSimpson ~ Compartment + Rootstock + Irrigation + Block, data = ALPHA_div_its_rare)
obso_its_fit <- lm(Observed ~ Compartment + Rootstock + Irrigation + Block, data = ALPHA_div_its_rare)
shan_its_fit <- lm(Shannon ~ Compartment + Rootstock + Irrigation + Block, data = ALPHA_div_its_rare)

# Print results, copied to Excel
as.data.frame(Anova(simpI_its_fit, type = "II"))
as.data.frame(Anova(obso_its_fit, type = "II"))
as.data.frame(Anova(shan_its_fit, type = "II"))

# Posthoc testing
pairs(emmeans(simpI_its_fit, ~ Compartment))
pairs(emmeans(simpI_its_fit, ~ Rootstock))
pairs(emmeans(simpI_its_fit, ~ Block))
pairs(emmeans(obso_its_fit, ~ Compartment))
pairs(emmeans(shan_its_fit, ~ Compartment))

# Reset contrasts
options(contrasts = c("contr.treatment", "contr.poly"))

##### 1.3) Beta Diversity #####

# 16s beta diversity calculations
# Calculate jaccard, weighted, and unweighted unifrac distances
out.wunifrac <- ordinate(rare_physeq_16s, method = "MDS", distance = "wunifrac")
out.uunifrac <- ordinate(rare_physeq_16s, method = "MDS", distance = "uunifrac")
out.jaccard  <- ordinate(rare_physeq_16s, method = "MDS", distance = "jaccard")

# Calculate distance matrixes for use in Vegan
out.dist.wunifrac <- phyloseq::distance(rare_physeq_16s, "wunifrac")
out.dist.uunifrac <- phyloseq::distance(rare_physeq_16s, "uunifrac")
out.dist.jaccard  <- phyloseq::distance(rare_physeq_16s, "jaccard", binary = TRUE)

# get variance explained by axes 1-3
sum(out.wunifrac$values$Eigenvalues[1:3])/sum(out.wunifrac$values$Eigenvalues) #71%
sum(out.uunifrac$values$Eigenvalues[1:3])/sum(out.uunifrac$values$Eigenvalues) #41%
sum(out.jaccard$values$Eigenvalues[1:3])/sum(out.jaccard$values$Eigenvalues) #39%

# Plot PCAs 1x2 and 1x3 for Jaccard, weighted, and unweighted unifrac
PCoA_1_2_wunif_16s <- PLOT_PCoA(physeq_16s, out.wunifrac, 1, 2)
PCoA_1_3_wunif_16s <- PLOT_PCoA(physeq_16s, out.wunifrac, 1, 3)
PCoA_1_2_uunif_16s <- PLOT_PCoA(physeq_16s, out.uunifrac, 1, 2)
PCoA_1_3_uunif_16s <- PLOT_PCoA(physeq_16s, out.uunifrac, 1, 3)
PCoA_1_2_jaccd_16s <- PLOT_PCoA(physeq_16s, out.jaccard, 1, 2)
PCoA_1_3_jaccd_16s <- PLOT_PCoA(physeq_16s, out.jaccard, 1, 3)

# Using this one in Figure 2, colored by compartment and not split by rootstock
PCoA_1_2_uunif_16s_c <- PLOT_PCoA(physeq_16s, out.uunifrac, 1, 2, split_by_rootstock = FALSE)

# Create sample data dataframe from phyloseq object to use in vegan/adonis
physeq_metadata_16s <- pssd2veg(rare_physeq_16s)

# Adonis testing on weighted and unweighted unifrac distance matrix generated using phyloseq::distance
# Full models
adonis2(out.dist.wunifrac ~ Compartment*Rootstock*Irrigation + Block, data = physeq_metadata_16s, permutations = 10000)
adonis2(out.dist.uunifrac ~ Compartment*Rootstock*Irrigation + Block, data = physeq_metadata_16s, permutations = 10000)
## Optimized model
adonis2(out.dist.wunifrac ~ Rootstock + Compartment + Irrigation + Block, data = physeq_metadata_16s, permutations = 10000)
adonis2(out.dist.uunifrac ~ Rootstock + Compartment + Irrigation + Block, data = physeq_metadata_16s, permutations = 10000)

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
# Using this one in Figure 2, colored by compartment and not split by rootstock
PCoA_1_2_bray_its_c <- PLOT_PCoA(physeq_its, out.bray_its, 1, 2, split_by_rootstock = FALSE)
  
# Create sample data dataframe from phyloseq object to use in vegan/adonis
physeq_metadata_its <- pssd2veg(rare_physeq_its)

# Adonis testing on weighted and unweighted unifrac distance matrix generated using
# Full models
adonis2(out.dist.bray_its ~ Compartment*Rootstock*Irrigation + Block, data = physeq_metadata_its, permutations = 10000)
adonis2(out.dist.jaccard_its ~ Compartment*Rootstock*Irrigation + Block, data = physeq_metadata_its, permutations = 10000)
## Optimized model (stepwise)
adonis2(out.dist.bray_its ~ Rootstock + Compartment + Irrigation + Block + Rootstock:Compartment, data = physeq_metadata_its, permutations = 10000)
adonis2(out.dist.jaccard_its ~ Rootstock + Compartment + Irrigation + Block, data = physeq_metadata_its, permutations = 10000)

##### 2.0) Taxonomic Barplots #####

##### 2.1) 16s #####

# Color palette for ggplot2 taxonomic barplots
taxapallete_16S <- c('#8c8fae', '#584563', '#3e2137', '#9a6348', '#d79b7d', '#f5edba', '#c0c741', '#647d34', '#e4943a', '#9d303b', '#d26471', '#582d63', '#7ec4c1', '#34859d', '#17434b', '#1f0e1c', '#213b25', '#c0d1cc', '#4da6ff', '#ffffff') # NA16 Palette, slightly modified (https://lospec.com/palette-list/na16)
taxapallete_ITS <- c('#ffffff', '#6df7c1', '#11adc1', '#606c81', '#393457', '#1e8875', '#5bb361', '#a1e55a', '#f7e476', '#f99252', '#cb4d68', '#6a3771', '#c92464', '#f48cb6', '#f7b69e', '#9b9c82') # Island joy 16 (https://lospec.com/palette-list/island-joy-16)

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
A <- cbind(Compartment=(taxa_BP_16s$Compartment), taxa_sumarized_16s, taxa_other)
B <- A[, -which(names(A) %in% taxa_other_list)]
C <- B %>% group_by(Compartment) %>% summarise_each(funs(sum))
D <- cbind(id = C[, 1], C[, -1]/rowSums(C[, -1]))
BP_16s <- D %>% gather(variable, value, -Compartment)

# Reseting factor level of the taxa names to alphabetic, this will keep it so each plot has the same order/legend order.
BP_16s$variable <- factor(BP_16s$variable, levels = c("D_0__Bacteria.D_1__Acidobacteria", "D_0__Bacteria.D_1__Actinobacteria", "D_0__Bacteria.D_1__Armatimonadetes", "D_0__Bacteria.D_1__Bacteroidetes", "D_0__Bacteria.D_1__Chloroflexi", "D_0__Bacteria.D_1__Deinococcus.Thermus", "D_0__Bacteria.D_1__Firmicutes", "D_0__Bacteria.D_1__Gemmatimonadetes", "D_0__Bacteria.D_1__Latescibacteria", "D_0__Bacteria.D_1__Nitrospirae", "D_0__Bacteria.D_1__Planctomycetes", "D_0__Bacteria.D_1__Proteobacteria", "D_0__Bacteria.D_1__Rokubacteria", "D_0__Bacteria.D_1__Verrucomicrobia", "taxa_other"))

# Barplot 16s by Compartment
taxa_barplot16s <- ggplot(BP_16s, aes(Compartment, y=value, fill = variable)) +
  geom_bar(stat = "identity", color = "black") +
  ylab("Relative abundance") +
  xlab("Compartment") +
  theme(legend.position = "right", legend.key.size = unit(0.75, "cm"), legend.text=element_text(size=12), legend.title=element_text(size=24)) +
  scale_fill_manual(name ="Bacterial phylum", values=taxapallete_16S, labels = c(expression(italic("Acidobacteria")), expression(italic("Actinobacteria")), expression(italic("Armatimonadetes")), expression(italic("Bacteroidetes")), expression(italic("Chloroflexi")), expression(italic("Deinococcus-Thermus")), expression(italic("Firmicutes")), expression(italic("Gemmatimonadetes")), expression(italic("Latescibacteria")), expression(italic("Nitrospirae")), expression(italic("Planctomycetes")), expression(italic("Proteobacteria")), expression(italic("Rokubacteria")), expression(italic("Verrucomicrobia")), "Low abundance taxa")) +
  theme(legend.text.align = 0)

# 16S by Compartment and rootstock
A<- cbind(Rootstock=(taxa_BP_16s$Rootstock), Compartment=(taxa_BP_16s$Compartment), taxa_sumarized_16s, taxa_other)
B<- A[, -which(names(A) %in% taxa_other_list)]
C<- B %>% group_by(Compartment, Rootstock) %>% summarise_each(funs(sum))

#1103P#
R_1103P <- C[c(1,5,9,13),]
R_1103P <- R_1103P[,-2]
R_1103P <- ungroup(R_1103P)
D <- cbind(id = R_1103P[, 1], R_1103P[, -1]/rowSums(R_1103P[, -1]))
R_1103P_BP <- D %>% gather(variable, value, -Compartment)
R_1103P_BP$variable <- factor(R_1103P_BP$variable, levels = c("D_0__Bacteria.D_1__Acidobacteria", "D_0__Bacteria.D_1__Actinobacteria", "D_0__Bacteria.D_1__Armatimonadetes", "D_0__Bacteria.D_1__Bacteroidetes", "D_0__Bacteria.D_1__Chloroflexi", "D_0__Bacteria.D_1__Deinococcus.Thermus", "D_0__Bacteria.D_1__Firmicutes", "D_0__Bacteria.D_1__Gemmatimonadetes", "D_0__Bacteria.D_1__Latescibacteria", "D_0__Bacteria.D_1__Nitrospirae", "D_0__Bacteria.D_1__Planctomycetes", "D_0__Bacteria.D_1__Proteobacteria", "D_0__Bacteria.D_1__Rokubacteria", "D_0__Bacteria.D_1__Verrucomicrobia", "taxa_other"))
#PLOT
taxa_barplot16s_1103P <- PLOT_BP_16S_ROOTSTOCK(R_1103P_BP)

#3309C#
R_3309C <- C[c(2,6,10,14),]
R_3309C <- R_3309C[,-2]
R_3309C <- ungroup(R_3309C)
D <- cbind(id = R_3309C[, 1], R_3309C[, -1]/rowSums(R_3309C[, -1]))
R_3309C_BP <- D %>% gather(variable, value, -Compartment)
R_3309C_BP$variable <- factor(R_3309C_BP$variable, levels = c("D_0__Bacteria.D_1__Acidobacteria", "D_0__Bacteria.D_1__Actinobacteria", "D_0__Bacteria.D_1__Armatimonadetes", "D_0__Bacteria.D_1__Bacteroidetes", "D_0__Bacteria.D_1__Chloroflexi", "D_0__Bacteria.D_1__Deinococcus.Thermus", "D_0__Bacteria.D_1__Firmicutes", "D_0__Bacteria.D_1__Gemmatimonadetes", "D_0__Bacteria.D_1__Latescibacteria", "D_0__Bacteria.D_1__Nitrospirae", "D_0__Bacteria.D_1__Planctomycetes", "D_0__Bacteria.D_1__Proteobacteria", "D_0__Bacteria.D_1__Rokubacteria", "D_0__Bacteria.D_1__Verrucomicrobia", "taxa_other"))
#PLOT
taxa_barplot16s_3309C <- PLOT_BP_16S_ROOTSTOCK(R_3309C_BP)

#S04#
R_S04 <- C[c(4,8,12,16),]
R_S04 <- R_S04[,-2]
R_S04 <- ungroup(R_S04)
D <- cbind(id = R_S04[, 1], R_S04[, -1]/rowSums(R_S04[, -1]))
R_S04_BP <- D %>% gather(variable, value, -Compartment)
R_S04_BP$variable <- factor(R_S04_BP$variable, levels = c("D_0__Bacteria.D_1__Acidobacteria", "D_0__Bacteria.D_1__Actinobacteria", "D_0__Bacteria.D_1__Armatimonadetes", "D_0__Bacteria.D_1__Bacteroidetes", "D_0__Bacteria.D_1__Chloroflexi", "D_0__Bacteria.D_1__Deinococcus.Thermus", "D_0__Bacteria.D_1__Firmicutes", "D_0__Bacteria.D_1__Gemmatimonadetes", "D_0__Bacteria.D_1__Latescibacteria", "D_0__Bacteria.D_1__Nitrospirae", "D_0__Bacteria.D_1__Planctomycetes", "D_0__Bacteria.D_1__Proteobacteria", "D_0__Bacteria.D_1__Rokubacteria", "D_0__Bacteria.D_1__Verrucomicrobia", "taxa_other"))
#PLOT
taxa_barplot16s_S04 <- PLOT_BP_16S_ROOTSTOCK(R_S04_BP)

#Ungrafted#
R_Ungrafted <- C[c(3,7,11,15),]
R_Ungrafted <- R_Ungrafted[,-2]
R_Ungrafted <- ungroup(R_Ungrafted)
D <- cbind(id = R_Ungrafted[, 1], R_Ungrafted[, -1]/rowSums(R_Ungrafted[, -1]))
R_Ungrafted_BP <- D %>% gather(variable, value, -Compartment)
R_Ungrafted_BP$variable <- factor(R_Ungrafted_BP$variable, levels = c("D_0__Bacteria.D_1__Acidobacteria", "D_0__Bacteria.D_1__Actinobacteria", "D_0__Bacteria.D_1__Armatimonadetes", "D_0__Bacteria.D_1__Bacteroidetes", "D_0__Bacteria.D_1__Chloroflexi", "D_0__Bacteria.D_1__Deinococcus.Thermus", "D_0__Bacteria.D_1__Firmicutes", "D_0__Bacteria.D_1__Gemmatimonadetes", "D_0__Bacteria.D_1__Latescibacteria", "D_0__Bacteria.D_1__Nitrospirae", "D_0__Bacteria.D_1__Planctomycetes", "D_0__Bacteria.D_1__Proteobacteria", "D_0__Bacteria.D_1__Rokubacteria", "D_0__Bacteria.D_1__Verrucomicrobia", "taxa_other"))
#PLOT
taxa_barplot16s_Ungrafted <- PLOT_BP_16S_ROOTSTOCK(R_Ungrafted_BP)

##### 2.2) ITS #####

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

# Taxa to remove from dataframe (low abundance)
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
A <- cbind(Compartment=(taxa_BP_its$Compartment), taxa_sumarized_its, taxa_other_its)
B <- A[, -which(names(A) %in% taxa_other_list_its)]
C <- B %>% group_by(Compartment) %>% summarise_each(funs(sum))
D <- cbind(id = C[, 1], C[, -1]/rowSums(C[, -1]))
BP_its <- D %>% gather(variable, value, -Compartment)

# Reseting factor level of the taxa names to alphabetic, this will keep it so each plot has the same order/legend order.
BP_its$variable <- factor(BP_its$variable, levels = c("k__Fungi.p__Basidiomycota.c__Agaricomycetes", "k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes", "k__Fungi.p__Ascomycota.c__Dothideomycetes", "k__Fungi.p__Glomeromycota.c__Glomeromycetes", "k__Fungi.p__Ascomycota.c__Leotiomycetes", "k__Fungi.p__Basidiomycota.c__Microbotryomycetes", "k__Fungi.p__Mortierellomycota.c__Mortierellomycetes", "k__Fungi.p__Glomeromycota.c__Paraglomeromycetes", "k__Fungi.p__Ascomycota.c__Pezizomycetes", "k__Fungi.p__Ascomycota.c__Saccharomycetes", "k__Fungi.p__Ascomycota.c__Sordariomycetes", "k__Fungi.p__Basidiomycota.c__Tremellomycetes", "taxa_other_its"))

### ITS barplot ###
taxa_barplotits <- ggplot(BP_its, aes(Compartment, y=value, fill = variable)) +
  geom_bar(stat = "identity", color = "black") +
  ylab("Relative abundance") +
  xlab("Compartment") +
  theme(legend.position = "right", legend.key.size = unit(0.75, "cm"), legend.text=element_text(size=12), legend.title=element_text(size=24)) +
  scale_fill_manual(name ="Fungal class", values=taxapallete_ITS, labels = c(expression(italic("Agaricomycetes")), expression(italic("Cystobasidiomycetes")), expression(italic("Dothideomycetes")), expression(italic("Glomeromycetes")), expression(italic("Leotiomycetes")), expression(italic("Microbotryomycetes")), expression(italic("Mortierellomycetes")), expression(italic("Paraglomeromycetes")), expression(italic("Pezizomycetes")), expression(italic("Saccharomycetes")), expression(italic("Sordariomycetes")), expression(italic("Tremellomycetes")), "Low abundance taxa")) +
  theme(legend.text.align = 0)
                                                                         
# ITS by Compartment and rootstock
A <- cbind(Rootstock=(taxa_BP_its$Rootstock), Compartment=(taxa_BP_its$Compartment), taxa_sumarized_its, taxa_other_its)
B <- A[, -which(names(A) %in% taxa_other_list_its)]
C <- B %>% group_by(Compartment, Rootstock) %>% summarise_each(funs(sum))

#1103P#
R_1103P_its <- C[c(1,5,9,13),]
R_1103P_its <- R_1103P_its[,-2]
R_1103P_its <- ungroup(R_1103P_its)
D <- cbind(id = R_1103P_its[, 1], R_1103P_its[, -1]/rowSums(R_1103P_its[, -1]))
R_1103P_its_BP <- D %>% gather(variable, value, -Compartment)
R_1103P_its_BP$variable <- factor(R_1103P_its_BP$variable, levels = c("k__Fungi.p__Basidiomycota.c__Agaricomycetes", "k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes", "k__Fungi.p__Ascomycota.c__Dothideomycetes", "k__Fungi.p__Glomeromycota.c__Glomeromycetes", "k__Fungi.p__Ascomycota.c__Leotiomycetes", "k__Fungi.p__Basidiomycota.c__Microbotryomycetes", "k__Fungi.p__Mortierellomycota.c__Mortierellomycetes", "k__Fungi.p__Glomeromycota.c__Paraglomeromycetes", "k__Fungi.p__Ascomycota.c__Pezizomycetes", "k__Fungi.p__Ascomycota.c__Saccharomycetes", "k__Fungi.p__Ascomycota.c__Sordariomycetes", "k__Fungi.p__Basidiomycota.c__Tremellomycetes", "taxa_other_its"))
#Plot
taxa_barplotits_1103P <- PLOT_BP_ITS_ROOTSTOCK(R_1103P_its_BP)

#3309C#
R_3309C_its <- C[c(2,6,10,14),]
R_3309C_its <- R_3309C_its[,-2]
R_3309C_its <- ungroup(R_3309C_its)
D <- cbind(id = R_3309C_its[, 1], R_3309C_its[, -1]/rowSums(R_3309C_its[, -1]))
R_3309C_its_BP <- D %>% gather(variable, value, -Compartment)
R_3309C_its_BP$variable <- factor(R_3309C_its_BP$variable, levels = c("k__Fungi.p__Basidiomycota.c__Agaricomycetes", "k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes", "k__Fungi.p__Ascomycota.c__Dothideomycetes", "k__Fungi.p__Glomeromycota.c__Glomeromycetes", "k__Fungi.p__Ascomycota.c__Leotiomycetes", "k__Fungi.p__Basidiomycota.c__Microbotryomycetes", "k__Fungi.p__Mortierellomycota.c__Mortierellomycetes", "k__Fungi.p__Glomeromycota.c__Paraglomeromycetes", "k__Fungi.p__Ascomycota.c__Pezizomycetes", "k__Fungi.p__Ascomycota.c__Saccharomycetes", "k__Fungi.p__Ascomycota.c__Sordariomycetes", "k__Fungi.p__Basidiomycota.c__Tremellomycetes", "taxa_other_its"))
#PLOT
taxa_barplotits_3309C <- PLOT_BP_ITS_ROOTSTOCK(R_3309C_its_BP)

#S04#
R_S04_its <- C[c(4,8,12,16),]
R_S04_its <- R_S04_its[,-2]
R_S04_its <- ungroup(R_S04_its)
D <- cbind(id = R_S04_its[, 1], R_S04_its[, -1]/rowSums(R_S04_its[, -1]))
R_S04_its_BP <- D %>% gather(variable, value, -Compartment)
R_S04_its_BP$variable <- factor(R_S04_its_BP$variable, levels = c("k__Fungi.p__Basidiomycota.c__Agaricomycetes", "k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes", "k__Fungi.p__Ascomycota.c__Dothideomycetes", "k__Fungi.p__Glomeromycota.c__Glomeromycetes", "k__Fungi.p__Ascomycota.c__Leotiomycetes", "k__Fungi.p__Basidiomycota.c__Microbotryomycetes", "k__Fungi.p__Mortierellomycota.c__Mortierellomycetes", "k__Fungi.p__Glomeromycota.c__Paraglomeromycetes", "k__Fungi.p__Ascomycota.c__Pezizomycetes", "k__Fungi.p__Ascomycota.c__Saccharomycetes", "k__Fungi.p__Ascomycota.c__Sordariomycetes", "k__Fungi.p__Basidiomycota.c__Tremellomycetes", "taxa_other_its"))
#PLOT
taxa_barplotits_S04 <- PLOT_BP_ITS_ROOTSTOCK(R_S04_its_BP)

#Ungrafted#
R_Ungrafted_its <- C[c(3,7,11,15),]
R_Ungrafted_its <- R_Ungrafted_its[,-2]
R_Ungrafted_its <- ungroup(R_Ungrafted_its)
D <- cbind(id = R_Ungrafted_its[, 1], R_Ungrafted_its[, -1]/rowSums(R_Ungrafted_its[, -1]))
R_Ungrafted_its_BP <- D %>% gather(variable, value, -Compartment)
R_Ungrafted_its_BP$variable <- factor(R_Ungrafted_its_BP$variable, levels = c("k__Fungi.p__Basidiomycota.c__Agaricomycetes", "k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes", "k__Fungi.p__Ascomycota.c__Dothideomycetes", "k__Fungi.p__Glomeromycota.c__Glomeromycetes", "k__Fungi.p__Ascomycota.c__Leotiomycetes", "k__Fungi.p__Basidiomycota.c__Microbotryomycetes", "k__Fungi.p__Mortierellomycota.c__Mortierellomycetes", "k__Fungi.p__Glomeromycota.c__Paraglomeromycetes", "k__Fungi.p__Ascomycota.c__Pezizomycetes", "k__Fungi.p__Ascomycota.c__Saccharomycetes", "k__Fungi.p__Ascomycota.c__Sordariomycetes", "k__Fungi.p__Basidiomycota.c__Tremellomycetes", "taxa_other_its"))
#PLOT
taxa_barplotits_Ungrafted <- PLOT_BP_ITS_ROOTSTOCK(R_Ungrafted_its_BP)

# Remove unneeded objects
rm("A" , "B" , "C" , "D")

##### 3.0) Differential abundance analysis  #####
##### 3.1) 16s #####
#DESEQ2
physeq_16s <- qza_to_phyloseq('16s/ALL/filtered-table-no-mitochondria-no-chloroplast.qza','16s/ALL/rooted-tree.qza','16s/ALL/taxonomy.qza','16s/ALL/16s_noMockorPosorNeg_metadata.tsv', tmp="C:/tmp")

# recode factors OWN to Ungrafted
physeq_16s@sam_data$Rootstock <- recode_factor(physeq_16s@sam_data$Rootstock , OWN = "Ungrafted")

# Fix issue with Silva taxonomy strings before conducting DESEQ2 object.
physeq_16s@tax_table[,1]  <- gsub('D_.__', '', physeq_16s@tax_table[,1])
tax_table(physeq_16s) <- as.matrix(separate(as.data.frame(physeq_16s@tax_table[,1]), col = Kingdom, sep = ";", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")))

# Remove taxa that are not found greater than 25 times in 10% of the samples
physeq_16s <- filter_taxa(physeq_16s, function(x) sum(x > 25) > (0.10*length(x)), TRUE)

# In order for DESEQ2 to run there can be no zeros in the OTU matrix, so I add 1 to each count
otu_table(physeq_16s) <- otu_table(physeq_16s) + 1

# Create Deseq2 object with study desgin formula
DS2_16s <- phyloseq_to_deseq2(physeq_16s, ~ Compartment*Rootstock + Irrigation + Irrigation:Compartment + Irrigation:Rootstock + Block)
DS2_16s <- DESeq(DS2_16s, test="Wald", fitType = "local") # Local fitype captures dispersion trend better than parametric.

# Results table contrasts
resultsNames(DS2_16s)

# Contrast list
Contrast_list <- c("Compartment_Leaf_vs_Berry", "Compartment_Root_vs_Berry", "Compartment_Soil_vs_Berry", "Rootstock_1103P_vs_Ungrafted", "Rootstock_3309C_vs_Ungrafted", "Rootstock_SO4_vs_Ungrafted", "Irrigation_none_vs_full", "Irrigation_rdi_vs_full", "Block_B_vs_A", "Block_C_vs_A", "CompartmentLeaf.Rootstock1103P", "CompartmentRoot.Rootstock1103P", "CompartmentSoil.Rootstock1103P", "CompartmentLeaf.Rootstock3309C", "CompartmentRoot.Rootstock3309C", "CompartmentSoil.Rootstock3309C", "CompartmentLeaf.RootstockSO4", "CompartmentRoot.RootstockSO4", "CompartmentSoil.RootstockSO4", "CompartmentLeaf.Irrigationnone", "CompartmentRoot.Irrigationnone", "CompartmentSoil.Irrigationnone", "CompartmentLeaf.Irrigationrdi", "CompartmentRoot.Irrigationrdi", "CompartmentSoil.Irrigationrdi", "Rootstock1103P.Irrigationnone", "Rootstock3309C.Irrigationnone", "RootstockSO4.Irrigationnone", "Rootstock1103P.Irrigationrdi", "Rootstock3309C.Irrigationrdi", "RootstockSO4.Irrigationrdi")

# Extracting significantly different abundant ASVs per factor
p <- c()
for (i in Contrast_list){
   p[[i]] <- as.data.frame(results(DS2_16s, cooksCutoff = FALSE, name = i))
   p[[i]] <- p[[i]][which(p[[i]]$padj < 0.05), ]
}

# Main effects Compartment, Rootstock, Irrigation, Block
T_num <- length(unique(c(rownames(p[["Compartment_Leaf_vs_Berry"]]), rownames(p[["Compartment_Root_vs_Berry"]]), rownames(p[["Compartment_Soil_vs_Berry"]]))))
R_num <- length(unique(c(rownames(p[["Rootstock_1103P_vs_Ungrafted"]]), rownames(p[["Rootstock_3309C_vs_Ungrafted"]]), rownames(p[["Rootstock_SO4_vs_Ungrafted"]]))))
I_num <- length(unique(c(rownames(p[["Irrigation_none_vs_full"]]), rownames(p[["Irrigation_rdi_vs_full"]]))))
B_num <- length(unique(c(rownames(p[["Block_B_vs_A"]]), rownames(p[["Block_C_vs_A"]]))))

# Interactions
RxT_num <- length(unique(c(rownames(p[["CompartmentLeaf.Rootstock1103P"]]), rownames(p[["CompartmentRoot.Rootstock1103P"]]), rownames(p[["CompartmentSoil.Rootstock1103P"]]), rownames(p[["CompartmentLeaf.Rootstock3309C"]]), rownames(p[["CompartmentRoot.Rootstock3309C"]]), rownames(p[["CompartmentSoil.Rootstock3309C"]]), rownames(p[["CompartmentLeaf.RootstockSO4"]]), rownames(p[["CompartmentRoot.RootstockSO4"]]), rownames(p[["CompartmentSoil.RootstockSO4"]]))))
TxI_num <- length(unique(c(rownames(p[["CompartmentLeaf.Irrigationnone"]]), rownames(p[["CompartmentRoot.Irrigationnone"]]), rownames(p[["CompartmentSoil.Irrigationnone"]]), rownames(p[["CompartmentLeaf.Irrigationrdi"]]), rownames(p[["CompartmentRoot.Irrigationrdi"]]), rownames(p[["CompartmentSoil.Irrigationrdi"]]))))
RxI_num <- length(unique(c(rownames(p[["Rootstock1103P.Irrigationnone"]]), rownames(p[["Rootstock3309C.Irrigationnone"]]), rownames(p[["RootstockSO4.Irrigationnone"]]), rownames(p[["Rootstock1103P.Irrigationrdi"]]), rownames(p[["Rootstock3309C.Irrigationrdi"]]), rownames(p[["RootstockSO4.Irrigationrdi"]]))))

Number_of_DiffAbund_ASVs_16s <- as.data.frame(c(T_num, R_num, I_num, B_num, RxT_num, TxI_num, RxI_num))
rownames(Number_of_DiffAbund_ASVs_16s) <- c( "Compartment", "Rootstock", "Irrigation", "Block", "Rootstock x Compartment", "Compartment x Irrigation", "Rootstock x Irrigation")
colnames(Number_of_DiffAbund_ASVs_16s) <- "num_ASVs"
Number_of_DiffAbund_ASVs_16s <- cbind(Number_of_DiffAbund_ASVs_16s, Factor = c("Compartment", "Rootstock", "Irrigation", "Block", "Rootstock x Compartment", "Compartment x Irrigation", "Rootstock x Irrigation"))
Number_of_DiffAbund_ASVs_16s$relabun <- (Number_of_DiffAbund_ASVs_16s$num_ASVs / 757)# This 757 number comes from the number of unique taxa from filtering the phyloseq object above

# Get Taxonomy of differentailly abundant ASVs per factor
# Main effects Compartment, Rootstock, Irrigation, Block
T_ASV_list <- as.data.frame(unique(c(rownames(p[["Compartment_Leaf_vs_Berry"]]), rownames(p[["Compartment_Root_vs_Berry"]]), rownames(p[["Compartment_Soil_vs_Berry"]]))))
R_ASV_list <- as.data.frame(unique(c(rownames(p[["Rootstock_1103P_vs_Ungrafted"]]), rownames(p[["Rootstock_3309C_vs_Ungrafted"]]), rownames(p[["Rootstock_SO4_vs_Ungrafted"]]))))
I_ASV_list <- as.data.frame(unique(c(rownames(p[["Irrigation_none_vs_full"]]), rownames(p[["Irrigation_rdi_vs_full"]]))))
B_ASV_list <- as.data.frame(unique(c(rownames(p[["Block_B_vs_A"]]), rownames(p[["Block_C_vs_A"]]))))

# Interactions
RxT_ASV_list <- as.data.frame(unique(c(rownames(p[["CompartmentLeaf.Rootstock1103P"]]), rownames(p[["CompartmentRoot.Rootstock1103P"]]), rownames(p[["CompartmentSoil.Rootstock1103P"]]), rownames(p[["CompartmentLeaf.Rootstock3309C"]]), rownames(p[["CompartmentRoot.Rootstock3309C"]]), rownames(p[["CompartmentSoil.Rootstock3309C"]]), rownames(p[["CompartmentLeaf.RootstockSO4"]]), rownames(p[["CompartmentRoot.RootstockSO4"]]), rownames(p[["CompartmentSoil.RootstockSO4"]]))))
TxI_ASV_list <- as.data.frame(unique(c(rownames(p[["CompartmentLeaf.Irrigationnone"]]), rownames(p[["CompartmentRoot.Irrigationnone"]]), rownames(p[["CompartmentSoil.Irrigationnone"]]), rownames(p[["CompartmentLeaf.Irrigationrdi"]]), rownames(p[["CompartmentRoot.Irrigationrdi"]]), rownames(p[["CompartmentSoil.Irrigationrdi"]]))))
RxI_ASV_list <- as.data.frame(unique(c(rownames(p[["Rootstock1103P.Irrigationnone"]]), rownames(p[["Rootstock3309C.Irrigationnone"]]), rownames(p[["RootstockSO4.Irrigationnone"]]), rownames(p[["Rootstock1103P.Irrigationrdi"]]), rownames(p[["Rootstock3309C.Irrigationrdi"]]), rownames(p[["RootstockSO4.Irrigationrdi"]]))))

# Get dataframe of all ASV hash number for use in next step
x <- as.data.frame(tax_table(physeq_16s))
x$ASV_hash <- rownames(x)

# Main effects Compartment, Rootstock, Irrigation, Block
T_ASV_list <- x[x$ASV_hash %in% T_ASV_list$uniq, ]
R_ASV_list <- x[x$ASV_hash %in% R_ASV_list$uniq, ]
I_ASV_list <- x[x$ASV_hash %in% I_ASV_list$uniq, ]
B_ASV_list <- x[x$ASV_hash %in% B_ASV_list$uniq, ]

# Interaction
RxT_ASV_list <- x[x$ASV_hash %in% RxT_ASV_list$uniq, ]
TxI_ASV_list <- x[x$ASV_hash %in% TxI_ASV_list$uniq, ]
RxI_ASV_list <- x[x$ASV_hash %in% RxI_ASV_list$uniq, ]

# Extract taxonomy
Class_ASV_df <- as.data.frame(bind_rows(summary(T_ASV_list$Phylum), summary(R_ASV_list$Phylum), summary(I_ASV_list$Phylum), summary(B_ASV_list$Phylum), summary(RxT_ASV_list$Phylum), summary(TxI_ASV_list$Phylum), summary(RxI_ASV_list$Phylum)))
Class_ASV_df <- Class_ASV_df[, colSums(Class_ASV_df != 0) > 0]

# Merge relative proportion of ASVs df with taxonomy df
num_tax_ASVs_16s <- cbind(Number_of_DiffAbund_ASVs_16s, Class_ASV_df)

# Covert taxonomy counts into proportions
num_tax_ASVs_16s[1,c(seq(4,22,1))] <- num_tax_ASVs_16s[1,3] * (num_tax_ASVs_16s[1,c(seq(4,22,1))] / sum(num_tax_ASVs_16s[1,c(seq(4,22,1))]))
num_tax_ASVs_16s[2,c(seq(4,22,1))] <- num_tax_ASVs_16s[2,3] * (num_tax_ASVs_16s[2,c(seq(4,22,1))] / sum(num_tax_ASVs_16s[2,c(seq(4,22,1))]))
num_tax_ASVs_16s[3,c(seq(4,22,1))] <- num_tax_ASVs_16s[3,3] * (num_tax_ASVs_16s[3,c(seq(4,22,1))] / sum(num_tax_ASVs_16s[3,c(seq(4,22,1))]))
num_tax_ASVs_16s[4,c(seq(4,22,1))] <- num_tax_ASVs_16s[4,3] * (num_tax_ASVs_16s[4,c(seq(4,22,1))] / sum(num_tax_ASVs_16s[4,c(seq(4,22,1))]))
num_tax_ASVs_16s[5,c(seq(4,22,1))] <- num_tax_ASVs_16s[5,3] * (num_tax_ASVs_16s[5,c(seq(4,22,1))] / sum(num_tax_ASVs_16s[5,c(seq(4,22,1))]))
num_tax_ASVs_16s[6,c(seq(4,22,1))] <- num_tax_ASVs_16s[6,3] * (num_tax_ASVs_16s[6,c(seq(4,22,1))] / sum(num_tax_ASVs_16s[6,c(seq(4,22,1))]))
num_tax_ASVs_16s[7,c(seq(4,22,1))] <- num_tax_ASVs_16s[7,3] * (num_tax_ASVs_16s[7,c(seq(4,22,1))] / sum(num_tax_ASVs_16s[7,c(seq(4,22,1))]))

# Gather columns into a single column for use in ggplot
num_tax_ASVs_16s <- gather(num_tax_ASVs_16s, key = "Taxon", value = "Freq", -c(Factor,relabun, num_ASVs))
num_tax_ASVs_16s$Taxon <- factor(num_tax_ASVs_16s$Taxon)

# Step the factor order of the taxon names
num_tax_ASVs_16s$Taxon <- factor(num_tax_ASVs_16s$Taxon, levels =  c("Acidobacteria", "Actinobacteria", "Armatimonadetes" , "Bacteroidetes", "Chloroflexi", "Deinococcus-Thermus", "Firmicutes", "Gemmatimonadetes", "Latescibacteria", "Nitrospirae", "Planctomycetes", "Proteobacteria", "Rokubacteria", "Verrucomicrobia", "Cyanobacteria", "Elusimicrobia", "Entotheonellaeota", "Fibrobacteres", "Spirochaetes", "Thaumarchaeota", "NA's"))

# Plotting of log fold changes by factor
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
T_lf2c["Factor"] = "Compartment"
B_lf2c["Factor"] = "Block"
I_lf2c["Factor"] = "Irrigation"
TxR_lf2c["Factor"] = "Rootstock x Compartment"
TxI_lf2c["Factor"] = "Compartment x Irrigation"
RxI_lf2c["Factor"] = "Rootstock x Irrigation"

# Bind by row all dfs
Log2fold_by_factor_16s <- rbind(R_lf2c, T_lf2c, B_lf2c, I_lf2c, TxR_lf2c, TxI_lf2c, RxI_lf2c)

# Tukey test for significant difference in means 
TukeyHSD(aov(Log2FoldChange ~ Factor, data = Log2fold_by_factor_16s), conf.level = 0.95)

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
# In order for DESEQ2 to run there can be no zeros in the OTU matrix, so I add 1 to each count
otu_table(physeq_its) <- otu_table(physeq_its) + 1
# Create Deseq2 object with study desgin formula
DS2_its <- phyloseq_to_deseq2(physeq_its, ~ Compartment*Rootstock + Irrigation + Irrigation:Compartment + Irrigation:Rootstock + Block)
DS2_its <- DESeq(DS2_its, test="Wald", fitType = "local")

# Results table
resultsNames(DS2_its)
# Contrast list
Contrast_list <- c("Compartment_Leaf_vs_Berry", "Compartment_Root_vs_Berry", "Compartment_Soil_vs_Berry", "Rootstock_1103P_vs_Ungrafted", "Rootstock_3309C_vs_Ungrafted", "Rootstock_SO4_vs_Ungrafted", "Irrigation_none_vs_full", "Irrigation_rdi_vs_full", "Block_B_vs_A", "Block_C_vs_A", "CompartmentLeaf.Rootstock1103P", "CompartmentRoot.Rootstock1103P", "CompartmentSoil.Rootstock1103P", "CompartmentLeaf.Rootstock3309C", "CompartmentRoot.Rootstock3309C", "CompartmentSoil.Rootstock3309C", "CompartmentLeaf.RootstockSO4", "CompartmentRoot.RootstockSO4", "CompartmentSoil.RootstockSO4", "CompartmentLeaf.Irrigationnone", "CompartmentRoot.Irrigationnone", "CompartmentSoil.Irrigationnone", "CompartmentLeaf.Irrigationrdi", "CompartmentRoot.Irrigationrdi", "CompartmentSoil.Irrigationrdi", "Rootstock1103P.Irrigationnone", "Rootstock3309C.Irrigationnone", "RootstockSO4.Irrigationnone", "Rootstock1103P.Irrigationrdi", "Rootstock3309C.Irrigationrdi", "RootstockSO4.Irrigationrdi")
# Rough plotting of the number of significantly different abundant ASVs per factor
p <- c()
for (i in Contrast_list){
  p[[i]] <- as.data.frame(results(DS2_its, cooksCutoff = FALSE, name = i))
  p[[i]] <- p[[i]][which(p[[i]]$padj < 0.05), ]
}

# Get Number of differentailly abundant ASVs per factor
# Main effects Compartment, Rootstock, Irrigation, Block
T_num <- length(unique(c(rownames(p[["Compartment_Leaf_vs_Berry"]]), rownames(p[["Compartment_Root_vs_Berry"]]), rownames(p[["Compartment_Soil_vs_Berry"]]))))
R_num <- length(unique(c(rownames(p[["Rootstock_1103P_vs_Ungrafted"]]), rownames(p[["Rootstock_3309C_vs_Ungrafted"]]), rownames(p[["Rootstock_SO4_vs_Ungrafted"]]))))
I_num <- length(unique(c(rownames(p[["Irrigation_none_vs_full"]]), rownames(p[["Irrigation_rdi_vs_full"]]))))
B_num <- length(unique(c(rownames(p[["Block_B_vs_A"]]), rownames(p[["Block_C_vs_A"]]))))
# Interactions
RxT_num <- length(unique(c(rownames(p[["CompartmentLeaf.Rootstock1103P"]]), rownames(p[["CompartmentRoot.Rootstock1103P"]]), rownames(p[["CompartmentSoil.Rootstock1103P"]]), rownames(p[["CompartmentLeaf.Rootstock3309C"]]), rownames(p[["CompartmentRoot.Rootstock3309C"]]), rownames(p[["CompartmentSoil.Rootstock3309C"]]), rownames(p[["CompartmentLeaf.RootstockSO4"]]), rownames(p[["CompartmentRoot.RootstockSO4"]]), rownames(p[["CompartmentSoil.RootstockSO4"]]))))
TxI_num <- length(unique(c(rownames(p[["CompartmentLeaf.Irrigationnone"]]), rownames(p[["CompartmentRoot.Irrigationnone"]]), rownames(p[["CompartmentSoil.Irrigationnone"]]), rownames(p[["CompartmentLeaf.Irrigationrdi"]]), rownames(p[["CompartmentRoot.Irrigationrdi"]]), rownames(p[["CompartmentSoil.Irrigationrdi"]]))))
RxI_num <- length(unique(c(rownames(p[["Rootstock1103P.Irrigationnone"]]), rownames(p[["Rootstock3309C.Irrigationnone"]]), rownames(p[["RootstockSO4.Irrigationnone"]]), rownames(p[["Rootstock1103P.Irrigationrdi"]]), rownames(p[["Rootstock3309C.Irrigationrdi"]]), rownames(p[["RootstockSO4.Irrigationrdi"]]))))

Number_of_DiffAbund_ASVs_its <- as.data.frame(c(T_num, R_num, I_num, B_num, RxT_num, TxI_num, RxI_num))
rownames(Number_of_DiffAbund_ASVs_its) <- c( "Compartment", "Rootstock", "Irrigation", "Block", "Rootstock x Compartment", "Compartment x Irrigation", "Rootstock x Irrigation")
colnames(Number_of_DiffAbund_ASVs_its) <- "num_ASVs"
Number_of_DiffAbund_ASVs_its <- cbind(Number_of_DiffAbund_ASVs_its, Factor = c("Compartment", "Rootstock", "Irrigation", "Block", "Rootstock x Compartment", "Compartment x Irrigation", "Rootstock x Irrigation"))
Number_of_DiffAbund_ASVs_its$relabun <- (Number_of_DiffAbund_ASVs_its$num_ASVs / 111) # This 111 number comes from the number of unique taxa from filtering the phyloseq object above

# Get Taxonomy of differentailly abundant ASVs per factor
# Main effects Compartment, Rootstock, Irrigation, Block
T_ASV_list <- as.data.frame(unique(c(rownames(p[["Compartment_Leaf_vs_Berry"]]), rownames(p[["Compartment_Root_vs_Berry"]]), rownames(p[["Compartment_Soil_vs_Berry"]]))))
R_ASV_list <- as.data.frame(unique(c(rownames(p[["Rootstock_1103P_vs_Ungrafted"]]), rownames(p[["Rootstock_3309C_vs_Ungrafted"]]), rownames(p[["Rootstock_SO4_vs_Ungrafted"]]))))
I_ASV_list <- as.data.frame(unique(c(rownames(p[["Irrigation_none_vs_full"]]), rownames(p[["Irrigation_rdi_vs_full"]]))))
B_ASV_list <- as.data.frame(unique(c(rownames(p[["Block_B_vs_A"]]), rownames(p[["Block_C_vs_A"]]))))
# Interactions
RxT_ASV_list <- as.data.frame(unique(c(rownames(p[["CompartmentLeaf.Rootstock1103P"]]), rownames(p[["CompartmentRoot.Rootstock1103P"]]), rownames(p[["CompartmentSoil.Rootstock1103P"]]), rownames(p[["CompartmentLeaf.Rootstock3309C"]]), rownames(p[["CompartmentRoot.Rootstock3309C"]]), rownames(p[["CompartmentSoil.Rootstock3309C"]]), rownames(p[["CompartmentLeaf.RootstockSO4"]]), rownames(p[["CompartmentRoot.RootstockSO4"]]), rownames(p[["CompartmentSoil.RootstockSO4"]]))))
TxI_ASV_list <- as.data.frame(unique(c(rownames(p[["CompartmentLeaf.Irrigationnone"]]), rownames(p[["CompartmentRoot.Irrigationnone"]]), rownames(p[["CompartmentSoil.Irrigationnone"]]), rownames(p[["CompartmentLeaf.Irrigationrdi"]]), rownames(p[["CompartmentRoot.Irrigationrdi"]]), rownames(p[["CompartmentSoil.Irrigationrdi"]]))))
RxI_ASV_list <- as.data.frame(unique(c(rownames(p[["Rootstock1103P.Irrigationnone"]]), rownames(p[["Rootstock3309C.Irrigationnone"]]), rownames(p[["RootstockSO4.Irrigationnone"]]), rownames(p[["Rootstock1103P.Irrigationrdi"]]), rownames(p[["Rootstock3309C.Irrigationrdi"]]), rownames(p[["RootstockSO4.Irrigationrdi"]]))))
# Get dataframe of all ASV hash number for use in next step
x <- as.data.frame(tax_table(physeq_its))
x$ASV_hash <- rownames(x)
# Main effects Compartment, Rootstock, Irrigation, Block
T_ASV_list <- x[x$ASV_hash %in% T_ASV_list$uniq, ]
R_ASV_list <- x[x$ASV_hash %in% R_ASV_list$uniq, ]
I_ASV_list <- x[x$ASV_hash %in% I_ASV_list$uniq, ]
B_ASV_list <- x[x$ASV_hash %in% B_ASV_list$uniq, ]
# Interaction
RxT_ASV_list <- x[x$ASV_hash %in% RxT_ASV_list$uniq, ]
TxI_ASV_list <- x[x$ASV_hash %in% TxI_ASV_list$uniq, ]
RxI_ASV_list <- x[x$ASV_hash %in% RxI_ASV_list$uniq, ]
# Extracting taxonomy at the class level
Class_ASV_df <- as.data.frame(bind_rows(summary(T_ASV_list$Class), summary(R_ASV_list$Class), summary(I_ASV_list$Class), summary(B_ASV_list$Class), summary(RxT_ASV_list$Class), summary(TxI_ASV_list$Class), summary(RxI_ASV_list$Class)))
Class_ASV_df[is.na(Class_ASV_df)] <- 0
# Merge relative proportion of ASVs df with taxonomy df
num_tax_ASVs_its <- cbind(Number_of_DiffAbund_ASVs_its, Class_ASV_df)
# Covert taxonomy counts into proportions
dim(num_tax_ASVs_its)
num_tax_ASVs_its[1,c(seq(4,16,1))] <- num_tax_ASVs_its[1,3] * (num_tax_ASVs_its[1,c(seq(4,16,1))] / sum(num_tax_ASVs_its[1,c(seq(4,16,1))]))
num_tax_ASVs_its[2,c(seq(4,16,1))] <- num_tax_ASVs_its[2,3] * (num_tax_ASVs_its[2,c(seq(4,16,1))] / sum(num_tax_ASVs_its[2,c(seq(4,16,1))]))
num_tax_ASVs_its[3,c(seq(4,16,1))] <- num_tax_ASVs_its[3,3] * (num_tax_ASVs_its[3,c(seq(4,16,1))] / sum(num_tax_ASVs_its[3,c(seq(4,16,1))]))
num_tax_ASVs_its[4,c(seq(4,16,1))] <- num_tax_ASVs_its[4,3] * (num_tax_ASVs_its[4,c(seq(4,16,1))] / sum(num_tax_ASVs_its[4,c(seq(4,16,1))]))
num_tax_ASVs_its[5,c(seq(4,16,1))] <- num_tax_ASVs_its[5,3] * (num_tax_ASVs_its[5,c(seq(4,16,1))] / sum(num_tax_ASVs_its[5,c(seq(4,16,1))]))
num_tax_ASVs_its[6,c(seq(4,16,1))] <- num_tax_ASVs_its[6,3] * (num_tax_ASVs_its[6,c(seq(4,16,1))] / sum(num_tax_ASVs_its[6,c(seq(4,16,1))]))
num_tax_ASVs_its[7,c(seq(4,16,1))] <- num_tax_ASVs_its[7,3] * (num_tax_ASVs_its[7,c(seq(4,16,1))] / sum(num_tax_ASVs_its[7,c(seq(4,16,1))]))
# Gather columns into a single column for use in ggplot
num_tax_ASVs_its <- gather(num_tax_ASVs_its, key = "Taxon", value = "Freq", -c(Factor,relabun, num_ASVs))
# Set the factor order of the taxon names
num_tax_ASVs_its$Taxon <- factor(num_tax_ASVs_its$Taxon, levels = c("Agaricomycetes", "Cystobasidiomycetes", "Dothideomycetes", "Glomeromycetes", "Leotiomycetes", "Microbotryomycetes", "Mortierellomycetes", "Paraglomeromycetes", "Pezizomycetes", "Eurotiomycetes", "Saccharomycetes", "Sordariomycetes", "Tremellomycetes", "NA's"))



# Plotting of log fold changes by factor
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
T_lf2c["Factor"] = "Compartment"
B_lf2c["Factor"] = "Block"
I_lf2c["Factor"] = "Irrigation"
TxR_lf2c["Factor"] = "Rootstock x Compartment"
TxI_lf2c["Factor"] = "Compartment x Irrigation"
RxI_lf2c["Factor"] = "Rootstock x Irrigation"
# Bind by row all dfs
Log2fold_by_factor_its <- rbind(R_lf2c, T_lf2c, B_lf2c, I_lf2c, TxR_lf2c, TxI_lf2c, RxI_lf2c)

# Combination bar plot
num_tax_ASVs_16s$Factor <- factor(num_tax_ASVs_16s$Factor, levels = c("Rootstock", "Compartment", "Irrigation", "Rootstock x Compartment", "Rootstock x Irrigation", "Compartment x Irrigation", "Block"))
num_tax_ASVs_its$Factor <- factor(num_tax_ASVs_its$Factor, levels = c("Rootstock", "Compartment", "Irrigation", "Rootstock x Compartment", "Rootstock x Irrigation", "Compartment x Irrigation", "Block"))
Numb_ASVs_16s_plot <- ggplot(data = num_tax_ASVs_16s, aes(x=Factor, y=Freq, fill=Taxon)) + geom_bar(stat="identity", colour="black", lwd=0.2) + scale_fill_manual(name ="Bacterial Phyla", values=taxapallete_16S) + ylab("Proportion of ASVs") + ylim(0,1) 
Numb_ASVs_its_plot <- ggplot(data = num_tax_ASVs_its, aes(x=Factor, y=Freq, fill=Taxon)) + geom_bar(stat="identity", colour="black", lwd=0.2) + scale_fill_manual(name ="Fungal Class", values=taxapallete_ITS) + ylab("Proportion of ASVs") + ylim(0,1) 

# Combination violin plot
Log2fold_by_factor_16s["Marker"] = "16S"
Log2fold_by_factor_its["Marker"] = "ITS"
Log2fold_by_factor_both <- rbind(Log2fold_by_factor_16s, Log2fold_by_factor_its)
Log2fold_by_factor_both$Factor <- as.factor(Log2fold_by_factor_both$Factor)
Log2fold_by_factor_both$Factor <- factor(Log2fold_by_factor_both$Factor, levels = c("Rootstock", "Compartment", "Irrigation", "Rootstock x Compartment", "Rootstock x Irrigation", "Compartment x Irrigation", "Block"))
Log2fold_both_plot <- ggplot(Log2fold_by_factor_both, aes(x=Factor, y=Log2FoldChange, fill=Marker)) + geom_violin(trim = FALSE, scale = "width") + stat_summary(fun.data=mean_sdl, geom="pointrange", position = position_dodge(width = 0.9))+ xlab("Factor") + scale_fill_manual(name = "Marker", values=c("gray36", "gray80"))  + labs(y=expression(paste(Log[2]," fold change")), x=("Source of variation"))

Bar_plot_and_log2foldchange_plot <- (Numb_ASVs_16s_plot | Numb_ASVs_its_plot) / Log2fold_both_plot
ggsave(filename = "Bar_plot_and_log2foldchange_plot.pdf", plot = Bar_plot_and_log2foldchange_plot, units = "in", height = 12, width = 24, path = "Figures")


# Mean and Stdev for log2 fold changes for bacteria and fungi
summary(Log2fold_by_factor_16s)
sd(Log2fold_by_factor_16s$Log2FoldChange)
summary(Log2fold_by_factor_its)
sd(Log2fold_by_factor_its$Log2FoldChange)

# Mean and sd by each factor and interactions 
tapply(Log2fold_by_factor_16s$Log2FoldChange, Log2fold_by_factor_16s$Factor, mean)
tapply(Log2fold_by_factor_16s$Log2FoldChange, Log2fold_by_factor_16s$Factor, sd)
tapply(Log2fold_by_factor_its$Log2FoldChange, Log2fold_by_factor_its$Factor, mean)
tapply(Log2fold_by_factor_its$Log2FoldChange, Log2fold_by_factor_its$Factor, sd)

# Testing for significant difference in factor log2foldchange means
TukeyHSD(aov(Log2FoldChange ~ Factor, data = Log2fold_by_factor_16s))
TukeyHSD(aov(Log2FoldChange ~ Factor, data = Log2fold_by_factor_its))


##### 4.0) Specific associations Acetobacterales and Saccharomycetes #####
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
## Removing 7 biological samples

## Fix taxonomy strings 
# Bacteria
rare_physeq_16s@tax_table[,1]  <- gsub('D_.__', '', rare_physeq_16s@tax_table[,1])
tax_table(rare_physeq_16s) <- as.matrix(separate(as.data.frame(rare_physeq_16s@tax_table[,1]), col = Kingdom, sep = ";", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")))
# Fungi
rare_physeq_its@tax_table[,1] <- gsub('.__', '', rare_physeq_its@tax_table[,1])
tax_table(rare_physeq_its) <- as.matrix(separate(as.data.frame(rare_physeq_its@tax_table[,1]), col = Kingdom, sep = ";", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")))

# Bacteria: Acetobacterales
# Subset to just Acetobacterales
acetobacterales_physeq16s <- subset_taxa(rare_physeq_16s, Order == "Acetobacterales")
plot_bar(acetobacterales_physeq16s, "ExtractionNum", "Abundance", "Genus", facet_grid="Rootstock~.")
# Check proportion of the abudance is driven by Gluconobacter and Acetobacter
# sum(rowSums(t(as.data.frame(otu_table(acetobacterales_physeq16s))))) #6809
# Gluconobacter_physeq16s <- subset_taxa(acetobacterales_physeq16s, Genus == "Gluconobacter")
# sum(rowSums(t(as.data.frame(otu_table(Gluconobacter_physeq16s))))) #5678
# Acetobacter_physeq16s <- subset_taxa(acetobacterales_physeq16s, Genus == "Acetobacter")
# sum(rowSums(t(as.data.frame(otu_table(Acetobacter_physeq16s))))) #417

# Collapse reads for members of Acetobacterales into a single sum
x <- t(as.data.frame(otu_table(acetobacterales_physeq16s)))
x <- as.data.frame(rowSums(x))
Aceto_abund <- cbind(x, sample_data(rare_physeq_16s))

# Rename rowsum column to Acetobacterales
colnames(Aceto_abund)[1] <- "Acetobacterales_abundance"

# Make abundance relative by dividing by the rarefaction amount (1500)
Aceto_abund$Acetobacterales_abundance <- Aceto_abund$Acetobacterales_abundance/1500

# Linear model full
lm_acetobacterales_full <- lm(Acetobacterales_abundance ~ Compartment*Rootstock*Irrigation + Block, data = Aceto_abund)
Anova(lm_acetobacterales_full, type = "III")

# Post-hoc testing
pairs(emmeans(lm_acetobacterales_full, ~ Rootstock | Compartment| Irrigation))

# Plots to investigate significant results of lm's
# Three way interaction from full model
Acetobacterales_abundance.graph <- ggplot(Aceto_abund, aes(Compartment, Acetobacterales_abundance, fill = Rootstock)) + geom_boxplot(outlier.shape = NA) + xlab ("Compartment") + ylab("Acetobacterales relative abundance") + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette) + theme(legend.position = "right")
Acetobacterales_abundance.graph <- Acetobacterales_abundance.graph + facet_wrap("Irrigation")
ggsave("Acetobacterales_abundance.png", Acetobacterales_abundance.graph, units = "in", height = 6, width = 12)

# Fungi: Saccharomycetes
# Subset to just Saccharomycetes
Saccharomycetes_physeqits <- subset_taxa(rare_physeq_its, Class == "Saccharomycetes")
# Quick plots
plot_bar(Saccharomycetes_physeqits, "ExtractionNum", "Abundance", "Family", facet_grid="Rootstock~.")
plot_bar(Saccharomycetes_physeqits, "ExtractionNum", "Abundance", "Genus", facet_grid="Rootstock~.")
#Candida_physeqits <- subset_taxa(Saccharomycetes_physeqits, Genus == "Candida")
#sum(rowSums(t(as.data.frame(otu_table(Candida_physeqits))))) #744
#Hanseniaspora_physeqits <- subset_taxa(Saccharomycetes_physeqits, Genus == "Hanseniaspora")
#sum(rowSums(t(as.data.frame(otu_table(Hanseniaspora_physeqits))))) #13562
#Pichia_physeqits <- subset_taxa(Saccharomycetes_physeqits, Genus == "Pichia")
#sum(rowSums(t(as.data.frame(otu_table(Pichia_physeqits))))) #10578

# Collapse reads for members of Acetobacterales into a single sum
x <- t(as.data.frame(otu_table(Saccharomycetes_physeqits)))
x <- as.data.frame(rowSums(x))
Sacchro_abund <- cbind(x, sample_data(rare_physeq_its))

# Rename rowsum column to Acetobacterales
colnames(Sacchro_abund)[1] <- "Saccharomycetes_abundance"

# Make abundance relative by dividing by the rarefaction amount (1500)
Sacchro_abund$Saccharomycetes_abundance <- Sacchro_abund$Saccharomycetes_abundance/5000

# Linear model full
lm_saccharomycetes_full <- lm(Saccharomycetes_abundance ~ Compartment*Rootstock*Irrigation + Block, data = Sacchro_abund)
Anova(lm_saccharomycetes_full, type = "III")

# Post-hoc testing
pairs(emmeans(lm_saccharomycetes_full, ~ Rootstock | Compartment))
pairs(emmeans(lm_saccharomycetes_full, ~ Rootstock))

# Plot
Saccharomycetes_abundance.graph <- ggplot(Sacchro_abund, aes(Compartment, Saccharomycetes_abundance, fill = Rootstock)) + geom_boxplot(outlier.shape = NA) + xlab ("Compartment") + ylab("Saccharomycetes relative abundance") + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette) + theme(legend.position = "right")
Saccharomycetes_abundance.graph <- Saccharomycetes_abundance.graph + facet_wrap("Irrigation")

ggplot(Aceto_abund, aes(Compartment, Acetobacterales_abundance, fill = Rootstock)) + geom_boxplot(outlier.shape = NA) + xlab ("Compartment") + ylab("Acetobacterales relative abundance") + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette) + theme(legend.position = "right")
ggplot(Sacchro_abund, aes(Irrigation, Saccharomycetes_abundance, fill = Rootstock)) + geom_boxplot(outlier.shape = NA) + xlab ("Rootstock") + ylab("Saccharomycetes relative abundance") + geom_point(position = position_jitterdodge(jitter.width = 0.2)) + scale_fill_manual(name = "Rootstock", values=zoe_palette) + theme(legend.position = "right")

Aceto_sacchro.graph <- ggarrange(Acetobacterales_abundance.graph, Saccharomycetes_abundance.graph, nrow = 2, labels ="AUTO")
ggsave("Acetobacterales_Saccharomycetes_byRandI.png", Aceto_sacchro.graph, units = "in", height = 8, width = 12)

# Looking at both aceto and sacchro
X <- data.frame("Sample" = rownames(Sacchro_abund), "Saccharomycetes_abund" = Sacchro_abund$Saccharomycetes_abundance, "Rootstock" = Sacchro_abund$Rootstock, "Compartment" = Sacchro_abund$Compartment)
Y <- data.frame("Sample" = rownames(Aceto_abund), "Acetobacterales_abund" = Aceto_abund$Acetobacterales_abundance, "Rootstock" = Aceto_abund$Rootstock, "Compartment" = Aceto_abund$Compartment)

# Check which samples are not in both lists
setdiff(X$Sample, Y$Sample)

# Merge only samples present within both lists 
combo_aceto_sacchro <- merge(X, Y, by = "Sample")

Compartment.abundance.Ace.Sac <- ggplot(combo_aceto_sacchro, aes(Saccharomycetes_abund, Acetobacterales_abund, color = Compartment.x, fill = Compartment.x)) +
  geom_smooth(size = 2, method = lm) + stat_cor(method = "spearman") +
  geom_point(size = 3, shape = 21, color = "Black") +
  labs(fill = "Compartment", color = "Compartment", x = "Saccharomycetes relative abundance", y = "Acetobacterales relative abundance")

Compartment_by_rootstock.abundance.Ace.Sac <- ggplot(combo_aceto_sacchro, aes(Saccharomycetes_abund, Acetobacterales_abund, color = Compartment.x, fill= Compartment.x)) +
  geom_smooth(size = 2, method = lm) + 
  stat_cor(method = "spearman") +
  geom_point(size = 3, shape = 21, color = "Black") + facet_wrap("Rootstock.x") +
  labs(fill = "Compartment", color = "Compartment", x = "Saccharomycetes relative abundance", y = "Acetobacterales relative abundance")

##### 5.0) Combining plots and saving #####

#Alpha diversity plots#
#16s by Compartment and rootsock
alpha_16s <- ggarrange(fphd_byR_T_16s, simpI_byR_T_16s, obso_byR_T_16s, shan_byR_T_16s, labels = c("A","B","C","D"), ncol = 2 , nrow = 2, common.legend = TRUE, legend = "right")
ggsave(filename = "16s_alphadiv.pdf", plot = alpha_16s, units = "in", height = 12, width = 12, path = "Figures")

#its by Compartment and rootstock
alpha_its <- ggarrange(simpI_byR_T_its, obso_byR_T_its, shan_byR_T_its, labels = c("A","B","C"), common.legend = TRUE, legend = "right")
ggsave(filename = "its_alphadiv.pdf", plot = alpha_its, units = "in", height = 12, width = 12, path = "Figures")

#Combine by Compartment and Rootstock
ggarrange(simpI_byR_T_16s, obso_byR_T_16s, shan_byR_T_16s, fphd_byR_T_16s, simpI_byR_T_its, obso_byR_T_its, shan_byR_T_its, labels = c("A","B","C","D","E","F","G"), ncol = 4 , nrow = 2, common.legend = TRUE, legend = "right")

#Combine by Compartment
alpha_16s_its <- ggarrange(shan_16s, obso_16s, fphd_16s, shan_its, obso_its, labels = c("A","B","C","D","E"), ncol = 3, nrow = 2)
ggsave(filename = "16s-its_alphadiv.pdf", plot = alpha_16s_its, units = "in", height = 12, width = 12, path = "Figures")

#Taxonomic bar plots#
#16s and ITS by Compartment
taxabp_16s_its <- ggarrange(taxa_barplot16s, taxa_barplotits, labels = c("A", "B"), ncol = 2)
ggsave(filename = "16s-its_taxonomic_barplot.pdf", plot = taxabp_16s_its, units = "in", height = 6, width = 12, path = "Figures")

#16s by Compartment and rootstock
taxabp_16s <- ggarrange(taxa_barplot16s_S04, taxa_barplot16s_3309C, taxa_barplot16s_1103P, taxa_barplot16s_Ungrafted, labels = c("   S04", "3309C", "1103P", "  Ungrafted"), label.x = 0.15, common.legend = TRUE, legend = "right")
ggsave(filename = "16s_taxonomic_barplot.pdf", plot = taxabp_16s, units = "in", height = 12, width = 12, path = "Figures")

#ITS by Compartment and rootstock
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
## Compartment combo figure, simpson diversity, PCoAs, and Taxonomic barplots
Compartment_16s_combo <- ggarrange(simpI_16s, PCoA_1_2_uunif_16s_c, taxa_barplot16s, widths = c(1,1.5,1), ncol = 3, labels = c("A", "B", "C"))
Compartment_its_combo <- ggarrange(simpI_its, PCoA_1_2_bray_its_c, taxa_barplotits,  widths = c(1,1.5,1), ncol = 3, labels = c("D", "E", "F"))
Compartment_16s_its <- ggarrange(Compartment_16s_combo, Compartment_its_combo, nrow = 2)
ggsave(filename = "Compartment_16s_its_simp_PCoA_taxabarplot.png", plot = Compartment_16s_its, units = "in", height = 18, width = 24, path = "Figures")
ggsave(filename = "Compartment_16s_its_simp_PCoA_taxabarplot.pdf", plot = Compartment_16s_its, units = "in", height = 18, width = 24, path = "Figures")

# Figure 3
## Rootstock combo figure, simpson diversity with individual rootstock barplots by tissue
#16s by Compartment and rootstock
taxabp_16s <- ggarrange(taxa_barplot16s_S04, taxa_barplot16s_3309C, taxa_barplot16s_1103P, taxa_barplot16s_Ungrafted, labels = c("   S04", "3309C", "1103P", "  Ungrafted"), label.x = 0.15, common.legend = TRUE, legend = "right", ncol = 4, nrow = 1)
taxabp_16s <- ggarrange(simpI_byR_T_16s, taxabp_16s, labels = c("A", "B"), widths = c(0.5,2))
#ggsave(filename = "16s_taxonomic_barplot_w_simpson.pdf", plot = taxabp_16s, units = "in", height = 12, width = 24, path = "Figures")
#ITS by Compartment and rootstock
taxabp_its <- ggarrange(taxa_barplotits_S04, taxa_barplotits_3309C, taxa_barplotits_1103P, taxa_barplotits_Ungrafted, labels = c("   S04", "3309C", "1103P", "  Ungrafted"), label.x = 0.15, common.legend = TRUE, legend = "right", ncol = 4, nrow = 1)
taxabp_its <- ggarrange(simpI_byR_T_its, taxabp_its, labels = c("C", "D"),  widths = c(0.5,2))
#ggsave(filename = "its_taxonomic_barplot.pdf", plot = taxabp_its, units = "in", height = 12, width = 24, path = "Figures")
#16s/ITS Together
taxabp_16s_its <- ggarrange(taxabp_16s, taxabp_its, nrow = 2)
ggsave(filename = "Rootstock_simp_taxonomicBP.pdf", plot = taxabp_16s_its, units = "in", height = 18, width = 24, path = "Figures")

# Figure 4
# Acetobacterales
Combo_correlation_plot <- ggarrange(Compartment.abundance.Ace.Sac, common.legend = TRUE, legend = "right")
ggsave("Correlation_plots_ACETO_SACCHRO_by_compartment.pdf", Combo_correlation_plot, width = 10, height = 8)