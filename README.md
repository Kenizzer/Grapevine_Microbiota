# Grapevine Microbiome

## Abridged abstract

An important question is the extent to which unique root system and shoot system genotypes, when grafted together, influence the microbiota of the graft partner. Our study sought to answer this question by utilizing an experimental vineyard composed of ‘Chambourcin’ vines growing ungrafted and grafted to three different rootstocks, replicated across three irrigation treatments. We characterized bacterial and fungal communities in roots, leaves, and berries, as well as surrounding soil. Our objectives were to **(1)** characterize the microbiota of compartments within the root system (roots and adjacent soil) and the shoot system (leaves and berries), **(2)** determine the influence of rootstock genotypes, irrigation, and their interaction on the microbiota of aboveground and belowground compartments, and **(3)** investigate the distribution of microorganisms implicated in the late-season grapevine bunch rot disease sour rot (Acetobacterales and Saccharomycetes).

Table of markers used within the study.

| Marker   | F Primer| F Primer Sequence       | R Primer| R Primer Sequence    | Size (bp) | Citation                                       |
|:--------:|:-------:|:-----------------------:|:-------:|:--------------------:|:---------:|:----------------------------------------------:|
| 16s      | 515F    | GTGYCAGCMGCCGCGGTAA     | 806R    | GGACTACNVGGGTWTCTAAT | 390       | Parada *et al.* 2016 and Apprill *et al.* 2015 |
| trnL     | ITS1f   | CTTGGTCATTTAGAGGAAGTAA  | ITS2    | GCTGCGTTCTTCATCGATGC | 250–600   | Smith and Peay 2014 and White *et al.* 1990    |

## Experimental and Sampling Design

**A)** Vineyard layout at the University of Missouri Southwest Research Center in Mount Vernon, MO consisting of 288 vines grafted to one of three rootstocks (‘SO4’, Selection Oppenheim 4; ‘3309C’, 3309 Courdec; ‘1103P’, 1103 Paulson) or ungrafted (UG). Each colored cell represents four vines; irrigation treatments (F, Full replacement of evapotranspiration; I, reduced replacement; N, unirrigated) and experimental blocks are listed along the bottom and top of the grid, respectively. **B)** Depiction of a grafted grapevine with sampling comparments highlighted, numbers correspond to legend (left).

![Image of experimental and sampling design](https://github.com/Kenizzer/Grapevine_microbiome/blob/master/Experimenetal_design_image/experimental_design_and_sampling_figure.png)

## Analysis 

### QIIME2 processing
***See Qiime_analysis.md***

Data processing was conducted in QIIME2 (Bolyen *et al.* 2018). Samples were demultiplexed (*qiime demux emp-paired*) according to barcode sequence (see metadata files). For 16S rRNA gene sequences, the QIIME2 plugin **DADA2** (Callahan *et al.* 2016) was used to denoise, dereplicate, and filter chimeric sequences. The first 13 nt of each sequence (forward and reverse) was trimmed and truncated at 150 nt to remove lower quality bases. For ITS sequences, **Cutadapt** (Martin 2011) was used to remove primer sequences from sequences prior to using DADA2. The first 12 bp of the 3′ end of each sequence were trimmed, the 5′ end of sequences were not truncated to preserve biologically relevant length variation (Schoch *et al.* 2012), and a max expected error rate of 3 was used. 

Taxonomic classification of ASVs was conducted with a naive Bayes classifier trained on either the **SILVA** (v.132; trimmed to the V4 region; Yilmaz *et al.* 2014) or **UNITE** database (v.8.0; Unite Community 2018) for bacteria and fungi, respectively. Bacterial and fungal ASVs not assigned to a phylum were removed along with bacterial ASVs assigned to mitochondria or chloroplasts. ASVs with less than 0.1% of the total filtered reads were removed.

### R analysis
***See R_analysis.r***

Samples were rarefied to 1,500 or 5,000 sequences for bacteria and fungi, respectively. Alpha-diversity and beta-diversity metrics were calculated in **phyloseq** (McMurdie and Holmes 2013).Principle coordinates analysis was used to visualize sample relationships. Venn diagrams were created by determining the intersection of ASV lists per compartment (URL: http://bioinformatics.psb.ugent.be/webtools/Venn/).

Linear models and PERMANOVA tests were run using the formula: response variable ~ Rootstock genotype × Compartment × Irrigation + Block. For each alpha-diversity metric, a linear model was fit using lm and assessed using type-III ANOVA with the **car** package (v3.0-3; Fox and Weisberg 2019). Beta-diversity metrics were subjected to PERMANOVA tests with adonis in the **vegan** package (Anderson, M. J. 2001; Oksanen *et al.* 2019).

For differential abundance analysis, we used unrarefied reads and removed ASVs that were not represented by a depth of least 25 reads in more than 10% of the samples to conduct differential abundance modeling with **DESeq2** (v.1.24.0; Love *et al.* 2014). DESeq2 was fit using the following model: response variable ~ Rootstock genotype (R) + Compartment (C) + Irrigation (I) + R × C + I × R + I × C + Block (B). For each factor we extracted the number of ASVs that showed a significant pattern of differential abundance as well as their fold changes (Log2) to generate summary plots. 

### Machine Learning
***See Machine_learing.r***

For the machine learning we used **ranger’s** implementation of random forest (v.0.11.2; Wright and Ziegler 2017) and tuned hyperparameters (number of trees, minimum node size, and number of features available at each node) with **Caret** (v.6.0-84; Kuhn 2008) on 80% of the dataset (20% withheld for testing). The optimal hyperparameters were selected by iteratively assessing the performance on out-of-bag samples for a given parameter set. We then used this final model to predict the label, either rootstock or compartment or both, on the withheld testing data assessing the prediction accuracy. We tested both regions (16S rRNA gene and ITS) using the above framework both separately and together. We found that the accuracy was similar for each region and the combined dataset, so we chose to report the results of the combined dataset. Visualizations were made using tile plots from the output confusion matrix. 

**Citations**
---
  * Anderson, M. J. 2001. A new method for non-parametric multivariate analysis of variance. Austral Ecology 26: 32–46. https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x
  * Apprill, A., McNally, S., Parsons, R., & Weber, L. (2015). Minor revision to V4 region SSU rRNA 806R gene primer greatly increases detection of SAR11 bacterioplankton. Aquatic Microbial Ecology, 75(2), 129–137. http://doi.org/10.3354/ame01753
  * Bolyen E, Rideout JR, Dillon MR, Bokulich NA, Abnet CC, Al-Ghalith GA, Alexander H, Alm EJ, Arumugam M, Asnicar F, Bai Y, Bisanz JE, Bittinger K, Brejnrod A, Brislawn CJ, Brown CT, Callahan BJ, Caraballo-Rodríguez AM, Chase J, Cope EK, Da Silva R, Diener C, Dorrestein PC, Douglas GM, Durall DM, Duvallet C, Edwardson CF, Ernst M, Estaki M, Fouquier J, Gauglitz JM, Gibbons SM, Gibson DL, Gonzalez A, Gorlick K, Guo J, Hillmann B, Holmes S, Holste H, Huttenhower C, Huttley GA, Janssen S, Jarmusch AK, Jiang L, Kaehler BD, Kang KB, Keefe CR, Keim P, Kelley ST, Knights D, Koester I, Kosciolek T, Kreps J, Langille MGI, Lee J, Ley R, Liu YX, Loftfield E, Lozupone C, Maher M, Marotz C, Martin BD, McDonald D, McIver LJ, Melnik AV, Metcalf JL, Morgan SC, Morton JT, Naimey AT, Navas-Molina JA, Nothias LF, Orchanian SB, Pearson T, Peoples SL, Petras D, Preuss ML, Pruesse E, Rasmussen LB, Rivers A, Robeson MS, Rosenthal P, Segata N, Shaffer M, Shiffer A, Sinha R, Song SJ, Spear JR, Swafford AD, Thompson LR, Torres PJ, Trinh P, Tripathi A, Turnbaugh PJ, Ul-Hasan S, van der Hooft JJJ, Vargas F, Vázquez-Baeza Y, Vogtmann E, von Hippel M, Walters W, Wan Y, Wang M, Warren J, Weber KC, Williamson CHD, Willis AD, Xu ZZ, Zaneveld JR, Zhang Y, Zhu Q, Knight R, and Caporaso JG. (2019). Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology 37: 852–857. https://doi.org/10.1038/s41587-019-0209-9
  * Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: high-resolution sample inference from Illumina amplicon data. Nature methods, 13(7), 581-583. https://doi.org/10.1038/nmeth.3869
  * Fox J, Weisberg S. 2019. An R Companion to Applied Regression. Thousands Oaks, CA: Sage.
  * Kuhn M. 2008. Building Predictive Models in R Using the caret Package. Journal of Statistical Software 28: 159–160.
  * Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome biology 15: 550. https://doi.org/10.1186/s13059-014-0550-8
  * Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. journal, 17(1), 10-12. https://doi.org/10.14806/ej.17.1.200
  * McMurdie PJ, Holmes S. 2013. Phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLoS ONE 8. https://doi.org/10.1371/journal.pone.0061217
  * Nilsson RH, Larsson K-H, Taylor AFS, Bengtsson-Palme J, Jeppesen TS, Schigel D, Kennedy P, Picard K, Glöckner FO, Tedersoo L, Saar I, Kõljalg U, Abarenkov K. (2018). The UNITE database for molecular identification of fungi: handling dark taxa and parallel taxonomic classifications. Nucleic Acids Research, https://doi.org/10.1093/nar/gky1022
  * Oksanen, J., Blanchet, F.G., Kindt, R., Legendre, P., Minchin, P.R., O’hara, R.B., Simpson, G.L., Solymos, P., Stevens, M.H.H., and Wagner, H., (2019). vegan: Community Ecology Package.
  * Schoch, C.L., Seifert, K.A., Huhndorf, S., Robert, V., Spouge, J.L., Levesque, C.A., Chen, W. and Fungal Barcoding Consortium, (2012). Nuclear ribosomal internal transcribed spacer (ITS) region as a universal DNA barcode marker for Fungi. Proceedings of the National Academy of Sciences, 109(16), 6241-6246. https://doi.org/10.1073/pnas.1117018109
  * Smith, D. P., & Peay, K. G. (2014). Sequence depth, not PCR replication, improves ecological inference from next generation DNA sequencing. PLoS ONE, 9(2), e90234–e90234. http://doi.org/10.1371/journal.pone.0090234
  * White, T. J., Bruns, T., Lee, S., & Taylor, J. (1990). Amplification and direct sequencing of fungal ribosomal RNA genes for phylogenetics. In PCR protocols: a guide to methods and applications (pp. 315–322). New York: Academic Press.
  * Wright MN, Ziegler A. 2017. Ranger: A fast implementation of random forests for high dimensional data in C++ and R. Journal of Statistical Software 77. https://arxiv.org/ct?url=https%3A%2F%2Fdx.doi.org%2F10.18637%2Fjss.v077.i01&v=c650da44
  * Yilmaz, P., Parfrey, L.W., Yarza, P., Gerken, J., Pruesse, E., Quast, C., Schweer, T., Peplies, J., Ludwig, W. and Glöckner, F.O., (2014). The SILVA and “all-species living tree project (LTP)” taxonomic frameworks. Nucleic acids research, 42(D1), D643-D648. https://doi.org/10.1093/nar/gkt1209
