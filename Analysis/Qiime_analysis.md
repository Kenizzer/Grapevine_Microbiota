# Targeted amplicon sequencing of samples from Mount Vernon Missouri funded by MGWI

#### MD5 Sums

##### 16s
NAME	Undetermined_S0_L001_I1_001.fastq
MD5	3b15e5b086627646b10637b045c0129a
SIZE	1.3 GB

NAME	Undetermined_S0_L001_R1_001.fastq
MD5	3a69d4eae812ce9125b2cd924eb8d067
SIZE	5.6 GB

NAME	Undetermined_S0_L001_R2_001.fastq
MD5	f8deaabb48560f84dae72c952b74a556
SIZE	5.6 GB

##### ITS
NAME	Undetermined_S0_L001_I1_001.fastq
MD5	74d9c5660d5af16e6258296c39865a32
SIZE	1.3 GB

NAME	Undetermined_S0_L001_R1_001_its.fastq
MD5	5416471ed1cd58ab8f09865718caa73e
SIZE	9.0 GB

NAME	Undetermined_S0_L001_R2_001.fastq
MD5	070a2c3afc2a6d5d3309c8cdbe073119
SIZE	9.0 GB


### 16S processing

##### Importing data into a Qiime2 artifact (.qza file)

Files were renamed to fit the Qiime2 EMPPairedEndSequences, for example forward.fastq.gz, reverse.fastq.gz, and barcode.fastq.gz

```bash
qiime tools import --type EMPPairedEndSequences --input-path /home/joel.swift/data/argonne_mtvernon/16s_data/ --output-path 16s-paired-end-sequences.qza
```

##### Initial processing of data

```bash
# I had to assign a new directory as temp to get Qiime demux to run
export TMPDIR='/xfs2/swiftjf/argonne_mtvernon/temp_dir_new/'

# Demultiplexing
qiime demux emp-paired \
  --m-barcodes-file 16s_metadata.tsv \
  --m-barcodes-column BarcodeSequence \
  --i-seqs 16s-paired-end-sequences.qza \
  --o-per-sample-sequences demux.qza \
  --p-rev-comp-mapping-barcodes

qiime demux summarize --i-data demux.qza --o-visualization demux.qzv

# Denoising and merging pair-end reads
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 13 \
  --p-trim-left-r 13 \
  --p-trunc-len-f 150 \
  --p-trunc-len-r 150 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza \
  --verbose \
  --p-n-threads 0

# Removing mock, positive, and negative control samples
qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file 16s_noMockorPosorNeg_metadata.tsv \
  --o-filtered-table filtered-table.qza

# Visualization of sample totals after removing mock, positive, and negative control samples
qiime feature-table summarize \
  --i-table filtered-table.qza \
  --o-visualization filtered-table.qzv \
  --m-sample-metadata-file 16s_metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

# MORE HERE
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```

##### Taxonomy analyses

```bash
# Getting the pretrained classifier
wget https://data.qiime2.org/2019.1/common/silva-132-99-515-806-nb-classifier.qza

# verifing md5sum 9f50514214ffb6fee9d2f87a47a51076
md5sum silva-132-99-515-806-nb-classifier.qza

# Obtaining taxonomy for the 16s sequences
qiime feature-classifier classify-sklearn \
  --i-classifier silva-132-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

# Received a small error when trying to run the command to create a visual from the a trailing whitespace character in the taxonomy.qza file.

# CategoricalMetadataColumn does not support values with leading or trailing whitespace characters. Column 'Taxon' has the following value: D_0__Bacteria;D_1__Bacteroidetes;D_2__Ignavibacteria;D_3__OPB56;D_4__uncultured bacterium 

# To remove this error I removed the trailing whitespace from the data file in the taxnomoy.qza zip file. This regex (\s\t[0-9.]*) was helpful in finding the offending spaces.

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# Determining the number of mitochondrial sequences removed
qiime taxa filter-table \
  --i-table filtered-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria \
  --o-filtered-table table-no-mitochondria.qza

qiime feature-table summarize \
  --i-table table-no-mitochondria.qza \
  --o-visualization table-no-mitochondria.qzv

#Total sequences: 10647004
#After filtering: 9670480
#Difference: 976,524
#% Diff: 9.172%

# Determining the number of Chloroplast sequences removed
qiime taxa filter-table \
  --i-table filtered-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude chloroplast \
  --o-filtered-table table-no-chloroplast.qza

qiime feature-table summarize \
  --i-table table-no-chloroplast.qza \
  --o-visualization table-no-chloroplast.qzv

#Total Sequences: 10647004
#After filtering: 7297866
#Difference: 3349138
#% Diff: 31.456%

# Filtering out unassigned sequences  
qiime taxa filter-table \
  --i-table filtered-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-include D_1__ \
  --o-filtered-table filtered_table-with-phyla.qza

# Filtering out mitochondria and chloroplast sequences
qiime taxa filter-table \
  --i-table filtered_table-with-phyla.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mitochondria-no-chloroplast.qza

qiime feature-table summarize \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --o-visualization table-no-mitochondria-no-chloroplast.qzv

# Removing ASVs that with less than 0.1% mean sample depth
# Mean sample read depth (excluding controls and wine) = 33624
# Min frequency = 34
# Number of ASVs (features) = 20602
# Removed ASVs (features) = 12320
qiime feature-table filter-features --i-table table-no-mitochondria-no-chloroplast.qza \
  --p-min-frequency 34 \
  --p-min-samples 1 \
  --o-filtered-table filtered-table-no-mitochondria-no-chloroplast.qza

qiime feature-table summarize \
  --i-table filtered-table-no-mitochondria-no-chloroplast.qza \
  --o-visualization filtered-table-no-mitochondria-no-chloroplast.qzv

# Creating bar graph
qiime taxa barplot \
  --i-table filtered-table-no-mitochondria-no-chloroplast.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file 16s_noMockorPosorNeg_metadata.tsv \
  --o-visualization taxa-bar-plots-filtered-no-mito-no-chloro.qzv

```

##### Diversity metrics ALL

```bash
# Testing rarefaction sampling depths

# Generation of rarefaction curves
qiime diversity alpha-rarefaction \
  --i-table filtered-table-no-mitochondria-no-chloroplast.qza \
  --p-max-depth 10000 \
  --p-steps 10 \
  --m-metadata-file 16s_noMockorPosorNeg_metadata.tsv \
  --o-visualization rarefaction_curves.qzv

# Samples rarefied to a sampling depth of 2000 (Removing ten biological samples)
# Sample Readcount
# MV16.3309C.14-15.B 1988
# MV15.1103P.10-11.B 1795
# MV09.3309C.10.B 1608
# MV08.OWN.14.B 923
# MV08.1103P.2.L 854
# MV09.OWN.2.B 846
# MV11.SO4.15.B 741
# MV16.OWN.6-7.B 647
# MV09.OWN.3.B 393
# MV08.1103P.2.B 89

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table filtered-table-no-mitochondria-no-chloroplast.qza \
  --p-sampling-depth 2000 \
  --m-metadata-file 16s_noMockorPosorNeg_metadata.tsv \
  --output-dir core-metrics-results-rare2000

# Samples rarefied to a sampling depth of 1500 (Removing seven biological samples)
# Sample Readcount
# MV08.OWN.14.B 923
# MV08.1103P.2.L 854
# MV09.OWN.2.B 846
# MV11.SO4.15.B 741
# MV16.OWN.6-7.B 647
# MV09.OWN.3.B 393
# MV08.1103P.2.B 89

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table filtered-table-no-mitochondria-no-chloroplast.qza \
  --p-sampling-depth 1500 \
  --m-metadata-file 16s_noMockorPosorNeg_metadata.tsv \
  --output-dir core-metrics-results-rare1500

# After examining PCAs of each of the rarefaction depths, there are no visual disparities between. Thus, I plan to use the 1500 as the sampling depth to preserve more of the samples.
```

### ITS processing

##### Importing data into a Qiime2 artifact (.qza file)

Files were renamed to fit the Qiime2 EMPPairedEndSequences, for example forward.fastq.gz, reverse.fastq.gz, and barcode.fastq.gz

```bash
qiime tools import \
  --type EMPPairedEndSequences \
  --input-path /home/joel.swift/data/ITS/ \
  --output-path its-paired-end-sequences.qza
```

##### Initial processing of data

```bash
# I had to assign a new directory as temp to get Qiime demux to run

# Demultiplexing
qiime demux emp-paired \
  --m-barcodes-file its_metadata.tsv \
  --m-barcodes-column BarcodeSequence \
  --i-seqs its-paired-end-sequences.qza \
  --o-per-sample-sequences demux.qza \
  --p-rev-comp-mapping-barcodes
qiime demux summarize --i-data demux.qza --o-visualization demux.qzv

# Using grep searches of the forward and reverse reads using the ITS primers (F-CTTGGTCATTTAGAGGAAGTAA and R-GCTGCGTTCTTCATCGATGC) I found that many of the reads show signs that the reads were sequenced through. This means that the reverse completement of the primers above was found at the tail end of the sequnecing reads. 

# Cutadapt
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux.qza \
  --p-adapter-f GCATCGATGAAGAACGCAGC \
  --p-front-f CTTGGTCATTTAGAGGAAGTAA \
  --p-adapter-r TTACTTCCTCTAAATGACCAAG \
  --p-front-r GCTGCGTTCTTCATCGATGC \
  --o-trimmed-sequences demux-trimmed.qza
qiime demux summarize --i-data demux-trimmed.qza --o-visualization demux-trimmed.qzv

# Denoising and merging pair-end reads, reads were not truncted in order to preserve biologically relevent length varation. --p-trunc-q 2 was used to remove very small length sequences since there was no truncation used (in this case I had to specify a --p-trunc-len-l of 0 to get around the fact these settings are required in QIIME2).

# I ran multiple dada2 runs to figure out an appropraite MaxEE value as many sequences were being lost to filtering.

## Per B. Callahan on QIIME2 forum
## Where are sequences being lost
## If filtering: Try truncating sooner, especially before quality crashes, or raise maxEE.
## If merging: Make sure your reads still overlap (by 20nt + amplicon-length-variation) after truncation.
## If chimera removal: You probably need to remove primers from the reads

#Dada2 with MAXee at 3 after testing with mock samples are taxa bar graphs
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-trimmed.qza \
  --p-trim-left-f 12 \
  --p-trim-left-r 12 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-max-ee 3 \
  --p-trunc-q 2 \
  --o-table table-trimmed.qza \
  --o-representative-sequences rep-seqs-trimmed.qza \
  --o-denoising-stats denoising-stats-trimmed.qza \
  --verbose \
  --p-n-threads 0

#Visualize dada2 results
qiime metadata tabulate \
  --m-input-file denoising-stats-trimmed.qza \
  --o-visualization denoising-stats-trimmed.qzv

# Removing mock, positive, and negative control samples
qiime feature-table filter-samples \
  --i-table table-trimmed.qza \
  --m-metadata-file its_metadata_noControlsorWine.tsv \
  --o-filtered-table filtered-trimmed-table.qza

# Removing ASVs that with less than 0.1% mean sample depth
# Mean sample read depth (excluding controls and wine) = 30576
# Min frequency = 31
# Removed 1046 ASVs (features)
qiime feature-table filter-features --i-table filtered-trimmed-table.qza \
  --p-min-frequency 31 \
  --p-min-samples 1 \
  --o-filtered-table filtered-nocontrol-trimmed-table.qza
```

##### Taxonomy analyses

```bash
# Getting UNITE database
#https://doi.org/10.15156/BIO/786334

wget https://files.plutof.ut.ee/public/orig/51/6F/516F387FC543287E1E2B04BA4654443082FE3D7050E92F5D53BA0702E4E77CD4.zip

unzip 516F387FC543287E1E2B04BA4654443082FE3D7050E92F5D53BA0702E4E77CD4.zip

#remove lowercase
awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' developer/sh_refs_qiime_ver8_99_02.02.2019_dev.fasta > developer/sh_refs_qiime_ver8_99_02.02.2019_dev_uppercase.fasta

# Import unite reference sequences into qiime object
qiime tools import \
  --type FeatureData[Sequence] \
  --input-path developer/sh_refs_qiime_ver8_99_02.02.2019_dev_uppercase.fasta \
  --output-path unite-ver8-99-seqs-02.02.2019.qza

# Import the taxonomy of reference seqeunces into qiime object
qiime tools import \
 --type FeatureData[Taxonomy] \
 --input-path developer/sh_taxonomy_qiime_ver8_99_02.02.2019_dev.txt \
 --output-path unite-ver8-99-tax-02.02.2019.qza \
 --input-format HeaderlessTSVTaxonomyFormat

# Train naive bayes classifier
qiime feature-classifier fit-classifier-naive-bayes \
 --i-reference-reads unite-ver8-99-seqs-02.02.2019.qza \
 --i-reference-taxonomy unite-ver8-99-tax-02.02.2019.qza \
 --o-classifier unite-ver8-99-classifier-02.02.2019.qza

# Obtaining taxonomy for the its sequences
qiime feature-classifier classify-sklearn \
  --i-classifier unite-ver8-99-classifier-02.02.2019.qza \
  --i-reads rep-seqs-trimmed.qza \
  --o-classification taxonomy.qza

# Filtering out unassigned sequences  
qiime taxa filter-table \
  --i-table filtered-nocontrol-trimmed-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-include p__ \
  --o-filtered-table filtered-nocontrol-trimmed-table.qza

# Barplot
qiime taxa barplot \
  --i-table filtered-nocontrol-trimmed-table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file its_metadata_noControlsorWine.tsv \
  --o-visualization taxa-bar-plots.qzv
```
##### Diversity analyses

```bash
# Testing rarefaction sampling depths

# Generation of rarefaction curves
qiime diversity alpha-rarefaction \
  --i-table filtered-nocontrol-trimmed-table.qza \
  --p-max-depth 10000 \
  --p-steps 10 \
  --m-metadata-file its_metadata_noControlsorWine.tsv \
  --o-visualization rarefaction_curves.qzv

# Looking at rarefaction curves I see that the values do not quite plateau at 5000 sequences per sample but higher values removed more samples than I would like. Thus I choose to error on the side of retaining samples at the cost of including more sequences per sample.

# At 5000 rarefaction 11 samples were removed from the dataset
# Sample
# MV11.3309C.2-3.S
# MV14.3309C.2-3.R
# MV16.3309C.14-15.S
# MV15.1103P.10-11.R
# MV09.3309C.10.R
# MV09.OWN.3.R
# MV09.SO4.15.R
# MV10.SO4.7.R
# MV15.1103P.10-11.S
# MV12.SO4.10-11.S
# MV10.1103P.14.R

qiime diversity core-metrics \
  --i-table filtered-nocontrol-trimmed-table.qza \
  --p-sampling-depth 5000 \
  --m-metadata-file its_metadata_noControlsorWine.tsv \
  --output-dir core-metrics-results-rare5000
```