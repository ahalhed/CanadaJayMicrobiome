# Author: Alicia Halhed
# Species: Canada (Grey) Jay
# Sample: Oral Swabs
# 16S rRNA
# working with qiime2-2020.2

# script starts here
# ________________________________________

# working in a SHARCNET folder for the jays
cd /home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome
# request interactive session for testing
# salloc --mem-per-cpu=4G --account=def-cottenie --time=0-01:00:00
# data transfer to graham (not github - for large files)
# rsync -avz --no-g --no-p <source> <destination>
# load miniconda
module load miniconda3
# Activate QIIME2
conda activate qiime2-2020.8

# import the sequence data
# the raw files directory includes the "first jay" that Mason sent
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path input/rawFiles \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza

# below was a quick run, so ran in login node
# downloaded to local to look at the plot
qiime demux summarize \
  --i-data demux-paired-end.qza \
  --o-visualization demux-paired-end.qzv

# going to trim at 250 bp, don't really drop much based on plot
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux-paired-end.qza \
  --p-trim-left 0 \
  --p-trunc-len 250 \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-table table-dada2.qza \
  --o-denoising-stats stats-dada2.qza

qiime metadata tabulate \
  --m-input-file stats-dada2.qza \
  --o-visualization stats-dada2.qzv


# in the original metadata spreadsheet, KOOLTOSR-2017 and KOOLTOSR-2018 were labelled as G19
# LOYLOOSR-2017 was labelled as G21, but fell between the two G19's
# I am going to assume there was mislabelling here and keep KOOLTOSR-2017 as G19
# KOOLTOSR-2018 is moved to G21 (21st entry) and LOYLOOSR-2017 as G20 (20th entry)
# need to replace the degree symbol with ' for qiime to take the metadata
# tabulate the metadata
qiime metadata tabulate \
  --m-input-file input/jay-met.tsv \
  --o-visualization tabulated-metadata.qzv
# to look at a visualization
qiime tools view tabulated-metadata.qzv

# Going to work with the DADA2 results
# will use this figure to make sampling depth decision for
qiime feature-table summarize \
  --i-table table-dada2.qza \
  --o-visualization table-dada2.qzv \
  --m-sample-metadata-file input/jay-met.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-dada2.qza \
  --o-visualization rep-seqs-dada2.qzv

# drop the singletons from the OTU table
qiime feature-table filter-features \
  --i-table table-dada2.qza \
  --p-min-samples 2 \
  --o-filtered-table filtered-table-no-singletons.qza
# filter representative sequences for this group
qiime feature-table filter-seqs \
  --i-data rep-seqs-dada2.qza \
  --i-table filtered-table-no-singletons.qza \
  --o-filtered-data rep-seqs-no-singletons.qza

# create the phylogenetic tree
# fasttree pipeline
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-no-singletons.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree trees/unrooted-tree.qza \
  --o-rooted-tree trees/rooted-tree.qza

# ASSIGN TAXONOMY (need more than 4G ram - 8G works for GG)
# Accessing the pretrained classifiers from here
# https://docs.qiime2.org/2020.8/data-resources/
# Download the  Green genes reference database (smaller, runs faster)
wget -O "gg-13-8-99-nb-classifier.qza" "https://data.qiime2.org/2020.8/common/gg-13-8-99-nb-classifier.qza"
# assignment
qiime feature-classifier classify-sklearn \
 --i-classifier references/gg-13-8-99-nb-classifier.qza \
 --p-n-jobs 16 \
 --i-reads rep-seqs-no-singletons.qza \
 --o-classification taxonomy/GG-taxonomy.qza
# Generating taxonomy visualization
qiime taxa barplot \
  --i-table filtered-table-no-singletons.qza \
  --i-taxonomy taxonomy/GG-taxonomy.qza \
  --m-metadata-file input/jay-met.tsv \
  --o-visualization taxonomy/GG-taxa-bar-plots.qzv
# Extracting Taxonomic Clasification
# Phylum
qiime taxa collapse \
  --i-table filtered-table-no-singletons.qza \
  --i-taxonomy taxonomy/GG-taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table taxonomy/GG-table-l2.qza
# Class
qiime taxa collapse \
  --i-table filtered-table-no-singletons.qza \
  --i-taxonomy taxonomy/GG-taxonomy.qza \
  --p-level 3 \
  --o-collapsed-table taxonomy/GG-table-l3.qza
# Order
qiime taxa collapse \
  --i-table filtered-table-no-singletons.qza \
  --i-taxonomy taxonomy/GG-taxonomy.qza \
  --p-level 4 \
  --o-collapsed-table taxonomy/GG-table-l4.qza
# Family
qiime taxa collapse \
  --i-table filtered-table-no-singletons.qza \
  --i-taxonomy taxonomy/GG-taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table taxonomy/GG-table-l5.qza
# Genus
qiime taxa collapse \
  --i-table filtered-table-no-singletons.qza \
  --i-taxonomy taxonomy/GG-taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table taxonomy/GG-table-l6.qza
# Species
qiime taxa collapse \
  --i-table filtered-table-no-singletons.qza \
  --i-taxonomy taxonomy/GG-taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table taxonomy/GG-table-l7.qza

# Obtaining SILVA reference database (much larger database, will likely do a better job at classifying)
wget -O "silva-132-99-nb-classifier.qza" "https://data.qiime2.org/2020.8/common/silva-138-99-nb-classifier.qza"
# This has a high memory requirement (mem=128G,ntasks=16), but runs relatively quick (<30 min)
# Classifying taxonomies
qiime feature-classifier classify-sklearn \
  --i-classifier references/silva-132-99-nb-classifier.qza \
  --p-n-jobs 16 \
  --i-reads rep-seqs-no-singletons.qza \
  --o-classification taxonomy/SILVA-taxonomy.qza
# Generating taxonomy visualization
qiime taxa barplot \
  --i-table filtered-table-no-singletons.qza \
  --i-taxonomy taxonomy/SILVA-taxonomy.qza \
  --m-metadata-file input/jay-met.tsv \
  --o-visualization taxonomy/SILVA-dada2-taxa-bar-plots.qzv
# Extracting Taxonomic Clasification
# Phylum
qiime taxa collapse \
  --i-table filtered-table-no-singletons.qza \
  --i-taxonomy taxonomy/SILVA-taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table taxonomy/SILVA-table-l2.qza
# Class
qiime taxa collapse \
  --i-table filtered-table-no-singletons.qza \
  --i-taxonomy taxonomy/SILVA-taxonomy.qza \
  --p-level 3 \
  --o-collapsed-table taxonomy/SILVA-table-l3.qza
# Order
qiime taxa collapse \
  --i-table filtered-table-no-singletons.qza \
  --i-taxonomy taxonomy/SILVA-taxonomy.qza \
  --p-level 4 \
  --o-collapsed-table taxonomy/SILVA-table-l4.qza
# Family
qiime taxa collapse \
  --i-table filtered-table-no-singletons.qza \
  --i-taxonomy taxonomy/SILVA-taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table taxonomy/SILVA-table-l5.qza
# Genus
qiime taxa collapse \
  --i-table filtered-table-no-singletons.qza \
  --i-taxonomy taxonomy/SILVA-taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table taxonomy/SILVA-table-l6.qza
# Species
qiime taxa collapse \
  --i-table filtered-table-no-singletons.qza \
  --i-taxonomy taxonomy/SILVA-taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table taxonomy/SILVA-table-l7.qza

# drop mitochondrial and chloroplast sequences
# these aren't necessarily informative of bacterial diversity
qiime taxa filter-table \
  --i-table filtered-table-no-singletons.qza \
  --i-taxonomy taxonomy/SILVA-taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table filtered-table-no-singletons-mitochondria-chloroplast.qza

# compositional data analysis
# compute aitchison distance matrix
qiime deicode rpca \
    --i-table filtered-table-no-singletons-mitochondria-chloroplast.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot aitchison-ordination.qza \
    --o-distance-matrix aitchison-distance.qza

qiime emperor biplot \
    --i-biplot aitchison-ordination.qza \
    --m-sample-metadata-file input/jay-met.tsv \
    --m-feature-metadata-file taxonomy/SILVA-taxonomy.qza \
    --o-visualization aitchison-biplot.qzv \
    --p-number-of-features 8


# multivariate comparison
# PERMANOVA & ANOSIM
# by Territory
qiime diversity beta-group-significance \
  --i-distance-matrix aitchison-distance.qza \
  --m-metadata-file input/jay-met.tsv  \
  --m-metadata-column Territory \
  --o-visualization beta/Territory-permanova-aitchison.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix aitchison-distance.qza \
  --m-metadata-file input/jay-met.tsv  \
  --m-metadata-column Territory \
  --p-method anosim \
  --o-visualization beta/Territory-anosim-aitchison.qzv

# by TerritoryQuality
qiime diversity beta-group-significance \
  --i-distance-matrix aitchison-distance.qza \
  --m-metadata-file input/jay-met.tsv \
  --m-metadata-column TerritoryQuality \
  --o-visualization beta/TerritoryQuality-permanova-aitchison.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix aitchison-distance.qza \
  --m-metadata-file input/jay-met.tsv \
  --m-metadata-column TerritoryQuality \
  --p-method anosim \
  --o-visualization beta/TerritoryQuality-anosim-aitchison.qzv

# by Food supplement
qiime diversity beta-group-significance \
  --i-distance-matrix aitchison-distance.qza \
  --m-metadata-file input/jay-met.tsv \
  --m-metadata-column FoodSupplement \
  --o-visualization beta/FoodSupplement-permanova-aitchison.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix aitchison-distance.qza \
  --m-metadata-file input/jay-met.tsv \
  --m-metadata-column FoodSupplement \
  --p-method anosim \
  --o-visualization beta/FoodSupplement-anosim-aitchison.qzv

# by jayID
qiime diversity beta-group-significance \
  --i-distance-matrix aitchison-distance.qza \
  --m-metadata-file input/jay-met.tsv \
  --m-metadata-column JayID \
  --o-visualization beta/JayID-permanova-aitchison.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix aitchison-distance.qza \
  --m-metadata-file input/jay-met.tsv \
  --m-metadata-column JayID \
  --p-method anosim \
  --o-visualization beta/JayID-anosim-aitchison.qzv

# by CollectionSeason
qiime diversity beta-group-significance \
  --i-distance-matrix aitchison-distance.qza \
  --m-metadata-file input/jay-met.tsv \
  --m-metadata-column CollectionSeason \
  --o-visualization beta/CollectionSeason-permanova-aitchison.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix aitchison-distance.qza \
  --m-metadata-file input/jay-met.tsv \
  --m-metadata-column CollectionSeason \
  --p-method anosim \
  --o-visualization beta/CollectionSeason-anosim-aitchison.qzv

# by MeanTempC (average temperature on day of collection in degrees Celcius)
# doesn't work b/c is numeric, not categorical

# bioenv
# long running (>2 hours)
qiime diversity bioenv \
  --i-distance-matrix aitchison-distance.qza \
  --m-metadata-file input/jay-met.tsv \
  --o-visualization bioenv-aitchison.qzv

# correlational analysis
# with SparCC
# attempting installation with: https://github.com/shafferm/q2-SCNIC
qiime SCNIC sparcc-filter \
  --i-table filtered-table-no-singletons-mitochondria-chloroplast.qza \
  --o-table-filtered filtered-table-SparCC.qza
# build a correlation table
qiime SCNIC calculate-correlations \
  --i-table filtered-table-SparCC.qza \
  --p-method sparcc \
  --o-correlation-table SparCC-correls.qza
# build the correlation table
qiime SCNIC build-correlation-network-r \
  --i-correlation-table SparCC-correls.qza \
  --p-min-val .35 \
  --o-correlation-network SparCC-net.qza
# summmarize results
qiime SCNIC make-modules-on-correlation-table \
  --i-correlation-table SparCC-correls.qza \
  --i-feature-table filtered-table-no-singletons-mitochondria-chloroplast.qza \
  --p-min-r .35 \
  --o-collapsed-table fake_data.collapsed.qza \
  --o-correlation-network SparCC-net-modules.qza \
  --o-module-membership SparCC-membership.qza
# visualize results
qiime meta tabulate \
  --m-input-file SparCC-membership.qza \
  --o-visualization SparCC-membership.qzv


# differential abundance with ancom
qiime composition add-pseudocount \
  --i-table filtered-table-no-singletons-mitochondria-chloroplast.qza \
  --o-composition-table filtered-table-pseudo.qza

# need to make sure all rows have values
# by Territory
qiime composition ancom \
  --i-table filtered-table-pseudo.qza \
  --m-metadata-file input/jay-met.tsv \
  --m-metadata-column Territory \
  --o-visualization ancom/Territory-ancom.qzv

# by TerritoryQuality
# unknown = U
qiime composition ancom \
  --i-table filtered-table-pseudo.qza \
  --m-metadata-file input/jay-met.tsv \
  --m-metadata-column TerritoryQuality \
  --o-visualization ancom/TerritoryQuality-ancom.qzv

# by FoodSupplement
# unknown = U
qiime composition ancom \
  --i-table filtered-table-pseudo.qza \
  --m-metadata-file input/jay-met.tsv \
  --m-metadata-column FoodSupplement \
  --o-visualization ancom/FoodSupplement-ancom.qzv

# by JayID
qiime composition ancom \
  --i-table filtered-table-pseudo.qza \
  --m-metadata-file input/jay-met.tsv \
  --m-metadata-column JayID \
  --o-visualization ancom/JayID-ancom.qzv

# by CollectionSeason
qiime composition ancom \
  --i-table filtered-table-pseudo.qza \
  --m-metadata-file input/jay-met.tsv \
  --m-metadata-column CollectionSeason \
  --o-visualization ancom/CollectionSeason-ancom.qzv

# by MeanTempC
# doesn't work - is numeric (must be categorical)
# could try forcing to categorical if needed

# ALDEx2 - not necessary because of ANCOM
# not working because of package/R version issues
# https://library.qiime2.org/plugins/q2-aldex2/24/

# Close QIIME2
conda deactivate