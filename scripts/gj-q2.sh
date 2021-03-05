# Author: Alicia Halhed
# Species: Canada (Grey) Jay
# Sample: Oral Swabs
# 16S rRNA
# working with qiime2-2020.8

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
conda activate qiime2-2020.11

# import the sequence data
# the raw files directory includes sequences from all samples
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path input/rawFiles \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza

qiime demux summarize \
  --i-data demux-paired-end.qza \
  --o-visualization demux-paired-end.qzv


# going to trim at 250 bp
# this took just over two hours to run
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

# Obtaining SILVA reference database 
# need this to do closed reference otu picking
wget -O "references/silva-138-99-seqs.qza" \
  "https://data.qiime2.org/2020.11/common/silva-138-99-seqs.qza"
# switch 2020.8 to 2020.11 if updating q2 version

# using closed reference clustering to account for the two different runs
qiime vsearch cluster-features-closed-reference \
  --i-table table-dada2.qza \
  --i-sequences rep-seqs-dada2.qza \
  --i-reference-sequences references/silva-138-99-seqs.qza \
  --p-perc-identity 0.99 \
  --o-clustered-table table-cr-99.qza \
  --o-clustered-sequences rep-seqs-cr-99.qza \
  --o-unmatched-sequences unmatched-cr-99.qza

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

# Going to work with the closed reference results
# will use this figure to make sampling depth decision (for core selection only)
qiime feature-table summarize \
  --i-table table-cr-99.qza \
  --o-visualization table-cr-99.qzv \
  --m-sample-metadata-file input/jay-met.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-cr-99.qza \
  --o-visualization rep-seqs-cr-99.qzv

# drop the singletons from the OTU tables
qiime feature-table filter-features \
  --i-table table-cr-99.qza \
  --p-min-samples 2 \
  --o-filtered-table filtered-table-no-singletons.qza

# filter representative sequences for all samples
qiime feature-table filter-seqs \
  --i-data rep-seqs-cr-99.qza \
  --i-table filtered-table-no-singletons.qza \
  --o-filtered-data rep-seqs-no-singletons.qza


# ASSIGN TAXONOMY (need more than 4G ram)
# Accessing the pretrained classifiers from here
# https://docs.qiime2.org/2020.8/data-resources/
# Obtaining SILVA reference database 
# need this to do closed reference otu picking
wget -O "references/silva-138-99-nb-classifier.qza" \
  "https://data.qiime2.org/2020.11/common/silva-138-99-nb-classifier.qza"
# classifier ref
# Bokulich, N.A., Robeson, M., Dillon, M.R. bokulich-lab/RESCRIPt. Zenodo. http://doi.org/10.5281/zenodo.3891931
# Bokulich, N.A., Kaehler, B.D., Rideout, J.R. et al. Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2’s q2-feature-classifier plugin. Microbiome 6, 90 (2018). https://doi.org/10.1186/s40168-018-0470-z

# This has a high memory requirement (mem=128G,ntasks=16), but runs relatively quick (<30 min)
# Classifying taxonomies (see 5taxonomyCR for closed reference)
qiime feature-classifier classify-sklearn \
  --i-classifier references/silva-138-99-nb-classifier.qza \
  --p-n-jobs 16 \
  --i-reads rep-seqs-no-singletons.qza \
  --o-classification taxonomy/SILVA-taxonomy.qza
# Generating taxonomy visualization
qiime taxa barplot \
  --i-table filtered-table-no-singletons.qza \
  --i-taxonomy taxonomy/SILVA-taxonomy.qza \
  --m-metadata-file input/jay-met.tsv \
  --o-visualization taxonomy/SILVA-cr99-taxa-bar-plots.qzv
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
qiime taxa filter-seqs \
  --i-sequences rep-seqs-no-singletons.qza \
  --i-taxonomy taxonomy/SILVA-taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-sequences rep-seqs-no-singletons-mitochondria-chloroplast.qza

# make phylogenetic tree (wanted to do insertion but ref was dif version so nah)
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-no-singletons-mitochondria-chloroplast.qza \
  --p-n-threads 'auto' \
  --output-dir trees

# filtering out blanks
qiime feature-table filter-samples \
    --i-table filtered-table-no-singletons-mitochondria-chloroplast.qza \
    --m-metadata-file input/jay-met.tsv \
    --p-exclude-ids 'TRUE' \
    --p-where "[JayID]='BLANK'" \
    --o-filtered-table filtered-table-no-blanks.qza

# rarefying to feed into core definition in R (needed nowhere else)
# 344 retains all samples
qiime feature-table rarefy \
    --i-table filtered-table-no-blanks.qza \
    --p-sampling-depth 344 \
    --p-no-with-replacement \
    --o-rarefied-table rarefied-table
# ran gj-core.R here to produce a list of core and rare features

# Hypothesis 1
# Prediction 1A+B - only breeders without supplementation
qiime feature-table filter-samples \
  --i-table filtered-table-no-blanks.qza \
  --m-metadata-file input/jay-met.tsv \
  --p-where "[BreedingStatus]='Breeder' AND [FoodSupplement]='N'" \
  --o-filtered-table P1AB-filtered-table.qza

# Prediction 1A - territory groups (all breeding statuses, without food supplementation)
qiime deicode rpca \
    --i-table P1AB-filtered-table.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot P1A-aitchison-ordination.qza \
    --o-distance-matrix P1A-aitchison-distance.qza
  
# Prediction 1C - only breeders without supplementation and with territory quality information
qiime feature-table filter-samples \
  --i-table P1AB-filtered-table.qza \
  --m-metadata-file input/jay-met.tsv \
  --p-where "[TerritoryQuality] IN ('H', 'M', 'L')" \
  --o-filtered-table P1C-filtered-table.qza
qiime deicode rpca \
    --i-table P1C-filtered-table.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot P1C-aitchison-ordination.qza \
    --o-distance-matrix P1C-aitchison-distance.qza

# differential abundance testing
qiime gneiss gradient-clustering \
  --i-table P1C-filtered-table.qza \
  --m-gradient-file input/jay-met.tsv \
  --m-gradient-column ProportionSpruceOnTerritory \
  --o-clustering P1C-gradient-hierarchy.qza
qiime gneiss dendrogram-heatmap \
  --i-table P1C-filtered-table.qza \
  --i-tree P1C-gradient-hierarchy.qza \
  --m-metadata-file input/jay-met.tsv \
  --m-metadata-column TerritoryQuality \
  --p-color-map viridis \
  --o-visualization P1C-heatmap.qzv
  

# Hypothesis 2
# Prediction 2A - host associated factors (all samples)
qiime deicode rpca \
    --i-table filtered-table-no-blanks.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot P2A-aitchison-ordination.qza \
    --o-distance-matrix P2A-aitchison-distance.qza


# Hypothesis 3
# Prediction 3A  - food caching (winter and spring only)
qiime feature-table filter-samples \
  --i-table filtered-table-no-blanks.qza \
  --m-metadata-file input/jay-met.tsv \
  --p-where "[CollectionSeason] IN ('Winter', 'Spring')" \
  --o-filtered-table P3A-filtered-table.qza

# Prediction 3B/C - food supplementation
qiime feature-table filter-samples \
  --i-table filtered-table-no-blanks.qza \
  --m-metadata-file input/jay-met.tsv \
  --p-where "[CollectionYear] IN ('2017', '2018') AND [BreedingStatus]='Breeder'" \
  --o-filtered-table P3BC-filtered-table.qza
# Prediction 3D - food supplementation
qiime deicode rpca \
    --i-table P3BC-filtered-table.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot P3C-aitchison-ordination.qza \
    --o-distance-matrix P3C-aitchison-distance.qza

# Hypothesis 4 - parental care
# all predictions
qiime feature-table filter-samples \
  --i-table filtered-table-no-blanks.qza \
  --m-metadata-file input/jay-met.tsv \
  --p-where "[CollectionSeason]='Spring' AND [CollectionYear]='2020'" \
  --o-filtered-table P4ABCD-filtered-table.qza

qiime deicode rpca \
    --i-table P4ABCD-filtered-table.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot P4ABCD-aitchison-ordination.qza \
    --o-distance-matrix P4ABCD-aitchison-distance.qza

# Hypothesis 5 - only samples with origin data
# H5-samples.tsv is a list of sampleid's to keep
# no nestlings from S20 were resampled in F20
qiime feature-table filter-samples \
  --i-table filtered-table-no-blanks.qza \
  --m-metadata-file CanadaJayMicrobiome/data/H5-samples.tsv \
  --o-filtered-table H5-filtered-table.qza

# P5A - breeders with origin information
qiime feature-table filter-samples \
  --i-table H5-filtered-table.qza \
  --m-metadata-file input/jay-met.tsv \
  --p-where "[BreedingStatus]='Breeder'" \
  --o-filtered-table P5A-filtered-table.qza

# P5B - non-breeders with origin information
qiime feature-table filter-samples \
  --i-table H5-filtered-table.qza \
  --m-metadata-file input/jay-met.tsv \
  --p-where "[BreedingStatus]='Non-breeder'" \
  --o-filtered-table P5B-filtered-table.qza


# Close QIIME2
conda deactivate