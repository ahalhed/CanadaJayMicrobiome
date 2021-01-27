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
conda activate qiime2-2020.8

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
  "https://data.qiime2.org/2018.4/common/silva-128-99-seqs.qza"
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

# create the phylogenetic tree
# fragment insertion tree (best for tree building on different 16S regions)
# have to use different versions of SILVA (128 here and 138 for OTUs) since 
# the QIIME2 version (2020.8) only has the recent (138) SILVA for OTU picking
# but the plugin (fragment-insertion) only has the older (128) SILVA version
wget -O "references/sepp-refs-silva-128.qza" \
  "https://data.qiime2.org/2020.8/common/sepp-refs-silva-128.qza"

qiime fragment-insertion sepp \
  --i-representative-sequences rep-seqs-no-singletons.qza \
  --i-reference-database references/sepp-refs-silva-128.qza \
  --o-tree trees/insertion-tree.qza \
  --o-placements trees/insertion-placements.qza
#Plugin error from fragment-insertion:
#Command '['run-sepp.sh', '/var/folders/wd/s68yr8892s137xwylsfl43nw0000gn/T/qiime2-archive-w__w93y1/540cc42f-8cdd-4386-ae0c-a0cc195c848f/data/dna-sequences.fasta', 'q2-fragment-insertion', '-x', '1', '-A', '1000', '-P', '5000', '-a', '/var/folders/wd/s68yr8892s137xwylsfl43nw0000gn/T/qiime2-archive-ixa4qg0f/e44b5e78-31e5-4a0f-9041-494bc3ca2df2/data/aligned-dna-sequences.fasta', '-t', '/var/folders/wd/s68yr8892s137xwylsfl43nw0000gn/T/qiime2-archive-ixa4qg0f/e44b5e78-31e5-4a0f-9041-494bc3ca2df2/data/tree.nwk', '-r', '/var/folders/wd/s68yr8892s137xwylsfl43nw0000gn/T/qiime2-archive-ixa4qg0f/e44b5e78-31e5-4a0f-9041-494bc3ca2df2/data/raxml-info.txt']' returned non-zero exit status 1.
# Debug info has been saved to /var/folders/wd/s68yr8892s137xwylsfl43nw0000gn/T/qiime2-q2cli-err-z1hbwmce.log

# filter otu table based on tree
qiime fragment-insertion filter-features \
  --i-table table.qza \
  --i-tree insertion-tree.qza \
  --o-filtered-table filtered-tree-table.qza \
  --o-removed-table removed-table.qza

# fasttree pipeline (denovo - will remove if the insertion works)
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-no-singletons.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree trees/unrooted-tree.qza \
  --o-rooted-tree trees/rooted-tree.qza


# ASSIGN TAXONOMY (need more than 4G ram)
# Accessing the pretrained classifiers from here
# https://docs.qiime2.org/2020.8/data-resources/
# Obtaining SILVA reference database 
# need this to do closed reference otu picking
wget -O "references/silva-138-99-nb-classifier.qza" \
  "https://data.qiime2.org/2020.8/common/silva-138-99-nb-classifier.qza"
# if switching to 2020.11, use below link
#https://data.qiime2.org/2020.11/common/silva-138-99-nb-classifier.qza
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

# filtering out blanks (for H1 & H2)
qiime feature-table filter-samples \
    --i-table filtered-table-no-singletons-mitochondria-chloroplast.qza \
    --m-metadata-file input/jay-met.tsv \
    --p-exclude-ids 'TRUE' \
    --p-where "[JayID]='BLANK'" \
    --o-filtered-table filtered-table-no-blanks.qza

# Hypothesis 1 - rarefaction
# rarefying to feed into core definition in R (needed nowhere else)
# 344 retains all samples
qiime feature-table rarefy \
    --i-table filtered-table-no-blanks.qza \
    --p-sampling-depth 344 \
    --p-no-with-replacement \
    --o-rarefied-table rarefied-table
# ran gj-core.R here, then produced a list of core and rare features
# filter for core variants
qiime feature-table filter-features \
  --i-table filtered-table-no-blanks.qza \
  --m-metadata-file CanadaJayMicrobiome/data/coreFeatures.tsv \
  --o-filtered-table filtered-table-core.qza
# filter for rare variants
qiime feature-table filter-features \
  --i-table filtered-table-no-blanks.qza \
  --m-metadata-file CanadaJayMicrobiome/data/coreFeatures.tsv \
  --p-exclude-ids TRUE \
  --o-filtered-table filtered-table-rare.qza
# Hypothesis 1 & 2 - full dataset
# compute aitchison distance matrix 
qiime deicode rpca \
    --i-table filtered-table-no-blanks.qza \
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
# repeat above for core variants
qiime deicode rpca \
    --i-table filtered-table-core.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot aitchison-ordination-core.qza \
    --o-distance-matrix aitchison-distance-core.qza

qiime emperor biplot \
    --i-biplot aitchison-ordination-core.qza \
    --m-sample-metadata-file input/jay-met.tsv \
    --m-feature-metadata-file taxonomy/SILVA-taxonomy.qza \
    --o-visualization aitchison-biplot-core.qzv \
    --p-number-of-features 8
# repeat above for rare variants
qiime deicode rpca \
    --i-table filtered-table-rare.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot aitchison-ordination-rare.qza \
    --o-distance-matrix aitchison-distance-rare.qza

qiime emperor biplot \
    --i-biplot aitchison-ordination-rare.qza \
    --m-sample-metadata-file input/jay-met.tsv \
    --m-feature-metadata-file taxonomy/SILVA-taxonomy.qza \
    --o-visualization aitchison-biplot-rare.qzv \
    --p-number-of-features 8
# Hypothesis 3 - spring 2020 samples (nest groups)
# using this i-table since we're filtering, so blanks will be dropped anyways
qiime feature-table filter-samples \
  --i-table filtered-table-no-singletons-mitochondria-chloroplast.qza \
  --m-metadata-file input/jay-met.tsv \
  --p-where "[CollectionSeason]='Spring' AND [CollectionYear]='2020'" \
  --o-filtered-table H3-filtered-table.qza

qiime deicode rpca \
    --i-table H3-filtered-table.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot H3-aitchison-ordination.qza \
    --o-distance-matrix H3-aitchison-distance.qza
# creates a sample tree
qiime diversity beta-rarefaction \
    --i-table H3-filtered-table.qza \
    --i-phylogeny trees/rooted-tree.qza \
    --p-metric 'aitchison' \
    --p-clustering-method 'nj' \
    --m-metadata-file input/jay-met.tsv \
    --p-sampling-depth 500 \
    --p-iterations 100 \
    --p-correlation-method 'spearman' \
    --p-color-scheme 'PuOr_r' \
    --o-visualization H3-aitchison-beta-rarefaction

# Hypothesis 4 - only samples with origin data
# H4-samples.tsv is a list of sampleid's to keep
qiime feature-table filter-samples \
  --i-table filtered-table-no-singletons-mitochondria-chloroplast.qza \
  --m-metadata-file CanadaJayMicrobiome/data/H4-samples.tsv \
  --o-filtered-table H4-filtered-table.qza

qiime deicode rpca \
    --i-table H4-filtered-table.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot H4-aitchison-ordination.qza \
    --o-distance-matrix H4-aitchison-distance.qza

# Close QIIME2
conda deactivate