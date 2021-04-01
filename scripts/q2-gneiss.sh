# Author: Alicia Halhed
# Species: Canada (Grey) Jay
# Sample: Oral Swabs
# 16S rRNA
# working with qiime2-2020.11

# script starts here
# ________________________________________

# differential abundance testing with gneiss
# load env on graham
module load nixpkgs/16.09 miniconda3
conda activate qiime2-2020.11

# make folder
mkdir differential
# move into folder
cd differential

# Hypothesis 1
# Prediction 1A - Territories
# fall 2017
qiime feature-table filter-samples \
  --i-table ../P1AB-filtered-table.qza \
  --m-metadata-file ../input/jay-met.tsv \
  --p-where "[CollectionYear]='2017' AND [CollectionSeason]='Fall'" \
  --o-filtered-table P1A-table-F17.qza

qiime gneiss correlation-clustering \
  --i-table P1A-table-F17.qza \
  --o-clustering P1A-hierarchy-F17.qza

qiime gneiss dendrogram-heatmap \
  --i-table P1A-table-F17.qza  \
  --i-tree P1A-hierarchy-F17.qza \
  --m-metadata-file ../input/jay-met.tsv  \
  --m-metadata-column Territory \
  --p-color-map viridis \
  --o-visualization P1A-heatmap-F17.qzv

# fall 2018
qiime feature-table filter-samples \
  --i-table ../P1AB-filtered-table.qza \
  --m-metadata-file ../input/jay-met.tsv \
  --p-where "[CollectionYear]='2018' AND [CollectionSeason]='Fall'" \
  --o-filtered-table P1A-table-F18.qza

qiime gneiss correlation-clustering \
  --i-table P1A-table-F18.qza \
  --o-clustering P1A-hierarchy-F18.qza

qiime gneiss dendrogram-heatmap \
  --i-table P1A-table-F18.qza  \
  --i-tree P1A-hierarchy-F18.qza \
  --m-metadata-file ../input/jay-met.tsv  \
  --m-metadata-column Territory \
  --p-color-map viridis \
  --o-visualization P1A-heatmap-F18.qzv

# fall 2020
qiime feature-table filter-samples \
  --i-table ../P1AB-filtered-table.qza \
  --m-metadata-file ../input/jay-met.tsv \
  --p-where "[CollectionYear]='2020' AND [CollectionSeason]='Fall'" \
  --o-filtered-table P1A-table-F20.qza

qiime gneiss correlation-clustering \
  --i-table P1A-table-F20.qza \
  --o-clustering P1A-hierarchy-F20.qza

qiime gneiss dendrogram-heatmap \
  --i-table P1A-table-F20.qza  \
  --i-tree P1A-hierarchy-F20.qza \
  --m-metadata-file ../input/jay-met.tsv  \
  --m-metadata-column Territory \
  --p-color-map viridis \
  --o-visualization P1A-heatmap-F20.qzv

# spring 2020
qiime feature-table filter-samples \
  --i-table ../P1AB-filtered-table.qza \
  --m-metadata-file ../input/jay-met.tsv \
  --p-where "[CollectionYear]='2020' AND [CollectionSeason]='Spring'" \
  --o-filtered-table P1A-table-S20.qza

qiime gneiss correlation-clustering \
  --i-table P1A-table-S20.qza \
  --o-clustering P1A-hierarchy-S20.qza

qiime gneiss dendrogram-heatmap \
  --i-table P1A-table-S20.qza  \
  --i-tree P1A-hierarchy-S20.qza \
  --m-metadata-file ../input/jay-met.tsv  \
  --m-metadata-column Territory \
  --p-color-map viridis \
  --o-visualization P1A-heatmap-S20.qzv

# Prediction 1C - only breeders without supplementation and with territory quality information
# fall 2017
qiime feature-table filter-samples \
  --i-table ../P1C-filtered-table.qza \
  --m-metadata-file ../input/jay-met.tsv \
  --p-where "[CollectionYear]='2017' AND [CollectionSeason]='Fall'" \
  --o-filtered-table P1C-table-F17.qza

qiime gneiss correlation-clustering \
  --i-table P1C-table-F17.qza \
  --o-clustering P1C-hierarchy-F17.qza

qiime gneiss dendrogram-heatmap \
  --i-table P1C-table-F17.qza  \
  --i-tree P1C-hierarchy-F17.qza \
  --m-metadata-file ../input/jay-met.tsv  \
  --m-metadata-column Territory \
  --p-color-map viridis \
  --o-visualization P1C-heatmap-F17.qzv

# fall 2018
qiime feature-table filter-samples \
  --i-table ../P1C-filtered-table.qza \
  --m-metadata-file ../input/jay-met.tsv \
  --p-where "[CollectionYear]='2018' AND [CollectionSeason]='Fall'" \
  --o-filtered-table P1C-table-F18.qza

qiime gneiss correlation-clustering \
  --i-table P1C-table-F18.qza \
  --o-clustering P1C-hierarchy-F18.qza

qiime gneiss dendrogram-heatmap \
  --i-table P1C-table-F18.qza  \
  --i-tree P1C-hierarchy-F18.qza \
  --m-metadata-file ../input/jay-met.tsv  \
  --m-metadata-column Territory \
  --p-color-map viridis \
  --o-visualization P1C-heatmap-F18.qzv

# fall 2020
qiime feature-table filter-samples \
  --i-table ../P1C-filtered-table.qza \
  --m-metadata-file ../input/jay-met.tsv \
  --p-where "[CollectionYear]='2020' AND [CollectionSeason]='Fall'" \
  --o-filtered-table P1C-table-F20.qza

qiime gneiss correlation-clustering \
  --i-table P1C-table-F20.qza \
  --o-clustering P1C-hierarchy-F20.qza

qiime gneiss dendrogram-heatmap \
  --i-table P1C-table-F20.qza  \
  --i-tree P1C-hierarchy-F20.qza \
  --m-metadata-file ../input/jay-met.tsv  \
  --m-metadata-column Territory \
  --p-color-map viridis \
  --o-visualization P1C-heatmap-F20.qzv

# spring 2020
qiime feature-table filter-samples \
  --i-table ../P1C-filtered-table.qza \
  --m-metadata-file ../input/jay-met.tsv \
  --p-where "[CollectionYear]='2020' AND [CollectionSeason]='Spring'" \
  --o-filtered-table P1C-table-S20.qza

qiime gneiss correlation-clustering \
  --i-table P1C-table-S20.qza \
  --o-clustering P1C-hierarchy-S20.qza

qiime gneiss dendrogram-heatmap \
  --i-table P1C-table-S20.qza  \
  --i-tree P1C-hierarchy-S20.qza \
  --m-metadata-file ../input/jay-met.tsv  \
  --m-metadata-column Territory \
  --p-color-map viridis \
  --o-visualization P1C-heatmap-S20.qzv

# Hypothesis 2 - host factors
# fall 2017
qiime feature-table filter-samples \
  --i-table ../filtered-table-no-blanks.qza \
  --m-metadata-file ../input/jay-met.tsv \
  --p-where "[CollectionYear]='2017' AND [CollectionSeason]='Fall'" \
  --o-filtered-table P2A-table-F17.qza

qiime gneiss correlation-clustering \
  --i-table P2A-table-F17.qza \
  --o-clustering P2A-hierarchy-F17.qza

qiime gneiss dendrogram-heatmap \
  --i-table P2A-table-F17.qza  \
  --i-tree P2A-hierarchy-F17.qza \
  --m-metadata-file ../input/jay-met.tsv  \
  --m-metadata-column Sex \
  --p-color-map viridis \
  --o-visualization P2A-heatmap-F17.qzv

# fall 2018
qiime feature-table filter-samples \
  --i-table ../filtered-table-no-blanks.qza \
  --m-metadata-file ../input/jay-met.tsv \
  --p-where "[CollectionYear]='2018' AND [CollectionSeason]='Fall'" \
  --o-filtered-table P2A-table-F18.qza

qiime gneiss correlation-clustering \
  --i-table P2A-table-F18.qza \
  --o-clustering P2A-hierarchy-F18.qza

qiime gneiss dendrogram-heatmap \
  --i-table P2A-table-F18.qza  \
  --i-tree P2A-hierarchy-F18.qza \
  --m-metadata-file ../input/jay-met.tsv  \
  --m-metadata-column Sex \
  --p-color-map viridis \
  --o-visualization P2A-heatmap-F18.qzv

# fall 2020
qiime feature-table filter-samples \
  --i-table ../filtered-table-no-blanks.qza \
  --m-metadata-file ../input/jay-met.tsv \
  --p-where "[CollectionYear]='2020' AND [CollectionSeason]='Fall'" \
  --o-filtered-table P2A-table-F20.qza

qiime gneiss correlation-clustering \
  --i-table P2A-table-F20.qza \
  --o-clustering P2A-hierarchy-F20.qza

qiime gneiss dendrogram-heatmap \
  --i-table P2A-table-F20.qza  \
  --i-tree P2A-hierarchy-F20.qza \
  --m-metadata-file ../input/jay-met.tsv  \
  --m-metadata-column Territory \
  --p-color-map viridis \
  --o-visualization P2A-heatmap-F20.qzv

# spring 2020
qiime feature-table filter-samples \
  --i-table ../filtered-table-no-blanks.qza \
  --m-metadata-file ../input/jay-met.tsv \
  --p-where "[CollectionYear]='2020' AND [CollectionSeason]='Spring'" \
  --o-filtered-table P2A-table-S20.qza

qiime gneiss correlation-clustering \
  --i-table P2A-table-S20.qza \
  --o-clustering P2A-hierarchy-S20.qza

qiime gneiss dendrogram-heatmap \
  --i-table P2A-table-S20.qza  \
  --i-tree P2A-hierarchy-S20.qza \
  --m-metadata-file ../input/jay-met.tsv  \
  --m-metadata-column Territory \
  --p-color-map viridis \
  --o-visualization P2A-heatmap-S20.qzv

# Hypothesis 3 - diet
# fall 2017
qiime feature-table filter-samples \
  --i-table ../P3C-filtered-table.qza \
  --m-metadata-file ../input/jay-met.tsv \
  --p-where "[CollectionYear]='2017' AND [CollectionSeason]='Fall'" \
  --o-filtered-table P3C-table-F17.qza

qiime gneiss correlation-clustering \
  --i-table P3C-table-F17.qza \
  --o-clustering P3C-hierarchy-F17.qza
# need to run figure
qiime gneiss dendrogram-heatmap \
  --i-table P3C-table-F17.qza  \
  --i-tree P3C-hierarchy-F17.qza \
  --m-metadata-file ../input/jay-met.tsv  \
  --m-metadata-column FoodSupplementation \
  --p-color-map viridis \
  --o-visualization P3C-heatmap-F17.qzv

# fall 2018
qiime feature-table filter-samples \
  --i-table ../P3C-filtered-table.qza \
  --m-metadata-file ../input/jay-met.tsv \
  --p-where "[CollectionYear]='2018' AND [CollectionSeason]='Fall'" \
  --o-filtered-table P3C-table-F18.qza

qiime gneiss correlation-clustering \
  --i-table P3C-table-F18.qza \
  --o-clustering P3C-hierarchy-F18.qza
# need to run figure
qiime gneiss dendrogram-heatmap \
  --i-table P3C-table-F18.qza  \
  --i-tree P3C-hierarchy-F18.qza \
  --m-metadata-file ../input/jay-met.tsv  \
  --m-metadata-column FoodSupplementation \
  --p-color-map viridis \
  --o-visualization P3C-heatmap-F18.qzv

# Hypothesis 4 - parental care
qiime gneiss correlation-clustering \
  --i-table ../H4-filtered-table.qza \
  --o-clustering H4-hierarchy.qza

qiime gneiss dendrogram-heatmap \
  --i-table ../H4-filtered-table.qza  \
  --i-tree H4-hierarchy.qza \
  --m-metadata-file ../input/jay-met.tsv  \
  --m-metadata-column BreedingStatus \
  --p-color-map viridis \
  --o-visualization H4-heatmap.qzv

# unload qiime environment
conda deactivate
