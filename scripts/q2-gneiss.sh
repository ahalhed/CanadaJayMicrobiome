# Author: Alicia Halhed
# Species: Canada (Grey) Jay
# Sample: Oral Swabs
# 16S rRNA
# working with qiime2-2020.11

# script starts here
# ________________________________________

# differential abundance testing with gneiss
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
  --i-table ../P1C-filtered-table.qza \
  --m-metadata-file ../input/jay-met.tsv \
  --p-where "[CollectionYear]='2018' AND [CollectionSeason]='Fall'" \
  --o-filtered-table P1C-table-F18.qza

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
  --i-table ../P1C-filtered-table.qza \
  --m-metadata-file ../input/jay-met.tsv \
  --p-where "[CollectionYear]='2020' AND [CollectionSeason]='Fall'" \
  --o-filtered-table P1C-table-F20.qza

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
  --i-table ../P1C-filtered-table.qza \
  --m-metadata-file ../input/jay-met.tsv \
  --p-where "[CollectionYear]='2020' AND [CollectionSeason]='Spring'" \
  --o-filtered-table P1C-table-S20.qza

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





