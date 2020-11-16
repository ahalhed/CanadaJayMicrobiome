# to open a qiime feature table from command line
# first time, install the biom dependancies
pip install numpy
pip install biom-format

# unzip the qza file (it's just a zip archive)
unzip filtered-table-no-singletons-mitochondria-chloroplast.qza

# change directories into the newly unzipped folder
# access the data folder
cd ./1d92e26b-e4f5-480a-bc3a-13f8843551bf/data/

# convert the biom file into a summary text file
biom summarize-table -i feature-table.biom --qualitative -o feature-table-summary.txt

# need to finish this
# http://www.metagenomics.wiki/tools/16s/qiime/otu-biom-table