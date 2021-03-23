# this script is for the environmental analysis
# to load R on interactive graham
# module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")

library(qiime2R)
library(vegan)
library(tidyverse)
# set theme
theme_set(theme_bw())

# read in the required data
nReads <- 344
otu <- read_qza("rarefied-table.qza")$data
map <- read_tsv("input/jay-met.tsv")[-1,]

otu_PA <- 1*((otu>0)==1)                                               # presence-absence data
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean)     # relative abundance 
occ_abun <- rownames_to_column(as.data.frame(cbind(otu_occ, otu_rel)),var="otu") # combining occupancy and abundance

# Ranking OTUs based on their occupancy
# For caluclating raking index we included following conditions:
#   - time-specific occupancy (sumF) = frequency of detection within time point (genotype or site)
#   - replication consistency (sumG) = has occupancy of 1 in at least one time point (genotype or site) (1 if occupancy 1, else 0)

PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>% 
  gather(sampleid, abun, -otu) %>%
  left_join(map, by = 'sampleid') %>%
  group_by(otu, Territory) %>%
  summarise(plot_freq=sum(abun>0)/length(abun),        # frequency of detection between time points
            coreSite=ifelse(plot_freq == 1, 1, 0), # 1 only if occupancy 1 with specific site, 0 if not
            detect=ifelse(plot_freq > 0, 1, 0)) %>%    # 1 if detected and 0 if not detected with specific site
  group_by(otu) %>%
  summarise(sumF=sum(plot_freq),
            sumG=sum(coreSite),
            nS=length(Territory)*2,
            Index=(sumF+sumG)/nS) # calculating weighting Index based on number of time points detected and 

otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by='otu') %>%
  transmute(otu=otu,
            rank=Index) %>%
  arrange(desc(rank))

# Calculating the contribution of ranked OTUs to the BC similarity
BCaddition <- NULL

otu_start <- otu_ranked$otu[1]
start_matrix <- as.matrix(otu[otu_start,])
start_matrix <- t(start_matrix)
x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(2*nReads))
x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
df_s <- data.frame(x_names,x)
names(df_s)[2] <- 1 
BCaddition <- rbind(BCaddition,df_s)
# calculating BC dissimilarity based on additon of ranked OTUs from 2nd to 644th (entire length).
for(i in 2:nrow(occ_abun)){
  otu_add <- otu_ranked$otu[i]
  add_matrix <- as.matrix(otu[otu_add,])
  add_matrix <- t(add_matrix)
  start_matrix <- rbind(start_matrix, add_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_a <- data.frame(x_names,x)
  names(df_a)[2] <- i 
  BCaddition <- left_join(BCaddition, df_a, by=c('x_names'))  
}

x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
df_full <- data.frame(x_names,x)
names(df_full)[2] <- length(rownames(otu))
BCfull <- left_join(BCaddition,df_full, by='x_names')

rownames(BCfull) <- BCfull$x_names
temp_BC <- BCfull
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC=mean(BC)) %>%
  arrange(-desc(MeanBC)) %>%
  mutate(proportionBC=MeanBC/max(MeanBC))

Increase <- BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
# moved the 0 in c to after Increase
increaseDF <- data.frame(IncreaseBC=c((Increase),0), rank=factor(c(1:(length(Increase)+1))))
BC_ranked <- left_join(BC_ranked, increaseDF, by = "rank")
BC_ranked <- BC_ranked[-nrow(BC_ranked),]


#Creating occupancy abundance plot
occ_abun$fill <- 'Rare'
# correcting the NA issue with the cbinds
occ_abun$fill[occ_abun$otu %in% cbind(BC_ranked, otu_ranked)$otu[cbind(BC_ranked, otu_ranked)$IncreaseBC>=1.02]] <- 'Core'

# add taxonomy information to occ_abun
tax <- read_qza("taxonomy/SILVA-taxonomy.qza")$data
# clean up/separatee taxonomy labels
tax$Taxon <- tax$Taxon %>%
  str_replace_all("d__", "") %>%
  str_replace_all("p__", "") %>%
  str_replace_all("c__", "") %>%
  str_replace_all("o__", "") %>%
  str_replace_all("f__", "") %>%
  str_replace_all("g__", "") %>%
  str_replace_all("s__", "")
occ_abunT <- tax %>% 
  mutate("Taxon" = gsub(".*__", "", .$Taxon),
         "Kingdom" = word(.$Taxon, 1, sep = ";"),
         "Phylum" = word(.$Taxon, 2, sep = ";"),
         "Class" = word(.$Taxon, 3, sep = ";"),
         "Order" = word(.$Taxon, 4, sep = ";"),
         "Family" = word(.$Taxon, 5, sep = ";"),
         "Genus" = word(.$Taxon, 6, sep = ";"),
         "Species" = word(.$Taxon, 7, sep = ";")) %>%
  # join to core labels
  right_join(occ_abun, by = c("Feature.ID" = "otu" ))

# filter for core names only
core <- occ_abunT[occ_abunT$fill == "Core",] %>% 
  rename(c("featureid" = "Feature.ID")) %>%
  select(featureid)

# save information to file
write.table(occ_abunT, file = "CanadaJayMicrobiome/data/coreJay.csv", sep = ",", quote = F, row.names = F)
write.table(core, file = "CanadaJayMicrobiome/data/coreFeatures.tsv", sep = "\t", quote = F, row.names = F)
