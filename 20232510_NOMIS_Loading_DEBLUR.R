## Loading data & Processing  - Filtering, averaging and rarefaction

library(phyloseq)
library(vegan)
library(data.table)
library(speedyseq)
library(phyloseqCompanion)
library(dplyr)
library(tidyverse)
setwd("~/Leilas_Drive/Projet SBER/NOMIS/NOMIS_16S/NOMIS_August2023/DEBLUR")
asv<-fread(file="081623_NOMIS_16S_merged_deblur_table.csv",header=TRUE, sep=",")
newnames <- colnames(asv)
asv_df<-as.data.frame(asv)
colnames(asv_df) <- newnames
rownames(asv_df) <- asv_df$Feature_ID
asv_df$Feature_ID <- NULL

tax<-read.csv(file="081623_NOMIS_16S_merged_deblur_taxonomy.csv",sep=",",header=TRUE,row.names=1)

## Import QIIME pipeline map-file
metadata_NOMIS="/Users/ezzat/Leilas_Drive/Projet SBER/NOMIS/NOMIS_16S/NOMIS_August2023/DEBLUR/20230508_NOMIS_metadata_deblur.txt"
metadata_NOMIS<-import_qiime_sample_data(metadata_NOMIS)

OTU_NOMIS <- otu_table((asv_df), taxa_are_rows=TRUE)
tax_NOMIS <- tax_table(as.matrix(tax))

## Merge everything (mapping file, tree and OTUs)
merged_NOMIS_DEBLUR <- merge_phyloseq(OTU_NOMIS, tax_NOMIS, metadata_NOMIS)

# ## Check the nb of reads and sequences
# posi = c("Blank","Mock")
# blank_mock <- subset_samples(merged_NOMIS_DEBLUR, sample_type %in% posi)
# sediment <- subset_samples(merged_NOMIS_DEBLUR, sample_type == "sed")
# sediment_up_debut <- subset_samples(sediment, Position == "UP")
# merge_phylo_objects <- merge_phyloseq(blank_mock, sediment_up_debut)
# 
# sum(sample_sums(merge_phylo_objects))
# 
# asv_table_raw <- otu_table(merge_phylo_objects, taxa_are_rows=T)
# asv_table_raw <- asv_table_raw[rowSums(asv_table_raw[])>0,] ## remove rows containing 0 values 

# tax_table_raw <- tax_table(merged_NOMIS_DEBLUR)
# tax_table_raw <- subset(tax_table_raw, rownames(tax_table_raw) %in% rownames(asv_table_raw))
# 
# write.csv(asv_table_raw, "20231708_NOMIS_raw_deblur_table.csv")
# write.csv(tax_table_raw, "20231708_NOMIS_raw_deblur_taxonomy.csv")

## Check Blank
# merged_NOMIS_blank_check <- subset_samples(merged_NOMIS_DEBLUR, sample_type %in% "Blank")
# asv_blank_check <- otu_table(merged_NOMIS_blank_check)
# asv_blank_check <- as.data.frame(asv_blank_check)
# tax_blank_check <- as.data.frame(tax_table(merged_NOMIS_blank_check))
# 
# write.csv(asv_blank_check, "20231708_NOMIS_blank_check_deblur_table.csv")
# write.csv(tax_blank_check, "20231708_NOMIS_blank_check_deblur_taxonomy.csv")
# 
## Check MOCK
# merged_NOMIS_mock_check <- subset_samples(merged_NOMIS_DEBLUR, sample_type %in% "Mock")
# asv_mock_check <- otu_table(merged_NOMIS_mock_check)
# asv_mock_check <- as.data.frame(asv_mock_check)
# tax_mock_check <- as.data.frame(tax_table(merged_NOMIS_mock_check))
# 
# write.csv(asv_mock_check, "20231708_NOMIS_mock_deblur_table.csv")
# write.csv(tax_mock_check, "20231708_NOMIS_mock_deblur_taxonomy.csv")
# 
## metadata
nomis_meta <- sample.data.frame(merged_NOMIS_DEBLUR)

## Subsetting taxa ## Careful there are some spaces within the next code lines that are necessary to filter out specific taxa
merged_NOMIS_sub <- subset_taxa(merged_NOMIS_DEBLUR, (Kingdom!="d__Eukaryota") | is.na(Kingdom)) 
merged_NOMIS_sub <- subset_taxa(merged_NOMIS_sub, (Kingdom!="d__Archaea") | is.na(Kingdom))
merged_NOMIS_sub <- subset_taxa(merged_NOMIS_sub, (Order!=" o__Chloroplast") )
merged_NOMIS_sub <- subset_taxa(merged_NOMIS_sub, (Family!=" f__Mitochondria"))

## subset blank 
merged_NOMIS_blank <- subset_samples(merged_NOMIS_sub, sample_type %in% "Blank")
asv_blank <- otu_table(merged_NOMIS_blank)
asv_blank <- as.data.frame(asv_blank)

asv_blank <- asv_blank[rowSums(asv_blank[])>0,]

pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

## Dataset with no blank reads
merged_NOMIS_noblank <- pop_taxa(merged_NOMIS_sub, rownames(asv_blank))

# posi = c("Blank","Mock")
# blank_mock <- subset_samples(merged_NOMIS_noblank, sample_type %in% posi)
# sediment <- subset_samples(merged_NOMIS_noblank, sample_type == "sed")
# sediment_up_debut <- subset_samples(sediment, Position == "UP")
# merge_phylo_objects <- merge_phyloseq(blank_mock, sediment_up_debut)
# asv_sedup_unfiltered <- otu_table(merge_phylo_objects)
# asv_sedup_unfiltered <- as.data.frame(asv_sedup_unfiltered)
# asv_sedup_unfiltered <- asv_sedup_unfiltered[rowSums(asv_sedup_unfiltered[])>0,]

## Select Substrate and Position
site = c("sed")
sediment <- subset_samples(merged_NOMIS_noblank, sample_type %in% site)
posi = c("UP")
sediment_up <- subset_samples(sediment, Position %in% posi)

## ASV table
asv_sedup_unfiltered <- otu_table(sediment_up)
asv_sedup_unfiltered <- as.data.frame(asv_sedup_unfiltered)
asv_sedup_unfiltered <- asv_sedup_unfiltered[rowSums(asv_sedup_unfiltered[])>0,]

## filtering step - removes ASVs that were not found in 50% of biological replicates

set.seed(132)
filtering_br <- function(phyloseq_obj){
  meta_sediment <- as.data.frame(sample_data(phyloseq_obj))
  asv_sedimentup <- as.data.frame(otu_table(phyloseq_obj, taxa_are_rows = T))
  
  new_asv <- c()
  for (Sample in meta_sediment$gl_code) {
    sub_asvsed <- asv_sedimentup[,startsWith(colnames(asv_sedimentup), paste0(Sample,'_',collapse = ''))]
    if (is.null(sub_asvsed)){}
    else{
      if (!(is.vector(sub_asvsed))){
        sub_asv_preval <- rowMeans(sub_asvsed > 0)
        new_asv = c(new_asv, rownames(sub_asvsed)[sub_asv_preval >= 0.5])}
      else {
        new_asv = c(new_asv, names(sub_asvsed)[sub_asvsed > 0])}
    }}
  new_asv = unique(new_asv)
  return(new_asv)
}

## filtering
filtered_asvs = filtering_br(sediment_up)

## subset ASVs that were filtered from ASV table.
filtered_table = asv_sedup_unfiltered[rownames(asv_sedup_unfiltered) %in% filtered_asvs, ]
filtered_asv_table<- otu_table(filtered_table, taxa_are_rows=T)

##How many sample per patch?
# Get column names
col_names <- colnames(filtered_asv_table)

# Create a function to count prefixes
count_prefixes <- function(prefix, col_names) {
  count <- sum(grepl(paste0("^", prefix, "_"), col_names))
  return(count)
}

# List of unique prefixes
unique_prefixes <- unique(sub("_.*", "", col_names))

# Initialize a data frame to store results
result_df <- data.frame(Prefix = unique_prefixes, Count = numeric(length(unique_prefixes)))

# Iterate through prefixes and count
for (i in 1:length(unique_prefixes)) {
  prefix <- unique_prefixes[i]
  count <- count_prefixes(prefix, col_names)
  result_df[i, "Count"] <- count
}

# Display the result
print(result_df)

##subset to show how many times we found 3/2/1 samples per patch

results_df_patch_3 <- result_df %>%
  summarise(sumcout=sum(Count==3))
#134
results_df_patch_2 <- result_df %>%
  summarise(sumcout=sum(Count==2))
#14
results_df_patch_1 <- result_df %>%
  summarise(sumcout=sum(Count==1))
#3

## melt ASV_table
melt_asv <- melt(filtered_asv_table)
code_glaciers <- sample.data.frame(nomis_meta[,c("sample","gl_code")])
rownames(code_glaciers)<-NULL

## merge ASV_table and gl_code
melt_merge_as <- merge(melt_asv, code_glaciers, by.x="Var2", by.y="sample")
#saveRDS(melt_merge_as,"20231908_melt_merge_as.RData")
#melt_merge_as <- readRDS("20231908_melt_merge_as.RData")
## take the average and transform as integer
melt_merge_asv <- melt_merge_as %>%
  group_by(gl_code, Var1) %>% 
  summarise(average=mean(value))

## re-build ASV table
dcast_asv <- dcast(melt_merge_asv, formula= Var1 ~ gl_code)
rownames(dcast_asv) <- dcast_asv$Var1
dcast_asv$Var1 <- NULL

min_non_zero = min(dcast_asv[dcast_asv > 0])
dcast_asv = as.data.frame(dcast_asv)*1/min_non_zero
dcast_asv_n = as.data.frame(sapply(dcast_asv, round))
dcast_asv_n = as.data.frame(sapply(dcast_asv_n, as.integer))
rownames(dcast_asv_n) <- rownames(dcast_asv)
asv_table_dcast <-(otu_table(dcast_asv_n, taxa_are_rows=T))

## re-build phyloseq object
metadata_glaciers="/Users/ezzat/Leilas_Drive/Projet SBER/NOMIS/NOMIS_16S/NOMIS_August2023/DEBLUR/202308_NOMIS_metadata_GFS.tsv"
metadata_glaciers<-import_qiime_sample_data(metadata_glaciers)
tax_sed <- as.matrix(tax_table(sediment_up))
merge_GFS_full <- merge_phyloseq(asv_table_dcast, tax_sed, metadata_glaciers)

## prune GFS sample which shows less than 25k reads (2)
prune_GFS<- prune_samples(sample_sums(merge_GFS_full)>60000,merge_GFS_full)
## save filtered tables
asv_filtered <- otu_table(prune_GFS, taxa_are_rows=T)
# tax_table_filtered <- tax_table(prune_GFS)
# write.csv(asv_filtered, "20230825_NOMIS_filtered_deblur_table.csv")
# write.csv(tax_table_filtered, "20230825_NOMIS_filtered_deblur_taxonomy.csv")
# 
## rarefaction
NOMIS_rarefied <- rarefy_even_depth(prune_GFS, sample.size=min(sample_sums(prune_GFS)), rngseed=678, replace=F, trimOTUs=TRUE, verbose=TRUE)

# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 39434 taxa and 149 samples ]:
# sample_data() Sample Data:        [ 149 samples by 4 sample variables ]:
# tax_table()   Taxonomy Table:     [ 39434 taxa by 7 taxonomic ranks ]:
# taxa are rows

## save rarefied tables
# asv_rarefied <- otu_table(NOMIS_rarefied, taxa_are_rows=T)
# tax_table_rarefied <- tax_table(NOMIS_rarefied)
# write.csv(asv_rarefied, "20230825_NOMIS_rarefied_deblur_table.csv")
# write.csv(tax_table_rarefied, "20230825_NOMIS_rarefied_deblur_taxonomy.csv")

#NOMIS_R <- saveRDS(NOMIS_rarefied,"20230825_NOMIS_rarefied.RData")
NOMIS_R <- readRDS("20230825_NOMIS_rarefied.RData")
## remove Uganda sample
uganda=c("Uganda")
prune_Uganda <- subset_samples(NOMIS_R, !site_c %in% uganda)
prune_Uganda_df <- as.matrix(otu_table(prune_Uganda, taxa_are_rows=T))

metadata_nomis <- sample.data.frame(prune_Uganda)
metadata_glaciers <- sample.data.frame(prune_Uganda)
