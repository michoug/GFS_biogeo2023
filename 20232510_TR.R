## Taxonomic Tree to represent GFS microbiome
library(phyloseq)
library(phyloseqCompanion)
library(vegan)
library(metacoder)
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(agricolae)
library(ape)

NOMIS_R <- readRDS("20230825_NOMIS_rarefied.RData")

# Build the object that could then be used by metacoder including sample from Uganda!
asv_taxo_metacoder <- read.csv("asv_taxo_metacoder082023.csv",header=T)
asv_metacoder_tibble <- as.tibble(asv_taxo_metacoder)

# load metadata
meta_glaciers_tibble <- as.tibble(sample.data.frame(NOMIS_R))

# parse
obj_tree_nomis <- parse_tax_data(asv_metacoder_tibble,
                      class_cols = "lineage",
                      class_sep = "; ", ##careful here because there is a space also!
                      class_regex = "^(.+)__(.+)$", 
                      class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))

# we did not choose to filter more ASVs here
no_reads <- rowSums(obj_tree_nomis$data$tax_data[, meta_glaciers_tibble$sample]) == 0
sum(no_reads)

# remove taxa with no names
obj_tree_nomis <- filter_taxa(obj_tree_nomis, taxon_names != "")
obj_tree_nomis <- filter_taxa(obj_tree_nomis, taxon_names != "uncultured")
obj_tree_nomis <- filter_taxa(obj_tree_nomis, taxon_names != "Bacteria;;;;;")


obj_tree_noms <- obj_tree_nomis %>% 
filter_taxa(taxon_ranks == "g", supertaxa=T)# subset to the order rank
  
# Getting per-taxon information
# We do have info regarding the abundance of each ASV. But to get to know more about the taxa,
# we can sum the abundance per-taxon and add the results to the tax-map
obj_tree_nomis$data$tax_abund <- calc_taxon_abund(obj_tree_nomis, "tax_data",
                                       cols = meta_glaciers_tibble$sample)


# Calculate the nb of samples that have reads for each taxon
obj_tree_nomis$data$tax_occ <- calc_n_samples(obj_tree_nomis, "tax_abund", cols = meta_glaciers_tibble$sample)
#
# Plot taxonomic tree
set.seed(3) # This makes the plot appear the same each time it is run 
heat_tree_nomis <- heat_tree(obj_tree_nomis, 
                            node_color=n_obs,
                            node_size=n_obs,
                            node_label = taxon_names,
                            edge_color_range = c("#CC8394FF", "#AC563BFF", "#CDA97CFF","#7C8EC5FF","#2B3C51FF"),
                            node_color_range = c("#CC8394FF", "#AC563BFF", "#CDA97CFF","#7C8EC5FF","#2B3C51FF"),
                            edge_color=n_samples,
                            initial_layout = "re", layout = "da")
heat_tree_nomis
