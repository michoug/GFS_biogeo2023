### Comparisons with other samples from Cryospheric Ecosystems from Massimo's paper

library(phyloseq)
library(phyloseqCompanion)
library(data.table)
library(tidyverse)
library(Rmisc)
library(vegan)
library(wesanderson)
library(paletteer)
library(FSA)

## Read in cryobiome table 
asv_table_NOMISCRYO = fread("cryobiome_table.txt", header=T)

## Load metadata CRYO
metadata_cryo="PP1_metadata_date.tsv"
metadata_cryo<-import_qiime_sample_data(metadata_cryo)

## We start here with an ASV table that takes into account the 149 samples and the ASVs that were filtered -- We keep samples from Uganda
merge_nomis_cryo <- merge(asv_filtered, asv_table_NOMISCRYO, by.x="row.names", by.y="Feature_ID", all.y=T, all.x=T)

## Now we need to keep only data from Cryospheric ecosystems
## Load metadata CRYO + NOMIS
metadata_merge_cryonom="metadata_merge_nomis_cryo_deblur.txt"
metadata_cryo_nomis <-import_qiime_sample_data(metadata_merge_cryonom)

## Replace NAs by 0 in merge_nomis_cryo because those NA are taxa that were absent from the GFS
merge_nomis_cryo[is.na(merge_nomis_cryo)] <- 0

#saveRDS(merge_nomis_cryo,"20230828_merge_nomis_cryo.RDS")
#merge_nomis_cryo <- readRDS("merge_nomis_cryo.RDS")
asv_table_NOMISCRYO = as.data.frame(merge_nomis_cryo)

unAsIs <- function(X) {
   if("AsIs" %in% class(X)) {
     class(X) <- class(X)[-match("AsIs", class(X))]
   }
  X
}

asv_table_NOMISCRYO$Row.names<-unAsIs(asv_table_NOMISCRYO$Row.names)
rownames(asv_table_NOMISCRYO) <- asv_table_NOMISCRYO$Row.names
asv_table_NOMISCRYO$Row.names <- NULL
asv_table_NOMISCRYO_m <- as.matrix(asv_table_NOMISCRYO)
## transform as otu table --> yeahhh
asv_tableNOMISCRYO_final <- otu_table(asv_table_NOMISCRYO_m, taxa_are_rows=T)

##now merge asv_table and metadata
merge_asvnomiscryo_metadata <- merge_phyloseq(asv_tableNOMISCRYO_final, metadata_cryo_nomis)

##Remove everything else that is not cryospheric!!!
eco = c("Yes")
Cryosphere_sample <- subset_samples(merge_asvnomiscryo_metadata, Cryosphere %in% eco)

# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 65520 taxa and 428 samples ]:
#   sample_data() Sample Data:        [ 428 samples by 9 sample variables ]:
#   taxa are rows

library(R.filesets)
#saveRDS(Cryosphere_sample,"20232808_Cryosphere_sample.RDS")
Cryosphere_sample <- readRDS("20232808_Cryosphere_sample.RDS")
prune_nomiscryo <- prune_samples(sample_sums(Cryosphere_sample) >10000,Cryosphere_sample)
NOMIS_CRYOS<- rarefy_even_depth(prune_nomiscryo, sample.size=min(sample_sums(prune_nomiscryo)), rngseed=678, replace=F, trimOTUs=TRUE, verbose=TRUE)

# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 57499 taxa and 327 samples ]:
#   sample_data() Sample Data:      [ 327 samples by 9 sample variables ]:
#   taxa are rows

#saveRDS(NOMIS_CRYOS,file="20232808_NOMISCRYOS.RData")
NOMIS_CRYOS<-readRDS("20232808_NOMISCRYOS.RData")
#load("NOMISCROYS.RData")
meta <- sample.data.frame(NOMIS_CRYOS)

## Dans la table à la base j'ai remplacé ce qui est permafrost par "terrestrial"
NOMIS_CRYOS_sub <- subset_samples(NOMIS_CRYOS, Habitat != "Rock")
NOMIS_CRYOS_sub <- subset_samples(NOMIS_CRYOS_sub, Habitat != "Ice/Snow")
NOMIS_CRYOS_sub <- subset_samples(NOMIS_CRYOS_sub, Habitat !="Water_lake")
NOMIS_CRYOS_sub <- subset_samples(NOMIS_CRYOS_sub, Habitat !="Water")
NOMIS_CRYOS_sub <- subset_samples(NOMIS_CRYOS_sub, Habitat !="Sediment")

## Calculate ASV richness
diversity_cryosphere <- plot_richness(NOMIS_CRYOS_sub, measures=c("Observed","Shannon"), color="Habitat", shape="Ecosystem")
alphadt_cryosphere <-data.table(diversity_cryosphere$data)

OR_cryo <- subset(alphadt_cryosphere, variable == "Observed")
ASVrichness_cryo <- OR_cryo %>% 
  group_by(Habitat) %>% 
  summarise(average=mean(value), std=sd(value))

# cryosphere_richness <- ggplot(OR_cryo,aes(x=Habitat,y=value, color=Habitat)) + 
#   geom_boxplot() + 
#   geom_jitter(aes(group=Habitat), position=position_jitterdodge()) + 
#   scale_fill_viridis(discrete = TRUE) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           panel.background = element_blank(), axis.line = element_line(colour = "black"))

median_OR_cryosphere <- OR_cryo %>%
  group_by(Habitat) %>% 
  summarize(median_x=median(value)) 

# GFS_cryoconite <- 1048/90
# [1] 11.64444
# GFS_glacier <- 1048/303
# [1] 3.458746
# GFS_Microbial_mat <- 1048/211
# [1] 4.966825
# GFS_Permafrost <- 1048/114
# [1] 9.192982
# GFS_snow <- 1048/164
# [1] 6.390244
# Terrestrial_GFS <- 1859/1048
# [1] 1.773855


Shannon_cryo <- subset(alphadt_cryosphere, variable == "Shannon")
Shannondiv_cryo <- Shannon_cryo %>% 
  group_by(Habitat) %>% 
  summarize(median=median(value))

# Shannondiv_cryo
# # A tibble: 7 × 2
# Habitat       median
# <chr>          <dbl>
# 1 Cryoconite      2.84
# 2 GFS             5.40  5.4/2.84 (about 2 times)
# 3 Glacier         4.26  5.4/4.26
# 4 Microbial mat   4.93  5.4/4.93
# 5 Permafrost      4.29
# 6 Snow            2.22
# 7 Terrestrial     6.61  6.61/5.4 (1.2 times)


# cryosphere_shannon <- ggplot(Shannon_cryo,aes(x=Habitat,y=value, color=Habitat)) + 
#   geom_boxplot()+ 
#   geom_jitter(aes(group=Habitat), position=position_jitterdodge()) + 
#   scale_fill_viridis(discrete = TRUE) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))

## NMDS cryospheric ecosystems
metadata_nmds_cryo <- sample.data.frame(sample_data(NOMIS_CRYOS_sub))
ord <- ordinate(NOMIS_CRYOS_sub, method = "NMDS", distance = "bray", trymax = 999, k=2)
stressplot(ord)
ord$stress
# [1] 0.1717098
data.scores = as.data.frame(scores(ord)$sites)

#add columns to data frame 
data.scores$Sample = metadata_nmds_cryo$Sample
data.scores$Site = metadata_nmds_cryo$Habitat
head(data.scores)

##plot NMDS
##k=2
nmds_bc_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 6, aes(color = Site), alpha=0.6) + 
  #stat_ellipse(aes(x=NMDS1, y=NMDS2,color=Site),type = "norm")+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Biomes", y = "NMDS2", shape = "Type") 
# 
# ## Statistical analyses to compare diversity across different biomes
# library(performance)
# library(multcomp)
# habitats_richness <-glm(sqrt(value) ~ as.factor(Habitat), data = OR_cryo)
# summary(habitats_richness)
# aov_region <-anova(habitats_richness)
# check_model(habitats_richness)
# check_normality(habitats_richness)
# check_heteroscedasticity(habitats_richness)
# 
# ##ollowing normality // Homoscedasticity of variances
# habitats_richness_np <- aov(value ~ Habitat, data = OR_cryo)
# summary(habitats_richness_np)
# 
# OR_tukey <- TukeyHSD(habitats_richness_np, "Habitat",ordered=T, method="holm")
# 
# ##Shannon diversity
# habitats_shannon <-glm(log(value) ~ as.factor(Habitat), data = Shannon_cryo)
# summary(habitats_shannon)
# aov_region <-anova(habitats_shannon)
# check_model(habitats_shannon)
# check_normality(habitats_shannon)
# check_heteroscedasticity(habitats_shannon)
# 
# ## using Kruskal-Wallis
# 
# habitats_shannon_np <- kruskal.test(Shannon_cryo$value, Shannon_cryo$Habitat)
# habitats_shannon_np
# 
# shannon_dunn <- dunnTest(value ~ Habitat,
#                     data=Shannon_cryo,
#                     method="holm")
# 

## Testing the effects of the habitat on community structure

vegan_otu <- function(physeq){
  OTU <- otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

vegan_matrix_cryosphere <- vegan_otu(NOMIS_CRYOS_sub)
NOMISCRYO_bray <- vegdist(vegan_matrix_cryosphere, method="bray")
betadisp_cryo <- betadisper(NOMISCRYO_bray, metadata_nmds_cryo$Habitat, type="centroid")
betadisp_cryo_perm <- permutest(betadisp_cryo, permutations=999)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
# Groups      6 0.49414 0.082357 19.634    999  0.001 ***
#   Residuals 258 1.08224 0.004195  

ado_habitat <- adonis2(NOMISCRYO_bray ~ Habitat, data=metadata_nmds_cryo, method="bray", permutations=999)

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = NOMISCRYO_bray ~ Habitat, data = metadata_nmds_cryo, permutations = 999, method = "bray")
# Df SumOfSqs      R2      F Pr(>F)    
# Habitat    6   19.372 0.15974 8.1746  0.001 ***
#   Residual 258  101.899 0.84026                  
# Total    264  121.271 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## Pairwise adonis
library(pairwiseAdonis)
pairwise_habitat <- pairwise.adonis(NOMISCRYO_bray, metadata_nmds_cryo$Habitat, p.adjust.m="holm")

# pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1          GFS vs Microbial mat  1 0.7936735  1.930326 0.01270534   0.003      0.021   .
# 2             GFS vs Permafrost  1 3.2581419  7.974434 0.04610181   0.001      0.021   .
# 3            GFS vs Terrestrial  1 6.2956915 15.610868 0.07667011   0.001      0.021   .
# 4                   GFS vs Snow  1 1.7102691  4.165223 0.02616918   0.001      0.021   .
# 5                GFS vs Glacier  1 3.2794546  8.257997 0.04966977   0.001      0.021   .
# 6             GFS vs Cryoconite  1 5.5548641 13.620510 0.06962721   0.001      0.021   .
# 7   Microbial mat vs Permafrost  1 0.7872700  1.957969 0.09342361   0.002      0.021   .
# 8  Microbial mat vs Terrestrial  1 0.8650931  2.264736 0.05116344   0.005      0.021   .
# 9         Microbial mat vs Snow  1 0.6740520  1.567946 0.14836811   0.023      0.028   .
# 10     Microbial mat vs Glacier  1 1.1620984  4.694240 0.28118920   0.006      0.021   .
# 11  Microbial mat vs Cryoconite  1 0.8085533  2.014360 0.05298944   0.014      0.028   .
# 12    Permafrost vs Terrestrial  1 3.0656546  8.021632 0.12336867   0.001      0.021   .
# 13           Permafrost vs Snow  1 1.4021998  3.501994 0.12733599   0.001      0.021   .
# 14        Permafrost vs Glacier  1 2.8258117  8.758428 0.24493325   0.001      0.021   .
# 15     Permafrost vs Cryoconite  1 2.8700200  7.249404 0.12445455   0.001      0.021   .
# 16          Terrestrial vs Snow  1 1.7721111  4.623756 0.08956644   0.001      0.021   .
# 17       Terrestrial vs Glacier  1 3.8005846 11.103100 0.18171091   0.001      0.021   .
# 18    Terrestrial vs Cryoconite  1 4.5423261 11.755140 0.13707796   0.001      0.021   .
# 19              Snow vs Glacier  1 2.0471032  7.043695 0.29295392   0.001      0.021   .
# 20           Snow vs Cryoconite  1 1.6103754  4.020984 0.08931356   0.001      0.021   .
# 21        Glacier vs Cryoconite  1 3.4575313  9.801109 0.18217299   0.001      0.021   .


## Taxonomy taxo
## Read in cryobiome table 
asv_table_NOMISCRYO = as.data.frame(fread("cryobiome_table.txt", header=T))
rownames(asv_table_NOMISCRYO)<-asv_table_NOMISCRYO$Feature_ID

asv_table_NOMISCRYO$Feature_ID <- NULL
asv_table_NOMISCRYO <- as.matrix(asv_table_NOMISCRYO)

## Load metadata CRYO
metadata_cryo="PP1_metadata_date.tsv"
metadata_cryo<-import_qiime_sample_data(metadata_cryo)

tax_cryo <- as.data.frame(read.csv("cryobiome_taxonomy.csv"))
rownames(tax_cryo)<-tax_cryo$Feature.ID
tax_cryo$Feature.ID <- NULL
tax_cryo_m <- tax_table(as.matrix(tax_cryo))

asv_cryo <- otu_table(asv_table_NOMISCRYO, taxa_are_rows=T)
merge_cryo_taxo <- merge_phyloseq(asv_cryo, tax_cryo_m, metadata_cryo)

merge_cryo_taxo_sub <- subset_samples(merge_cryo_taxo, Habitat != "Rock")
merge_cryo_taxo_sub <- subset_samples(merge_cryo_taxo_sub, Habitat != "Ice/Snow")
merge_cryo_taxo_sub <- subset_samples(merge_cryo_taxo_sub, Habitat !="Water_lake")
merge_cryo_taxo_sub <- subset_samples(merge_cryo_taxo_sub, Habitat !="Water")
merge_cryo_taxo_sub <- subset_samples(merge_cryo_taxo_sub, Habitat !="Sediment")
merge_cryo_taxo_sub <- subset_samples(merge_cryo_taxo_sub, Habitat !="Mat")

cryo_taxglom <- tax_glom(merge_cryo_taxo_sub, taxrank=rank_names(merge_cryo_taxo_sub)[6], NArm=F)
transf_cryo <- transform_sample_counts(cryo_taxglom, function(x) x / sum(x))
sample_merge_habitat <- merge_samples(transf_cryo, "Habitat")
region_cryo <- transform_sample_counts(sample_merge_habitat, function(x) x / sum(x))

TopASV_f <- names(sort(taxa_sums(region_cryo), TRUE)[1:25])
top15_NOMIS_f <- prune_species(TopASV_f, region_cryo)
top15_NOMIS_f <- prune_taxa(taxa_sums(top15_NOMIS_f)>0, top15_NOMIS_f)
top_genus<-as.data.frame(tax_table(top15_NOMIS_f))

##Turn all ASVs into Genus counts
cryo_df <- psmelt(region_cryo) # create dataframe from phyloseq object
cryo_df$Genus<- as.character(cryo_df$Genus) #convert to character

## On va choisir de mettre en Other uniquement les ASVs qui ne font pas partie des Familles les plus abondantes
cryo_df$Genus[!(cryo_df$Genus %in% top_genus$Genus)] <- "Other"
cryo_df$Genus[(cryo_df$Genus == "")] <- "Other"
cryo_df$Genus[(cryo_df$Genus == "g__uncultured")] <- "Other"

##Generating n colors via rbrewer
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

barplot_biogeo <- ggplot(data=cryo_df, aes(x=Sample, y=Abundance, fill=Genus))
barplot_biogeo + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#1B9E77", "#D95F02", "#7570B3", "#E7298A","#66A61E",
                               "#E6AB02","#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6"
                              )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(ncol=3)) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




