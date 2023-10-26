##Heatmap
library(circlize)
library(ComplexHeatmap)
library(reshape2)
library(phyloseqCompanion)

character_data <- as.data.frame(sample_data(prune_Uganda))
asv_count_level <- otu_table(prune_Uganda, taxa_are_rows=T)
asv_count_levelf <- asv_count_level[rowSums(asv_count_level[])>0,]

asv_character_ab  = transform_sample_counts(prune_Uganda, function(x) x / sum(x))
asv_character <- otu_table(asv_character_ab, taxa_are_rows=T)
asv_characterf <- asv_character[rowSums(asv_character[])>0,]

##Relative abundance of the character taxa
##asv_character -> table abundance relative 
merge_character_abondance <- asv_characterf[row.names(asv_characterf) %in% row.names(character_table)]
merge_character_ab <- otu_table(merge_character_abondance, taxa_are_rows=T)

## Merge everything (mapping file, tree and OTUs) -- the metadatafile comes from the 1_loadingdata_biogeo_betadiv.R
merged_NOMIS_character_ab<- merge_phyloseq(merge_character_ab, tax_NOMIS, character_data)
metadata_character <- sample.data.frame(merged_NOMIS_character_ab)
#taxglom_biogeo_barplot <- tax_glom(physeq=NOMIS_FR, taxrank=rank_names(NOMIS_FR)[5], NArm=F)

merged_NOMIS_character_ab_otu <- (as.matrix(otu_table(merged_NOMIS_character_ab, taxa_are_rows=T)))

melt_asv_character <- melt(merged_NOMIS_character_ab_otu)
merge_asv_character <- merge(as.data.frame(melt_asv_character),as.matrix(metadata_nomis), by.x="Var2",by.y="sample")

##Here we would need to divide by the number of glaciers per mountain ranges. 
sum_mr_character <- merge_asv_character %>% group_by(site_c)%>% summarize(summrr=sum(value), n=n_distinct(Var2))%>% summarize(ar_mr=summrr/n, site_c)

median_indicator_ab <- sum_mr_character %>% 
  summarise(median=median(ar_mr), x = quantile(ar_mr, c(0.25, 0.5, 0.75)))
# # A tibble: 3 × 2
# median     x
# <dbl> <dbl>
# 1  0.210 0.170
# 2  0.210 0.210
# 3  0.210 0.420

##Lets get the most abundant ASVs - les premiers 500
TopASV_f <- names(sort(taxa_sums(merged_NOMIS_character_ab), TRUE)[1:500])
top500_NOMIS_f <- prune_species(TopASV_f, merged_NOMIS_character_ab)
top500_NOMIS_f <- prune_taxa(taxa_sums(top500_NOMIS_f)>0, top500_NOMIS_f)
top_ASV<-as.data.frame(tax_table(top500_NOMIS_f))

merge_character_ab_sub <- merge_character_ab[rownames(merge_character_ab) %in% rownames(top_ASV),]

hm_mat = (scale(log1p(as.matrix(merge_character_ab_sub))))

col_fun = colorRamp2(c(0, 1, 2), c("steelblue", "white", "darkred"))

ht <- ComplexHeatmap::Heatmap((hm_mat), use_raster = T, col = col_fun, column_split = metadata_character$site_c,
             show_column_names = F, show_row_names = F, name = 'Normalised abundance',
             row_dend_reorder = TRUE, show_row_dend = F, raster_quality = 5)

dev.off()
pdf('heatmap_indicator_2023.pdf', width = 9, height = 10)
draw(ht, auto_adjust = FALSE)
dev.off()


