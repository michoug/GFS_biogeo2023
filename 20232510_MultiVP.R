## Multivariate analysis and Variance Partitioning
## We start with rarefied data -- 148 glaciers -- Updated August 23
library(phyloseq)
library(phyloseqCompanion)
library(tidyverse)


##Community Matrix## 
library(R.filesets)
NOMIS_FR <- readRDS("20230825_NOMIS_rarefied.RData")
metadata_nomis<- sample.data.frame(NOMIS_FR)

uganda=c("Uganda")
prune_Uganda <- subset_samples(NOMIS_FR, !site_c %in% uganda)
comm_table <- as.data.frame(otu_table(prune_Uganda, taxa_are_rows=T))
metadata_nomis_wu <- sample.data.frame(prune_Uganda)

## Environmental data
## We need to include chla and climatic data
metadata_nomis_db <- read.csv("nomis_metadata_database_2023.csv")
minerals_nomis <- read.csv("minerals_nomis_2023.csv")
chla_nomis <- read.csv("chla_nomis_2023.csv")
chla_mean <- chla_nomis %>% 
  group_by(gl_code)%>% summarize(mean_chla=mean(chla))

chla_mean_df <- as.data.frame(chla_mean)
climate_nomis <- read.csv("climate_nomis_2023.csv")

metadata_merge <- merge(metadata_nomis_db, metadata_nomis_wu, by="sample")
metadata_merge <- merge(metadata_merge, minerals_nomis, by="sample")
metadata_merge <- merge(metadata_merge, chla_mean, by.x="sample", by.y="gl_code")
metadata_merge <- merge(metadata_merge, climate_nomis, by.x="sample", by.y="Sample")

## Keep the data you want (to include ions and minerals) -- We do not include Clays as this is collinear with another mineral!
multi_data_subset <- metadata_merge[c("sample","site_c","water_temp", "pH", "cond", "turb","doc","srp", "nh4","no3","no2", "mean_chla","gl_sa.x","gl_cov.x","lat.x","lon.x","ele","Feldspar","Calcite","Quartz","pr","tas","scd")]

multi_data_subset<- multi_data_subset  %>% 
  mutate(no2 = coalesce(no2,0),no3=coalesce(no3,0),nh4=coalesce(nh4,0))

#Sum inorganic nitrogen into DIN
rowsum_nut <- rowSums(multi_data_subset[, c("no2", "no3", "nh4")])
multi_data_subset$DIN <- rowsum_nut 

## Convert col to numeric
convert_columns_to_numeric <- function(data, columns) {
  for (col in columns) {
    data[[col]] <- as.numeric(data[[col]])
  }
  return(data)
}
multi_data_num <- convert_columns_to_numeric(multi_data_subset, c("water_temp", "pH", "cond", "turb","doc","srp", "mean_chla","gl_sa.x","gl_cov.x","lat.x","lon.x","ele","Feldspar","Calcite","Quartz","pr","scd","DIN"))

## Remove duplicate lines
multi_data_num_u <- distinct(multi_data_num)
multi_data_num_u <- multi_data_num_u[c("sample","site_c","water_temp", "pH", "cond", "turb","doc","srp", "mean_chla","gl_sa.x","gl_cov.x","lat.x","lon.x","ele","Feldspar","Calcite","Quartz","pr","scd","DIN")]

## Remove NA values -> 148 to 132 GFS
hist(multi_data_na$water_temp)

## Now we have to deal with 0 values!
## There are 0 values in water_temp, mean_chla, scd
## We add a minimum value for columns which show 0

multi_data_na <-(na.omit(multi_data_num_u))

add_const <- function(x) {
  min_nonzero <- min(x[which(x > 0)])  # Find the minimum nonzero value
  return((x + (min_nonzero/2)))
}

# Specify the column names you want to modify
columns_to_modify <- c("water_temp", "scd", "mean_chla")

# Apply add_const only to the specified columns
multi_data_trans <- multi_data_na %>%
  mutate_at(vars(all_of(columns_to_modify)), add_const)

print(multi_data_trans)

## gl_code are the rownames
comm_table_t <- as.data.frame(t(comm_table))

## rename the first column (là en fait on ajoute une column)
comm_table_t$sample<- rownames(comm_table_t)

## Geographic Distance Matrix
metadata_distance <- as.data.frame(multi_data_trans %>% select(lon.x, lat.x))
metadata_dist_df <- as.data.frame(multi_data_trans %>% select(lon.x, lat.x, ele))
metadata_dist_y <- metadata_dist_df
metadata_dist_df <- as.data.frame(sapply(metadata_dist_df, as.numeric))
library(geosphere)
dist_geo<-distm(metadata_distance, fun=distGeo)

## Environmental data without geographic coordinates and altitude
multi_data_trans_ming <-subset(multi_data_trans, select=-c(lat.x, lon.x, ele))

## Community MATRIX with filtered samples -- Keep only the samples that do not show any NA
asv_community_trim  <- comm_table_t%>% 
  filter(sample %in% multi_data_trans$sample)
## We have to remove the last column, which contains a "sample"
row.names(asv_community_trim) == multi_data_trans$sample
colnames(asv_community_trim)[39435] ## sample
Y.com <- asv_community_trim[,1:(dim(asv_community_trim)[2]-1)]
## Convert into numeric values
Y.com <- sapply(Y.com, as.numeric)
## Saving Files 
save(Y.com,file="comY_20232309.Rdata")
save(multi_data_trans_ming, file="Env_20232309.Rdata")
save(dist_geo,file= "Geo_20232309.Rdata")
save(metadata_dist_df,file="geocoordinates_20232309.Rdata")

####################################################################################################

## Multivariate analyses core microbiomes :) 
core_table <-read.csv("core_list_2023.csv")

merge_core_count <- merge(comm_table, core_table, by.x="row.names", by.y="ASV" )

# Specify the column after which you want to remove all columns
column_to_keep <- "GL99"
# Find the index of the specified column
column_index <- which(names(merge_core_count) == column_to_keep)
# Select columns up to the specified column
merge_core_count_f <- merge_core_count[, 1:column_index]
# View the resulting data frame
print(merge_core_count_f)

## Formate row names
rownames(merge_core_count_f)<- merge_core_count_f$Row.names
merge_core_count_f$Row.names <- NULL

merge_core_count_f <- otu_table(merge_core_count_f, taxa_are_rows=T)
tax_core <-as.matrix(tax_table(NOMIS_FR))
metadata_glaciers="/Users/ezzat/Leilas_Drive/Projet SBER/NOMIS/NOMIS_16S/NOMIS_August2023/DEBLUR/202308_NOMIS_metadata_GFS.tsv"
metadata_glaciers<-import_qiime_sample_data(metadata_glaciers)
## Merge everything (mapping file, tree and OTUs) -- the metadatafile comes from the 1_loadingdata_biogeo_betadiv.R
core_for_multivariate<- merge_phyloseq(merge_core_count_f, tax_core, metadata_glaciers)
metadata_core <- sample.data.frame(core_for_multivariate)
core_fm <- otu_table(core_for_multivariate, taxa_are_rows=T)

## Environmental data
## We need to include chla and climatic data
metadata_nomis_db <- read.csv("nomis_metadata_database_2023.csv")
minerals_nomis <- read.csv("minerals_nomis_2023.csv")
chla_nomis <- read.csv("chla_nomis_2023.csv")
chla_mean <- chla_nomis %>% 
  group_by(gl_code)%>% summarize(mean_chla=mean(chla))

chla_mean_df <- as.data.frame(chla_mean)
climate_nomis <- read.csv("climate_nomis_2023.csv")

metadata_core_merge <- merge(metadata_nomis_db, metadata_core, by="sample")
metadata_core_merge <- merge(metadata_core_merge, minerals_nomis, by="sample")
metadata_core_merge <- merge(metadata_core_merge, chla_mean, by.x="sample", by.y="gl_code")
metadata_core_merge <- merge(metadata_core_merge, climate_nomis, by.x="sample", by.y="Sample")

## Keep the data you want (to include ions and minerals) -- We do not include Clays as this is collinear with another mineral!
multi_data_core_subset <- metadata_core_merge[c("sample","site_c","water_temp", "pH", "cond", "turb","doc","srp", "nh4","no3","no2", "mean_chla","gl_sa.x","gl_cov.x","lat.x","lon.x","ele","Feldspar","Calcite","Quartz","pr","tas","scd")]

multi_data_core_subset<- multi_data_core_subset  %>% 
  mutate(no2 = coalesce(no2,0),no3=coalesce(no3,0),nh4=coalesce(nh4,0))

#Sum inorganic nitrogen into DIN
rowsum_nut <- rowSums(multi_data_core_subset[, c("no2", "no3", "nh4")])
multi_data_core_subset$DIN <- rowsum_nut 

## Convert col to numeric
convert_columns_to_numeric <- function(data, columns) {
  for (col in columns) {
    data[[col]] <- as.numeric(data[[col]])
  }
  return(data)
}
multi_data_core_num <- convert_columns_to_numeric(multi_data_core_subset, c("water_temp", "pH", "cond", "turb","doc","srp", "mean_chla","gl_sa.x","gl_cov.x","lat.x","lon.x","ele","Feldspar","Calcite","Quartz","pr","scd","DIN"))

## Remove duplicate lines
multi_data_num_core_u <- distinct(multi_data_core_num)
multi_data_num_core_u <- multi_data_num_core_u[c("sample","site_c","water_temp", "pH", "cond", "turb","doc","srp", "mean_chla","gl_sa.x","gl_cov.x","lat.x","lon.x","ele","Feldspar","Calcite","Quartz","pr","scd","DIN")]

## Now we have to deal with 0 values!
## There are 0 values in water_temp, mean_chla, scd
## We add a minimum value for columns which show 0

multi_data_core_na <-(na.omit(multi_data_num_core_u))

add_const <- function(x) {
  min_nonzero <- min(x[which(x > 0)])  # Find the minimum nonzero value
  return((x + (min_nonzero/2)))
}

# Specify the column names you want to modify
columns_to_modify <- c("water_temp", "scd", "mean_chla")

# Apply add_const only to the specified columns
multi_data_core_trans <- multi_data_core_na %>%
  mutate_at(vars(all_of(columns_to_modify)), add_const)

print(multi_data_core_trans)

## gl_code are the rownames
comm_core_table_t <- as.data.frame(t(core_fm))

## rename the first column (là en fait on ajoute une column)
comm_core_table_t$sample<- rownames(comm_core_table_t)

## Geographic Distance Matrix
metadata_distance_core <- as.data.frame(multi_data_core_trans %>% select(lon.x, lat.x))
metadata_dist_core_df <- as.data.frame(multi_data_core_trans %>% select(lon.x, lat.x, ele))
metadata_dist_core_y <- metadata_dist_core_df
metadata_dist_core_df <- as.data.frame(sapply(metadata_dist_core_df, as.numeric))
library(geosphere)
dist_geo_core<-distm(metadata_distance_core, fun=distGeo)

## Environmental data without geographic coordinates and altitude
multi_data_trans_core_ming <-subset(multi_data_core_trans, select=-c(lat.x, lon.x, ele))

## Community MATRIX with filtered samples -- Keep only the samples that do not show any NA
asv_community_trim_core  <- comm_core_table_t%>% 
  filter(sample %in% multi_data_core_trans$sample)
## We have to remove the last column, which contains a "sample"
row.names(asv_community_trim_core) == multi_data_core_trans$sample
colnames(asv_community_trim_core)[165] ## sample
Y.com_core <- asv_community_trim_core[,1:(dim(asv_community_trim_core)[2]-1)]
## Convert into numeric values
Y.com_core <- sapply(Y.com_core, as.numeric)
## Saving Files 
save(Y.com_core,file="comY_core_20232309.Rdata")
save(multi_data_trans_core_ming, file="Env_core_20232309.Rdata")
save(dist_geo_core,file= "Geo_core_20232309.Rdata")
save(metadata_dist_core_df,file="geocoordinates_core_20232309.Rdata")

# 
# ################## Multivariate analysis endemic ASVs
# endemic_list <- read.csv("endemic_list_2023.csv")
# 
# merge_endemic_count <- merge(comm_table, endemic_list, by.x="row.names", by.y="Var2" )
# 
# # Specify the column after which you want to remove all columns
# column_to_keep <- "GL99"
# # Find the index of the specified column
# column_index <- which(names(merge_endemic_count) == column_to_keep)
# # Select columns up to the specified column
# merge_endemic_count_f <- merge_endemic_count[, 1:column_index]
# # View the resulting data frame
# print(merge_endemic_count_f)
# 
# ## Formate row names
# rownames(merge_endemic_count_f)<- merge_endemic_count_f$Row.names
# merge_endemic_count_f$Row.names <- NULL
# 
# merge_endemic_count_f <- otu_table(merge_endemic_count_f, taxa_are_rows=T)
# tax_endemic <-as.matrix(tax_table(NOMIS_FR))
# metadata_glaciers="/Users/ezzat/Leilas_Drive/Projet SBER/NOMIS/NOMIS_16S/NOMIS_August2023/DEBLUR/202308_NOMIS_metadata_GFS.tsv"
# metadata_glaciers<-import_qiime_sample_data(metadata_glaciers)
# ## Merge everything (mapping file, tree and OTUs) -- the metadatafile comes from the 1_loadingdata_biogeo_betadiv.R
# endemic_for_multivariate<- merge_phyloseq(merge_endemic_count_f, tax_endemic, metadata_glaciers)
# metadata_endemic <- sample.data.frame(endemic_for_multivariate)
# endemic_fm <- otu_table(endemic_for_multivariate, taxa_are_rows=T)
# 
# ## Environmental data
# ## We need to include chla and climatic data
# metadata_nomis_db <- read.csv("nomis_metadata_database_2023.csv")
# minerals_nomis <- read.csv("minerals_nomis_2023.csv")
# chla_nomis <- read.csv("chla_nomis_2023.csv")
# chla_mean <- chla_nomis %>% 
#   group_by(gl_code)%>% summarize(mean_chla=mean(chla))
# 
# chla_mean_df <- as.data.frame(chla_mean)
# climate_nomis <- read.csv("climate_nomis_2023.csv")
# 
# metadata_endemic_merge <- merge(metadata_nomis_db, metadata_endemic, by="sample")
# metadata_endemic_merge <- merge(metadata_endemic_merge, minerals_nomis, by="sample")
# metadata_endemic_merge <- merge(metadata_endemic_merge, chla_mean, by.x="sample", by.y="gl_code")
# metadata_endemic_merge <- merge(metadata_endemic_merge, climate_nomis, by.x="sample", by.y="Sample")
# 
# ## Keep the data you want (to include ions and minerals) -- We do not include Clays as this is collinear with another mineral!
# multi_data_endemic_subset <- metadata_endemic_merge[c("sample","site_c","water_temp", "pH", "cond", "turb","doc","srp", "nh4","no3","no2", "mean_chla","gl_sa.x","gl_cov.x","lat.x","lon.x","ele","Feldspar","Calcite","Quartz","pr","tas","scd")]
# 
# multi_data_endemic_subset<- multi_data_endemic_subset  %>% 
#   mutate(no2 = coalesce(no2,0),no3=coalesce(no3,0),nh4=coalesce(nh4,0))
# 
# #Sum inorganic nitrogen into DIN
# rowsum_nut <- rowSums(multi_data_endemic_subset[, c("no2", "no3", "nh4")])
# multi_data_endemic_subset$DIN <- rowsum_nut 
# 
# ## Convert col to numeric
# convert_columns_to_numeric <- function(data, columns) {
#   for (col in columns) {
#     data[[col]] <- as.numeric(data[[col]])
#   }
#   return(data)
# }
# multi_data_endemic_num <- convert_columns_to_numeric(multi_data_endemic_subset, c("water_temp", "pH", "cond", "turb","doc","srp", "mean_chla","gl_sa.x","gl_cov.x","lat.x","lon.x","ele","Feldspar","Calcite","Quartz","pr","scd","DIN"))
# 
# ## Remove duplicate lines
# multi_data_num_endemic_u <- distinct(multi_data_endemic_num)
# multi_data_num_endemic_u <- multi_data_num_endemic_u[c("sample","site_c","water_temp", "pH", "cond", "turb","doc","srp", "mean_chla","gl_sa.x","gl_cov.x","lat.x","lon.x","ele","Feldspar","Calcite","Quartz","pr","scd","DIN")]
# 
# ## Now we have to deal with 0 values!
# ## There are 0 values in water_temp, mean_chla, scd
# ## We add a minimum value for columns which show 0
# 
# multi_data_endemic_na <-(na.omit(multi_data_num_endemic_u))
# 
# add_const <- function(x) {
#   min_nonzero <- min(x[which(x > 0)])  # Find the minimum nonzero value
#   return((x + (min_nonzero/2)))
# }
# 
# # Specify the column names you want to modify
# columns_to_modify <- c("water_temp", "scd", "mean_chla")
# 
# # Apply add_const only to the specified columns
# multi_data_endemic_trans <- multi_data_endemic_na %>%
#   mutate_at(vars(all_of(columns_to_modify)), add_const)
# 
# print(multi_data_endemic_trans)
# 
# ## gl_code are the rownames
# comm_endemic_table_t <- as.data.frame(t(endemic_fm))
# 
# ## rename the first column (là en fait on ajoute une column)
# comm_endemic_table_t$sample<- rownames(comm_endemic_table_t)
# 
# ## Geographic Distance Matrix
# metadata_distance_endemic <- as.data.frame(multi_data_endemic_trans %>% select(lon.x, lat.x))
# metadata_dist_endemic_df <- as.data.frame(multi_data_endemic_trans %>% select(lon.x, lat.x, ele))
# metadata_dist_endemic_y <- metadata_dist_endemic_df
# metadata_dist_endemic_df <- as.data.frame(sapply(metadata_dist_endemic_df, as.numeric))
# library(geosphere)
# dist_geo_endemic<-distm(metadata_distance_endemic, fun=distGeo)
# 
# ## Environmental data without geographic coordinates and altitude
# multi_data_trans_endemic_ming <-subset(multi_data_endemic_trans, select=-c(lat.x, lon.x, ele))
# 
# ## Community MATRIX with filtered samples -- Keep only the samples that do not show any NA
# asv_community_trim_endemic  <- comm_endemic_table_t%>% 
#   filter(sample %in% multi_data_endemic_trans$sample)
# ## We have to remove the last column, which contains a "sample"
# row.names(asv_community_trim_endemic) == multi_data_endemic_trans$sample
# colnames(asv_community_trim_endemic)[22774] ## sample
# Y.com_endemic <- asv_community_trim_endemic[,1:(dim(asv_community_trim_endemic)[2]-1)]
# ## Convert into numeric values
# Y.com_endemic <- sapply(Y.com_endemic, as.numeric)
# ## Saving Files 
# save(Y.com_endemic,file="comY_endemic_20232309.Rdata")
# save(multi_data_trans_endemic_ming, file="Env_endemic_20232309.Rdata")
# save(dist_geo_endemic,file= "Geo_endemic_20232309.Rdata")
# save(metadata_dist_endemic_df,file="geocoordinates_endemic_20232309.Rdata")
# 
