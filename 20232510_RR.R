## Radar per region
#Computing radarplots from environmental dataset
require(ggplot2)
require(reshape)
require(scales)
library(scales)
library(ggthemes)
library(ggpubr)
library(phyloseq)
library(phyloseqCompanion)
library(fmsb)
library(ggvanced)
library(fishualize)


prune_Uganda_df <- as.matrix(otu_table(prune_Uganda_df, taxa_are_rows=T))
metadata_radar <- sample.data.frame(prune_Uganda)

## Merging the different environmental datasets - nomis db, minerals, chl and climatic variables
metadata_nomis_db <- read.csv("nomis_metadata_database_2023.csv")
minerals_nomis <- read.csv("minerals_nomis_2023.csv")
chla_nomis <- read.csv("chla_nomis_2023.csv")

#lets take the mean for chla values per patch!
chla_mean <- chla_nomis %>% 
  group_by(gl_code)%>% summarize(mean_chla=mean(chla))

chla_mean_df <- as.data.frame(chla_mean)

metadata_radar_db <- merge(metadata_radar, metadata_nomis_db, by="sample", no.dups = T)
metadata_radar_db <- merge(metadata_radar_db, minerals_nomis, by="sample", no.dups=T)
metadata_radar_db <- merge(metadata_radar_db, chla_mean_df, by.x="sample", by.y="gl_code",no.dups=T)

## replace NAs by 0 in nutrient colums so that we can sum NH4, NO3 and NO2 together!
metadata_radar_db_sub <- metadata_radar_db  %>% 
  mutate(no2 = coalesce(no2,0),no3=coalesce(no3,0),nh4=coalesce(nh4,0))

#Sum inorganic nitrogen into DIN
rowsum_nut <- rowSums(metadata_radar_db_sub[, c("no2", "no3", "nh4")])
metadata_radar_db_sub$DIN <- rowsum_nut 

## Convert col to numeric
convert_columns_to_numeric <- function(data, columns) {
  for (col in columns) {
    data[[col]] <- as.numeric(data[[col]])
  }
  return(data)
}

metadata_radar <- convert_columns_to_numeric(metadata_radar_db_sub, c('DIN', 'srp', 'turb', 'doc', 'mean_chla', 'Clays', 'Feldspar', 'Quartz', 'Calcite','water_temp','pH'))
metadata_radar_subset <- metadata_radar[c('sample','site_c', 'lat.x','lon.x','water_temp', 'pH', 'cond', 'turb', 'gl_sa.x', 'gl_cov.x', 'srp', 'DIN', 'doc', 'mean_chla', 'Clays', 'Feldspar', 'Quartz', 'Calcite','sn_sp_dist')]

#remove duplicate lines
metadata_radar_sub_u <- distinct(metadata_radar_subset)
metadata_radar_sub_u[is.na(metadata_radar_sub_u)] <- 0

metadata_radar_sub_u <- metadata_radar_sub_u[c('site_c','water_temp', 'pH', 'cond', 'turb', 'gl_sa.x', 'gl_cov.x', 'srp', 'DIN','doc')]
metadata_radar_sub_u_table <- metadata_radar_sub_u[c('site_c','water_temp', 'pH', 'cond', 'turb', 'gl_sa.x', 'gl_cov.x', 'srp', 'DIN','doc','Feldspar','Calcite','Quartz','Clays','mean_chla','sn_sp_dist')]

metadata_median <- as.data.frame(metadata_radar_sub_u %>% group_by(site_c) %>% summarize_if(is.numeric, list("q25" = ~quantile(., 0.25,type=2),"q75" = ~quantile(., 0.75,type=2),"q50" = ~quantile(., 0.50,type=2))))
median_only <- as.data.frame(metadata_radar_sub_u %>% group_by(site_c) %>% summarize_if(is.numeric, list("q50" = ~quantile(., 0.50,type=2))))
metadata_median_table<- as.data.frame(metadata_radar_sub_u_table %>% group_by(site_c) %>% summarize_if(is.numeric, list("q75" = ~quantile(., 0.75,type=2))))


process_site_data <- function(data, site_name) {
  # Choose a region
  metadata_site <- subset(data, site_c == site_name)
  melt_site <- melt(metadata_site)
  
  # Filter data for q25, q50, and q75
  q25_data <- subset(melt_site, grepl("_q25$", variable))
  q50_data <- subset(melt_site, grepl("_q50$", variable))
  q75_data <- subset(melt_site, grepl("_q75$", variable))
  
  # Pivot the data
  q25_pivot <- pivot_wider(q25_data, names_from = variable, values_from = value)
  q50_pivot <- pivot_wider(q50_data, names_from = variable, values_from = value)
  q75_pivot <- pivot_wider(q75_data, names_from = variable, values_from = value)
  
  # Rename the columns
  # colnames(q25_pivot) <- c("group","water_temp","pH","cond","turb","gl_sa","gl_cov","srp","DIN")
  # colnames(q50_pivot) <- c("group","water_temp","pH","cond","turb","gl_sa","gl_cov","srp","DIN")
  # colnames(q75_pivot) <- c("group","water_temp","pH","cond","turb","gl_sa","gl_cov","srp","DIN")
  # 
  colnames(q25_pivot) <- c("group","water_temp","pH","cond","turb","gl_sa","gl_cov","srp","DIN","Feldspar","Calcite","Quartz","Clays","mean_chla","sn_sp_dist")
  colnames(q50_pivot) <- c("group","water_temp","pH","cond","turb","gl_sa","gl_cov","srp","DIN","Feldspar","Calcite","Quartz","Clays","mean_chla","sn_sp_dist")
  colnames(q75_pivot) <- c("group","water_temp","pH","cond","turb","gl_sa","gl_cov","srp","DIN","Feldspar","Calcite","Quartz","Clays","mean_chla","sn_sp_dist")
  
  
  # Combine q25, q50, and q75 data
  min <- data.frame(matrix(rep(c(0), 9), nrow = 1))
  colnames(min)<-c("group","water_temp","pH","cond","turb","gl_sa","gl_cov","srp","DIN","Feldspar","Calcite","Quartz","Clays","mean_chla","sn_sp_dist")
  #colnames(min)<-c("group","water_temp","pH","cond","turb","gl_sa","gl_cov","srp","DIN")
  
  combined_data <- rbind(min, q75_pivot, q50_pivot)
  
  # Add a new column "group1"
  combined_data$group1 <- c("min", "max", "Q50")
  combined_data$group <- NULL
  
  combined_data <- combined_data %>%
    select(group1, everything())
  
  return(combined_data)
}



generate_and_save_spider_plot <- function(data, site_name) {
  gg <- ggspider(data, axis_name_offset = 0.1, background_color = "white", fill_opacity = 0.4) +
    scale_fill_manual(values = c("Q25" = "black", "Q50" = "#3F459BFF", "Q75" = "white")) +
    scale_color_manual(values = c("Q25" = "black", "Q50" = "#3F459BFF", "Q75" = "#3F459BFF"))
  
  # Start a PDF device
  pdf(paste0(site_name, "_spider_plot.pdf"), width = 8, height = 8)
  
  # Print the plot
  print(gg)
  
  # Close the PDF device
  dev.off()
}

##Compute dataframe
alaska_data <- process_site_data(metadata_median, "Alaska")
caucasus_data <- process_site_data(metadata_median, "Caucasus")
norway_data <- process_site_data(metadata_median, "Norway")
nepal_data <- process_site_data(metadata_median, "Nepal")
alps_data <- process_site_data(metadata_median, "Alps")
greenland_data <- process_site_data(metadata_median, "Greenland")
ecuador_data <- process_site_data(metadata_median, "Ecuador")
chile_data <- process_site_data(metadata_median, "Chile")
kh_data <- process_site_data(metadata_median, "Kirghizistan")
nz_data <- process_site_data(metadata_median, "New_Zealand")

# Generate and save the plots!
generate_and_save_spider_plot(alaska_data, "Alaska")
generate_and_save_spider_plot(caucasus_data, "Caucasus")
generate_and_save_spider_plot(norway_data, "Norway")
generate_and_save_spider_plot(nepal_data, "Nepal")
generate_and_save_spider_plot(alps_data, "Alps")
generate_and_save_spider_plot(greenland_data, "Greenland")
generate_and_save_spider_plot(ecuador_data, "Ecuador")
generate_and_save_spider_plot(chile_data, "Chile")
generate_and_save_spider_plot(nz_data, "New_Zealand")
generate_and_save_spider_plot(kh_data, "Kirghizistan")

#### This is what we keep for the paper
## We keep the function create beautiful radarchart and re-scale btw 0-100
print(median_only)

#         site_c  water_temp_q50 pH_q50 cond_q50  turb_q50 gl_sa.x_q50 gl_cov.x_q50 srp_q50 DIN_q50  doc_q50
# 1        Alaska           0.70 8.8430    91.50 134.30000       4.910        0.530   3.540  29.880 151.5567
# 2          Alps           0.45 7.6950    63.25  69.15000       3.175        0.570   2.015 182.555 158.1667
# 3      Caucasus           1.00 7.8300    81.60  30.30000       3.040        0.510   1.160 263.030 184.4333
# 4         Chile           2.80 7.8935    81.55 127.49500       6.925        0.580   4.335  49.595 108.7500
# 5       Ecuador           0.90 5.0800     6.10 759.70000       0.420        0.880  16.440  38.260 262.7667
# 6     Greenland           0.35 6.0300     4.70  69.75000       0.805        0.650   5.020  16.810 126.0000
# 7  Kirghizistan           0.40 8.5520    86.00  23.76585       1.480        0.620   1.080 447.220 249.7767
# 8         Nepal           0.70 8.3130   117.05  19.54256       4.460        0.455   1.150 117.110 109.1667
# 9   New_Zealand           2.95 7.7050    33.05   8.85000       0.710        0.540   6.800  46.400  91.5500
# 10       Norway           0.80 7.1450     6.25  32.90000       3.765        0.665   2.130  59.735 108.2083


# Function to scale values between 0 and 100 for each column !!!
scale_column_to_0_100 <- function(x) {
  rescaled <- scale(x, center = FALSE, scale = max(x))
  scaled <- rescaled * 100
  return(scaled)
}

# Apply scaling to each numeric column
for (col_name in colnames(median_only)[-1]) {
  median_only[, col_name] <- scale_column_to_0_100(median_only[, col_name])
}

median_only_rescaled <-median_only

## Adding two rows, one for the min=0 and one for the max=100
## Now we need to include two additional rows max and min to each of the variable
min <- data.frame(
  site_c = "Min",
  water_temp_q50 = 0,
  pH_q50 = 0,
  cond_q50 = 0,
  turb_q50 = 0,
  gl_sa.x_q50 = 0,
  gl_cov.x_q50 =0,
  srp_q50 = 0,
  DIN_q50 = 0,
  doc_q50=0
)

max <- data.frame(
  site_c = "max",
  water_temp_q50 = 100,
  pH_q50 = 100,
  cond_q50 = 100,
  turb_q50 = 100,
  gl_sa.x_q50 = 100,
  gl_cov.x_q50 =100,
  srp_q50 = 100,
  DIN_q50 = 100,
  doc_q50=100
)

## Now bind the different columns together!
bind_columns <- rbind(max,min,median_only_rescaled)
rownames(bind_columns) <- bind_columns$site_c
bind_columns$site_c <- NULL


# Split the screen in 4 parts
colors <- c("#2E2A2BFF", "#CF4E9CFF","#8C57A2FF",
            "#3EBCB6","#82581FFF","#2F509EFF",
            "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183")
titles <- c("Alaska", "Alps", "Caucasus","Chile","Ecuador","Greenland","Kirghizistan","Nepal","New_Zealand","Norway")

op <- par(mar = c(1, 1, 1, 1))# Adjust the margins as needed
par(mfrow = c(1, 2))# Create a single row with 2 columns

# Create the radar chart

create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}

op <- par(mar = c(1, 1, 1, 1))# Adjust the margins as needed
par(mfrow = c(1, 2))# Create a single row with 2 columns

for(i in 9:10){
  create_beautiful_radarchart(
    data = bind_columns[c(1, 2, i+2), ], caxislabels =  c(0, 25, 50, 75, 100),
    color = colors[i], title = titles[i]
  )
}
par(op)

