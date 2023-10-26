## Beta diversity - Community Composition - NMDS and statistics
library(vegan)
library(ggplot2)
library(betadisper)
library(adonis)
library(ecodist)
library(phyloseq)

vegan_otu <- function(physeq){
  OTU <- otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

## be sure to compute this with the correct object - either prune_Uganda for statistical purposes, or NOMIS_rarefied to show the NMDS including 
## all samples

metadata_nmds <- sample.data.frame(prune_Uganda)

metadata_nmds$lat_attribute <- "A"
metadata_nmds$lat_attribute <- ifelse((metadata_nmds$lat) < 1, metadata_nmds$lat_attribute == "B", metadata_nmds$lat_attribute == "A")
metadata_nmds$lat_attribute <- as.factor(metadata_nmds$lat_attribute)

##NMDS plot
asv_table_nmds <- as.matrix((otu_table(prune_Uganda, taxa_are_rows=T)))
asv_table_nmds_f <- asv_table_nmds[rowSums(asv_table_nmds[])>0,]
nmds_bc_nomis <- metaMDS(t(log1p(asv_table_nmds_f)), distance = "bray", k = 2, trymax=999)

stressplot(nmds_bc_nomis)
data.scores = as.data.frame(scores(nmds_bc_nomis)$sites)
#add columns to data frame 
data.scores$Sample = metadata_nmds$Sample
data.scores$Site = metadata_nmds$site_c
head(data.scores)

##plot nmds without Uganda
##stress=0.1565905
##k=2

##plot nmds with Uganda
##stress=0.1561976
##k=2

nmds_bc_GFS_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 6, aes(colour = Site))+ 
 #stat_ellipse(aes(x=NMDS1, y=NMDS2,color=Site),type = "norm")+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Region", y = "NMDS2", shape = "Type") 

nmds_bc_GFS_plot

##Include Smooth line to show the latitude
##First use
colors<-c("#2E2A2BFF", "#CF4E9CFF","#8C57A2FF",
          "#3EBCB6","#82581FFF","#2F509EFF",
          "#E5614CFF","#97A1A7FF","#bee183","#DC9445FF","#EDD03E")

plot(x=data.scores$NMDS1, y=data.scores$NMDS2, type="n", xlim = c(-2.5, 2.5), ylim = c(-1.6, 1.6))
points(nmds_bc_nomis, display = "sites", cex = 2.3, pch=19, col=alpha(colors[factor(data.scores$Site)], 0.8))
ordisurf(nmds_bc_nomis, metadata_nmds$lat, add = TRUE, col="blue", labcex=1)

merge_data <- merge_data[merge_data$sample != "GL140", ]


climate_nomis <- read.csv("climate_nomis_2023.csv")
climate_nomis <- climate_nomis[climate_nomis$Sample %in% metadata_nmds$sample, ]


plot(x=data.scores$NMDS1, y=data.scores$NMDS2, type="n", xlim = c(-2, 2), ylim = c(-1.3, 1.3))
points(nmds_bc_nomis, display = "sites", cex = 2.3, pch=19, col=alpha(colors[factor(data.scores$Site)], 0.8))
ordisurf(nmds_bc_nomis, climate_nomis$pr, add = TRUE, col="blue", labcex=1)


## Statistical analyses
vegan_matrix_GFS<- vegan_otu(prune_Uganda)
GFS_bray <- vegdist(log1p(vegan_matrix_GFS), method="bray")
GFS_ado <- adonis2(GFS_bray ~ site_c, permutations = 999, method = "bray", data=metadata_nmds)
GFS_ado_latitude <- adonis2(GFS_bray ~ lat_attribute, permutations = 999, method = "bray", data=metadata_nmds)

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = GFS_bray ~ site_c, data = metadata_nmds, permutations = 999, method = "bray")
# Df SumOfSqs      R2      F Pr(>F)    
# site_c     9   18.271 0.33814 7.8336  0.001 ***
# Residual 138   35.763 0.66186                  
# Total    147   54.034 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = GFS_bray ~ lat_attribute, data = metadata_nmds, permutations = 999, method = "bray")
# Df SumOfSqs      R2      F Pr(>F)    
# lat_attribute   1    4.798 0.08879 14.226  0.001 ***
#   Residual      146   49.236 0.91121                  
# Total         147   54.034 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


library(pairwiseAdonis)
pairwise_GFS <- pairwise.adonis(GFS_bray, metadata_nmds$site_c, p.adjust.m="holm")

# pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1          New_Zealand vs Alps  1 3.7010688 14.196883 0.24394576   0.001      0.045   .
# 2         New_Zealand vs Nepal  1 3.0216817 11.739246 0.25665587   0.001      0.045   .
# 3  New_Zealand vs Kirghizistan  1 3.8971191 15.775461 0.31069064   0.001      0.045   .
# 4         New_Zealand vs Chile  1 1.4052987  4.887220 0.14860545   0.001      0.045   .
# 5        New_Zealand vs Alaska  1 2.9179171 11.394913 0.25667159   0.001      0.045   .
# 6     New_Zealand vs Greenland  1 1.8074541  7.087916 0.22799585   0.001      0.045   .
# 7      New_Zealand vs Caucasus  1 3.5055006 14.956796 0.28786986   0.001      0.045   .
# 8       New_Zealand vs Ecuador  1 2.1678221  8.603130 0.24163970   0.001      0.045   .
# 9        New_Zealand vs Norway  1 2.5938468 10.360769 0.27008763   0.001      0.045   .
# 10               Alps vs Nepal  1 1.5057511  5.645262 0.12367684   0.001      0.045   .
# 11        Alps vs Kirghizistan  1 2.6082074 10.122901 0.19801108   0.001      0.045   .
# 12               Alps vs Chile  1 1.7598314  6.002148 0.15004564   0.001      0.045   .
# 13              Alps vs Alaska  1 1.0050966  3.780780 0.08837567   0.001      0.045   .
# 14           Alps vs Greenland  1 0.9984483  3.726641 0.11049547   0.002      0.045   .
# 15            Alps vs Caucasus  1 0.9632885  3.911587 0.08338211   0.001      0.045   .
# 16             Alps vs Ecuador  1 2.3638205  8.944615 0.21324824   0.001      0.045   .
# 17              Alps vs Norway  1 1.6396696  6.244675 0.15516773   0.001      0.045   .
# 18       Nepal vs Kirghizistan  1 2.0219194  7.989941 0.20492314   0.001      0.045   .
# 19              Nepal vs Chile  1 1.4987439  4.961503 0.17131371   0.001      0.045   .
# 20             Nepal vs Alaska  1 1.2727947  4.825676 0.14266310   0.001      0.045   .
# 21          Nepal vs Greenland  1 1.5086932  5.673236 0.22097860   0.002      0.045   .
# 22           Nepal vs Caucasus  1 1.7597377  7.378456 0.18273249   0.001      0.045   .
# 23            Nepal vs Ecuador  1 2.4047408  9.215116 0.28604944   0.001      0.045   .
# 24             Nepal vs Norway  1 2.3096696  8.928574 0.27114974   0.001      0.045   .
# 25       Kirghizistan vs Chile  1 1.9216739  6.724332 0.21196134   0.001      0.045   .
# 26      Kirghizistan vs Alaska  1 2.2573283  8.977156 0.23031840   0.001      0.045   .
# 27   Kirghizistan vs Greenland  1 1.8014496  7.256494 0.25680801   0.001      0.045   .
# 28    Kirghizistan vs Caucasus  1 2.7601635 12.085587 0.26224223   0.001      0.045   .
# 29     Kirghizistan vs Ecuador  1 2.6576048 10.816684 0.31067531   0.001      0.045   .
# 30      Kirghizistan vs Norway  1 2.7421212 11.232534 0.31001237   0.001      0.045   .
# 31             Chile vs Alaska  1 1.4853508  4.916565 0.17611641   0.001      0.045   .
# 32          Chile vs Greenland  1 1.2100953  3.668296 0.20762026   0.001      0.045   .
# 33           Chile vs Caucasus  1 1.9226018  7.239852 0.21144520   0.001      0.045   .
# 34            Chile vs Ecuador  1 1.3165249  4.221487 0.19892513   0.001      0.045   .
# 35             Chile vs Norway  1 1.6243629  5.308311 0.22774326   0.001      0.045   .
# 36         Alaska vs Greenland  1 1.2748657  4.827682 0.20260812   0.001      0.045   .
# 37          Alaska vs Caucasus  1 1.6055550  6.787770 0.17499769   0.001      0.045   .
# 38           Alaska vs Ecuador  1 2.2822830  8.807631 0.28589120   0.001      0.045   .
# 39            Alaska vs Norway  1 1.8668571  7.268753 0.24014049   0.001      0.045   .
# 40       Greenland vs Caucasus  1 0.9505404  4.173085 0.15357419   0.001      0.045   .
# 41        Greenland vs Ecuador  1 1.2300924  4.744443 0.26737628   0.001      0.045   .
# 42         Greenland vs Norway  1 0.3540393  1.385707 0.09006457   0.082      0.082    
# 43         Caucasus vs Ecuador  1 2.4176600 10.614031 0.28988971   0.001      0.045   .
# 44          Caucasus vs Norway  1 1.6781060  7.392936 0.21495507   0.001      0.045   .
# 45           Ecuador vs Norway  1 1.5907424  6.347634 0.27187484   0.001      0.045   .

pairwise_GFS_lat <- pairwise.adonis(GFS_bray, metadata_nmds$lat_attribute, p.adjust.m="holm")
# pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
# 1 FALSE vs TRUE  1  4.797646 14.22639 0.08878929   0.001      0.001  **


