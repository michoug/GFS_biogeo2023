library(vegan)
library(RColorBrewer)
library(reshape2)

ASVs<-read.csv("20230825_NOMIS_rarefied_deblur_table.csv", head=T, row.names = 1)
ASVs_uganda <- ASVs %>% select(-GL140)
summary(colnames(ASVs_uganda)==metadata_glaciers$sample)

ind.ASVs<-read.csv("indicataor_taxa_2023_NMDS.csv", head=T)

ind.ASVs.long<-melt(ind.ASVs, id.vars="X")
ind.ASVs.long<-na.omit(ind.ASVs.long)
ind.ASVs.long<-ind.ASVs.long[order(ind.ASVs.long$X),]

cols=c("#2E2A2BFF","#CF4E9CFF","#8C57A2FF",
  "#3EBCB6","#82581FFF","#2F509EFF",
"#E5614CFF","#97A1A7FF","#bee183","#DC9445FF","#EDD03E","#000000")

nMDS<-metaMDS(t(log1p(ASVs_uganda)), distance = "bray", k = 2, trymax=999)
plot(nMDS, type="n")
ordispider(nMDS, groups=metadata_glaciers$site_c, lwd=2)
points(nMDS, disp="sites", pch=21, bg=cols[factor(metadata_glaciers$site_c)], cex=1.5)

spec.scores<-scores(nMDS, display="species")
spec.scores.db.MANOVA<-data.frame(spec.scores[rownames(spec.scores) %in% ind.ASVs$X,])
spec.scores.db.MANOVA<-spec.scores.db.MANOVA[order(row.names(spec.scores.db.MANOVA)),]
summary(ind.ASVs.long$X==rownames(spec.scores.db.MANOVA))

#spec.scores.db.MANOVA<-data.frame(spec.scores[rownames(spec.scores) %in% X$ASV,])
points(spec.scores.db.MANOVA$NMDS1,spec.scores.db.MANOVA$NMDS2, pch=16, cex=0.2, col=cols[factor(ind.ASVs.long$variable)])


