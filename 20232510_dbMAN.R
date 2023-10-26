ASVs<-read.csv("20230825_NOMIS_rarefied_deblur_table.csv", head=T, row.names = 1)
samp<-read.csv("samples.csv")

# indicator ASVs ----------------------------------------------------------
library(foreach)
library(doParallel)
library(adiv)


#remove Uganda
samp2<-samp[samp$region !="Uganda",]
asv_table<-ASVs[colnames(ASVs) %in% samp2$ID]
asv_table<-asv_table[rowSums(asv_table)>0,]
asv_table<-asv_table[rowSums(asv_table>0)>1,] #remove unique ASVs


#hellinger-transformation
asv_table<-decostand(asv_table, method="hellinger", MARGIN=2)

n.cores<-detectCores()
registerDoParallel(cl <- makeCluster(n.cores-1))

res.f.all<-foreach(f = levels(factor(samp2$region)), .packages = c("adiv")) %dopar% {
  samp2.f<-samp2
  samp2.f$region[!grepl(f, samp2$region)]<-"other"
  Q.f<-dbMANOVAspecies(t(asv_table),samp2.f$region, nrep=999, global=FALSE,species=TRUE, padj="BH")
  Q.f.pairwise_adj <- dbMANOVAspecies_pairwise(Q.f)
  c(summary(Q.f.pairwise_adj))
}
save.image("db.MANOVA.res.RData")

# indicator NMDS

ind.ASVs<-read.csv("indicator_taxa.csv", head=T)

ind.ASVs.long<-melt(ind.ASVs, id.vars="X")
ind.ASVs.long<-na.omit(ind.ASVs.long)
ind.ASVs.long<-ind.ASVs.long[order(ind.ASVs.long$X),]

cols<-brewer.pal(10, "Spectral")

nMDS<-metaMDS(t(ASVs))
plot(nMDS, type="n")
ordispider(nMDS, groups=samp$region, col=cols, lwd=2)
points(nMDS, disp="sites", pch=21, bg=cols[factor(samp$region)], cex=1.5)

spec.scores<-scores(nMDS, display="species")
spec.scores.db.MANOVA<-data.frame(spec.scores[rownames(spec.scores) %in% ind.ASVs$X,])
spec.scores.db.MANOVA<-spec.scores.db.MANOVA[order(row.names(spec.scores.db.MANOVA)),]
summary(ind.ASVs.long$X==rownames(spec.scores.db.MANOVA))

spec.scores.db.MANOVA<-data.frame(spec.scores[rownames(spec.scores) %in% X$ASV,])
points(spec.scores.db.MANOVA$NMDS1,spec.scores.db.MANOVA$NMDS2, pch=16, cex=0.2, col=cols[factor(ind.ASVs.long$variable)])


