library(phyloseqCompanion)

## Rank abundance Curves, all Glaciers! So we have to take into account the biogeo betadiv 
comm<-as(otu_table(prune_Uganda), "matrix")
env<-as(sample.data.frame(prune_Uganda), "matrix")

comm<-as.data.frame(comm)
env<-as.data.frame(env)

### Applying this to the full dataset
comm_melt <- reshape2::melt(as.matrix(comm), na.rm=TRUE)
colnames(comm_melt)<-c("Var1","sample","Abundance")
head(comm_melt)
comm_melt<-as.data.frame(comm_melt)
comm_melt_scaling<-comm_melt[comm_melt$Abundance!=0,] ##Removing taxa with 0 abundance
head(comm_melt_scaling)

comm_melt_RA<-comm_melt_scaling %>%
  group_by(sample) %>%
  mutate(Scale_RA=Abundance/sum(Abundance)) ##Takes the relative abundance for each sample
head(comm_melt_RA)

##Sanity check
comm_melt_RA %>%
  group_by(sample)%>% summarize(somme=sum(Scale_RA))

comm_melt_RANK<-comm_melt_RA %>%
  group_by(sample) %>%
  mutate(my_ranks = order(order(Scale_RA,decreasing=T))) #Ranks

## Pour comprendre comment order fonctionne, en fait -- on prendre notre vecteur
## et on regarde à quelle position les chiffres/nombres sont placés et on les replace
## en ordre!
Rank_merge<-merge(comm_melt_RANK, env, by="sample")

Rank_merge_mod <- expand.grid(site_c = unique(Rank_merge$site_c),
                              Rank = 1:max(Rank_merge$my_ranks))

##Adding the average for each mountain range
Rank_merge_mod$Scale_RA = vapply(1:nrow(Rank_merge_mod), 
                                 function(i) mean(Rank_merge$Scale_RA[(Rank_merge$my_ranks == Rank_merge_mod$Rank[i]) & 
                                                                 (Rank_merge$site_c == Rank_merge_mod$site_c[i])]),FUN.VALUE = numeric(1))
colors<-c("#2E2A2BFF", "#CF4E9CFF", "#8C57A2FF",
          "#3EBCB6","#82581FFF","#2F509EFF",
          "#E5614CFF","#97A1A7FF","#bee183","#DC9445FF","#EDD03E")
rare_curve <- ggplot() +
  geom_line(data=Rank_merge, aes(x=my_ranks, y=log(Scale_RA), group=sample),size=0.2,colour="black",alpha=0.2)+
  geom_line(data=Rank_merge_mod, aes(x=Rank, y=log(Scale_RA), colour=site_c))+
  scale_color_manual(values= c("Alaska"="#2E2A2BFF","Alps"="#CF4E9CFF","Caucasus"="#8C57A2FF",
                               "Chile"="#3EBCB6","Ecuador"="#82581FFF","Greenland"="#2F509EFF",
                               "Kirghizistan"="#E5614CFF","Nepal"="#97A1A7FF","New_Zealand"="#bee183","Norway"="#DC9445FF"))+
  xlab("Rank")+
  facet_wrap(~site_c, nrow=3)+
  theme(panel.background = element_rect(fill = 'white', color="grey60"))

## Frequence relative des ASVs
taxa_names(prune_Uganda) <- paste("ASV", 1:ntaxa(prune_Uganda), sep="_")
sort(sample_sums(prune_Uganda))
comm_count <- as.matrix(t(otu_table(prune_Uganda, taxa_are_rows=T)))

spe_table <- data.frame(species=colnames(comm_count),
                        freq.rel=apply(comm_count,2,sum)/sum(comm_count)*100)

spe_table_10 <- subset(spe_table[1:15000,])

setwd(here("Figures"))
pdf("Species_regional_occurence.pdf", width=8, height=15)
ggplot(spe_table, aes(x=reorder(species, freq.rel), y=freq.rel)) +
  geom_bar(position="dodge", stat="identity", fill="steelblue2", colour="grey") +
  xlab("") + ylab("Relative species frequency") +
  coord_flip()+
  theme(axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()) #remove y axis ticks
dev.off()
