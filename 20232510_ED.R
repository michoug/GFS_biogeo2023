library(phyloseq)
library(phyloseqCompanion)
library(vegan)
library(tidyverse)

##Some metrics regarding the Number of ASVs found across the dataset
##Removing the Uganda sample!
NOMIS_R <- readRDS("20230825_NOMIS_rarefied.RData")
uganda=c("Uganda")
prune_Uganda <- subset_samples(NOMIS_R, !site_c %in% uganda)
prune_Uganda_df <- as.matrix(otu_table(prune_Uganda, taxa_are_rows=T))
prune_Uganda_df <- prune_Uganda_df[rowSums(prune_Uganda_df)>0,]


metadata_nomis <- sample.data.frame(prune_Uganda)
asv_df <- as.data.frame(t(otu_table(prune_Uganda_df, taxa_are_rows=T)))

## Here we are investigating the prevalence of ASV across the dataset
howmanyasv<- as.data.frame(colSums(asv_df != 0))
colnames(howmanyasv) <- c("Count_nb")
howmanyasv$ASV <- rownames(howmanyasv)
rownames(howmanyasv) <- NULL

## Melt asv table
asvdfmelt <- melt(as.matrix(asv_df))
## Keep only values that are > 0
asvdfmelt <- asvdfmelt[asvdfmelt$value >0,]

##### Endemism
## Liste des ASVs qui ont une prévalence de 1 : 11231 soit 28.8% -- comment sont-ils répartis par mountain ranges?
endemic_oneGFS <- howmanyasv %>% 
  filter(Count_nb == 1) 

endemic_oneGFS_n <- howmanyasv %>% 
  filter(Count_nb == 1) %>%
  summarize(number=n())
#n
#1 11231

## total number of endemic ASVs found in one single GFS
endemic_oneGFS_n/length(unique(row.names(prune_Uganda_df)))
# number
# 1 0.2880999

## On prend notre table d'ASV et on filtre pour avoir le nom des échantillons ## value c'est le nb de count, et count_nb c'est le nb d'échantillons
merge_equalone <- merge(endemic_oneGFS, asvdfmelt, by.x="ASV",by.y="Var2")
write.csv(merge_equalone, "unique_asv.csv")
## assign the name of the mountain range
merge_equalone$MR <- vapply((merge_equalone$Var1), function(x) metadata_nomis$site_c[metadata_nomis$sample == x], FUN.VALUE = character(1))
## rename column
colnames(merge_equalone) <- c("ASV","prev","glname","nb_count","MR")

## Mtn on va obtenir le nombre d'ASVs qui sont endemiques par région -- prenant en compte une prévalence de 1 (Observé dans 1 GFS)
endemic_MR_prev <- merge_equalone %>% group_by(MR) %>% summarize(prev=sum(prev)) 
endemic_oneGFS_plot <- merge_equalone %>% group_by(MR) ## here question, nb_count?
endemic_oneGFS_plot<- endemic_oneGFS_plot[c("MR","ASV","prev")]
endemic_oneGFS_plot$Color <- "C_uniquetoone"


## Calculate unique ASVs per mountain range. This means that these ASVs are found in only one GFS within their mountain range.
prop_unique_MR<- endemic_MR_prev%>%
  group_by(MR)%>% 
  summarize(prop=prev/endemic_oneGFS_n)
colnames(prop_unique_MR) <- c("mountain_range","prop_unique")

# 
# A tibble: 10 × 2
# MR           prop$n
# <chr>         <dbl>
# 1 Alaska       0.130 
# 2 Alps         0.0841
# 3 Caucasus     0.0705
# 4 Chile        0.138 
# 5 Ecuador      0.126 
# 6 Greenland    0.0212
# 7 Kirghizistan 0.0971
# 8 Nepal        0.126 
# 9 New_Zealand  0.118 
# 10 Norway      0.0897

## Barplot 
prop_unique_MR$mountain_range<- factor(prop_unique_MR$mountain_range, levels = prop_unique_MR$mountain_range[order(prop_unique_MR$prop_unique$n)])
p_unique<-ggplot(data=prop_unique_MR, aes(x=mountain_range, y=prop_unique$n)) +geom_bar(stat="identity")
p_unique + coord_flip() + theme_minimal()


## Filter the table of prevalence >1 & <10 GFS
twoandnine<- howmanyasv %>% filter(Count_nb > 1 & Count_nb < 10) ## 20468 ASVs 0.5250494 -- 52.5% des ASVs
twoandnine_n <- twoandnine %>% summarize(n())
twoandnine_n/length(unique(row.names(prune_Uganda_df)))

## On merge twoandnine avec la table d'ASV original pour obtenir le nom des échantillons puis pour chaque échantillon on lui attribue un nom de MR!
merge_twoandnine <- merge(twoandnine, asvdfmelt, by.x="ASV",by.y="Var2")
merge_twoandnine$MR <- vapply((merge_twoandnine$Var1), function(x) metadata_nomis$site_c[metadata_nomis$sample == x], FUN.VALUE = character(1))

colnames(merge_twoandnine) <- c("ASV","prev_GFS","glname","nb_count","MR")

## Maintenant on va compter le nombre de fois que l'ASV est présent dans la region, donc on va summer par GL code!(et on somme la value)
## en gros on compte le nombre de lignes par ASV et Glcode
twoandnine_end<- merge_twoandnine %>% group_by(MR,ASV) %>% summarize(prev=n()) 
twoandnine_end$Color <- "B_twoandnine"

## Maintenant on veut qu'il soit unique dans la mountain range! On  mutate pour pouvoir compter le nombre de fois ou l'ASV est observé (dans combien de MR), puis on filtre seulement les ASVs qui sont observés 1X dans la MR
## Après cela on va grouper par MR et on va faire la sum des N pour savoir combien il y a d'ASVs dans chaque MR
endemism_twoandnine <- twoandnine_end %>% group_by(ASV) %>% mutate(n=n()) %>% filter(n==1)%>%
  ungroup()%>% group_by(MR) %>%summarize(number=n()) ## 
# 
# A tibble: 10 × 2
# MR           number
# <chr>         <int>
# 1 Alaska          709
# 2 Alps            868
# 3 Caucasus        684
# 4 Chile          1009
# 5 Ecuador        1207
# 6 Greenland        77
# 7 Kirghizistan   1061
# 8 Nepal          1643
# 9 New_Zealand    3316
# 10 Norway         639

endemism_twoandnine_plot <- twoandnine_end %>% group_by(ASV) %>% mutate(n=n()) %>% filter(n==1)

# Number of ASVs that are endemic (present in 2-9 GFSs)
sum_endemic_twoandnine <- endemism_twoandnine %>% summarize(sum=sum(number))
# sum
# <int>
#   1 11213

# 11213/22773
# 49.2% are endemic 

sum_endemic_twoandnine/length(unique(row.names(prune_Uganda_df)))
# sum
#1 0.2876382

prop_unique_twonine<- endemism_twoandnine%>%
  group_by(MR)%>% 
  summarize(prop=number/twoandnine_n)

colnames(prop_unique_twonine)<-c("mountain_range", "prop_unique")

# # A tibble: 10 × 2
# MR           prop$`n()`
# <chr>             <dbl>
# 1 Alaska        0.0346 
# 2 Alps            0.0424 
# 3 Caucasus        0.0334 
# 4 Chile           0.0493 
# 5 Ecuador         0.0590 
# 6 Greenland       0.00376
# 7 Kirghizistan    0.0518 
# 8 Nepal           0.0803 
# 9 New_Zealand     0.162  
# 10 Norway          0.0312 


## We can do the same with the last section
tenandmore<- howmanyasv %>% filter(Count_nb >= 10) ## 7284 -> 18.7% of all ASVs
tenandmore_n <- tenandmore %>% summarize(n())
tenandmore_n/length(unique(row.names(prune_Uganda_df)))
# n()
# 1 0.1847137

## On merge tenandmore avec la table d'ASV original pour obtenir le nom des échantillons puis pour chaque échantillon on lui attribue un nom de MR!
merge_tenandmore <- merge(tenandmore, asvdfmelt, by.x="ASV",by.y="Var2")
merge_tenandmore$MR <- vapply((merge_tenandmore$Var1), function(x) metadata_nomis$site_c[metadata_nomis$sample == x], FUN.VALUE = character(1))

colnames(merge_tenandmore) <- c("ASV","prev_GFS","glname","nb_count","MR")

## Là on va sommer les ASVs (donner la prévalence de cet ASV dans cette chaine de montagne) lorsque le nombre de count est >0
tenandmore_end<- merge_tenandmore %>% group_by(MR,ASV) %>% summarize(prev=sum(nb_count>0))

tenandmore_end$Color <-"A_tenandemore"
## Maintenant on veut qu'il soit unique dans la mountain range! On  mutate pour pouvoir compte le nombre de fois ou l'ASV est observé, puis on filtre seulement les ASVs qui sont observés 1X dans la MR
## Après cela on va grouper par MR et on va faire la sum des N pour savoir combien il y a d'ASVs dans chaque MR
endemism_tenandmore <- tenandmore_end %>% group_by(ASV) %>% mutate(n=n()) %>% filter(n==1)%>%
  ungroup()%>% group_by(MR) %>%summarize(number=n()) ## 329/38983=0.8% des ASVs totaux sont endemiques ds cette section (>=10 GFSs).

endemism_tenandmore_plot <- tenandmore_end %>% group_by(ASV) %>% mutate(n=n()) %>% filter(n==1)

sum(endemism_tenandmore$number)
# [1] 329

# MR           number
# <chr>         <int>
# 1 Alaska            6
# 2 Alps             24
# 3 Caucasus         12
# 4 Kirghizistan     46
# 5 Nepal             6
# 6 New_Zealand     235

## Control Sanity check!
control_asv <- asvdfmelt
control_asv$MR <- vapply((control_asv$Var1), function(x) metadata_nomis$site_c[metadata_nomis$sample == x], FUN.VALUE = character(1))
## We want to know how many ASVs are endemic, and found in only one mountain range! We start from the original ASV table
## First we group by MR and ASV, and we count the nb of times the ASV appears for a given mountain ranges (in how many GFSs are they present?)
## Then we ungroup and we want to know how many ASVs are endemic so we just group by ASV and count the nb of times an ASV is observed throughout the dataset 
## meaning in how many MR does it appear? We would like only one MR, so we apply the filter(n==1)..
## and then we count the nb of lines -- meaning the nb of ASVs that are endemic (unique to one MR)
controlasv_end<- control_asv %>% group_by(MR,Var2) %>% summarize(prev=n()) %>%
  ungroup()%>% group_by(Var2) %>% mutate(n=n())%>% filter(n==1)%>%
  ungroup()%>% group_by(MR) %>%summarize(number=n())

# # A tibble: 10 × 2
# MR           number
# <chr>         <int>
# 1 Alaska        2179
# 2 Alps           1837
# 3 Caucasus       1488
# 4 Chile          2558
# 5 Ecuador        2619
# 6 Greenland       315
# 7 Kirghizistan   2197
# 8 Nepal          3062
# 9 New_Zealand    4872
# 10 Norway       1646


endemic_total<- control_asv %>% group_by(MR,Var2) %>% summarize(prev=n()) %>%
  ungroup()%>% group_by(Var2) %>% mutate(n=n())%>% filter(n==1)
write.csv(endemic_total,"endemic_list_2023.csv")
#22773 which should be the sum of 11213+329+11231=22773 its alright!!
##22773/38983=58.4%

## to get to know how many ASVs are present per mountain range (Observed)
controlasv_total<- control_asv %>% group_by(MR,Var2) %>% summarize(sumi=sum(n()))%>%
  ungroup()%>% group_by(MR)%>%summarize(sumii=sum(n()))

# MR           sumii
# <chr>        <int>
# 1 Alaska      8279
# 2 Alps         10971
# 3 Caucasus     10867
# 4 Chile         6980
# 5 Ecuador       5827
# 6 Greenland     6337
# 7 Kirghizistan  7519
# 8 Nepal        10204
# 9 New_Zealand   9771
# 10 Norway        6464

##Now we would like to know the proportion of endemic ASVs within mountain range --> what we report!!
endemic_mr <- merge(controlasv_total, controlasv_end, by="MR")
endemic_mr <- endemic_mr %>%
  group_by(MR) %>% summarize(proportion_endemic=number/sumii)

# A tibble: 10 × 2
# MR           proportion_endemic
# <chr>                     <dbl>
# 1 Alaska                  0.263 
# 2 Alps                     0.167 
# 3 Caucasus                 0.137 
# 4 Chile                    0.366 
# 5 Ecuador                  0.449 
# 6 Greenland                0.0497
# 7 Kirghizistan             0.292 
# 8 Nepal                    0.300 
# 9 New_Zealand              0.499 
# 10 Norway                   0.255 


##Proportion of unique ASVs -- 49.3% of endemic ASVs are unique to one GFS!
endemic_oneGFS_n/sum(controlasv_end$number)
# n
# 1 0.4931717

## % d'abondance relative des ASVS.
#1) prendre l'abondance relative de toute la table d'ASVs
#2) merger la table avec la liste des endemiques
#3) faire un object phyloseq...et filtrer

## Lets do some graphs
## We rbind the 3 different dataframes to plot the ASVs that are endemic to the different mountain ranges
endemic_oneGFS_plot$n <- 1
df_full <- rbind(endemic_oneGFS_plot, endemism_tenandmore_plot, endemism_twoandnine_plot)
niveaux <- c("New_Zealand","Nepal","Alps","Ecuador","Chile","Kirghizistan","Caucasus","Norway","Alaska","Greenland")
df_full$MR<- factor(df_full$MR, levels = niveaux)
df_full$MR <- fct_rev(df_full$MR)

prune_uganda_abondance= transform_sample_counts(prune_Uganda, function(x) x / sum(x))
prune_uganda_asv <- otu_table(prune_uganda_abondance, taxa_are_rows=T)
prune_uganda_asv <- prune_uganda_asv[rowSums(prune_uganda_asv[])>0,]

merge_endemic_abondance <- merge(df_full, prune_uganda_asv, by.x="ASV", by.y="row.names")

##remove the columns that we dont need anymore
merge_endemic_abondance<- as.data.frame(merge_endemic_abondance[,-c(2:5)])
##create phyloseq object
row.names(merge_endemic_abondance)<-merge_endemic_abondance$ASV
merge_endemic_abondance$ASV <- NULL
endemic_table_abondance <- otu_table(merge_endemic_abondance, taxa_are_rows=T)
merge_endemic_abondance_phylo <- merge_phyloseq(endemic_table_abondance, tax_NOMIS,metadata_nomis)

end_ab_phylo_table <- (as.matrix(otu_table(merge_endemic_abondance_phylo, taxa_are_rows=T)))

melt_asv <- melt(end_ab_phylo_table)
merge_asv_data <- merge(as.data.frame(melt_asv), as.matrix(metadata_nomis), by.x="Var2",by.y="sample")

##Here we would need to divide by the number of glaciers per mountain ranges. 
sum_mr <- merge_asv_data %>% group_by(site_c)%>% summarize(summrr=sum(value), n=n_distinct(Var2))%>% summarize(ar_mr=summrr/n, site_c)
## This is the proportion of endemic within the full asv table, not what we report! But we report the median 
# # A tibble: 10 × 2
# ar_mr site_c      
# <dbl> <chr>       
# 1 0.0920 Alaska      
# 2 0.0341 Alps        
# 3 0.0280 Caucasus    
# 4 0.275  Chile       
# 5 0.313  Ecuador     
# 6 0.0219 Greenland   
# 7 0.146  Kirghizistan
# 8 0.0954 Nepal       
# 9 0.291  New_Zealand 
# 10 0.0506 Norway     
# 
mean(sum_mr$ar_mr)
# [1] 0.1347194
sd(sum_mr$ar_mr)
# [1] 0.1157538

median_endemism_ab <- sum_mr %>% 
  summarise(median=median(ar_mr), x = quantile(ar_mr, c(0.25, 0.5, 0.75)))
# A tibble: 3 × 2
# median      x
# <dbl>  <dbl>
# 1 0.0937 0.0383
# 2 0.0937 0.0937
# 3 0.0937 0.243 

## Plot Stackplot
ggplot(df_full, aes(fill=Color, y=MR)) + 
  geom_bar(position="stack", stat="count")+
  scale_fill_manual(values = c("#D69C4E","#ECCBAE", "#046C9A"))+
  theme_minimal()

