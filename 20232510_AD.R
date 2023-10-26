###Computing the average ASV richness from rarefied data
library(dplyr)
library(tidyverse)
library(breakaway)
library(data.table)
library(phyloseq)
library(microViz)
library(viridis)
library(hrbrthemes)
library(paletteer)

NOMIS_R <- readRDS("20230825_NOMIS_rarefied.RData")

uganda=c("Uganda")
prune_Uganda <- subset_samples(NOMIS_R, !site_c %in% uganda)
prune_Uganda_df <- as.matrix(otu_table(prune_Uganda, taxa_are_rows=T))
prune_Uganda_df <- prune_Uganda_df[rowSums(prune_Uganda_df[])>0,]


## Calculate diversity metrics
diversity_nomis <- plot_richness(prune_Uganda, measures=c("Observed"), color="site_c")
alphadt_nomis<- data.table(diversity_nomis$data)

## Compute the average observed richness per mountain range
## Filter asv richness
OR <- subset(alphadt_nomis, variable == "Observed")
hist(OR$value) ##normally distributed
ASVrichness_GFS <- OR %>% 
  group_by(site_c) %>% 
  summarise(average=mean(value), std=sd(value))
write.csv(alphadt_nomis, "table_diversite_nomis2023.csv")
# 
# site_c       average   std
# <chr>          <dbl> <dbl>
# 1 Alaska         1641.  389.
# 2 Alps           1894.  609.
# 3 Caucasus       2489.  500.
# 4 Chile          1326.  546.
# 5 Ecuador        1466.  637.
# 6 Greenland      1678. 1044.
# 7 Kirghizistan   1454.  544.
# 8 Nepal          2129.  463.
# 9 New_Zealand    2049.  400.
# 10 Norway        1585.  490.

average_richness <- alphadt_nomis %>% 
  summarise(average=mean(value), std=sd(value))
# average      std
# 1 1846.561 630.8855

median_richness <- alphadt_nomis %>% 
  summarise(median=median(value), x = quantile(value, c(0.25, 0.5, 0.75)))
# 
# > median_richness
# median       x
# 1   1801 1359.25
# 2   1801 1801.00
# 3   1801 2260.75

## Plot ASV richness per region using violin plot and include jitter 
plot_div <- ggplot(OR,aes(x=alphadt_nomis$site_c,y=value, color=site_c)) + 
  geom_violin(width=1)  + 
  geom_boxplot(width=0.1, color="black", alpha=0.2) +
  geom_jitter(aes(group=site_c), position=position_jitterdodge()) + 
  theme(axis.text.x = element_text(angle = 45, margin=margin(0.5, unit="cm"))) + theme_bw()

plot_div + scale_colour_manual(values=c("#2E2A2BFF","#CF4E9CFF","#8C57A2FF",
                                        "#3EBCB6","#82581FFF","#2F509EFF",
                                        "#E5614CFF","#97A1A7FF","#bee183","#DC9445FF","#EDD03E"))


ggplot( aes(x=myaxis, y=value, fill=name)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="black", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A Violin wrapping a boxplot") +
  xlab("")


Shannon <- subset(alphadt_nomis, variable == "Shannon")
Shannon_GFS <- Shannon %>% 
  group_by(site_c) %>% 
  summarise(average=mean(value), std=sd(value))

# # A tibble: 10 × 3
# site_c       average   std
# <chr>          <dbl> <dbl>
# 1 Alaska          5.19 0.816
# 2 Alps            5.53 0.565
# 3 Caucasus        5.52 0.478
# 4 Chile           4.92 0.947
# 5 Ecuador         5.27 0.678
# 6 Greenland       5.44 0.832
# 7 Kirghizistan    5.29 0.561
# 8 Nepal           5.78 0.446
# 9 New_Zealand     5.44 0.579
# 10 Norway         5.28 0.724

median_shannon <- Shannon %>% 
  summarise(median=median(value), x = quantile(value, c(0.25, 0.5, 0.75)))
# median        x
# 1 5.47767 4.986756
# 2 5.47767 5.477670
# 3 5.47767 5.870378

## Plot ASV richness per region using violin plot and include jitter 
plot_shannon <- ggplot(Shannon,aes(x=site_c,y=value, color=site_c)) + 
  geom_violin(scale="width") + 
  geom_jitter(aes(group=site_c), position=position_jitterdodge()) + 
  stat_summary(fun.y="mean",geom="crossbar", mapping=aes(ymin=..y.., ymax=..y..), 
               width=1, position=position_dodge(),show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 45, margin=margin(0.5, unit="cm"))) + theme_bw()
plot_shannon + scale_colour_manual(values=c("#2E2A2BFF","#CF4E9CFF","#8C57A2FF",
                                        "#3EBCB6","#82581FFF","#2F509EFF",
                                        "#E5614CFF","#97A1A7FF","#bee183","#DC9445FF","#EDD03E"))

# Evenness with Bulla index
bulla_estimate <- microbiome::evenness(prune_Uganda, index="all")
alphadt_nomis$even <-bulla_estimate$bulla
hist(alphadt_nomis$even)
# Compute the average evenness per mountain range
evenness_GFS <- alphadt_nomis %>% 
  group_by(site_c) %>% 
  summarise(average=mean(even), std=sd(even))

# # A tibble: 10 × 3
# site_c       average    std
# <chr>          <dbl>  <dbl>
# 1 Alaska         0.290 0.0748
# 2 Alps           0.300 0.0391
# 3 Caucasus       0.280 0.0417
# 4 Chile          0.289 0.0803
# 5 Ecuador        0.318 0.0525
# 6 Greenland      0.314 0.0747
# 7 Kirghizistan   0.293 0.0512
# 8 Nepal          0.325 0.0504
# 9 New_Zealand    0.291 0.0471
# 10 Norway        0.291 0.0684


evenness_GFS_median <- alphadt_nomis %>% 
  #summarise(average=mean(even), std=sd(even)) %>% 
  summarise(median=median(even), x = quantile(even, c(0.25, 0.5, 0.75)))

# > evenness_GFS_median
# median         x
# 1 0.3021268 0.2561128
# 2 0.3021268 0.3021268
# 3 0.3021268 0.3340257

evenness_GFS_total <- alphadt_nomis %>% 
  summarise(average=mean(even), std=sd(even))

# ##lets do some stat
# library(performance)
# library(multcomp)
# region_effect_anova_even <-aov(even ~ site_c, data = alphadt_nomis)
# summary(region_effect_anova_even)
# aov_region_even <-anova(region_effect_anova_even)
# check_model(region_effect_anova_even)
# check_normality(region_effect_anova_even)
# check_heteroscedasticity(region_effect_anova_even)
# 
# region_effect_kruskal <-kruskal.test(even ~ site_c, data = alphadt_nomis)
# summary(region_effect_kruskal)
# # Kruskal-Wallis rank sum test
# # 
# # data:  even by site_c
# # Kruskal-Wallis chi-squared = 7.7787, df = 9, p-value = 0.5566
# 
# ## Richness
# region_effect_anova_richness <-aov(value ~ site_c, data = OR)
# summary(region_effect_anova_richness)
# aov_region_richness <-anova(region_effect_anova_richness)
# check_model(region_effect_anova_richness)
# check_normality(region_effect_anova_richness)
# check_heteroscedasticity(region_effect_anova_richness)
# 
# tukey_richness <- TukeyHSD(region_effect_anova_richness)
# ## Shannon
# region_effect_anova_shannon <-aov(sqrt(value) ~ site_c, data = Shannon)
# summary(region_effect_anova_shannon)
# aov_region_shannon <-anova(region_effect_anova_shannon)
# check_model(region_effect_anova_shannon)
# check_normality(region_effect_anova_shannon)
# check_heteroscedasticity(region_effect_anova_shannon)
# ##Ok we have to use non-parametric test
# 
# kruskal_shannon <-kruskal.test(value ~ site_c, data = Shannon)
# summary(kruskal_shannon)
# library(FSA)
# ##Posthoc
# DT_shannon = dunnTest(value ~ site_c, data = Shannon,
#               method="holm")     



##Results from INEXT - observed and estimated richness
results_inext = read.csv("results_INEXT_2022111.csv",sep=",",header=T)
results_inext <- as.data.frame(results_inext)

names(results_inext) <- c("Site","diversity","asv_obs","estimated","SE","LCL","UCL")

inext_res <- results_inext |> pivot_longer(cols = -c(Site, diversity, LCL, UCL), names_to = "Type") 
head(inext_res)

merged_pivot <- merge(results_inext, inext_res, by.x="Site",by.y="Site")

##Assign SD of Percent to 0
merged_pivot$SE[merged_pivot$Type == 'Observed'] = 0

merged_pivot <- as.data.frame(merged_pivot)
merged_pivot$SE <- as.numeric(merged_pivot$SE)
merged_pivot$value <- as.numeric(merged_pivot$value)

##Plot everything!
regional_div <- ggplot(merged_pivot, aes(x=Site, y = value, fill= Type)) + 
  geom_col(position="dodge") + 
  labs(y="Average")+ geom_errorbar(aes(ymin=(value - SE), ymax=(value + SE)),
                                   width=.2, colour="red", 
                                   position=position_dodge(.9))

regional_div + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

##With Facet
facet_regional <- ggplot(merged_pivot, aes( x = Site, y = value, fill=Type) ) + 
  geom_col(position="dodge") + labs(y="Number of ASVs")+ 
  geom_errorbar(aes(ymin=(value - SE), ymax=(value + SE)), width=.2, colour="red", 
                position=position_dodge(.9))+
  facet_wrap( ~ Site ) + 
  xlab("Mountain Ranges")  
p1 + theme_bw() + 
  theme( axis.text.x = element_text( angle = 90,  hjust = 1 ) )


## Effect of altitude on species richness
## On check comment la latitude et l'altitude sont lilées entre elles 
## de la latitude

lat_ele <- ggplot(data=alphadt_nomis,aes(x=lat,y=ele_sp))+
   geom_point()+
   geom_smooth(method="gam", formula = y ~ s(x, bs = 'tp'))

 lat_ele + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# We test the effect of latitude on altitude via a spline
 library(mgcv)
 mod_lat <- gam(data=alphadt_nomis, formula = ele_sp ~ s(lat, bs='tp'))
 summary(mod_lat)
 alphadt_nomis$ele_resids = mod_lat$residuals
 # 
 # Family: gaussian 
 # Link function: identity 
 # 
 # Formula:
 #   ele_sp ~ s(lat, bs = "tp")
 # 
 # Parametric coefficients:
 #   Estimate Std. Error t value Pr(>|t|)    
 # (Intercept)   2453.5       36.3   67.59   <2e-16 ***
 #   ---
 #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 # 
 # Approximate significance of smooth terms:
 #   edf Ref.df   F p-value    
 # s(lat) 7.517  8.202 168  <2e-16 ***
 #   ---
 #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 # 
 # R-sq.(adj) =  0.903   Deviance explained = 90.8%
 # GCV = 2.069e+05  Scale est. = 1.9499e+05  n = 148


# ## We try now the inverse and investigate how altitude influences latitude
# mod_ele <- gam(data=alphadt_nomis, formula = lat_sp ~ s(ele_sp, bs='tp'))
# summary(mod_ele)
# alphadt_nomis$lat_resids = mod_ele$residuals

# ## On detrend en fonction de la latitude -- Robin advice
# ele_sp <- lm(data=alphadt_nomis, formula = value ~ ele_sp + lat_sp)
# summary(ele_sp)
# check_model(ele_sp)
# check_normality(ele_sp)
# check_heteroscedasticity(ele_sp)

# null_model_elevation <- glm(value ~ 1, data=alphadt_nomis)
# anova(ele_sp, null_model_elevation, test="F")
# 
# spe_ele_plot <- ggplot(alphadt_nomis, aes(x=ele_sp, y=value)) + geom_point() + geom_smooth(method='lm', color="black")
# spe_ele_plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 
# ele_sp$coefficients
# (Intercept)        ele_sp   abs(lat_sp) 
# 1921.50376520   -0.01085261   -7.87143852 

# ## We retrieve the residuals of the models and investigate the correlations between altitude or latitude and species richness
pred <- predict(lm(data=OR, formula = value ~ ele_resids), 
                 se.fit = TRUE, interval = "confidence")
limits <- as.data.frame(pred$fit)
spe_ele <- ggplot(alphadt_nomis, aes(x=ele_resids, y=value)) + geom_point(size=3, alpha=0.4,color="#3F459BFF") + 
geom_smooth(method='gam', se=T, formula = y ~ s(x, bs='tp'), fill="blue", color="black", alpha=0.2, span=0.3)  +
geom_line(aes(x = ele_resids, y = limits$lwr), 
          linetype = 2) +
  geom_line(aes(x = ele_resids, y = limits$upr), 
            linetype = 2) 
spe_ele + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
 
model_elevation<-gam(data=alphadt_nomis, formula = value ~ s(ele_resids, bs='tp'))
summary(model_elevation)
model_elevation$coefficients

# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   value ~ s(ele_resids, bs = "tp")
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1846.56      47.21   39.11   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F  p-value    
# s(ele_resids)   1      1 31.35 5.76e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.171   Deviance explained = 17.7%
# GCV = 3.3443e+05  Scale est. = 3.2991e+05  n = 148

##There is no relationship between the %of glacier coverage and the altitude
cov_ele <- ggplot(OR,aes(x=(gl_cov),y=ele_sp))+
  geom_point()+
  geom_smooth(method="gam", formula = y ~ s(x, bs = 'tp'))

cov_ele + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# We test the effect % of glacier coverage on species diversity

pred_cov <- predict(lm(OR, formula = value ~ gl_cov), 
                se.fit = TRUE, interval = "confidence")
limits_cov<- as.data.frame(pred_cov$fit)

spe_cov <- ggplot(OR, aes(x=(gl_cov), y=value)) + geom_point(size=3, alpha=0.4,color="#3F459BFF") + 
  geom_smooth(method='gam', se=T, formula = y ~ s(x, bs='tp'), fill="blue", color="black", alpha=0.2, span=0.3)  +
  geom_line(aes(x = gl_cov, y = limits_cov$lwr), 
            linetype = 2) +
  geom_line(aes(x = gl_cov, y = limits_cov$upr), 
            linetype = 2) 
spe_cov + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

model_cov<-gam(data=OR, formula = value ~ s(gl_cov, bs='tp'))
summary(model_cov)
model_cov$coefficients

# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   value ~ s(gl_cov, bs = "tp")
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   1846.6       49.7   37.15   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value    
# s(gl_cov)   1      1 14.02 0.00026 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0813   Deviance explained = 8.76%
# GCV = 3.7065e+05  Scale est. = 3.6564e+05  n = 148


