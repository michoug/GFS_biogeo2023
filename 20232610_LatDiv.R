library(ggplot2)
library(mgcv)
library(tidyverse)
library(ggpubr)
library(purrr)
library(reshape2)
library(vegan)

setwd('~/Desktop/leila_data')
alphadt_nomis = read.csv('table_diversite_nomis2023.csv', sep=',')
lcbd = read.csv('LCBD.res.csv')
kegg_data = read.csv('kegg_alpha_div_all.csv')

lcbd$X == alphadt_nomis$sample
alphadt_nomis$LCBD = lcbd$LCBD
kegg_data = kegg_data %>% filter(!startsWith(X, 'GLR'), endsWith(clean_name, 'UP'))
kegg_data$clean_name = map_chr(kegg_data$clean_name, function(x) strsplit(x, '_')[[1]][1])
kegg_data = kegg_data %>% select(clean_name, ACE, Chao1, Shannon, Reads) %>% group_by(clean_name) %>% summarise(Chao1=mean(Chao1), Shannon=mean(Shannon), ACE=mean(ACE))

kegg_counts = read.csv('NOMIS_merged_KEGG_counts.csv')
kegg_counts[, colnames(kegg_counts) != 'Geneid'] = sweep(kegg_counts[, colnames(kegg_counts) != 'Geneid'], 
                                                         2, colSums(kegg_counts[, colnames(kegg_counts) != 'Geneid']), '/')
colnames(kegg_counts) = map_chr(colnames(kegg_counts), function(x) gsub('Up', 'UP', x))
kegg_counts = kegg_counts[,(grepl('UP', colnames(kegg_counts)) & !(grepl('GLR', colnames(kegg_counts)))) | (colnames(kegg_counts) == 'Geneid')]
kegg_counts = kegg_counts[rowSums(kegg_counts[,colnames(kegg_counts) != 'Geneid']) > 0,]
kegg_pathways = read.csv('keggPathwayGood.txt', sep='\t')
kegg_counts$KEGG_cat = map_chr(kegg_counts$Geneid, function(x) paste0(kegg_pathways$Path_2[kegg_pathways$ko == x], collapse = ','))
colnames(kegg_counts) = map_chr(colnames(kegg_counts), function(x) strsplit(x, '_')[[1]][1])
kegg_counts_melt = melt(kegg_counts, value.name = 'Count', variable.name = 'sample')
kegg_counts_melt = kegg_counts_melt %>% separate_longer_delim(KEGG, delim = ',')
kegg_counts_melt = kegg_counts_melt %>% filter(Count > 0) %>% group_by(KEGG, sample) %>% summarise(Shannon = diversity(Count))
kegg_counts_melt = kegg_counts_melt %>% filter(KEGG != '')
kegg_counts_melt$KEGG = map_chr(kegg_counts_melt$KEGG, function(x) paste0(strsplit(x,' ')[[1]], collapse = '_'))

kegg_cat_list = c()
for (kegg_cat in unique(kegg_counts_melt$KEGG)){
  kegg_cat_list = c(kegg_cat_list, kegg_cat)
  for (samp in unique(kegg_counts_melt$sample)){
  alphadt_nomis[alphadt_nomis$sample == samp,kegg_cat] = ifelse(kegg_counts_melt %>% filter(sample == samp, KEGG == kegg_cat) %>% nrow() > 0, 
                                                                kegg_counts_melt %>% filter(sample == samp, KEGG == kegg_cat) %>% pull(Shannon),
                                                                0)
}}

alphadt_nomis$Shannon = map_dbl(alphadt_nomis$sample, function(x) ifelse(x %in% kegg_data$clean_name, kegg_data$Shannon[kegg_data$clean_name == x], NA))
alphadt_nomis$Chao1 = map_dbl(alphadt_nomis$sample, function(x) ifelse(x %in% kegg_data$clean_name, kegg_data$Chao1[kegg_data$clean_name == x], NA))
alphadt_nomis$ACE = map_dbl(alphadt_nomis$sample, function(x) ifelse(x %in% kegg_data$clean_name, kegg_data$ACE[kegg_data$clean_name == x], NA))

#########################################################################################################################################################################
## Effect of altitude on species richness
## On check comment la latitude et l'altitude sont lil??es entre elles
## de la latitude
lat_ele <- ggplot(data=alphadt_nomis,aes(x=lat,y=ele_sp))+
  geom_point()+
  geom_smooth(method="gam", formula = y ~ s(x, bs = 'tp', k=3))
lat_ele + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

mod_lat <- gam(data=alphadt_nomis, formula = ele_sp ~ s(lat, bs='tp', k=3))
summary(mod_lat)
alphadt_nomis$ele_resids = mod_lat$residuals

cor_tab = data.frame()
plot_tab = data.frame()

for (var in c('value', 'LCBD', 'Shannon', 'Chao1', 'ACE')){
  ctest = cor.test(alphadt_nomis[,var], alphadt_nomis$ele_resids, method='spearman')
  cor_tab = rbind(cor_tab, data.frame(Variable = var, p = ctest$p.value, statistic = ctest$statistic, rho=ctest$estimate))
  plot_tab = rbind(plot_tab, data.frame(Variable=rep(var, 148), value=alphadt_nomis[,var], norm_elevation=alphadt_nomis$ele_resids))}

for (kegg_cat in kegg_cat_list){
  ctest = cor.test(alphadt_nomis[,kegg_cat], alphadt_nomis$ele_resids)
  cor_tab = rbind(cor_tab, data.frame(Variable = kegg_cat, p = ctest$p.value, statistic = ctest$statistic, rho=ctest$estimate))
  plot_tab = rbind(plot_tab, data.frame(Variable=rep(kegg_cat, 148), value=alphadt_nomis[,kegg_cat], norm_elevation=alphadt_nomis$ele_resids))}

cor_tab$padj = p.adjust(cor_tab$p, method = 'bonferroni')
cor_tab$Variable[cor_tab$Variable == 'value'] = 'ASV_number'
plot_tab$Variable[plot_tab$Variable == 'value'] = 'ASV_number'

write.csv(cor_tab, file='correlations_elevation.csv')

p = ggplot(plot_tab, aes(x=norm_elevation, y=value)) + geom_point() + geom_smooth(method='loess') + facet_wrap(Variable~., scales = 'free', ncol = 3) +
  xlab('Normalised elevation [m]') + ylab('Value')
ggsave('loess_elevation.pdf', p, width = 12, height = 18)


cor_tab %>% filter(padj < 0.05)
#      Variable            p statistic        rho         padj
#                                      ASV_number 9.056167e-10 7.978429e+05 -0.4767376 2.264042e-08
#                                            LCBD 7.530328e-04 3.917580e+05  0.2748901 1.882582e-02
#  09109_Metabolism_of_terpenoids_and_polyketides 2.621585e-06 5.096751e+00  0.5123026 6.553962e-05



##################################################################################################################################################################
cor_tab = data.frame()
plot_tab = data.frame()

for (var in c('value', 'LCBD', 'Shannon', 'Chao1', 'ACE')){
  ctest = cor.test(alphadt_nomis[,var], alphadt_nomis$gl_cov, method='spearman')
  cor_tab = rbind(cor_tab, data.frame(Variable = var, p = ctest$p.value, statistic = ctest$statistic, rho=ctest$estimate))
  plot_tab = rbind(plot_tab, data.frame(Variable=rep(var, 148), value=alphadt_nomis[,var], gl_cov=alphadt_nomis$gl_cov))}

for (kegg_cat in kegg_cat_list){
  ctest = cor.test(alphadt_nomis[,kegg_cat], alphadt_nomis$gl_cov)
  cor_tab = rbind(cor_tab, data.frame(Variable = kegg_cat, p = ctest$p.value, statistic = ctest$statistic, rho=ctest$estimate))
  plot_tab = rbind(plot_tab, data.frame(Variable=rep(kegg_cat, 148), value=alphadt_nomis[,kegg_cat], gl_cov=alphadt_nomis$gl_cov))}

cor_tab$padj = p.adjust(cor_tab$p, method = 'bonferroni')
cor_tab$Variable[cor_tab$Variable == 'value'] = 'ASV_number'
plot_tab$Variable[plot_tab$Variable == 'value'] = 'ASV_number'

write.csv(cor_tab, file='correlations_glcov.csv')

p = ggplot(plot_tab, aes(x=gl_cov, y=value)) + geom_point() + geom_smooth(method='loess') + facet_wrap(Variable~., scales = 'free', ncol = 3) +
  xlab('Normalised elevation [m]') + ylab('Value')
ggsave('loess_glcov.pdf', p, width = 12, height = 18)


cor_tab %>% filter(padj < 0.05)

















summary(gam(data=alphadt_nomis, formula = value ~ gl_cov))

summary(lm(data=alphadt_nomis, formula = LCBD ~ ele_resids:site_c + site_c))

summary(lm(data=alphadt_nomis, formula = Shannon ~ ele_resids:site_c + site_c))






p1 = ggplot(alphadt_nomis, aes(x=ele_resids, y=value, colour=site_c, group=site_c)) + xlab('Elevation model residuals') + ylab('Observed richness') +
  geom_point() + geom_smooth(method='lm', se=F, formula = y ~ x) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + labs(colour="Mountain range")

p2 = ggplot(alphadt_nomis, aes(x=gl_cov, y=value, colour=site_c, group=site_c)) + xlab('Glacier coverage') + ylab('Observed richness') +
  geom_point() + geom_smooth(method='lm', se=F, formula = y ~ x) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + labs(colour="Mountain range")

p = ggarrange(p1, p2, nrow = 2)

ggsave('latitude_glacier_diversity.pdf', p, width = 8, height = 8)


p = ggplot() + xlab('Normalised elevation') + ylab('LCBD') +
  geom_point(data=alphadt_nomis, aes(x=ele_resids, y=LCBD, colour=site_c, group=site_c)) + 
  geom_smooth(data=alphadt_nomis, aes(x=ele_resids, y=LCBD, colour=site_c, group=site_c), method='lm', se=F) +  
  geom_smooth(data=alphadt_nomis, aes(x=ele_resids, y=LCBD), method='lm', colour='black',  linewidth=1.5) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + labs(colour="Mountain range")

ggsave('lcbd_elevation.pdf', p, width = 8, height = 8)





# ## We retrieve the residuals of the models and investigate the correlations between altitude or latitude and species richness
pred <- predict(lm(data=OR, formula = value ~ ele_resids:site_c),
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
