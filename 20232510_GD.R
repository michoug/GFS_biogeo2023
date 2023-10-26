library(phyloseq)
library(biomformat)
library(ggplot2)
library(vegan)
library(knitr)
library(dplyr)
library(data.table)
library(plyr)
library(DescTools)
library(devtools)
library(phyloseqCompanion)
library(iNEXT)


NOMIS_R <- readRDS("20230825_NOMIS_rarefied.RData")
prune_Uganda <- subset_samples(NOMIS_R, !site_c %in% "Uganda")
metadata_nomis_inext <- sample.data.frame(prune_Uganda)

sites <- c("Greenland", "Norway", "Alps", "New_Zealand", "Ecuador", "Caucasus", "Kirghizistan", "Chile", "Alaska", "Nepal")

sample_lists <- list()
asv_lists <- list()
vec_lists <- list()

# Loop through each site and perform the subset and merge operations
for (site in sites) {
   subset_result <- subset_samples(prune_Uganda, site_c == site)
   merged_result <- merge_samples(subset_result,"sample")

   # ## create asv table
   asv_table_site <- otu_table(merged_result, taxa_are_rows=T)
   asv_table_site_t <- t(asv_table_site)
   asv_lists[[site]] <- asv_table_site_t
   #
   # ## calculate rowSums and
   asv_table_df <- asv_table_site_t
   sumrow_site <- unname(rowSums(asv_table_df>0))
   sort_site<- sort(sumrow_site, decreasing=T)
   vec_site <- sort_site[sort_site >0]
   vec_lists[[site]] <- vec_site
}

list_exped_all <- list(alps=c(ncol(asv_lists$Alps),vec_lists$Alps),nz=c(ncol(asv_lists$New_Zealand),vec_lists$New_Zealand),
                       caucasus=c(ncol(asv_lists$Caucasus),vec_lists$Caucasus),
                       kh=c(ncol(asv_lists$Kirghizistan),vec_lists$Kirghizistan),
                       greenland=c(ncol(asv_lists$Greenland),vec_lists$Greenland),
                       norway=c(ncol(asv_lists$Norway),vec_lists$Norway), 
                       ecuador=c(ncol(asv_lists$Ecuador),vec_lists$Ecuador),
                       nepal=c(ncol(asv_lists$Nepal),vec_lists$Nepal),
                       alaska=c(ncol(asv_lists$Alaska),vec_lists$Alaska),
                       chile=c(ncol(asv_lists$Chile),vec_lists$Chile))
                       

out_all_exped <- iNEXT(list_exped_all, q=0, datatype="incidence_freq", se=T, conf=0.95, nboot=99)

df <- fortify(out_all_exped, type =1)

df.point <- df[which(df$Method=="Observed"),]
df.line <- df[which(df$Method!="Observed"),]
df.line$Method <- factor(df.line$Method, 
                         c("Rarefaction", "Extrapolation"),
                       )

df.asympote <- data.frame(y = c(24,10),
                          Asymptote = c("alps","nz","caucasus","kh","greenland","norway","ecuador","nepal","alaska","chile"))


ggplot(df, aes(x=x, y=y, colour=Assemblage)) + 
  #geom_point(aes(shape=Assemblage), size=5, data=df.point) +
  geom_line(aes(linetype= Method), lwd=1.5, data=df.line) +
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
                  fill=Assemblage, colour=NULL), alpha=0.2) +
  labs(x="Number of GFS", y="Species diversity") +
scale_fill_manual(values=c("#2E2A2BFF", "#CF4E9CFF","#8C57A2FF",
                              "#3EBCB6","#82581FFF","#2F509EFF",
                              "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183")
)+
  scale_color_manual(values=c("#2E2A2BFF", "#CF4E9CFF","#8C57A2FF",
                             "#3EBCB6","#82581FFF","#2F509EFF",
                             "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183")
  )+
  scale_linetype_discrete(name ="Method")+
theme_bw() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


# $DataInfo
# Assemblage  T     U S.obs     SC   Q1   Q2   Q3   Q4  Q5  Q6  Q7  Q8  Q9 Q10
# 1        alps 26 49239 10971 0.9368 3273 2118 1229  772 630 495 372 310 269 243
# 2          nz 20 40976  9771 0.9509 2218 2162 1339  975 648 480 404 302 265 198
# 3    caucasus 19 47292 10867 0.9439 2867 2074 1334  996 715 529 411 332 262 228
# 4          kh 16 34069 10204 0.9213 2964 2337 1611 1008 630 438 293 219 168 149
# 5   greenland  6 10071  6337 0.5978 4404  961  461  271 162  78   0   0   0   0
# 6      norway 10 15852  6464 0.8400 2825 1445  876  484 280 176 127  99  84  68
# 7     ecuador  9 13191  5827 0.8239 2650 1493  649  419 210 145 124  89  48   0
# 8       nepal 17 24725  7519 0.8886 2932 1502  859  530 361 277 217 182 143 133
# 9      alaska 15 24617  8279 0.8755 3286 1668 1074  646 451 300 221 160 117  92
# 10      chile 15 13259  6980 0.7600 3457 2095  761  324 162  97  46  28   5   5
# 
# $iNextEst
# $iNextEst$size_based
# Assemblage  t        Method Order.q         qD     qD.LCL     qD.UCL        SC    SC.LCL    SC.UCL
# 1        alps  1   Rarefaction       0  1893.8077  1881.8662  1905.7492 0.3277889 0.3249398 0.3306381
# 2        alps  6   Rarefaction       0  6093.4138  6057.5013  6129.3262 0.7443547 0.7417868 0.7469227
# 3        alps 11   Rarefaction       0  8067.8378  8017.7100  8117.9657 0.8439177 0.8419025 0.8459330
# 4        alps 16   Rarefaction       0  9348.3824  9287.8070  9408.9577 0.8904722 0.8888031 0.8921413
# 5        alps 21   Rarefaction       0 10271.4130 10202.7126 10340.1134 0.9181519 0.9165838 0.9197199
# 6        alps 26      Observed       0 10971.0000 10894.7536 11047.2464 0.9368001 0.9350645 0.9385357
# 7          nz  1   Rarefaction       0  2048.8000  2034.6676  2062.9324 0.3405845 0.3367691 0.3443998
# 8          nz  6   Rarefaction       0  6403.3054  6358.5854  6448.0255 0.7705798 0.7677646 0.7733950
# 9          nz 11   Rarefaction       0  8232.7214  8174.4225  8291.0204 0.8780386 0.8760277 0.8800495
# 10         nz 16   Rarefaction       0  9254.2268  9189.5898  9318.8638 0.9273908 0.9254103 0.9293713
# 11         nz 19   Rarefaction       0  9660.1000  9591.5925  9728.6075 0.9458708 0.9437255 0.9480160
# 12         nz 20      Observed       0  9771.0000  9700.8811  9841.1189 0.9509079 0.9487073 0.9531084
# 13         nz 21 Extrapolation       0  9871.5799  9799.6525  9943.5074 0.9554763 0.9532350 0.9577176
# 14         nz 26 Extrapolation       0 10250.3328 10166.7135 10333.9520 0.9726794 0.9704417 0.9749170
# 15   caucasus  1   Rarefaction       0  2489.0526  2472.2575  2505.8478 0.3897723 0.3868574 0.3926872
# 16   caucasus  6   Rarefaction       0  7246.5907  7203.9500  7289.2315 0.7979407 0.7958204 0.8000611
# 17   caucasus 11   Rarefaction       0  9221.0993  9166.5633  9275.6352 0.8889849 0.8874713 0.8904985
# 18   caucasus 16   Rarefaction       0 10376.5531 10315.0176 10438.0887 0.9290779 0.9277070 0.9304488
# 19   caucasus 18   Rarefaction       0 10716.1053 10652.2758 10779.9347 0.9393766 0.9379827 0.9407706
# 20   caucasus 19      Observed       0 10867.0000 10802.0346 10931.9654 0.9438869 0.9424732 0.9453006
# 21   caucasus 20 Extrapolation       0 11006.6684 10940.5427 11072.7941 0.9480616 0.9466298 0.9494935
# 22   caucasus 21 Extrapolation       0 11135.9458 11068.6160 11203.2756 0.9519258 0.9504798 0.9533717
# 23   caucasus 26 Extrapolation       0 11651.5986 11577.2512 11725.9460 0.9673387 0.9659063 0.9687711
# 24         kh  1   Rarefaction       0  2129.3125  2113.4120  2145.2130 0.3177943 0.3145192 0.3210693
# 25         kh  6   Rarefaction       0  6961.7662  6926.0985  6997.4340 0.7470584 0.7441202 0.7499966
# 26         kh 11   Rarefaction       0  9051.3187  9005.5415  9097.0959 0.8672011 0.8645032 0.8698990
# 27         kh 15   Rarefaction       0 10018.7500  9962.0343 10075.4657 0.9130001 0.9102868 0.9157134
# 28         kh 16      Observed       0 10204.0000 10144.1077 10263.8923 0.9212762 0.9185334 0.9240189
# 29         kh 17 Extrapolation       0 10371.6276 10308.3701 10434.8851 0.9287650 0.9260047 0.9315253
# 30         kh 21 Extrapolation       0 10897.1392 10818.8666 10975.4118 0.9522424 0.9495625 0.9549222
# 31         kh 26 Extrapolation       0 11317.6306 11219.3576 11415.9036 0.9710279 0.9687008 0.9733550
# 32  greenland  1   Rarefaction       0  1678.5000  1654.7321  1702.2679 0.2684937 0.2612170 0.2757703
# 33  greenland  5   Rarefaction       0  5603.0000  5504.6087  5701.3913 0.5627048 0.5526726 0.5727369
# 34  greenland  6      Observed       0  6337.0000  6221.8740  6452.1260 0.5978097 0.5878263 0.6077931
# 35  greenland  7 Extrapolation       0  7012.0764  6881.0450  7143.1079 0.6300965 0.6201131 0.6400798
# 36  greenland 11 Extrapolation       0  9212.2356  9023.0948  9401.3763 0.7353231 0.7256768 0.7449694
# 37  greenland 16 Extrapolation       0 11104.3956 10851.9442 11356.8471 0.8258191 0.8176229 0.8340154
# 38  greenland 21 Extrapolation       0 12349.6048 12045.4269 12653.7827 0.8853735 0.8790430 0.8917041
# 39  greenland 26 Extrapolation       0 13169.0629 12825.5062 13512.6197 0.9245656 0.9199602 0.9291710
# 40     norway  1   Rarefaction       0  1585.2000  1562.2965  1608.1035 0.3356043 0.3301562 0.3410524
# 41     norway  6   Rarefaction       0  5109.8286  5049.7636  5169.8936 0.7457494 0.7401995 0.7512992
# 42     norway  9   Rarefaction       0  6181.5000  6106.3525  6256.6475 0.8217890 0.8162189 0.8273592
# 43     norway 10      Observed       0  6464.0000  6383.6982  6544.3018 0.8399783 0.8342848 0.8456719
# 44     norway 11 Extrapolation       0  6717.6663  6632.0156  6803.3171 0.8563111 0.8505335 0.8620888
# 45     norway 16 Extrapolation       0  7646.6052  7531.5066  7761.7037 0.9161226 0.9107308 0.9215143
# 46     norway 21 Extrapolation       0  8188.8670  8044.5529  8333.1812 0.9510371 0.9467843 0.9552899
# 47     norway 26 Extrapolation       0  8505.4088  8337.2272  8673.5904 0.9714182 0.9683596 0.9744769
# 48    ecuador  1   Rarefaction       0  1465.6667  1444.9331  1486.4003 0.3231938 0.3161163 0.3302713
# 49    ecuador  6   Rarefaction       0  4811.5238  4750.8823  4872.1653 0.7372424 0.7309549 0.7435299
# 50    ecuador  8   Rarefaction       0  5532.5556  5461.5001  5603.6110 0.7991055 0.7924271 0.8057838
# 51    ecuador  9      Observed       0  5827.0000  5751.0612  5902.9388 0.8239079 0.8169330 0.8308827
# 52    ecuador 10 Extrapolation       0  6085.0924  6004.1144  6166.0704 0.8456482 0.8384680 0.8528284
# 53    ecuador 11 Extrapolation       0  6311.3207  6224.9731  6397.6683 0.8647044 0.8574343 0.8719746
# 54    ecuador 16 Extrapolation       0  7086.3986  6968.6474  7204.1499 0.9299928 0.9236005 0.9363851
# 55    ecuador 21 Extrapolation       0  7487.4541  7338.8741  7636.0341 0.9637756 0.9591061 0.9684450
# 56    ecuador 26 Extrapolation       0  7694.9758  7522.8570  7867.0946 0.9812561 0.9781503 0.9843618
# 57      nepal  1   Rarefaction       0  1454.4118  1440.6937  1468.1298 0.3388018 0.3344954 0.3431082
# 58      nepal  6   Rarefaction       0  4689.8589  4651.8515  4727.8663 0.7332572 0.7294899 0.7370245
# 59      nepal 11   Rarefaction       0  6289.5372  6233.9096  6345.1649 0.8330153 0.8297397 0.8362909
# 60      nepal 16   Rarefaction       0  7346.5294  7274.8439  7418.2150 0.8814156 0.8779010 0.8849302
# 61      nepal 17      Observed       0  7519.0000  7443.9638  7594.0362 0.8885521 0.8849019 0.8922023
# 62      nepal 18 Extrapolation       0  7681.0911  7602.5777  7759.6046 0.8952592 0.8914853 0.8990330
# 63      nepal 21 Extrapolation       0  8111.1484  8021.3158  8200.9811 0.9130542 0.9090231 0.9170854
# 64      nepal 26 Extrapolation       0  8671.7650  8560.4548  8783.0753 0.9362516 0.9321764 0.9403268
# 65     alaska  1   Rarefaction       0  1641.1333  1625.9735  1656.2932 0.3185255 0.3139031 0.3231479
# 66     alaska  6   Rarefaction       0  5452.1650  5403.8090  5500.5211 0.7234834 0.7198410 0.7271258
# 67     alaska 11   Rarefaction       0  7297.5040  7228.7249  7366.2832 0.8328726 0.8297903 0.8359548
# 68     alaska 14   Rarefaction       0  8059.9333  7982.1488  8137.7178 0.8665150 0.8633198 0.8697102
# 69     alaska 15      Observed       0  8279.0000  8198.6134  8359.3866 0.8755403 0.8722653 0.8788152
# 70     alaska 16 Extrapolation       0  8483.2550  8400.3719  8566.1382 0.8839553 0.8806041 0.8873065
# 71     alaska 21 Extrapolation       0  9315.1305  9219.5600  9410.7009 0.9182274 0.9147247 0.9217302
# 72     alaska 26 Extrapolation       0  9901.3234  9791.1624 10011.4843 0.9423778 0.9391026 0.9456530
# 73      chile  1   Rarefaction       0   883.9333   870.6759   897.1907 0.1244653 0.1205765 0.1283542
# 74      chile  6   Rarefaction       0  4008.3916  3960.7097  4056.0735 0.4849088 0.4772996 0.4925179
# 75      chile 11   Rarefaction       0  5931.4916  5867.5404  5995.4428 0.6656097 0.6579175 0.6733019
# 76      chile 14   Rarefaction       0  6749.5333  6676.1696  6822.8971 0.7392714 0.7317294 0.7468135
# 77      chile 15      Observed       0  6980.0000  6903.2513  7056.7487 0.7600452 0.7524891 0.7676014
# 78      chile 16 Extrapolation       0  7192.1040  7111.8190  7272.3890 0.7791639 0.7715900 0.7867378
# 79      chile 21 Extrapolation       0  8024.5014  7924.3454  8124.6573 0.8541945 0.8468088 0.8615802
# 80      chile 26 Extrapolation       0  8574.0859  8452.0569  8696.1148 0.9037329 0.8970722 0.9103936
# 
# $iNextEst$coverage_based
# Assemblage        SC  t        Method Order.q         qD     qD.LCL     qD.UCL
# 1        alps 0.3277889  1   Rarefaction       0  1893.8077  1881.8662  1905.7492
# 2        alps 0.7443541  6   Rarefaction       0  6093.4025  6039.6117  6147.1933
# 3        alps 0.8439180 11   Rarefaction       0  8067.8431  7992.9190  8142.7671
# 4        alps 0.8904721 16   Rarefaction       0  9348.3811  9260.7019  9436.0603
# 5        alps 0.9181516 21   Rarefaction       0 10271.4024 10176.2001 10366.6048
# 6        alps 0.9368001 26      Observed       0 10971.0000 10864.0622 11077.9378
# 7          nz 0.3405845  1   Rarefaction       0  2048.8000  2034.6676  2062.9324
# 8          nz 0.7705796  6   Rarefaction       0  6403.3019  6339.9697  6466.6340
# 9          nz 0.8780385 11   Rarefaction       0  8232.7189  8158.5555  8306.8822
# 10         nz 0.9273903 16   Rarefaction       0  9254.2168  9172.9049  9335.5287
# 11         nz 0.9458708 19   Rarefaction       0  9622.3755  9488.8946  9755.8565
# 12         nz 0.9509079 20      Observed       0  9771.0000  9674.5338  9867.4662
# 13         nz 0.9554763 21 Extrapolation       0  9871.5799  9663.9978 10079.1621
# 14         nz 0.9726794 26 Extrapolation       0 10250.3328 10124.3680 10376.2976
# 15   caucasus 0.3897723  1   Rarefaction       0  2489.0526  2472.2575  2505.8478
# 16   caucasus 0.7979404  6   Rarefaction       0  7246.5849  7191.6982  7301.4715
# 17   caucasus 0.8889848 11   Rarefaction       0  9221.0969  9151.6606  9290.5332
# 18   caucasus 0.9290779 16   Rarefaction       0 10376.5505 10300.0589 10453.0422
# 19   caucasus 0.9393766 17   Rarefaction       0 10592.6858 10518.7900 10666.5816
# 20   caucasus 0.9438869 19      Observed       0 10867.0000 10776.2743 10957.7257
# 21   caucasus 0.9480616 20 Extrapolation       0 11006.6684 10903.3232 11110.0137
# 22   caucasus 0.9519258 21 Extrapolation       0 11135.9458 10888.3017 11383.5899
# 23   caucasus 0.9673387 26 Extrapolation       0 11651.5986 11552.9696 11750.2275
# 24         kh 0.3177943  1   Rarefaction       0  2129.3125  2113.4120  2145.2130
# 25         kh 0.7470579  6   Rarefaction       0  6961.7593  6916.8918  7006.6267
# 26         kh 0.8672009 11   Rarefaction       0  9051.3147  8978.3608  9124.2686
# 27         kh 0.9130001 14   Rarefaction       0  9892.4432  9775.4958 10009.3906
# 28         kh 0.9212762 16      Observed       0 10204.0000 10060.3202 10347.6798
# 29         kh 0.9287650 17 Extrapolation       0 10371.6276 10232.9162 10510.3391
# 30         kh 0.9522424 21 Extrapolation       0 10897.1392 10771.8470 11022.4314
# 31         kh 0.9710279 26 Extrapolation       0 11317.6306 11174.7506 11460.5107
# 32  greenland 0.2684937  1   Rarefaction       0  1678.5000  1654.7321  1702.2679
# 33  greenland 0.5627048  4   Rarefaction       0  5083.2927  4847.8148  5318.7706
# 34  greenland 0.5978097  6      Observed       0  6337.0000  6025.4366  6648.5634
# 35  greenland 0.6300965  7 Extrapolation       0  7012.0764  6753.2565  7270.8963
# 36  greenland 0.7353231 11 Extrapolation       0  9212.2356  8934.7662  9489.7049
# 37  greenland 0.8258191 16 Extrapolation       0 11104.3956 10779.4232 11429.3681
# 38  greenland 0.8853735 21 Extrapolation       0 12349.6048 11990.4039 12708.8057
# 39  greenland 0.9245656 26 Extrapolation       0 13169.0629 12786.4114 13551.7144
# 40     norway 0.3356043  1   Rarefaction       0  1585.2000  1562.2965  1608.1035
# 41     norway 0.7457484  6   Rarefaction       0  5109.8159  5026.6197  5193.0120
# 42     norway 0.8217890  8   Rarefaction       0  5933.6592  5812.8068  6054.5117
# 43     norway 0.8399783 10      Observed       0  6464.0000  6327.2333  6600.7667
# 44     norway 0.8563111 11 Extrapolation       0  6717.6663  6612.1051  6823.2276
# 45     norway 0.9161226 16 Extrapolation       0  7646.6052  7486.5681  7806.6422
# 46     norway 0.9510371 21 Extrapolation       0  8188.8670  8005.6259  8372.1082
# 47     norway 0.9714182 26 Extrapolation       0  8505.4088  8307.8838  8702.9339
# 48    ecuador 0.3231938  1   Rarefaction       0  1465.6667  1444.9331  1486.4003
# 49    ecuador 0.7372401  6   Rarefaction       0  4811.4952  4722.6590  4900.3313
# 50    ecuador 0.7991055  8   Rarefaction       0  5433.3196  5317.2469  5549.3922
# 51    ecuador 0.8239079  9      Observed       0  5827.0000  5639.2927  6014.7073
# 52    ecuador 0.8456482 10 Extrapolation       0  6085.0924  6006.0322  6164.1525
# 53    ecuador 0.8647044 11 Extrapolation       0  6311.3207  6128.2523  6494.3891
# 54    ecuador 0.9299928 16 Extrapolation       0  7086.3986  6921.0276  7251.7696
# 55    ecuador 0.9637756 21 Extrapolation       0  7487.4541  7300.0377  7674.8705
# 56    ecuador 0.9812561 26 Extrapolation       0  7694.9758  7495.6294  7894.3222
# 57      nepal 0.3388018  1   Rarefaction       0  1454.4118  1440.6937  1468.1298
# 58      nepal 0.7332568  6   Rarefaction       0  4689.8532  4631.1765  4748.5299
# 59      nepal 0.8330153 11   Rarefaction       0  6289.5371  6202.7224  6376.3518
# 60      nepal 0.8814156 15   Rarefaction       0  7213.5153  7100.1230  7326.9077
# 61      nepal 0.8885521 17      Observed       0  7519.0000  7401.4386  7636.5614
# 62      nepal 0.8952592 18 Extrapolation       0  7681.0911  7566.2749  7795.9074
# 63      nepal 0.9130542 21 Extrapolation       0  8111.1484  7968.3624  8253.9344
# 64      nepal 0.9362516 26 Extrapolation       0  8671.7650  8501.9451  8841.5849
# 65     alaska 0.3185255  1   Rarefaction       0  1641.1333  1625.9735  1656.2932
# 66     alaska 0.7234833  6   Rarefaction       0  5452.1641  5382.1754  5522.1528
# 67     alaska 0.8328724 11   Rarefaction       0  7297.4992  7200.0968  7394.9016
# 68     alaska 0.8665150 13   Rarefaction       0  7902.2041  7795.5164  8008.8918
# 69     alaska 0.8755403 15      Observed       0  8279.0000  8160.0451  8397.9549
# 70     alaska 0.8839553 16 Extrapolation       0  8483.2550  8308.9436  8657.5665
# 71     alaska 0.9182274 21 Extrapolation       0  9315.1305  9182.1494  9448.1116
# 72     alaska 0.9423778 26 Extrapolation       0  9901.3234  9752.1084 10050.5383
# 73      chile 0.1244657  1   Rarefaction       0   883.9357   870.6783   897.1931
# 74      chile 0.4849081  6   Rarefaction       0  4008.3843  3949.9017  4066.8669
# 75      chile 0.6656082 11   Rarefaction       0  5931.4745  5841.2210  6021.7279
# 76      chile 0.7392714 13   Rarefaction       0  6581.3014  6473.4706  6689.1322
# 77      chile 0.7600452 15      Observed       0  6980.0000  6869.1185  7090.8815
# 78      chile 0.7791639 16 Extrapolation       0  7192.1040  7065.7080  7318.5000
# 79      chile 0.8541945 21 Extrapolation       0  8024.5014  7873.1847  8175.8180
# 80      chile 0.9037329 26 Extrapolation       0  8574.0859  8399.2233  8748.9484
# 
# 
# $AsyEst
# Assemblage         Diversity  Observed Estimator      s.e.       LCL       UCL
# 1      alaska  Species richness  8279.000 11299.966 112.93000 11078.627 11521.305
# 2      alaska Shannon diversity  5971.696  7486.666  39.67549  7408.903  7564.428
# 3      alaska Simpson diversity  4509.139  5152.283  29.46771  5094.527  5210.038
# 4        alps  Species richness 10971.000 13402.659 102.13119 13202.486 13602.833
# 5        alps Shannon diversity  7224.810  8253.602  33.08169  8188.763  8318.441
# 6        alps Simpson diversity  5355.136  5777.521  24.75980  5728.993  5826.049
# 7    caucasus  Species richness 10867.000 12744.308  93.31084 12561.422 12927.194
# 8    caucasus Shannon diversity  7627.022  8633.000  30.89630  8572.445  8693.556
# 9    caucasus Simpson diversity  5899.773  6385.914  24.96278  6336.988  6434.840
# 10      chile  Species richness  6980.000  9642.083 130.24675  9386.804  9897.361
# 11      chile Shannon diversity  5844.307  8629.151  65.31320  8501.139  8757.162
# 12      chile Simpson diversity  4834.614  7101.843  64.47879  6975.467  7228.219
# 13    ecuador  Species richness  5827.000  7917.496 111.00976  7699.921  8135.072
# 14    ecuador Shannon diversity  4598.855  6149.780  50.68496  6050.440  6249.121
# 15    ecuador Simpson diversity  3678.933  4534.947  36.06043  4464.269  4605.624
# 16  greenland  Species richness  6337.000 14746.303 238.65467 14278.548 15214.057
# 17  greenland Shannon diversity  5283.485 10366.624 130.30364 10111.234 10622.014
# 18  greenland Simpson diversity  4299.311  6251.544  73.95375  6106.597  6396.490
# 19         kh  Species richness 10204.000 11966.134  93.97231 11781.952 12150.316
# 20         kh Shannon diversity  7608.094  8987.334  33.75327  8921.179  9053.489
# 21         kh Simpson diversity  5907.665  6700.286  28.27566  6644.866  6755.705
# 22      nepal  Species richness  7519.000 10212.389 118.23874  9980.645 10444.133
# 23      nepal Shannon diversity  5185.046  6355.740  37.36464  6282.507  6428.974
# 24      nepal Simpson diversity  3850.749  4292.810  24.74401  4244.313  4341.307
# 25     norway  Species richness  6464.000  8949.316  98.84172  8755.590  9143.042
# 26     norway Shannon diversity  4989.834  6552.162  41.01421  6471.776  6632.549
# 27     norway Simpson diversity  3942.853  4723.420  31.46197  4661.755  4785.084
# 28         nz  Species richness  9771.000 10851.839  73.64223 10707.503 10996.175
# 29         nz Shannon diversity  7067.157  7983.304  28.97589  7926.512  8040.096
# 30         nz Simpson diversity  5484.597  6015.541  27.74868  5961.154  6069.927


inext_freq_results <- out_all_exped$AsyEst  
inext_freq_results$prop <- inext_freq_results$Observed/inext_freq_results$Estimator
inext_freq_results<-inext_freq_results[inext_freq_results$Diversity == 'Species richness',]
median_GD_freq <-inext_freq_results %>% 
  summarise(median=median(prop), x = quantile(prop, c(0.25, 0.5, 0.75)))
# > median_GD_freq
# median         x
# 1 0.7361138 0.7260968
# 2 0.7361138 0.7361138
# 3 0.7361138 0.8441630


## Here we investigate gamma diversity but through abundance 
t <- seq(1, 30, by=5)
out_all_exped_abundance <- iNEXT(list_exped_all, q=0, datatype="abundance", size=t)

# $AsyEst: asymptotic diversity estimates along with related statistics.
# Assemblage         Diversity  Observed Estimator      s.e.       LCL       UCL
# 1      alaska  Species richness  8280.000 11516.618 127.78777 11266.159 11767.077
# 2      alaska Shannon diversity  5970.637  7762.031  51.19454  7661.692  7862.371
# 3      alaska Simpson diversity  4507.090  5516.255  52.42421  5413.505  5619.004
# 4        alps  Species richness 10972.000 13500.874 108.71263 13287.801 13713.947
# 5        alps Shannon diversity  7223.519  8418.669  38.98782  8342.254  8495.084
# 6        alps Simpson diversity  5352.801  6005.174  38.95408  5928.825  6081.522
# 7    caucasus  Species richness 10868.000 12849.561 100.37407 12652.831 13046.291
# 8    caucasus Shannon diversity  7626.656  8871.819  37.27515  8798.761  8944.877
# 9    caucasus Simpson diversity  5898.897  6739.016  38.11113  6664.320  6813.713
# 10      chile  Species richness  6981.000  9833.016 118.06626  9601.611 10064.422
# 11      chile Shannon diversity  5838.443  8956.558  65.29938  8828.574  9084.543
# 12      chile Simpson diversity  4815.761  7557.082  78.41582  7403.389  7710.774
# 13    ecuador  Species richness  5828.000  8179.630 118.25742  7947.850  8411.411
# 14    ecuador Shannon diversity  4598.406  6550.202  61.95108  6428.781  6671.624
# 15    ecuador Simpson diversity  3677.656  5097.630  63.37547  4973.416  5221.844
# 16  greenland  Species richness  6338.000 16428.162 318.86093 15803.206 17053.118
# 17  greenland Shannon diversity  5283.025 11644.769 173.00971 11305.677 11983.862
# 18  greenland Simpson diversity  4297.877  7493.421 140.37926  7218.283  7768.559
# 19         kh  Species richness 10205.000 12084.555 108.60319 11871.696 12297.413
# 20         kh Shannon diversity  7607.118  9286.595  50.59403  9187.433  9385.757
# 21         kh Simpson diversity  5905.521  7142.920  52.85124  7039.334  7246.507
# 22      nepal  Species richness  7520.000 10381.610 115.95424 10154.344 10608.876
# 23      nepal Shannon diversity  5184.081  6559.086  48.57340  6463.884  6654.288
# 24      nepal Simpson diversity  3849.039  4557.950  44.03363  4471.646  4644.254
# 25     norway  Species richness  6465.000  9226.288 130.34508  8970.816  9481.760
# 26     norway Shannon diversity  4989.373  6930.627  56.12707  6820.619  7040.634
# 27     norway Simpson diversity  3941.644  5244.678  55.67814  5135.550  5353.805
# 28         nz  Species richness  9772.000 10909.698  78.50996 10755.821 11063.574
# 29         nz Shannon diversity  7066.336  8190.699  34.57596  8122.932  8258.467
# 30         nz Simpson diversity  5482.788  6329.107  38.43849  6253.769  6404.445


inext_abun_results <- out_all_exped_abundance$AsyEst  
inext_abun_results$prop <- inext_abun_results$Observed/inext_abun_results$Estimator
inext_abun_results<-inext_abun_results[inext_abun_results$Diversity == 'Species richness',]
median_GD_abun <-inext_abun_results %>% 
  summarise(median=median(prop), x = quantile(prop, c(0.25, 0.5, 0.75)))
# > median_GD_abun
# median         x
# 1 0.7216594 0.7105917
# 2 0.7216594 0.7216594
# 3 0.7216594 0.8365218

colnames(inext_freq_results) <- c("Assemblage","Diversity","a_observed","b_estimated","se","LCL","UCL","prop")
res_inext <- inext_freq_results |> pivot_longer(cols = -c(Assemblage, Diversity, LCL, UCL), names_to = "Type") 
head(res_inext)
merged_pivot <- merge(inext_freq_results, res_inext, by.x="Assemblage",by.y="Assemblage")

##Assign SD of Percent to 0
merged_pivot$se[merged_pivot$Type == 'a_observed'] = 0

merged_pivot <- as.data.frame(merged_pivot)
merged_pivot$s.e. <- as.numeric(merged_pivot$se)
merged_pivot$value <- as.numeric(merged_pivot$value)

merged_pivot<-merged_pivot[merged_pivot$Type != 'prop',]
merged_pivot<-merged_pivot[merged_pivot$Type != 'se',]

##Plot everything!
regional_div <- ggplot(merged_pivot, aes(x=Assemblage, y = value, fill= Type)) + 
  geom_col(position="dodge") + 
  labs(y="Average")+ geom_errorbar(aes(ymin=(value - s.e.), ymax=(value + s.e.)),
                                   width=.2, colour="red", 
                                   position=position_dodge(.9))

regional_div + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



