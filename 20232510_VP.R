##Variance partition 

#### Setting R ####
# load packages
library(here)
library(ggplot2)
library(tidyverse)
library(geosphere)   # distm()
library(vegan)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(svd)
library(propack.svd)
# load functions (source Legendre 2016) - not currently used (need to install package AEM)
'PCNM' <-
  function(matdist, thresh=NULL, dbMEM=FALSE, moran=NULL, all=FALSE, include.zero=FALSE, silent=FALSE)
    #
    # Compute the PCNM or dbMEM eigenfunctions corresponding to
    # all eigenvalues (+, 0, -).
    #    In PCNM computation, the diagonal of D = 0.
    #    In dbMEM, the diagonal of D = 4*threshh.
    #    Distance-based MEM are described in Dray et al. 2006.
    #    The name was abbreviated to db-MEM by PPN & PL (subm.)
    # Input file: distance matrix produced by the function "dist".
    # Computation of the threshold requires a function of the library "ape".
    #
    # Original PCNM function: Stephane Dray, November 11, 2004
# The present version: Pierre Legendre, August 2007, January and March 2009
  {
    require(vegan)
    epsilon <- sqrt(.Machine$double.eps)
    a <- system.time({
      if(is.null(moran)) {
        if(dbMEM) { moran=FALSE } else { moran=TRUE }
      }
      single <- FALSE
      if(moran) {
        # cat("The site coordinates were computed from 'matdist'.",'\n')
        pcoa.xy <- pcoa.all(matdist)
        
        if(is.na(pcoa.xy$values[2]) | (pcoa.xy$values[2] < epsilon)) {
          if(!silent) cat("The sites form a straight line on the map.",'\n')
          xy <- pcoa.xy$vectors
          single <- TRUE
        } else {
          xy <- pcoa.xy$vectors[,1:2]
        }
      }
      
      matdist <- as.matrix(matdist)
      n <- nrow(matdist)
      
      # Truncation of distance matrix
      if(is.null(thresh)) {
        spanning <- vegan::spantree(as.dist(matdist))
        threshh <- max(spanning$dist)
        if(!silent) cat("Truncation level =",threshh+0.000001,'\n')
      } else {
        threshh = thresh
        if(!silent) cat("User-provided truncation threshold =",thresh,'\n')
      }
      matdist[matdist > threshh] <- 4*threshh
      
      if(dbMEM==FALSE) { diagonal <- 0 } else { diagonal <- 4*threshh }
      
      mypcnm.all <- pcoa.all(matdist, diagonal=diagonal, all=all, include.zero=include.zero, rn=rownames(matdist))
      
      # Compute Moran's I
      if(moran) {
        require(AEM)
        if(single) {
          nb <- dnearneigh(matrix(c(xy,rep(0,n)),n,2), 0, (threshh + epsilon))
        } else {
          nb <- dnearneigh(xy, 0, (threshh + epsilon))
        }
        fr.to.pcnm2 <- as.matrix(listw2sn(nb2listw(nb))[,1:2])
        weight.dist.coord.mat <- as.matrix(1-(as.dist(matdist)/(4*threshh))^2)
        weight <- weight.dist.coord.mat[fr.to.pcnm2]
        res <- moran.I.multi(mypcnm.all$vectors, link=fr.to.pcnm2, weight=weight)
        Moran <- res$res.mat[,1:2]
        positive <- rep(FALSE,length(mypcnm.all$values))
        positive[which(Moran[,1] > res$expected)] <- TRUE
        Moran <- cbind(as.data.frame(Moran), positive)
        colnames(Moran) <- c("Moran","p.value","Positive")
      }
    })
    a[3] <- sprintf("%2f",a[3])
    if(!silent) cat("Time to compute PCNMs =",a[3]," sec",'\n')
    if(is.null(thresh)) {
      if(moran) {
        res <- list(values=mypcnm.all$values, vectors=mypcnm.all$vectors, Moran_I=Moran, expected_Moran=res$expected, spanning=spanning, thresh=threshh+0.000001)
      } else {
        res <- list(values=mypcnm.all$values, vectors=mypcnm.all$vectors, spanning=spanning, thresh=threshh+0.000001)
      }
    } else {
      if(moran) {
        res <- list(values=mypcnm.all$values, vectors=mypcnm.all$vectors, Moran_I=Moran, expected_Moran=res$expected, thresh=thresh)
      } else {
        res <- list(values=mypcnm.all$values, vectors=mypcnm.all$vectors, thresh=threshh+0.000001)
      }
    }
    res
  }
'pcoa.all' <- function(D, diagonal=0, all=FALSE, include.zero=FALSE, rn=NULL)
  # Principal coordinate decomposition of a square distance matrix D
  # Get the eigenvectors corresponding to all eigenvalues, positive and negative
  # Pierre Legendre, 2005, 2007
  #
  # D : A distance matrix of class 'dist' or 'matrix'.
  # all : If TRUE, the eigenvectors corresponding to all eigenvalues, positive and negative, are shown in the output list.
  # include.zero : If FALSE (default value), the zero eigenvalues as well as their eigenvectors are excluded from the output list.
  # rn : An optional vector of row names, of length n, for the objects.
{
  epsilon <- sqrt(.Machine$double.eps)
  # replace by:     epsilon <- .Machine$double.eps * 10^2
  D <- as.matrix(D)
  n <- nrow(D)
  D <- D + diag(rep(diagonal,n))
  
  # Gower centring, matrix formula
  One <- matrix(1,n,n)
  mat <- diag(n) - One/n
  Dpr2 <- -0.5 * mat %*% (D^2) %*% mat
  trace <- sum(diag(Dpr2))
  
  # Eigenvalue decomposition
  D.eig <- eigen(Dpr2, symmetric=TRUE)
  rel.values <- D.eig$values/trace
  rel.cum <- cumsum(rel.values)
  if(length(rn)!=0) {
    rownames(D.eig$vectors) <- rn
  } else {
    rownames(D.eig$vectors) <- rownames(D)
  }
  
  # Output the results: k eigenvalues and eigenvectors
  if(all) {
    select <- 1:n
    if(!include.zero) {
      exclude <- which(abs(D.eig$values) < epsilon)
      select <- select[-exclude]
    }
    k <- length(select)
    res <- list(values=D.eig$values[select], rel.values=rel.values[select], rel.cum.values=rel.cum[select], vectors=D.eig$vectors[,select], trace=trace)
    # cat("k =",k,"Select =",select,'\n')
    
  } else {
    
    k <- length(which(D.eig$values > epsilon))        
    weight <- sqrt(D.eig$values[1:k])
    if(k == 1) {
      vectors <- D.eig$vectors[,1]*sqrt(D.eig$values[1])
    } else {
      vectors <- D.eig$vectors[,1:k]%*%diag(weight)
    }
    res <- list(values=D.eig$values[1:k], rel.values=rel.values[1:k], rel.cum.values=rel.cum[1:k], vectors=vectors, trace=trace)
  }
  res
}

# load R data
load("comY_20232309.Rdata") # asv_community_trim
load("Env_20232309.Rdata")   # multi_data_na_filter
load("Geo_20232309.Rdata")   # dist_geo
load("geocoordinates_20232309.Rdata")   # metadata_dist_df

##log transformation
bio.b <- log1p(Y.com)

# rename environmental data for ease of use
## we checked several times and we had to remove Clays + Feldspar for multicollinearity purpose!
env.data <- multi_data_trans_ming
env.data <- multi_data_trans_ming %>%
  select(-gl_sa.x,-Feldspar) ## We remove them for multicollinearity !

env.names <- c("water_temp", "pH","cond", "turb","doc", "srp", "mean_chla", "gl_cov.x","Calcite","Quartz","pr","scd","DIN")
env.mat <- multi_data_trans_ming[,env.names]
env.mat <- as.data.frame(sapply(env.mat, as.numeric))

# rename geographic coordinates for ease of us
xyz.dat <- metadata_dist_df

# check data structure ## nb of glaciers
with(env.data, table(site_c))
# 
# site_c
# Alaska         Alps     Caucasus        Chile      Ecuador    Greenland Kirghizistan        Nepal 
# 15           21           15            9            9            6           17           16 
# New_Zealand       Norway 
# 18            9 

#### 1) Linear trend in community structure #### On regarde s'il y a un effet de la position géographique sur la composition (ça pourrait expliquer 5%)
rda.linear <- dbrda(bio.b~., data=xyz.dat, distance="bray")
anova(rda.linear, step=1000)

# Permutation test for dbrda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: dbrda(formula = bio.b ~ lon.x + lat.x + ele, data = xyz.dat, distance = "bray")
# Df SumOfSqs      F Pr(>F)    
# Model      3    8.048 8.4505  0.001 ***
#   Residual 131   41.587                  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 

RsquareAdj(rda.linear)$adj.r.squared 
# > RsquareAdj(rda.linear)$adj.r.squared 
# [1] 0.1429563


#### Nested spatial model ####  modified December 7th 2022
## Analyse of microbial spatial variation among and within regions by means of a two-level spatial model
## The among-region component is modeled by a set of dummy variables (N-1 variables with N the number of regions)
## The within-region component is modeled by a set of db-MEM variables for each region
## The db-MEM variables were arranged in blocks corresponding to each region
## within each block, all sites belonging to other regions received the value 0
## See Borcard et al. 2011 (Num ecol with R), Declerck et al. (2011) and function create.MEM.model()

# creating data.frame to store the dbMEMs
var.data.hier <- env.data
var.data.hier <- as.tibble(var.data.hier)
library(geosphere)
##create matrix of geographic coordinates
xy.dat = xyz.dat[,1:2]

## ****************************** ##
## SUBSET = Alps                  ##
## ****************************** ##
i = "Alps"
var.data.temp <- subset(env.data, site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
dbMEM.region <- cbind(data.frame(sample=var.data.temp$sample), dbMEM.region)
var.data.hier <- merge(var.data.hier, dbMEM.region, by="sample", sort=F, all.x=T, nomatch=0)


## ****************************** ##
## SUBSET = Greenland             ##
## ****************************** ##
i = "Greenland"
var.data.temp <- subset(env.data, site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(sample=var.data.temp$sample), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="sample", sort=F, all.x=T, nomatch=0)


# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region),distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared  
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call
# 
# R2a <- RsquareAdj(rda.region)$adj.r.squared   # R2a = 0.169
# eigenvals(rda.region)[1:2]/sum(eigenvals(rda.region))*RsquareAdj(rda.region)$adj.r.squared/RsquareAdj(rda.region)$r.squared
# axes.region <- scores(rda.region, choices=c(1:2), display="lc", scaling=1)
# pdf(here("Figures", "Figure_Greeland.pdf"), width=12, height=12)
# PCNM.region(axes.region, 1, "Axe1 (16,9%%)")
# dev.off()

## ****************************** ##
## SUBSET = Caucasus              ##
## ****************************** ##
i = "Caucasus"
var.data.temp <- subset(env.data, site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(sample=var.data.temp$sample), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="sample", sort=F, all.x=T, nomatch=0)


# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region),distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared  
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call
# 
# eigenvals(rda.region)[1:2]/sum(eigenvals(rda.region))*RsquareAdj(rda.region)$adj.r.squared/RsquareAdj(rda.region)$r.squared
# axes.region <- scores(rda.region, choices=c(1:2), display="lc", scaling=1)
# pdf(here("Figures", "Figure_Caucasus"), width=12, height=12)
# PCNM.region(axes.region, 1, "Axe1 (5.1%)")
# dev.off()


## ****************************** ##
## SUBSET = Chile                 ##
## ****************************** ##
i = "Chile"
var.data.temp <- subset(env.data, site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(sample=var.data.temp$sample), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="sample", sort=F, all.x=T, nomatch=0)


# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared   # R2a = [1] -0.1062698
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call
# 
# 
# eigenvals(rda.region)[1:2]/sum(eigenvals(rda.region))*RsquareAdj(rda.region)$adj.r.squared/RsquareAdj(rda.region)$r.squared
# axes.region <- scores(rda.region, choices=c(1:2), display="lc", scaling=1)
# pdf(here("Figures", "Figure_Caucasus"), width=12, height=12)
# PCNM.region(axes.region, 1, "Axe1 (5.1%)")
# dev.off()
# 

## ****************************** ##
## SUBSET = Norway                ##
## ****************************** ##
i = "Norway"
var.data.temp <- subset(env.data, site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
(sel.pos <- which(dbMEM.temp$values > 0))
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(sample=var.data.temp$sample), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="sample", sort=F, all.x=T, nomatch=0)

# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared   # R2a = [1] -0.1062698
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call
# 


## ****************************** ##
## SUBSET = Nepal                 ##
## ****************************** ##
i = "Nepal"
var.data.temp <- subset(env.data, site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(sample=var.data.temp$sample), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="sample", sort=F, all.x=T, nomatch=0)

# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared  
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call


## ****************************** ##
## SUBSET = Kirghizistan          ##
## ****************************** ##
i = "Kirghizistan"
var.data.temp <- subset(env.data, site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(sample=var.data.temp$sample), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="sample", sort=F, all.x=T, nomatch=0)

# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared 
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call
# 


## ****************************** ##
## SUBSET = Ecuador               ##
## ****************************** ##
i = "Ecuador"
var.data.temp <- subset(env.data, site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(sample=var.data.temp$sample), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="sample", sort=F, all.x=T, nomatch=0)


# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared 
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call


## ****************************** ##
## SUBSET = New_Zealand           ##
## ****************************** ##
i = "New_Zealand"
var.data.temp <- subset(env.data, site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
(sel.pos <- which(dbMEM.temp$values > 0))
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(sample=var.data.temp$sample), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="sample", sort=F, all.x=T, nomatch=0)

# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared 
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call


## ****************************** ##
## SUBSET = Alaska                ##
## ****************************** ##
i = "Alaska"
var.data.temp <- subset(env.data, site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
(sel.pos <- which(dbMEM.temp$values > 0))
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(sample=var.data.temp$sample), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="sample", sort=F, all.x=T, nomatch=0)

# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared 
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call
# 


## ****************************** ##
## Editing variables              ##
## ****************************** ##
# Replace all NA in db-MEM by 0
var.data.hier[is.na(var.data.hier)] <- 0
# reorder to match community matrix
var.data.hier$sample == env.data$sample
var.data.hier2 <- merge(env.data, var.data.hier[,!(colnames(var.data.hier)%in%env.names)], by=c("sample","site_c"), sort=F, all.x=T, nomatch=0)
var.data.hier2$sample == env.data$sample
# create a region dummy variable (N-1, with N the number of regions)
region.dum <- as.matrix(model.matrix(~-1+var.data.hier2$site_c))
region.dum <- region.dum[,1:(dim(region.dum)[2]-1)]
# db-MEM only dataset  
region.dbMEM <- var.data.hier2[,16:dim(var.data.hier2)[2]]



#### 3) Models ####
# 3a) variation among regions at the global scale
# Spatial model
rda1.S <-dbrda(bio.b ~., data=as.data.frame(region.dum), distance="bray")
# anova(rda1.S, step=1000)
# Permutation test for dbrda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: dbrda(formula = bio.b ~ `var.data.hier2$site_cAlaska` + `var.data.hier2$site_cAlps` + `var.data.hier2$site_cCaucasus` + `var.data.hier2$site_cChile` + `var.data.hier2$site_cEcuador` + `var.data.hier2$site_cGreenland` + `var.data.hier2$site_cKirghizistan` + `var.data.hier2$site_cNepal` + `var.data.hier2$site_cNew_Zealand`, data = as.data.frame(region.dum), distance = "bray")
#             Df SumOfSqs      F Pr(>F)    
#   Model      9   17.434 7.5196  0.001 ***
#   Residual 125   32.201                  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


R2a <- RsquareAdj(rda1.S)$adj.r.squared #[1] 0.3045347


# Environmental model
rda1.E <- dbrda(bio.b~., data=env.mat, distance="bray")
vif.cca(rda1.E)
# water_temp         pH       cond       turb        doc        srp  mean_chla   gl_cov.x    Calcite 
# 1.334795   1.309891   2.070996   1.665702   1.233989   1.434614   1.196329   1.376692   1.771330 
# Quartz         pr        scd        DIN 
# 1.218167   1.390589   1.400890   1.435401

# library(corrplot)
# library(RColorBrewer)
# library(PerformanceAnalytics)
# chart.Correlation(as.data.frame(env.mat))

eigenvals(rda1.E)[1:2]/sum(eigenvals(rda1.E))*RsquareAdj(rda1.E)$adj.r.squared/RsquareAdj(rda1.E)$r.squared
# dbRDA1     dbRDA2 
# 0.06280750 0.03730144 
anova(rda1.E, step=1000) 
# Permutation test for dbrda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: dbrda(formula = bio.b ~ water_temp + pH + cond + turb + doc + srp + mean_chla + gl_cov.x + Calcite + Quartz + pr + scd + DIN, data = env.mat, distance = "bray")
# Df SumOfSqs      F Pr(>F)    
# Model     13   14.053 3.6761  0.001 ***
#   Residual 121   35.582                  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

R2a <- RsquareAdj(rda1.E)$adj.r.squared # [1] 0.2061079
##Compute forward selection
mod0 <- dbrda(bio.b ~ 1, data=env.mat, distance="bray")
rda1.E.fwd <- ordistep(mod0, scope=formula(rda1.E), direction="forward", perm.max=200)
rda1.E.fwd$call

rda2 = dbrda(formula = bio.b ~ pH + pr + DIN + turb + cond + mean_chla + 
               water_temp + scd + Calcite + srp + Quartz, data = env.mat, 
             distance = "bray")

anova(rda2, step=1000)


# Permutation test for dbrda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: dbrda(formula = bio.b ~ pH + pr + DIN + turb + cond + mean_chla + water_temp + scd + Calcite + srp + Quartz, data = env.mat, distance = "bray")
# Df SumOfSqs      F Pr(>F)    
# Model     11   13.323 4.1025  0.001 ***
#   Residual 123   36.313                  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



anova(rda2, by="terms")


# Permutation test for dbrda under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# Model: dbrda(formula = bio.b ~ pH + pr + DIN + turb + cond + mean_chla + water_temp + scd + Calcite + srp + Quartz, data = env.mat, distance = "bray")
# Df SumOfSqs       F Pr(>F)    
#   pH           1    3.131 10.6039  0.001 ***
#   pr           1    2.403  8.1386  0.001 ***
#   DIN          1    1.466  4.9654  0.001 ***
#   turb         1    1.239  4.1978  0.001 ***
#   cond         1    1.236  4.1853  0.001 ***
#   mean_chla    1    0.888  3.0065  0.001 ***
#   water_temp   1    0.804  2.7220  0.001 ***
#   scd          1    0.756  2.5592  0.002 ** 
#   Calcite      1    0.554  1.8771  0.004 ** 
#   srp          1    0.420  1.4225  0.058 .  
#   Quartz       1    0.428  1.4495  0.052 .  
# Residual   123   36.313  

anova(rda2, by="axis")

# Permutation test for dbrda under reduced model
# Forward tests for axes
# Permutation: free
# Number of permutations: 999
# 
# Model: dbrda(formula = bio.b ~ pH + pr + DIN + turb + cond + mean_chla + water_temp + scd + Calcite + srp + Quartz, data = env.mat, distance = "bray")
# Df SumOfSqs       F Pr(>F)    
# dbRDA1     1    4.265 14.4471  0.001 ***
# dbRDA2     1    2.479  8.3966  0.001 ***
# dbRDA3     1    1.791  6.0667  0.001 ***
# dbRDA4     1    1.486  5.0320  0.001 ***
# dbRDA5     1    0.888  3.0073  0.002 ** 
# dbRDA6     1    0.801  2.7132  0.001 ***
# dbRDA7     1    0.437  1.4810  0.481    
# dbRDA8     1    0.393  1.3312  0.565    
# dbRDA9     1    0.320  1.0856  0.828    
# dbRDA10    1    0.262  0.8889  0.947    
# dbRDA11    1    0.200  0.6783  0.985    
# Residual 123   36.313  

R2a <- RsquareAdj(rda2)$adj.r.squared #[1] 0.2029876
eigenvals(rda2)[1:2]/sum(eigenvals(rda2))*RsquareAdj(rda2)$adj.r.squared/RsquareAdj(rda2)$r.squared
# dbRDA1     dbRDA2 
# 0.06498396 0.03776829 

##then subset environmental variables to enter it in the variance partitioning
env.mat.adj <- subset(env.mat, select = -c(doc,gl_cov.x))

vp1 <- varpart(vegdist(bio.b), xyz.dat, region.dum, region.dbMEM, env.mat.adj)
vp1;plot(vp1)


# Partition of squared Bray distance in dbRDA 
# 
# Call: varpart(Y = vegdist(bio.b), X = xyz.dat, region.dum, region.dbMEM, env.mat.adj)
# 
# Explanatory tables:
#   X1:  xyz.dat
# X2:  region.dum
# X3:  region.dbMEM
# X4:  env.mat.adj 
# 
# No. of explanatory tables: 4 
# Total variation (SS): 49.636 
# No. of observations: 135 
# 
# Partition table:
#   Df R.square Adj.R.square Testable
# [aeghklno] = X1              3  0.16214      0.14296     TRUE
# [befiklmo] = X2              9  0.35125      0.30453     TRUE
# [cfgjlmno] = X3             60  0.39756     -0.09091     TRUE
# [dhijkmno] = X4             11  0.26841      0.20299     TRUE
# [abefghiklmno] = X1+X2      12  0.38839      0.32823     TRUE
# [acefghjklmno] = X1+X3      63  0.55904      0.16776     TRUE
# [adeghijklmno] = X1+X4      14  0.35298      0.27749     TRUE
# [bcefgijklmno] = X2+X3      69  0.74880      0.48214     TRUE
# [bdefhijklmno] = X2+X4      20  0.45268      0.35666     TRUE
# [cdfghijklmno] = X3+X4      71  0.66084      0.27861     TRUE
# [abcefghijklmno] = X1+X2+X3 72  0.76432      0.49062     TRUE
# [abdefghijklmno] = X1+X2+X4 23  0.48332      0.37626     TRUE
# [acdefghijklmno] = X1+X3+X4 74  0.71694      0.36784     TRUE
# [bcdefghijklmno] = X2+X3+X4 80  0.79623      0.49434     TRUE
# [abcdefghijklmno] = All     83  0.81150      0.50472     TRUE
# Individual fractions                                         
# [a] = X1 | X2+X3+X4          3               0.01038     TRUE
# [b] = X2 | X1+X3+X4          9               0.13688     TRUE
# [c] = X3 | X1+X2+X4         60               0.12846     TRUE
# [d] = X4 | X1+X2+X3         11               0.01410     TRUE
# [e]                          0               0.07885    FALSE
# [f]                          0              -0.03812    FALSE
# [g]                          0               0.00921    FALSE
# [h]                          0              -0.00191    FALSE
# [i]                          0               0.18598    FALSE
# [j]                          0               0.03392    FALSE
# [k]                          0               0.17135    FALSE
# [l]                          0              -0.02394    FALSE
# [m]                          0              -0.09946    FALSE
# [n]                          0               0.00601    FALSE
# [o]                          0              -0.10700    FALSE
# [p] = Residuals              0               0.49528    FALSE
# Controlling 2 tables X                                       
# [ae] = X1 | X3+X4            3               0.08923     TRUE
# [ag] = X1 | X2+X4            3               0.01959     TRUE
# [ah] = X1 | X2+X3            3               0.00847     TRUE
# [be] = X2 | X3+X4            9               0.21574     TRUE
# [bf] = X2 | X1+X4            9               0.09876     TRUE
# [bi] = X2 | X1+X3            9               0.32286     TRUE
# [cf] = X3 | X1+X4           60               0.09035     TRUE
# [cg] = X3 | X2+X4           60               0.13768     TRUE
# [cj] = X3 | X1+X2           60               0.16238     TRUE
# [dh] = X4 | X2+X3           11               0.01220     TRUE
# [di] = X4 | X1+X3           11               0.20008     TRUE
# [dj] = X4 | X1+X2           11               0.04802     TRUE
# Controlling 1 table X                                        
# [aghn] = X1 | X2             3               0.02370     TRUE
# [aehk] = X1 | X3             3               0.25867     TRUE
# [aegl] = X1 | X4             3               0.07451     TRUE
# [bfim] = X2 | X1             9               0.18528     TRUE
# [beik] = X2 | X3             9               0.57306     TRUE
# [befl] = X2 | X4             9               0.15368     TRUE
# [cfjm] = X3 | X1            60               0.02481     TRUE
# [cgjn] = X3 | X2            60               0.17761     TRUE
# [cfgl] = X3 | X4            60               0.07562     TRUE
# [dijm] = X4 | X1            11               0.13454     TRUE
# [dhjn] = X4 | X2            11               0.05213     TRUE
# [dhik] = X4 | X3            11               0.36952     TRUE


# Test of individual fractions by mean of RDA
rda1.Sx <- dbrda(bio.b ~. +Condition(as.matrix(xyz.dat)), data=as.data.frame(region.dum), distance="bray")
anova(rda1.Sx, step=1000)
# Permutation test for dbrda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: dbrda(formula = bio.b ~ `var.data.hier2$site_cAlaska` + `var.data.hier2$site_cAlps` + `var.data.hier2$site_cCaucasus` + `var.data.hier2$site_cChile` + `var.data.hier2$site_cEcuador` + `var.data.hier2$site_cGreenland` + `var.data.hier2$site_cKirghizistan` + `var.data.hier2$site_cNepal` + `var.data.hier2$site_cNew_Zealand` + Condition(as.matrix(xyz.dat)), data = as.data.frame(region.dum), distance = "bray")
# Df SumOfSqs      F Pr(>F)    
# Model      9   11.230 5.0145  0.001 ***
#   Residual 122   30.358                  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1               
R2a <- RsquareAdj(rda1.Sx)$adj.r.squared #  [1] 0.1852785

rda1.Ex <- dbrda(bio.b ~. + Condition(as.matrix(region.dum)) + Condition(as.matrix(region.dbMEM)) +Condition(as.matrix(xyz.dat)), data=as.data.frame(env.mat.adj),distance="bray")
anova(rda1.Ex, step=1000)
(R2a <- RsquareAdj(rda1.Ex)$adj.r.squared) # 

# 3b) Within region scale
rda2.S <- dbrda(bio.b ~. + Condition(as.matrix(region.dum))+Condition(as.matrix(xyz.dat)),data=as.data.frame(region.dbMEM), distance="bray")
anova(rda2.S, step=1000)
# Df SumOfSqs      F Pr(>F)    
# Model    60   18.659 1.6482  0.001 ***
#   Residual 62   11.698                  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

R2a <- RsquareAdj(rda2.S)$adj.r.squared #[1]0.1623823


