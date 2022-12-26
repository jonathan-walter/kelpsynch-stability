rm(list=ls())


library(ncdf4)
library(ncf)
library(nlme)
library(wsyn)

source("/Users/jonathanwalter/Documents/Research/kelpStability/partialSpearman.R")

dat.raw <- nc_open("/Users/jonathanwalter/Documents/Research/kelpStability/CAkelpCanopyEnv_2021.nc")
print(dat.raw)

# xcoord <- c(ncvar_get(dat.raw, "utm_x"))
# ycoord <- c(ncvar_get(dat.raw, "utm_y"))
latitude <-c(ncvar_get(dat.raw, "lat"))
longitude <-c(ncvar_get(dat.raw, "lon"))
year <- c(ncvar_get(dat.raw, "year"))
quarter <- c(ncvar_get(dat.raw, "quarter"))
biomass <- ncvar_get(dat.raw, "biomass")
#cover <- ncvar_get(dat.raw, "area")/900
nitrate <- ncvar_get(dat.raw, "nitrate")
waves <- ncvar_get(dat.raw, "hsmax")

#rm(dat.raw)


## filter data

meetsNAcrit <- apply(biomass, 1, function(x){mean(is.na(x)) < 0.3})
meetsPresCrit <-apply(biomass, 1, function(x){mean(x>0, na.rm=T) > 0.3})
meetsPresCrit[is.na(meetsPresCrit)] <- FALSE

biomass <- biomass[meetsNAcrit & meetsPresCrit,]
nitrate <- nitrate[meetsNAcrit & meetsPresCrit,]
waves <- waves[meetsNAcrit & meetsPresCrit,]
lon <- longitude[meetsNAcrit & meetsPresCrit]
lat <- latitude[meetsNAcrit & meetsPresCrit]

#point conception is lat 34.44881; san francisco bay north edge  is lat 37.815483

#norcal <- lat > 37.815483
sbchannel <- lat < 34.44881 & lat > 34.00000
cencal <- lat >= 34.44881 & lat <= 37.815483
socal <- lat < 34.00000

year2 <- unique(year)

biomass.ann <- t(apply(biomass, 1, function(x){aggregate(x, by=list(year), FUN=max, na.rm=TRUE)$x}))
biomass.ann[is.infinite(biomass.ann)] <- NA
waves.ann <- t(apply(waves, 1, function(x){aggregate(x, by=list(year), FUN=max, na.rm=TRUE)$x}))
nitrate.ann <- t(apply(nitrate, 1, function(x){aggregate(x, by=list(year), FUN=mean, na.rm=TRUE)$x}))

#clean form of data for subsequent wavelet analyses
biomass.ann.cln <- biomass.ann
#fill missing values with medians
if(any(is.na(biomass.ann.cln))){
  for(aa in 1:nrow(biomass.ann.cln)){
    for(bb in 1:ncol(biomass.ann.cln)){
      if(is.na(biomass.ann.cln[aa,bb])){
        biomass.ann.cln[aa,bb] <- median(biomass.ann.cln[aa,], na.rm=T)
      }
    }
  }
}
biomass.ann.cln <- cleandat(biomass.ann.cln, times=1:38, clev=5)$cdat




#measure synchrony, mean exit time, mean recovery time, CV of biomass
biomassCV <- apply(biomass.ann, 1, function(x){sd(x, na.rm=T)/mean(x, na.rm=T)})
wavesCV <- apply(waves.ann, 1, function(x){sd(x, na.rm=T)/mean(x, na.rm=T)})
nitrateCV <- apply(nitrate.ann, 1, function(x){sd(x, na.rm=T)/mean(x, na.rm=T)})

occupancy <- apply(biomass.ann, 1, function(x){mean(x > 0, na.rm=T)})

exitTime <- function(x){
  #remove leading and trailing NAs
  tmp <- rle(is.na(x))
  if(tmp$values[1]){
    x <- x[-c(1:tmp$lengths[1])]
  }
  tmp <- rle(is.na(x))
  if(tmp$values[length(tmp$values)]){
    x <- x[-c(length(x)-tmp$lengths[length(tmp$values)]+1:length(x))]
  }
  
  #consider NAs present if surrounding years are present; otherwise, absent
  x1 <- x > 0
  if(any(is.na(x1))){
    nainds <- which(is.na(x1))
    for(ii in nainds){
      if(any(is.na(c(x1[ii-1], x1[ii+1])))){x1[ii] <- FALSE}
      else if(x1[ii-1] & x1[ii+1]){x1[ii] <- TRUE}
      else{x1[ii] <- FALSE}
    }
  }
  
  tmp <- rle(x1)
  return(mean(tmp$lengths[tmp$values]))
}


kelpExitTime <- apply(biomass.ann, 1, exitTime)

recoveryTime <- function(x){
  #remove leading and trailing NAs
  tmp <- rle(is.na(x))
  if(tmp$values[1]){
    x <- x[-c(1:tmp$lengths[1])]
  }
  tmp <- rle(is.na(x))
  if(tmp$values[length(tmp$values)]){
    x <- x[-c(length(x)-tmp$lengths[length(tmp$values)]+1:length(x))]
  }
  
  x1 <- x > 0
  #consider NAs present if surrounding years are present; otherwise, absent
  if(any(is.na(x1))){
    nainds <- which(is.na(x1))
    for(ii in nainds){
      if(any(is.na(c(x1[ii-1], x1[ii+1])))){x1[ii] <- FALSE}
      else if(x1[ii-1] & x1[ii+1]){x1[ii] <- TRUE}
      else{x1[ii] <- FALSE}
    }
  }
  
  tmp <- rle(x1)
  return(mean(tmp$lengths[!tmp$values]))
  
}

kelpRecoveryTime <- apply(biomass.ann, 1, recoveryTime)


kelpCor <- rep(NA, nrow(biomass.ann))

nsample <- 200 #used for sub-sampling the data to conserve memory and reduce computations
dmax <- 250 #meters; range for quantifying synchrony

for(ii in 1:nrow(biomass.ann)){
  
  toget <- (ii - nsample/2):(ii + nsample/2)
  toget <- toget[toget >= 1 & toget <= nrow(biomass.ann)]
  coords <- cbind(lon[toget],lat[toget])
  foc <- which(toget==ii)
  
  dd <- rep(NA, length(toget))
  for(jj in 1:length(dd)){
    dd[jj] <- gcdist(x=c(coords[foc,1],coords[jj,1]), y=c(coords[foc,2],coords[jj,2]))[2,1]*1000
  }

  #check if nsample is large enough to encompass dmax
  if(max(dd) < dmax){stop("nsample is not large enough to encompass dmax")}
  
  biomass.ii <- biomass.ann[toget,]
  sync <- rep(NA, length(toget))
  for(jj in 1:length(sync)){
    if(jj == foc){next}
    sync[jj] <- cor(biomass.ii[jj,], biomass.ii[foc,], use="pairwise.complete.obs")
  }
  
  sync <- sync[dd < dmax]
  kelpCor[ii] <- mean(sync, na.rm=T)
  
}



kelpXWTst <- rep(NA, nrow(biomass.ann))

for(ii in 1:nrow(biomass.ann.cln)){
  
  toget <- (ii - nsample/2):(ii + nsample/2)
  toget <- toget[toget >= 1 & toget <= nrow(biomass.ann)]
  coords <- cbind(lon[toget],lat[toget])
  foc <- which(toget==ii)
  
  dd <- rep(NA, length(toget))
  for(jj in 1:length(dd)){
    dd[jj] <- gcdist(x=c(coords[foc,1],coords[jj,1]), y=c(coords[foc,2],coords[jj,2]))[2,1]*1000
  }
  
  #check if nsample is large enough to encompass dmax
  if(max(dd) < dmax){stop("nsample is not large enough to encompass dmax")}
  
  biomass.ii <- biomass.ann.cln[toget,]
  
  sync <- rep(NA, length(toget))
  for(jj in 1:length(sync)){
    if(jj == foc){next}
    sync[jj] <- synmat(rbind(biomass.ii[foc,], biomass.ii[jj,]), times=1:38, method="ReXWT", tsrange=c(2,4))[2,1]
  }
  
  sync <- sync[dd < dmax]
  kelpXWTst[ii] <- mean(sync, na.rm=T)
  
}


kelpXWTlt <- rep(NA, nrow(biomass.ann))

for(ii in 1:nrow(biomass.ann.cln)){
  
  toget <- (ii - nsample/2):(ii + nsample/2)
  toget <- toget[toget >= 1 & toget <= nrow(biomass.ann)]
  coords <- cbind(lon[toget],lat[toget])
  foc <- which(toget==ii)
  
  dd <- rep(NA, length(toget))
  for(jj in 1:length(dd)){
    dd[jj] <- gcdist(x=c(coords[foc,1],coords[jj,1]), y=c(coords[foc,2],coords[jj,2]))[2,1]*1000
  }
  
  #check if nsample is large enough to encompass dmax
  if(max(dd) < dmax){stop("nsample is not large enough to encompass dmax")}
  
  biomass.ii <- biomass.ann.cln[toget,]
  
  sync <- rep(NA, length(toget))
  for(jj in 1:length(sync)){
    if(jj == foc){next}
    sync[jj] <- synmat(rbind(biomass.ii[foc,], biomass.ii[jj,]), times=1:38, method="ReXWT", tsrange=c(4,Inf))[2,1]
  }
  
  sync <- sync[dd < dmax]
  kelpXWTlt[ii] <- mean(sync, na.rm=T)
  
}


kelpLTA <- rep(NA, nrow(biomass.ann))
kelpUTA <- rep(NA, nrow(biomass.ann))
kelpTailedness <- rep(NA, nrow(biomass.ann))

for(ii in 1:nrow(biomass.ann)){
  
  toget <- (ii - nsample/2):(ii + nsample/2)
  toget <- toget[toget >= 1 & toget <= nrow(biomass.ann)]
  coords <- cbind(lon[toget],lat[toget])
  foc <- which(toget==ii)
  
  dd <- rep(NA, length(toget))
  for(jj in 1:length(dd)){
    dd[jj] <- gcdist(x=c(coords[foc,1],coords[jj,1]), y=c(coords[foc,2],coords[jj,2]))[2,1]*1000
  }
  
  #check if nsample is large enough to encompass dmax
  if(max(dd) < dmax){stop("nsample is not large enough to encompass dmax")}
  
  biomass.ii <- biomass.ann[toget,]
  
  lta <- rep(NA, length(toget))
  uta <- rep(NA, length(toget))
  for(jj in 1:length(sync)){
    if(jj == foc){next}
    lta[jj] <- partialSpearman(biomass.ii[foc,], biomass.ii[jj,], c(0,0.5))
    uta[jj] <- partialSpearman(biomass.ii[foc,], biomass.ii[jj,], c(0.5,1))
  }
  
  lta <- lta[dd < dmax]
  uta <- uta[dd < dmax]
  tailedness <- uta - lta
  
  kelpLTA[ii] <- mean(lta, na.rm=T)
  kelpUTA[ii] <- mean(uta, na.rm=T)
  kelpTailedness[ii] <- mean(tailedness, na.rm=T)
  
}


## simple figures

png("/Users/jonathanwalter/Documents/Research/kelpStability/stability_vs_pearson.png")
par(mfrow=c(2,2), mar=c(3.1,3.1,1.1,1.1), mgp=c(1.8,0.8,0))
plot(kelpCor, occupancy, xlab="Pearson correlation", ylab="Occupancy", pch=16, cex=0.5)
mtext(paste0("r = ", round(cor(kelpCor, occupancy, use="pairwise.complete.obs"), 3)), line=0.1)
plot(kelpCor, kelpRecoveryTime, xlab="Pearson correlation", ylab="Recovery time", pch=16, cex=0.5)
mtext(paste0("r = ", round(cor(kelpCor, kelpRecoveryTime, use="pairwise.complete.obs"), 3)), line=0.1)
plot(kelpCor, kelpExitTime, xlab="Pearson correlation", ylab="Exit time", pch=16, cex=0.5)
mtext(paste0("r = ", round(cor(kelpCor, kelpExitTime, use="pairwise.complete.obs"), 3)), line=0.1)
plot(kelpCor, biomassCV, xlab="Pearson correlation", ylab="CV(biomass)", pch=16, cex=0.5)
mtext(paste0("r = ", round(cor(kelpCor, biomassCV, use="pairwise.complete.obs"), 3)), line=0.1)
dev.off()


png("/Users/jonathanwalter/Documents/Research/kelpStability/stability_vs_wavst.png")
par(mfrow=c(2,2), mar=c(3.1,3.1,1.1,1.1), mgp=c(1.8,0.8,0))
plot(kelpXWTst, occupancy, xlab="Short timescales", ylab="Occupancy", pch=16, cex=0.5)
mtext(paste0("r = ", round(cor(kelpXWTst, occupancy, use="pairwise.complete.obs"), 3)), line=0.1)
plot(kelpXWTst, kelpRecoveryTime, xlab="Short timescales", ylab="Recovery time", pch=16, cex=0.5)
mtext(paste0("r = ", round(cor(kelpXWTst, kelpRecoveryTime, use="pairwise.complete.obs"), 3)), line=0.1)
plot(kelpXWTst, kelpExitTime, xlab="Short timescales", ylab="Exit time", pch=16, cex=0.5)
mtext(paste0("r = ", round(cor(kelpXWTst, kelpExitTime, use="pairwise.complete.obs"), 3)), line=0.1)
plot(kelpXWTst, biomassCV, xlab="Short timescales", ylab="CV(biomass)", pch=16, cex=0.5)
mtext(paste0("r = ", round(cor(kelpXWTst, biomassCV, use="pairwise.complete.obs"), 3)), line=0.1)
dev.off()


png("/Users/jonathanwalter/Documents/Research/kelpStability/stability_vs_wavlt.png")
par(mfrow=c(2,2), mar=c(3.1,3.1,1.1,1.1), mgp=c(1.8,0.8,0))
plot(kelpXWTlt, occupancy, xlab="Long timescales", ylab="Occupancy", pch=16, cex=0.5)
mtext(paste0("r = ", round(cor(kelpXWTlt, occupancy, use="pairwise.complete.obs"), 3)), line=0.1)
plot(kelpXWTlt, kelpRecoveryTime, xlab="Long timescales", ylab="Recovery time", pch=16, cex=0.5)
mtext(paste0("r = ", round(cor(kelpXWTlt, kelpRecoveryTime, use="pairwise.complete.obs"), 3)), line=0.1)
plot(kelpXWTlt, kelpExitTime, xlab="Long timescales", ylab="Exit time", pch=16, cex=0.5)
mtext(paste0("r = ", round(cor(kelpXWTlt, kelpExitTime, use="pairwise.complete.obs"), 3)), line=0.1)
plot(kelpXWTlt, biomassCV, xlab="Long timescales", ylab="CV(biomass)", pch=16, cex=0.5)
mtext(paste0("r = ", round(cor(kelpXWTlt, biomassCV, use="pairwise.complete.obs"), 3)), line=0.1)
dev.off()

png("/Users/jonathanwalter/Documents/Research/kelpStability/stability_vs_lta.png")
par(mfrow=c(2,2), mar=c(3.1,3.1,1.1,1.1), mgp=c(1.8,0.8,0))
plot(kelpLTA, occupancy, xlab="Lower-tail synchrony", ylab="Occupancy", pch=16, cex=0.5)
mtext(paste0("r = ", round(cor(kelpLTA, occupancy, use="pairwise.complete.obs"), 3)), line=0.1)
plot(kelpLTA, kelpRecoveryTime, xlab="Lower-tail synchrony", ylab="Recovery time", pch=16, cex=0.5)
mtext(paste0("r = ", round(cor(kelpLTA, kelpRecoveryTime, use="pairwise.complete.obs"), 3)), line=0.1)
plot(kelpLTA, kelpExitTime, xlab="Lower-tail synchrony", ylab="Exit time", pch=16, cex=0.5)
mtext(paste0("r = ", round(cor(kelpLTA, kelpExitTime, use="pairwise.complete.obs"), 3)), line=0.1)
plot(kelpLTA, biomassCV, xlab="Lower-tail synchrony", ylab="CV(biomass)", pch=16, cex=0.5)
mtext(paste0("r = ", round(cor(kelpLTA, biomassCV, use="pairwise.complete.obs"), 3)), line=0.1)
dev.off()

# quartz()
# par(mfrow=c(2,2), mar=c(4.1,4.1,3.1,1.1))
# plot(kelpSynchrony, occupancy, xlab="Synchrony", ylab="Occupancy")
# mtext("All California", line=2)
# mtext(paste0("r = ", round(cor(kelpSynchrony, occupancy, use="pairwise.complete.obs"), 2)), line=1)
# 
# plot(kelpSynchrony[cencal], occupancy[cencal], xlab="Synchrony", ylab="Occupancy")
# mtext("Central California", line=2)
# mtext(paste0("r = ", round(cor(kelpSynchrony[cencal], occupancy[cencal], use="pairwise.complete.obs"), 2)), line=1)
# 
# plot(kelpSynchrony[sbchannel], occupancy[sbchannel], xlab="Synchrony", ylab="Occupancy")
# mtext("Santa Barbara Channel", line=2)
# mtext(paste0("r = ", round(cor(kelpSynchrony[sbchannel], occupancy[sbchannel], use="pairwise.complete.obs"), 2)), line=1)
# 
# plot(kelpSynchrony[socal], occupancy[socal], xlab="Synchrony", ylab="Occupancy")
# mtext("Southern California", line=2)
# mtext(paste0("r = ", round(cor(kelpSynchrony[socal], occupancy[socal], use="pairwise.complete.obs"), 2)), line=1)
# 
# 
# 
# quartz()
# par(mfrow=c(2,2), mar=c(4.1,4.1,3.1,1.1))
# plot(kelpSynchrony, biomassCV, xlab="Synchrony", ylab="CV(Biomass)")
# mtext("All California", line=2)
# mtext(paste0("r = ", round(cor(kelpSynchrony, biomassCV, use="pairwise.complete.obs"), 2)), line=1)
# 
# plot(kelpSynchrony[cencal], biomassCV[cencal], xlab="Synchrony", ylab="CV(Biomass)")
# mtext("Central California", line=2)
# mtext(paste0("r = ", round(cor(kelpSynchrony[cencal], biomassCV[cencal], use="pairwise.complete.obs"), 2)), line=1)
# 
# plot(kelpSynchrony[sbchannel], biomassCV[sbchannel], xlab="Synchrony", ylab="CV(Biomass)")
# mtext("Santa Barbara Channel", line=2)
# mtext(paste0("r = ", round(cor(kelpSynchrony[sbchannel], biomassCV[sbchannel], use="pairwise.complete.obs"), 2)), line=1)
# 
# plot(kelpSynchrony[socal], biomassCV[socal], xlab="Synchrony", ylab="CV(Biomass)")
# mtext("Southern California", line=2)
# mtext(paste0("r = ", round(cor(kelpSynchrony[socal], biomassCV[socal], use="pairwise.complete.obs"), 2)), line=1)
# 
# 
# 
# quartz()
# par(mfrow=c(2,2), mar=c(4.1,4.1,3.1,1.1))
# plot(kelpSynchrony, kelpRecoveryTime, xlab="Synchrony", ylab="Recovery time")
# mtext("All California", line=2)
# mtext(paste0("r = ", round(cor(kelpSynchrony, kelpRecoveryTime, use="pairwise.complete.obs"), 2)), line=1)
# 
# plot(kelpSynchrony[cencal], kelpRecoveryTime[cencal], xlab="Synchrony", ylab="Recovery time")
# mtext("Central California", line=2)
# mtext(paste0("r = ", round(cor(kelpSynchrony[cencal], kelpRecoveryTime[cencal], use="pairwise.complete.obs"), 2)), line=1)
# 
# plot(kelpSynchrony[sbchannel], kelpRecoveryTime[sbchannel], xlab="Synchrony", ylab="Recovery time")
# mtext("Santa Barbara Channel", line=2)
# mtext(paste0("r = ", round(cor(kelpSynchrony[sbchannel], kelpRecoveryTime[sbchannel], use="pairwise.complete.obs"), 2)), line=1)
# 
# plot(kelpSynchrony[socal], kelpRecoveryTime[socal], xlab="Synchrony", ylab="Recovery time")
# mtext("Southern California", line=2)
# mtext(paste0("r = ", round(cor(kelpSynchrony[socal], kelpRecoveryTime[socal], use="pairwise.complete.obs"), 2)), line=1)
# 
# 
# 
# par(mfrow=c(2,2), mar=c(4.1,4.1,3.1,1.1))
# plot(kelpSynchrony, kelpExitTime, xlab="Synchrony", ylab="Exit time")
# mtext("All California", line=2)
# mtext(paste0("r = ", round(cor(kelpSynchrony, kelpExitTime, use="pairwise.complete.obs"), 2)), line=1)
# 
# plot(kelpSynchrony[cencal], kelpExitTime[cencal], xlab="Synchrony", ylab="Exit time")
# mtext("Central California", line=2)
# mtext(paste0("r = ", round(cor(kelpSynchrony[cencal], kelpExitTime[cencal], use="pairwise.complete.obs"), 2)), line=1)
# 
# plot(kelpSynchrony[sbchannel], kelpExitTime[sbchannel], xlab="Synchrony", ylab="Exit time")
# mtext("Santa Barbara Channel", line=2)
# mtext(paste0("r = ", round(cor(kelpSynchrony[sbchannel], kelpExitTime[sbchannel], use="pairwise.complete.obs"), 2)), line=1)
# 
# plot(kelpSynchrony[socal], kelpExitTime[socal], xlab="Synchrony", ylab="Exit time")
# mtext("Southern California", line=2)
# mtext(paste0("r = ", round(cor(kelpSynchrony[socal], kelpExitTime[socal], use="pairwise.complete.obs"), 2)), line=1)



## statistical models also testing for effects of variability in environmental conditions

moddf <- data.frame(occupancy = occupancy,
                    biomassCV = biomassCV,
                    exitTime = kelpExitTime,
                    recoveryTime = kelpRecoveryTime,
                    sync.pearson = scale(kelpCor), #note predictor variables are scaled 
                    sync.st = scale(kelpXWTst),
                    sync.lt = scale(kelpXWTlt),
                    lta = scale(kelpLTA),
                    uta = scale(kelpUTA),
                    tailedness = scale(kelpTailedness),
                    nitrateCV = scale(nitrateCV),
                    wavesCV = scale(wavesCV),
                    lat = lat,
                    lon = lon
                    )

write.csv(moddf, "/Users/jonathanwalter/Documents/Research/kelpStability/modelvars.csv", row.names=FALSE)

moddf <- moddf[complete.cases(moddf),]

## pearson correlation
mod.occ <- lm(occupancy ~ sync.pearson + nitrateCV + wavesCV, data=moddf)
summary(mod.occ)

mod.biomassCV <- lm(biomassCV ~ sync.pearson + nitrateCV + wavesCV, data=moddf)
summary(mod.biomassCV)

mod.exitTime <- lm(exitTime ~ sync.pearson + nitrateCV + wavesCV, data=moddf)
summary(mod.exitTime)

mod.recoveryTime <- lm(recoveryTime ~ sync.pearson + nitrateCV + wavesCV, data=moddf)
summary(mod.recoveryTime)


## short timescale synchrony
mod.occ <- lm(occupancy ~ sync.st + nitrateCV + wavesCV, data=moddf)
summary(mod.occ)

mod.biomassCV <- lm(biomassCV ~ sync.st + nitrateCV + wavesCV, data=moddf)
summary(mod.biomassCV)

mod.exitTime <- lm(exitTime ~ sync.st + nitrateCV + wavesCV, data=moddf)
summary(mod.exitTime)

mod.recoveryTime <- lm(recoveryTime ~ sync.st + nitrateCV + wavesCV, data=moddf)
summary(mod.recoveryTime)


#long timescale synchrony
mod.occ <- lm(occupancy ~ sync.lt + nitrateCV + wavesCV, data=moddf)
summary(mod.occ)

mod.biomassCV <- lm(biomassCV ~ sync.lt + nitrateCV + wavesCV, data=moddf)
summary(mod.biomassCV)

mod.exitTime <- lm(exitTime ~ sync.lt + nitrateCV + wavesCV, data=moddf)
summary(mod.exitTime)

mod.recoveryTime <- lm(recoveryTime ~ sync.lt + nitrateCV + wavesCV, data=moddf)
summary(mod.recoveryTime)


#lta
mod.occ <- lm(occupancy ~ lta + nitrateCV + wavesCV, data=moddf)
summary(mod.occ)

mod.biomassCV <- lm(biomassCV ~ lta + nitrateCV + wavesCV, data=moddf)
summary(mod.biomassCV)

mod.exitTime <- lm(exitTime ~ lta + nitrateCV + wavesCV, data=moddf)
summary(mod.exitTime)

mod.recoveryTime <- lm(recoveryTime ~ lta + nitrateCV + wavesCV, data=moddf)
summary(mod.recoveryTime)


## does tailedness improve models?
mod.occ <- lm(occupancy ~ sync.pearson + nitrateCV + wavesCV, data=moddf)
mod.occ2 <- lm(occupancy ~ sync.pearson + tailedness + nitrateCV + wavesCV, data=moddf)
anova(mod.occ, mod.occ2)
summary(mod.occ2)

mod.biomassCV <- lm(biomassCV ~ sync.pearson + nitrateCV + wavesCV, data=moddf)
mod.biomassCV2 <- lm(biomassCV ~ sync.pearson + tailedness + nitrateCV + wavesCV, data=moddf)
anova(mod.biomassCV, mod.biomassCV2)
summary(mod.biomassCV2)

mod.exitTime <- lm(exitTime ~ sync.pearson + nitrateCV + wavesCV, data=moddf)
mod.exitTime2 <- lm(exitTime ~ sync.pearson + tailedness + nitrateCV + wavesCV, data=moddf)
anova(mod.exitTime, mod.exitTime2)
summary(mod.exitTime2)

mod.recoveryTime <- lm(recoveryTime ~ sync.pearson + nitrateCV + wavesCV, data=moddf)
mod.recoveryTime2 <- lm(recoveryTime ~ sync.pearson + tailedness + nitrateCV + wavesCV, data=moddf)
anova(mod.recoveryTime, mod.recoveryTime2)


save.image("/Users/jonathanwalter/Documents/Research/kelpStability/workspace.RData")
