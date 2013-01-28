require(sp)
require(raster)
require(dismo)
require(msm)

# Set initial conditions
##########################
BNG<- CRS("+init=epsg:27700")

potential <- raster("studyareadem.asc", crs=BNG)

rcl.pot <-  matrix(c(-1, 99, NA, 100, 650, 0, 651, 2500, NA), nrow=3, ncol=3, byrow=T)

potential <- reclass(potential, rcl.pot)

potential <- potential + raster("maxent.asc")

writeRaster(potential, "potential.asc", overwrite=T)

rm(potential)

meta.runs <- 100 #Number of repeat simulations

leks <- 500 #Initial number of leks

years <- 50 #Length of simulation

K <- 40 #Carrying capacity per lek

fec <- rtnorm(1000, mean=1.7, sd=0.25, lower=0)

# Note that survial is a data frame rather than ditribution so that the distribution for each lek site can be scaled by MaxEnt score
surv <- rtnorm(n=1000, mean=0.66, sd=0.25, lower=0, upper=1)

juv.surv <- rtnorm(1000, mean=0.56, sd=0.25, lower=0, upper=1)

max.mean <- mean(extract(raster("maxent.asc"), read.csv("BK_pres.csv")), na.rm=T)

cp <- read.csv("changespred.csv") #Leks in 1994
cp <- cp$count #Lek size in 1994

lek.dist <- rtnorm(1000, mean=mean(cp), sd=sd(cp), lower=0) #Create intial distribution from lek sizes in 1994
###########################

# Define metapopulation function
################################
meta <- function(N0, d0, f, surv, js, years, xy, maxent, p.map, K){
  popn <- matrix(0, nrow=length(N0), ncol=years)
  dispersal <- matrix(0, nrow=length(N0), ncol=years)
  fecundity <- matrix(0, nrow=length(N0), ncol=years)
  chicks <- matrix(0, nrow=length(N0), ncol=years)
  left <- matrix(0, nrow=length(N0), ncol=years)
  disp.prop <- matrix(0, nrow=length(N0), ncol=years)
  max.avg <- numeric(years)
popn[,1] <- N0
dispersal[,1] <- d0
if(length(maxent[which(popn[,1] > 0)]) >= 1){
		max.avg[1] <- mean(maxent[which(popn[,1] > 0)])
	} else {
		max.avg[1] <- 0
	}
for(i in 2:years-1){
  N <- popn[,i]
  for(j in 1:length(N0)){
	N[j] <- N[j] * sample(surv, 1)
	N[j] <- N[j] + dispersal[j,i] * sample(js, 1)
	lek.fec <- sample(f,1)/2
   	d.j <- N[j] * lek.fec
    	popn[j,i + 1] <- round(N[j], 0) * abs(1 - N[i]/K)
    	chicks[j,i] <- round(d.j, 0)
    	fecundity[j,i] <- lek.fec
	if(length(maxent[which(popn[,i] > 0)]) >= 1){
		max.avg[i] <- mean(maxent[which(popn[,i] > 0)])
	} else {
		max.avg[i] <- 0
	}
}
for (k in 2:length(popn[,1])) {
		disp.run <- rep(0, length(popn[,1]))
		p <- SpatialPoints(xy[k,], proj4string=BNG)
		d <- spDists(p.map, p,longlat=F) 
		d.which <- which(d < 15000)
  choose.disp <- 1:length(d.which)
  d.which <- d.which[sample(choose.disp, (length(choose.disp)*0.4))] # Arbitrary number of leks to disperse to
		d.dist <- d[d.which]+1
		dist.cont <- sum(d.dist)/d.dist
		dist.scale <- 1/dist.cont
		d.count <- popn[d.which,i]
		d.males <- chicks[d.which,i]
		d.count <- d.count + d.males
		d.count <- d.count * (1-d.count/K)
		for (m in 1:length(d.count)){
		if (d.count[m] < 0) {
			d.count[m] <- 0
			}
		}
		d.count <- sqrt(d.count) + 0.001
		d.count <- d.count/sum(d.count)
		d.fec <- chicks[k,i]
		disp <- dist.scale + d.count
		disp <- disp/sum(disp)
		disp2 <- round(disp * d.fec , 0)
    while (sum(disp2) > d.fec){
      cor.dist <- 1:length(disp2[disp2>0])
      cor.which <- which(disp2 > 0)
      x <- sample(cor.dist, 1)
      x <- cor.which[x]
      disp2[x] <- disp2[x] -1
      }
    if (d.fec-sum(disp2)> 0){
      dispersal[k, i+1] <- dispersal[k, i+1] + (d.fec - sum(disp2))
    }
    left[k, i] <- d.fec- sum(disp2)
		disp.run[d.which] <- disp2
		dispersal[,i + 1] <- dispersal[,i + 1] + disp.run
}
}
  max.avg[years] <- mean(maxent[which(popn[,years] > 0)])
  tot.pop <- apply(popn, 2, sum)
  tot.chicks <- apply(chicks, 2, sum)
  n.l <- apply(popn > 0, 2, which)
  n.leks <- as.numeric(summary(n.l)[,1])
  lek.size <- tot.pop/n.leks
  n.chicks <- tot.chicks/tot.pop
  return(list(popn=popn, fecundity=fecundity, dispersal=dispersal, chicks=chicks, tot.pop=tot.pop, lek.size=lek.size, n.leks = n.leks, n.chicks=n.chicks, left=left, max.avg=max.avg))
}
############################
#Plot function
####################
plot.meta <- function(meta){
  op <- par()
  par(mfrow=c(2,2), ask=T)
  plot(1:ncol(meta$popn), meta$tot.pop, type="b", , xlab="time", ylab="Population size")
  plot(1:ncol(meta$popn), meta$n.leks, type="b",xlab="time", ylab="Number of leks")
  plot(1:ncol(meta$popn), meta$lek.size, type="b",xlab="time", ylab="Mean lek size")
  plot(1:ncol(meta$popn), meta$n.chicks, type="b",xlab="time", ylab="Mean brood size")
  par(mfrow=c(1,1))
  matplot(1:ncol(meta$popn), t(meta$dispersal), type="l", ylab="Number of incoming juveniles", xlab="time")
  matplot(1:ncol(meta$popn), t(meta$chicks), type="l", ylab="Number of chicks produced", xlab="time")
  matplot(1:ncol(meta$popn), t(meta$popn), type="l", ylab="Lek size", xlab="time")
  par(mfrow=c(1,1), ask=F)
  }
###########################

source("meta_max_nomax.R", echo=T)

