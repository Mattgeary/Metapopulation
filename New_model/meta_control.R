require(sp)
require(raster)
require(dismo)
require(msm)

# Set initial conditions
##########################
BNG<- CRS("+init=epsg:27700")

meta.runs <- 500 #Number of repeat simulations

leks <- 316 #Initial number of leks

years <- 50 #Length of simulation

K <- 40 #Carrying capacity per lek

fec <- rtnorm(1000, mean=1.7, sd=0.25, lower=0)

surv <- rtnorm(1000, mean=0.66, sd=0.25, lower=0, upper=1)

juv.surv <- rtnorm(1000, mean=0.56, sd=0.25, lower=0, upper=1)

cp <- read.csv("changespred.csv") #Leks in 1994
cp <- cp$count #Lek size in 1994

lek.dist <- rtnorm(1000, mean=mean(cp), sd=sd(cp), lower=0) #Create intial distribution from lek sizes in 1994
###########################

# Define metapopulation function
################################
meta <- function(N0, d0, f, s, js, years, xy, maxent, p.map, K){
  popn <- matrix(0, nrow=length(N0), ncol=years)
  dispersal <- matrix(0, nrow=length(N0), ncol=years)
  fecundity <- matrix(0, nrow=length(N0), ncol=years)
  chicks <- matrix(0, nrow=length(N0), ncol=years)
  left <- matrix(0, nrow=length(N0), ncol=years)
  disp.prop <- matrix(0, nrow=length(N0), ncol=years)
popn[,1] <- N0
dispersal[,1] <- d0
for(i in 2:years-1){
  N <- popn[,i]
  for(j in 1:length(N0)){
	N[j] <- N[j] * (sample(s, 1) * (0.5 + maxent[j]))
	N[j] <- N[j] + dispersal[j,i] * sample(js, 1)
	lek.fec <- sample(f,1)/2
   	d.j <- N[j] * lek.fec
    	popn[j,i + 1] <- round(N[j], 0) * (1 - N[i]/K)
    	chicks[j,i] <- round(d.j, 0)
    	fecundity[j,i] <- lek.fec	
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
		d.max <- maxent[d.which]
                d.fec <- chicks[k,i]
		disp <- dist.scale  * d.max + d.count
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
  tot.pop <- apply(popn, 2, sum)
  tot.chicks <- apply(chicks, 2, sum)
  n.l <- apply(popn > 0, 2, which)
  n.leks <- as.numeric(summary(n.l)[,1])
  lek.size <- tot.pop/n.leks
  n.chicks <- tot.chicks/tot.pop
  return(list(popn=popn, fecundity=fecundity, dispersal=dispersal, chicks=chicks, tot.pop=tot.pop, lek.size=lek.size, n.leks = n.leks, n.chicks=n.chicks, left=left))
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


