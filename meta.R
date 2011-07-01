library(sp)
library(raster)
library(dismo)
library(msm)

leks <- 250

fec <- rtnorm(1000, mean=2.5, sd=0, lower=0)

surv <- rtnorm(1000, mean=0.66, sd=0, lower=0, upper=1)

juv.surv <- rtnorm(1000, mean=0.56, sd=0, lower=0, upper=1)

years <- 10

K <- 200

BNG<- CRS("+init=epsg:27700")

map <- read.asciigrid("maxent.asc", proj4=BNG)
map <- raster(map)

xy <- randomPoints(map, leks)
xy <- as.data.frame(xy)
p.map <- SpatialPoints(xy, proj4string=BNG)

maxent <- extract(map, xy)

cp <- read.csv("changespred.csv")
cp <- cp$count

lek.dist <- rtnorm(1000, mean=mean(cp), sd=sd(cp), lower=0)

leks.rnd <- 1:leks
which.leks <- sample(leks.rnd, round(leks/3, 0))
start.leks <- rep(0, leks)
start.leks[which.leks] <- round(sample(lek.dist, length(which.leks)),0)

disp.mean <- (mean(cp)*mean(fec)) * mean(juv.surv)
disp.dist <- rtnorm(1000, mean=disp.mean, sd=sd(cp), lower=0)

disp.rnd <- 1:leks
which.disp <- sample(disp.rnd, round(leks/3, 0))
start.disp <- rep(0, leks)
start.disp[which.disp] <- round(sample(disp.dist, length(which.leks)), 0)

meta <- function(N0, d0, f, s, js, years, xy, maxent, p.map, K){
  popn <- matrix(0, nrow=length(N0), ncol=years)
  dispersal <- matrix(0, nrow=length(N0), ncol=years)
  fecundity <- matrix(0, nrow=length(N0), ncol=years)
  chicks <- matrix(0, nrow=length(N0), ncol=years)
popn[,1] <- N0
dispersal[,1] <- d0
for(i in 2:years-1){
  N <- popn[,i]
  for(j in 1:length(N0)){
	N[j] <- N[j] * (sample(s, 1) * (0.5 + maxent[j]))
	N[j] <- N[j] + dispersal[j,i] * sample(js, 1)
	lek.fec <- (sample(f,1)/2)
   	d.j <- N[j] * lek.fec 
    	popn[j,i + 1] <- round(N[j], 0)
    	chicks[j,i] <- round(d.j, 0)
    	fecundity[j,i] <- lek.fec	
}
for (k in 2:length(popn[,1])) {
		disp.run <- rep(0, length(popn[,1]))
		p <- SpatialPoints(xy[k,], proj4string=BNG)
		d <- spDists(p.map, p,longlat=F) 
		d.which <- which(d < 30000)
		d.dist <- d[d < 30000]
		dist.cont <- sum(d.dist)/d.dist
		dist.scale <- 1/dist.cont
		d.count <- popn[d.which,i]
		d.males <- chicks[d.which,i]
		d.count <- sqrt(d.count/2) + 0.001 + d.males
		d.count <- d.count/sum(d.count)
		d.max <- maxent[d.which]
                d.fec <- chicks[k,i]
		disp <- dist.scale  * d.max + d.count
		disp <- disp/sum(disp)
		disp2 <- round(disp * d.fec , 0)
		disp2[min(disp2)] <- d.fec- sum(disp2)
		disp.run[d.which] <- disp2
		dispersal[,i + 1] <- dispersal[,i + 1] + disp.run
}
}
  tot.pop <- apply(popn, 2, sum)
  tot.chicks <- apply(chicks, 2, sum)
  n.l <- apply(popn > 0, 2, which)
  n.leks <- as.numeric(summary(n.l)[,1])
  lek.size <- tot.pop/n.leks
  n.chicks <- (tot.chicks/n.leks)/tot.pop
  return(list(popn=popn, fecundity=fecundity, dispersal=dispersal, chicks=chicks, tot.pop=tot.pop, lek.size=lek.size, n.leks = n.leks, n.chicks=n.chicks))
}

meta.1 <- meta(N0=start.leks, d0=start.disp, f=fec, s=surv, js=juv.surv, maxent=maxent, p.map=p.map, years=years, xy=xy, K=K)

plot.meta <- function(meta){
op <- par()
par(mfrow=c(2,2), ask=T)
plot(1:ncol(meta$popn), meta$tot.pop, type="b", , xlab="time", ylab="Population size")
plot(1:ncol(meta$popn), meta$n.leks, type="b",xlab="time", ylab="Number of leks")
plot(1:ncol(meta$popn), meta$lek.size, type="b",xlab="time", ylab="Mean lek size")
plot(1:ncol(meta$popn), meta$n.chicks, type="b",xlab="time", ylab="Mean brood size")
par(mfrow=c(1,1))
matplot(1:ncol(meta$popn), t(meta.1$dispersal), type="l", ylab="Number of incoming juveniles", xlab="time")
matplot(1:ncol(meta$popn), t(meta.1$chicks), type="l", ylab="Number of chicks produced", xlab="time")
matplot(1:ncol(meta$popn), t(meta.1$popn), type="l", ylab="Lek size", xlab="time")
par(mfrow=c(1,1), ask=F)
}

plot.meta(meta.1)
