meta.94_pop <- data.frame(rep(list(run=numeric(years)), meta.runs)) 
meta.94_lek.size <- data.frame(rep(list(run=numeric(years)), meta.runs)) 
meta.94_n.leks <- data.frame(rep(list(run=numeric(years)), meta.runs)) 
meta.94_max <- data.frame(rep(list(run=numeric(years)), meta.runs))
meta.94_max_all <- data.frame(rep(list(run=numeric(years)), meta.runs)) 

for (i in 1:meta.runs){

map <- raster("maxent.asc", crs=BNG)

potential <- raster("potential.asc", crs=BNG)

xy <- randomPoints(potential, leks)
xy <- as.data.frame(xy)
p.map <- SpatialPoints(xy, proj4string=BNG)

maxent <- extract(map, xy)

leks.rnd <- 1:leks
which.leks <- sample(leks.rnd, round(leks/6, 0))
start.leks <- rep(0, leks)
start.leks[which.leks] <- round(sample(lek.dist, length(which.leks)),0)

disp.mean <- (mean(cp)*mean(fec)) * mean(juv.surv)
disp.dist <- rtnorm(1000, mean=disp.mean, sd=sd(cp), lower=0)

disp.rnd <- 1:leks
which.disp <- sample(disp.rnd, round(leks/6, 0))
start.disp <- rep(0, leks)
start.disp[which.disp] <- round(sample(disp.dist, length(which.leks)), 0)

meta.1 <- meta(N0=start.leks, d0=start.disp, f=fec, surv=surv, js=juv.surv, maxent=maxent, p.map=p.map, years=years, xy=xy, K=K)
  
#plot.meta(meta.1)

#meta.94_full[[i]] <- meta.1

meta.94_pop[,i] <- meta.1$tot.pop

meta.94_lek.size[,i] <- meta.1$lek.size

meta.94_n.leks[,i] <- meta.1$n.leks

meta.94_max[,i] <- meta.1$max.avg

meta.94_max_all[,i] <- meta.1$maxent

print(i)
}

write.csv(meta.94_pop, "P94_pop.csv",row.names=F)
write.csv(meta.94_lek.size, "P94_lek_size.csv",row.names=F)
write.csv(meta.94_n.leks, "P94_leks.csv",row.names=F)
write.csv(meta.94_max, "P94_max.csv", row.names=F)
write.csv(meta.94_max_all, "P94_max_all.csv", row.names=F)
save.image("P94.RData")

