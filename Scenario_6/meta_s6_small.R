meta.s6_pop <- data.frame(rep(list(run=numeric(years)), meta.runs)) 
meta.s6_lek.size <- data.frame(rep(list(run=numeric(years)), meta.runs)) 
meta.s6_n.leks <- data.frame(rep(list(run=numeric(years)), meta.runs)) 
meta.s6_max <- data.frame(rep(list(run=numeric(years)), meta.runs)) 

maps <- c(1:5)

for (i in 1:meta.runs){

sm <- sample(maps,1)

map <- raster(paste("s6_large_map_",sm, ".asc", sep=""), proj4=BNG)

potential <- raster("potential.asc", crs=BNG)

xy <- randomPoints(potential, leks)
xy <- as.data.frame(xy)
p.map <- SpatialPoints(xy, proj4string=BNG)

maxent <- extract(map, xy)

leks.rnd <- 1:leks
which.leks <- sample(leks.rnd, round(leks/5, 0))
start.leks <- rep(0, leks)
start.leks[which.leks] <- round(sample(lek.dist, length(which.leks)),0)

disp.mean <- (mean(cp)*mean(fec)) * mean(juv.surv)
disp.dist <- rtnorm(1000, mean=disp.mean, sd=sd(cp), lower=0)

disp.rnd <- 1:leks
which.disp <- sample(disp.rnd, round(leks/5, 0))
start.disp <- rep(0, leks)
start.disp[which.disp] <- round(sample(disp.dist, length(which.leks)), 0)

meta.1 <- meta(N0=start.leks, d0=start.disp, f=fec, surv=surv, js=juv.surv, maxent=maxent, p.map=p.map, years=years, xy=xy, K=K)
  
#plot.meta(meta.1)

#meta.s6_full[[i]] <- meta.1

meta.s6_pop[,i] <- meta.1$tot.pop

meta.s6_lek.size[,i] <- meta.1$lek.size

meta.s6_n.leks[,i] <- meta.1$n.leks

meta.s6_max[,i] <- meta.1$max.avg

print(i)
}

write.csv(meta.s6_pop, "Scenario_6_pop_small.csv", row.names=F)
write.csv(meta.s6_lek.size, "Scenario_6_lek_size_small.csv", row.names=F)
write.csv(meta.s6_n.leks, "Scenario_6_leks_small.csv", row.names=F)
save.image("Scenario_6_small.RData")
