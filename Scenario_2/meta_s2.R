meta.s2_pop <- data.frame(rep(list(run=numeric(years)), meta.runs)) 
meta.s2_lek.size <- data.frame(rep(list(run=numeric(years)), meta.runs)) 
meta.s2_n.leks <- data.frame(rep(list(run=numeric(years)), meta.runs))
meta.s2_max <- data.frame(rep(list(run=numeric(years)), meta.runs))
meta.s2_max_all <- data.frame(rep(list(run=numeric(years)), meta.runs)) 

maps <- c(1:5)

for (i in 1:meta.runs){

sm <- sample(maps,1)

map <- raster(paste("s2_large_map_",sm, ".asc", sep=""), crs=BNG)

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

#meta.s2_full[[i]] <- meta.1

meta.s2_pop[,i] <- meta.1$tot.pop

meta.s2_lek.size[,i] <- meta.1$lek.size

meta.s2_n.leks[,i] <- meta.1$n.leks

meta.s2_max[,i] <- meta.1$max.avg

meta.s2_max_all[,i] <- meta.1$maxent

print(i)
}

write.csv(meta.s2_pop, "Scenario_2_pop.csv",row.names=F)
write.csv(meta.s2_lek.size, "Scenario_2_lek_size.csv",row.names=F)
write.csv(meta.s2_n.leks, "Scenario_2_leks.csv",row.names=F)
write.csv(meta.s1_max, "Scenario_2_max.csv", row.names=F)
write.csv(meta.s1_max_all, "Scenario_2_max_all.csv", row.names=F)
save.image("Scenario_2.RData")
