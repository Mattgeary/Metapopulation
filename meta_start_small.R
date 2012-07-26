require(sp)
require(raster)
require(dismo)
require(msm)

work <- "/home/alan/Work/Metapopulation/New_model"

source("meta_control_small.R", echo=T)

source("meta_94_small.R", echo=T)

setwd(paste(work, "Scenario_1", sep="/"))

source("meta_s1_small.R", echo=T)

setwd(work)

setwd(paste(work, "Scenario_2", sep="/"))

source("meta_s2_small.R", echo=T)

setwd(work)

setwd(paste(work, "Scenario_6", sep="/"))

source("meta_s6_small.R", echo=T)

setwd(work)

save.image("All_small.RData")

rm(list=ls())
