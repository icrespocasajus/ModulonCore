library(usethis)
network.LCMV = readRDS(file="./network.LCMV.Rds")
network.LCMV = unique(network.LCMV)
modulons.LCMV=readRDS(file="./TF.AUC.clusters.LCMV.subset.Rds")
usethis::use_data(network.LCMV, overwrite = TRUE)
usethis::use_data(modulons.LCMV, overwrite = TRUE)