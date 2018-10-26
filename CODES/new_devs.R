

for (i in 46:50) {
  setwd("C:\\Users\\marshalle\\Dropbox\\2017_2020_PHD_MELBOURNE\\CHAPTER_2_PVASIMULATION\\PVA_1_SquirrelGlider\\CODES")
  raster_dev <- raster(paste0("../../LARGEDEVELOPMENTS/NEW_DEVELOPMENTS_", i, ".asc"), sep="")
new_devs <- raster_dev * SQG_SDM
writeRaster(new_devs, paste0("SQG_LRGEDEVS_2/SQG_DEVELOPMENTS_NEW_",i,".asc"), sep="_",format = "ascii",overwrite=TRUE)
}


developments <- raster("../../LARGEDEVELOPMENTS/NEW_DEVELOPMENTS_46.asc")

plot(developments)

dev_new <- developments*SQG_SDM
plot(dev_new)

new <- raster("SQG_LRGEDEVS_2/SQG_DEVELOPMENTS_NEW_13.asc")
plot(new)