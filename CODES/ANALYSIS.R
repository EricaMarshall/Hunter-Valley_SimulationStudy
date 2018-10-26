source("OFFSET_function.R")
source("VEGETATION_function_3.R")
source("DEVELOPMENT_function.R")
source("DEVELOPMENT_THRESHOLD_function.R")
source("LOAD_DATA_function.R")
par(mfrow=c(1,1))
SQG_SDM <- load_base("Petaurus_norfolcensis_pre1750_SDM_GH.tif")
ConLayer <- load_weight("ConditionLayer_masked.asc") # Our weightRaster which is here a condition layer reflecting current landuse in our landscape
#Dev <-load_devs("SQG_LRGDEVS/SQG_DEVELOPMENT_", n = 50)

for (i in 1:10) {raster_dev <- Developing(baseRaster = SQG_SDM, 
                                          weightRaster = ConLayer, 
                                          devArea=10000, 
                                          repeats=10,
                                          bufZone=30000,
                                          maxvalue = TRUE)
writeRaster(raster_dev, paste0("SQG_LRGEDEVS_4/SQG_DEVELOPMENT_MAX_",i,".asc"), sep="_",format = "ascii",overwrite=TRUE)
}

Test <- raster("SQG_LRGEDEVS_3/SQG_DEVELOPMENT_TSH_3.asc")
NAs1 <- which(is.na(values(SQG_SDM)))
NAs2 <- which(is.na(values(Test)))
developedArea <- length(NAs2) - length(NAs1)
message(paste("The number of development cells is", developedArea))



##### Run function on each development file #######
for (i in 1:50) {raster_off <- offsetting(baseRaster = SQG_SDM, 
                                          weightRaster = ConLayer, 
                                          devRaster = Dev[[i]], 
                                          bufferZone = 40000, 
                                          offPatches = 10, 
                                          type = "area", 
                                          threshold = 1,
                                          offsetMultiplier= 2,
                                          incrementCell = 40)
writeRaster(raster_off, paste0("SQG_LRGEOFFS/SQG_AREAXSDM/SQG_LRG_AREAXSDM_OFFSET",i,".asc"), sep="_",format = "ascii",overwrite=TRUE)
}

for (i in 1:50){raster_off <- Vegetation_offsetting(baseRaster=SQG_SDM , 
                                                    weightRaster1=ConLayer, 
                                                    weightRaster2=vegelayer, 
                                                    devRaster=Dev[[i]], 
                                                    bufferZone=40000, 
                                                    type="value", 
                                                    threshold=0.5,
                                                    offsetMultiplier=1, 
                                                    incrementCell=40, 
                                                    offPatches = 10)
writeRaster(raster_off, paste0("SQG_LRGEOFFS/SQG_VEGECONXSDM/SQG_LRG_VEGCONXSDM_OFFSET", i ,".asc"), sep="_", format="ascii", overwrite=TRUE)
}
