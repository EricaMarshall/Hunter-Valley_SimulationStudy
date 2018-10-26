# ----- Load data -----# 

# Load base raster layer 

load_base <- function(baseRaster){
  require(raster)
  mydata <- raster(baseRaster)
  mydata <- mydata/1000
  plot(mydata)
  return(mydata)
}

# Load weight raster layer

load_weight <- function(weightRaster){
  require(raster)
  weight_layer <- raster(weightRaster)
  plot(weight_layer)
  return(weight_layer)
}


# load list of development rasters 
load_devs <- function(devRaster, n){
  require(raster)
  raster_list <- list()
  for (i in c(1:n)) {
    raster_list[[i]] <- raster(print(paste0(devRaster,i,".asc", sep = "")))
  }
  Dev <- stack(raster_list)
  return(Dev)
}


# load list of offset rasters 

load_offs <- function(offRaster, n){
  require(raster)
  raster_list_off <- list()
  for (i in c(1:n)) {
    raster_list_off[[i]] <- raster(print(paste0(offRaster,i,".asc", sep = "")))
  }
  off <- stack(raster_list_off)
  return(off)
}


