
library(raster)
library(reshape)
Vegetation_offsetting <- function(baseRaster , weightRaster1, weightRaster2, devRaster, bufferZone=15000, type="value", threshold=0.5,
                                  offsetMultiplier=1, incrementCell=3, offPatches = 10){
  require(dplyr)
  require(sf)
  require(dismo)
  require(fasterize)
  require(data.table)
  if(is.null(weightRaster1)){
    newRaster <- baseRaster
  } else{
    newRaster <- baseRaster * weightRaster1
  }
  NAs1 <- which(is.na(values(newRaster))) #finding the NAs for newraster 
  NAs2 <- which(is.na(values(devRaster)))
  developedArea <- length(NAs2) - length(NAs1)
  message(paste("The number of development cells is", developedArea))
  values(devRaster)[NAs1] <- 0.5
  NAs <- which(!is.na(values(devRaster)))
  devRaster <- newRaster
  values(devRaster)[NAs] <- NA
  devRaster <- reclassify(devRaster, c(-Inf, Inf, 1))
  dev <- rasterToPolygons(devRaster, dissolve = TRUE)
  devBuff <- buffer(dev, bufferZone) %>% 
    crop(newRaster)
  values(newRaster)[which(!is.na(values(devRaster)))] <- NA
  par(mfrow=c(1,2))
  plot(newRaster)
  plot(dev, add=TRUE)
  plot(devBuff, add=T, border="red")
  plot(newRaster)
  plot(devBuff, add=T, border="red")
  devBuff <- devBuff - dev
  buffRaster <- fasterize(sf::st_as_sf(devBuff), devRaster)
  Alldevs <- seq_len(nlayers(devRaster))
  if(type == "area"){
    vals <- hist(weightRaster2 - devRaster)
    VegValue <- data.frame(class = 1:9,
                           NumCells = vals$counts)
    VegValue <- VegValue[VegValue$NumCells > 0,]
    offsetvalueTotal <- sum(VegValue[,2]) * offsetMultiplier ### Multiply the number of cells by the offset multiplier so when used in type == area twice the amount of area is required
    conimpact <- offsetvalueTotal / offPatches
  } else if(type == "value"){
    for(t in Alldevs){
      offsetvalueTotal <- data.frame(class=1:9, offsetValue=NA)
      if(is.null(weightRaster2)){
        conimpact <- sum(values(newRaster)[which(!is.na(values(devRaster[[t]])))])
        offsetvalueTotal <- conimpact * offsetMultiplier 
      } else { for (vg in 1:9) { 
        temp <- weightRaster2 == vg
        offsetvalueTotal[vg,2] <- sum(values(temp*baseRaster)[which(!is.na(values(devRaster[[t]])))],na.rm = T)
      } 
        offsetvalueTotal <- offsetvalueTotal[offsetvalueTotal$offsetValue > 0,]
        offsetvalueTotal[, 2] <- offsetvalueTotal[, 2] * offsetMultiplier
        #### Calculate the number of developed cells belowing to each vegetation class
        vals <- hist(weightRaster2-devRaster)
        
        VegValue <- data.frame(class = 1:9,
                               NumCells = vals$counts)
        VegValue <- VegValue[VegValue$NumCells > 0,]
        VegValue[,2] <- VegValue[,2] * offsetMultiplier ### Multiply the number of cells by the offset multiplier so when used in type == area twice the amount of area is required
        OverallOffset <- merge((data.table(offsetvalueTotal, key = "class")),
                               (data.table(VegValue, key = "class")))
        offsetConValue <- OverallOffset[,1:2]
        offsetvegvalue <- OverallOffset[,c(1,3)]
        conimpact <- sum(offsetConValue[,2])
        offsetvalueTotal <- OverallOffset ### So that further down the first column for both area and type can be used in the for loop
      }
      message("Pre-processing is done!")    
    }
    oo <- c()
    gg <- c()
    fullRaster <- newRaster
    for (p in 1:offPatches) {
      of_indices <- c()
      highPoints <- c()
      offsetGain <- c() # to count the offset gain we need to add to the raster later
      e <- 0
      s <- 0
      n <- 0
      rfp <- NA
      while(is.na(rfp)){
        nn <- dismo::randomPoints(buffRaster, 1)
        rfp <- values(newRaster)[cellFromXY(newRaster, nn)]
      }
      # calculate minimum buffer
      if(type == "area"){
        resRaster <- res(newRaster)[1] ^ 2
        rr <- sqrt((resRaster * offsetvalueTotal / 3.141593))
      } else if(type == "value"){
        rr <- sum(offsetConValue[,2]) / 2
      }
      points(nn)
      for(class in offsetvalueTotal[,1]){
        Vegecover <- c()
        c <- 0
        # while (c < offsetvegvalue[,2]) {
        #   b <- buffer(SpatialPoints(nn), rr + e)
        #   e <- e + (res(newRaster)[1] * incrementCell)
        #   v <- unlist(cellFromPolygon(newRaster, b))
        #   v <- v[which(!is.na(values(newRaster)[v]))]
        # }
        while(s < conimpact){
          b <- buffer(SpatialPoints(nn), rr + e)
          # plot(b, add=TRUE)
          e <- e + (res(newRaster)[1] * incrementCell)
          v <- unlist(cellFromPolygon(newRaster, b))
          v <- v[which(!is.na(values(newRaster)[v]))]
          for(i in v){
            if(values(newRaster)[i] <= threshold){ # value of suitability or combined raster
              n <- n + 1
              of_indices[n] <- i
              if(is.null(weightRaster1)){
                offsetGain[n] <- threshold
                if(type == "value"){
                  s <- s + (threshold - values(newRaster)[i]) # condition value
                } else if(type == "area"){
                  s <- n
                }
              } else{
                if(values(baseRaster)[i] < threshold){ # values of the base raster
                  offsetGain[n] <- values(baseRaster)[i]
                  if(type == "value"){
                    s <- s + (values(baseRaster)[i] - values(newRaster)[i])
                  } else if(type == "area"){
                    s <- n
                  }
                } else{
                  offsetGain[n] <- threshold
                  if(type == "value"){
                    s <- s + (threshold - values(newRaster)[i]) # condition value
                  } else if(type == "area"){
                    s <- n
                  }
                }
              }
            } else{
              highPoints <- append(highPoints, i)
            }
            if(s >= conimpact){
              break
              if(c >= offsetvalueTotal[,1]){
                break
              }
            }
          }
          searchedPoints <- unique(append(of_indices, highPoints))
          points(xyFromCell(newRaster, of_indices), col="blue", cex=0.1)
          values(newRaster)[searchedPoints] <- NA
          values(baseRaster)[searchedPoints] <- NA
          values(buffRaster)[searchedPoints] <- NA
        }
      }
      # plot(devRaster)
      # points(xyFromCell(newRaster, of_indices), col="blue", cex=0.3)
      message("Biodiversity offsetting recovered the habitat quality! :)")
      if(length(of_indices) > 10 || !is.null(of_indices)){
        gg <- append(gg, offsetGain)
        oo <- append(oo, of_indices)
      } else{
        message(paste("The offset of repeat", p, "has not been taken"))
        message(paste("The offset value is", offsetvalue))
      }
      
    }
    
    message(paste("Repeat", p, "is done!", "Hooray! :))"))
    
  }
  fullRaster[oo] <- gg
  return(fullRaster)
}

