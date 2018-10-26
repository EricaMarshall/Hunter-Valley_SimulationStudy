# incrementCell: the number of cells that increase around the previous buffer each time - the higher the faster
# but too much high might make a patch disconnected
offsetting <- function(baseRaster, weightRaster, devRaster, bufferZone=15000, offPatches=10, type="value", threshold=0.5,
                       offsetMultiplier=1, incrementCell=3){
  require(dplyr)
  require(sf)
  require(dismo)
  require(fasterize)
  if(is.null(weightRaster)){
    newRaster <- baseRaster
  } else{
    newRaster <- baseRaster * weightRaster
  }
  NAs1 <- which(is.na(values(newRaster)))
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
  devValue <- sum(values(newRaster)[which(!is.na(values(devRaster)))])
  values(newRaster)[which(!is.na(values(devRaster)))] <- NA
  values(baseRaster)[which(!is.na(values(devRaster)))] <- NA
  par(mfrow=c(1,2))
  plot(newRaster)
  plot(dev, add=TRUE)
  plot(devBuff, add=T, border="red")
  plot(newRaster)
  plot(devBuff, add=T, border="red")
  devBuff <- devBuff - dev
  buffRaster <- fasterize(sf::st_as_sf(devBuff), devRaster)
  # plot(buffRaster)
  if(type == "area"){
    offsetvalueTotal <- developedArea * offsetMultiplier
  } else if(type == "value"){
    offsetvalueTotal <- devValue * offsetMultiplier
  }
  offsetvalue <- offsetvalueTotal / offPatches
  message("Pre-processing is done!")
  oo <- c()
  gg <- c()
  fullRaster <- newRaster
  for(t in 1:offPatches){
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
    points(nn)
    # calculate minimum buffer
    if(type == "area"){
      resRaster <- res(newRaster)[1] ^ 2
      rr <- sqrt((resRaster * offsetvalue) / 3.141593)
    } else if(type == "value"){
      rr <- offsetvalue / 2
    }
    while(s < offsetvalue){
      b <- buffer(SpatialPoints(nn), rr + e)
      # plot(b, add=TRUE)
      e <- e + (res(newRaster)[1] * incrementCell)
      v <- unlist(cellFromPolygon(newRaster, b))
      v <- v[which(!is.na(values(newRaster)[v]))]
      for(i in v){
        if(values(newRaster)[i] <= threshold){ # value of suitability or combined raster
          n <- n + 1
          of_indices[n] <- i
          if(is.null(weightRaster)){
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
        if(s >= offsetvalue){
          break
        }
      }
      searchedPoints <- unique(append(of_indices, highPoints))
      points(xyFromCell(newRaster, of_indices), col="blue", cex=0.1)
      values(newRaster)[searchedPoints] <- NA
      values(baseRaster)[searchedPoints] <- NA
      values(buffRaster)[searchedPoints] <- NA
    }
    # plot(devRaster)
    # points(xyFromCell(newRaster, of_indices), col="blue", cex=0.3)
    message("Biodiversity offsetting recovered the habitat quality! :)")
    if(length(of_indices) > 10 || !is.null(of_indices)){
      gg <- append(gg, offsetGain)
      oo <- append(oo, of_indices)
    } else{
      message(paste("The offset of repeat", t, "has not been taken"))
      message(paste("The offset value is", offsetvalue))
    }
    message(paste("Repeat", t, "is done!", "Hooray! :))"))
  }
  print(paste("The number of cells used for offsetting:", length(oo)))
  print(paste("Sum of the development value previously taken:", devValue))
  print(paste("Sum of the offset gain added:", sum(gg) - sum(values(fullRaster)[oo])))
  fullRaster[oo] <- gg
  return(fullRaster)
}

