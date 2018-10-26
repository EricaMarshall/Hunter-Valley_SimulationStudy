biOffset <- function(baseRaster, weightRaster=NULL, devArea=100, threshold=0.5, type="value", repeats=1, 
                     offsetMultiplier=2, bufZone=1500, maxvalue=FALSE){
  require(dplyr)
  require(dismo)
  require(spdep)
  if(is.null(weightRaster)){
    newRaster <- baseRaster
  } else{
    newRaster <- baseRaster * weightRaster
  }
  fullRaster <- newRaster
  message("Pre-processing is done!")
  dd2 <- c()
  oo2 <- c()
  gg <- c()
  for(t in 1:repeats){
    npixle <- 0
    while(npixle < devArea * 3){ # find buffers  with more than the minimum number of cells needed
      oldw <- getOption("warn")
      options(warn = -1)
      firstrand <- dismo::randomPoints(newRaster, 1, prob = TRUE)
      buf <- raster::buffer(SpatialPoints(firstrand), bufZone) # buffer size!
      samllRaster <- newRaster %>% 
        crop(buf) %>% 
        mask(buf)
      point2 <- rasterToPoints(samllRaster, spatial=TRUE)
      names(point2)[1] <- "suitability"
      if(!is.null(weightRaster)){
        point2$base <- raster::extract(baseRaster, point2)
      }
      npixle <- nrow(point2)
      options(warn = oldw)
    }
    par(mfrow=c(1,2))
    plot(newRaster, zlim=c(0,1), legend.width=1.2)
    plot(buf, add=T)
    plot(samllRaster, zlim=c(0,1), legend.width=1.2)
    neighbour <- spdep::knearneigh(point2, k=8)
    mp <- as.data.frame(firstrand)
    indices <- c()
    indices[1] <- which(neighbour$x[,1] == mp[1,1] & neighbour$x[,2] == mp[1,2])
    for(j in 2:devArea){
      nn <- neighbour$nn[indices,]
      nn <- nn[!nn %in% indices] # remove the previous indices from the neighbours
      if(maxvalue==TRUE){
        npoints <- point2[nn, ]
        npmax <- which(npoints@data[,"suitability"] == max(npoints$suitability))
        mp3 <- as.data.frame(npoints[npmax,]@coords)
        randomnn <- which(neighbour$x[,1] == mp3[1,1] & neighbour$x[,2] == mp3[1,2])
      } else{
        if(length(nn) > 1){ # when lenght(nn) < 1 it consider it 1:nn
          randomnn <- sample(nn, size=1, prob=point2$suitability[nn])
          indices[j] <- randomnn
          thePoint <- neighbour$x[randomnn,]
          points(thePoint[1], thePoint[2], pch=16, cex=0.4)
        } else if(length(nn) == 1){
          randomnn <- nn
          indices[j] <- randomnn
          thePoint <- neighbour$x[randomnn,]
          points(thePoint[1], thePoint[2], pch=16, cex=0.4)
        } else{
          message("Very strangly end up finding a NULL neighbourhood!")
        }
      }
    }
    message("Some development happened! :/")
    newPoint <- point2[-indices,]
    newNeighbour <-  spdep::knearneigh(newPoint, k=8)
    of_indices <- sample(1:nrow(newPoint), 1)
    if(type == "value"){
      offsetvalue <- sum(point2$suitability[indices]) * offsetMultiplier # the valuse should be added for offset
    } else if(type == "area"){
      offsetvalue <- length(indices) * offsetMultiplier
    }
    s <- 0 # offseting gain in each point
    n <- 0
    highPoints <- c()
    offsetGain <- c() # to count the offset gain we need to add to the raster later
    while(s < offsetvalue){
      highPoints <- unique(append(highPoints, of_indices))
      nn <- newNeighbour$nn[highPoints,]
      nn <- nn[!nn %in% highPoints]
      if(length(nn) > 1){ # when lenght(nn) < 1 it consider it 1:nn
        randomnn <- sample(nn, 1)
      } else if(length(nn) == 1){
        randomnn <- nn
      } else{
        message("Couldn't find enough points for offsetting. Please increase the buffer size on other trials!")
        bufZone <- bufZone * 1.1
        break
      }
      if(newPoint$suitability[randomnn] < threshold){ # condtion on the value of condition raster
        n <- n + 1
        of_indices[n] <- randomnn
        # thePoint <- newNeighbour$x[randomnn,]
        # points(thePoint[1], thePoint[2], pch=16, cex=0.4, col="red")
        if(is.null(weightRaster)){
          offsetGain[n] <- threshold
          if(type == "value"){
            s <- s + (newPoint$suitability[randomnn] - threshold)
          } else if(type == "area"){
            s <- n
          }
        } else{
          if(newPoint$base[randomnn] < threshold){
            offsetGain[n] <- newPoint$base[randomnn]
            if(type == "value"){
              s <- s + (newPoint$base[randomnn] - newPoint$suitability[randomnn])
            } else if(type == "area"){
              s <- n
            }
          } else{
            offsetGain[n] <- threshold
            if(type == "value"){
              s <- s + (threshold - newPoint$suitability[randomnn])
            } else if(type == "area"){
              s <- n
            }
          }
        }
      } else{
        highPoints <- append(highPoints, randomnn)
      }
      if(length(of_indices) %% devArea == 0){
        thePoint <- newNeighbour$x[of_indices,]
        points(thePoint, pch=16, cex=0.4, col="red")
      }
      if(length(unique(highPoints)) == (nrow(newPoint) - 2)){
        message("Couldn't find enough points for offsetting. Please increase the buffer size on other trials!")
        bufZone <- bufZone * 1.1
        break
      }
    }
    thePoint <- newNeighbour$x[of_indices,]
    points(thePoint, pch=16, cex=0.4, col="red")
    message("Biodiversity offsetting recovered the habitat quality! :)")
    if(length(of_indices) > 10 || !is.null(of_indices)){
      dd <- c()
      oo <- c()
      dd <- cellFromXY(newRaster, neighbour$x[indices, ])
      oo <- cellFromXY(newRaster, newNeighbour$x[of_indices, ])
      gg <- append(gg, offsetGain)
      dd2 <- append(dd2, dd)
      oo2 <- append(oo2, oo)
      if(repeats > 1 && t < repeats){
        # point <- point[-append(dd, oo), ] # surprisingly this takes time in sf objects, I don't know why!
        raster::values(newRaster)[append(dd, oo)] <- NA
      }
    } else{
      message(paste("The offset of repeat", t, "has not been taken"))
      message(paste("The offset value is", offsetvalue))
    }
    message(paste("Repeat", t, "is done!", "Hooray! :))"))
    raster::writeRaster(newRaster, "modelSave.asc", overwrite=TRUE)
    write.table(dd2, "development.txt")
    write.table(oo2, "Offset.txt")
    write.table(gg, "Gain.txt")
  }
  ls <- list(development=dd2, offset=oo2, gain=gg)
  print(paste("The number of deveopment cells removed:", length(dd2)))
  print(paste("The number of offsetting cells added:", length(oo2)))
  print(paste("The total amount of development removed:", sum(values(fullRaster)[dd2], na.rm=T)))
  print(paste("The total amount of offset gain:", sum(gg)))
  # add --coverting to raster-- here
  outDevRaster <- fullRaster
  raster::values(outDevRaster)[ls[["development"]]] <- NA
  outOffRaster <- outDevRaster
  raster::values(outOffRaster)[ls[["offset"]]] <- ls[["gain"]]
  outDevRaster <- stack(outDevRaster, outOffRaster)
  names(outDevRaster) <- c("Development", "Offset")
  return(outDevRaster)
}

##### Load Target SDMs #######
SQG_SDM <- raster("Petaurus_norfolcensis_pre1750_SDM_GH.tif") # Our baseRaster
SQG_SDM <- SQG_SDM/1000 # To change scale to between 0 and 1 because this happens to have a maxent output of 0-1000
plot(SQG_SDM) #check 

ConLayer <- raster("GH.Condition.tif") # Our weightRaster which is here a condition layer reflecting current landuse in our landscape
plot(ConLayer) #Check

