
###############################################################
#generating times-slice spatial plots of the reanalysis data
#to check its ability to recreate TC alby.

###############################################################


library(ncdf4);library(RNetCDF);library(raster); library(maptools)
library(lubridate);library(stringr); library(rasterVis); library(viridis); library(RColorBrewer);
library(colorRamps); library(tidyverse); library(animation)

setwd("C:/Users/u18343/Desktop/rotation_3/reanalysis_validation")

# Can call CDO operators from R using system() e.g.:
#   file_grb2 = "001.grb2"
#   file_ncdf ="001.nc"
#   system(paste("cd ~/DATA/prate; cdo -f nc copy ",file_grb2,file_ncdf,sep=(" ")))

#----functions for code-------------------------------------------------------------------
splitLevels_JRA55 = function(fileloc, varname, timeVar_name = "initial_time0_hours", offset = 0,
                             convert_tz = F, varType = "pressure", export_ncdf = F,
                             return_dates = F){
  #JRA55 exports 3hr and 6hr forecasts as different levels with
  #the same timestamp. This function seperates the levels, adds the
  #correct timestamp, and then remerges the files to create a single
  #timeseries of rasterstacks for plotting/exporting etc.
  
  # fileloc = string, 
  # varname = string, 
  # timeVar_name = name of the time dimnesion in your netcdf file
  # offset = Numeric. Offset for time dimension, in seconds.
  # convert_tz = Boolean, if True converts UTC to Perth local time
  # varType = string, options are 'wind' or 'pressure' (converts units)
  # export_ncdf = Boolen, if true export the file as a netcdf
  # return_dates = Boolean, if true return only the sequence of dates (and not the rasterstack, runs quickly)
  
  #read in files as seperate levels
  d = fileloc
  b_lev1 = brick(d, varname = varname, level=1)
  crs(b_lev1) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  
  d = fileloc
  b_lev2 = brick(d, varname = varname, level=2)
  crs(b_lev2) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  
  #get the time dimension from the netcdf file
  data.nc<- nc_open(fileloc)
  Zdim = ncvar_get(data.nc,varid = timeVar_name)
  
  #get the time units from netcdf
  tunits = ncatt_get(data.nc,timeVar_name, attname="units")
  tustr<-strsplit(tunits$value, " ")
  origin=paste(unlist(tustr)[3], unlist(tustr)[4], sep=" ")
  
  #turn numeric time values into R date object
  Zdim = as.POSIXct((Zdim*3600), tz='UTC', origin=origin) + offset 
  #create date objects for each level
  dates_lev1 = Zdim - 10800 #10800 is num of secs in 3hrs
  dates_lev2 = Zdim
  #bind these two dates to create an object with all the dates
  alldates= rbind(unname(dates_lev1), unname(dates_lev2))
  unix2POSIXct <- function (time){
    structure(time, class = c("POSIXt","POSIXct"))
  } 
  alldates = unix2POSIXct(alldates)
  alldates = with_tz(alldates, tzone="UTC")
  
  if (convert_tz == T) {
    #convert dates to Perth local time "AWST"
    dates_lev1 = with_tz(dates_lev1, tzone="Australia/Perth")
    dates_lev2 = with_tz(dates_lev2, tzone="Australia/Perth")
    alldates = with_tz(alldates, tzone = "Australia/Perth")
  }
  if (return_dates == F) {
    #Label the rasters by their dates
    b_lev1 = setNames(b_lev1, dates_lev1)
    b_lev2 = setNames(b_lev2, dates_lev2)
    
    #merge the stacks, then sort them by date, then set the Z variable
    b = stack(b_lev1, b_lev2) 
    b = subset(b, sort(names(b)))
    
    #convert units
    if (varType == "pressure"){
      b = b/100
      b = setZ(b, alldates)
      #export netcdf
      if (export_ncdf == T){
        writeRaster(b, filename="results/prmsl_3hr_flattened", format="CDF", 
                    varname = "prmsl", varunit = "hPa", longname = "Pressure reduced to mean sea level",
                    zname = "time",  zunit=tunits, overwrite=TRUE, force_v4=TRUE, compression=7)
      }
    }  else if (varType == "wind"){
      b = setZ(b, alldates)
      #export netcdf
      if (export_ncdf == T){
        writeRaster(b, filename="results/wsmax_surf_3hr_flattened", format="CDF", 
                    varname = "ws", varunit = "m/s", longname = "3hr maximum wind speed",
                    zname = "time", zunit=tunits, overwrite=TRUE, force_v4=TRUE, compression=7)
      } 
    }
    return(b)
  } else if (return_dates == T){
    return(alldates)
  }    
}

rasterStat_for_plot = function(rasterstack, stat = "minimum"){
  #find statistic for individual raster layers for adding to plot title
  
  #crop the stack so we can exclude the jetstream lows 
  tc_crop <- as(raster::extent(85, 140.0, -40 , 0.0), "SpatialPolygons")
  proj4string(tc_crop) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  stack_trim =crop(rasterstack, tc_crop)
  
  #find stat of TC
  var=list()
  for (layer in 1:nlayers(stack_trim)) {
    if (stat == "minimum"){
      tmp <- raster::minValue(stack_trim[[layer]])
      tmp = round(tmp, digits=0)
    } else if (stat == "maximum"){
      tmp <- raster::maxValue(stack_trim[[layer]])
      tmp = round(tmp, digits=1) 
    }  
    var[[layer]] = tmp
  }
  return(unlist(var))
}

##################################################################################################
#bring in data

#Windspeed. Using 3hr max wind speed from derived surface variables.
jra55_windMax = splitLevels_JRA55(fileloc="data/renanalysis/surface/wsmax_surf_3hr_originalFile.nc",
                                  varname = "WSMX_GDS4_HTGL", offset = 21600, varType= "wind")

#sea level pressure
jra55_prmsl = splitLevels_JRA55(fileloc="data/renanalysis/surface/fcst_surf_prmsl_3hr_originalFile.nc",
                                varname = "PRMSL_GDS4_MSL",offset = 21600, varType= "pressure")

#creating bounding box
  extent_aus <- as(raster::extent(85, 140.0, -62 , 0.0), "SpatialPolygons")
  proj4string(extent_aus) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  
#bring in world continents shapefile
  continents = readShapeSpatial("data/spatial/world_continents.shp")
  crs(continents) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

#bring in world TC track shapefile
  tcAlby_path = readShapeSpatial("data/spatial/ibtracks/tc_alby_path.shp")
  crs(tcAlby_path) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  #tcAlby_path = as(tcAlby_path, "SpatialPointsDataFrame")

#crop
  jra55_windMax =crop(jra55_windMax, extent_aus)
  jra55_prmsl =crop(jra55_prmsl, extent_aus)
  continents = crop(continents, extent_aus)
  tcAlby_path = crop(tcAlby_path, extent_aus)

#----creating stacks and labelling layers for plotting-------------------------------------------------  

#creating stacks of 6 time slices corresponding to timing of satellite images.
  
#surface pressure-----
  jra55_sp_stack = stack(jra55_prmsl[['X1978.03.30.00.00.00']], jra55_prmsl[['X1978.04.02.00.00.00']], 
                       jra55_prmsl[['X1978.04.03.00.00.00']], jra55_prmsl[['X1978.04.04.00.00.00']],
                       jra55_prmsl[['X1978.04.04.12.00.00']], jra55_prmsl[['X1978.04.04.18.00.00']])  
  #----------
  #finding the lowest pressure of the TC for each timeslice in the plot
  minpress= rasterStat_for_plot(jra55_sp_stack)

  #add the minimum pressure of TC to the layer title
  sp_stackNames = c(paste("30/3 00:00,", minpress[1], "hPa"), paste("02/4 00:00,", minpress[2], "hPa"),
                    paste("03/4 00:00,", minpress[3], "hPa"), paste("04/4 00:00,", minpress[4], "hPa"),
                    paste("04/4 12:00,", minpress[5], "hPa"), paste("04/4 18:00,", minpress[6], "hPa"))
  
  jra55_sp_stack = setNames(jra55_sp_stack, sp_stackNames)

#windspeed-----
  jra55_windMax_stack =  stack(jra55_windMax[['X1978.03.30.00.00.00']], jra55_windMax[['X1978.04.02.00.00.00']], 
                                jra55_windMax[['X1978.04.03.00.00.00']], jra55_windMax[['X1978.04.04.00.00.00']],
                                jra55_windMax[['X1978.04.04.12.00.00']], jra55_windMax[['X1978.04.04.18.00.00']])
  
  #----------
  #finding the highest windspeed of the TC for each timeslice in the plot
  maxwind = rasterStat_for_plot(jra55_windMax_stack, stat="maximum")
  
  #add the max windspeed to the layer title
  ws_stackNames = c(paste("30/3 00:00,", maxwind[1], "m/s"), paste("02/4 00:00,", maxwind[2], "m/s"),
                    paste("03/4 00:00,", maxwind[3], "m/s"), paste("04/4 00:00,", maxwind[4], "m/s"),
                    paste("04/4 12:00,", maxwind[5], "m/s"), paste("04/4 18:00,", maxwind[6], "m/s"))
  
  jra55_windMax_stack = setNames(jra55_windMax_stack, ws_stackNames)


#-------plotting---------------------------------------------------------------------------------------------

  p.strip <- list(cex=1.5, lines=1) 

  #--surface pressure--
  levelplot(jra55_sp_stack, col.regions = colorRampPalette(brewer.pal(8, 'RdBu')), contour=T, lwd=0.5,
          at=c(-Inf, seq(950, 1030, 10), Inf), names.attr= sp_stackNames, par.strip.text=p.strip,
          ylab = "", xlab="", pretty=TRUE) +
  latticeExtra::layer(sp.polygons(continents, col = "black", aplha=0.1, lwd = 2)) +
  spplot(tcAlby_path, zcol='Pcentre', lwd=3, at=c(-Inf, seq(950, 1030, 10), Inf), 
         col.regions = brewer.pal(8,'RdBu'), colorkey=F)
  

  #--wind 6hr---
  levelplot(jra55_windMax_stack, col.regions= colorRampPalette(c("white", "blue", "purple", "red")),
          at=c(seq(0, 35, 2.5), Inf),
          names.attr= ws_stackNames, par.strip.text=p.strip,
          ylab = "", xlab="", pretty=TRUE) +
  latticeExtra::layer(sp.polygons(continents, col = "black", lwd = 2))+  
  latticeExtra::layer(sp.polygons(tcAlby_path, col = "darkgreen",lwd = 3))


#----create animation to show cyclone path-------------------------------------------------------------------------

#clean out empty layers to prevent animation crashing, and clip to the time period we want
jra55_prmsl_animate = dropLayer(jra55_prmsl, c(1:90,138,180:nlayers(jra55_prmsl)))

animation::saveGIF({
  for(i in 1:nlayers(jra55_prmsl_animate)){
    l <- levelplot(jra55_prmsl_animate[[i]], col.regions = colorRampPalette(brewer.pal(8, 'RdBu')), contour=T, lwd=0.5,
                     main=substr(names(jra55_prmsl_animate)[i], 2,17),
                     names.attr= substr(names(jra55_prmsl_animate)[i], 2,17),
                     at=c(-Inf, seq(950, 1030, 10), Inf), ylab = "", xlab="", pretty=TRUE, margin=F) +
       latticeExtra::layer(sp.polygons(continents, col = "black", aplha=0.1, lwd = 2)) +
       spplot(tcAlby_path, zcol='Pcentre', lwd=3, at=c(-Inf, seq(950, 1030, 10), Inf), 
             col.regions = brewer.pal(8,'RdBu'), colorkey=F)
    plot(l)
  }
}, interval=0.15, movie.name="results/prmsl_animation.gif")



















