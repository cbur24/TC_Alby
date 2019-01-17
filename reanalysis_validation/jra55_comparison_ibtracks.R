


##########################################################################
# This script compares maxwind and lowest pressure data from IBtracks with
# JRA55 for TC ALBY
##########################################################################

library(ncdf4);library(RNetCDF);library(raster); library(maptools)
library(lubridate);library(stringr); library(rasterVis); library(viridis); library(RColorBrewer);
library(colorRamps); library(tidyverse)

setwd("C:/Users/u18343/Desktop/rotation_3/reanalysis_validation")

#--functions for code------------------------
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


#---bring in data----------------------------------------------------------------------------------


jra55_windMax = splitLevels_JRA55(fileloc="data/renanalysis/surface/wsmax_surf_3hr_originalFile.nc",
                                  varname = "WSMX_GDS4_HTGL", offset = 21600, varType= "wind")

#sea level pressure
jra55_prmsl = splitLevels_JRA55(fileloc="data/renanalysis/surface/fcst_surf_prmsl_3hr_originalFile.nc",
                                varname = "PRMSL_GDS4_MSL",offset = 21600, varType= "pressure")

dates_press = splitLevels_JRA55(fileloc="data/renanalysis/surface/fcst_surf_prmsl_3hr_originalFile.nc",
                                varname = "PRMSL_GDS4_MSL",offset = 21600, varType= "pressure",
                                return_dates =T)

dates_wind= splitLevels_JRA55(fileloc="data/renanalysis/surface/wsmax_surf_3hr_originalFile.nc",
                                  varname = "WSMX_GDS4_HTGL", offset = 21600,varType= "wind", 
                                  return_dates =T)

#import ibtracks data on tc_alby
tc_alby_data = read_csv('data/tc_alby.csv', col_types = cols(Year = col_character(),Month = col_character(),
                                                             Day = col_character(), Hour = col_character(),
                                                             Minute = col_character()))

#mutate to get dateTime column
tc_alby_data = tc_alby_data %>% 
  dplyr::mutate(date = paste(tc_alby_data$Year, tc_alby_data$Month, tc_alby_data$Day, sep="-")) %>%
  dplyr::mutate(time = paste(tc_alby_data$Hour, tc_alby_data$Minute, '00',sep=":")) %>%
  dplyr::mutate(dateTime = ymd_hms(paste(date, time, sep=" ")))


#---pressure---------------------------------------------------------------------------------

ibtracks_press = tc_alby_data %>%
                  dplyr::select(dateTime, 'pressure' = Pcentre)

#get vortex centre pressure from the jra55
jra55_vortex_press = rasterStat_for_plot(jra55_prmsl)
jra55_vortex_press = data_frame('dateTime' = as.character(dates_press), 'pressure' = jra55_vortex_press)
jra55_vortex_press$dateTime = ymd_hms(jra55_vortex_press$dateTime) #handling datetimes sucks
jra55_vortex_press = jra55_vortex_press %>%
                      filter(dateTime >= as.Date("1978-03-25") & dateTime <= as.Date("1978-04-06"))

#---wind-------------------------------------------------------------------------------------

ibtracks_wind = tc_alby_data %>%
  dplyr::select(dateTime, 'wind' = MaxWind)

#get vortex centre pressure from the jra55
jra55_maxwind = rasterStat_for_plot(jra55_windMax, stat="maximum")
jra55_maxwind = data_frame('dateTime' = as.character(dates_wind), 'wind' = jra55_maxwind)
jra55_maxwind$dateTime = ymd_hms(jra55_maxwind$dateTime)
jra55_maxwind = jra55_maxwind %>%
  filter(dateTime >= as.Date("1978-03-25") & dateTime <= as.Date("1978-04-06"))


#---plotting----------------------------------------------------------------------------------

ggplot(NULL, aes(y = pressure, x = dateTime)) + 
  geom_point(data = na.omit(jra55_vortex_press), size=1.5, aes(colour= "JRA55")) +  
  geom_point(data = na.omit(ibtracks_press), size=1.5, aes(colour="IBtracks")) +
  geom_line(data = na.omit(jra55_vortex_press), size=1, aes(colour="JRA55")) +  
  geom_line(data = na.omit(ibtracks_press), size=0.1, aes(colour="IBtracks")) +
  xlab("Date") + ylab("Pressure (hPa)")+
  scale_colour_manual("Model",values=c("red","darkgreen"))

ggplot(NULL, aes(y = wind, x = dateTime)) + 
  geom_point(data = na.omit(jra55_maxwind), size=1.5, aes(colour= "JRA55")) +  
  geom_point(data = na.omit(ibtracks_wind), size=1.5, aes(colour="IBtracks")) +
  geom_line(data = na.omit(jra55_maxwind), size=1, aes(colour="JRA55")) +  
  geom_line(data = na.omit(ibtracks_wind), size=0.1, aes(colour="IBtracks")) +
  xlab("Date") + ylab("Max. Wind Speed (m/s)")+
  scale_colour_manual("Model",values=c("red","darkgreen"))



































