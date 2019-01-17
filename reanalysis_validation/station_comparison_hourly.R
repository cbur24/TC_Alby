
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#None of this data goes back to 1979
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

library(raster);library(readr); library(dplyr); library(lubridate); library(ncdf4); library(ggplot2)


#-------------------import netcdf and spatial files and process--------------
#JRA55
#6hr wind
d = "data/renanalysis/surface/windspeed_6hr.nc"
jra55_wind_6hr = stack(d, varname = "ws")
crs(jra55_wind_6hr) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"


#------import weather station files and clean----------------------------------

#station information header file
stationInfo = read_csv('data/weather_stations/headers/HM01X_StnDet_999999999425091.txt', col_names = F)
stationInfo = stationInfo[,-c(1,3,5,6,9:13,15:22)]
col_names = c("StationID", "StationName", "Lat", "Long", "year_started")
colnames(stationInfo) = col_names

#filter out stations north of Geraldton and east of Albany 
stationInfo = stationInfo %>% filter(Lat < -28.0) %>% 
                              filter(Long < 118.5) %>% 
                              filter(year_started < 1979)

#create a vector of filenames that correspond to those we filtered for earlier
filesToBind = character()  
for (station in stationInfo$StationID){
  tmp = paste("HM01X_Data_",station ,  "_999999999425091.txt",sep="")
  filesToBind[[station]] <- tmp
}
filesToBind = unname(filesToBind)

#import individual weather station data
setwd("data/weather_stations/hourly_WA/") 
station_wind <- do.call(rbind, lapply(filesToBind, read_csv))
setwd("C:/Users/Chad/Desktop/reanalysis_validation") 

