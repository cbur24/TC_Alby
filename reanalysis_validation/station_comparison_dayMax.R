

#######################################################
# This script produces point-to-pixel scatterplots 
# of daily maximum wind gusts from Jan 1978 to Dec 1978
#######################################################

#Script is functioning, but is not robust - requires automating some sections:
#       - Add user inputs section for adding filepaths etc
#       - Maybe add more statistics functions

#written by Chad Burton, 2/11/2018.
library(RNetCDF); library(maptools);library(rasterVis); library(viridis); 
library(RColorBrewer);library(colorRamps); library(raster);library(tidyverse);
library(lubridate); library(ncdf4); library(hydroGOF)

setwd("C:/Users/u18343/Desktop/rotation_3/reanalysis_validation")


#--functions for script-------

FUN.df_wrangle <- function(x){           
  #this function does the data wrangling to get the 
  #extraction df into the format I want
  x =  x %>% t() %>% data.frame() 
  colnames(x) = as.character(unique(station_wind_join$StationID))
  x = x[-1,]                      
  x = stack(x)                   
  no_of_stations= length(unique(station_wind_join$StationID))
  dates_col <- rep(as.POSIXlt(seq(ymd('1978-01-01'),          
                                  ymd('1978-12-31'), by = 'days')), times=no_of_stations) #"times=" the number of stations extracted
  dates_col = as.Date(dates_col)
  x = data.frame(dateTime=dates_col, x)
  #colnames(x) = c("dateTime", "ws", "StationID") #no idea why this isnt working
}

mae_ <- function(obs, model){
  x <- abs(obs-model)
  mean(x, na.rm=T)
}

#-------------------import netcdf--------------
#JRA55
#Maximum wind gusts per day, from 3hr wsmax file.
# This file was run through the 'spatial_plots' script (splitting levels), exported, then taken to cdo and had 
# settaxis run, follwed by daymax.)
d = "data/renanalysis/surface/wsmax_surf_3hrFlat_dayMax.nc"
jra55_1978_windMax = stack(d, varname = "ws")
crs(jra55_1978_windMax) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"


#----import weather station files, clean, and wrangle into order----------------------------------

#station information header file
stationInfo = read_csv('data/weather_stations/headers/DC02D_StnDet_999999999425050.txt')
stationInfo = stationInfo %>% dplyr::select('Bureau of Meteorology Station Number', 
                                     'Station Name',
                                     'Latitude to 4 decimal places in decimal degrees', 
                                     'Longitude to 4 decimal places in decimal degrees',
                                     'State',
                                     'First year of data supplied in data file',
                                     'Last year of data supplied in data file') %>% 
                              dplyr::rename('StationID' = 'Bureau of Meteorology Station Number', 
                                     'StationName' ='Station Name',
                                     'Lat'='Latitude to 4 decimal places in decimal degrees', 
                                     'Long' = 'Longitude to 4 decimal places in decimal degrees',
                                     'yearStart' = 'First year of data supplied in data file',
                                     'yearEnd' = 'Last year of data supplied in data file')


#filter out stations north of Geraldton, east of Albany and by dates
stationInfo = stationInfo %>% filter(State == "WA") %>% 
                              filter(Lat < -28.0) %>% 
                              filter(Long < 118.5) %>% 
                              filter(yearStart <= 1978) %>% 
                              filter(yearEnd > 1978)


#create a vector of filenames that correspond to those we filtered for.
filesToBind = character()  
for (station in stationInfo$StationID){
  tmp = paste("DC02D_Data_",station ,  "_999999999425050.txt",sep="")
  filesToBind[[station]] <- tmp
}
filesToBind = unname(filesToBind)

#import individual weather station data
setwd("data/weather_stations/daily_max_wind_gust/") 
station_wind <- do.call(rbind, lapply(filesToBind, read_csv))
setwd("C:/Users/u18343/Desktop/rotation_3/reanalysis_validation") 

#get rid of the '#' column
station_wind = station_wind[,-c(12)]

#create a dateTime column, remove unnecesary columns, filter dates, rename columns
station_wind = station_wind %>% 
               dplyr::mutate(dateTime = ymd(paste(station_wind$Year, station_wind$Month, station_wind$Day, sep="-"))) %>% 
               dplyr::select(-Year, -Day, -Month, -dc) %>% 
               dplyr::filter(dateTime >= as.Date("1978-01-01") & dateTime <= as.Date("1978-12-31")) %>% 
               dplyr::rename('ws' = `Speed of maximum wind gust in km/h`,
                             'StationID' = 'Station Number')

#calculate the number of non NA observations per station, if < 20 remove station
stations_toRemove = station_wind %>% group_by(StationID) %>% 
                                     summarise(obs = sum(!is.na(ws))) #count observations
stations_toRemove = as.character((stations_toRemove %>% filter(obs < 20))[1]) #select out station IDs to remove

#get stationIDs into same class first
stationInfo$StationID = as.integer(stationInfo$StationID) 
stationInfo$StationID = as.character(stationInfo$StationID)
station_wind$StationID = as.character(station_wind$StationID)

#now remove stations with low observations
station_wind = station_wind %>% filter(as.character(StationID) != stations_toRemove)
stationInfo = stationInfo %>% filter(as.character(StationID) != stations_toRemove)

#adjust windspeed to m/s
station_wind$ws = station_wind$ws/3.6              

#join the weather stations with their metadata
station_wind_join = left_join(stationInfo, station_wind, by = c("StationID"))


#----extract pixel values from JRA55 corresponding with station locations-----------------

#create spatial points
extractionPoints = stationInfo %>% dplyr::select(StationID, Lat, Long)
projection = crs(jra55_1978_windMax)
extractionPoints_spatial <- SpatialPointsDataFrame(extractionPoints[,3:2],
                                                   extractionPoints, 
                                                   proj4string = projection) 
#extract into dataframe
windMax_jra55_df = raster::extract(jra55_1978_windMax, extractionPoints_spatial,
                            method = 'simple', df = TRUE) 



windMax_jra55_df <- FUN.df_wrangle(windMax_jra55_df) #apply the function
colnames(windMax_jra55_df) = c("dateTime", "ws", "StationID") #change the colnames

#join the station dataframe with the model dataframe
jra55_windMax_Validation = left_join(station_wind_join, windMax_jra55_df, by = c("dateTime", "StationID"), 
                                 suffix=c("_station", "_model"))



#---calculate statistics--------------------
# Mean Absolute Error

stats_df = jra55_windMax_Validation %>% group_by(StationID) %>% 
                                        summarise(mae = mae_(ws_station,ws_model), #MAE
                                        d1 = hydroGOF::md (ws_model,ws_station))   #Index of agreement

#join stats to the validation dataframe
jra55_windMax_Validation <- left_join(jra55_windMax_Validation, stats_df, by="StationID")

#filter out just the dates around the TC event so we can plot those as well
jra55_windMax_Validation_marApril = jra55_windMax_Validation  %>% 
                                    dplyr::filter(dateTime >= as.Date("1978-03-15") & dateTime <= as.Date("1978-04-10"))


#----plotting-----------------------------
#plot all of 1978 data
ggplot(data=jra55_windMax_Validation, aes(x=ws_station, y=ws_model, colour=StationName)) + 
  facet_wrap(~StationName,scales="free_x",ncol=3)+
  geom_point(alpha=1, size=1.5, fill=NA, shape=21) +
  stat_smooth(method = "lm", linetype="longdash", col="black", size=0.6)+
  geom_abline(intercept = 0,slope=1)+
  theme_bw()+
  theme(aspect.ratio = 1)+
  xlab("Station wind speed (m/s)")+
  ylab("Model wind speed (m/s")+
  theme(legend.position="none")+
  theme(text = element_text(size=20))+
  theme(strip.background =element_rect(fill="peachpuff2")) +
  geom_text(aes(label=paste("D1 = ", round(d1,digits=2), sep = "")),col="grey30",
             x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, size=5)+
  geom_text(aes(label=paste("MAE = ", round(mae,digits=2), sep = "")),col="grey30",
          x=-Inf, y=Inf, hjust=-0.2, vjust=2.4, size=5)

#plot just the march april data
ggplot(data=jra55_windMax_Validation_marApril, aes(x=ws_station, y=ws_model, colour=StationName)) + 
  facet_wrap(~StationName,scales="free_x",ncol=3)+
  geom_point(alpha=1, size=1.5, fill=NA, shape=21) +
  stat_smooth(method = "lm", linetype="longdash", col="black", size=0.6)+
  geom_abline(intercept = 0,slope=1)+
  theme_bw()+
  theme(aspect.ratio = 1)+
  xlab("Station wind speed (m/s)")+
  ylab("Model wind speed (m/s")+
  theme(legend.position="none")+
  theme(strip.background =element_rect(fill="peachpuff2"))
#   geom_text(aes(label=paste("D1 = ", round(d1,digits=2), sep = "")),col="grey30",
#             x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, size=3.5)+
#   geom_text(aes(label=paste("MAE = ", round(mae,digits=2), sep = "")),col="grey30",
#             x=-Inf, y=Inf, hjust=-0.2, vjust=2.4, size=3.5)

