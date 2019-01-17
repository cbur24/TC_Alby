# TC_Alby

The Advanced Weather Research and Forecasting (WRF) model was used to simulate Tropical Cyclone Alby. This repo contains the scripts, namelists, and Vtables used to run WRF and post-process the results.  

*namelists_etc* : conatins namelists, Vtables, and shell scripts used for running WPS and WRF.

reanalysis_validation : contains scripts used to validate JRA55 as a sensible product for initializing and providing boundary conditions                           to the WRF model.

CreateTracks.ncl : generates a plot that compares IBTrACS cyclone paths against the simulated cyclone path.

plot_wrf_results.ipynb :  Jupyter Notebook that takes as an input a wrfout_do1* file and produces plots demostrating the evolution of the                           simulated cyclone. 

track_comparisons.ipynb : Jupyter Notebook that ompares IBTrACS cyclone paths against the simulated cyclone path. An alternative to using                           the CreateTracks.ncl script
