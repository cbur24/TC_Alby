# TC_Alby

The Advanced Weather Research and Forecasting (WRF) model was used to simulate Tropical Cyclone Alby. This repo contains the scripts, namelists, and Vtables used to run WRF and post-process the results.  

_namelists_etc_ : contains namelists, Vtables, and shell scripts for running WPS and WRF.

_reanalysis_validation_ : contains scripts for validating JRA55 as a sensible product for initializing and providing boundary conditions                           to the WRF model.

_CreateTracks.ncl_ : generates a plot that compares IBTrACS cyclone paths against the simulated cyclone path.

_plot_wrf_results.ipynb_ :  Jupyter Notebook that takes as an input a wrfout_do1* file and produces plots demonstrating the evolution of the                           simulated cyclone. 

_track_comparisons.ipynb_ : Jupyter Notebook that compares IBTrACS cyclone paths against the simulated cyclone path. An alternative to using                           the CreateTracks.ncl script
