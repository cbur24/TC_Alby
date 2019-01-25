# TC_Alby

The Advanced Weather Research and Forecasting (WRF) model was used to simulate Tropical Cyclone Alby. This repo contains the scripts, namelists, and Vtables used to run WRF and post-process the results.  The jupyter notebooks contain instructions for their use at the top of the script

_namelists_etc_ : contains namelists, Vtables, and shell scripts for running WPS and WRF, and for initiating the .py scripts

_reanalysis_validation_ : contains scripts for validating JRA55 as a sensible product for initializing and providing boundary conditions                           to the WRF model.
SRC/
        _CreateTracks.ncl_ : generates a plot that compares IBTrACS cyclone paths against the simulated cyclone path.

        _plot_wrfout_d01.ipynb_ :  Jupyter Notebook that takes as an input a wrfout_do1* file and produces plots demonstrating the evolution of the                           simulated cyclone. 

        _plot_wrfout_static.py_ :  This is a python file that performs the same operations as the jupyter notebook _plot_wrfout_d01.ipynb_, only it can be run from the command line where it will prompt you for the user inputs. This is offered as an alternative to the notebook for instances in which the wrfout files are too large to be processed by a desktop computer; the .py file can be run on a HPC platform. Check out the 'run_plot_wrfout.py.sh' file in namelists_etc for intsructions on running this file

        _plot_wrfout_movingNests.ipynb_ :  Jupyter Notebook that takes as an input a nested domain (ie wrfout_d02, wrfout_d03 etc) file and produces plots demonstrating the evolution of the                           simulated cyclone.

        _plot_wrfout_movingNests.py_ :  This is a python file that performs the same operations as the jupyter notebook of the same name, only it can be run from the command line where it will prompt you for the user inputs. This is offered as an alternative to the notebook for instances in which the wrfout files are too large to be processed by a desktop computer; the .py file can be run on a HPC platform. Check out the 'run_plot_wrfout.py.sh' file in namelists_etc for intsructions on running this file

        _track_comparisons.ipynb_ : Jupyter Notebook that compares IBTrACS cyclone paths against the simulated cyclone path. An alternative to using                           the CreateTracks.ncl script
