#notes on compiling and running WRFV3.7.1

#IS THIS A BETTER OPTION FOR THE FUTURE? :  https://github.com/GIS4WRF/gis4wrf
#IT RUNS WPS/WRF FROM QGIS AND PLOTS RESULTS

#################################################################################
#  									Compiling
#################################################################################

#---module load-----
module load dot
module load intel-fc/12.1.9.293
module load intel-cc/12.1.9.293
module load openmpi/1.6.3
module load netcdf/4.3.3.1
export JASPERINC=/usr/include
export JASPERLIB=/usr/lib64
export WRFIO_NCD_LARGE_FILE_SUPPORT=1

#--clone model to your workspace------
cd /short/w85/cb3058/
git clone /projects/WRF/WRFV_3.7.1

#--Build WRF (ARW)--------
cd /short/w85/cb3058/WRFV_3.7.1/WRFV3/
./clean -a
./run_compile
#select intel compiler: choose option 3; dmpar (distrubuted memory)
#select nesting type: choose option 3; vortext following

##scrpit will then run, and it'll submit a job (WRF_compile) to the express queue.

#--Build WPS------------- 
cd /short/w85/cb3058/WRF/WRFV_3.7.1/WPS/
clean -a
./run_compile
#Select option 4 (dmpar_NOGRIB2) for the build methods.

#scrpit will then run, and it'll submit a job (WPS_compile) to the express queue.

#################################################################################
#    									WPS
#################################################################################

#--geogrid---

	#Use the ncl code util/plotgrids_new.ncl to plot up your domains before running
	#the geogrid program.

#link the namelist to WPS folder, and the right geogrid table to geogrid.tbl
ln -s /g/data/w85/cb3058/wrf/namelists/namelist_wps_tcAlby.txt namelist.wps
ln –s geogrid/GEOGRID.TBL.ARW GEOGRID.TBL
#then run geogrid
./geogrid.exe


#--Ungrib----

# link the Vtable, link the filenames with link_csh, run ungrib (change prefix to jra55_plev
ln -sf /g/data/w85/cb3058/wrf/initialization_files/Vtable.JMAGSM_plev Vtable
./link_grib.csh /g/data1a/w85/cb3058/wrf/initialization_files/jra55/grib/all/jra55*
./ungrib.exe >& ungrib.output_jra55

#NOT RUN. Rerun this process but for the fixed data (change prefix in namelist to 'jra55_fixed',
# and change dates in &share to 1981-01-01 00:00:00)
./link_grib.csh /g/data1a/w85/cb3058/wrf/initialization_files/jra55/grib/fixed/jra55*
./ungrib.exe >& ungrib.output_fixed

#Rerun this process but for the SST data (change prefix in namelist to 'sst', change back dates)
ln -sf /g/data/w85/cb3058/wrf/initialization_files/Vtable.E20C Vtable
./link_grib.csh /g/data1a/w85/cb3058/wrf/initialization_files/sst/sst*
./ungrib.exe >& ungrib.output_sst


	#--IF USING ECMWF-P
	#copy ecmwf_coeffs plaintext file to the working directory, this file was created
	#by using the coefficients in the JRA55 technical guide.
	#Link the .exe file to the working directory, relink JRA vtable, link coefficents table to dir	
	ln -sf util/calc_ecmwf_p.exe calc_ecmwf_p.exe
	ln -sf /g/data/w85/cb3058/wrf/initialization_files/Vtable.JMAGSM Vtable
	ln -sf /g/data1a/w85/cb3058/wrf/initialization_files/ecmwf_coeffs_original ecmwf_coeffs
	#In the working directory run:
	./calc_ecmwf_p.exe >& calc_ecmwf_p_output
	# if this runs succesfully, add 'PRES' to fg_name= before running metgrid.

#--Metgrid----
ln –s metgrid/METGRID.TBL.ARW METGRID.TBL
./metgrid.exe >& metgrid_output

#shift the newly created met_em files into a folder called 'met_em_files' in the
# WPS directory

#################################################################################
#    									WRF
#################################################################################

#now move to the WRF dir
cd /short/w85/cb3058/WRFV_3.7.1/WRFV3/test/em_real/
#link all of the em_met files to the folder (the "." at the end uses the file name)
ln -sf /short/w85/cb3058/WRFV_3.7.1/WPS/met_em_files/met_em.d* .

#link or copy your namelist.input file

	# IF BOGUSSING
	#edit the namelist record &tc to reflect the storm dimensions you wish to insert,
	#make sure to change the namelist to just the first datetime, and maxdom=1.
	#Need to set fine_input_stream = 0,2,2 so the nests use the met data from domain 1.
	#Need to set input_from_file = T,F,F otherwise get weird artefacts
	
	#First edit the namelist to REMOVE the storm.
	qsub run_tc.sh	#THIS NEEDS TO BE SERIALLY COMPILED
	#NOW RENAME OUTPUT FILE from "auxinput_do1 etc' to the met_em.d01<date> file that we just overwrote
	# now repeat above to insert the new storm.

#run real.exe to vertically interpolate the layers to the model layers
qsub run_real.sh 

#Once the output for real has been created (wrfbdy) and wrfinput) then we can run wrf
#the run_wrf.sh has mpirun before the "./wrf.exe" command to ensure it runs on multiple
#nodes. 
qsub run_wrf.sh 



