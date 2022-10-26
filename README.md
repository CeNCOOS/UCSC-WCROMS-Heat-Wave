# UCSC-WCROMS-Heat-Wave
This is code to use the hobday heatwave criteria as implemented by Eric. C. J. Oliver using a 31 year ROMS climatology with WCOFS model output.
The idea is to use the WCOFS ocean model to create a 3-day forecast of heatwave conditions off the US West Coast.  This particular code uses the SST data.

Note data files have not been upload to github except the shapefile below.  The yearly files containing daily values can be created using several of the programs
below.

# File and descriptions:
``` USMaritimeLimitsNBoundaries.shp ``` A shape file used to truncate model domain data to the US West coast

``` compute_30_yr_heat_content_heatwave.py ``` The file to generate yearly output from the 31 year UCSC ROMS reanalysis product.  This code outputs the upper ocean 100m heat content

``` compute_30_yr_heat_content_heatwave_new.py ``` Do the same as the previous file but only use the last 4 days of each update for the model run

``` compute_30_yr_sst_new.py ``` Do the same as the file above but only extract the SST and don't compute heat content.

``` grab_sst_for_output.py ``` Grab all the SST data from the aggregated UCSC ROMS NRT data and output as a NetCDF file.  This is used by the heatwave computation.

``` heat-content.py ``` Patrick Daniel's upper 100m ocean heat content code for use with the UCSC ROMS NRT output.

``` heat_wave_plot_UCSC_ROMSsst.py ``` Code to read the output files from the heatwave computation and plot them as daily png files.  Forecast plots have a red title date.

``` heat_wave_run_with_climate_sst_WCOFS.py ``` Code to compute MHW from the UCSC ROMS 31 year reanalysis run, with the UCSC ROMS NRT model run and the WCOFS runs for 2021 to the present with the forecast from WCOFS.  The output is stored as python pickle files to seperate steps and prevent having to run the code again if there was
an error in the plotting routines.

``` heatwave_to_netcdf.ipynb ``` Python notebook to take heatwave output computed from upper ocean heat content pickle files and write a NetCDF file.

``` heatwave_to_netcdf_sst.ipynb ``` Python notebook to take heatwave computed from SST and write and NetCDF file.

``` load_wcofs_for_heatwave.py ``` Load the current WCOFS output up to the present and then also load the forecast data for use with the heatwave code.

``` marineHeatWave.py ``` Eric C. J. Oliver's heatwave computation code base on the Hobday definition.  See below for link to his code

``` merge_heat.ipynb ``` Python notebook exploring how to merge data from 31 year UCSC ROMS and UCSC ROMS NRT.  The data are the upper ocean heat content.
The data are used to compute the marine heatwave and then output into pickle files.

``` merge_sst.ipynb ``` Python notebook to do the same as above but using sst data instead of heat content data.

``` set_depth.py ``` Python code from ROMS for computing depths from sigma coordinates

``` setup_heatwave.py ``` Eric C.J. Oliver's code for setting up the marine heatwave code

``` stretching.py ``` Python code from ROMS using in the depth computation from sigma coordinates.

``` truncate_wcofs.py ``` Code to truncate WCOFS to the domain of interest, used only on archived data files.

``` wcof_heatcontent.py ``` Code only using WCOFS data to compute 100m upper ocean heat content

``` wcof_load.py ``` Initial code to load WCOFS data but doesn't take into account weird timing with WCOFS output files.

``` wcof_merge_data.py ``` Code to merge test output from python pickle files and create a NetCDF file.

``` wcof_sst_output.py ``` Code to test merging files and outputing not used by the above.

``` wcofs_compare.py ``` Code for comparing some of the WCOFS files from nowcast vs forecast.  Not used above.

``` wcofs_lonlat_2_xy.py ``` Code from Alexander Kuropov (author of WCOFS) to regrid local coordinats.  Ultimately wasn't used with the heatwave code above.

### References
Hobday, A.J. et al. (2016), A hierarchical approach to defining marine heatwaves, Progress in Oceanography, 141, pp. 227-238, doi: 10.1016/j.pocean.2015.12.014 pdf

### Heatwave computation code:
Eric C. J. Oliver's implementation of the Hobday heatwave criteria:

https://github.com/ecjoliver/marineHeatWaves

