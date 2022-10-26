# UCSC-WCROMS-Heat-Wave
This is code to use the hobday heatwave criteria as implemented by Eric. C. J. Oliver using a 31 year ROMS climatology with WCOFS model output.
The idea is to use the WCOFS ocean model to create a 3-day forecast of heatwave conditions off the US West Coast.  This particular code uses the SST data.

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







### References
Hobday, A.J. et al. (2016), A hierarchical approach to defining marine heatwaves, Progress in Oceanography, 141, pp. 227-238, doi: 10.1016/j.pocean.2015.12.014 pdf

### Heatwave computation code:
https://github.com/ecjoliver/marineHeatWaves

