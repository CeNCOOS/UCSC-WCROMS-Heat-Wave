# UCSC-WCROMS-Heat-Wave
This is code to use the hobday heatwave criteria as implemented by Eric. C. J. Oliver using a 31 year ROMS climatology with WCOFS model output.
The idea is to use the WCOFS ocean model to create a 3-day forecast of heatwave conditions off the US West Coast.  This particular code uses the SST data.

# File and descriptions:
``` USMaritimeLimitsNBoundaries.shp ``` A shape file used to truncate model domain data to the US West coast

``` compute_30_yr_heat_content_heatwave.py ``` The file to generate yearly output from the 31 year UCSC ROMS reanalysis product.  This code outputs the upper ocean 100m heat content

``` compute_30_yr_heat_content_heatwave_new.py ``` Do the same as the previous file but only use the last 4 days of each update for the model run







### References
Hobday, A.J. et al. (2016), A hierarchical approach to defining marine heatwaves, Progress in Oceanography, 141, pp. 227-238, doi: 10.1016/j.pocean.2015.12.014 pdf

### Heatwave computation code:
https://github.com/ecjoliver/marineHeatWaves

