import sys
import numpy as np
import xarray as xr
import scipy.io
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy import integrate
import datetime as dt
import time, scipy.io
# Weird Bug here, but make sure to run this first
import salem
from shapely import geometry
sanctuary_outline= salem.read_shapefile('/home/flbahr/heat_content/shapefiles/USMaritimeLimitsAndBoundariesSHP/USMaritimeLimitsNBoundaries.shp')
sanctuary_outline.crs = 'epsg:4326'
eez_shape = sanctuary_outline[(sanctuary_outline['REGION'] == "Pacific Coast") & (sanctuary_outline['EEZ'])]
bbox = ((eez_shape['max_y'].max(),eez_shape['min_y'].min()),(eez_shape['max_x'].max(),eez_shape['min_x'].min()))
geom = geometry.box(minx=bbox[1][1],maxx=bbox[1][0],miny=bbox[0][1],maxy=bbox[0][0])
#
currtime=dt.datetime.now()
#
theyear=2021
#if currtime.month==1:
#    theyear=currtime.year-1
#else:
#    theyear=currtime.year
url='/home/flbahr/heat_content/WCOFS_SST_'+str(theyear)+'.nc'
wcofs=xr.open_dataset(url,decode_cf=True)
wcofs_temp=wcofs['wcofs_temperature']
wcofs_time=wcofs['time']
wcofs_latitude=wcofs['lat']
wcofs_longitude=wcofs['lon']
westcoast_temp = wcofs_temp.sel(lat=slice(bbox[0][1], bbox[0][0]), lon=slice(bbox[1][1], bbox[1][0]))
longitude=westcoast_temp.lon.values
latitude=westcoast_temp.lat.values
data=wcofs_temp
time=wcofs_time
dims=['time','lat','lon']
#ds=xr.Dataset({'temperature_2023_surface':(dims,westcoast_temp)},
#              coords={'time':time,
#                      'lon':longitude,
#                      'lat':latitude,
#                      })
#ds.attrs['title'] = "West Coast Upper SST 2023 - 0 meters"
ds=xr.Dataset({'temperature_'+str(theyear)+'_surface':(dims,westcoast_temp)},
              coords={'time':time,
                      'lon':longitude,
                      'lat':latitude,
                      })
ds.attrs['title'] = "West Coast Upper SST "+str(theyear)+" - 0 meters"
ds.attrs['notes'] = "Created on "+dt.datetime.today().strftime("%Y-%m-%d") + " by flbahr"
fname = "/home/flbahr/heat_content/West_Coast_Temperature_"+str(theyear)+"_WCOFS.nc"
ds.to_netcdf(path=fname, mode='w')
