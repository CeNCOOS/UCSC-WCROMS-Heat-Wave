import sys
# sys.path.append("..")
from set_depth import set_depth
import numpy as np
import xarray as xr
# from shapely.geometry import Polygon, Point
import scipy.io
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy import integrate
import pandas as pd
import datetime as dt
import time, cmocean, scipy.io
from tqdm import tqdm
import seaborn as sns

# Weird Bug here, but make sure to run this first
import salem
from shapely import geometry
sanctuary_outline= salem.read_shapefile('./shapefiles/USMaritimeLimitsAndBoundariesSHP/USMaritimeLimitsNBoundaries.shp')
sanctuary_outline.crs = 'epsg:4326'

eez_shape = sanctuary_outline[(sanctuary_outline['REGION'] == "Pacific Coast") & (sanctuary_outline['EEZ'])]
bbox = ((eez_shape['max_y'].max(),eez_shape['min_y'].min()),(eez_shape['max_x'].max(),eez_shape['min_x'].min()))
geom = geometry.box(minx=bbox[1][1],maxx=bbox[1][0],miny=bbox[0][1],maxy=bbox[0][0])

url='WCOFS_SST_2022.nc'
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
ds=xr.Dataset({'temperature_2022_surface':(dims,westcoast_temp)},
              coords={'time':time,
                      'lon':longitude,
                      'lat':latitude,
                      })
ds.attrs['title'] = "West Coast Upper SST 2021 - 0 meters"
ds.attrs['notes'] = "Created on "+dt.datetime.today().strftime("%Y-%m-%d") + " by flbahr"
fname = "West_Coast_Temperature_2022_WCOFS.nc"
ds.to_netcdf(path=fname)
