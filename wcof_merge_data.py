# file to merge data and output?
import xarray as xr
import numpy as np
from set_depth import set_depth
import scipy.io
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy import integrate
from scipy.interpolate import griddata
import pandas as pd
import datetime as dt
import time, cmocean, scipy.io
from tqdm import tqdm
import seaborn as sns
#from wcofs_lonlat_2_xy import wcofs_lonlat_2_xy
from wcof_load import wcof_load
import pickle
import pdb
from netCDF4 import Dataset, date2num, num2date
#
#
#
#mons=['march','april','may','june','july','august','september','october','november','december']
mons=['january','february','march','april','may','june','july','august','september','october']
imon=np.arange(0,len(mons))
abigtime=[]
abigtemp=[]
for i in imon:
    fname='wcofs_'+mons[i]+'_2022.p'
    thefile=open(fname,'rb')
    [time,temp]=pickle.load(thefile)
    #close(fname)
    if i==0:
        abigtime=time
        abigtemp=temp
        #pdb.set_trace()
    else:
        abigtime=np.append(abigtime,time)
        abigtemp=np.vstack((abigtemp,temp))
#
#ncout=Dataset('WCOFS_Temperature_2021.nc','w','NETCDF4')
#nc.createDimension

url='https://oceanmodeling.ucsc.edu:8443/thredds/dodsC/ccsra_2016a_phys_agg_slevs/fmrc/CCSRA_2016a_Phys_ROMS_Sigma-level_Aggregation_best.ncd'
roms_ds=xr.open_dataset(url,decode_cf=True)
lat=roms_ds['lat_rho'].values[:,0]
lon=roms_ds['lon_rho'].values[0,:]
#pdb.set_trace()
#roms_ds=xr.open_dataset('West_Coast_Temperature_2019_last4days.nc',decode_cf=True)
#lat=roms_ds['lat']
#lon=roms_ds['lon']
latitude=lat
longitude=lon
dims=['time','lat','lon']
ds=xr.Dataset({'wcofs_temperature':(dims,abigtemp)},
              coords={'time':abigtime,
                     'lon':longitude,
                     'lat':latitude,
                     })
ds.attrs['title']='WCOFS SST UCSC ROMS grid 2022'
ds.attrs['notes']="Created on "+dt.datetime.today().strftime("%&-%m-%d")+" by flbahr, interpolated onto UCSC grid"
fname="WCOFS_SST_2022.nc"
ds.to_netcdf(path=fname)
