# file to merge data and output?
import xarray as xr
import numpy as np
import scipy.io
from scipy import integrate
from scipy.interpolate import griddata
import datetime as dt
import time, scipy.io
import pandas as pd
import pickle
import pdb
from netCDF4 import Dataset, date2num, num2date
#
# Need to get current time
atime=dt.datetime.now()
#
#current_year=2022
current_year=atime.year
#current_month=12
current_month=atime.month
#
# deal with case of january
if current_month==1:
    current_year=current_year-1
    lastmonth=13
else:
    lastmonth=current_month
#
#nummonth=np.arange(0,lastmonth-1)
nummonth=np.arange(0,lastmonth)
#
mons=['january','february','march','april','may','june','july','august','september','october','november','december']
abigtime=[]
abigtemp=[]
for i in nummonth:
    fname='/home/flbahr/heat_content/wcofs_'+mons[i]+'_'+str(current_year)+'.p'
    thefile=open(fname,'rb')
    [time,temp]=pickle.load(thefile)
    if i==0:
        abigtime=time
        abigtemp=temp
    else:
        abigtime=np.append(abigtime,time)
        abigtemp=np.vstack((abigtemp,temp))
#
# this is just to load the latitude and longitude arrays so year can be hard coded.
roms_ds=xr.open_dataset('/home/flbahr/heat_content/WCOFS_SST_2022.nc',decode_cf=True)
lat=roms_ds['lat']
lon=roms_ds['lon']
latitude=lat
longitude=lon
#
string_year=str(current_year)
dims=['time','lat','lon']
ds=xr.Dataset({'wcofs_temperature':(dims,abigtemp)},
              coords={'time':abigtime,
                     'lon':longitude,
                     'lat':latitude,
                     })
ds.attrs['title']='WCOFS SST UCSC ROMS grid '+string_year
ds.attrs['notes']="Created on "+dt.datetime.today().strftime("%&-%m-%d")+" by flbahr, interpolated onto UCSC grid"
fname="/home/flbahr/heat_content/WCOFS_SST_"+string_year+".nc"
ds.to_netcdf(path=fname)
