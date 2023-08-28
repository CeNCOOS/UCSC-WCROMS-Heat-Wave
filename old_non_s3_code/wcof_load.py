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
import pdb
#
# atime is now, thedays
def wcof_load(atime,offset,numbs,mx,my):
    io=0
    # this goes upto 3 days behind current time
    urlpre='https://ncei.noaa.gov/thredds/dodsC/model-wcofs-files/'
    # this contains forecast but doesn't go too far back in time
    urlpre1='https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/WCOFS/MODELS/'
    filestr='nos.wcofs.2ds.'
    filepost='.t03z.nc'
    # we have an issue that nos.wcofs.2ds.n001.20220913.t03z.nc is not at time 001 but at time 2022-09-13T04:00:00
    #
    thedate=atime-dt.timedelta(days=int(offset))
    # compute the parts so we can construct the url and filename
    year=thedate.year
    month=thedate.month
    day=thedate.day
    print(str(year)+'/'+str(month)+'/'+str(day))
    # okay now we have year month and day of our offset date so we want to find the hour 1 to hour 24 values
    datetoget=dt.datetime(year,month,day)
    # so datetoget is the date we want to get but to get the actual hours for that day we need to get 3 hours earlier
    astart_date=datetoget-dt.timedelta(hours=3)
    # set up appropiate loop?
    if numbs > 24:
        thehours=np.arange(0,72,1)
        forns='f'
    else:
        thehours=np.arange(0,24,1)
        forns='n'
    
    for zhours in thehours:
        # compute offset from start
        timetoget=astart_date+dt.timedelta(hours=int(zhours))
        year=timetoget.year
        month=timetoget.month
        day=timetoget.day
        hour=timetoget.hour
        # case to get 24 hr values and not have n000 or f000 since we have n024 and f024
        if hour==0:
            timetoget=astart_date+dt.timedelta(hours=int(zhours-1))
            year=timetoget.year
            month=timetoget.month
            day=timetoget.day
            hour=24
        # end if
            
        syear=str(year)
        if month < 10:
            smonth='0'+str(month)
        else:
            smonth=str(month)
        if day < 10:
            sday='0'+str(day)
        else:
            sday=str(day)
        if hour < 10:
            forehour=forns+'00'+str(hour)
        else:
            forehour=forns+'0'+str(hour)
# This is for the NCEI site
        urlstr=urlpre+syear+'/'+smonth+'/'+filestr+forehour+'.'+syear+smonth+sday+filepost
# this works for CO-OPS site
#        urlstr=urlpre+syear+'/'+smonth+'/'+sday+'/'+filestr+forehour+'.'+syear+smonth+sday+filepost
        print(urlstr)
    # ocean time is relative to 2016-01-01
    # It turns out we can have missing nowcast files at ncei so we need to trap for that
        try:
            dataset=xr.open_dataset(urlstr)
    #dataset=xr.open_dataset("https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/WCOFS/MODELS/2022/08/23/nos.wcofs.fields.n024.20220823.t03z.nc")
    #dataset=xr.open_dataset("https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/WCOFS/MODELS/2022/08/23/nos.wcofs.2ds.n024.20220823.t03z.nc")
    #dataset=open_url("https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/WCOFS/MODELS/2022/08/23/nos.wcofs.fields.n024.20220823.t03z.nc")
#
            lat=dataset['lat_rho']
            lon=dataset['lon_rho']
#
            otime=dataset['ocean_time']
            #print(urlstr)
            #print(otime)
            roms_temp = dataset['temp_sur']
            roms_time = dataset['ocean_time']
            roms_temp_latitude = roms_temp['lat_rho'].values[:,0]
            roms_temp_longitude = roms_temp['lon_rho'].values[0,:]
            roms_temp = roms_temp.assign_coords({"eta_rho":roms_temp_latitude, "xi_rho":roms_temp_longitude})
            roms_temp = roms_temp.rename({"eta_rho":"latitude", "xi_rho":"longitude"})
            roms_temp = roms_temp.drop('lat_rho', errors='ignore')
            roms_temp = roms_temp.drop('lon_rho', errors='ignore')
            roms_temp = roms_temp.drop('time_run', errors='ignore')
            roms_temp = roms_temp.transpose('ocean_time','latitude','longitude')
            roms_temp=np.squeeze(roms_temp)
#
            t=np.ravel(roms_temp)
            xl=np.ravel(lon)
            yl=np.ravel(lat)
            zgrid=griddata((xl,yl),t,(mx,my),method='linear')
            if io==0:
                #testtemp=np.expand_dims(roms_temp,axis=0)
                bigtemp=np.expand_dims(zgrid,axis=0)
                bigtime=otime
                io=1
            else:
                bigtime=np.append(bigtime,otime)
                btmp=np.expand_dims(zgrid,axis=0)
                #ttmp=np.expand_dims(roms_temp,axis=0)
                bigtemp=np.vstack((bigtemp,btmp))
                #testtemp=np.vstack((testtemp,ttmp))
        except:
            print(urlstr+' is missing from set')
    # average in time
    meantemp1=np.mean(bigtemp,axis=0)
    meantime=pd.Series(bigtime).mean()
    # No need to keep big array and average as differnce between interpolated data and avraged big array that is
    # downsampled is in the 1e-6 range.  So in the noise of the computation.
    #meantemp2=np.mean(testtemp,axis=0)
    #tt=np.ravel(meantemp2)
    #qgrid=griddata((xl,yl),tt,(mx,my),method='linear')
    #pdb.set_trace()
    return [meantime, meantemp1]
