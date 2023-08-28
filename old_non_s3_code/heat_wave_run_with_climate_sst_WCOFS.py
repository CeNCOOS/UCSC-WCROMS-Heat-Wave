# Load supporting libraries we need
#
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
from load_wcofs_for_heatwave import wcof_load
import pdb
#
#
import salem
from shapely import geometry
sanctuary_outline= salem.read_shapefile('/home/flbahr/heat_content/shapefiles/USMaritimeLimitsAndBoundariesSHP/USMaritimeLimitsNBoundaries.shp')
sanctuary_outline.crs = 'epsg:4326'

eez_shape = sanctuary_outline[(sanctuary_outline['REGION'] == "Pacific Coast") & (sanctuary_outline['EEZ'])]
bbox = ((eez_shape['max_y'].max(),eez_shape['min_y'].min()),(eez_shape['max_x'].max(),eez_shape['min_x'].min()))
geom = geometry.box(minx=bbox[1][1],maxx=bbox[1][0],miny=bbox[0][1],maxy=bbox[0][0])
#
# get current time
#
now=dt.datetime.now()
#
# find 30 years prior
#
cyear=now.year-32
cstart=dt.datetime(cyear,now.month,now.day,now.hour,now.minute,now.second)
#
# load the files from cstart to present
#
ny=np.arange(cyear,2020)
#
import marineHeatWave as marineHeatWave
#
ic=0
for i in ny:
    file='/home/flbahr/heat_content/West_Coast_Temperature_'+str(i)+'_last4days.nc'
    ds=xr.open_dataset(file)
    print(file)
    tz=ds['time']
    varname='temperature_'+str(i)+'_surface'
    th=ds[varname]
    if ic==0:
        bigarray=th
        bigtime=tz
        ic=1;
    else:
        bigarray=np.dstack((bigarray,th))
        bigtime=np.append(bigtime,tz)

#
# load current aggregation of temperatue data for West Coast ROMS
# 
thefile='/home/flbahr/heat_content/West_Coast_Temperature_surface_present.nc'
#thefile='/home/pdaniel/heat-content/data/West_Coast_upper_heat_content.nc'
ds1=xr.open_dataset(thefile)

th1=ds1['temperature_surface']
tz1=ds1['time']
# remove values that are already in the 30 year reanalysis
# this still contains data that overlaps with the WCOFS data below
ind=np.where(tz1 > bigtime[-1])
tz1=tz1[ind]
th1=th1[ind[0],:,:]
# load the WCOFS data
wcof1='/home/flbahr/heat_content/West_Coast_Temperature_2021_WCOFS.nc'
wcof2='/home/flbahr/heat_content/West_Coast_Temperature_2022_WCOFS.nc'
wcf1=xr.open_dataset(wcof1)
wcf2=xr.open_dataset(wcof2)

wh1=wcf1['temperature_2021_surface']
wt1=wcf1['time']
wh2=wcf2['temperature_2022_surface']
wt2=wcf2['time']
# remove overlap (so we just have cover the gap between 2019 and 2021
ind2=np.where(tz1 < wt1[0])
tz1=tz1[ind2]
th1=th1[ind2[0],:,:]

tz1=np.append(tz1,wt1)
tz1=np.append(tz1,wt2)

th1=np.vstack((th1,wh1))
th1=np.vstack((th1,wh2))

lat=ds1['lat'].values
lon=ds1['lon'].values
[mx,my]=np.meshgrid(lon,lat)

[timearray,temparray]=wcof_load(tz1,mx,my)
# find common times
# need to convert from timestamp to datetime64
lt=len(timearray)
a=np.arange(0,lt)
newt=[]
for i in a:
    xxx=timearray[i].to_datetime64()
    if i==0:
        newt=xxx
    else:
        newt=np.append(newt,xxx)

ind3=np.where(newt > tz1[-1])
newt=newt[ind3]
#timearray=timearray[ind3]
temparray=temparray[ind3[0],:,:]
th1=np.vstack((th1,temparray))
tz1=np.append(tz1,newt)

#thefile2='/home/pdaniel/heat-content/data/West_Coast_upper_heat_content_recent.nc'
#ds2=xr.open_dataset(thefile)

#th2=ds2['heat_content_100_meters']
#tz2=ds2['time']
# remove duplicate times
#values,indext,count=np.unique(tz1,return_index=True,return_counts=True)
#newt=tz1[indext]
#newh=th1[indext,:,:].values
#pdb.set_trace()

# change shape to append 
newh=th1.transpose((1,2,0))
#newh=newh.transpose((1,2,0))
oldarray=bigarray
oldtime=bigtime
bigarray=np.dstack((bigarray,newh))
bigtime=np.append(bigtime,tz1)


[nlo,nla,nt]=bigarray.shape
alat=ds['lat']
alon=ds['lon']
#
# convert time to ordinals for heatwave code
#
time=bigtime
newtime=[]
it=np.arange(0,len(time))
for t in it:    
    tx=dt.datetime.utcfromtimestamp(time[t].tolist()/1e9)
    tx=dt.datetime.toordinal(tx)
    newtime=np.append(newtime,tx)
newtime=newtime.astype('int')
#pdb.set_trace()
#
# set up arrays for storing heatwave computation results
#
heatwavesnumber=np.zeros((175,117))
hwv=np.zeros((175,117))
hwevents=np.zeros((175,117))
hwmoderate=np.zeros((175,117))
hwstrong=np.zeros((175,117))
hwsevere=np.zeros((175,117))
hwextreme=np.zeros((175,117))
startindexs=[]
stopindexs=[]
heatvalueindex=np.zeros((nt,nlo,nla))
heatmoderate=np.zeros((nt,nlo,nla))
heatstrong=np.zeros((nt,nlo,nla))
heatsevere=np.zeros((nt,nlo,nla))
heatextreme=np.zeros((nt,nlo,nla))

# Compute the heat-wave!!

ix=np.arange(0,175)
iy=np.arange(0,117)
# first actual data at i=69, j=59
for i in ix:
    print(i)
    for j in iy:
        xtest=bigarray[i,j,:]
        #xtest=test.values        
        # how to check for NaN
        if np.isnan(xtest[0])==False:
            # we have data to run the code on
            # here is where we can run the marine Heat Wave code
            mhw=marineHeatWave.detect(newtime,xtest)
            stats=mhw[0]
            other=mhw[1]
            indexs=stats['index_start']
            indexe=stats['index_end']
            idiff=(np.array(indexe)-np.array(indexs))+1
            heatwavesnumber[i,j]=sum(idiff)
            hwv[i,j]=sum(stats['duration'])
            hwmoderate[i,j]=sum(stats['duration_moderate'])
            hwstrong[i,j]=sum(stats['duration_strong'])
            hwsevere[i,j]=sum(stats['duration_severe'])
            hwextreme[i,j]=sum(stats['duration_extreme'])
            hwevents[i,j]=stats['n_events']
            startindexs.append(mhw[0]['index_start'])
            stopindexs.append(mhw[0]['index_end'])
            # keep track of when a heat wave occurs.  Need to know type of heat wave also?
            category=stats['category']
            lid=len(indexs)
            for k in np.arange(0,lid):
                heatvalueindex[indexs[k]:indexe[k]+1,i,j]=1
                if category[k]=='Moderate':
                    heatmoderate[indexs[k]:indexe[k]+1,i,j]=1
                if category[k]=='Strong':
                    heatstrong[indexs[k]:indexe[k]+1,i,j]=1
                if category[k]=='Severe':
                    heatsevere[indexs[k]:indexe[k]+1,i,j]=1
                if category[k]=='Extreme':
                    heatextreme[indexs[k]:indexe[k]+1,i,j]=1

#
# Do we want to save here incase of problems?
#
import pickle
#pickle.dump([bigarray],open("/home/flbahr/heat_content/sstarray.p","wb"))
del oldarray
del bigarray
pickle.dump([time,newtime],open("/home/flbahr/heat_content/heat_wave_westcoastsst_time_fore.p","wb"))
pickle.dump([heatextreme],open("/home/flbahr/heat_content/heat_wavesst_extreme_fore.p","wb"))
pickle.dump([heatsevere],open("/home/flbahr/heat_content/heat_wavesst_severe_fore.p","wb"))
pickle.dump([heatstrong],open("/home/flbahr/heat_content/heat_wavesst_strong_fore.p","wb"))
pickle.dump([heatmoderate],open("/home/flbahr/heat_content/heat_wavesst_moderate_fore.p","wb"))
pickle.dump([heatwavesnumber,hwv,hwmoderate,hwstrong,hwsevere,hwextreme,hwevents,startindexs,stopindexs,heatvalueindex], open("/home/flbahr/heat_content/heat_wavesst_meanvalues_fore.p","wb"))            
                                 

#ind=np.where(np.diff(tz1)==0)
#td=newh1.values
#td=td.transpose((1,2,0))
#oldarray=bigarray;
#oldtime=bigtime;
#bigarray=np.dstack((bigarray,td))
#bigtime=np.append(bigtime,newz1)
