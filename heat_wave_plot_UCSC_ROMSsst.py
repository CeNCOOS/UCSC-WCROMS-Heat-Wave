import sys
sys.path.append("..")
import numpy as np
import xarray as xr
# from shapely.geometry import Polygon, Point
import scipy.io
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter)
from scipy import integrate
import pandas as pd
import datetime as dt
import time, cmocean, scipy.io
from tqdm import tqdm
import seaborn as sns
import geopandas as gpd
from statsmodels.tsa.seasonal import seasonal_decompose
import pocean
import copy
import pdb
ds=xr.open_dataset('/home/flbahr/heat_content/West_Coast_upper_heat_content_2019_last4days.nc')
#ds=xr.open_dataset('c:/Users/flbahr/heatcontent_heatwave_2sst.nc')
#ds=xr.open_dataset('c:/Users/flbahr/heatcontent_heatwave_20220808.nc')
#dsst=xr.open_dataset('x:/heat_content/heat_wave_data/heatcontent_heatwave_2sst_20220809.nc')
import pickle
tfile=open('/home/flbahr/heat_content/heat_wave_westcoastsst_time_fore.p','rb')
[time,newtime]=pickle.load(tfile)
del newtime
efile=open('/home/flbahr/heat_content/heat_wavesst_extreme_fore.p','rb')
sefile=open('/home/flbahr/heat_content/heat_wavesst_severe_fore.p','rb')
stfile=open('/home/flbahr/heat_content/heat_wavesst_strong_fore.p','rb')
mfile=open('/home/flbahr/heat_content/heat_wavesst_moderate_fore.p','rb')
[heatextreme]=pickle.load(efile)
[heatsevere]=pickle.load(sefile)
[heatstrong]=pickle.load(stfile)
[heatmoderate]=pickle.load(mfile)
#pdb.set_trace()
#mod=ds['heatmoderate'][:,:,:]
#sto=ds['heatstrong'][:,:,:]
#sev=ds['heatsevere'][:,:,:]
#ext=ds['heatextreme'][:,:,:]

zmod=np.mean(heatmoderate,axis=(1,2))
zsto=np.mean(heatstrong,axis=(1,2))
zsev=np.mean(heatsevere,axis=(1,2))
zext=np.mean(heatextreme,axis=(1,2))

lat=ds['lat']
lon=ds['lon']
#time=ds['time']

from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import cm
import matplotlib.ticker as mticker

amap=plt.cm.YlOrRd
amaplist=[amap(i) for i in range(amap.N)]
amaplist[0]=(1,1,1,1)
#bounds=[0,1,2,3,4,5]
bounds=[0,0.5,1.5,2.5,3.5,4.5]
cmap=mpl.colors.LinearSegmentedColormap.from_list('Custom cmap',amaplist,amap.N)
bounds
norm=mpl.colors.BoundaryNorm(bounds,cmap.N)
[n,il1,il2]=heatmoderate.shape
#pdb.set_trace()
# suppress too many plots open warning
plt.rcParams.update({'figure.max_open_warning': 0})
#
#
plt.rcParams['figure.figsize']=[12,10]
iz=np.arange(n-60,n) # we used to have this as 180
#iz=np.arange(11857,n)
for i in iz:
    zzz=str(i)+' '+str(n)
    print(zzz)
    fig=plt.figure()
    strong=heatstrong[i,:,:]
    moderate=heatmoderate[i,:,:]
    severe=heatsevere[i,:,:]
    extreme=heatextreme[i,:,:]
#    strong=ds['heatstrong'][i,:,:]
#    moderate=ds['heatmoderate'][i,:,:]
#    severe=ds['heatsevere'][i,:,:]
#    extreme=ds['heatextreme'][i,:,:]
    thewave=moderate+strong*2+severe*3+extreme*4
    #ax0=plt.subplot2grid((5,5),(0,0),colspan=4,rowspan=4,projection=ccrs.PlateCarree())
    ax0=plt.subplot2grid((6,5),(0,0),colspan=4,rowspan=4,projection=ccrs.PlateCarree())
    g1=ax0.gridlines(draw_labels=True)
    g1.top_labels=False
    #g1.xlabels_top=False
    g1.right_labels=False
    #g1.ylabels_right=False
    ax0.set_xlim(-129.1, -117.5)
    g1.xlocator=mticker.FixedLocator([-129,-127,-125,-123,-121,-119,-117])
    g1.xlabel_style={'size':15}
    g1.ylabel_style={'size':15}
    #g1.ylabel('Longitude')
    #g1.xlabel('Latitude')
    ax0.set_ylim(30.6,48)
    #ax0.yaxis.set_major_formatter(test)
    #ax0.set_ylabel('Latitude')
    #t=time[i].values
    t=time[i]
    tt=pd.to_datetime(str(t))
    ts=tt.strftime('%d %b %Y')
    if i < n-2: # changed from n-2 to n-1 for when to use the red color
        plt.title(ts, fontsize=18)
    else:
        plt.title(ts, fontsize=18, color='red')
    im=ax0.pcolor(lon,lat,np.where(thewave==0, np.nan, thewave),cmap=cmap,norm=norm,shading='auto')
    #im=ax0.pcolor(lon,lat,np.where(thewave==0, np.nan, thewave),cmap=cmap,norm=norm,vmin=0,vmax=5,shading='auto')
    #im=ax0.pcolor(lon,lat,np.where(thewave==0, np.nan, thewave),cmap=cmap,norm=norm,vmin=0,vmax=5)
    ax0.coastlines('10m')
    ax0.add_feature(cfeature.NaturalEarthFeature('physical','land','10m', edgecolor='face',facecolor='#808080'))
    #g1.xformatter=LongitudeFormatter()
    #g1.yformatter=LatitudeFormatter()

    cbar=plt.colorbar(im,ticks=[0,1,2,3,4,5])
    #cbar=plt.colorbar(im,ticks=[0,0.5,1.5,2.5,3.5,4.5])
    #cbar.ax.set_xticks([0.5,1.5,2.5,3.5,4.5], minor=True)
    cbar.ax.set_yticklabels(['None','Moderate','Strong','Severe','Extreme'])
    cbar.ax.tick_params(labelsize=15)
    #ax1=plt.subplot2grid((5,5),(4,1),colspan=3,rowspan=1)
    ax1=plt.subplot2grid((6,5),(5,1),colspan=3,rowspan=1)
    #ax1.plot(ds['time'],zmod+zsto+zsev+zext)
    ax1.plot(time,zmod+zsto+zsev+zext)
    ax1.plot(time[i],zmod[i]+zsto[i]+zsev[i]+zext[i],'r*')
    #ax1.plot(ds['time'][i],zmod[i]+zsto[i]+zsev[i]+zext[i],'r*')
    ax1.tick_params(axis='x', labelsize=15)
    ax1.tick_params(axis='y', labelsize=15)
    #ax1.set_ylabel('Heatwave Intensity',fontsize=15)
    #ax1.left_labels=False
    ax1.get_yaxis().set_ticks([])
    ax1.set_ylabel('Area Intensity',fontsize=15)
    plt.subplots_adjust(hspace=-.4)
    #ax1.ylabel_style={'size':15}
    #plt.savefig('x:/heat_content/heat_contnet_sst_'+str(i)+'_refine.png',dpi=300,bbox_inches='tight', pad_inches=0.25)
    #plt.savefig('x:/heat_content/heat_wave_plots/heat_contnet_sst_'+str(i)+'_refine.png',dpi=300,bbox_inches='tight', pad_inches=0.25)
    plt.savefig('/home/flbahr/heat_content/heat_wave_sst_'+str(i)+'_refine_wcofs.png',dpi=300,bbox_inches='tight', pad_inches=0.25)
    plt.clf()
    #plt.savefig('heat_test_to_refine.png',dpi=300,bbox_inches='tight', pad_inches=0.25)

