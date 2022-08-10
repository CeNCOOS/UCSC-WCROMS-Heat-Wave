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

# date span of data files
# need code to comute these numbers
startdate1=dt.datetime(2019,1,12)
enddate1=startdate1+dt.timedelta(days=8)
ic=0
#
filemid=43471
for i in np.arange(0,27): # had been 92, 67
    if i==0:
        startdate=startdate1
        enddate=enddate1
    else:
        startdate=startdate+dt.timedelta(days=4)
        enddate=startdate+dt.timedelta(days=8)
    filestr=startdate.strftime('%Y%m%d')+'_'+enddate.strftime('%Y%m%d')
    #thefilename='wc12_ccsra31_his_'+filestr+'.nc'
    filemid=filemid+4
    thefilename='fwd_'+str(filemid)+'_001_'+filestr+'.nc'
    #print(thefilename)
    # open the data set
    url='https://oceanmodeling.ucsc.edu:8443/thredds/dodsC/wc12.0_ccsra31_ext_01/posterior/'+thefilename
    
    #url = 'https://oceanmodeling.ucsc.edu:8443/thredds/dodsC/wc12.0_ccsra31_01/posterior/'+thefilename
    roms_ds = xr.open_dataset(url+'#fillmismatch', decode_cf=True) 
    #roms_ds = xr.open_dataset(url, decode_cf=True) 
    # set up the stretching coordinates for the sigma coordinate system
    s_coor_streching_rho = roms_ds['Cs_r'].values # S-coordinate stretching curves at RHO-points
    s_coor_streching_w = roms_ds['Cs_w'].values # S-coordinate stretching curves at w-points
    vertical_stretching_function = roms_ds['Vstretching'] # vertical terrain following stretching function
    vertical_transformation_function = roms_ds['Vtransform'] # vertical terrain following tranformation equation
    ##
    h = roms_ds['h'] # bathymetry at rho points
    h_critical = roms_ds['hc'] # S-coordinate parameter, critical depth
    theta_bottom = roms_ds['theta_b'].values # S-coordinate bottom control parameter
    theta_surface = roms_ds['theta_s'].values # S-coordinate surface control parameter
    zeta = roms_ds['zeta'] # free surface
    N_levels=42 # number of levels, note w has 43 levels
    zeta = np.squeeze(zeta[0,:,:]) # reduce dimensions to    
    ## Calculate z-depths for the rho points   
    ##%%time
    igrid=1 # density point grid, for (T,S) use igrid=3, for u and igrid=4 for v # use igrid=5 for w
    z_rho = set_depth(vertical_transformation_function, vertical_stretching_function, theta_surface, theta_bottom, h_critical, N_levels, igrid, h, zeta, report=0)
    ## Calculate the z-depths for the stream function points
    ##%%time
    igrid=2 # density point grid, for (T,S) use igrid=3, for u and igrid=4 for v # use igrid=5 for w
    z_stream = set_depth(vertical_transformation_function, vertical_stretching_function, theta_surface, theta_bottom, h_critical, N_levels, igrid, h, zeta, report=0)    
    #
    roms_ds = xr.open_dataset(url+'#fillmismatch',decode_cf=True)
    roms_temp = roms_ds['temp']
    roms_time = roms_ds['ocean_time']
    roms_temp_latitude = roms_temp['lat_rho'].values[:,0]
    roms_temp_longitude = roms_temp['lon_rho'].values[0,:]
    roms_temp = roms_temp.assign_coords({"eta_rho":roms_temp_latitude, "xi_rho":roms_temp_longitude})
    roms_temp = roms_temp.rename({"eta_rho":"latitude", "xi_rho":"longitude"})
    roms_temp = roms_temp.drop('lat_rho', errors='ignore')
    roms_temp = roms_temp.drop('lon_rho', errors='ignore')
    roms_temp = roms_temp.drop('time_run', errors='ignore')
    ##roms_temp = roms_temp.transpose('time','latitude','longitude','s_rho')
    roms_temp = roms_temp.transpose('ocean_time','latitude','longitude','s_rho')
    roms_temp['z_rho'] = (('latitude','longitude','s_rho'),z_rho) # Add the depth values back to the Dataset    
    #
    #%%time
    #westcoast_temp = roms_temp.sel(latitude=slice(bbox[0][1], bbox[0][0]), longitude=slice(bbox[1][1], bbox[1][0]), time=slice(last_time,today))
    westcoast_temp = roms_temp.sel(latitude=slice(bbox[0][1], bbox[0][0]), longitude=slice(bbox[1][1], bbox[1][0]))
    #westcoast_heat_content = np.zeros(shape=westcoast_temp[:,:,:,0].shape)

    #for i in tqdm(range(len(westcoast_temp)):
    #    westcoast_temp[i,:,:]=roms_temp[]
    #    upper_heat_content = calculate_upper_heat_content(westcoast_temp.isel(ocean_time=i), upper_depth=-100)
    #    westcoast_heat_content[i,:,:] = upper_heat_content
    #for i in tqdm(range(len(westcoast_heat_content))):
    #    upper_heat_content = calculate_upper_heat_content(westcoast_temp.isel(ocean_time=i), upper_depth=-100)
    #    westcoast_heat_content[i,:,:] = upper_heat_content

    # okay we now have the heat content for a file we only want daily output for this so average
    # since assimilation occurs at the start of each run only need first 4 days?
    if i==66:
        startavg=np.arange(0,32,4)
    else:
        startavg=np.arange(0,16,4)
    for j in startavg:
        avgheat=np.nanmean(westcoast_temp[j:j+4,:,:,-1],axis=0)
    #    # grab a time associated with this 
        datest=np.datetime_as_string(roms_time[j],unit='D')
        parts=datest.split('-')
        centertime=dt.datetime(int(parts[0]),int(parts[1]),int(parts[2]),12,0,0)
    ##     
    ##avgheat=np.nanmean(westcoast_heat_content[0:4,:,:],axis=0)
    ## add third axis to stack along
        avgheat=np.expand_dims(avgheat,axis=2)
        if ic==0:
            the_heat=avgheat
            the_time=centertime
            ic=1
        else:
            the_heat=np.dstack((the_heat,avgheat))
            the_time=np.append(the_time,centertime)
    
    # 0:4, 4:8, 8:12, 12:16, 16:20
    #  9   10   11     12     13
    
#19800109_19800117
#19800113_19800121
#19800117_19800125
#19800121_19800129
#19800125_19800202
#19800129_19800206
#19800202_19800210
#19800206_19800214
#19800210_19800218

data = the_heat
time = the_time
longitude = westcoast_temp.longitude.values
latitude = westcoast_temp.latitude.values
dims = ['lat', 'lon','time']
ds = xr.Dataset({'temperature_2019_surface': (dims, the_heat)},
                coords={"time":time,
                        "lon": longitude,
                        "lat": latitude,
                        })
ds.attrs['title'] = "West Coast Upper SST 2019 - 0 meters"
ds.attrs['notes'] = "Created on "+dt.datetime.today().strftime("%Y-%m-%d") + " by flbahr"

fname = "West_Coast_Temperature_2019_last4days.nc"
ds.to_netcdf(path=fname)

print(the_time[-1])
print(thefilename)
