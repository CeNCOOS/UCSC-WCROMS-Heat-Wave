#from pydap.client import open_url
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
from wcofs_lonlat_2_xy import wcofs_lonlat_2_xy
import pdb
#
#
# Weird Bug here, but make sure to run this first
import salem
from shapely import geometry
sanctuary_outline= salem.read_shapefile('./shapefiles/USMaritimeLimitsAndBoundariesSHP/USMaritimeLimitsNBoundaries.shp')
sanctuary_outline.crs = 'epsg:4326'
eez_shape = sanctuary_outline[(sanctuary_outline['REGION'] == "Pacific Coast") & (sanctuary_outline['EEZ'])]
bbox = ((eez_shape['max_y'].max(),eez_shape['min_y'].min()),(eez_shape['max_x'].max(),eez_shape['min_x'].min()))
geom = geometry.box(minx=bbox[1][1],maxx=bbox[1][0],miny=bbox[0][1],maxy=bbox[0][0])
#
#
url='https://oceanmodeling.ucsc.edu:8443/thredds/dodsC/ccsra_2016a_phys_agg_slevs/fmrc/CCSRA_2016a_Phys_ROMS_Sigma-level_Aggregation_best.ncd'
roms_ds = xr.open_dataset(url, decode_cf=True) 
roms_tempUCSC = roms_ds['temp'][:,-1,:,:]
roms_temp_latitudeUCSC = roms_tempUCSC['lat_rho'].values[:,0]
roms_temp_longitudeUCSC = roms_tempUCSC['lon_rho'].values[0,:]

[mx,my]=np.meshgrid(roms_temp_longitudeUCSC,roms_temp_latitudeUCSC)

#
atime=dt.datetime.now()
#
# move this outside loop so we only have to compute this once and use multiple times
# forecast hour for 3d fields
#fhour=np.arange(3,75,3)
#nhour=np.arange(3,27,3)
# for 2d surface fields
fshour=np.arange(1,73,1)
nshour=np.arange(1,25,1)
urlpre='https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/WCOFS/MODELS/'
#filestr='nos.wcofs.fields.'
filestr='nos.wcofs.2ds.'
filepost='.t03z.nc'
# define a loop in time
#
dataset1=xr.open_dataset("https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/WCOFS/MODELS/2022/09/11/nos.wcofs.2ds.n001.20220911.t03z.nc")
dataset2=xr.open_dataset("https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/WCOFS/MODELS/2022/09/11/nos.wcofs.2ds.f001.20220911.t03z.nc")
dataset3=xr.open_dataset("https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/WCOFS/MODELS/2022/09/10/nos.wcofs.2ds.f025.20220910.t03z.nc")
dataset4=xr.open_dataset("https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/WCOFS/MODELS/2022/09/09/nos.wcofs.2ds.f049.20220909.t03z.nc")

#dataset1=xr.open_dataset("https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/WCOFS/MODELS/2022/08/23/nos.wcofs.2ds.n024.20220823.t03z.nc")
#dataset2=xr.open_dataset("https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/WCOFS/MODELS/2022/08/22/nos.wcofs.2ds.f048.20220822.t03z.nc")
#dataset3=xr.open_dataset("https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/WCOFS/MODELS/2022/08/24/nos.wcofs.2ds.n001.20220824.t03z.nc")
#
#
t1=dataset1['temp_sur']
t2=dataset2['temp_sur']
t3=dataset3['temp_sur']
t4=dataset4['temp_sur']
lat=dataset1['lat_rho']
lon=dataset2['lon_rho']
t1=np.squeeze(t1)
t2=np.squeeze(t2)
t3=np.squeeze(t3)
t4=np.squeeze(t4)
pdb.set_trace()
thedays=np.arange(2,-1,-1)
io=0
for offset in thedays:
    print(offset)
    thedate=atime-dt.timedelta(days=int(offset))
#starttime=atime-dt.timedelta(days=10)
    # compute the parts so we can construct the url and filename
    year=thedate.year
    month=thedate.month
    day=thedate.day
    syear=str(year)
    if month < 10:
        smonth='0'+str(month)
    else:
        smonth=str(month)
    if day < 10:
        sday='0'+str(day)
    else:
        sday=str(day)
    #
    # hour output?
    # n=nowcast? f=forecast
    # n024, n021, n018, n015, n012, n009, n006, n003
    # f072, f069, f066, f063, f060, f057, f054, f051, f048, f045, f042, f039, f036, f033, f030, f027, f024, f021, f018 ,f015, f012, f009, f006, f003
    #
    # fields
    # lat_rho, lon_rho, ocean_time, mask_rho, temp, salt

    # variable names 2ds
    # temp_sur, salt_sur, u_sur, v_sur, lat_rho, lon_rho, ocean_time, mask_rho
    for zhours in fshour:
        if zhours < 10:
            forehour='f00'+str(zhours)
        else:
            forehour='f0'+str(zhours)
        urlstr=urlpre+syear+'/'+smonth+'/'+sday+'/'+filestr+forehour+'.'+syear+smonth+sday+filepost
    # ocean time is relative to 2016-01-01
        dataset=xr.open_dataset(urlstr)
    #dataset=xr.open_dataset("https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/WCOFS/MODELS/2022/08/23/nos.wcofs.fields.n024.20220823.t03z.nc")
    #dataset=xr.open_dataset("https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/WCOFS/MODELS/2022/08/23/nos.wcofs.2ds.n024.20220823.t03z.nc")
    #dataset=open_url("https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/WCOFS/MODELS/2022/08/23/nos.wcofs.fields.n024.20220823.t03z.nc")
#
#
#
        lat=dataset['lat_rho']
        lon=dataset['lon_rho']
    #xy=wcofs_lonlat_2_xy(lon,lat,1);
#
#
#
#pdb.set_trace()
#
        otime=dataset['ocean_time']

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
#
#y=np.ravel(xy['x'])
#x=np.ravel(xy['y'])
#x=-x
        t=np.ravel(roms_temp)
#xgrid=griddata((x,y),t,(mx,my),method='linear')
        xl=np.ravel(lon)
        yl=np.ravel(lat)
#ygrid=griddata((xl,yl),t,(mx,my),method='linear')
        zgrid=griddata((xl,yl),t,(mx,my),method='nearest')
        if io==0:
            bigtemp=np.expand_dims(zgrid,axis=0)
            bigtime=otime
            io=1
        else:
            bigtime=np.append(bigtime,otime)
            btmp=np.expand_dims(zgrid,axis=0)
            bigtemp=np.vstack((bigtemp,btmp))
#levels=np.arange(10,26,2)
#plt.subplot(131)
#plt.contourf(mx,my,ygrid,levels=levels)
##plt.clim(7.5, 27.5)
#plt.colorbar()
#plt.title('Linear')
#plt.subplot(132)
#plt.contourf(mx,my,zgrid,levels=levels)
##plt.clim(7.5,27.5)
#plt.colorbar()
#plt.title('Nearest Neighbor')
#plt.subplot(133)
#plt.contourf(lon,lat,roms_temp,levels=levels)
#plt.xlim([-134, -115.6])
#plt.ylim([30,48])
#plt.title('Original')
##plt.clim(10,24)
#plt.colorbar()
##plt.show()
#plt.savefig('wcofs_test.png')
pdb.set_trace()
#xxgrid=griddata((xy['x'],xy['y']),roms_temp,(mx,my),method='linear')
#
# I don't think this "slice" really works properly for this grid
#
#westcoast_temp = roms_temp.sel(latitude=slice(bbox[0][1], bbox[0][0]), longitude=slice(bbox[1][1], bbox[1][0]))
#westcoast_temp=np.squeeze(westcoast_temp)
pdb.set_trace()

#
#
#
def calculate_upper_heat_content(temperature, upper_depth=-100 ):
    """ 
    Calculate upper ocean heat content from ROMS
    
    Calculate the upper ocean heat (J/m^2) content from the CA ROMS model output
    Density is assumed constant at 1025 kg/m^3
    specific heat is 3850 J/(kg C)
    Args:
        temperature: xarray dataset with dimensions (latitude,longitude,s_sho) 
        z_rho: array mapping sigma depth bins to depths in meters with dimensions (latitude,longitude,s_sho) 
        upper_depth: depth to integrate to from the surface, values below the surface are negative

    Returns:
        upper_heat_content: array of upper heat content values in the shape of the temperature input dimensions (latitude, longitude)
        upper_heat_content_b_bins: array of the number of depth bins used for integrated upper heat content in the shape of the temperature input dimensions (latitude, longitude)
    """
    z_rho = temperature.z_rho.values
    cp = 3850
    density = 1025
    
    if len(z_rho.shape) == 1:
        '''Calculating at single point'''
        upper_heat_content = np.zeros(shape=(1))
        roms_temperature_flattened = temperature.values.reshape(-1, temperature.shape[-1])
        
    else:
        upper_heat_content = np.zeros(shape=(z_rho.shape[0]*z_rho.shape[1]))
        roms_temperature_flattened = temperature.values.reshape(-1, temperature.shape[-1])
        
    for i, grid_point in enumerate(z_rho.reshape(-1, z_rho.shape[-1])):
        #Find index of depth less than 100 meters
        upper_ix = np.where(grid_point >= upper_depth)[0]
        # Integrate
        if not upper_ix.size == 0:
            
            interp_top = np.interp(0,grid_point,roms_temperature_flattened[i,:])
            interp_bottom = np.interp(-100, grid_point,roms_temperature_flattened[i,:])
            # Add interpolated values to temperature array for integrtaion
            upper_values = roms_temperature_flattened[i,upper_ix]
            upper_values = np.append(interp_bottom,upper_values)
            upper_values = np.append(upper_values,interp_top)
            # Add interpolated depths to depth array for integrtaion
            upper_depths = grid_point[upper_ix]
            upper_depths = np.append(-100,upper_depths)
            upper_depths = np.append(upper_depths,0)
            integrated_temp = integrate.simps(upper_values,upper_depths)                    
            upper_heat_content[i] = density * cp * integrated_temp
            
        else:
            upper_heat_content[i] = np.nan
    
    if len(z_rho.shape) == 1:
        return upper_heat_content
    
    else:
        upper_heat_content = upper_heat_content.reshape(z_rho.shape[:2])
        return upper_heat_content
#
#
#
s_coor_streching_rho = dataset['Cs_r'].values # S-coordinate stretching curves at RHO-points
s_coor_streching_w = dataset['Cs_w'].values # S-coordinate stretching curves at w-points
vertical_stretching_function = dataset['Vstretching'] # vertical terrain following stretching function
vertical_transformation_function = dataset['Vtransform'] # vertical terrain following tranformation equation
    #
h = dataset['h'] # bathymetry at rho points
h_critical = dataset['hc'] # S-coordinate parameter, critical depth
theta_bottom = dataset['theta_b'].values # S-coordinate bottom control parameter
theta_surface = dataset['theta_s'].values # S-coordinate surface control parameter
zeta = dataset['zeta'] # free surface
N_levels=40 # number of levels, note w has 43 levels
zeta = np.squeeze(zeta[0,:,:]) # reduce dimensions to    
# Calculate z-depths for the rho points   
#%%time
igrid=1 # density point grid, for (T,S) use igrid=3, for u and igrid=4 for v # use igrid=5 for w
z_rho = set_depth(vertical_transformation_function, vertical_stretching_function, theta_surface, theta_bottom, h_critical, N_levels, igrid, h, zeta, report=0)
roms_temp = dataset['temp']
roms_time = dataset['ocean_time']
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

westcoast_temp = roms_temp.sel(latitude=slice(bbox[0][1], bbox[0][0]), longitude=slice(bbox[1][1], bbox[1][0]))
westcoast_heat_content = np.zeros(shape=westcoast_temp[:,:,:,0].shape)
for i in tqdm(range(len(westcoast_heat_content))):
    upper_heat_content = calculate_upper_heat_content(westcoast_temp.isel(ocean_time=i), upper_depth=-100)
    westcoast_heat_content[i,:,:] = upper_heat_content
