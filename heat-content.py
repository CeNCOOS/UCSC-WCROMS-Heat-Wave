from set_depth import set_depth
import numpy as np
import xarray as xr
from scipy import integrate
import pandas as pd
import datetime as dt
from tqdm import tqdm
import seaborn as sns
# Weird Bug here, but make sure to run this first
import salem, os, logging
import pdb
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

logging.basicConfig(filename='/home/flbahr/heat_content/hc-log.log', level=logging.DEBUG)

FILENAME = "/home/flbahr/heat_content/data/West_Coast_upper_heat_content.nc"

def load_eez_shapefile():
    """ EEZ shapefile will be used as a mask for the model to reduce the size of the stored data and focus on areas of interest. 

    Return the bounding box of the EEZ shape.
    """
    
    eez_outline = salem.read_shapefile('/home/flbahr/heat_content/shapefiles/USMaritimeLimitsAndBoundariesSHP/USMaritimeLimitsNBoundaries.shp')
    eez_outline.crs = 'epsg:4326'
    eez_shape = eez_outline[(eez_outline['REGION'] == "Pacific Coast") & (eez_outline['EEZ'])]
    bbox = ((eez_shape['max_y'].max(),eez_shape['min_y'].min()),(eez_shape['max_x'].max(),eez_shape['min_x'].min()))
    return bbox


def get_last_time():
    """ Check if the generated heat content data is recent (within the past day). If so, there is no need to run the model. If not, we want to run the model and also know when the heat content dataset ends.

    Returns:
        boolean: True if the model is not recent. False if the model is up-to-date.
        tuple: tuple of the last model datetime (DateTime) and today (DateTime).
    """
    try:
        ds = xr.open_dataset(FILENAME)
        last_time = pd.to_datetime(ds['time'].isel(time=-1).values)
        today = pd.to_datetime("today")
        if (today - last_time > dt.timedelta(days=1)):
            return True, (last_time, today)
        else:
            return False, None
    except:
        return True, None

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

def get_depth(ds):
    """Calulate descrete depth bins from sigma-coordinate.
    Args:
        ds (xarray.Dataset): 4DVAR roms model data from UCSC THREDDS server
    """
    vertical_stretching_function = ds['Vstretching'] # vertical terrain following stretching function
    vertical_transformation_function = ds['Vtransform'] # vertical terrain following tranformation equation
    h = ds['h'] # bathymetry at rho points
    h_critical = ds['hc'] # S-coordinate parameter, critical depth
    theta_bottom = ds['theta_b'].values # S-coordinate bottom control parameter
    theta_surface = ds['theta_s'].values # S-coordinate surface control parameter
    zeta = ds['zeta'] # free surface
    N_levels=42 # number of levels, note w has 43 levels
    zeta = np.squeeze(zeta[0,:,:]) # reduce dimensions to
    igrid=1 # density point grid, for (T,S) use igrid=3, for u and igrid=4 for v # use igrid=5 for w
    return  set_depth(vertical_transformation_function, vertical_stretching_function, theta_surface, theta_bottom, h_critical, N_levels, igrid, h, zeta, report=0)


def load_model_data():
    """
    Load model data from the UCSC THREDDS servers. Temperature data is subsetted and discrete depths are calcualted from the sigma coordinates. 

    Returns:
        xarray.DataArray: Model temperature.
    """
    #url = 'https://oceanmodeling.ucsc.edu:8443/thredds/dodsC/ccsra_2016a_phys_agg_slevs/fmrc/CCSRA_2016a_Phys_ROMS_Sigma-level_Aggregation_best.ncd'
    url = 'https://oceanmodeling.ucsc.edu:8443/thredds/dodsC/ccsra_2016a_phys_agg_slevs/fmrc/CCSRA_2016a_Phys_ROMS_Sigma-level_Aggregation_best.ncd'
    roms_ds = xr.open_dataset(url)
    z_rho = get_depth(roms_ds)
    roms_temp = roms_ds['temp']
    roms_temp_latitude = roms_temp['lat_rho'].values[:,0]
    roms_temp_longitude = roms_temp['lon_rho'].values[0,:]
    roms_temp = roms_temp.assign_coords({"eta_rho":roms_temp_latitude, "xi_rho":roms_temp_longitude})
    roms_temp = roms_temp.rename({"eta_rho":"latitude", "xi_rho":"longitude"})
    roms_temp = roms_temp.drop('lat_rho', errors='ignore')
    roms_temp = roms_temp.drop('lon_rho', errors='ignore')
    roms_temp = roms_temp.drop('time_run', errors='ignore')
    roms_temp = roms_temp.transpose('time','latitude','longitude','s_rho')    
    roms_temp['z_rho'] = (('latitude','longitude','s_rho'), z_rho) # Add the depth values back to the Dataset    
    return roms_temp


def subset_dataset(ds, bbox):
    """
    Subset the model data spatially within the bounds of the bounding box of the EEZ shapefile.
    Args:
        ds (xarray.DataArray): Model temperature data
        bbox (tuple): Bounding box with coordinates of the range of the EEZ

    Returns:
        xarray.DataArray: spatially subsetting model dataset.
    """
    return ds.sel(latitude=slice(bbox[0][1], bbox[0][0]), longitude=slice(bbox[1][1], bbox[1][0]))


def format_output(ds, heat_content_data):
    """Create and format a new xarray DataData with the calcuated heatcontent data. This data will be used to append the new original file. Writes that Dataset to a netcdf file.

    Args:
        ds (xarray.DataArray): Model temperature data used for the calculation
        heat_content_data (np.array): heat content calculation array. 

    """
    time = ds.time.values
    longitude = ds.longitude.values
    latitude = ds.latitude.values
    dims = ['time', 'lat', 'lon']
    ds_output = xr.Dataset({'heat_content_100_meters': (dims, heat_content_data)},
                    coords={"time":time,
                            "lon": longitude,
                            "lat": latitude,
                            })
    ds_output.attrs['title'] = "West Coast Upper Ocean Heat Content - 0-100 meters"
    ds_output.attrs['notes'] = "Created on "+dt.datetime.today().strftime("%Y-%m-%d") + " by pdaniel. Data is generated daily on the concave.shore.mbari."
    fname = "/home/flbahr/heat_content/data/West_Coast_wide_upper_heat_content_recent.nc"
    ds_output.to_netcdf(path=fname)
 
    
def append_datasets():
    """
    Open heatcontent file, append with newly generated data and overwrite the original file.
    """
    fname_recent = "/home/flbahr/heat_content/data/West_Coast_wide_upper_heat_content_recent.nc"
    fname_original = FILENAME
    # THis is a little hacky. 
    try:
        ds = xr.open_mfdataset([fname_original,fname_recent])
        if os.path.exists(fname_original):
            os.remove(fname_original)

    except Exception as e:
        logging.debug("Failed to open both new and old data. Error: {}".format(e))
    try:
        # Overwrite original
        ds.to_netcdf(fname_original)
        logging.info("Overwriting data file.")
    except Exception as e:
        logging.debug("Failed to write new netCDF. Error: {}".format(e))


##def copy_file_to_webserver(fname):
##    """Copy images from model runs to webserver where they can be viewed publically."""
##    try:
##        os.system('scp -i /etc/ssh/keys/pdaniel/scp_rsa {}  spyglass.mbari.org:/opt/apache-tomcat/content/thredds/public/data/models'.format(fname))
##        logging.info('MOVING: {} to spyglass-thredds'.format(fname))
##    except:
##        logging.debug('Unabled to move {} to spyglass'.format(fname))


def estimate_Heat_Content_WestCoast():
    isStale, times = get_last_time()
    # what do we do if times is None?
    if times==None:
        x=dt.datetime(2022,1,1)
        y=dt.datetime(2022,4,25)
        times=[x,y]
    if isStale:
        bbox = load_eez_shapefile()
        ds = load_model_data()
        #pdb.set_trace()
        ds = ds.sel(time=slice(times[0],times[1]))
        ds = subset_dataset(ds, bbox)
        westcoast_heat_content = np.zeros(shape=ds[:,:,:,0].shape)
        #pdb.set_trace()
        for i in tqdm(range(len(westcoast_heat_content))):
            upper_heat_content = calculate_upper_heat_content(ds.isel(time=i), upper_depth=-100)
            westcoast_heat_content[i,:,:] = upper_heat_content
        #pdb.set_trace()
        format_output(ds, westcoast_heat_content)
        append_datasets()
        #copy_file_to_webserver(FILENAME)

    else:
        logging.info("Model is up-to-date")


if __name__ == "__main__":
    estimate_Heat_Content_WestCoast()
