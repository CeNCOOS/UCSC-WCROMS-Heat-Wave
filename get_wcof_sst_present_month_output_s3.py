import xarray as xr
import numpy as np
import scipy.io
from scipy import integrate
from scipy.interpolate import griddata
import datetime as dt
import time, scipy.io
from wcof_load_s3 import wcof_load_s3
import pickle
import pdb
#
# This code is to try and deal with issues arrising from NOAA server errors.  I'm having too many to keep the code up and running.
# These errors are out of my control.  So I have to build code to handle these faults and try to make up for them.
#
# Weird Bug here, but make sure to run this first
#import salem
#from shapely import geometry
# grab some boundary information for possible use later?
#sanctuary_outline= salem.read_shapefile('/home/flbahr/heat_content/shapefiles/USMaritimeLimitsAndBoundariesSHP/USMaritimeLimitsNBoundaries.shp')
#sanctuary_outline.crs = 'epsg:4326'
#eez_shape = sanctuary_outline[(sanctuary_outline['REGION'] == "Pacific Coast") & (sanctuary_outline['EEZ'])]
#bbox = ((eez_shape['max_y'].max(),eez_shape['min_y'].min()),(eez_shape['max_x'].max(),eez_shape['min_x'].min()))
#geom = geometry.box(minx=bbox[1][1],maxx=bbox[1][0],miny=bbox[0][1],maxy=bbox[0][0])
#
# Get latitude and longitude that defines the model domain.  Since the THREDDS server at UCSC is down for the moment working to get at the code.
#
roms_ds=xr.open_dataset('/home/flbahr/heat_content/WCOFS_SST_2022.nc',decode_cf=True) #different size lat and lon, this may be what we need before truncating?
roms_temp_latitudeUCSC =roms_ds['lat']
roms_temp_longitudeUCSC=roms_ds['lon']
[mx,my]=np.meshgrid(roms_temp_longitudeUCSC,roms_temp_latitudeUCSC)
#
# get the current date so we can get the current month
#
atime=dt.datetime.now()
# current month
current_month=atime.month
# so since we are in the next month the previous month should all be there
current_year=atime.year
#
# book keeping on months and number of days per month in a standard year
dayspermonth=[31,28,31,30,31,30,31,31,30,31,30,31]
monthstring=['january',
             'february',
             'march',
             'april',
             'may',
             'june',
             'july',
             'august',
             'september',
             'october',
             'november',
             'december']
#
# check if it is a leap year
#
isleap=current_year%4
if isleap==0:
    dayspermonth[1]=29
#
string_year=str(current_year)
if current_month > 9:
    string_month=str(current_month)
else:
    string_month='0'+str(current_month)
#
# define a loop in time
#
thedays=np.arange(atime.day,-1,-1)
io=0
# define start time so we can see how long this takes to run
tic=time.time()
#
atime=dt.datetime(current_year,current_month,atime.day)
#
# see if we have the current month file?
# NOTE CODE BELOW IS INCOMPLETE 
try:
    apfile=open('/home/flbahr/heat_content/wcofs_'+monthstring[current_month-1]+'_'+string_year+'.p','rb')
    [bigtime,bigtemp]=pickle.load(apfile)
    apfile.close()
    lastday=bigtime[-1].day
    lastmonth=bigtime[-1].month
    lastyear=bigtime[-1].year
    lasttime=dt.datetime(lastyear,lastmonth,lastday)
    n=(atime-lasttime).days-1
    thedays=np.arange(n,-1,-1)
#   Need to check and make sure we don't duplicate data already in the file.  We don't want to append if we already have the point
#   this would really mess things up...
# Need the time
    io=1
except:
    thedays=np.arange(atime.day,-1,-1)
    io=0
#pdb.set_trace()
for offset in thedays:
    print(offset)
    try:
        [timearray,temparray]=wcof_load_s3(atime,offset,24,mx,my)
        if io==0:
            bigtime=timearray
            bigtemp=np.expand_dims(temparray,axis=0)
            io=1
            #pdb.set_trace()
        else:
            bigtime=np.append(bigtime,timearray)
            tmptmp=np.expand_dims(temparray,axis=0)
            bigtemp=np.vstack((bigtemp,tmptmp))
    except:
        print('failed for this date\n')
#        pass # failed to get data for this date
#pdb.set_trace()
toc=time.time()
print(toc-tic)
pfile=open('/home/flbahr/heat_content/wcofs_'+monthstring[current_month-1]+'_'+string_year+'.p','wb')
pobj=[bigtime,bigtemp]
pickle.dump(pobj,pfile)
# Have to close pickle file for it to write EOF and thus be readable.  If we don't do this
# the file has the data but is unreadable and just takes up space
pfile.close()
