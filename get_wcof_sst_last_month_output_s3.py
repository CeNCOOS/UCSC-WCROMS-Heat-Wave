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
#
#
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
last_month=current_month-1
current_year=atime.year
# look at edge case of current time is january 1
# then last_month would be 0
if last_month==0:
    last_month=12
    current_year=atime.year-1
#
# will this work for getting january so current_month=2, last_month=1, current_year is correct so no need to change
#
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
# check if it is a leap year only if last_month==2
#
if last_month==2:
    isleap=current_year%4
    if isleap==0:
        dayspermonth[1]=29
#
string_year=str(current_year)
if last_month > 9:
    string_month=str(last_month)
else:
    string_month='0'+str(last_month)
#
# Need to set up the edge cases of January/December?
#
#
# define a loop in time
#
thedays=np.arange(dayspermonth[last_month-1]-1,-1,-1)
#thedays=np.arange(30,-1,-1)
io=0
tic=time.time()
#
atime=dt.datetime(current_year,last_month,dayspermonth[last_month-1])
for offset in thedays:
    print(offset)
    try:
        [timearray,temparray]=wcof_load_s3(atime,offset,24,mx,my)
        if io==0:
            bigtime=timearray
            bigtemp=np.expand_dims(temparray,axis=0)
            io=1
        else:
            bigtime=np.append(bigtime,timearray)
            tmptmp=np.expand_dims(temparray,axis=0)
            bigtemp=np.vstack((bigtemp,tmptmp))
    except:
        print('failed for this date\n')
#        pass # failed to get data for this date
toc=time.time()
print(toc-tic)
pfile=open('/home/flbahr/heat_content/wcofs_'+monthstring[last_month-1]+'_'+string_year+'_20230605.p','wb')
pobj=[bigtime,bigtemp]
pickle.dump(pobj,pfile)
# Have to close pickle file for it to write EOF and thus be readable.  If we don't do this
# the file has the data but is unreadable and just takes up space
pfile.close()
