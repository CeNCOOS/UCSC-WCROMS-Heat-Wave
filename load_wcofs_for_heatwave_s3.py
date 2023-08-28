import xarray as xr
import numpy as np
import scipy.io
from scipy import integrate
from scipy.interpolate import griddata
import pandas as pd
import datetime as dt
import time, scipy.io
import s3fs
import pdb
#
fs=s3fs.S3FileSystem(anon=True)
#
# atime is last time in wcofs file
# 
# Need to get todays date to figure out when to look for forecast data
#
def wcof_load(atime,mx,my):
    ioerr=0
    temparray=[]
    timearray=[]
    # grab current time
    current=dt.datetime.now()
    # get last time in wcofs file
    wlast=pd.to_datetime(atime[-1])
    # compute number of days to load of wcofs data from co-ops
    timedelta=current-wlast
    #
    offset=timedelta.days
    daystoload=np.arange(offset,-1,-1)
    
    io=0
    # this goes upto 3 days behind current time
    urlpre='s3://noaa-nos-ofs-pds/wcofs/netcdf/'    
    #urlpre1='https://ncei.noaa.gov/thredds/dodsC/model-wcofs-files/'
    # this contains forecast but doesn't go too far back in time
    #urlpre2='https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/WCOFS/MODELS/'
    #urlpre3='/home/flbahr/heat_content/data/'
    filestr='nos.wcofs.2ds.'
    filepost='.t03z.nc'
    # we have an issue that nos.wcofs.2ds.n001.20220913.t03z.nc is not at time 001 but at time 2022-09-13T04:00:00
    #
    for daystoget in daystoload:
        # compute the date based upon the offset
        thedate=current-dt.timedelta(days=int(daystoget))
        # get the components
        year=thedate.year
        month=thedate.month
        day=thedate.day
        #if year==2023 and month==2 and day==27:
        #    urlpre=urlpre3
        #if year==2023 and month==2 and day==28:
        #    urlpre=urlpre3
        # create an datetime to get
        datetoget=dt.datetime(year,month,day)
        # since files start on hour 3 we need to go back 3 hours from the 00 file
        astart_date=datetoget-dt.timedelta(hours=3)
        thehours=np.arange(0,24,1)
        forns='n'
        for zhours in thehours:
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
            urlstr=urlpre+syear+smonth+'/'+filestr+forehour+'.'+syear+smonth+sday+filepost
            print(urlstr)
            print('the hours= '+str(zhours))
            #pdb.set_trace()
            # try and load the data file
            fileid=fs.open(urlstr)
            dataset=xr.open_dataset(fileid)
            fileid.close()
                #
            lat=dataset['lat_rho']
            lon=dataset['lon_rho']
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
            t=np.ravel(roms_temp)
            xl=np.ravel(lon)
            yl=np.ravel(lat)
            zgrid=griddata((xl,yl),t,(mx,my),method='linear')

            if io==0:
                bigtemp=np.expand_dims(zgrid,axis=0)
                bigtime=otime
                io=1
            else:
                bigtime=np.append(bigtime,otime)
                btmp=np.expand_dims(zgrid,axis=0)
                bigtemp=np.vstack((bigtemp,btmp))
# need a check to see if bigtemp exists
        meantemp1=np.mean(bigtemp,axis=0)
        meantime=pd.Series(bigtime).mean()
        meantemp1=np.expand_dims(meantemp1,axis=0)
        io=0
        if len(temparray)==0:
            temparray=meantemp1
        else:
            temparray=np.vstack((temparray,meantemp1))
        timearray=np.append(timearray,meantime)
    # now to get the forecast data...
    # something odd file time is listed  as 20221011 20th hour but ocean time is 2022-10-10T23:00:00 so almost a day off
    # not sure I understand what is going on here...
    thehours=np.arange(24,28,1)
    for zhours in thehours:
        timetoget=astart_date+dt.timedelta(hours=int(zhours))
        print(timetoget)
        year=timetoget.year
        month=timetoget.month
        day=timetoget.day
        hour=timetoget.hour
        if hour==0:
            timetoget=astart_date+dt.timedelta(hours=int(zhours-1))
            year=timetoget.year
            month=timetoget.month
            day=timetoget.day
            hour=24        
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
        urlstr=urlpre+syear+smonth+'/'+filestr+forehour+'.'+syear+smonth+sday+filepost
        #urlstr=urlpre+syear+'/'+smonth+'/'+sday+'/'+filestr+forehour+'.'+syear+smonth+sday+filepost
        print(urlstr)
        print('the hours 24 to 28 '+str(zhours))
#        try:
        fileid=fs.open(urlstr)
        dataset=xr.open_dataset(fileid)
        fileid.close()
#            ioerr=0
#        except HTTPError as err:
#            if (err.code==502) or (err.code==503):
#                try:
#                    dataset=xr.open_dataset(urlstr)
#                    ioerr=0
#                except Exception as e:
#                    print('Exception '+str(e))
#                    ioerr=1
#            else:
#                print('Some other error occurred')
#                ioerr=1
#            #
        if ioerr==0:
            lat=dataset['lat_rho']
            lon=dataset['lon_rho']
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
            t=np.ravel(roms_temp)
            xl=np.ravel(lon)
            yl=np.ravel(lat)
            zgrid=griddata((xl,yl),t,(mx,my),method='linear')
            if io==0:
                bigtemp=np.expand_dims(zgrid,axis=0)
                bigtime=otime
                io=1
            else:
                bigtime=np.append(bigtime,otime)
                btmp=np.expand_dims(zgrid,axis=0)
                bigtemp=np.vstack((bigtemp,btmp))
#        except:
#            # deal with missing data files
#            print(urlstr+' is missing from set')
#    forec=np.arange(27,48)
    #forec=np.arange(1,69)
    forec=np.arange(1,21)
    forns='f'
#    pdb.set_trace()
    for zhours in forec:
        forns='f'
        #print(zhours)
        timetoget=astart_date+dt.timedelta(hours=int(26))
        #pdb.set_trace()
        print(timetoget)
        year=timetoget.year
        month=timetoget.month
        day=timetoget.day
        hour=zhours
#        hour=timetoget.hour
        if hour==0:
            timetoget=astart_date+dt.timedelta(hours=int(zhours-1))
            year=timetoget.year
            month=timetoget.month
            day=timetoget.day
            hour=24
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
        urlstr=urlpre+syear+smonth+'/'+filestr+forehour+'.'+syear+smonth+sday+filepost
        #urlstr=urlpre+syear+'/'+smonth+'/'+sday+'/'+filestr+forehour+'.'+syear+smonth+sday+filepost
        print(urlstr)
        print('The hours 1 to 21 '+str(zhours))
#            try:
        fileid=fs.open(urlstr)
        dataset=xr.open_dataset(fileid)
        fileid.close()
#                ioerr=0
#            except HTTPError as err:
#                if (err.code==502) or (err.code==503):
#                    try:
#                        dataset=xr.open_dataset(urlstr)
#                        ioerr=0
#                    except Exception as e:
#                        print('Exception '+str(e))
#                else:
#                    print('Some other error occurred')
#                    ioerr=1
        #
        lat=dataset['lat_rho']
        lon=dataset['lon_rho']
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
        t=np.ravel(roms_temp)
        xl=np.ravel(lon)
        yl=np.ravel(lat)
        zgrid=griddata((xl,yl),t,(mx,my),method='linear')
        if io==0:
            bigtemp=np.expand_dims(zgrid,axis=0)
            bigtime=otime
            io=1
        else:
            bigtime=np.append(bigtime,otime)
            btmp=np.expand_dims(zgrid,axis=0)
            bigtemp=np.vstack((bigtemp,btmp))
#        except:
#            # deal with missing data files
#            print(urlstr+' is missing from set')
#            #pdb.set_trace()
    meantemp1=np.mean(bigtemp,axis=0)
    meantime=pd.Series(bigtime).mean()
    meantemp1=np.expand_dims(meantemp1,axis=0)
    #pdb.set_trace()
    io=0
    if len(temparray)==0:
        temparray=meantemp1
    else:
        temparray=np.vstack((temparray,meantemp1))
    timearray=np.append(timearray,meantime)
    #
    #pdb.set_trace()
    forec=np.arange(21,45,1)
    forns='f'
    for zhours in forec:
        #print(zhours)
        timetoget=astart_date+dt.timedelta(hours=int(26))
        year=timetoget.year
        month=timetoget.month
        day=timetoget.day
        hour=zhours
        #hour=timetoget.hour
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
        urlstr=urlpre+syear+smonth+'/'+filestr+forehour+'.'+syear+smonth+sday+filepost
#        urlstr=urlpre+syear+'/'+smonth+'/'+sday+'/'+filestr+forehour+'.'+syear+smonth+sday+filepost
        print(urlstr)
        print('The hours 21 to 45 '+str(zhours))
#        try:
        fileid=fs.open(urlstr)
        dataset=xr.open_dataset(fileid)
        fileid.close()
#            ioerr=0
#        except HTTPError as err:
#            if (err.code==502) or (err.code==503):
#                try:
#                    dataset=xr.open_dataset(urlstr)
#                    ioerr=0
#                except Exception as e:
#                    print('Exception '+str(e))
#                    ioerr=1
#            else:
#                print('Some other error occurred')
#                ioerr=1
#        if ioerr==0:
        #
        lat=dataset['lat_rho']
        lon=dataset['lon_rho']
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
        t=np.ravel(roms_temp)
        xl=np.ravel(lon)
        yl=np.ravel(lat)
        zgrid=griddata((xl,yl),t,(mx,my),method='linear')
        if io==0:
            bigtemp=np.expand_dims(zgrid,axis=0)
            bigtime=otime
            io=1
        else:
            bigtime=np.append(bigtime,otime)
            btmp=np.expand_dims(zgrid,axis=0)
            bigtemp=np.vstack((bigtemp,btmp))
#       except:
 #           # deal with missing data files
 #           print(urlstr+' is missing from set')
    meantemp1=np.mean(bigtemp,axis=0)
    meantime=pd.Series(bigtime).mean()
    meantemp1=np.expand_dims(meantemp1,axis=0)
    io=0
    if len(temparray)==0:
        temparray=meantemp1
    else:
        temparray=np.vstack((temparray,meantemp1))
    timearray=np.append(timearray,meantime)
    #pdb.set_trace()
    forec=np.arange(45,69,1)
    forns='f'
    for zhours in forec:
        timetoget=astart_date+dt.timedelta(hours=int(26))
        year=timetoget.year
        month=timetoget.month
        day=timetoget.day
        hour=zhours
        #hour=timetoget.hour
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
        urlstr=urlpre+syear+smonth+'/'+filestr+forehour+'.'+syear+smonth+sday+filepost
#        urlstr=urlpre+syear+'/'+smonth+'/'+sday+'/'+filestr+forehour+'.'+syear+smonth+sday+filepost
        print(urlstr)
#        try:
        fileid=fs.open(urlstr)
        dataset=xr.open_dataset(fileid)
        fileid.close()
#            ioerr=0
#        except HTTPError as err:
#            if (err.code==502) or (err.code==503):
#                try:
#                    dataset=xr.open_dataset(urlstr)
#                    ioerr=0
#                except Exception as e:
#                    print('Exception '+str(e))
#                    ioerr=1
#            else:
#                print('Some other error occurred')
#                ioerr=1
#        if ioerr==0:
        #
        lat=dataset['lat_rho']
        lon=dataset['lon_rho']
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
        t=np.ravel(roms_temp)
        xl=np.ravel(lon)
        yl=np.ravel(lat)
        zgrid=griddata((xl,yl),t,(mx,my),method='linear')
        if io==0:
            bigtemp=np.expand_dims(zgrid,axis=0)
            bigtime=otime
            io=1
        else:
            bigtime=np.append(bigtime,otime)
            btmp=np.expand_dims(zgrid,axis=0)
            bigtemp=np.vstack((bigtemp,btmp))
#        except:
#            # deal with missing data files
#            print(urlstr+' is missing from set')
    meantemp1=np.mean(bigtemp,axis=0)
    meantime=pd.Series(bigtime).mean()
    meantemp1=np.expand_dims(meantemp1,axis=0)
    io=0
    if len(temparray)==0:
        temparray=meantemp1
    else:
        temparray=np.vstack((temparray,meantemp1))
    timearray=np.append(timearray,meantime)

    # need to append to array?
    return [timearray, temparray]


        
##    thedate=atime-dt.timedelta(days=int(offset))
##    # compute the parts so we can construct the url and filename
##    year=thedate.year
##    month=thedate.month
##    day=thedate.day
##    print(str(year)+'/'+str(month)+'/'+str(day))
##    # okay now we have year month and day of our offset date so we want to find the hour 1 to hour 24 values
##    datetoget=dt.datetime(year,month,day)
##    # so datetoget is the date we want to get but to get the actual hours for that day we need to get 3 hours earlier
##    astart_date=datetoget-dt.timedelta(hours=3)
##    # set up appropiate loop?
##    if numbs > 24:
##        thehours=np.arange(0,72,1)
##        forns='f'
##    else:
##        thehours=np.arange(0,24,1)
##        forns='n'
##    
##    for zhours in thehours:
##        # compute offset from start
##        timetoget=astart_date+dt.timedelta(hours=int(zhours))
##        year=timetoget.year
##        month=timetoget.month
##        day=timetoget.day
##        hour=timetoget.hour
##        # case to get 24 hr values and not have n000 or f000 since we have n024 and f024
##        if hour==0:
##            timetoget=astart_date+dt.timedelta(hours=int(zhours-1))
##            year=timetoget.year
##            month=timetoget.month
##            day=timetoget.day
##            hour=24
##        # end if
##            
##        syear=str(year)
##        if month < 10:
##            smonth='0'+str(month)
##        else:
##            smonth=str(month)
##        if day < 10:
##            sday='0'+str(day)
##        else:
##            sday=str(day)
##        if hour < 10:
##            forehour=forns+'00'+str(hour)
##        else:
##            forehour=forns+'0'+str(hour)
### This is for the NCEI site
##        urlstr=urlpre+syear+'/'+smonth+'/'+filestr+forehour+'.'+syear+smonth+sday+filepost
### this works for CO-OPS site
###        urlstr=urlpre+syear+'/'+smonth+'/'+sday+'/'+filestr+forehour+'.'+syear+smonth+sday+filepost
##        print(urlstr)
##    # ocean time is relative to 2016-01-01
##    # It turns out we can have missing nowcast files at ncei so we need to trap for that
##        try:
##            dataset=xr.open_dataset(urlstr)
##    #dataset=xr.open_dataset("https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/WCOFS/MODELS/2022/08/23/nos.wcofs.fields.n024.20220823.t03z.nc")
##    #dataset=xr.open_dataset("https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/WCOFS/MODELS/2022/08/23/nos.wcofs.2ds.n024.20220823.t03z.nc")
##    #dataset=open_url("https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/WCOFS/MODELS/2022/08/23/nos.wcofs.fields.n024.20220823.t03z.nc")
###
##            lat=dataset['lat_rho']
##            lon=dataset['lon_rho']
###
##            otime=dataset['ocean_time']
##            #print(urlstr)
##            #print(otime)
##            roms_temp = dataset['temp_sur']
##            roms_time = dataset['ocean_time']
##            roms_temp_latitude = roms_temp['lat_rho'].values[:,0]
##            roms_temp_longitude = roms_temp['lon_rho'].values[0,:]
##            roms_temp = roms_temp.assign_coords({"eta_rho":roms_temp_latitude, "xi_rho":roms_temp_longitude})
##            roms_temp = roms_temp.rename({"eta_rho":"latitude", "xi_rho":"longitude"})
##            roms_temp = roms_temp.drop('lat_rho', errors='ignore')
##            roms_temp = roms_temp.drop('lon_rho', errors='ignore')
##            roms_temp = roms_temp.drop('time_run', errors='ignore')
##            roms_temp = roms_temp.transpose('ocean_time','latitude','longitude')
##            roms_temp=np.squeeze(roms_temp)
###
##            t=np.ravel(roms_temp)
##            xl=np.ravel(lon)
##            yl=np.ravel(lat)
##            zgrid=griddata((xl,yl),t,(mx,my),method='linear')
##            if io==0:
##                #testtemp=np.expand_dims(roms_temp,axis=0)
##                bigtemp=np.expand_dims(zgrid,axis=0)
##                bigtime=otime
##                io=1
##            else:
##                bigtime=np.append(bigtime,otime)
##                btmp=np.expand_dims(zgrid,axis=0)
##                #ttmp=np.expand_dims(roms_temp,axis=0)
##                bigtemp=np.vstack((bigtemp,btmp))
##                #testtemp=np.vstack((testtemp,ttmp))
##        except:
##            print(urlstr+' is missing from set')
##    # average in time
##    meantemp1=np.mean(bigtemp,axis=0)
##    meantime=pd.Series(bigtime).mean()
##    # No need to keep big array and average as differnce between interpolated data and avraged big array that is
##    # downsampled is in the 1e-6 range.  So in the noise of the computation.
##    #meantemp2=np.mean(testtemp,axis=0)
##    #tt=np.ravel(meantemp2)
##    #qgrid=griddata((xl,yl),tt,(mx,my),method='linear')
##    #pdb.set_trace()
##    return [meantime, meantemp1]
