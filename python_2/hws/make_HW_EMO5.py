from netCDF4 import Dataset
from netCDF4 import num2date
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import pandas as pd
import xarray as xr
import scipy
from scipy import special
from scipy.stats import norm
from scipy.stats import ks_2samp as ks_2samp
import sys
import cdo
from cdo import *   # python version
cdo = Cdo()
import os

root_dir = '/work3/lorenzosangelantoni/scratch75/gcm_driven_experiment_from_nird' 

DIR_OBS     = f'{root_dir}/EMO-5'
DIR_OBS_p90 = f'{root_dir}/EMO-5/p90'
DIR_OUT     = f'{root_dir}/OUTPUT_HW'

PERIOD = 'historical'
VAR = 'tx'

#YSTART = int(sys.argv[1])
#YEND=YSTART+1

YSTART = 1996
YEND   = 2006

filename     = 'EMO-5-tx_1996_2005_JJA.nc'
filename_p90 = 'EMO-5-tx_1996_2005.nc'

f_out = f'/work3/lorenzosangelantoni/scratch75/gcm_driven_experiment_from_nird/scripts/python_2/hws/outputs_post_prague/historical/{filename}' 

# Domain mask
boundaries = np.array([1,17.0,40,50]) ; lon_min,lon_max,lat_min,lat_max = boundaries

tasmax_ifile=f'{filename}'
tasmax_p90_ifile=f'{filename_p90}' 

# Load variables 
cdo.sellonlatbox(lon_min,lon_max,lat_min,lat_max, input=f'{DIR_OBS}/{filename}', output=tasmax_ifile)
cdo.sellonlatbox(lon_min,lon_max,lat_min,lat_max, input=f'{DIR_OBS_p90}/{filename_p90}', output=tasmax_p90_ifile)

dataset = Dataset(f'{tasmax_ifile}', mode='r')
lon = dataset.variables['lon'][:]
lat = dataset.variables['lat'][:]
tim = dataset.variables['time'][:]
date = num2date(tim[:],dataset['time'].units)
units = dataset['time'].units
tasmax = dataset.variables['tx'][:]
dataset.close()

dataset = Dataset(f'{tasmax_p90_ifile}', mode='r')
tasmax_p90 = dataset.variables['tx'][:]
dataset.close()

YEARS=np.arange(YSTART,YEND,1)


ntim = 10 #tasmax.shape[0]
nlat = tasmax.shape[1]
nlon = tasmax.shape[2]

HW_persistence       = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_index_event_start = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_index_event_end   = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_HWMI              = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_mean_tmax         = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_max_tmax          = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_number            = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan


cc = -1 

for Y in range(YSTART,YEND,1) : 
    
    print(f'Y')
    
    cc = cc + 1 
    
    # Count heatwave  
    tasmax_y  = cdo.selyear(Y,input=f'{tasmax_ifile}',returnCdf=True).variables['tx'][:] 
    tasmax_y = tasmax_y.filled(np.nan)
    
    for ii in range(0,tasmax.shape[1]) :  #(30,31) : 
    
        print("Current progress: " + str((ii)*100/len(lon[1])) + "%")
        
        for jj in range(0,tasmax.shape[2]) : #(20,21) : 
            
            if np.isnan(tasmax[0,ii,jj]) == False :
                
                pacchetti         = [] 
                dummy_index       = [] 
                index_event_start = np.empty(92)*np.nan
                index_event_end   = np.empty(92)*np.nan
                index_events      = [] 
                HWMI              = np.empty(92)*np.nan 

                array  = tasmax_y[:,ii,jj]
                thresh = tasmax_p90[:,ii,jj]
                
                Count = -1     
                i = 0 
                cval = 3 

                IQR = scipy.stats.iqr(tasmax[:,ii,jj],rng=(25,75))
                p25 = np.percentile(tasmax[:,ii,jj],25)               
                
                
                # ------------------------------------------------
                while i < len(array)-1 : 


                    if array[i] >= thresh[i] : 

                        i = i + 1;

                        if i >= len(array) : 

                            break

                        c = 1;

                        while array[i] >= thresh[i] : # Start counting (c) time steps

                            c = c + 1
                            i = i + 1

                            if i >= len(array) :

                                break
                                
                        if c >= cval : # Start counting (Count) when "c" is above 3 consecutive days

                                Count = Count + 1
                                pacchetti.append(c)
                                dummy_index.append(i)
                                index_events.append(list(range(int(dummy_index[-1] - pacchetti[-1]), int(dummy_index[-1]))))

                    else :
                        i = i + 1

                # ------------------------------------------------
                
                if 'pacchetti' in locals() and len(pacchetti) > 0: 
                    if len(pacchetti) > 1 : 
                        #amplitude = [np.sum(array[index_events[p]] - thresh[index_events[p]]) for p in range(len(pacchetti))]
                        
                        amplitude = np.array([np.sum((array[index_events[p]] - p25) / IQR)  for p in range(len(pacchetti))])
                        
                        # Select the most intense event in the i year
                        max_amplitude_index = np.argmax(amplitude)
                        index_most_intense_event = index_events[max_amplitude_index]
                        HW_persistence[cc,ii,jj] = np.nanmax(pacchetti) 
                        HW_index_event_start[cc,ii,jj] = np.min(index_most_intense_event) 
                        HW_index_event_end[cc,ii,jj]   = np.max(index_most_intense_event) 
                        HW_mean_tmax[cc,ii,jj] = np.nanmean(array[index_most_intense_event])  
                        HW_max_tmax[cc,ii,jj]  = np.nanmax(array[index_most_intense_event])  		
                        HW_number[cc,ii,jj]    = len(pacchetti) 
                        HW_HWMI[cc, ii, jj]         = np.sum((array[index_most_intense_event] - p25) / IQR) 
                        
                    else :
                        index_event = index_events
                        HW_persistence[cc,ii,jj] = np.array(pacchetti)
                        HW_index_event_start[cc,ii,jj] = np.min(index_event)
                        HW_index_event_end[cc,ii,jj]   = np.max(index_event)
                        HW_mean_tmax[cc,ii,jj] = np.nanmean(array[index_event])
                        HW_max_tmax[cc,ii,jj]  = np.nanmax(array[index_event])
                        HW_number[cc,ii,jj] = len(pacchetti)
                        HW_HWMI[cc, ii, jj] = np.sum((array[index_event] - p25) / IQR)
                        
                        
os.remove(tasmax_ifile)
os.remove(tasmax_p90_ifile)
           
# *** define output file
out_nc1 = Dataset(f_out, "w", format="NETCDF4")


# *** define time to store in the netcdf file 
cdo.yearmean(input = f'{DIR_OBS}/{filename}',    output = 'ofile.nc', options = '-f nc') 
dataset = Dataset('ofile.nc', mode='r')
time_period = dataset.variables['time'][:]
date = num2date(time_period,dataset['time'].units)
dataset.close()

# **** define dimensions
lats = out_nc1.createDimension('y', nlat)
lons = out_nc1.createDimension('x', nlon)
time = out_nc1.createDimension('time',ntim)

# **** define varaibles
latitudes = out_nc1.createVariable('lat', 'f4', ('y','x'))
longitudes= out_nc1.createVariable('lon', 'f4', ('y','x'))
times     = out_nc1.createVariable('time', 'f8', ('time',))
v1 = out_nc1.createVariable('HW_persistence', 'f4', ('time', 'y', 'x'))
v2 = out_nc1.createVariable('HW_HWMI',        'f4', ('time', 'y', 'x'))
v3 = out_nc1.createVariable('HW_index_event_start','f4', ('time', 'y', 'x'))
v4 = out_nc1.createVariable('HW_index_event_end',  'f4', ('time', 'y', 'x'))
v5 = out_nc1.createVariable('HW_mean_tmax',        'f4', ('time', 'y', 'x'))
v6 = out_nc1.createVariable('HW_max_tmax',         'f4', ('time', 'y', 'x')) 
v7 = out_nc1.createVariable('HW_number',           'f4', ('time', 'y', 'x'))

# **** define variables attribute
latitudes.long_name = 'latitude'
latitudes.units = 'degree_north'

longitudes.long_name = 'longitude'
longitudes.units = 'degree_east'

times.long_name = 'time'
times.units = units

v1.long_name = 'HW persistence'
v1.units = 'DAYS'
v2.long_name = 'Heatwave magnitude index daily'
v2.units = 'HWMId'
v3.long_name = 'HW start time index'
v3.units = 'summer day'
v4.long_name = 'HW end time index'
v4.units = 'summer day'
v5.long_name = 'HW days mean maximum temperature'
v5.units = '°C'
v6.long_name = 'HW days max maximum temperature'
v6.units = '°C'
v7.long_name = 'HW number'
v7.units = 'N.'

# **** assign values to variables
latitudes[:,:] = lat
longitudes[:,:] = lon
times[:] = time_period.data[:]  
v1[:,:,:] = HW_persistence
v2[:,:,:] = HW_HWMI
v3[:,:,:] = HW_index_event_start
v4[:,:,:] = HW_index_event_end
v5[:,:,:] = HW_mean_tmax
v6[:,:,:]= HW_max_tmax 
v7[:,:,:]= HW_number 

out_nc1.close()

# Open output
dataset = Dataset(f'{f_out}', mode='r')

ds = []
ds = xr.Dataset(

    data_vars=dict(

         HW_persistence      = (["time","x", "y"], dataset.variables['HW_persistence'][:]),
         HW_HWMI             = (["time","x", "y"], dataset.variables['HW_HWMI'][:]),
         HW_index_event_start= (["time","x", "y"], dataset.variables['HW_index_event_start'][:]),         
         HW_index_event_end  = (["time","x", "y"], dataset.variables['HW_index_event_end'][:]),
         HW_mean_tmax        = (["time","x", "y"], dataset.variables['HW_mean_tmax'][:]-273.16),
         HW_max_tmax         = (["time","x", "y"], dataset.variables['HW_max_tmax'][:]-273.16),
         HW_number           = (["time","x", "y"], dataset.variables['HW_number'][:]),
),

    coords=dict(
        lon=(["x", "y"], dataset.variables['lon'][:]),
        lat=(["x", "y"], dataset.variables['lat'][:]),
        time=num2date(dataset['time'][:],dataset['time'].units)  #dataset.variables['time'][:] date = num2date(tim[:],dataset['time'].units)

        ),
)

# Map 
# -------------------------------------------------------------------
#plt.clf()
#levels = np.linspace(10,40,30)
my_cmap = 'jet'
plt.figure(figsize=(12,4))
data_crs = ccrs.PlateCarree() # Define transform
ax = plt.axes(projection=ccrs.PlateCarree()) # Define projection 
ds.HW_mean_tmax.mean('time').plot(ax=ax, x='lon', y='lat',transform=ccrs.PlateCarree(),cmap=my_cmap) 
ax.coastlines(linewidth = 1)
ax.gridlines(linewidth = .5)
ax.set_facecolor('white')
#plt.show()
plt.savefig(f'./figures/EMO5.png')
# ---------------

