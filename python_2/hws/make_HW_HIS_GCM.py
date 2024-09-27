#!/usr/bin/env python
# coding: utf-8

from netCDF4 import Dataset
from netCDF4 import num2date
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import pandas as pd
import xarray as xr
import scipy
from scipy import special
from scipy.stats import norm
from scipy.stats import ks_2samp as ks_2samp
import xesmf as xe
import sys
import cdo
from cdo import *   # python version
cdo = Cdo()
import os
import xesmf as xe 


root_dir = '/work3/lorenzosangelantoni/scratch75/gcm_driven_experiment_from_nird/'

YSTART = 1996
YEND = 2006

boundaries = np.array([-30,30,30,60])
boundaries = np.array([-20,40,20,80]) ; lon_min,lon_max,lat_min,lat_max = boundaries

# Choose if intrpolating GCMs into a common grid 
interp = 1 

M = int(sys.argv[1])
#M=0

MODELS = [ 'CNRM-CM5_r1i1p1', 'EC-EARTH_r12i1p1', 'HadGEM2-ES_r1i1p1', 'MPI-ESM-LR_r1i1p1', 'NorESM1-M_r1i1p1']
RES='GCM'

DIR_GCM     = f'{root_dir}/GCMs/{MODELS[M]}/surface'
DIR_GCM_out = f'{root_dir}/GCMs/{MODELS[0]}/surface'
DIR_GCM_p90 = f'{root_dir}/GCMs/{MODELS[M]}/surface/p90'

filename = [ \
f'tasmax_day_CNRM-CM5_historical_r1i1p1_box_1996_2005.nc', \
f'tasmax_day_EC-EARTH_historical_r12i1p1_1996_2005_box.nc', \
f'tasmax_day_HadGEM2-ES_historical_r1i1p1_box_1996_2005.nc', \
f'tasmax_day_MPI-ESM-LR_historical_r1i1p1_box_1996_2005.nc', \
f'tasmax_day_NorESM1-M_historical_r1i1p1_box_1996_2005.nc', \
]

filename_p90 = [ \
f'tasmax_day_CNRM-CM5_historical_r1i1p1_box_1996_2005.nc', \
f'tasmax_day_EC-EARTH_historical_r12i1p1_1996_2005_box.nc', \
f'tasmax_day_HadGEM2-ES_historical_r1i1p1_box_1996_2005.nc', \
f'tasmax_day_MPI-ESM-LR_historical_r1i1p1_box_1996_2005.nc', \
f'tasmax_day_NorESM1-M_historical_r1i1p1_box_1996_2005.nc', \
]

filename_LSM = [ f'sftlf_fx_CNRM-CM5_historical_r0i0p0_box.nc', \
f'sftlf_fx_CNRM-CM5_historical_r0i0p0.nc', \
f'sftlf_fx_HadGEM2-ES_historical_r0i0p0.nc', \
f'sftlf_fx_MPI-ESM-LR_historical_r0i0p0.nc', \
f'sftlf_fx_NorESM1-M_historical_r0i0p0.nc', \
]

LSM_GCM = f'{root_dir}/GCMs/LSM/{filename_LSM[0]}'

print(f'MODELS: {MODELS[M]}')
print(f'FILE LSM: {LSM_GCM}')

tasmax_ifile     = f'tasmax_{MODELS[M]}_{RES}.nc'
tasmax_p90_ifile = f'tasmax_p90_{MODELS[M]}_{RES}.nc' 
lsmask_ifile     = f'lsmask_{MODELS[M]}_{RES}.nc'

#os.remove(f'/scratch/lorenzosangelantoni/gcm_driven_experiment_from_nird/scripts/python_2/hws/outputs/{filename[M]}') 

cdo.sellonlatbox(lon_min,lon_max,lat_min,lat_max, input=f'{DIR_GCM}/{filename[M]}', output=tasmax_ifile)
cdo.sellonlatbox(lon_min,lon_max,lat_min,lat_max, input=f'{DIR_GCM_p90}/{filename_p90[M]}', output=tasmax_p90_ifile)
cdo.sellonlatbox(lon_min,lon_max,lat_min,lat_max, input=f'{LSM_GCM}', output=lsmask_ifile)

#f_out = f'./outputs/{filename[M]}'
f_out = f'/work3/lorenzosangelantoni/scratch75/gcm_driven_experiment_from_nird/scripts/python_2/hws/outputs_post_prague//historical/GCM/{filename[M]}'

# tasmax 
dataset = Dataset(f'{tasmax_ifile}', mode='r')
lon = dataset.variables['lon'][:]
lat = dataset.variables['lat'][:]
if np.min(lon)>=0 : 
    lon_lsm = lon - 180  
tim = dataset.variables['time'][:]
date = num2date(tim[:],dataset['time'].units)
units = dataset['time'].units
tasmax = dataset.variables['tasmax'][:]
dataset.close()
print(f'SIZE TASMAX:  {tasmax.shape}')

# tasmax p90 
dataset = Dataset(f'{tasmax_p90_ifile}', mode='r')
tasmax_p90 = dataset.variables['tasmax'][:]
dataset.close()


# land-sea mask    
# numpy
dataset = Dataset(f'{LSM_GCM}', mode='r')
lon_lsm = dataset.variables['lon'][:]
lat_lsm = dataset.variables['lat'][:]
lsm=dataset.variables['sftlf'][:]
lsm.fill_value=0
if np.min(lon_lsm)>=0 : 
    lon_lsm = lon_lsm - 180  
dataset.close()
MASK = lsm <=50 
# xarray
ds_mask = xr.open_dataset(f'{LSM_GCM}').drop_dims('bnds')
da_mask = ds_mask.sftlf


# Define destination grid for the HW outputs
dataset = Dataset(f'{DIR_GCM_out}/{filename[0]}', mode='r')
lon_out = dataset.variables['lon'][:]
lat_out = dataset.variables['lat'][:]
if np.min(lon_out)>=0 : 
    lon_lsm = lon - 180  
dataset.close()

if interp > 0 : 
    ntim = 10 #tasmax.shape[0]
    nlat = lat_out.shape[0]
    nlon = lon_out.shape[0]
else : 
    ntim = 10 #tasmax.shape[0]
    nlat = lat.shape[0]
    nlon = lon.shape[0]

YEARS=np.arange(YSTART,YEND,1)

class structtype():
    pass

HW_persistence       = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_intensity         = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_index_event_start = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_index_event_end   = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_HWMI              = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_mean_tmax         = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_max_tmax          = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_number            = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan


cc = -1 
for Y in range(YSTART,YEND,1) : 
    
    print("MODEL:",MODELS[M])
    print("YEARS: ",Y)

    cc = cc + 1 
    # Count heatwave  
    #tasmax_y  = cdo.selyear(Y,input=f'{tasmax_ifile}',returnCdf=True).variables['tasmax'][:] 
    tasmax_y  = cdo.selyear(Y,input='-selmonth,6/8 '+f'{tasmax_ifile}',returnCdf=True).variables['tasmax'][:]     
    for ii in range(0,tasmax.shape[1]) :  #(30,31) : 
    
        #print("Current progress: " + str((ii)*100/len(lon)) + "%")
        
        for jj in range(0,tasmax.shape[2]) : #(20,21) : 
            
            if np.isnan(tasmax[0,ii,jj]) == False :
                
                pacchetti         = [] #np.empty(92)*np.nan
                dummy_index       = [] #np.empty(92)*np.nan
                index_event_start = [] #np.empty(92)*np.nan
                index_event_end   = [] #np.empty(92)*np.nan
                index_events      = [] # np.empty(92) * np.nan # It considers all the HWs in one year
                HWMI              = [] #np.empty(92)*np.nan 
                hfls_deficit      = [] #np.empty(92)*np.nan
            
            
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
os.remove(lsmask_ifile)


# Interp to a common GCM grid (CNRM) the HW metrics 
grid_in = {"lon": lon, "lat": lat}
grid_out = {"lon": lon_out, "lat": lat_out}
regridder = xe.Regridder(grid_in, grid_out, "bilinear")


# *** define output file
out_nc1 = Dataset(f_out, "w", format="NETCDF4")

# *** define time to store in the netcdf file 
cdo.yearmean(input = f'{root_dir}/GCMs/{MODELS[0]}/surface/{filename[0]}',output = 'ofile.nc',options = '-f nc') # Attention to the time_period from one model only
dataset = Dataset('ofile.nc', mode='r')
time_period = dataset.variables['time'][:]
date = num2date(time_period,dataset['time'].units)
units = dataset['time'].units
dataset.close()

# **** define dimensions
lats = out_nc1.createDimension('y', nlat)
lons = out_nc1.createDimension('x', nlon)
time = out_nc1.createDimension('time',ntim)
model_dim = out_nc1.createDimension('model', 1)

# **** define varaibles
latitudes = out_nc1.createVariable('lat', 'f4', 'y')
longitudes= out_nc1.createVariable('lon', 'f4', 'x')
times     = out_nc1.createVariable('time', 'f8', ('time',))
v1 = out_nc1.createVariable('HW_persistence', 'f4', ('time', 'y', 'x'))
v2 = out_nc1.createVariable('HW_HWMI',        'f4', ('time', 'y', 'x'))
v3 = out_nc1.createVariable('HW_index_event_start','f4', ('time', 'y', 'x'))
v4 = out_nc1.createVariable('HW_index_event_end',  'f4', ('time', 'y', 'x'))
v5 = out_nc1.createVariable('HW_mean_tmax',        'f4', ('time', 'y', 'x'))
v6 = out_nc1.createVariable('HW_max_tmax',         'f4', ('time', 'y', 'x')) 
v7 = out_nc1.createVariable('HW_number',           'f4', ('time', 'y', 'x'))
model_var = out_nc1.createVariable('model', 'S10', ('model',))

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
if interp > 0 : 
    latitudes[:]  = lat_out
    longitudes[:] = lon_out
    times[:] = time_period.data[:]  
else : 
    latitudes[:] = lat
    longitudes[:] = lon
    times[:] = time_period.data[:]  


v1[:,:,:] = regridder(HW_persistence)
v2[:,:,:] = regridder(HW_HWMI)
v3[:,:,:] = regridder(HW_index_event_start)
v4[:,:,:] = regridder(HW_index_event_end)
v5[:,:,:] = regridder(HW_mean_tmax)
v6[:,:,:] = regridder(HW_max_tmax)
v7[:,:,:] = regridder(HW_number)
model_var[0] = np.string_(MODELS[M])


out_nc1.close()

dataset = Dataset(f'{f_out}', mode='r')
var = dataset.variables['HW_mean_tmax'][:]

"""
fig, axes = plt.subplots(ncols=2,nrows=2, figsize=(12,7.2),
                             subplot_kw={'projection': ccrs.PlateCarree()})
    
    
P1 = axes[0,0].pcolormesh(
            lon_out,
            lat_out,
            np.nanmean(var,axis=0),vmin=290, vmax=np.nanmax(var),
            transform=ccrs.PlateCarree())    

interval=5
levels = np.arange(np.nanmin(var), np.nanmax(var) + interval, interval)
axes[0,0].coastlines(linewidth = 1)
axes[0, 0].set_extent([-10, 30, 35, 66], crs=ccrs.PlateCarree())
cb = fig.colorbar(P1, ax=(axes[0,0]), orientation='vertical',
                      aspect=15,shrink=.95,pad=.015, ticks=levels)
plt.show()      
""" 





