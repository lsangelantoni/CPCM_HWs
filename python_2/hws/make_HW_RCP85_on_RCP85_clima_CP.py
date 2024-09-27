#!/usr/bin/env python
# coding: utf-8

import os
os.environ["PATH"] = "/users_home/cmcc/ls21622/.conda/envs/shell/bin:" + os.environ["PATH"]
from netCDF4 import Dataset
from netCDF4 import num2date
import numpy as np
import pandas as pd
import xarray as xr
import scipy
import sys
import cdo
from cdo import *   # python version
cdo = Cdo()
import os

root_dir = '/data/cmcc/ls21622/gcm_driven_experiment_from_nird'
DIR_nonCP     = f'{root_dir}/FILES_CP/tasmax/rcp85'
DIR_nonCP_hfls= f'{root_dir}/FILES_CP/hfls/rcp85'
DIR_nonCP_p90 = f'{root_dir}/FILES_CP/tasmax/rcp85/p90'
tmp_cdo_dir  ='/work/cmcc/ls21622/tmp'

PERIOD = 'rcp85'
MODELS = [ 'BCCR-AUTH','BTU','CMCC','CNRM','ETHZ','FZJ-IDL','HCLIM','ICTP','KIT','KNMI','UKMO','JLU'] 
RES='CP'

# Choose model through M variable 
M = int(sys.argv[1])

YSTART = 2090
YEND   = 2100

if MODELS[M] == "UKMO" :
    YSTART=2096
    YEND = 2106
    
filename = [ \
'tasmax_ALP-3_NorESM1-ME_rcp85_r1i1p1_AUTH-MC-WRF381D_fpsconv-x1n2-v1_day_JJA.nc',\
'tasmax_ALP-3_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_CLMcom-BTU-CCLM5-0-14_fpsconv-x2yn2-v1_day_JJA.nc',\
'tasmax_ALP-3_ICHEC-EC-EARTH_rcp85_r12i1p1_CLMcom-CMCC-CCLM5-0-9_x2yn2v1_day_JJA.nc',\
'tasmax_ALP-3_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_CNRM-AROME41t1_fpsconv-x2yn2-v1_day_JJA.nc',\
'tasmax_ALP-3_MPI_rcp85_r1i1p1_COSMO-pompa_5.0_2019.1_day_JJA.nc',\
'tasmax_ALP-3_SMHI-EC-EARTH_rcp85_r12i1p1_FZJ-IDL-WRF381CA_fpsconv-x1n2-v1_day_JJA.nc',\
'tasmax_ALP-3_ICHEC-EC-EARTH_rcp85_r12i1p1_HCLIMcom-HCLIM38-AROME_fpsconv-x2yn2-v1_day_JJA.nc',\
'tasmax_ALP-3_MOHC-HadGEM2-ES_rcp85_r1i1p1_ICTP-RegCM4-7_fpsconv-x2yn2-v1_day_JJA.nc',\
'tasmax_ALP-3_MPI-M-MPI-ESM-LR_rcp85_r1i1p1_CLMcom-KIT-CCLM5-0-15_fpsconv-x2yn2-v1_day_JJA.nc',\
'tasmax_ALP-3_KNMI-EC-EARTH_rcp85_r04i1p1_KNMI-HCLIM38h1-AROME_fpsconv-x2yn2-v1_day_JJA.nc',\
'tasmax_ALP-3_HadGEM3-GC3.1-N512_rcp85_r1i1p1_HadREM3-RA-UM10.1_fpsconv-x0n1-v1_day_JJA.nc',\
'tasmax_ALP-3_MPI-M-MPI-ESM-LR_rcp85_r1i1p1_CLMcom-JLU-CCLM5-0-15_fpsconv-x0n1-v1_day_JJA.nc',\
]

filename_p90 = [ \
'tasmax_ALP-3_NorESM1-ME_rcp85_r1i1p1_AUTH-MC-WRF381D_fpsconv-x1n2-v1_day.nc',\
'tasmax_ALP-3_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_CLMcom-BTU-CCLM5-0-14_fpsconv-x2yn2-v1_day.nc',\
'tasmax_ALP-3_ICHEC-EC-EARTH_rcp85_r12i1p1_CLMcom-CMCC-CCLM5-0-9_x2yn2v1_day.nc',\
'tasmax_ALP-3_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_CNRM-AROME41t1_fpsconv-x2yn2-v1_day.nc',\
'tasmax_ALP-3_MPI_rcp85_r1i1p1_COSMO-pompa_5.0_2019.1_day.nc',\
'tasmax_ALP-3_SMHI-EC-EARTH_rcp85_r12i1p1_FZJ-IDL-WRF381CA_fpsconv-x1n2-v1_day.nc',\
'tasmax_ALP-3_ICHEC-EC-EARTH_rcp85_r12i1p1_HCLIMcom-HCLIM38-AROME_fpsconv-x2yn2-v1_day.nc',\
'tasmax_ALP-3_MOHC-HadGEM2-ES_rcp85_r1i1p1_ICTP-RegCM4-7_fpsconv-x2yn2-v1_day.nc',\
'tasmax_ALP-3_MPI-M-MPI-ESM-LR_rcp85_r1i1p1_CLMcom-KIT-CCLM5-0-15_fpsconv-x2yn2-v1_day.nc',\
'tasmax_ALP-3_KNMI-EC-EARTH_rcp85_r04i1p1_KNMI-HCLIM38h1-AROME_fpsconv-x2yn2-v1_day.nc',\
'tasmax_ALP-3_HadGEM3-GC3.1-N512_rcp85_r1i1p1_HadREM3-RA-UM10.1_fpsconv-x0n1-v1_day.nc',\
'tasmax_ALP-3_MPI-M-MPI-ESM-LR_rcp85_r1i1p1_CLMcom-JLU-CCLM5-0-15_fpsconv-x0n1-v1_day.nc', \
]

filename_hfls = [ \
'hfls_ALP-3_NorESM1-ME_rcp85_r1i1p1_AUTH-MC-WRF381D_fpsconv-x1n2-v1_day_JJA.nc',\
'hfls_ALP-3_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_CLMcom-BTU-CCLM5-0-14_fpsconv-x2yn2-v1_day_JJA.nc',\
'hfls_ALP-3_ICHEC-EC-EARTH_rcp85_r12i1p1_CLMcom-CMCC-CCLM5-0-9_x2yn2v1_day_JJA.nc',\
'hfls_ALP-3_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_CNRM-AROME41t1_fpsconv-x2yn2-v1_day_JJA.nc',\
'hfls_ALP-3_MPI_rcp85_r1i1p1_COSMO-pompa_5.0_2019.1_day_JJA.nc',\
'hfls_ALP-3_SMHI-EC-EARTH_rcp85_r12i1p1_FZJ-IDL-WRF381CA_fpsconv-x1n2-v1_day_JJA.nc',\
'hfls_ALP-3_ICHEC-EC-EARTH_rcp85_r12i1p1_HCLIMcom-HCLIM38-AROME_fpsconv-x2yn2-v1_day_JJA.nc',\
'hfls_ALP-3_MOHC-HadGEM2-ES_rcp85_r1i1p1_ICTP-RegCM4-7_fpsconv-x2yn2-v1_day_JJA.nc',\
'hfls_ALP-3_MPI-M-MPI-ESM-LR_rcp85_r1i1p1_CLMcom-KIT-CCLM5-0-15_fpsconv-x2yn2-v1_day_JJA.nc',\
'hfls_ALP-3_KNMI-EC-EARTH_rcp85_r04i1p1_KNMI-HCLIM38h1-AROME_fpsconv-x2yn2-v1_day_JJA.nc',\
'hfls_ALP-3_HadGEM3-GC3.1-N512_rcp85_r1i1p1_HadREM3-RA-UM10.1_fpsconv-x0n1-v1_day_JJA.nc',\
'hfls_ALP-3_MPI-M-MPI-ESM-LR_rcp85_r1i1p1_CLMcom-JLU-CCLM5-0-15_fpsconv-x0n1-v1_day_JJA.nc',\
]

filename_tasmax_rcp85_out =  [ \
'tasmax_ALP-3_NorESM1-ME_rcp85_r1i1p1_AUTH-MC-WRF381D_fpsconv-x1n2-v1_day_on_rcp85_clima_JJA.nc',\
'tasmax_ALP-3_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_CLMcom-BTU-CCLM5-0-14_fpsconv-x2yn2-v1_day_on_rcp85_clima_JJA.nc',\
'tasmax_ALP-3_ICHEC-EC-EARTH_rcp85_r12i1p1_CLMcom-CMCC-CCLM5-0-9_x2yn2v1_day_on_rcp85_clima_JJA.nc',\
'tasmax_ALP-3_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_CNRM-AROME41t1_fpsconv-x2yn2-v1_day_on_rcp85_clima_JJA.nc',\
'tasmax_ALP-3_MPI_rcp85_r1i1p1_COSMO-pompa_5.0_2019.1_day_on_rcp85_clima_JJA.nc',\
'tasmax_ALP-3_SMHI-EC-EARTH_rcp85_r12i1p1_FZJ-IDL-WRF381CA_fpsconv-x1n2-v1_day_on_rcp85_clima_JJA.nc',\
'tasmax_ALP-3_ICHEC-EC-EARTH_rcp85_r12i1p1_HCLIMcom-HCLIM38-AROME_fpsconv-x2yn2-v1_day_on_rcp85_clima_JJA.nc',\
'tasmax_ALP-3_MOHC-HadGEM2-ES_rcp85_r1i1p1_ICTP-RegCM4-7_fpsconv-x2yn2-v1_on_rcp85_clima_JJA.nc',\
'tasmax_ALP-3_MPI-M-MPI-ESM-LR_rcp85_r1i1p1_CLMcom-KIT-CCLM5-0-15_fpsconv-x2yn2-v1_day_on_rcp85_clima_JJA.nc',\
'tasmax_ALP-3_KNMI-EC-EARTH_rcp85_r04i1p1_KNMI-HCLIM38h1-AROME_fpsconv-x2yn2-v1_day_on_rcp85_clima_JJA.nc',\
'tasmax_ALP-3_HadGEM3-GC3.1-N512_rcp85_r1i1p1_HadREM3-RA-UM10.1_fpsconv-x0n1-v1_day_on_rcp85_clima_JJA.nc',\
'tasmax_ALP-3_MPI-M-MPI-ESM-LR_rcp85_r1i1p1_CLMcom-JLU-CCLM5-0-15_fpsconv-x0n1-v1_day_on_rcp85_clima_JJA.nc',\
]

f_out = f'/work/{filename_tasmax_rcp85_out[M]}'

# Land - sea mask
mask_CP = f'{root_dir}/OROG_SFTLF_WRF/cordex_fps_domain/sftlf_ALP-3_ECMWF-ERAINT_evaluation_r1i1p1_UCAN-WRF381BI_x1n2_fx_on_BCCR_grid.nc'

# Domain mask
boundaries = np.array([1,17.0,40,50]) ; lon_min,lon_max,lat_min,lat_max = boundaries

tasmax_ifile     = f'{tmp_cdo_dir}/tasmax_{MODELS[M]}_{RES}.nc'
tasmax_p90_ifile = f'{tmp_cdo_dir}/tasmax_p90_{MODELS[M]}_{RES}.nc' 
lsmask_ifile     = f'{tmp_cdo_dir}/lsmask_{MODELS[M]}_{RES}.nc'
tasmax_y_ifile   = f'{tmp_cdo_dir}/tasmax_{MODELS[M]}_Y_{RES}.nc'
hfls_ifile       = f'{tmp_cdo_dir}/hfls_{MODELS[M]}_{RES}.nc'
hfls_y_ifile     = f'{tmp_cdo_dir}/hfls_{MODELS[M]}_Y_{RES}.nc'

# Load variables 
cdo.sellonlatbox(lon_min,lon_max,lat_min,lat_max, input=f'{DIR_nonCP}/{filename[M]}', output=tasmax_ifile)
cdo.sellonlatbox(lon_min,lon_max,lat_min,lat_max, input=f'{DIR_nonCP_p90}/{filename_p90[M]}', output=tasmax_p90_ifile)
cdo.sellonlatbox(lon_min,lon_max,lat_min,lat_max, input=f'{mask_CP}', output=lsmask_ifile)
cdo.sellonlatbox(lon_min,lon_max,lat_min,lat_max, input=f'{DIR_nonCP_hfls}/{filename_hfls[M]}', output = hfls_ifile)

dataset = Dataset(f'{tasmax_ifile}', mode='r')
lon = dataset.variables['lon'][:]
lat = dataset.variables['lat'][:]
tim = dataset.variables['time'][:]
date = num2date(tim[:],dataset['time'].units)
units = dataset['time'].units
tasmax = dataset.variables['tasmax'][:].squeeze()
dataset.close()

dataset = Dataset(f'{hfls_ifile}', mode='r')
hfls = dataset.variables['hfls'][:].squeeze()
dataset.close()

dataset = Dataset(f'{tasmax_p90_ifile}', mode='r')
tasmax_p90 = dataset.variables['tasmax'][:].squeeze()
dataset.close()

dataset = Dataset(f'{lsmask_ifile}', mode='r')
mask = dataset['sftlf'][:]
dataset.close()

ntim = 10 #tasmax.shape[0]
nlat = tasmax.shape[1]
nlon = tasmax.shape[2]
nmodels = len(MODELS[M])

# Apply land-sea mask
tasmax.mask = mask==0
tasmax = tasmax.filled(np.nan)
hfls.mask = mask==0
hfls = hfls.filled(np.nan)
tasmax_p90.mask = mask==0
tasmax_p90 = tasmax_p90.filled(np.nan)


YEARS=np.arange(YSTART,YEND,1)

class structtype():
    pass

HW_persistence       = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_index_event_start = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_index_event_end   = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_HWMI              = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_mean_tmax         = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_max_tmax          = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_number            = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_hfls_deficit      = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan
HW_mean_hfls         = np.empty([len(YEARS),tasmax.shape[1],tasmax.shape[2]])*np.nan


cc = -1 
for Y in range(YSTART,YEND,1) : 
    
    print("MODEL:",MODELS[M])
    print("YEARS: ",Y)

    cc = cc + 1 
    # Count heatwave  
    cdo.selyear(Y, input=f'{tasmax_ifile}', output=tasmax_y_ifile)
    dataset  = Dataset(f'{tasmax_y_ifile}', mode='r');
    tasmax_y = dataset.variables['tasmax'][:]; dataset.close(); os.remove(tasmax_y_ifile)
    tasmax_y.mask = mask==0
    tasmax_y = tasmax_y.filled(np.nan).squeeze()
      
    cdo.selyear(Y, input=f'{hfls_ifile}', output=hfls_y_ifile)
    dataset  = Dataset(f'{hfls_y_ifile}', mode='r');
    if MODELS[M] == 'ETHZ' :
        hfls_y = dataset.variables['hfls'][:] * -1; dataset.close(); os.remove(hfls_y_ifile)
    else :
        hfls_y = dataset.variables['hfls'][:]; dataset.close(); os.remove(hfls_y_ifile)
    hfls_y.mask = mask==0
    hfls_y = hfls_y.filled(np.nan).squeeze()
    
    
    for ii in range(0,tasmax.shape[1]) :  #(30,31) : 
    
        #print("Current progress: " + str((ii)*100/len(lon[1])) + "%")
        
        for jj in range(0,tasmax.shape[2]) : #(20,21) : 
            
            if np.isnan(tasmax.squeeze()[0,ii,jj]) == False :
                
                pacchetti         = [] #np.empty(92)*np.nan
                dummy_index       = [] #np.empty(92)*np.nan
                index_event_start = [] #np.empty(92)*np.nan
                index_event_end   = [] #np.empty(92)*np.nan
                index_events      = [] # np.empty(92) * np.nan # It considers all the HWs in one year
                HWMI              = [] #np.empty(92)*np.nan 
                hfls_deficit      = [] #np.empty(92)*np.nan
            
            
                array  = tasmax_y[:,ii,jj]
                thresh = tasmax_p90[:,ii,jj].squeeze()
                array_hfls = hfls_y[:,ii,jj]
            
                Count = -1     
                i = 0 
                cval = 3 
		    
                IQR = scipy.stats.iqr(tasmax[:,ii,jj],rng=(25,75))
                p25 = np.percentile(tasmax[:,ii,jj],25)               
		    
                IQR_hfls = scipy.stats.iqr(hfls[:,ii,jj],rng=(25,75))
                p75_hfls = np.percentile(hfls[:,ii,jj],75)
            
            
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
                        HW_hfls_deficit[cc, ii, jj] = np.sum((array_hfls[index_most_intense_event] - p75_hfls)/IQR_hfls) * -1 
                        HW_mean_hfls[cc,ii,jj]      = np.nanmean(array_hfls[index_most_intense_event])
                        
                    else :
                        index_event = index_events
                        HW_persistence[cc,ii,jj] = np.array(pacchetti)
                        HW_index_event_start[cc,ii,jj] = np.min(index_event)
                        HW_index_event_end[cc,ii,jj]   = np.max(index_event)
                        HW_mean_tmax[cc,ii,jj] = np.nanmean(array[index_event])
                        HW_max_tmax[cc,ii,jj]  = np.nanmax(array[index_event])
                        HW_number[cc,ii,jj] = len(pacchetti)
                        HW_HWMI[cc, ii, jj] = np.sum((array[index_event] - p25) / IQR)
                        HW_hfls_deficit[cc, ii, jj] = np.sum((array_hfls[index_event] - p75_hfls)/IQR_hfls) * -1 
                        HW_mean_hfls[cc,ii,jj]      = np.nanmean(array_hfls[index_event])
                        
os.remove(tasmax_ifile)
os.remove(tasmax_p90_ifile)
os.remove(lsmask_ifile)
os.remove(hfls_ifile)
HW_hfls_deficit = HW_hfls_deficit.squeeze()

f_out = f'/users_home/cmcc/ls21622/tmp/{filename_tasmax_rcp85_out[M]}'

# *** define output file
out_nc1 = Dataset(f_out, "w", format="NETCDF4")

# *** define time to store in the netcdf file 
cdo.yearmean(input = f'{DIR_nonCP}/{filename[0]}',    output = 'ofile.nc', options = '-f nc') # Attention to the time_period from one model only
dataset = Dataset('ofile.nc', mode='r')
if MODELS[M] == 'UKMO' : 
    time_period = dataset.variables['time'][:-1]
else :
    time_period = dataset.variables['time'][:]
date = num2date(time_period,dataset['time'].units)
units = dataset['time'].units
dataset.close()

# **** define dimensions
lats = out_nc1.createDimension('y', nlat)
lons = out_nc1.createDimension('x', nlon)

if MODELS[M] == 'UKMO' : 
    time = out_nc1.createDimension('time',ntim-0)
else : 
    time = out_nc1.createDimension('time',ntim)

model_dim = out_nc1.createDimension('model', 1)


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
v8 = out_nc1.createVariable('HW_hfls_deficit',      'f4', ('time', 'y', 'x'))
v9 = out_nc1.createVariable('HW_mean_hfls',         'f4', ('time', 'y', 'x'))

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
v8.long_name = 'HW hfls deficit'
v8.units = 'unitless'
v9.long_name = 'mean HW hfls'
v9.units = 'w/m^2'

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
v8[:,:,:]= HW_hfls_deficit 
v9[:,:,:]= HW_mean_hfls 

model_var[0] = np.string_(MODELS[M])

out_nc1.close()

# Open output
dataset = Dataset(f'{f_out}', mode='r')