#!/usr/bin/env python
# coding: utf-8


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
DIR_nonCP     = f'{root_dir}/FILES_nonCP/tasmax/historical'
DIR_nonCP_hfls= f'{root_dir}/FILES_nonCP/hfls/historical'
DIR_nonCP_p90 = f'{root_dir}/FILES_nonCP/tasmax/historical/p90'

#f_out = f'/work3/lorenzosangelantoni/scratch75/gcm_driven_experiment_from_nird/scripts/python_2/hws/outputs_matlab_like/historical/{RES}/tasmaxc{filename[M]}'

RES='nonCP'

YSTART = 1996
YEND   = 2006

# Choose model through M variable 
M = int(sys.argv[1])
#M = 0

MODELS = [ 'BCCR-AUTH','BTU','CMCC','CNRM','ETHZ','FZJ-IDL','HCLIM','ICTP','KIT','KNMI','UKMO'] 

filename = [ \
f'tasmax_EUR-15_NorESM1-ME_historical_r1i1p1_BCCR-WRF381DA_fpsconv-x1n2-v1_day_JJA.nc', \
f'tasmax_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CLMcom-CCLM4-8-17_cordex-x0n1-v1_day_JJA.nc', \
f'tasmax_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CMCC-CCLM5-0-9_v1_day_JJA.nc', \
f'tasmax_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CNRM-ALADIN63_v2_day_JJA.nc', \
f'tasmax_EUR-11_MPI_historical_r1i1p1_COSMO-pompa_5.0_2019.1_day_JJA.nc', \
f'tasmax_EUR-15_SMHI-EC-EARTH_historical_r12_FZJ-IBG3-WRF381CA_v1_day_JJA.nc', \
f'tasmax_ALP-12_ICHEC-EC-EARTH_historical_r12i1p1_HCLIMcom-HCLIM38-ALADIN_v1_day_JJA.nc', \
f'tasmax_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_ICTP-RegCM4-7_fpsconv-x2yn1-v1_day_JJA.nc', \
f'tasmax_EUR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_JJA.nc', \
f'tasmax_WCE-11_KNMI-EC-EARTH_historical_r14i1p1_KNMI-RACMO23E_v1_day_JJA.nc', \
f'tasmax_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_MOHC-HadREM3-GA7-05_v1_day_JJA.nc', \
]

filename_p90 = [ \
f'tasmax_EUR-15_NorESM1-ME_historical_r1i1p1_BCCR-WRF381DA_fpsconv-x1n2-v1_day.nc', \
f'tasmax_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CLMcom-CCLM4-8-17_cordex-x0n1-v1_day.nc', \
f'tasmax_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CMCC-CCLM5-0-9_v1_day.nc', \
f'tasmax_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CNRM-ALADIN63_v2_day.nc', \
f'tasmax_EUR-11_MPI_historical_r1i1p1_COSMO-pompa_5.0_2019.1_day.nc', \
f'tasmax_EUR-15_SMHI-EC-EARTH_historical_r12_FZJ-IBG3-WRF381CA_v1_day.nc', \
f'tasmax_ALP-12_ICHEC-EC-EARTH_historical_r12i1p1_HCLIMcom-HCLIM38-ALADIN_v1_day.nc', \
f'tasmax_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_ICTP-RegCM4-7_fpsconv-x2yn1-v1_day.nc' ,\
f'tasmax_EUR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day.nc', \
f'tasmax_WCE-11_KNMI-EC-EARTH_historical_r14i1p1_KNMI-RACMO23E_v1_day.nc', \
f'tasmax_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_MOHC-HadREM3-GA7-05_v1_day.nc', \
]

filename_hfls = [ \
f'hfls_EUR-15_NorESM1-ME_historical_r1i1p1_BCCR-WRF381DA_fpsconv-x1n2-v1_day_JJA.nc', \
f'hfls_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_JJA.nc', \
f'hfls_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CMCC-CCLM5-0-9_v1_day_JJA.nc', \
f'hfls_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CNRM-ALADIN63_v2_day_JJA.nc', \
f'hfls_EUR-11_MPI_historical_r1i1p1_COSMO-pompa_5.0_2019.1_day_JJA.nc', \
f'hfls_EUR-15_SMHI-EC-EARTH_historical_r12_FZJ-IBG3-WRF381CA_v1_day_JJA.nc', \
f'hfls_ALP-12_ICHEC-EC-EARTH_historical_r12i1p1_HCLIMcom-HCLIM38-ALADIN_v1_day_JJA.nc', \
f'hfls_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_ICTP-RegCM4-7_fpsconv-x2yn1-v1_day_JJA.nc', \
f'hfls_EUR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_JJA.nc', \
f'hfls_WCE-11_KNMI-EC-EARTH_historical_r14i1p1_KNMI-RACMO23E_v1_day_JJA.nc', \
f'hfls_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_MOHC-HadREM3-GA7-05_v1_day_JJA.nc', \
]


f_out = f'{root_dir}/scripts/python_2/hws/TMP_TO_TEST_SIMILARITY/{filename[M]}'

# Land - sea mask
mask_nonCP = f'{root_dir}/OROG_SFTLF_WRF/cordex_fps_domain/sftlf_EUR-15_ECMWF-ERAINT_evaluation_r1i1p1_UCAN-WRF381BI_x1n2_fx.nc'

# Domain mask
boundaries = np.array([1,17.0,40,50]) ; lon_min,lon_max,lat_min,lat_max = boundaries

tasmax_ifile     = f'tasmax_{MODELS[M]}_{RES}.nc'
tasmax_p90_ifile = f'tasmax_p90_{MODELS[M]}_{RES}.nc' 
lsmask_ifile     = f'lsmask_{MODELS[M]}_{RES}.nc'
tasmax_y_ifile   = f'tasmax_{MODELS[M]}_Y_{RES}.nc'
hfls_ifile       = f'hfls_{MODELS[M]}_{RES}.nc'
hfls_y_ifile     = f'hfls_{MODELS[M]}_Y_{RES}.nc'

# Load variables 
cdo.sellonlatbox(lon_min,lon_max,lat_min,lat_max, input=f'{DIR_nonCP}/{filename[M]}', output=tasmax_ifile)
cdo.sellonlatbox(lon_min,lon_max,lat_min,lat_max, input=f'{DIR_nonCP_p90}/{filename_p90[M]}', output=tasmax_p90_ifile)
cdo.sellonlatbox(lon_min,lon_max,lat_min,lat_max, input=f'{mask_nonCP}', output=lsmask_ifile)
cdo.sellonlatbox(lon_min,lon_max,lat_min,lat_max, input=f'{DIR_nonCP_hfls}/{filename_hfls[M]}', output = hfls_ifile)


dataset = Dataset(f'{tasmax_ifile}', mode='r')
lon = dataset.variables['lon'][:]
lat = dataset.variables['lat'][:]
tim = dataset.variables['time'][:]
date = num2date(tim[:],dataset['time'].units)
units = dataset['time'].units
tasmax = dataset.variables['tasmax'][:]
dataset.close()

dataset = Dataset(f'{hfls_ifile}', mode='r')
hfls = dataset.variables['hfls'][:]
dataset.close()

dataset = Dataset(f'{tasmax_p90_ifile}', mode='r')
tasmax_p90 = dataset.variables['tasmax'][:]
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
    tasmax_y = tasmax_y.filled(np.nan)
      
    cdo.selyear(Y, input=f'{hfls_ifile}', output=hfls_y_ifile)
    dataset  = Dataset(f'{hfls_y_ifile}', mode='r');
    if MODELS[M] == 'ETHZ' :
        hfls_y = dataset.variables['hfls'][:] * -1; dataset.close(); os.remove(hfls_y_ifile)
    else :
        hfls_y = dataset.variables['hfls'][:]; dataset.close(); os.remove(hfls_y_ifile)
    hfls_y.mask = mask==0
    hfls_y = hfls_y.filled(np.nan)
    
    
    for ii in range(0,tasmax.shape[1]) :  #(30,31) : 
    
        print("Current progress: " + str((ii)*100/len(lon[1])) + "%")
        
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


# *** define output file
out_nc1 = Dataset(f_out, "w", format="NETCDF4")


# *** define time to store in the netcdf file 
#time_period = pd.date_range(f'{YSTART}-01-01',f'{YEND}-01-01',freq="A-JAN")
cdo.yearmean(input = f'{DIR_nonCP}/{filename[0]}',    output = 'ofile.nc', options = '-f nc') #Attention to the time_period from one model only
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
