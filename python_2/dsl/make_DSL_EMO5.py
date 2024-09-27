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
import functions_for_dsl
import sys
import cdo
from cdo import *   # python version
cdo = Cdo()
import os
import xesmf as xe

MODELS = 'EMO5'

DIR = '/work3/lorenzosangelantoni/scratch75/gcm_driven_experiment_from_nird/EMO-5'

PERIOD = 'historical'
VAR = 'pr'

YSTART = 1996
YEND = 2006

filename = 'EMO-5-pr_1996_2005_JJA.nc'
f_out = f'./outputs/{filename}' 

# Domain mask
boundaries = np.array([1,17.0,40,50]) ; lon_min,lon_max,lat_min,lat_max = boundaries

pr_ifile     = f'pr_EMO5.nc'
pr_y_ifile   = f'pr_EMO5_Y.nc'

# Load variables 
cdo.sellonlatbox(lon_min,lon_max,lat_min,lat_max, input=f'{DIR}/{filename}', output = pr_ifile)
#cdo.sellonlatbox(lon_min,lon_max,lat_min,lat_max, input=f'{mask}', output = lsmask_ifile)

dataset = Dataset(f'{pr_ifile}', mode='r')
lon = dataset.variables['lon'][:]
lat = dataset.variables['lat'][:]
tim = dataset.variables['time'][:]
date = num2date(tim[:],dataset['time'].units)
units = dataset['time'].units
pr = dataset.variables['pr'][:]
dataset.close()

ntim = 10 #tasmax.shape[0]
nlat = lat.shape[0]
nlon = lon.shape[1]
nmodels = len(MODELS)

YEARS=np.arange(YSTART,YEND,1)

class structtype():
    pass

DSL_index_event_start = np.empty([len(YEARS),pr.shape[1],pr.shape[2]])*np.nan
DSL_index_event_end   = np.empty([len(YEARS),pr.shape[1],pr.shape[2]])*np.nan
DSL_mean              = np.empty([len(YEARS),pr.shape[1],pr.shape[2]])*np.nan
DSL_max               = np.empty([len(YEARS),pr.shape[1],pr.shape[2]])*np.nan
DSL_number            = np.empty([len(YEARS),pr.shape[1],pr.shape[2]])*np.nan


cc = -1 
for Y in range(YSTART,YEND,1) : 
    
    print("MODEL:",MODELS)
    print("YEARS: ",Y)

    cc = cc + 1 
    # Count no rain  
    cdo.selyear(Y, input=f'{pr_ifile}', output=pr_y_ifile)
    dataset  = Dataset(f'{pr_y_ifile}', mode='r');
    pr_y = dataset.variables['pr'][:]; dataset.close(); os.remove(pr_y_ifile)
    #pr_y.mask = mask==0
    #pr_y = pr_y.filled(np.nan)

    for ii in range(0,pr.shape[1]) : 
    #for ii in range(82,83) :    
    
        for jj in range(0,pr.shape[2]) :        
        #for jj in range(104,105) : 
    
            pacchetti         = np.empty(92)*np.nan
            dummy_index       = np.empty(92)*np.nan
            index_event_start = np.empty(92)*np.nan
            index_event_end   = np.empty(92)*np.nan
            index_events      = np.empty([92,92]) * np.nan # It considers all the DSLs in one year
            index_event = [ structtype() for i in range(92) ]
            
            array  = pr_y[:,ii,jj]
            thresh = 1 
            
            Count = -1     
            i = 0 
            cval = 2 # Attention here, the value from which we start counting DSL. HWS was three cons. days here two.
            
            if pd.isna(array[0]) != True :
            
                pacchetti, index_event, index_events, index_event_start, index_event_end = functions_for_dsl.count_dsl(array,thresh,Count,i,cval,pacchetti,dummy_index,index_event_start,index_event_end,index_event,index_events)  
                index_events = index_events[~np.isnan(index_events)]
                index_events = [ int(x) for x in index_events ]


            if pd.isna(pacchetti[0]) != True :  
                pacchetti = pacchetti[~np.isnan(pacchetti)] 
                
                # Find the most intense envent 
                index_event_max = (index_event[(np.argmax(pacchetti))])   
                
                DSL_index_event_start[cc,ii,jj] = np.min(index_event_max) 
                DSL_index_event_end[cc,ii,jj]   = np.max(index_event_max) 
                                               
                #HW_mean_tmax[cc,ii,jj] = np.nanmean(array[index_event_max].values)  
                #HW_max_tmax[cc,ii,jj]  = np.nanmax(array[index_event_max].values)  		
                
                DSL_mean[cc,ii,jj] = np.nanmean(pacchetti)  
                DSL_max[cc,ii,jj]  = np.nanmax(pacchetti)  

                DSL_number[cc,ii,jj] = len(pacchetti) 
                
os.remove(pr_ifile)


# *** define output file
out_nc1 = Dataset(f_out, "w", format="NETCDF4")


# *** define time to store in the netcdf file 
cdo.yearmean(input = f'{DIR}/{filename}',    output = 'ofile.nc', options = '-f nc') #Attention to the time_period from one model only
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
v1 = out_nc1.createVariable('DSL_mean', 'f4', ('time', 'y', 'x'))
v2 = out_nc1.createVariable('DSL_max',  'f4', ('time', 'y', 'x'))
v3 = out_nc1.createVariable('DSL_number','f4', ('time', 'y', 'x'))
v4 = out_nc1.createVariable('DSL_index_event_start',  'f4', ('time', 'y', 'x'))
v5 = out_nc1.createVariable('DSL_index_event_end',  'f4', ('time', 'y', 'x'))
model_var = out_nc1.createVariable('model', 'S10', ('model',))

# **** define variables attribute
latitudes.long_name = 'latitude'
latitudes.units = 'degree_north'

longitudes.long_name = 'longitude'
longitudes.units = 'degree_east'

times.long_name = 'time'
times.units = units

v1.long_name = 'mean DSL'
v1.units = 'DAYS'
v2.long_name = 'max DSL'
v2.units = 'DAYS'
v3.long_name = 'DSL number'
v3.units = 'N'
v3.long_name = 'DSL start time index'
v4.units = 'summer day'
v5.long_name = 'DSL end time index'
v5.units = 'summer day'

# **** assign values to variables
latitudes[:,:] = lat
longitudes[:,:] = lon
times[:] = time_period.data[:]  
v1[:,:,:] =  DSL_mean
v1[:,:,:] =  DSL_mean
v2[:,:,:] =  DSL_max
v3[:,:,:] =  DSL_number
v4[:,:,:] =  DSL_index_event_start
v5[:,:,:] =  DSL_index_event_end
model_var = np.string_(MODELS)

out_nc1.close()

# Open output
dataset = Dataset(f'{f_out}', mode='r')

ds = []
ds = xr.Dataset(

    data_vars=dict(

         DSL_mean              = (["time","y", "x"], dataset.variables['DSL_mean'][:]),
         DSL_max               = (["time","y", "x"], dataset.variables['DSL_max'][:]),
         DSL_count             = (["time","y", "x"], dataset.variables['DSL_number'][:]),
         DSL_index_event_start = (["time","y", "x"], dataset.variables['DSL_index_event_start'][:]),         
         DSL_index_event_end   = (["time","y", "x"], dataset.variables['DSL_index_event_end'][:]),
),

    coords=dict(
        lon=(["y", "x"], dataset.variables['lon'][:]),
        lat=(["y", "x"], dataset.variables['lat'][:]),
        time=num2date(dataset['time'][:],dataset['time'].units)
        ),
)


#levels = np.linspace(10,40,30)
my_cmap = 'jet'
fig = plt.figure(figsize=(12,4))
data_crs = ccrs.PlateCarree() # Define transform
ax = plt.axes(projection=ccrs.PlateCarree()) # Define projection 
P=plt.pcolor(lon,lat,ds.DSL_mean.mean('time'),transform=ccrs.PlateCarree(),cmap=my_cmap) 
ax.coastlines(linewidth = 1)
ax.gridlines(linewidth = .5)
ax.set_facecolor('white')
cb = fig.colorbar(P,ax=ax, orientation='vertical',aspect=20,shrink=.7,pad=.015)
#plt.show()
plt.savefig(f'./figures/{MODELS}_his_test.png')
# ---------------

