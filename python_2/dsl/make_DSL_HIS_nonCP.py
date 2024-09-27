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

# Choose model through M variable 
M = int(sys.argv[1])
#M = 1

RES='nonCP'
MODELS = [ 'BCCR-AUTH','BTU','CMCC','CNRM','ETHZ','FZJ-IDL','HCLIM','ICTP','KIT','KNMI','UKMO'] 

# Compute DSL metrics from RCM simulations grid node-by-grid node 
root_dir = '/work3/lorenzosangelantoni/scratch75/gcm_driven_experiment_from_nird'

YSTART=1996
YEND = 2006

#if MODELS[M] == "UKMO" :
#    YSTART=2091
    

DIR = f'{root_dir}/FILES_nonCP/pr/historical'
DIR_OUT = f'{root_dir}/OUTPUT_HW'

MODELS = [ 'BCCR-AUTH','BTU','CMCC','CNRM','ETHZ','FZJ-IDL','HCLIM','ICTP','KIT','KNMI','UKMO','JLU'] 

PERIOD = 'rcp85'
VAR = 'pr'

filename = [ \
f'pr_EUR-15_NorESM1-ME_historical_r1i1p1_BCCR-WRF381DA_fpsconv-x1n2-v1_day_JJA.nc', \
f'pr_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CLMcom-CCLM4-8-17_cordex-x0n1-v1_day_JJA.nc', \
f'pr_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CMCC-CCLM5-0-9_v1_day_JJA.nc', \
f'pr_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CNRM-ALADIN63_v2_day_JJA.nc', \
f'pr_EUR-11_MPI_historical_r1i1p1_COSMO-pompa_5.0_2019.1_day_JJA.nc', \
f'pr_EUR-15_SMHI-EC-EARTH_historical_r12i1p1_FZJ-IDL-WRF381DA_fpsconv-x1n2-v1_day_JJA.nc',\
f'pr_ALP-12_ICHEC-EC-EARTH_historical_r12i1p1_HCLIMcom-HCLIM38-ALADIN_v1_day_JJA.nc', \
f'pr_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_ICTP-RegCM4-7_fpsconv-x2yn1-v1_day_JJA.nc', \
f'pr_EUR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_JJA.nc', \
f'pr_WCE-11_KNMI-EC-EARTH_historical_r14i1p1_KNMI-RACMO23E_v1_day_JJA.nc', \
f'pr_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_MOHC-HadREM3-GA7-05_v1_day_JJA.nc', \
]

f_out = f'./outputs/{filename[M]}' 

# Land - sea mask
mask = f'{root_dir}/OROG_SFTLF_WRF/cordex_fps_domain/sftlf_EUR-15_ECMWF-ERAINT_evaluation_r1i1p1_UCAN-WRF381BI_x1n2_fx.nc'
#mask=f'{root_dir}/OROG_SFTLF_WRF/cordex_fps_domain/sftlf_ALP-3_ECMWF-ERAINT_evaluation_r1i1p1_UCAN-WRF381BI_x1n2_fx_on_BCCR_grid.nc'

# Domain mask
boundaries = np.array([1,17.0,40,50]) ; lon_min,lon_max,lat_min,lat_max = boundaries

pr_ifile     = f'pr_{MODELS[M]}_{RES}.nc'
tmp_ifile     = f'tmp_{MODELS[M]}_{RES}.nc'
pr_y_ifile   = f'pr_{MODELS[M]}_Y_{RES}.nc'
lsmask_ifile = f'lsmask_{MODELS[M]}_{RES}.nc'

# Load variables 
cdo.sellonlatbox(lon_min,lon_max,lat_min,lat_max, input=f'{DIR}/{filename[M]}', output = pr_ifile)
cdo.sellonlatbox(lon_min,lon_max,lat_min,lat_max, input=f'{mask}', output = lsmask_ifile)

dataset = Dataset(f'{pr_ifile}', mode='r')
if MODELS[M] != 'UKMO' : 
    lat = dataset.variables['lat'][:]
    lon = dataset.variables['lon'][:]
else : 
    lat = dataset.variables['latitude'][:]
    lon = dataset.variables['longitude'][:]

tim = dataset.variables['time'][:]
date = num2date(tim[:],dataset['time'].units)
units = dataset['time'].units
pr = dataset.variables['pr'][:]
dataset.close()

cdo.sellonlatbox(lon_min,lon_max,lat_min,lat_max, input=f'{DIR}/{filename[0]}', output = tmp_ifile)
dataset = Dataset(f'{tmp_ifile}', mode='r')
lat_out = dataset.variables['lat'][:]
lon_out = dataset.variables['lon'][:]

dataset = Dataset(f'{lsmask_ifile}', mode='r')
mask = dataset['sftlf'][:]
dataset.close()

ntim = 10 #tasmax.shape[0]
nlat = lat_out.shape[0]
nlon = lon_out.shape[1]
nmodels = len(MODELS[M])


if MODELS[M] != 'BCCR-AUTH' : 
    # Interp to a common GCM grid (CNRM) the HW metrics 
    grid_in = {"lon": lon, "lat": lat}
    grid_out = {"lon": lon_out, "lat": lat_out}
    regridder = xe.Regridder(grid_in, grid_out, "bilinear")


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
    
    print("MODEL:",MODELS[M])
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
os.remove(tmp_ifile)
os.remove(lsmask_ifile)


# *** define output file
out_nc1 = Dataset(f_out, "w", format="NETCDF4")


# *** define time to store in the netcdf file 
#time_period = pd.date_range(f'{YSTART}-01-01',f'{YEND}-01-01',freq="A-JAN")
cdo.yearmean(input = f'{DIR}/{filename[0]}',    output = 'ofile.nc', options = '-f nc') #Attention to the time_period from one model only
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
latitudes[:,:] = lat_out
longitudes[:,:] = lon_out
times[:] = time_period.data[:]  
if MODELS[M] != 'BCCR-AUTH' : 
    v1[:,:,:] =  xr.DataArray(regridder(DSL_mean)).where(mask!=0).rename({'dim_0':'time','dim_1':'lat','dim_2':'lon'})
    v2[:,:,:] =  xr.DataArray(regridder(DSL_max)).where(mask!=0).rename({'dim_0':'time','dim_1':'lat','dim_2':'lon'})
    v3[:,:,:] =  xr.DataArray(regridder(DSL_number)).where(mask!=0).rename({'dim_0':'time','dim_1':'lat','dim_2':'lon'})
    v4[:,:,:] =  xr.DataArray(regridder(DSL_index_event_start)).where(mask!=0).rename({'dim_0':'time','dim_1':'lat','dim_2':'lon'})
    v5[:,:,:] =  xr.DataArray(regridder(DSL_index_event_end)).where(mask!=0).rename({'dim_0':'time','dim_1':'lat','dim_2':'lon'})
else : 
    v1[:,:,:] =  xr.DataArray(DSL_mean).where(mask!=0).rename({'dim_0':'time','dim_1':'lat','dim_2':'lon'})
    v2[:,:,:] =  xr.DataArray(DSL_max).where(mask!=0).rename({'dim_0':'time','dim_1':'lat','dim_2':'lon'})
    v3[:,:,:] =  xr.DataArray(DSL_number).where(mask!=0).rename({'dim_0':'time','dim_1':'lat','dim_2':'lon'})
    v4[:,:,:] =  xr.DataArray(DSL_index_event_start).where(mask!=0).rename({'dim_0':'time','dim_1':'lat','dim_2':'lon'})
    v5[:,:,:] =  xr.DataArray(DSL_index_event_end).where(mask!=0).rename({'dim_0':'time','dim_1':'lat','dim_2':'lon'})

model_var[0] = np.string_(MODELS[M])

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
P = plt.pcolor(lon_out,lat_out,ds.DSL_mean.mean('time'),transform=ccrs.PlateCarree(),cmap=my_cmap) 
ax.coastlines(linewidth = 1)
ax.gridlines(linewidth = .5)
ax.set_facecolor('white')
cb = fig.colorbar(P,ax=ax, orientation='vertical',aspect=20,shrink=.7,pad=.015)
#plt.show()
plt.savefig(f'./figures/{MODELS[M]}_nonCP_his_test.png')
# ---------------
