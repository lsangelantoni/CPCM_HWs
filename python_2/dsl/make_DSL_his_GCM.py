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
import xesmf as xe

# Compute DSL metrics from RCM simulations grid node-by-grid node 

YSTART = 1996
YEND   = 2006

boundaries = np.array([-0,30,38,60])
boundaries = np.array([-20,40,20,80])

# Choose model through M variable 
M = int(sys.argv[1])
#M = 10

DIR_GCM = '/scratch/lorenzosangelantoni/GCMs'

DIR_OUT = '/scratch/lorenzosangelantoni/gcm_driven_experiment_from_nird/OUTPUT_HW'

MODELS = [ 'CNRM-CM5_r1i1p1', 'EC-EARTH_r12i1p1', 'HadGEM2-ES_r1i1p1', 'MPI-ESM-LR_r1i1p1',  'NorESM1-M_r1i1p1'] 

filename_LSM = [ f'sftlf_fx_CNRM-CM5_historical_r0i0p0.nc', \
f'sftlf_fx_CNRM-CM5_historical_r0i0p0.nc', \
f'sftlf_fx_HadGEM2-ES_historical_r0i0p0.nc', \
f'sftlf_fx_MPI-ESM-LR_historical_r0i0p0.nc', \
f'sftlf_fx_NorESM1-M_historical_r0i0p0.nc', \
]

PERIOD = 'historical'
VAR = 'pr'

filename = [ \
f'pr_day_CNRM-CM5_historical_r1i1p1_box_1996_2005.nc', \
f'pr_day_EC-EARTH_historical_r12i1p1_1996_2005_box.nc', \
f'pr_day_HadGEM2-ES_historical_r1i1p1_box_1996_2005.nc', \
f'pr_day_MPI-ESM-LR_historical_r1i1p1_box_1996_2005.nc', \
f'pr_day_NorESM1-M_historical_r1i1p1_box_1996_2005.nc', \
]

filename_pr = f'{DIR_GCM}/{MODELS[M]}/surface/{filename[M]}'
ds_pr = xr.open_dataset(f'{filename_pr}') 
ds_pr = ds_pr.assign_coords({'institution':MODELS[M]})
ds_pr['pr'] = ds_pr['pr']*86400  

LSM_GCM=f'/scratch/lorenzosangelantoni/GCMs/LSM/{filename_LSM[M]}'
ds_LSM = xr.open_dataset(LSM_GCM)
lon_min,lon_max,lat_min,lat_max = boundaries

if lon_min < 0:
    #change lon coord to -180,180
    with xr.set_options(keep_attrs=True):
        ds_LSM['lon'] = ((ds_LSM.lon + 180) % 360) - 180
        ds_LSM = ds_LSM.sortby(ds_LSM.lon)


# Cut the subdomain for the LSMASK
mask = np.logical_and(np.logical_and(boundaries[0] <= ds_LSM.lon, ds_LSM.lon <= boundaries[1]),
                      np.logical_and(boundaries[2] <= ds_LSM.lat, ds_LSM.lat <= boundaries[3]))
da_lsmask = ds_LSM.sftlf.where(mask, drop=True)


# First interpolation:
# Interpolate LSM to SIM grid
method = 'bilinear'
def regrid(ds,ds_out,method):
    regridder = xe.Regridder(ds,ds_out, method = method)
    ds = regridder(ds,ds_out,method)
    return ds
ds_out = xr.open_dataset(f'{filename_pr}')
da_lsmask_remap = regrid(da_lsmask,ds_out,method)

ds_pr['pr'] = ds_pr.pr.where(da_lsmask_remap.values > 50)

YEND = 2006
#YEND=1997

YEARS=np.arange(YSTART,YEND,1)

class structtype():
    pass


DSL_index_event_start = np.empty([len(YEARS),ds_pr.pr.shape[1],ds_pr.pr.shape[2]])*np.nan
DSL_index_event_end   = np.empty([len(YEARS),ds_pr.pr.shape[1],ds_pr.pr.shape[2]])*np.nan
DSL_mean              = np.empty([len(YEARS),ds_pr.pr.shape[1],ds_pr.pr.shape[2]])*np.nan
DSL_max               = np.empty([len(YEARS),ds_pr.pr.shape[1],ds_pr.pr.shape[2]])*np.nan
DSL_number            = np.empty([len(YEARS),ds_pr.pr.shape[1],ds_pr.pr.shape[2]])*np.nan


cc = -1 
for Y in range(YSTART,YEND,1) : 
    
    print("MODEL:",MODELS[M])
    print("YEARS: ",Y)

    cc = cc + 1 
    # Count heatwave  
    ds_pr_yr = ds_pr.sel(time=ds_pr.time.dt.year.isin(Y))
    ds_pr_yr = ds_pr_yr.sel(time=ds_pr_yr.time.dt.month.isin([6,7,8]))

    for ii in range(0,ds_pr.pr.shape[1]) : 
    #for ii in range(82,83) :    
    
        print("Current progress: " + str((ii)*100/len(ds_pr.lon)) + "%")
        
        for jj in range(0,ds_pr.pr.shape[2]) :        
        #for jj in range(104,105) : 
            pacchetti         = np.empty(92)*np.nan
            dummy_index       = np.empty(92)*np.nan
            index_event_start = np.empty(92)*np.nan
            index_event_end   = np.empty(92)*np.nan
            index_events      = np.empty([92,92]) * np.nan # It considers all the DSLs in one year
            index_event = [ structtype() for i in range(92) ]
            
            array  = ds_pr_yr.pr[:,ii,jj]
            thresh = 1 #ds_pr_p90.pr[:,ii,jj]
            
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
                index_event_max = (index_event[(np.argmax(pacchetti))].values)   
                
                DSL_index_event_start[cc,ii,jj] = np.min(index_event_max) 
                DSL_index_event_end[cc,ii,jj]   = np.max(index_event_max) 
                                               
                #HW_mean_tmax[cc,ii,jj] = np.nanmean(array[index_event_max].values)  
                #HW_max_tmax[cc,ii,jj]  = np.nanmax(array[index_event_max].values)  		
                
                DSL_mean[cc,ii,jj] = np.nanmean(pacchetti)  
                DSL_max[cc,ii,jj]  = np.nanmax(pacchetti)  

                DSL_number[cc,ii,jj] = len(pacchetti) 
                
#time = pd.date_range("1996-01-01","2006-01-01",freq="A-JAN")
time = pd.date_range(f'{YSTART}-01-01',f'{YEND}-01-01',freq="A-JAN")
ds = []
ds = xr.Dataset(

    data_vars=dict(
        
         DSL_index_event_start=(["time","y", "x"], DSL_index_event_start),
         DSL_index_event_end=(["time","y", "x"], DSL_index_event_end),
         
         DSL_mean=(["time","y", "x"], DSL_mean),
         DSL_max=(["time","y", "x"], DSL_max),
         DSL_number=(["time","y", "x"], DSL_number),
        
        ),

    coords=dict(
        lon=("x", ds_pr.lon.values),
        lat=("y", ds_pr.lat.values),
        time=time
        ),
)

ds = ds.assign_coords(model=MODELS[M])
ds = ds.expand_dims('model') 

# Second interpolation:
# Interpolate models to a common grid
method = 'bilinear'
def regrid(ds,ds_out,method):
    regridder = xe.Regridder(ds,ds_out, method = method)
    ds = regridder(ds,ds_out,method)
    return ds

if MODELS[M] != 'CNRM-CM5_r1i1p1' : #'HadGEM2-ES_r1i1p1' :
    ds_out = xr.open_dataset(f'/scratch/lorenzosangelantoni/GCMs/CNRM-CM5_r1i1p1/surface/pr_day_CNRM-CM5_historical_r1i1p1_box_1996_2005.nc')
    ds_out = xr.DataArray.to_dataset(ds_out.pr)
    ds_remap = regrid(ds,ds_out,method)
    # Write netcdf with the remapped outputs
    ds_remap.to_netcdf(path=f'{DIR_OUT}/{filename[M]}')
    print('SHAPE OF FILES SENT TO PRINT:')
    print(ds_remap)
else:
    ds = ds.rename({'x':'lon','y':'lat'})
    ds.to_netcdf(path=f'{DIR_OUT}/{filename[M]}')
    print('SHAPE OF FILES SENT TO PRINT:')
    print(ds)

# Open output
ds_out = []
ds_out = xr.open_dataset(f'{DIR_OUT}/{filename[M]}') 
print(ds_out) 

# Map 
# -------------------------------------------------------------------
plt.clf()
levels = np.linspace(2,15,14)
my_cmap = 'jet'
plt.figure(figsize=(12,4))
data_crs = ccrs.PlateCarree() # Define transform
ax = plt.axes(projection=ccrs.PlateCarree()) # Define projection 
plot = ds_out.DSL_mean.mean('time').plot(ax=ax, x='lon', y='lat',transform=ccrs.PlateCarree(),
                                                        cmap=my_cmap,levels=levels, cbar_kwargs={"label": 'days'})

ax.coastlines(linewidth = 1)
ax.gridlines(linewidth = .5)
ax.set_facecolor('white')
ax.set_extent([boundaries[0], boundaries[1],boundaries[2],boundaries[3]])
plt.savefig(f'{MODELS[M]}_{PERIOD}_test.png')
# ---------------




