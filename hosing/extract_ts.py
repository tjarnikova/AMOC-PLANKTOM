#sulini
import numpy as np
import glob
import xarray as xr

def make_yearlist(yrst, yrend, tr, dtype = 'grid_T'):
    yrs = np.arange(yrst,yrend+1,1)
    ylist = []
    for i in range(0,len(yrs)):
        ty = f'/gpfs/home/mep22dku/scratch/ModelRuns/TOM12_TJ_{tr}/ORCA2_1m_{yrs[i]}*{dtype}*.nc'
        t2 = glob.glob(ty)
        #print(t2)
        #print(t2)
        ylist.append(t2[0])
    return ylist


dats = ['R4B1','R4A1','R4C1','R4D1','RWB0']
for dat in dats:

    print(dat)
    
    ySRs = xr.open_mfdataset(make_yearlist(1950,2020,dat))
    single_var = ySRs['wfo']
    new_ds = single_var.to_dataset()
    try:
        new_ds.to_netcdf(f'./{dat}_wfo.nc')
    except:
        print('no wfo')

    single_var = ySRs['sos']
    new_ds = single_var.to_dataset()
    try:
        new_ds.to_netcdf(f'./{dat}_sos.nc')
    except:
        print('no sos')

dats = ['R4B0']
for dat in dats:
    
    print(dat)

    ySRs = xr.open_mfdataset(make_yearlist(1950,1958,dat))
    single_var = ySRs['wfo']
    new_ds = single_var.to_dataset()
    try:
        new_ds.to_netcdf(f'./{dat}_wfo.nc')
    except:
        print('no wfo')

    single_var = ySRs['sos']
    new_ds = single_var.to_dataset()
    try:
        new_ds.to_netcdf(f'./{dat}_sos.nc')
    except:
        print('no sos')