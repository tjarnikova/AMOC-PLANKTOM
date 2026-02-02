import xarray as xr
import numpy as np
import glob

BASE = 'TOM12_TJ_LA50'
HOS = 'TOM12_TJ_LAH3'
CCLIM = 'TOM12_TJ_LA50'
bdir = '/gpfs/data/greenocean/software/runs/clims/'



def make_yearlist(yrst, yrend, dtype, tr, baseDir):
    yrs = np.arange(yrst,yrend+1,1)
    ylist = []
    for i in range(0,len(yrs)):
        ty = f'{baseDir}/{tr}/ORCA2_1m_{yrs[i]}*{dtype}*.nc'
        t2 = glob.glob(ty)
        #print(t2)
        ylist.append(t2[0])
    return ylist

#grid_T, diad_T, ptrc_T, LoP_T
varty = ['grid_T', 'diad_T', 'ptrc_T', 'LoP_T']
#varty = ['grid_T']
mods = [BASE,HOS]#,CCLIM]
#mods = [BASE]
ys = 2010; ye = 2019
for mod in mods:
    for tvar in varty:
        try:
            w = xr.open_mfdataset(make_yearlist(ys,ye,tvar,mod,bdir))
            w = w.groupby('time_counter.month').mean('time_counter')
            w = w.rename({'month': 'time'})
            w.to_netcdf(f'{bdir}/{mod}/ORCA2_1m_clim_{ys}_{ye}_{tvar}.nc')
            print(f'YES {bdir}/{mod}/ORCA2_1m_clim_{ys}_{ye}_{tvar}.nc')
        except:
            print(f'--no for {bdir}/{mod}/ORCA2_1m_clim_{ys}_{ye}_{tvar}.nc')
