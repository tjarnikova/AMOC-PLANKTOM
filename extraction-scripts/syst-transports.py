#put python script here
import numpy as np
from cmocean import cm
import netCDF4 as nc
import matplotlib.pyplot as plt
import xarray as xr
from scipy.stats import pearsonr

import warnings
from datetime import datetime
warnings.filterwarnings('ignore')
from importlib import reload
import matplotlib.path as mpath
import glob
import pickle
import pandas as pd
import seawater
import time

def make_yearlist(yrst, yrend, dtype, tr, bd):
    yrs = np.arange(yrst,yrend+1,1)
    ylist = []
    for i in range(0,len(yrs)):
        ty = f'/gpfs/home/mep22dku/scratch/AMOC-PLANKTOM/transports/{tr}_{yrs[i]}_{dtype}_row{bd}*.nc'
        t2 = glob.glob(ty)
        #print(t2)
        ylist.append(t2[0])
    return ylist

def make_yearlist_FNAT(yrst, yrend, dtype, bd):
    yrs = np.arange(yrst,yrend+1,1)
    ylist = []
    for i in range(0,len(yrs)):
        ty = f'/gpfs/home/mep22dku/scratch/AMOC-PLANKTOM/transports/FNAT*_{yrs[i]}_{dtype}_row{bd}*.nc'
        t2 = glob.glob(ty)
        #print(t2)
        ylist.append(t2[0])
    return ylist


def get_cortable(amoc,yrst,yrend,trac,run,savenam, typ = 'pos', water = False, annual = False):
    depths = np.arange(0,24,1)
    rows = np.arange(73,110,1)

    # Initialize a DataFrame to hold the correlation coefficients
    correlation_table = pd.DataFrame(index=depths, columns=rows)
    
    if water:
        tvar = 'water_transports'
    else:
        tvar = 'tracer_transports'

    for L in rows:
        print(L)
        row60 = xr.open_mfdataset(make_yearlist(yrst,yrend,trac,run, bd = L))

        for d in depths:
            row60_h = row60[tvar].sel(time_counter=slice(f'{1950}-01-01', f'{2020}-12-31'))\
            .sel(depthlevel = d).sel(type = typ)
            if annual:
                row60_h = row60[tvar].sel(time_counter=slice(f'{1950}-01-01', f'{2020}-12-31'))\
                    .sel(depthlevel = d).sel(type = typ).groupby('time_counter.year').mean()
            variable = row60_h.values
            
            correlation_coefficient, pp = pearsonr(amoc, variable)
            if pp > 0.05:
                correlation_coefficient = 0
            correlation_table.at[d, L] = correlation_coefficient

    correlation_table.to_csv(savenam, index_label='Depth')


def get_cortable_FNAT(amoc,yrst,yrend,trac,savenam, typ = 'pos', water = False, annual = False):
    depths = np.arange(0,24,1)
    rows = np.arange(73,110,1)

    # Initialize a DataFrame to hold the correlation coefficients
    correlation_table = pd.DataFrame(index=depths, columns=rows)
    
    if water:
        tvar = 'water_transports'
    else:
        tvar = 'tracer_transports'

    for L in rows:
        print(L)
        row60 = xr.open_mfdataset(make_yearlist_FNAT(yrst,yrend,trac, bd = L))

        for d in depths:
            row60_h = row60[tvar].sel(time_counter=slice(f'{1950}-01-01', f'{2020}-12-31'))\
            .sel(depthlevel = d).sel(type = typ)
            if annual:
                row60_h = row60[tvar].sel(time_counter=slice(f'{1950}-01-01', f'{2020}-12-31'))\
                    .sel(depthlevel = d).sel(type = typ).groupby('time_counter.year').mean()
            variable = row60_h.values
            
            correlation_coefficient, pp = pearsonr(amoc, variable)
            if pp > 0.05:
                correlation_coefficient = 0
            correlation_table.at[d, L] = correlation_coefficient

    correlation_table.to_csv(savenam, index_label='Depth')


def get_transport_matrix(run,yrst,yrend,trac,typ = 'pos'):

    RVA0_amoc = xr.open_dataset('/gpfs/home/mep22dku/scratch/AMOC-PLANKTOM/data/AMOC_TOM12_TJ_RVA0_1945-2022.nc')
    amoc = RVA0_amoc.AMOC.sel(TIME=slice(f'{yrst}-01-01', f'{yrend}-12-31'))

    sn = '/gpfs/home/mep22dku/scratch/AMOC-PLANKTOM/transports/'
    snam = f'{run}-{typ}-{trac}-transport-vs-AMOC-{yrst}-{yrend}-monthly.csv'
    print(snam)
    savenam = f'{sn}{snam}'
    get_cortable(amoc,yrst,yrend,trac,run,savenam, typ, water = False, annual = False)
    
    snam = f'{run}-{typ}-water-transport-vs-AMOC-{yrst}-{yrend}-monthly.csv'
    print(snam)
    savenam = f'{sn}{snam}'
    get_cortable(amoc,yrst,yrend,trac,run,savenam, typ, water = True, annual = False)

    amoc = RVA0_amoc.AMOC.sel(TIME=slice(f'{yrst}-01-01', f'{yrend}-12-31')).groupby('TIME.year').mean()
    sn = '/gpfs/home/mep22dku/scratch/AMOC-PLANKTOM/transports/'
    snam = f'{run}-{typ}-{trac}-transport-vs-AMOC-{yrst}-{yrend}-annual.csv'
    print(snam)
    savenam = f'{sn}{snam}'
    get_cortable(amoc,yrst,yrend,trac,run,savenam, typ, water = False, annual = True)
    
    snam = f'{run}-{typ}-water-transport-vs-AMOC-{yrst}-{yrend}-annual.csv'
    print(snam)
    savenam = f'{sn}{snam}'
    get_cortable(amoc,yrst,yrend,trac,run,savenam, typ, water = True, annual = True)

    print(snam)

def get_transport_matrix_FNAT(yrst,yrend,typ = 'pos'):

    trac = 'DIC'
    RVA0_amoc = xr.open_dataset('/gpfs/home/mep22dku/scratch/AMOC-PLANKTOM/data/AMOC_TOM12_TJ_RVA0_1945-2022.nc')
    amoc = RVA0_amoc.AMOC.sel(TIME=slice(f'{yrst}-01-01', f'{yrend}-12-31'))

    sn = '/gpfs/home/mep22dku/scratch/AMOC-PLANKTOM/transports/'
    snam = f'FNAT-RVD0-RVB0-{typ}-{trac}-transport-vs-AMOC-{yrst}-{yrend}-monthly.csv'
    print(snam)
    savenam = f'{sn}{snam}'
    get_cortable_FNAT(amoc,yrst,yrend,trac,savenam, typ, water = False, annual = False)
    amoc = RVA0_amoc.AMOC.sel(TIME=slice(f'{yrst}-01-01', f'{yrend}-12-31')).groupby('TIME.year').mean()
    sn = '/gpfs/home/mep22dku/scratch/AMOC-PLANKTOM/transports/'
    snam = f'FNAT-RVD0-RVB0-{typ}-{trac}-transport-vs-AMOC-{yrst}-{yrend}-annual.csv'
    print(snam)
    savenam = f'{sn}{snam}'
    get_cortable_FNAT(amoc,yrst,yrend,trac,savenam, typ, water = False, annual = True)

    print(snam)

#get_transport_matrix('RVD0',1950,2020,'DIC','net')
# get_transport_matrix('RVD0',1950,2020,'NO3','net')

# get_transport_matrix('RVD0',1980,2020,'DIC','net')
# get_transport_matrix('RVD0',1980,2020,'NO3','net')

get_transport_matrix_FNAT(1980,2020,'net')
get_transport_matrix_FNAT(1950,2020,'net')