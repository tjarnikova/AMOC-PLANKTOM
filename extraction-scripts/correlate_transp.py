
from scipy.stats import pearsonr
import numpy as np
import xarray as xr
import pandas as pd
import glob

RVA0_amoc = xr.open_dataset('/gpfs/home/mep22dku/scratch/AMOC-PLANKTOM/data/AMOC_TOM12_TJ_RVA0_1945-2022.nc')

def make_yearlist(yrst, yrend, dtype, tr, bd):
    yrs = np.arange(yrst,yrend+1,1)
    ylist = []
    for i in range(0,len(yrs)):
        ty = f'/gpfs/home/mep22dku/scratch/AMOC-PLANKTOM/transports/{tr}_{yrs[i]}_{dtype}_row{bd}*.nc'
        t2 = glob.glob(ty)
        #print(t2)
        ylist.append(t2[0])
    return ylist

amoc = RVA0_amoc.AMOC.sel(TIME=slice(f'{1950}-01-01', f'{2020}-12-31'))

# Define the range of depths and latitudes
depths = np.arange(0,31,1)
rows = np.arange(60,110,1)

# Initialize a DataFrame to hold the correlation coefficients
correlation_table = pd.DataFrame(index=depths, columns=rows)
correlation_table2 = pd.DataFrame(index=depths, columns=rows)
correlation_table3 = pd.DataFrame(index=depths, columns=rows)
correlation_table4 = pd.DataFrame(index=depths, columns=rows)

# Populate the table with correlation coefficients
for L in rows:
    print(L)
    A = xr.open_mfdataset(make_yearlist(1950, 2021, 'DIC', 'RVA0', bd = L))
    B = xr.open_mfdataset(make_yearlist(1950, 2021, 'DIC', 'RVB0', bd = L))
    C = xr.open_mfdataset(make_yearlist(1950, 2021, 'DIC', 'RVC0', bd = L))
    D = xr.open_mfdataset(make_yearlist(1950, 2021, 'DIC', 'RVD0', bd = L))

    TNat_ns = D - B
    TAnt_ss = C - B
    TAnt_ns = A - B - TNat_ns - TAnt_ss
    TTot = A - B 


    for d in depths:
        tdat = TNat_ns.tracer_transports.sel(time_counter=slice(f'{1950}-01-01', f'{2020}-12-31'))\
        .sel(depthlevel = d).sel(type = 'net')
        variable = tdat.values
        correlation_coefficient, _ = pearsonr(amoc, variable)
        correlation_table.at[d, L] = correlation_coefficient

    for d in depths:
        tdat = TAnt_ss.tracer_transports.sel(time_counter=slice(f'{1950}-01-01', f'{2020}-12-31'))\
        .sel(depthlevel = d).sel(type = 'net')
        variable = tdat.values
        correlation_coefficient, _ = pearsonr(amoc, variable)
        correlation_table2.at[d, L] = correlation_coefficient

    for d in depths:
        tdat = TAnt_ns.tracer_transports.sel(time_counter=slice(f'{1950}-01-01', f'{2020}-12-31'))\
        .sel(depthlevel = d).sel(type = 'net')
        variable = tdat.values
        correlation_coefficient, _ = pearsonr(amoc, variable)
        correlation_table3.at[d, L] = correlation_coefficient            

    for d in depths:
        tdat = TTot.tracer_transports.sel(time_counter=slice(f'{1950}-01-01', f'{2020}-12-31'))\
        .sel(depthlevel = d).sel(type = 'net')
        variable = tdat.values
        correlation_coefficient, _ = pearsonr(amoc, variable)
        correlation_table4.at[d, L] = correlation_coefficient        

print(correlation_table)
correlation_table.to_csv('./transports/TNat_ns_net_DIC_transport_vs_amoc_1950_2020.csv', index_label='Depth')
correlation_table2.to_csv('./transports/TAnt_ss_net_DIC_transport_vs_amoc_1950_2020.csv', index_label='Depth')
correlation_table3.to_csv('./transports/TAnt_ns_net_DIC_transport_vs_amoc_1950_2020.csv', index_label='Depth')
correlation_table4.to_csv('./transports/TTot_net_DIC_transport_vs_amoc_1950_2020.csv', index_label='Depth')

