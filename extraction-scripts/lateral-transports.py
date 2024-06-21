#put python script here
import numpy as np
from cmocean import cm
import netCDF4 as nc
import matplotlib.pyplot as plt
import xarray as xr

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
plt.rcParams.update({'font.size': 12})
font = {'family' : 'normal',
'weight' : 'normal',
'size'   : 12}

plt.rc('font', **font)

N1_x1 = 115; N1_x2 = 117; N1_y = 125
N2_x1 = 126; N2_x2 = 129; N2_y = 122
N3_x1 = 134; N3_x2 = 141; N3_y = 121
S1_x1 = 100; S1_x2 = 131; S1_y = 89

def get_mask(x1,x2,y):
    tmesh = xr.open_dataset('/gpfs/data/greenocean/software/resources/regrid/mesh_mask3_6.nc')
    vmask = tmesh['vmask'][0,:,y,:]
    e3v = tmesh['e3v_0'][0,:,y,:]
    e1v = (tmesh.e1v[0,y,:])
    e1vb = np.zeros([31,182])
    for i in range(0,31):
       e1vb[i,:] =  e1v
    csize_N1 = vmask*e1vb*e3v
    
    csize_N1[:,0:x1] = 0
    csize_N1[:,x2+1:-1]=0
    csize_N1v = np.copy(csize_N1)
    csize_N1v[np.where(csize_N1 < 1)] = np.nan
    
    return csize_N1v

def get_transport(year=1984,model='RVA0',transectname='S1',x1=S1_x1,x2=S1_x2,y=S1_y,\
                  trac='water',northIsIn = True, showPlot = False, units = 'kmol/s'):
    w1 = time.time()
    savenam = f'./transports/{model}_{year}_{trac}_{transectname}.nc'

    tmask = get_mask(x1,x2,y)
    tmb = np.zeros([12,np.shape(tmask)[0],np.shape(tmask)[1]])

    for i in range(0,12):
        tmb[i,:,:] = tmask
    tmb1 = tmb[:,:,x1:x2+1]

    tdi = '/gpfs/data/greenocean/software/runs/TOM12_TJ_'
    tfi = xr.open_dataset(f'{tdi}{model}/ORCA2_1m_{year}0101_{year}1231_ptrc_T.nc')
    vfi = xr.open_dataset(f'{tdi}{model}/ORCA2_1m_{year}0101_{year}1231_grid_V.nc')

    vo = vfi['vomecrty'][:,:,y,x1:x2+1].values
    transp = vo * tmb1

    to1 =  tfi[trac][:,:,y,x1:x2+1].values
    to2 =  tfi[trac][:,:,y+1,x1:x2+1].values
    to = (to1+to2)/2
    totransp = transp * to
    
    if (northIsIn == False):
        transp = -1 * transp
        totransp = -1 * totransp

    if showPlot:
        fact = 1.1
        fig, axs = plt.subplots(1,5, figsize=(20*fact, 3*fact), facecolor='w', edgecolor='k')
        axs = axs.ravel()

        w = axs[0].pcolormesh(tmb1[0,:,:], cmap = plt.cm.Spectral)
        plt.colorbar(w, ax = axs[0])
        w = axs[1].pcolormesh(vo[0,:,:], cmap = plt.cm.Spectral)
        plt.colorbar(w, ax = axs[1])
        w = axs[2].pcolormesh(transp[0,:,:]/1e6, cmap = cm.balance, vmin = -1, vmax = 1)
        plt.colorbar(w, ax = axs[2])

        ttrac = to[0,:,:]
        ttrac[np.where(np.isnan(tmb1[0,:,:]))] = np.nan
        w = axs[3].pcolormesh(ttrac, cmap = plt.cm.Spectral)
        plt.colorbar(w, ax = axs[3])        

        w = axs[4].pcolormesh(totransp[0,:,:], cmap = plt.cm.Spectral)
        plt.colorbar(w, ax = axs[4])      

        tits = ['mask','vomecrty','transport (Sv)',trac,f'{trac} transport ({units})']
        for i in range(0,5):
            axs[i].invert_yaxis()
            axs[i].set_title(tits[i])

    #transp
    #totransp

    transp_res = np.zeros([12,5,3])
    totransp_res = np.zeros([12,5,3])

    for i in range(0,12):

        tt = transp[i,0:10,:]
        ttpos = (np.nansum(tt[np.where(tt>=0)]))
        ttneg = (np.nansum(tt[np.where(tt<0)]))
        ttnet = np.nansum(tt)
        transp_res[i,0,0] = ttpos
        transp_res[i,0,1] = ttneg
        transp_res[i,0,2] = ttnet

        tt = transp[i,0:16,:]
        ttpos = (np.nansum(tt[np.where(tt>=0)]))
        ttneg = (np.nansum(tt[np.where(tt<0)]))
        ttnet = np.nansum(tt)
        transp_res[i,1,0] = ttpos
        transp_res[i,1,1] = ttneg
        transp_res[i,1,2] = ttnet

        tt = transp[i,0:20,:]
        ttpos = (np.nansum(tt[np.where(tt>=0)]))
        ttneg = (np.nansum(tt[np.where(tt<0)]))
        ttnet = np.nansum(tt)
        transp_res[i,2,0] = ttpos
        transp_res[i,2,1] = ttneg
        transp_res[i,2,2] = ttnet

        tt = transp[i,0:22,:]
        ttpos = (np.nansum(tt[np.where(tt>=0)]))
        ttneg = (np.nansum(tt[np.where(tt<0)]))
        ttnet = np.nansum(tt)
        transp_res[i,3,0] = ttpos
        transp_res[i,3,1] = ttneg
        transp_res[i,3,2] = ttnet

        tt = transp[i,:,:]
        ttpos = (np.nansum(tt[np.where(tt>=0)]))
        ttneg = (np.nansum(tt[np.where(tt<0)]))
        ttnet = np.nansum(tt)
        transp_res[i,4,0] = ttpos
        transp_res[i,4,1] = ttneg
        transp_res[i,4,2] = ttnet

        ####################
        tt = totransp[i,0:10,:]
        ttpos = (np.nansum(tt[np.where(tt>=0)]))
        ttneg = (np.nansum(tt[np.where(tt<0)]))
        ttnet = np.nansum(tt)
        totransp_res[i,0,0] = ttpos
        totransp_res[i,0,1] = ttneg
        totransp_res[i,0,2] = ttnet

        tt = totransp[i,0:16,:]
        ttpos = (np.nansum(tt[np.where(tt>=0)]))
        ttneg = (np.nansum(tt[np.where(tt<0)]))
        ttnet = np.nansum(tt)
        totransp_res[i,1,0] = ttpos
        totransp_res[i,1,1] = ttneg
        totransp_res[i,1,2] = ttnet

        tt = totransp[i,0:20,:]
        ttpos = (np.nansum(tt[np.where(tt>=0)]))
        ttneg = (np.nansum(tt[np.where(tt<0)]))
        ttnet = np.nansum(tt)
        totransp_res[i,2,0] = ttpos
        totransp_res[i,2,1] = ttneg
        totransp_res[i,2,2] = ttnet

        tt = totransp[i,0:22,:]
        ttpos = (np.nansum(tt[np.where(tt>=0)]))
        ttneg = (np.nansum(tt[np.where(tt<0)]))
        ttnet = np.nansum(tt)
        totransp_res[i,3,0] = ttpos
        totransp_res[i,3,1] = ttneg
        totransp_res[i,3,2] = ttnet

        tt = totransp[i,:,:]
        ttpos = (np.nansum(tt[np.where(tt>=0)]))
        ttneg = (np.nansum(tt[np.where(tt<0)]))
        ttnet = np.nansum(tt)
        totransp_res[i,4,0] = ttpos
        totransp_res[i,4,1] = ttneg
        totransp_res[i,4,2] = ttnet

    times = pd.date_range(f"{year}/01/01",f"{year+1}/01/01",freq='MS',closed='left')
    #print(times)



    data_vars = {'tracer_transports':(['time_counter', 'depthlevel', 'type'], totransp_res,
    {'units': units,
    'long_name':trac}),
                 'water_transports':(['time_counter', 'depthlevel', 'type'], transp_res,
    {'units': 'm3/s',
    'long_name':'water only'}),
    }
    # define coordinates
    coords = {'time_counter':  times,
    'depthlevel': (['100','200','500','1000','tot']),
    'type': (['pos','neg','net']),\
             }
    # define global attributes
    attrs = {'made in':'scratch/AMOC-PLANKTOM/lateral-transports.ipynb',
    'desc': 'lateral transports'
    }
    ds = xr.Dataset(data_vars=data_vars,
    coords=coords,
    attrs=attrs)
    

    w2 = time.time()
    print(f'{savenam} {w2-w1}')# print(w2-w1)
    ds.to_netcdf(savenam)
        
        


for year in range(1950,2023):
    get_transport(year=year,model='RVA0',transectname='S1',x1=S1_x1,x2=S1_x2,y=S1_y,\
                     trac='DIC',northIsIn = True, showPlot = False)
for year in range(1950,2023):
    get_transport(year=year,model='RVA0',transectname='N1',x1=N1_x1,x2=N1_x2,y=N1_y,\
                     trac='DIC',northIsIn = True, showPlot = False)
for year in range(1950,2023):
    get_transport(year=year,model='RVA0',transectname='N2',x1=N2_x1,x2=N2_x2,y=N2_y,\
                     trac='DIC',northIsIn = True, showPlot = False)
for year in range(1950,2023):
    get_transport(year=year,model='RVA0',transectname='N3',x1=N3_x1,x2=N3_x2,y=N3_y,\
                     trac='DIC',northIsIn = True, showPlot = False)
    
for year in range(1950,2023):
    get_transport(year=year,model='RVB0',transectname='S1',x1=S1_x1,x2=S1_x2,y=S1_y,\
                     trac='DIC',northIsIn = True, showPlot = False)
for year in range(1950,2023):
    get_transport(year=year,model='RVB0',transectname='N1',x1=N1_x1,x2=N1_x2,y=N1_y,\
                     trac='DIC',northIsIn = True, showPlot = False)
for year in range(1950,2023):
    get_transport(year=year,model='RVB0',transectname='N2',x1=N2_x1,x2=N2_x2,y=N2_y,\
                     trac='DIC',northIsIn = True, showPlot = False)
for year in range(1950,2023):
    get_transport(year=year,model='RVB0',transectname='N3',x1=N3_x1,x2=N3_x2,y=N3_y,\
                     trac='DIC',northIsIn = True, showPlot = False)
    
#
for year in range(1950,2023):
    get_transport(year=year,model='RVC0',transectname='S1',x1=S1_x1,x2=S1_x2,y=S1_y,\
                     trac='DIC',northIsIn = True, showPlot = False)
for year in range(1950,2023):
    get_transport(year=year,model='RVC0',transectname='N1',x1=N1_x1,x2=N1_x2,y=N1_y,\
                     trac='DIC',northIsIn = True, showPlot = False)
for year in range(1950,2023):
    get_transport(year=year,model='RVC0',transectname='N2',x1=N2_x1,x2=N2_x2,y=N2_y,\
                     trac='DIC',northIsIn = True, showPlot = False)
for year in range(1950,2023):
    get_transport(year=year,model='RVC0',transectname='N3',x1=N3_x1,x2=N3_x2,y=N3_y,\
                     trac='DIC',northIsIn = True, showPlot = False)
    
    
for year in range(1950,2023):
    get_transport(year=year,model='RVD0',transectname='S1',x1=S1_x1,x2=S1_x2,y=S1_y,\
                     trac='DIC',northIsIn = True, showPlot = False)
for year in range(1950,2023):
    get_transport(year=year,model='RVD0',transectname='N1',x1=N1_x1,x2=N1_x2,y=N1_y,\
                     trac='DIC',northIsIn = True, showPlot = False)
for year in range(1950,2023):
    get_transport(year=year,model='RVD0',transectname='N2',x1=N2_x1,x2=N2_x2,y=N2_y,\
                     trac='DIC',northIsIn = True, showPlot = False)
for year in range(1950,2023):
    get_transport(year=year,model='RVD0',transectname='N3',x1=N3_x1,x2=N3_x2,y=N3_y,\
                     trac='DIC',northIsIn = True, showPlot = False)