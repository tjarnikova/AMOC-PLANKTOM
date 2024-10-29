import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from cmocean import cm
import time
import pandas as pd


def get_atlantic(array, row_index, target_index = 130):
    row = array[row_index]
    if row[target_index] != 1:
        return None, None

    # Find the start of the segment of 0s
    start_index = target_index
    while start_index > 0 and row[start_index - 1] == 1:
        start_index -= 1

    # Find the end of the segment of 0s
    end_index = target_index
    while end_index < len(row) - 1 and row[end_index + 1] == 1:
        end_index += 1

    if row_index in [91,92]:
        start_index = 93

    return start_index, end_index


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
    csize_N1[:,x2+1:190]=0
    csize_N1v = np.copy(csize_N1)
    csize_N1v[np.where(csize_N1 < 1)] = np.nan
    
    return csize_N1v


def get_transports(year,model,row,trac,northIsIn=True, showPlot=False, units = 'kmol/s'):
    w1 = time.time()
    y = row
    tdi = '/gpfs/data/greenocean/software/runs/TOM12_TJ_'
    tfi = xr.open_dataset(f'{tdi}{model}/ORCA2_1m_{year}0101_{year}1231_ptrc_T.nc')
    vfi = xr.open_dataset(f'{tdi}{model}/ORCA2_1m_{year}0101_{year}1231_grid_V.nc')

    tmesh = xr.open_dataset('/gpfs/data/greenocean/software/resources/regrid/mesh_mask3_6.nc')

    ##size mask for transect only
    x1, x2 = get_atlantic(tmesh['tmaskutil'][0,:,:].values,row)
    cmask = get_mask(x1,x2+1,row)
    tmb = np.zeros([12,np.shape(cmask)[0],np.shape(cmask)[1]])

    for i in range(0,12):
        tmb[i,:,:] = cmask
    tmb1 = tmb[:,:,x1:x2+1]

    savenam = f'/gpfs/home/mep22dku/scratch/AMOC-PLANKTOM/transports/{model}_{year}_{trac}_row{row}_x{x1}-{x2}.nc'

    vo = vfi['vomecrty'][:,:,y,:].values
    transp = vo * tmb

    to1 =  tfi[trac][:,:,y,:].values
    to2 =  tfi[trac][:,:,y+1,:].values

    to = (to1+to2)/2
    totransp = transp * to

    if (northIsIn == False):
        transp = -1 * transp
        totransp = -1 * totransp


    if showPlot:
        fact = 1.1
        fig, axs = plt.subplots(1,5, figsize=(20*fact, 3*fact), facecolor='w', edgecolor='k')
        axs = axs.ravel()

        w = axs[0].pcolormesh(tmb[0,:,:], cmap = plt.cm.Spectral)
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

        tits = ['mask','vomecrty','transport (Sv)',trac,f'{trac} transport (kmol/s)']
        for i in range(0,5):
            axs[i].invert_yaxis()
            axs[i].set_title(tits[i])


    transp_res = np.zeros([12,31,3])
    totransp_res = np.zeros([12,31,3])

    for i in range(0,12):
        for k in range(0,31):

            tt = transp[i,0:k+1,:]
            ttpos = (np.nansum(tt[np.where(tt>=0)]))
            ttneg = (np.nansum(tt[np.where(tt<0)]))
            ttnet = np.nansum(tt)
            transp_res[i,k,0] = ttpos
            transp_res[i,k,1] = ttneg
            transp_res[i,k,2] = ttnet

            tt = totransp[i,0:k+1,:]
            ttpos = (np.nansum(tt[np.where(tt>=0)]))
            ttneg = (np.nansum(tt[np.where(tt<0)]))
            ttnet = np.nansum(tt)
            totransp_res[i,k,0] = ttpos
            totransp_res[i,k,1] = ttneg
            totransp_res[i,k,2] = ttnet

    times = pd.date_range(f"{year}/01/01",f"{year+1}/01/01",freq='MS',closed='left')

    data_vars = {'tracer_transports':(['time_counter', 'depthlevel', 'type'], totransp_res,
    {'units': units,
    'long_name':trac}),
                 'water_transports':(['time_counter', 'depthlevel', 'type'], transp_res,
    {'units': 'm3/s',
    'long_name':'water only'}),
    }
    # define coordinates
    coords = {'time_counter':  times,
    'depthlevel': (np.arange(0,31,1)),
    'type': (['pos','neg','net']),\
             }
    # define global attributes
    attrs = {'made in':'scratch/AMOC-PLANKTOM/transprun.py',
    'desc': 'lateral transports'
    }
    ds = xr.Dataset(data_vars=data_vars,
    coords=coords,
    attrs=attrs)


    w2 = time.time()
    print(f'{savenam} {w2-w1}')
    ds.to_netcdf(savenam)


def get_transports_FNAT(year,row,trac,northIsIn=True, showPlot=False, units = 'kmol/s'):
    
    ##for FNAT, which is RVD0 - RVB0
    w1 = time.time()
    y = row
    tdi = '/gpfs/data/greenocean/software/runs/TOM12_TJ_'
    tfi = xr.open_dataset(f'{tdi}RVD0/ORCA2_1m_{year}0101_{year}1231_ptrc_T.nc')
    tfi2 = xr.open_dataset(f'{tdi}RVB0/ORCA2_1m_{year}0101_{year}1231_ptrc_T.nc')
    vfi = xr.open_dataset(f'{tdi}RVD0/ORCA2_1m_{year}0101_{year}1231_grid_V.nc')

    tmesh = xr.open_dataset('/gpfs/data/greenocean/software/resources/regrid/mesh_mask3_6.nc')

    ##size mask for transect only
    x1, x2 = get_atlantic(tmesh['tmaskutil'][0,:,:].values,row)
    cmask = get_mask(x1,x2+1,row)
    tmb = np.zeros([12,np.shape(cmask)[0],np.shape(cmask)[1]])

    for i in range(0,12):
        tmb[i,:,:] = cmask
    tmb1 = tmb[:,:,x1:x2+1]

    savenam = f'/gpfs/home/mep22dku/scratch/AMOC-PLANKTOM/transports/FNAT-RVD0-RVB0_{year}_{trac}_row{row}_x{x1}-{x2}.nc'

    vo = vfi['vomecrty'][:,:,y,:].values
    transp = vo * tmb

    to1 =  tfi[trac][:,:,y,:].values
    to2 =  tfi[trac][:,:,y+1,:].values
    to = (to1+to2)/2

    too1 =  tfi2[trac][:,:,y,:].values
    too2 =  tfi2[trac][:,:,y+1,:].values
    too = (too1+too2)/2
    tores = to-too
    
    totransp = transp * tores

    if (northIsIn == False):
        transp = -1 * transp
        totransp = -1 * totransp


    if showPlot:
        fact = 1.1
        fig, axs = plt.subplots(1,5, figsize=(20*fact, 3*fact), facecolor='w', edgecolor='k')
        axs = axs.ravel()

        w = axs[0].pcolormesh(tmb[0,:,:], cmap = plt.cm.Spectral)
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

        tits = ['mask','vomecrty','transport (Sv)',trac,f'{trac} transport (kmol/s)']
        for i in range(0,5):
            axs[i].invert_yaxis()
            axs[i].set_title(tits[i])


    transp_res = np.zeros([12,31,3])
    totransp_res = np.zeros([12,31,3])

    for i in range(0,12):
        for k in range(0,31):

            tt = transp[i,0:k+1,:]
            ttpos = (np.nansum(tt[np.where(tt>=0)]))
            ttneg = (np.nansum(tt[np.where(tt<0)]))
            ttnet = np.nansum(tt)
            transp_res[i,k,0] = ttpos
            transp_res[i,k,1] = ttneg
            transp_res[i,k,2] = ttnet

            tt = totransp[i,0:k+1,:]
            ttpos = (np.nansum(tt[np.where(tt>=0)]))
            ttneg = (np.nansum(tt[np.where(tt<0)]))
            ttnet = np.nansum(tt)
            totransp_res[i,k,0] = ttpos
            totransp_res[i,k,1] = ttneg
            totransp_res[i,k,2] = ttnet

    times = pd.date_range(f"{year}/01/01",f"{year+1}/01/01",freq='MS',closed='left')

    data_vars = {'tracer_transports':(['time_counter', 'depthlevel', 'type'], totransp_res,
    {'units': units,
    'long_name':trac}),
                 'water_transports':(['time_counter', 'depthlevel', 'type'], transp_res,
    {'units': 'm3/s',
    'long_name':'water only'}),
    }
    # define coordinates
    coords = {'time_counter':  times,
    'depthlevel': (np.arange(0,31,1)),
    'type': (['pos','neg','net']),\
             }
    # define global attributes
    attrs = {'made in':'/gpfs/home/mep22dku/scratch/AMOC-PLANKTOM/extraction-scripts/tranpart.py using transporter.py',
    'desc': 'lateral transports'
    }
    ds = xr.Dataset(data_vars=data_vars,
    coords=coords,
    attrs=attrs)


    w2 = time.time()
    print(f'{savenam} {w2-w1}')
    ds.to_netcdf(savenam)