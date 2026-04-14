import marimo

__generated_with = "0.23.1"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _():
    import xarray as xr
    import glob
    import matplotlib.pyplot as plt
    import numpy as np
    import warnings
    warnings.filterwarnings('ignore')
    return glob, np, plt, xr


@app.cell(hide_code=True)
def _(xr):
    mask = xr.open_dataset('/gpfs/home/mep22dku/scratch/SOZONE/UTILS/mesh_mask3pt6_nicedims.nc')
    return (mask,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## look at the ice concentration

    this is a standard model run from our restarts -- it looks like we do an aggressive hosing right at the beginning and have a jump in sea ice concentration. not clear why.
    """)
    return


@app.cell(hide_code=True)
def _(glob, mask, np, plt, xr):
    # Load data for both models
    files_T1 = sorted(glob.glob('/gpfs/afm/greenocean/software/runs/TOM12_TJ_50T1/ORCA2_1m_192*_icemod.nc') +
                      glob.glob('/gpfs/afm/greenocean/software/runs/TOM12_TJ_50T1/ORCA2_1m_193*_icemod.nc'))

    files_T1_gridT = sorted(glob.glob('/gpfs/afm/greenocean/software/runs/TOM12_TJ_50T1/ORCA2_1m_192*_grid_T.nc') +
                             glob.glob('/gpfs/afm/greenocean/software/runs/TOM12_TJ_50T1/ORCA2_1m_193*_grid_T.nc'))

    ds_T1 = xr.concat([xr.open_dataset(f) for f in files_T1], dim='time_counter')
    ds_T1_gridT = xr.concat([xr.open_dataset(f) for f in files_T1_gridT], dim='time_counter')

    # Weighted averages
    avg_T1 = (ds_T1.ice_pres * mask.csize).sum(dim=['y', 'x']) / mask.csize.sum(dim=['y', 'x'])
    avg_T1_sos = (ds_T1_gridT.sos * mask.csize).sum(dim=['y', 'x']) / mask.csize.sum(dim=['y', 'x'])

    years3 = 1920 + np.arange(len(avg_T1)) / 12

    # Plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

    ax1.plot(years3, avg_T1.values, label='TOM12_TJ_50T1')
    ax1.legend()
    ax1.set_title('Weighted average ice presence 1920-1940')
    ax1.set_ylabel('Ice presence')

    ax2.plot(years3, avg_T1_sos.values, label='TOM12_TJ_50T1', color='tab:blue')
    ax2.legend()
    ax2.set_title('Weighted average sea surface salinity 1920-1940')
    ax2.set_ylabel('SSS (psu)')
    ax2.set_xlabel('Year')

    plt.tight_layout()
    plt.show()
    return


if __name__ == "__main__":
    app.run()
