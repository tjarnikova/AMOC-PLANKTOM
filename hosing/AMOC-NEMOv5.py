import marimo

__generated_with = "0.23.1"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Visualization of AMOC and MLD in NEMOv5

    -
    - notes on TOM12_PT_KF51, monitor in `/gpfs/scratch/avd22gnu/ModelRuns/NEMO5/monitor/TOM12_PT_KF51/`

    To calculate streamfunction, need to rename dimensions, saved in `/gpfs/data/greenocean/users/mep22dku/clims/TOM12_PT_KF51`
    streamfunction calculation for NEMOv5:
    - interactive mode
    - `cd /gpfs/data/greenocean/software/resources/CDFTOOLS`
    `module add netcdf/4.7.4/gcc gcc/11.1.0`
    - check that meshmasks point appropriately, if not see `link_tommask_v5.sh`
    overturning call is in `calcover_nemo5test.sh` which calls cdfmoc
    """)
    return


@app.cell(hide_code=True)
def _():
    import marimo as mo
    import pandas as pd
    import xarray as xr
    import numpy as np
    import glob
    import os
    import matplotlib.pyplot as plt
    from pathlib import Path

    return Path, glob, mo, np, os, pd, plt, xr


@app.cell(hide_code=True)
def _(glob, os, xr):
    ex = False
    if ex:
        model_path = "/gpfs/scratch/avd22gnu/ModelRuns/NEMO5/TOM12_PT_KF51"
        model_name = os.path.basename(model_path)
        out_base = f"/gpfs/data/greenocean/users/mep22dku/clims/{model_name}"
        os.makedirs(out_base, exist_ok=True)
    
        for year in range(1750, 2101):
            pattern = os.path.join(model_path, f"ORCA2_7d_{year}0101_{year}1231_grid_V.nc")
            matches = glob.glob(pattern)
    
            if not matches:
                print(f"{year} file not found, skipping")
                continue
    
            fpath = matches[0]
            fname = os.path.basename(fpath)
    
            ds = xr.open_dataset(fpath)
    
            drop_vars = [v for v in ["x_grid_T", "y_grid_T", "nav_lat_grid_T", "nav_lon_grid_T", "vos", "tauvo"] if v in ds]
            ds = ds.drop_vars(drop_vars)
    
            rename_map = {}
            if "x_grid_V" in ds.dims: rename_map["x_grid_V"] = "x"
            if "y_grid_V" in ds.dims: rename_map["y_grid_V"] = "y"
            if "nav_lat_grid_V" in ds: rename_map["nav_lat_grid_V"] = "nav_lat"
            if "nav_lon_grid_V" in ds: rename_map["nav_lon_grid_V"] = "nav_lon"
            ds = ds.rename(rename_map)
    
            ds.attrs["history"] = "made in /gpfs/home/mep22dku/scratch/AMOC-PLANKTOM/hosing/AMOC-NEMOv5.py"
    
            out_path = os.path.join(out_base, fname)
            ds.to_netcdf(out_path)
            ds.close()
    
            print(f"saved {fname}")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## sample AMOC streamfunction -- declines in time
    """)
    return


@app.cell(hide_code=True)
def _(mo, np, plt, xr):
    def plot_moc(year):
        fpath = f"/gpfs/data/greenocean/software/resources/CDFTOOLS/MOCresults/TOM12_PT_KF51_7d_{year}0101_{year}1231_MOC.nc"
        ds = xr.open_dataset(fpath)
        moc = ds["zomsfatl"].mean(dim="time_counter").squeeze()
        depth = ds["depthw"].values
        lat = ds["nav_lat"].values
        if lat.ndim == 2:
            lat = lat[:, 0]

        fig, ax = plt.subplots(figsize=(10, 5))
        cf = ax.contourf(lat, depth, moc.values, levels=np.linspace(-20, 20, 41),
                         cmap="RdBu_r", vmin=-15, vmax=15, extend="both")
        ax.contour(lat, depth, moc.values, levels=np.linspace(-20, 20, 41),
                   colors="k", linewidths=0.3, alpha=0.4)
        plt.colorbar(cf, ax=ax, label="Sv")
        ax.set_xlabel("Latitude (°N)")
        ax.set_ylabel("Depth (m)")
        ax.set_title(f"Atlantic MOC streamfunction — TOM12_PT_KF51 ({year}, annual mean)")
        plt.tight_layout()
        return fig

    mo.vstack([plot_moc(y) for y in [1950, 1970, 2000]])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### sample sf at 26.1
    """)
    return


@app.cell(hide_code=True)
def _(mo, np, plt, xr):
    def plot_moc_26n(year):
        fpath = f"/gpfs/data/greenocean/software/resources/CDFTOOLS/MOCresults/TOM12_PT_KF51_7d_{year}0101_{year}1231_MOC.nc"
        ds = xr.open_dataset(fpath)
        moc = ds["zomsfatl"].mean(dim="time_counter").squeeze()
        depth = ds["depthw"].values
        lat = ds["nav_lat"].values.flatten()

        y_26n = int(np.abs(lat - 26.5).argmin())
        profile = moc.isel(y=y_26n).values

        fig, ax = plt.subplots(figsize=(4, 6))
        ax.plot(profile, depth)
        ax.axvline(0, color="k", linewidth=0.7)
        ax.set_xlabel("Sv")
        ax.set_ylabel("Depth (m)")
        ax.set_title(f"MOC at {lat[y_26n]:.1f}°N — {year}")

        plt.tight_layout()
        return fig


    mo.vstack([plot_moc_26n(y) for y in [1950, 1970, 2000]])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### SF timeseries - below 500m
    """)
    return


@app.cell(hide_code=True)
def _(pd, plt, xr):
    def plot_amoc(model="TOM12_PT_KF51", yrst=1920, yrend=2024):
        base = f"/gpfs/data/greenocean/users/mep22dku/clims/{model}"

        ds_26n   = xr.open_dataset(f"{base}/{model}_AMOC_26N_{yrst}_{yrend}.nc")
        ds_north = xr.open_dataset(f"{base}/{model}_AMOC_north_{yrst}_{yrend}.nc")

        t_26n   = ds_26n.time_counter.values
        t_north = ds_north.time_counter.values

        amoc_26n_raw   = pd.Series(ds_26n["AMOC_26N"].values,     index=t_26n)
        amoc_north_raw = pd.Series(ds_north["AMOC_north"].values,  index=t_north)
        amoc_26n_roll   = amoc_26n_raw.rolling(52, center=True).mean()
        amoc_north_roll = amoc_north_raw.rolling(52, center=True).mean()

        fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)

        for ax, raw, roll, label in zip(
            axes,
            [amoc_26n_raw, amoc_north_raw],
            [amoc_26n_roll, amoc_north_roll],
            ["26.5°N (below 500 m)", "Max north of 0°N (below 500 m)"]
        ):
            ax.plot(raw.index,  raw.values,  color="steelblue", linewidth=0.5, alpha=0.4, label="weekly")
            ax.plot(roll.index, roll.values, color="steelblue", linewidth=1.5, label="1-year rolling mean")
            ax.axhline(0, color="k", linewidth=0.5, linestyle="--")
            ax.set_ylim(0, 25)
            ax.set_ylabel("Sv")
            ax.set_title(label)
            ax.legend(fontsize=8)

        fig.suptitle(f"AMOC — {model} ({yrst}–{yrend})")
        plt.tight_layout()

        ds_26n.close()
        ds_north.close()
        return fig

    plot_amoc()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### annual-mean MLD in lab sea and nordic seas
    """)
    return


@app.cell
def _(mo, plt, xr):
    def plot_mld(year=1950):
        ds = xr.open_dataset(f"/gpfs/scratch/avd22gnu/ModelRuns/NEMO5/TOM12_PT_KF51/ORCA2_7d_{year}0101_{year}1231_grid_T.nc")
        mld = ds["mldr10_1"].mean(dim="time_counter")
        mld = mld.where(mld != 0)
    
        regions = {
            'Labrador': {'x': (110, 130), 'y': (105, 125), 'color': 'red'},
            'Nordic':   {'x': (130, 158), 'y': (120, 145), 'color': 'blue'},
        }

        fig, ax = plt.subplots(figsize=(10, 6))
        w = ax.pcolormesh(mld.values, cmap="viridis", vmin = 0, vmax = 500)

        for name, r in regions.items():
            x0, x1 = r['x']; y0, y1 = r['y']
            ax.add_patch(plt.Rectangle((x0, y0), x1-x0, y1-y0,
                         edgecolor=r['color'], facecolor='none', linewidth=2, label=name))
        

        ax.set_xlim(80, 165)
        ax.set_ylim(90, 148)
        ax.legend()
        ax.set_title(f"Mean MLD — {year}")
        plt.tight_layout()
        plt.colorbar(w)
        return fig



    mo.vstack([plot_mld(y) for y in [1950, 1970, 2000]])
    return


@app.cell(hide_code=True)
def _(Path, np, pd, xr):
    def get_mld_march(model="TOM12_PT_KF51"):
        base_in  = f"/gpfs/scratch/avd22gnu/ModelRuns/NEMO5/{model}"
        base_out = Path(f"/gpfs/data/greenocean/users/mep22dku/clims/{model}")
        base_out.mkdir(parents=True, exist_ok=True)

        regions = {
            'Labrador': {'x': (110, 130), 'y': (105, 125)},
            'Nordic':   {'x': (130, 158), 'y': (120, 145)},
        }

        results = {r: {'mean': [], 'max': [], 'time': []} for r in regions}

        for year in range(1750, 2101):
            fpath = f"{base_in}/ORCA2_7d_{year}0101_{year}1231_grid_T.nc"
            try:
                ds = xr.open_dataset(fpath)
            except FileNotFoundError:
                print(f"{year} not found, skipping")
                continue

            times = pd.to_datetime([pd.Timestamp(t.isoformat()) for t in ds.time_counter.values])
            march = np.where([t.month == 3 for t in times])[0]
            if len(march) == 0:
                print(f"{year} no March timesteps, skipping")
                ds.close()
                continue

            mld = ds["mldr10_1"].isel(time_counter=march).where(lambda x: x != 0)

            for name, r in regions.items():
                box = mld.isel(x_grid_T_2D_inner=slice(*r['x']), y_grid_T_2D_inner=slice(*r['y'])).mean(dim='time_counter')
                results[name]['mean'].append(float(box.mean()))
                results[name]['max'].append(float(box.max()))
                results[name]['time'].append(pd.Timestamp(f"{year}-03-15"))

            ds.close()
            print(f"processed {year}")

        for name in regions:
            if not results[name]['time']:
                continue
            r = results[name]
            ds_out = xr.Dataset({
                'mld_mean': xr.DataArray(r['mean'], coords={'time': r['time']}, dims='time',
                                         attrs={'long_name': f'Mean March MLD in {name} box', 'units': 'm'}),
                'mld_max':  xr.DataArray(r['max'],  coords={'time': r['time']}, dims='time',
                                         attrs={'long_name': f'Max March MLD in {name} box',  'units': 'm'}),
            }, attrs={'model': model, 'region': name, 'made_in': '/gpfs/home/mep22dku/scratch/AMOC-PLANKTOM/hosing/AMOC-NEMOv5.py'})
            outpath = base_out / f"{model}_MLD_march_{name}.nc"
            ds_out.to_netcdf(outpath)
            print(f"saved {outpath}")

    ext = False; 
    if ext: get_mld_march()
    return


@app.cell(hide_code=True)
def _(pd, plt, xr):
    def plot_mld_march(model="TOM12_PT_KF51"):
        base = f"/gpfs/data/greenocean/users/mep22dku/clims/{model}"
        regions = ['Labrador', 'Nordic']

        fig, axes = plt.subplots(2, 2, figsize=(14, 7), sharex=True)

        for col, name in enumerate(regions):
            ds = xr.open_dataset(f"{base}/{model}_MLD_march_{name}.nc")
            t  = ds.time.values
            for row, var in enumerate(['mld_mean', 'mld_max']):
                ax  = axes[row, col]
                raw  = pd.Series(ds[var].values, index=t)
                roll = raw.rolling(10, center=True).mean()
                ax.plot(t, raw.values,  color="steelblue", linewidth=0.5, alpha=0.4, label="yearly")
                ax.plot(t, roll.values, color="steelblue", linewidth=1.5, label="10-yr rolling")
                ax.set_title(f"{name} — {'Mean' if var == 'mld_mean' else 'Max'} March MLD")
                ax.set_ylabel("m")
                ax.legend(fontsize=8)
            ds.close()

        fig.suptitle(f"March MLD — {model}")
        plt.tight_layout()
        return fig

    plot_mld_march()
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
