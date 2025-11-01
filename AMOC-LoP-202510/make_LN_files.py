import xarray as xr
import numpy as np

def get_limiter(run='TOM12_TJ_LC00', year=1920):
    """
    Calculate limiting nutrients for different phytoplankton functional types.
    
    Nutrient codes:
        3 = Fe (Iron)
        4 = P (Phosphorus)
        5 = Si (Silica, diatoms only)
        6 = N (Nitrogen)
    
    Args:
        run: Model run identifier
        year: Year to process
    """
    # Load data
    tdir = f'/gpfs/data/greenocean/software/runs/{run}/'
    tfi = f'ORCA2_1m_{year}0101_{year}1231_diad_T.nc'
    w = xr.open_dataset(f'{tdir}/{tfi}')
    print(f'{run} {year}')
    
    # Load and prepare mesh mask
    tmesh = xr.open_dataset('/gpfs/data/greenocean/software/resources/regrid/mesh_mask3_6.nc')
    tm = tmesh.tmask.values
    
    # Broadcast mesh mask to match data dimensions [time, depth, lat, lon]
    tm_broad = np.broadcast_to(tm[0, :, :, :], (12, 31, 149, 182)).copy()
    
    # Process each phytoplankton functional type
    pfts = ['dia', 'mix', 'coc', 'pic', 'pha', 'fix']
    
    for pft in pfts:
        # Gather limitation variables
        limitation_vars = [
            w[f'lim3fe_{pft}'],   # Fe limitation
            w[f'lim4po4_{pft}'],  # P limitation
        ]
        
        # Diatoms also have silica limitation
        if pft == 'dia':
            limitation_vars.append(w[f'lim5si_{pft}'])  # Si limitation
        
        limitation_vars.append(w[f'lim6din_{pft}'])  # N limitation
        
        # Stack limitations and find minimum (most limiting nutrient)
        stacked = xr.concat(limitation_vars, dim='variable')
        stacked = stacked.where(stacked != 0, np.nan)
        stacked = stacked.fillna(np.inf)
        
        # Store minimum limitation value
        w[f'LV_{pft}'] = stacked.min(dim='variable', skipna=True)
        
        # Identify which nutrient is limiting
        min_indices = stacked.argmin(dim='variable', skipna=True)
        q2 = np.copy(min_indices).astype(float)
        q2[tm_broad == 0] = np.nan
        
        # Map indices to nutrient codes
        limiter_codes = (q2 + 1) * 10
        
        if pft == 'dia':
            # Diatoms: 0->Fe(3), 1->P(4), 2->Si(5), 3->N(6)
            mapping = {10: 3, 20: 4, 30: 5, 40: 6}
        else:
            # Other PFTs: 0->Fe(3), 1->P(4), 2->N(6)
            mapping = {10: 3, 20: 4, 30: 6}
        
        for old_val, new_val in mapping.items():
            limiter_codes[limiter_codes == old_val] = new_val
        
        limiter_codes[tm_broad == 0] = np.nan
        
        # Store limiting nutrient identity
        w[f'LN_{pft}'] = w[f'LV_{pft}']
        w[f'LN_{pft}'].data = limiter_codes
        
        # Add attributes
        w[f'LV_{pft}'].attrs['long_name'] = f'Minimum limitation value for {pft}'
        w[f'LN_{pft}'].attrs['long_name'] = f'Limiting nutrient for {pft}'
        w[f'LN_{pft}'].attrs['note'] = 'Nutrient codes: 3=Fe, 4=P, 5=Si (diatoms only), 6=N'
    
    # Select output variables
    output_vars = []
    for pft in pfts:
        output_vars.extend([f'LV_{pft}', f'LN_{pft}'])
    w_sel = w[output_vars]
    
    # Add global attributes
    w_sel.attrs['description'] = 'Limiting nutrient analysis for phytoplankton functional types'
    w_sel.attrs['nutrient_codes'] = '3=Fe, 4=P, 5=Si (diatoms only), 6=N'
    
    # Save to NetCDF
    output_file = f'{tdir}/ORCA2_1m_{year}0101_{year}1231_LN_T.nc'
    try:
        w_sel.to_netcdf(output_file)
        print(f'Saved {run} {year}:')
        print(f'{output_file}\n')
    except Exception as e:
        print(f'Failed to save {run} {year}: {e}\n')
    
    w_sel.close()
    w.close()


horse = True
if horse:
    
    runs = ['TOM12_TJ_LA50']

    for run in runs:
        for y in range(1940,2030):
            get_limiter(run,y)
            #except: print(f'no for {run} {y}')
