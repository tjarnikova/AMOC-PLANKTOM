#put python script here

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from cmocean import cm
import time
import pandas as pd

import transporter as tran
# model = 'RVA0'
# trac = 'DIC'
# for row in range(60,110):
#     for year in range(1950,2022):
#         tran.get_transports(year,model,row,trac,northIsIn=True, showPlot=False, units = 'kmol/s')
# model = 'RVC0'
# trac = 'DIC'
# for row in range(60,110):
#     for year in range(1950,2022):
#         tran.get_transports(year,model,row,trac,northIsIn=True, showPlot=False, units = 'kmol/s')


# model = 'RVB0'
# trac = 'DIC'
# for row in range(60,110):
#     for year in range(1950,2022):
#         tran.get_transports(year,model,row,trac,northIsIn=True, showPlot=False, units = 'kmol/s')

for row in range(73,110):
    for year in range(1950,2022):
        tran.get_transports_FNAT(year,row,'DIC',northIsIn=True, showPlot=False, units = 'kmol/s')