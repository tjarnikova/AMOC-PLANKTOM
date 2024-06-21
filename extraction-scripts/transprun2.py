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
# for row in rows:
#     for year in range(1950,2022):
#         tran.get_transports(year,model,row,trac,northIsIn=True, showPlot=False, units = 'kmol/s')
model = 'RVD0'
trac = 'DIC'
rows = [73,90,100,109]
for row in rows:
    for year in range(1950,2022):
        tran.get_transports(year,model,row,trac,northIsIn=True, showPlot=False, units = 'kmol/s')

model = 'RVD0'
trac = 'NO3'
for row in rows:
    for year in range(1950,2022):
        tran.get_transports(year,model,row,trac,northIsIn=True, showPlot=False, units = 'kmol/s')

model = 'RVD0'
trac = 'Fer'
for row in rows:
    for year in range(1950,2022):
        tran.get_transports(year,model,row,trac,northIsIn=True, showPlot=False, units = 'kmol/s')

model = 'RVD0'
trac = 'Si'
for row in rows:
    for year in range(1950,2022):
        tran.get_transports(year,model,row,trac,northIsIn=True, showPlot=False, units = 'kmol/s')
