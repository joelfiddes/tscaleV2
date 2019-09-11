#NCAR archive
#https://rda.ucar.edu/datasets/ds084.1/#!description

#Data
#https://nomads.ncep.noaa.gov:9090/dods/gfs_0p25

#tHIS ONE:
#https://stackoverflow.com/questions/52406826/load-selection-gfs-ensemble-opendap-data-into-memory-python

from netCDF4 import Dataset
import numpy as np
import pandas as pd
import xarray as xr
import time as tm

# Set time to download data from (this is always the 00UTC run of the present day)
time_year = str(tm.localtime()[0])
time_month = str(tm.localtime()[1])
time_day = str(tm.localtime()[2])

if len(time_month)== 1:
    time_month = '0' + time_month
if len(time_day)== 1:
    time_day = '0' + time_day

datestr = time_year + time_month + time_day
print('The run chosen is the 00 UTC run of ' + time_day + '-' + time_month + '-' + time_year)

# Define server information
#serverstring='http://nomads.ncep.noaa.gov:9090/dods/gens_bc/gens' + datestr + '/gep_all_00z'
 
 
serverstring='https://nomads.ncep.noaa.gov:9090/dods/gfs_0p25_1hr/gfs'+time_year+time_month+time_day+'/gfs_0p25_1hr_00z'
print(serverstring) # Load data
dataset = xr.open_dataset(serverstring)
time = dataset.variables['time']
lat = dataset.variables['lat'][:]
lon = dataset.variables['lon'][:] 
lev = dataset.variables['lev'][:] 

import myplot
tmax = dataset.variables['tmax2m'][:]
t =tmax[1,:,:]
myplot.main(t)



tmax = dataset.variables['tmax2m'][:]

ugrdprs 
vgrdprs
tmp2m
rhprs
pratesfc
tmp2m
