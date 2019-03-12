#!/usr/bin/env python

"""tscale2d_era5.py

2D Interpolation of ERA5 surface data (no surface effects)

This module interpolates ERA5 surface data ssrd, strd, p from coarse ERA5 grid 
to a finescale DEM or set of points (points.csv).


Example:
	python tscale2d_era5.py grid t

Attributes:
	ARG1: "mode": grid or points
	ARG2: "var": t2m (air temp) ,'d2m','ssrd', strd, p
	ARG3: test "True" limit of 20 timesteps "False"= all timesteps
Todo:
	
Notes:

points.csv:
id,lon,lat,ele,slp,asp,sky
1,6.988,46.48,1992,22.07,105.6,1
2,6.994,46.48,1716,3.42,347.1,0.97

only first 4 cols required by this script. Others required for surface 
variables (sw,lw).

"""

import recapp_era5 as rc
import time
import myplot
import os
from os import path, remove
import numpy as np
import pandas as pd
import sys

# Global stuff
#var="v"
#mode ='grid' #'grid'


def main(wdir, mode, var, starti, endi):

	start_time = time.time()
	sa   =  wdir+'/forcing/SURF.nc' #dataImport.plf_get()    # pressure level air temperature
	stationsfile= wdir+'/points.csv'


	if mode=="grid":
		dem_ncdf     = wdir+'/predictos/ele.nc'#dir_data+'/dem1km.nc'


#============ POINTS RUN =======================================================
	if mode=="points" or mode=="point":
		if os.path.isfile(stationsfile) == False:
			print("No points.csv found!")

		mystations=pd.read_csv(stationsfile)

		# station case
		ds = rc.tscale3dPl( sa = sa, pl=None, dem=None)
		timesteps = ds.sa.variables['time'][:].size
		out_xyz_dem, lats, lons, shape, names= ds.demGrid(stations=mystations)

		# init grid stack
		xdim=shape[0]
		sa_vec = np.zeros((xdim))


		for timestep in range(starti, endi):
			#print(str(round(float(timestep)/float(timesteps)*100,0))+ "% done")
                        sa_t = ds.surTaPoint(timestep, mystations, var)
			sa_vec=np.column_stack((sa_vec,sa_t))
                
		# drop init row
		sa_vec =sa_vec[:,1:]

		# export dataframe
		tofile='False'
		if tofile=="True":
			df=pd.DataFrame(sa_vec.T)
			df.to_csv("TSP_"+var+".csv")
		
		print(" %f minutes for total run" % round((time.time()/60 - start_time/60),2) )
		return sa_vec.T
	# Benchmark
	#38809 1km cells
	#6hr data
	#3.6 mins/year/variable

#============ GRID RUN =======================================================
	if mode=="grid":
		ds = rc.tscale3dPl( sa=sa, dem=dem_ncdf, pl=None)
		timesteps = ds.sa.variables['time'][:].size
		out_xyz_dem, lats, lons, shape = ds.demGrid()

		# init grid stack
		xdim=shape[0]
		ydim=shape[1]
		sa_out = np.zeros((xdim,ydim))
	
		for timestep in range(starti, endi):
			print(str(round(float(timestep)/float(timesteps)*100,0))+ "% done")
			sa_t = ds.surTaGrid(timestep, lats, lons, var)
			sa_grid = sa_t.reshape(xdim,ydim)
			sa_out = np.dstack((sa_out, sa_grid))
			
		# drop init blank layer
		a =sa_out[:,:,1:]

		print(" %f minutes for total run" % round((time.time()/60 - start_time/60),2) )

		# mean
		#l = a.mean(axis=(2))
		#plot
		#myplot.main(l)
		return a
#===============================================================================
#	Calling Main
#===============================================================================
if __name__ == '__main__':
	
	mode = sys.argv[1]
	var = sys.argv[2]
	test=sys.argv[3]
	main(mode, var, test)

