#!/usr/bin/env python

"""tscale3d_era5.py

Downscaling in 3d ERA5 pressure level data.

This module downscale ERA5 pl data air temp, rel hum, u, v and produces 2D grids
(based on high res DEM) or point based timeseries (based on csv file)

Example:
	python tscale3d_era5.py grid t

Attributes:
	ARG1: "mode": grid or points
	ARG2: "var": t (air temp) ,r (relhum),u (windU) or v (windV)

Todo:
    * remove dependency on dem when running "points"
	* input dem directly as netcdf or convert from tif
* make robust tif 	
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
#import myplot
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
	pl   = wdir+ '/forcing/PLEV.nc' #dataImport.plf_get()    # pressure level air temperature
	stationsfile= wdir+'/points.csv'
        
	if mode=="grid":
		dem_ncdf     =  wdir+"/predictors/ele.nc"

#============ POINTS RUN =======================================================
	if mode=="points" or mode=="point": # catches common typo
		if os.path.isfile(stationsfile) == False:
			print("No points.csv found!")

		mystations=pd.read_csv(stationsfile)

		# station case
		ds = rc.tscale3dPl( pl=pl)
		timesteps = ds.pl.variables['time'][:].size
		out_xyz_dem, lats, lons, shape, names= ds.demGrid(stations=mystations)

		# init grid stack
		xdim=shape[0]
		sa_vec = np.zeros((xdim))

		for timestep in range(starti, endi):
			#print(str(round(float(timestep)/float(timesteps)*100,0))+ "% done")
			gridT,gridZ,gridLat,gridLon=ds.gridValue(var,timestep)
			t_interp, z_interp = ds.inLevelInterp(gridT,gridZ, gridLat,gridLon,out_xyz_dem)
			pl_obs = ds.fast1d(t_interp, z_interp, out_xyz_dem)
			sa_vec=np.column_stack((sa_vec,pl_obs))
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
		ds = rc.tscale3dPl( pl=pl, dem=dem_ncdf)
		timesteps = ds.pl.variables['time'][:].size
		out_xyz_dem, lats, lons, shape = ds.demGrid()

		# init grid stack
		xdim=shape[0]
		ydim=shape[1]
		sa_out = np.zeros((xdim,ydim))

		for timestep in range(starti, endi):
			print(str(round(float(timestep)/float(timesteps)*100,0))+ "% done")
			gridT,gridZ,gridLat,gridLon=ds.gridValue(var,timestep)
			t_interp, z_interp = ds.inLevelInterp(gridT,gridZ, gridLat,gridLon,out_xyz_dem)
			pl_obs = ds.fast1d(t_interp, z_interp, out_xyz_dem)
			sa_grid = pl_obs.reshape(xdim,ydim)
			sa_grid = np.flip(sa_grid, 0) 
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

