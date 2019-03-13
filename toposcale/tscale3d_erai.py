#!/usr/bin/env python

"""tscale3d.py


Example:



Attributes:

Todo:
    * remove dependency on dem when running "points"
	* input dem directly as netcdf or convert from tif	



"""

import recapp_erai as rc
import time
#import myplot
import os
from os import path, remove
import numpy as np

start_time = time.time()

mode ='points' #'grid'


if mode=="grid":
	# input filenames
	dir_data = './data'
	dataImport = rc.rawData(dir_data)

	dem_ascii = dir_data +'/dem1km.asc'
	dem_ncdf     = dir_data+'/dem1km.nc'
	if path.exists(dem_ncdf):
	    	remove(dem_ncdf)
	#convert ASCII DEM to netcdf format.
	dataImport.ascii2ncdf(dem_ascii, dem_ncdf)





pl   = './data/ecmwf_erai_pl_m_151201_151205.nc' #dataImport.plf_get()    # pressure level air temperature

# Format of stations for which to provide time series
# ['name':'siteName','lat':latNumber, 'lon':lonNumber, 'ele':eleNumber]             
mystations=[{'name':'COV','lat': 46.41801198, 'lon': 9.821232448, 'ele': 3350.5},
          {'name':'SAM','lat': 46.52639523, 'lon': 9.878944266, 'ele': 1756.2}]






var="Temperature"

if mode=="points":
	# station case
	ds = rc.tscale3dPl( pl)




	# point example
	out_xyz_dem, lats, lons, shape, names= ds.demGrid(stations=mystations)


	#out_xyz_ori = ds.geoGrid()
	#out_xyz_sur = ds.surGrid(lats, lons,stations=None)
	#surTa = ds.surTa(0, out_xyz_sur)


	# init grid stack
	xdim=shape[0]
	sa_vec = np.zeros((xdim))

	for timestep in range(0,20):
		print(timestep)
		gridT,gridZ,gridLat,gridLon=ds.gridValue(var,timestep)
		t_interp, z_interp = ds.inLevelInterp(gridT,gridZ, gridLat,gridLon,out_xyz_dem)
		#pl_sa = ds.fast1d(t_interp, z_interp, out_xyz_sur)
		pl_obs = ds.fast1d(t_interp, z_interp, out_xyz_dem)

		
		sa_vec=np.column_stack((sa_vec,pl_obs))
			




	print(" %f minutes for total run" % round((time.time()/60 - start_time/60),2) )





# Benchmark
#38809 1km cells
#6hr data
#3.6 mins/year/variable

#==============================================================================
if mode=="grid":
	# station case
	ds = rc.tscale3dPl( pl, dem_ncdf)


	out_xyz_dem, lats, lons, shape = ds.demGrid()


	#out_xyz_ori = ds.geoGrid()
	#out_xyz_sur = ds.surGrid(lats, lons,stations=None)
	#surTa = ds.surTa(0, out_xyz_sur)


	# init grid stack
	xdim=shape[0]
	ydim=shape[1]
	sa_out = np.zeros((xdim,ydim))

	for timestep in range(0,20):
		print(timestep)
		gridT,gridZ,gridLat,gridLon=ds.gridValue(var,timestep)
		t_interp, z_interp = ds.inLevelInterp(gridT,gridZ, gridLat,gridLon,out_xyz_dem)
		#pl_sa = ds.fast1d(t_interp, z_interp, out_xyz_sur)
		pl_obs = ds.fast1d(t_interp, z_interp, out_xyz_dem)

		sa_grid = pl_obs.reshape(xdim,ydim)
		sa_out = np.dstack((sa_out, sa_grid))

			
	# drop init blank layer

	a =sa_out[:,:,1:]



	print(" %f minutes for total run" % round((time.time()/60 - start_time/60),2) )


	# mean
	l = a.mean(axis=(2))

	#plot
	myplot.main(l)


