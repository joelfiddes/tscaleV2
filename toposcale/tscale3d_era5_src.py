

"""tscale3d_era5_EDA.py

Downscaling in 3d ERA5 EDA pressure level data.

This module downscale ERA5 pl data air temp, rel hum, u, v and produces 2D grids
(based on high res DEM) or point based timeseries (based on csv file). It is called 
by tscale3D.py.

Example:
	python tscale3d_era5_src.py '/home/joel/sim/cci_perm/cci_test_ensemble/' 't' 'points' 0 10 1 "EDA"

Attributes:
	ARG1: "wdir": workdir	
	ARG2: "mode": grid or points
	ARG3: "var": t (air temp) ,r (relhum),u (windU) or v (windV)
	ARG4: "starti": start index (numeric index for nc file)
	ARG5: "endi": end index (numeric index for nc file)
	ARG6: dataset: "HRES" (ERA5) or "EDA" (ensemble ERA5)
	ARG7: ensemble member [1-10] or not declared (HRES)

	
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
import os
from os import path, remove
import numpy as np
import pandas as pd
import sys
#from osgeo import gdal
#from osgeo import gdal_array
#from osgeo import osr


## Test
#wdir='/home/joel/sim/cci_perm/cci_test_ensemble/'
#var="t"
#mode ='points'
#starti=0
#endi=10
#member=1 
#dataset="EDA" 

#wdir='/home/joel/sim/cci_perm/cci_test_ensemble/'
#var="t"
#mode ='grid'
#starti=0
#endi=10
#member=1 
#dataset="EDA" 

#wdir='/home/joel/sim/cci_perm/cci_test/'
#var="t"
#mode ='points'
#starti=0
#endi=10
#member=1 
#dataset="HRES" 

#wdir='/home/joel/sim/cci_perm/cci_test/'
#var="t"
#mode ='grid'
#starti=0
#endi=10
#member=1 
#dataset="HRES" 

def main(wdir, mode, var, starti, endi,dataset, member=None):
	start_time = time.time()
	


	pl   = wdir+ '/forcing/PLEV.nc' #dataImport.plf_get()    # pressure level air temperature
	stationsfile= wdir+'/listpoints.txt'
        
	if mode=="grid":
		dem_ncdf     =  wdir+"/predictors/ele.nc"

#============ POINTS RUN =======================================================
	if mode=="points" or mode=="point": # catches common typo
		if os.path.isfile(stationsfile) == False:
			print("No points.csv found!")

		mystations=pd.read_csv(stationsfile)

		# Set up class object
		if dataset=="HRES":
			ds = rc.t3d( pl=pl)
		if dataset=="EDA":
			ds = rc.t3d_eda( pl=pl)


		timesteps = ds.pl.variables['time'][:].size
		out_xyz_dem, lats, lons, shape, names= ds.demGrid(stations=mystations)

		# init grid stack
		xdim=shape[0]
		sa_vec = np.zeros(xdim)

		for timestep in range(starti, endi):
			#print(str(round(float(timestep)/float(endi-starti)*100,0))+ "% done")
			gridT,gridZ,gridLat,gridLon=ds.gridValue(var,timestep)

			if dataset=="HRES":
				t_interp, z_interp = ds.inLevelInterp(gridT,gridZ, gridLat,gridLon,out_xyz_dem)
			if dataset=="EDA":
				t_interp, z_interp = ds.inLevelInterp(gridT,gridZ, gridLat,gridLon,out_xyz_dem,member)
#			if dataset=="EDA":
#				t_interp, z_interp = ds.inLevelInterp(gridT[member,:,:,:],gridZ[member,:,:,:], gridLat,gridLon,out_xyz_dem)
			pl_obs = ds.fast1d(t_interp, z_interp, out_xyz_dem)
			sa_vec=np.column_stack((sa_vec,pl_obs))
                # drop init row
		sa_vec =sa_vec[:,1:]
		
		# export dataframe
		tofile='False'
		if tofile=="True":
			df=pd.DataFrame(sa_vec.T)
			df.to_csv(wdir+"points_"+var+".csv")

		print(" %f minutes for total run of " + var % round((time.time()/60 - start_time/60),2) )
		return sa_vec.T

	# Benchmark
	#38809 1km cells
	#6hr data
	#3.6 mins/year/variable

#============ GRID RUN =======================================================
	if mode=="grid":

		# Set up class object
		if dataset=="HRES":
			ds = rc.t3d( pl=pl, dem=dem_ncdf)
		if dataset=="EDA":
			ds = rc.t3d_eda( pl=pl, dem=dem_ncdf)


		timesteps = ds.pl.variables['time'][:].size
		out_xyz_dem, lats, lons, shape = ds.demGrid()

		# init grid stack
		xdim=shape[0]
		ydim=shape[1]
		sa_out = np.zeros((xdim,ydim))

		for timestep in range(starti, endi):
			#print(str(round(float(timestep)/float(endi-starti)*100,0))+ "% done")
			gridT,gridZ,gridLat,gridLon=ds.gridValue(var,timestep)

			if dataset=="HRES":
				t_interp, z_interp = ds.inLevelInterp(gridT,gridZ, gridLat,gridLon,out_xyz_dem)
			if dataset=="EDA":
				t_interp, z_interp = ds.inLevelInterp(gridT,gridZ, gridLat,gridLon,out_xyz_dem,member)

			pl_obs = ds.fast1d(t_interp, z_interp, out_xyz_dem)
			sa_grid = pl_obs.reshape(xdim,ydim)
			sa_grid = np.flip(sa_grid, 0) 
			sa_out = np.dstack((sa_out, sa_grid))
			
		# drop init blank layer
		a =sa_out[:,:,1:]

		print(" %f minutes for run" % round((time.time()/60 - start_time/60),2) )

		# for testing
		plot='True'
		if plot=="True":
			import myplot
			l = a.mean(axis=(2))
			myplot.main(l)

		writegrid='False'
		if writegrid=="True":
		# https://gis.stackexchange.com/questions/37238/writing-numpy-array-to-raster-file
			

			

			array = l    
			lat = out_xyz_dem[:,0].reshape(l.shape)
			lon = out_xyz_dem[:,1].reshape(l.shape)

			xmin,ymin,xmax,ymax = [lon.min(),lat.min(),lon.max(),lat.max()]
			nrows,ncols = np.shape(array)
			xres = (xmax-xmin)/float(ncols)
			yres = (ymax-ymin)/float(nrows)
			geotransform=(xmin,xres,0,ymax,0, -yres)   

			output_raster = gdal.GetDriverByName('GTiff').Create('myraster.tif',ncols, nrows, 1 ,gdal.GDT_Float32)# Open the file
			output_raster.GetRasterBand(1).WriteArray( array )  # Writes my array to the raster
			output_raster.SetGeoTransform(geotransform)# Specify its coordinates
			srs = osr.SpatialReference()# Establish its coordinate encoding
			srs.ImportFromEPSG(4326)   # This one specifies WGS84 lat long.
			output_raster.SetProjection(srs.ExportToWkt())# Exports the coordinate system 
			output_raster = None

	
		return a
#===============================================================================
#	Calling Main
#===============================================================================
if __name__ == '__main__':
	
	wdir = sys.argv[1]
	mode = sys.argv[2]
	var=sys.argv[3]
	starti=sys.argv[4]
	endi=sys.argv[5]
	dataset=sys.argv[6]
	member=sys.argv[7]
	main(wdir, mode,var, starti, endi, dataset, member)




