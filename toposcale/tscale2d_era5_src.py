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
import os
from os import path, remove
import numpy as np
import pandas as pd
import sys
from osgeo import gdal
from osgeo import gdal_array
from osgeo import osr
		
# Global stuff
#var="v"
#mode ='grid' #'grid'

# Args
#wdir='/home/joel/sim/cci_perm/cci_test_ensemble/'
#var="t2m"
#mode ='points'
#starti=0
#endi=10
#member=1 
#dataset="EDA" 

#wdir='/home/joel/sim/cci_perm/cci_test_ensemble/'
#var="t2m"
#mode ='grid'
#starti=0
#endi=10
#member=1 
#dataset="EDA" 

#wdir='/home/joel/sim/cci_perm/cci_test/'
#var="t2m"
#mode ='points'
#starti=0
#endi=10
#member=1 
#dataset="HRES" 

#wdir='/home/joel/sim/cci_perm/cci_test/'
#var="t2m"
#mode ='grid'
#starti=0
#endi=10
#member=1 
#dataset="HRES"


def main(wdir, mode, var, starti, endi, dataset, member=None):

	start_time = time.time()
	sa   =  wdir+'/forcing/SURF.nc' #dataImport.plf_get()    # pressure level air temperature
	stationsfile= wdir+'/points.csv'


	if mode=="grid":
		dem_ncdf     = wdir+'/predictors/ele.nc'#dir_data+'/dem1km.nc'


#============ POINTS RUN =======================================================
	if mode=="points" or mode=="point":
		if os.path.isfile(stationsfile) == False:
			print("No points.csv found!")

		mystations=pd.read_csv(stationsfile)


		# Set up class object
		if dataset=="HRES":
			ds = rc.t3d( sa = sa, pl=None, dem=None)
		if dataset=="EDA":
			ds = rc.t3d_eda( sa = sa, pl=None, dem=None)


		timesteps = ds.sa.variables['time'][:].size
		out_xyz_dem, lats, lons, shape, names= ds.demGrid(stations=mystations)

		# init grid stack
		xdim=shape[0]
		sa_vec = np.zeros((xdim))


		for timestep in range(starti, endi):
			#print(str(round(float(timestep)/float(endi-starti)*100,0))+ "% done")
			if dataset=="HRES":
                       		sa_t = ds.surVarPoint(timestep, mystations, var)
			if dataset=="EDA":
                       		sa_t = ds.surVarPoint(timestep, mystations, var, member)
			sa_vec=np.column_stack((sa_vec,sa_t))
                
		# drop init row
		sa_vec =sa_vec[:,1:]

		# export dataframe
		tofile='False'
		if tofile=="True":
			df=pd.DataFrame(sa_vec.T)
			df.to_csv("points_"+var+".csv")
		
		print(" %f minutes for total run" % round((time.time()/60 - start_time/60),2) )
		return sa_vec.T
	# Benchmark
	#38809 1km cells
	#6hr data
	#3.6 mins/year/variable

#============ GRID RUN =======================================================
	if mode=="grid":
	
		# Set up class object
		if dataset=="HRES":
			ds = rc.t3d( sa = sa, pl=None, dem=dem_ncdf)
		if dataset=="EDA":
			ds = rc.t3d_eda( sa = sa, pl=None, dem=dem_ncdf)

		timesteps = ds.sa.variables['time'][:].size
		out_xyz_dem, lats, lons, shape = ds.demGrid()

		# init grid stack
		xdim=shape[0]
		ydim=shape[1]
		sa_out = np.zeros((xdim,ydim))
	
		for timestep in range(starti, endi):
			#print(str(round(float(timestep)/float(timesteps)*100,0))+ "% done")

			if dataset=="HRES":
                       		sa_t = ds.surVarGrid(timestep, lats, lons, var)
			if dataset=="EDA":
                       		sa_t = ds.surVarGrid(timestep, lats, lons, var, member)

			sa_grid = sa_t.reshape(xdim,ydim)
			sa_out = np.dstack((sa_out, sa_grid))
			
		# drop init blank layer
		a =sa_out[:,:,1:]

		print(" %f minutes for total run" % round((time.time()/60 - start_time/60),2) )

		# for testing
		plot='False'
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

	main(wdir, mode, var, starti, endi, dataset, member)

