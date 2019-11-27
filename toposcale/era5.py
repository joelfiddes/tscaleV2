"""ERA5.py

ERA5 IO plugin

This module rquires input files:

	* PLEV.nc: 
		all pressure level variables for entire domain (single/multiple 
		CGCs) [time X grid cell X  level X variable]
	* SURF.nc: 
		all surface variables for entire domain (single/multiple 
		CGCs) [time X grid cell X variable]

This plugin contains methods to:

	* classes for both pressurelevel and surface objects 
	* extract coarse grid cell (CGC) based on longitude latitude of station 
		or CGC centre
	* convert values to toposcale standard
	* writes single parameter files
	* TopoSCALE standard input:

	   	- air temperature - K
		- precipitation - mmh*1
		- shortwave Wm**2
		- longwave Wm**2
		- wind - U and V vectors
		- time - iso
		- relative humidity

Example:
	Initialise new era5 instance::

	p=era5.Plev(fp, stat.lat, stat.lon)

Attributes:

Todo:
 
"""



#		*
#		*

import subprocess
import numpy as np
import netCDF4 as nc
import pandas as pd

#1. concat monthly files to single using CDO 
#2. conversions
#2. write out single var files


class Plev(object):
	"""
	Makes a plev object which is array of all pressure level variables 
	processed to standard toposcale units
	
	Args:
		fp: filepath to concatenated PLEVEL.nc file
		mylat: latitude of station or grid centre 
		mylon: longitude of station or gride centre
	Example:
		p=Plev(fp)
		varnames=p.varnames()
	"""
	def __init__(self,fp, mylat, mylon):
		self.fp = fp
		self.varnames = []
		self.mylat=mylat
		self.mylon=mylon



	def getVarNames(self):
		"""	# returns list of variables excluding dimension names time, 
		lon,lat, level"""
		f = nc.Dataset(self.fp)
		
		for v in f.variables:
				if (v != ( u'time') 
					and v != ( u'latitude') 
					and v != ( u'longitude') 
					and v != ( u'number') 
					and v != ( u'level')):
			 		self.varnames.append(v)
		#return self.varnames

	def getVar(self, var):
		"""extract variables (remains an nc object)"""
		f = nc.Dataset(self.fp)
		self.myvar = f.variables[var]
		#return myvar

	def extractCgc(self,var, startIndex, endIndex):
		"""extract variable and cgc (now np array)
		dimension order: time, level, lat,lon
		"""
		f = nc.Dataset(self.fp)

		latbounds = [ self.mylat , self.mylat ]
		lonbounds = [ self.mylon , self.mylon ] # degrees east ? 
		lats = f.variables['latitude'][:] 
		lons = f.variables['longitude'][:]

		# latitude lower  index
		latli = np.argmin( np.abs( lats - latbounds[0] ) )

		# longitude lower index
		lonli = np.argmin( np.abs( lons - lonbounds[0] ) )

		# subset
		self.var = f.variables[var][ startIndex:endIndex,: ,latli , lonli] 
		#return mysub

	def extractCgc5d(self,var, member, startIndex, endIndex):
		"""extract variable, ensemble memeber and cgc (now np array) from 5d
		nc files eg era5 ensemble (extra dimension ensemble 'number') 

		dimension order: time, ensemble, level, lat,lon
		"""
		f = nc.Dataset(self.fp)

		latbounds = [ self.mylat , self.mylat ]
		lonbounds = [ self.mylon , self.mylon ] # degrees east ? 
		lats = f.variables['latitude'][:] 
		lons = f.variables['longitude'][:]

		# latitude lower  index
		latli = np.argmin( np.abs( lats - latbounds[0] ) )

		# longitude lower index
		lonli = np.argmin( np.abs( lons - lonbounds[0] ) )
		
		# subset : memeber -1 convert from index 1-10 (input)to 0-9 (python)
		self.var = f.variables[var][ startIndex:endIndex,member-1 ,:,latli , lonli] 
		#return mysub


	def addVar(self,varname,dat):
		""" rename attribute"""
		setattr(self, varname, dat)

	def addTime(self):
		""" add time vector and convert to ISO 
			Return datetime objects given numeric time values. 
			The units of the numeric time values are described by the units 
			argument and the calendar keyword. The returned datetime objects 
			represent UTC with no time-zone offset, even if the specified 
			units contain a time-zone offset.

			calender options defined here:
			http://unidata.github.io/netcdf4-python/#netCDF4.num2date
			
		"""
		f = nc.Dataset(self.fp)
		self.nctime = f.variables['time']
		self.dtime = pd.to_datetime(nc.num2date(self.nctime[:],self.nctime.units, calendar="standard"))

	def addShape(self):
		""" adds two dimensions of time and levels """
		self.myshape = self.var.shape
	# def convZ(self):
	# 	""" create elevation (m) from geopotential """
	# 	self.z = self.z/self.g

	def plevels(self):
		f = nc.Dataset(self.fp)
		self.levels = f.variables['level'][:]

	


class Plev_interp(Plev):
	"""
	Makes a plev object which is array of all variables processed to 
	standard toposcale units and interpolated to x, y, z 
	
	Args:
		fp: filepath to concatenated PLEVEL.nc file
	Example:
		p=Plev(fp)
		varnames=p.varnames()
	"""

	def interpCgc(self,var):
		"""	# interp variable and cgc (now np array) !!NOT FINISHED!! 
		Perhaps a  different class as here includes z interpolation"""
		f = nc.Dataset(self.fp)

		latbounds = [ self.mylat , self.mylat ]
		lonbounds = [ self.mylon , self.mylon ] # degrees east ? 
		lats = f.variables['latitude'][:] 
		lons = f.variables['longitude'][:]

		# latitude index
		latli = np.argmin( np.abs( lats - latbounds[0] ) )

		# longitudeindex
		lonli = np.argmin( np.abs( lons - lonbounds[0] ) )

		# find neighbours 4 or 8?
		# dissaggregate to high res grid
		# interpolate in x and y

		"""
		Evaluate a simple example function on the points of a 3D grid:

		from scipy.interpolate import RegularGridInterpolator
		def f(x,y,z):
		    return 2 * x**3 + 3 * y**2 - z
		x = np.linspace(1, 4, 11)
		y = np.linspace(4, 7, 22)
		z = np.linspace(7, 9, 33)
		data = f(*np.meshgrid(x, y, z, indexing='ij', sparse=True))
		#data is now a 3D array with data[i,j,k] = f(x[i], y[j], z[k]). Next, define an interpolating function from this data:

		my_interpolating_function = RegularGridInterpolator((x, y, z), data)
		#Evaluate the interpolating function at the two points (x,y,z) = (2.1, 6.2, 8.3) and (3.3, 5.2, 7.1):

		pts = np.array([[2.1, 6.2, 8.3], [3.3, 5.2, 7.1]])
		my_interpolating_function(pts)

		#which is indeed a close approximation to [f(2.1, 6.2, 8.3), f(3.3, 5.2, 7.1)].
		"""

		# subset
		self.var = f.variables[var][ :,: ,latli , lonli] 
		#return mysub

class Surf(Plev):
	"""
	Makes a plev object which is array of all surface variables 
	processed to standard toposcale units
	
	Args:
		fp: filepath to concatenated PLEVEL.nc file
	Example:
		p=Plev(fp)
		varnames=p.varnames()
	"""
	
	def extractCgc(self,var,startIndex, endIndex):
		"""extract variable and cgc (now np array)"""
		f = nc.Dataset(self.fp)

		latbounds = [ self.mylat , self.mylat ]
		lonbounds = [ self.mylon , self.mylon ] # degrees east ? 
		lats = f.variables['latitude'][:] 
		lons = f.variables['longitude'][:]

		# latitude lower  index
		latli = np.argmin( np.abs( lats - latbounds[0] ) )

		# longitude lower index
		lonli = np.argmin( np.abs( lons - lonbounds[0] ) )

		# subset
		self.var = f.variables[var][ startIndex:endIndex ,latli , lonli] 
		#return mysub

	def extractCgc4d(self,var, member,startIndex, endIndex):
		"""extract variable, ensemble memeber and cgc (now np array)from 4d
		nc files eg era5 ensemble (extra dimension ensemble 'number') 

		dimension order: time, ensemble, lat,lon"""
		f = nc.Dataset(self.fp)

		latbounds = [ self.mylat , self.mylat ]
		lonbounds = [ self.mylon , self.mylon ] # degrees east ? 
		lats = f.variables['latitude'][:] 
		lons = f.variables['longitude'][:]

		# latitude lower  index
		latli = np.argmin( np.abs( lats - latbounds[0] ) )

		# longitude lower index
		lonli = np.argmin( np.abs( lons - lonbounds[0] ) )

		# subset : memeber -1 convert from index 1-10 (input)to 0-9 (python)
		self.var = f.variables[var][ startIndex:endIndex ,member-1, latli , lonli] 
		#return mysub

	def instRad(self, step):
		""" Convert SWin from accumulated quantities in J/m2 to 
		instantaneous W/m2 see: 
		https://confluence.ecmwf.int/pages/viewpage.action?pageId=104241513

		Args:
			step: timstep in seconds (era5=3600, ensemble=10800)

		Note: both EDA (ensemble 3h) and HRES (1h) are accumulated over the timestep
		and therefore treated here the same.
		https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
				"""
		self.strd = self.strd/step  
		self.ssrd = self.ssrd/step
		self.tisr = self.tisr/step 

	def tp2rate(self, step):
		""" convert tp from m/timestep (total accumulation over timestep) to rate in mm/h 

			Args:
				step: timstep in seconds (era5=3600, ensemble=10800)

			Note: both EDA (ensemble 3h) and HRES (1h) are accumulated over the timestep
			and therefore treated here the same.
			https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
		"""
		stepinhr=step/(60*60)
		self.prate = self.tp*1000 # convert m per timestep -> mm/hour [PRATE] use approx scaling method until monthly correction fixed. 
		# Although this is actually pretty accurate. tp is accumulated over hour so prate IS tp at coarser timstep (eg 6h) we introduce uncertainty
		self.psum = self.tp* stepinhr	*1000 # mm/ timestep  [PSUM] scale by stepinhr
	
	def gridEle(self):
		""" compute surface elevation of coarse grid"""
		self.gridEle = self.z/9.80665












