# TopoSCALE: ERA5 plugin
#
# === DESCRIPTION ==============================================================
#
#	This plugin 
#		* reads  monthly ERA5 data files (MARS efficiency)
#		* concats to single file timeseries
#		* converts values to toposcale standard
#		* writes single parameter files
#		*
#		*
#		*
#		*
# === COPYRIGHT AND LICENCE ====================================================
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# === CONTRIBUTIONS ====================================================
#
# === NOTES ============================================================

# === REQUIREMENTS =====================================================

# === DEPENDENCIES =====================================================
#	
#	* CDO: sudo apt install cdo
# ======================================================================
import subprocess
import numpy as np
import netCDF4 as nc

#1. concat monthly files to single using CDO 
#2. conversions
#2. write out single var files

class EraCat(object):
	"""
	Concats monthly era (interim or 5) files by some keyword *grepStr*

	Args:
		wd: directory of era monthly datafiles
		grepStr: "PLEVEL" or "SURF"
	Example:
		wd=/home/eraDat/
		grepStr= "PLEVEL"
		eraCat("/home/joel/mnt/myserver/sim/wfj_era5/eraDat/", "PLEVEL")
	"""

	def __init__(self, wd, grepStr):
		cmd     = ("cdo -b F64 -f nc2 mergetime " 
				+ wd 
				+  grepStr
				+ "* " 
				+ wd 
				+"/"
				+grepStr
				+".nc")
		subprocess.check_output(cmd, shell = "TRUE")

class Plev(object):
	"""
	Makes a plev object which is array of all pressure level variables 
	processed to standard toposcale units
	
	Args:
		fp: filepath to concatenated PLEVEL.nc file
	Example:
		p=Plev(fp)
		vars=p.vars()
	"""

	def __init__(self,fp, mylat, mylon):
	  	self.fp = fp
	  	self.vars = []
	  	self.mylat=mylat
	  	self.mylon=mylon
	  	self.g          = 9.80665 #m s-2
        self.absZero    = 273.15


	def getVarNames(self):
		"""	# returns list of variables excluding dimension names time, 
		lon,lat, level"""
		f = nc.Dataset(self.fp)
		
		for v in f.variables:
				if (v != ( u'time') 
					and v != ( u'latitude') 
					and v != ( u'longitude') 
					and v != ( u'level')):
			 		self.vars.append(v)
		#return self.vars

	
	def getVar(self, var):
		"""extract variables (remains an nc object)"""
		f = nc.Dataset(self.fp)
		self.myvar = f.variables[var]
		#return myvar

	
	def extractCgc(self,var):
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
		self.var = f.variables[var][ latli , lonli,:,: ] 
		#return mysub

		
	def addVar(self,varname,dat):
		""" rename attribute"""
		setattr(self, varname, dat)

	# create elevation (m) from geopotential	
	def convZ(self):
		self.z = self.z/self.g
		
	# creat wind speed and wind direction from u and V vectors	
	def convUV(self):	
class Plev_interp(Plev):
	"""
	Makes a plev object which is array of all variables processed to 
	standard toposcale units and interpolated to x, y, z 
	
	Args:
		fp: filepath to concatenated PLEVEL.nc file
	Example:
		p=Plev(fp)
		vars=p.vars()
	"""


	# interp variable and cgc (now np array) !!NOT FINISHED!! Perhaps a 
	# different class as here includes z interpolation
	def interpCgc(self,var):
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
		self.var = f.variables[var][ latli , lonli,:,: ] 
		#return mysub

class Surf(Plev):
	"""
	Makes a plev object which is array of all surface variables 
	processed to standard toposcale units
	
	Args:
		fp: filepath to concatenated PLEVEL.nc file
	Example:
		p=Plev(fp)
		vars=p.vars()
	"""
	
	def extractCgc(self,var):
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
		self.var = f.variables[var][ latli , lonli,: ] 
		#return mysub




#https://confluence.ecmwf.int/pages/viewpage.action?pageId=104241513





