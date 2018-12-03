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
# === CONTRIBUTIONS ============================================================
#
# === NOTES ====================================================================

# === REQUIREMENTS =============================================================

# === DEPENDENCIES =============================================================
#	
#	* CDO: sudo apt install cdo
# ==============================================================================
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
		cmd     = "cdo -b F64 -f nc2 mergetime " + wd +  grepStr+ "* " + wd +"/"+grepStr+".nc"
		subprocess.check_output(cmd, shell = "TRUE")

class Plev(object):
	"""
	Makes a plev object which is array of all variables processed to standard toposcale units
	
	Args:
		fp: filepath to concatenated PLEVEL.nc file
	Example:
		p=Plev(fp)
		vars=p.vars()
	"""

	def __init__(self,fp):
	  	self.fp = fp

	# returns list of variables excluding dimension names time, lon,lat, level
	def getVarNames(self):
		f = nc.Dataset(self.fp)
		self.vars = []
		for v in f.variables:
				if v != ( u'time') and v != ( u'latitude') and v != ( u'longitude') and v != ( u'level'):
			 		self.vars.append(v)
		#return self.vars

	# extract variables	(still an nc object)
	def getVar(self, var):
		f = nc.Dataset(self.fp)
		self.myvar = f.variables[var]
		#return myvar

	# extract variable and cgc (now np array)
	def extractCgc(self,var, lon, lat):
		f = nc.Dataset(fp)

		latbounds = [ mylat , mylat ]
		lonbounds = [ mylon , mylon ] # degrees east ? 
		lats = f.variables['latitude'][:] 
		lons = f.variables['longitude'][:]

		# latitude lower and upper index
		latli = np.argmin( np.abs( lats - latbounds[0] ) )
		latui = np.argmin( np.abs( lats - latbounds[1] ) ) 

		# longitude lower and upper index
		lonli = np.argmin( np.abs( lons - lonbounds[0] ) )
		lonui = np.argmin( np.abs( lons - lonbounds[1] ) )  

		# Air (time, latitude, longitude) 
		mysub = f.variables[var][ latli , lonli,:,: ] 
		return mysub

	def addVar(self,varname,dat):
		#self.varname = dat
		setattr(self, varname, dat)




#class Surf(object)


#https://confluence.ecmwf.int/pages/viewpage.action?pageId=104241513





