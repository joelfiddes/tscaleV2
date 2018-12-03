
# TopoSCALE: Topograhic based downscaling of atmospheric datasets
#
# === DESCRIPTION ==============================================================

# === COPYRIGHT AND LICENCE ====================================================
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# === CONTRIBUTIONS ============================================================
#
# === NOTES ====================================================================

# === REQUIREMENTS =============================================================
#	* PLEV.nc - all pressure level variables for entire domain (single/multiple 
#		CGCs) [time X grid cell X  level X variable]
#	* SURF.nc - all surface variables for entire domain (single/multiple 
#		CGCs) [time X grid cell X variable]

import numpy as np
import netCDF4 as nc


class cgcExtract(object):
	"""
	Return extracted cgc of interest

	Args:
		lon: cgc centre lon
		lat: cgc centre lat
		nc: data netcdf

	Example:
		mylon=9.3
		mylat=46.2
		nc= /home/joel/PLEV.nc

		cgcExtract(mylon, mylat, ncfile)

	""" 


	def showVars(ncfile):

		f = nc.Dataset(ncfile)
	def getVar(mylon, mylat, ncfile, var):

		f = nc.Dataset(ncfile)

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



class cgcInterp(object):
	"""
	Return interpolated cgc to point of interest

	Args:
		lon: POI lon
		lat: POI lat
		nc: data netcdf

	Example:
		lon=7.5
		lat=48.5
		nc= /home/joel/PLEV.nc

		cgcExtract(lon, lat, nc)
		
	""" 

class downscale(object):
    """
    Return object for downscaling that has methods for interpolating
	pressure level data based on a given elevation. If a 
    
    Args:
        dem: A required fine-scale dem in netcdf format
    
    Example:
        dem  = 'example_alps.nc'
        geop = 'alps_geop.nc'
        sa   = 'alps_sa_79_15.nc'
        pl   = 'alps_pl_79_15.nc'
        
        downscaling = downscaling(dem, geop, sa, pl)
        
    """
        def __init__(self, geop, sa, pl, dem = None):
        self.g    = 9.80665 #Gravitational acceleration [m/s2]
        self.geop = nc.Dataset(geop)
        self.sa   = nc.Dataset(sa)
        self.pl   = nc.Dataset(pl)
        if not (dem is None):
            self.dem  = nc.Dataset(dem)