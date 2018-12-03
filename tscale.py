
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

class cgcExtract(object):
	"""
	Return extracted cgc of interest

	Args:
		lon: cgc centre lon
		lat: cgc centre lat
		nc: data netcdf

	Example:
		lon=7.5
		lat=48.5
		nc= /home/joel/PLEV.nc

		cgcExtract(lon, lat, nc)

	""" 

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
    Return object for downscaling that has methods for interpolationg
    upper-air temperature and surface influences at surface level
    based on disaggregating coarse-grid reanalysis and dem.
    
    Args:
        dem: A required fine-scale dem in netcdf format
    
    Example:
        dem  = 'example_alps.nc'
        geop = 'alps_geop.nc'
        sa   = 'alps_sa_79_15.nc'
        pl   = 'alps_pl_79_15.nc'
        
        downscaling = downscaling(dem, geop, sa, pl)
        
    """