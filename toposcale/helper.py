"""helper methods """

#from osgeo import ogr
import os



class Bunch:
	""" Class to create a Bunch
	whenever you want to group a few variables:
	http://code.activestate.com/recipes/52308-the-simple-but-handy-collector-of-a-bunch-of-named/?in=user-97991 """
	def __init__(self, **kwds):
		self.__dict__.update(kwds)


#def get_shp_extent(shpPath):
#	"""Get extent of shapefile. 

#	Args: 
#		shpPath: full path to shapefile

#	Returns:
#		Extent of shapefile as tuple [lonW,lonE, latS, latN]
#		
#	"""
#	driver = ogr.GetDriverByName("ESRI Shapefile")
#	dataSource = driver.Open(shpPath, 0)
#	layer = dataSource.GetLayer()
#	return(layer.GetExtent())



#def get_nc_bbox(nc):
#	f = nc.Dataset(nc)
#	lats = f.variables['latitude'][:] 
#	lons = f.variables['longitude'][:]
#	cellRes= (lats.max()-lats.min())/ (lats.size-1)
#	lats.min() - cellRes



