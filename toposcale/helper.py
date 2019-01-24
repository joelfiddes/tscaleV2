"""helper methods """

from osgeo import ogr
import os



class Bunch:
	""" Class to create a Bunch
	whenever you want to group a few variables:
	http://code.activestate.com/recipes/52308-the-simple-but-handy-collector-of-a-bunch-of-named/?in=user-97991 """
	def __init__(self, **kwds):
		self.__dict__.update(kwds)


def get_shp_extent(shpPath):
	"""Get extent of shapefile. 

	Args: 
		shpPath: full path to shapefile

	Returns:
		Extent of shapefile as tuple [lonW,lonE, latS, latN]
		
	"""
	driver = ogr.GetDriverByName("ESRI Shapefile")
	dataSource = driver.Open(shpPath, 0)
	layer = dataSource.GetLayer()
	return(layer.GetExtent())

	#file = ogr.Open(shpPath)
	#shape = file.GetLayer(0)
	##first feature of the shapefile
	#feature = shape.GetFeature(0)
	#first = feature.ExportToJson()
	#print first # (GeoJSON format)
	#{"geometry": {"type": "LineString", "coordinates": [[0.0, 0.0], [25.0, 10.0], [50.0, 50.0]]}, "type": "Feature", "properties": {"FID": 0.0}, "id": 0}



