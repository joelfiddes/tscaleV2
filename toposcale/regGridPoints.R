# creates regular grid of spatial points to tscale contained by shp
# this is used for example as input to SWAT
# https://stackoverflow.com/questions/41787313/how-to-create-a-grid-of-spatial-points

library(sp)
library(rgdal)
library(raster)

resol=0.1# cellsize in map units!
myshape <- shapefile('/home/joel/data/zisp/inputs/zarafshanBasin_cut.shp')
clip=TRUE

# check the CRS to know which map units are used
#proj4string(myshape)
# "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Create a grid of points within the bbox of the SpatialPolygonsDataFrame 
# myshape with decimal degrees as map units
grid <- makegrid(myshape, cellsize = resol) # cellsize in map units!

# grid is a data.frame. To change it to a spatial data set we have to
grid <- SpatialPoints(grid, proj4string = CRS(proj4string(myshape)))

if (clip == TRUE){
grid <- grid[myshape, ] #To extract only the points within the polygon, use `[` to subset the points based on the location like this:
}
plot(myshape)
plot(grid, pch = ".", add = T)
shapefile(grid , "swat_shp.shp")
