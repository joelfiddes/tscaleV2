""" fetch_srtm

Fetch SRTM DEM based on input shpefile domain (points or polygon). Create 
merged DEM.

Account required https://urs.earthdata.nasa.gov/profile. Authentication 
controlled by file ~/.netrc

~/.netrc content:
	machine urs.earthdata.nasa.gov
	login USER
	password PWD


"""
#import subprocess
import netCDF4 as nc
#====================================================================
# PARAMETERS/ARGS
#====================================================================
#args = commandArgs(trailingOnly=TRUE)
wd=config["main"]["wd"]
 # args[1] #"/home/joel/sim/testDem" #args[1]
src = config["main"]["srcdir"]
#= args[2] #"/home/joel/data/DEM/srtm" 
rst = "era5grid30.tif" #"/home/joel/src/topoMAPP/dat/eraigrid75.tif"
shp= shpPath #"/home/joel/data/GCOS/wfj_poly.shp"

#====================================================================
# Credentials/ access
#====================================================================
SERVICE=unlist(strsplit(readLines("~/.netrc")[[1]]," "))[2]
print(paste0('using credentials for: ', SERVICE))
USER=unlist(strsplit(readLines("~/.netrc")[[2]]," "))[2]
PWD=unlist(strsplit(readLines("~/.netrc")[[3]]," "))[2]


demDir = wd + '/spatial/'
forcing_grid = wd + '/forcing/PLEVEL.nc'

f = nc.Dataset(forcing_grid)
#subprocess.check_output("gdalwarp -of GTiff -cutline " + shp + " -cl area_of_interest  -crop_to_cutline " +forcing_grid+  " out.tiff")


