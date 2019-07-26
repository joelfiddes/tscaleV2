args = commandArgs(trailingOnly=TRUE)
require(raster)
require(horizon) # https://cran.r-project.org/web/packages/horizon/horizon.pdf
				# R implementation of Dozier and Frew 1990 algorithms

# This function computes asp,slp, svf  grids used by SW routine in tscale3D mode=grids only and writes all as nc files
wdir =args[1] #"/home/joel/sim/tscale3Dbig" #args[2]
angles = as.numeric(args[2]) #number of sectors 8 good shown by sensitivity test
dist = as.numeric(args[3])	#in meteres 5000 is good shown by sensitivity test
dem=raster(paste0(wdir,"/predictors/ele.tif"))


if(!file.exists(paste0(wdir, "/predictors/svf.tif"))){

	s <- svf(dem, nAngles=angles, maxDist=dist, ll=TRUE)
	}

if(file.exists(paste0(wdir, "/predictors/svf.tif"))){

	s <- raster(paste0(wdir, "/predictors/svf.tif"))
	}
#====================================================================
# EXTRACT SLP/ASP
#================================================================= ==
slp=terrain(dem, opt="slope", unit="degrees", neighbors=8, filename='')
asp=terrain(dem, opt="aspect", unit="degrees", neighbors=8, filename='')

#====================================================================
# WRITE OUTPUTS
#====================================================================
setwd(paste0(wdir,'/predictors'))
writeRaster(round(dem,0), "ele.nc", overwrite=TRUE) #write and reduce precision
writeRaster(round(slp,0), "slp.nc", overwrite=TRUE) #write and reduce precision
writeRaster(round(asp,0), "asp.nc", overwrite=TRUE) #write and reduce precision
writeRaster(round(s,2), "svf.nc", overwrite=TRUE)


