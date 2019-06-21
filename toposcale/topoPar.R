
require(raster)
require(horizon)

# This function computes asp,slp, svf  grids used by SW routine in tscale3D mode=grids only
wdir ="/home/joel/sim/tscale3Dbig" #args[2]
angles = as.numeric(args[2])
dist = as.numeric(args[3])
dem=raster(paste0(wdir,"/predictors/ele.nc"))


if(!file.exists(paste0(wdir, "/predictors/svf.tif"))){

	s <- svf(dem, nAngles=angles, maxDist=dist, ll=TRUE)
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
writeRaster(round(slp,0), "slp.nc", overwrite=TRUE) #write and reduce precision
writeRaster(round(asp,0), "asp.nc", overwrite=TRUE) #write and reduce precision
writeRaster(round(s,2), "svf.nc", overwrite=TRUE)


