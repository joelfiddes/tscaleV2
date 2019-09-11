#===============================================================================

require(raster)
require(ncdf4)
require(rgeos)

args = commandArgs(trailingOnly=TRUE)
wd=args[1]


setwd(paste0(wd,"/forcing/"))
surfnc = paste0(wd,"/forcing/SURF.nc")
plevnc = paste0(wd,"/forcing/PLEV.nc")

basin=shapefile(paste0(wd,"/basins/basins.shp"))
centroids = gCentroid(basin,byid=T)
mylons=centroids@coords[,1]
mylats=centroids@coords[,2]

for (i in 1:length(mylons)){
	mylon=mylons[i]
	mylat=mylats[i]
	outfile=paste0("g",i,"_surf.nc")
	cmd = paste0("cdo intgridbil,lon=",mylon,"_lat=",mylat," ",surfnc," "  ,outfile)
	system(cmd)
}


for (i in 1:length(mylons)){
	mylon=mylons[i]
	mylat=mylats[i]
	outfile=paste0("g",i,"_plev.nc")
	cmd = paste0("cdo intgridbil,lon=",mylon,"_lat=",mylat," ",plevnc," "  ,outfile)
	system(cmd)
}
