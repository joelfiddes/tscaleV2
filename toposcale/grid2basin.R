#domain should be bigger than shape - actually not
require(raster)
require(ncdf4)
require(rgeos)
args = commandArgs(trailingOnly=TRUE)
wd=args[1]

basinshp=paste0(wd,'basins/basins.shp')
surffile=paste0(wd,"/forcing/SURF.nc")
plevfile=paste0(wd,"/forcing/PLEV.nc")
#surffile='SURF.nc'
#fp=inDir+"/PLEV.nc"
#lpfile=home + "/listpoints.txt"

s = stack(surffile)
fp = s[[1]]

values(fp)<-1:ncell(s)
shp=shapefile(basinshp)

# ids of era5 grids that belong to each polygon
basins =extract(fp,shp)



#=============================================
#	SURFACE VARIABLES
#=============================================

#read in meteo
nc = nc_open(surffile)

#=== t2m ===
t2m_nc =ncvar_get(nc, 't2m')
z_nc =ncvar_get(nc, 'z')
d2m_nc =ncvar_get(nc, 'd2m')
strd_nc =ncvar_get(nc, 'strd')
ssrd_nc =ncvar_get(nc, 'ssrd')
tp_nc =ncvar_get(nc, 'tp')
tisr_nc =ncvar_get(nc, 'tisr')
time_nc =ncvar_get(nc, 'time')

#loop through basins
t2m_mat=c()
z_mat=c()
d2m_mat=c()
strd_mat=c()
ssrd_mat=c()
tp_mat=c()
tisr_mat=c()

for ( b in basins){

# get grids that exist in basin
print(b)
cellIndex=na.omit(unlist(b))
#cellIndex=b

# get row/col of each grid
myrow = ceiling(cellIndex/ncol(s))
mycol = cellIndex-((myrow-1)*33)

# extract subset of variables corresponding to those grids
t2m_sub=t2m_nc[mycol,myrow,]
z_sub=z_nc[mycol,myrow,]
d2m_sub=d2m_nc[mycol,myrow,]
strd_sub=strd_nc[mycol,myrow,]
ssrd_sub=ssrd_nc[mycol,myrow,]
tp_sub=tp_nc[mycol,myrow,]
tisr_sub=tisr_nc[mycol,myrow,]



if ( length(dim(t2m_sub))){
	# the mean timeseries for basin
	t2m_sub_mean = apply(t2m_sub, MARGIN=c(3), FUN='mean')
	z_sub_mean = apply(z_sub, MARGIN=c(3), FUN='mean')
	d2m_sub_mean = apply(d2m_sub, MARGIN=c(3), FUN='mean')
	strd_sub_mean = apply(strd_sub, MARGIN=c(3), FUN='mean')
	ssrd_sub_mean = apply(ssrd_sub, MARGIN=c(3), FUN='mean')
	tp_sub_mean = apply(tp_sub, MARGIN=c(3), FUN='mean')
	tisr_sub_mean = apply(tisr_sub, MARGIN=c(3), FUN='mean')

	# get some stats on the distribution
	#t2m_sub_sd = apply(t2m_sub, MARGIN=c(3), FUN='sd')

	#compute p lapse rate here
	
	}else{
t2m_sub_mean <- t2m_sub
z_sub_mean <- z_sub
d2m_sub_mean <- d2m_sub
strd_sub_mean <- strd_sub
ssrd_sub_mean <- ssrd_sub
tp_sub_mean <- tp_sub
tisr_sub_mean <- tisr_sub
}


t2m_mat=c(t2m_mat, t2m_sub_mean)
z_mat=c(z_mat, z_sub_mean)
d2m_mat=c(d2m_mat, d2m_sub_mean)
strd_mat=c(strd_mat, strd_sub_mean)
ssrd_mat=c(ssrd_mat, ssrd_sub_mean)
tp_mat=c(tp_mat, tp_sub_mean)
tisr_mat=c(tisr_mat, tisr_sub_mean)


}
rm(t2m_nc)
gc()
rm(z_nc)
gc()
rm(d2m_nc)
gc()
rm(strd_nc)
gc()
rm(ssrd_nc)
gc()
rm(tp_nc)
gc()
rm(tisr_nc)
gc()
rm(t2m_sub)
gc()
rm(z_sub)
gc()
rm(d2m_sub)
gc()
rm(strd_sub)
gc()
rm(ssrd_sub)
gc()
rm(tp_sub)
gc()
rm(tisr_sub)
gc()

t2m = matrix(t2m_mat, nrow=length(basins), ncol=length(t2m_sub_mean), byrow=T)
zs = matrix(z_mat, nrow=length(basins), ncol=length(z_sub_mean), byrow=T)
d2m = matrix(d2m_mat, nrow=length(basins), ncol=length(d2m_sub_mean), byrow=T)
strd = matrix(strd_mat, nrow=length(basins), ncol=length(strd_sub_mean), byrow=T)
ssrd = matrix(ssrd_mat, nrow=length(basins), ncol=length(ssrd_sub_mean), byrow=T)
tp = matrix(tp_mat, nrow=length(basins), ncol=length(tp_sub_mean), byrow=T)
tisr = matrix(tisr_mat, nrow=length(basins), ncol=length(tisr_sub_mean), byrow=T)

rm(t2m_mat)
gc()
rm(z_mat)
gc()
rm(d2m_mat)
gc()
rm(strd_mat)
gc()
rm(ssrd_mat)
gc()
rm(tp_mat)
gc()
rm(tisr_mat)
gc()
# get centroids
centroids = gCentroid(shp,byid=T)
mylon=centroids@coords[,1]
mylat=centroids@coords[,2]




#=============================================
#	writeNCDF
#=============================================

writeNcdf =function(varname,units,ncfname,title,dlname,data,time_nc,wd){
	fillvalue <- -9999
	basindim<- ncdim_def("basin", "degrees_north", as.double(1:length(basins)))
	tdim <- ncdim_def("time", nc$dim$time$units, as.double(time_nc))
	# define variables
	tmp.def <- ncvar_def(varname, units, list(basindim, tdim), fillvalue, 
		             dlname, prec = "single")
	ncout <- nc_create(paste0(wd,"/forcing/",ncfname), list(tmp.def), force_v4 = TRUE)
	# put the array
	ncvar_put(ncout, tmp.def, data)
	ncatt_put(ncout, "basin", "axis", "Y")
	# add global attributes
	ncatt_put(ncout, 0, "title", title)
	# close the file, writing data to disk
	nc_close(ncout)
	}


varname="t2m"
units="k"
ncfname <- "t2m.nc"
title <- "t2m at basin centroids"
dlname <- "2m Surface temperature"
data=t2m
writeNcdf(varname,units,ncfname,title,dlname,data,time_nc,wd)
rm(t2m)
gc()
varname="z"
units="m**2 s**-2"
ncfname <- "zs.nc"
title <- "geopotential at basin centroids"
dlname <- "Surface geopotential"
data=zs
writeNcdf(varname,units,ncfname,title,dlname,data,time_nc,wd)
rm(z)
gc()
varname="d2m"
units="K"
ncfname <- "d2m.nc"
title <- "2 metre dewpoint temperature at basin centroids"
dlname <- "2 metre dewpoint temperature"
data=d2m
writeNcdf(varname,units,ncfname,title,dlname,data,time_nc,wd)
rm(d2m)
gc()
varname="strd"
units="J m**-2"
ncfname <- "strd.nc"
title <- "Surface thermal radiation downwards at basin centroids"
dlname <- "Surface thermal radiation downwards"
data=strd
writeNcdf(varname,units,ncfname,title,dlname,data,time_nc,wd)
rm(strd)
gc()
varname="ssrd"
units="J m**-2"
ncfname <- "ssrd.nc"
title <- "Surface solar radiation downwards at basin centroids"
dlname <- "Surface solar radiation downwards"
data=ssrd
writeNcdf(varname,units,ncfname,title,dlname,data,time_nc,wd)
rm(ssrd)
gc()

tp[tp<0]<-0
tp[which(is.infinite(as.vector(tp)))]<-0 # finds and remove Inf's why these exist?! must be TP correct method

varname="tp"
units="m"
ncfname <- "tp.nc"
title <- "total precipitation at basin centroids"
dlname <- "total precipitation"
data=tp
writeNcdf(varname,units,ncfname,title,dlname,data,time_nc,wd)
rm(tp)
gc()
varname="tisr"
units="J m**-2"
ncfname <- "tisr.nc"
title <- "TOA incident solar radiation at basin centroids"
dlname <- "TOA incident solar radiation"
data=tisr
writeNcdf(varname,units,ncfname,title,dlname,data,time_nc,wd)
rm(tisr)
gc()



#=============================================
#	PRESSURE LEVELS VARIABLES
#=============================================

#read in meteo
nc = nc_open(plevfile)
lev =ncvar_get(nc, 'level')

#=== t2m ===
z_nc =ncvar_get(nc, 'z')
t_nc =ncvar_get(nc, 't')
u_nc =ncvar_get(nc, 'u')
v_nc =ncvar_get(nc, 'v')
r_nc =ncvar_get(nc, 'r')



#loop through basins
z_mat=c()
t_mat=c()
u_mat=c()
v_mat=c()
r_mat=c()


for ( b in basins){

# get grids that exist in basin
print(b)
cellIndex=b

# get row/col of each grid
myrow = ceiling(cellIndex/ncol(s))
mycol = cellIndex-((myrow-1)*33)

# extract subset of variables corresponding to those grids
z_sub=z_nc[mycol,myrow,,]
t_sub=t_nc[mycol,myrow,,]
u_sub=u_nc[mycol,myrow,,]
v_sub=v_nc[mycol,myrow,,]
r_sub=r_nc[mycol,myrow,,]



if ( length(dim(z_sub))>2){
	# the mean timeseries for basin
	z_sub_mean = as.vector(t(apply(z_sub, MARGIN=c(3,4), FUN='mean')))
	t_sub_mean = as.vector(t(apply(t_sub, MARGIN=c(3,4), FUN='mean')))
	u_sub_mean = as.vector(t(apply(u_sub, MARGIN=c(3,4), FUN='mean')))
	v_sub_mean = as.vector(t(apply(v_sub, MARGIN=c(3,4), FUN='mean')))
	r_sub_mean = as.vector(t(apply(r_sub, MARGIN=c(3,4), FUN='mean')))


	# get some stats on the distribution
	#t2m_sub_sd = apply(t2m_sub, MARGIN=c(3), FUN='sd')

	#compute p lapse rate here
	
	}else{
z_sub_mean <- as.vector(z_sub)
t_sub_mean <- as.vector(t_sub)
u_sub_mean <- as.vector(u_sub)
v_sub_mean <- as.vector(v_sub)
r_sub_mean <- as.vector(r_sub)

}



z_mat=c(z_mat, z_sub_mean)
t_mat=c(t_mat, t_sub_mean)
u_mat=c(u_mat, u_sub_mean)
v_mat=c(v_mat, v_sub_mean)
r_mat=c(r_mat, r_sub_mean)

}

rm(z_nc)
gc()
rm(t_nc)
gc()
rm(u_nc)
gc()
rm(v_nc)
gc()
rm(r_nc)
gc()
rm(z_sub)
gc()
rm(t_sub)
gc()
rm(u_sub)
gc()
rm(v_sub)
gc()
rm(r_sub)
gc()


zp = matrix(z_mat, nrow=length(basins), ncol=length(z_sub_mean), byrow=T)
t = matrix(t_mat, nrow=length(basins), ncol=length(t_sub_mean), byrow=T)
u = matrix(u_mat, nrow=length(basins), ncol=length(u_sub_mean), byrow=T)
v = matrix(v_mat, nrow=length(basins), ncol=length(v_sub_mean), byrow=T)
r = matrix(r_mat, nrow=length(basins), ncol=length(r_sub_mean), byrow=T)

# reshape to 3D
zp_array = array(zp, dim=c(length(basins), length(time_nc), length(lev)))
rm(zp)
gc()

t_array = array(t, dim=c(length(basins), length(time_nc), length(lev)))
rm(t)
gc()
u_array = array(u, dim=c(length(basins), length(time_nc), length(lev)))
rm(u)
gc()
v_array = array(v, dim=c(length(basins), length(time_nc), length(lev)))
rm(v)
gc()
r_array = array(r, dim=c(length(basins), length(time_nc), length(lev)))
rm(r)
gc()

# get centroids
centroids = gCentroid(shp,byid=T)
mylon=centroids@coords[,1]
mylat=centroids@coords[,2]




#=============================================
#	writeNCDF - plevel
#=============================================
writeNcdf =function(varname,units,ncfname,title,dlname,data,time_nc,wd,lev){
	fillvalue <- -9999
	basindim<- ncdim_def("basin", "degrees_north", as.double(1:length(basins)))
	plevdim<- ncdim_def("levels", "mb", as.double(lev))
	tdim <- ncdim_def("time", nc$dim$time$units, as.double(time_nc))
	# define variables
	tmp.def <- ncvar_def(varname, units, list(basindim, tdim, plevdim), fillvalue, 
		             dlname, prec = "single")
	ncout <- nc_create(paste0(wd,"/forcing/",ncfname), list(tmp.def), force_v4 = TRUE)
	# put the array
	ncvar_put(ncout, tmp.def, data)
	ncatt_put(ncout, "basin", "axis", "Y")
	ncatt_put(ncout, "levels", "axis", "Z")
	# add global attributes
	ncatt_put(ncout, 0, "title", title)
	# close the file, writing data to disk
	nc_close(ncout)
	}


varname="z"
units="m**2 s**-2"
ncfname <- "zp.nc"
title <- "geopotential at basin centroids"
dlname <- "geopotential"
data=zp_array
writeNcdf(varname,units,ncfname,title,dlname,data,time_nc,wd,lev)

varname="t"
units="k"
ncfname <- "t.nc"
title <- "t at basin centroids"
dlname <- "temperature"
data=t_array
writeNcdf(varname,units,ncfname,title,dlname,data,time_nc,wd,lev)

varname="u"
units=" m s**-1"
ncfname <- "u.nc"
title <- "u"
dlname <- "u"
data=u_array
writeNcdf(varname,units,ncfname,title,dlname,data,time_nc,wd,lev)

varname="v"
units=" m s**-1"
ncfname <- "v.nc"
title <- "v"
dlname <- "v"
data=v_array
writeNcdf(varname,units,ncfname,title,dlname,data,time_nc,wd,lev)

varname="r"
units="%"
ncfname <- "r.nc"
title <- "Relative humidity"
dlname <- "Relative humidity"
data=r_array
writeNcdf(varname,units,ncfname,title,dlname,data,time_nc,wd,lev)




