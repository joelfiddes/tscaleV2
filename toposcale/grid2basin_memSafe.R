#domain should be bigger than shape - actually not
require(raster)
require(ncdf4)
require(rgeos)
args = commandArgs(trailingOnly=TRUE)
wd=args[1]


writeNcdfSurf =function(varname,units,ncfname,title,dlname,data,time_nc,wd){
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

writeNcdfPlev =function(varname,units,ncfname,title,dlname,data,time_nc,wd,lev){
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
print ("Compute surface vars")
#read in meteo
nc = nc_open(surffile)
time_nc =ncvar_get(nc, 'time')
vars=c("t2m", 'z', 'd2m', 'strd', 'ssrd', 'tp', 'tisr')

for (var in vars){
	print(var)
t2m_nc =ncvar_get(nc, var)
t2m_mat=c()

#loop through basins
for ( b in basins){

# get grids that exist in basin
print(b)
cellIndex=na.omit(unlist(b))
#cellIndex=b

# get row/col of each grid
myrow = ceiling(cellIndex/ncol(s))
mycol = cellIndex-((myrow-1)*ncol(s))

# extract subset of variables corresponding to those grids
t2m_sub=t2m_nc[mycol,myrow,]

# compute grid ele - and lapse rate
# if (var=='tp'){
# z =ncvar_get(nc, 'z')
# z_sub=z[mycol,myrow,]
# x =apply(z_sub,MARGIN=c(1,2),FUN=mean)/9.81
# gridEle = diag(x)

# tp =ncvar_get(nc, 'tp')
# tp_sub=tp[mycol,myrow,]
# p =apply(tp_sub,MARGIN=c(1,2),FUN=sum)
# p2=diag(p)
# plot(p2/39,gridEle)
# }

if ( length(dim(t2m_sub))){
	# the mean timeseries for basin
	t2m_sub_mean = apply(t2m_sub, MARGIN=c(3), FUN='mean')

	# get some stats on the distribution
	#t2m_sub_sd = apply(t2m_sub, MARGIN=c(3), FUN='sd')
	#compute p lapse rate here
	
	}else{ # case of only one gridbox in basin
t2m_sub_mean <- t2m_sub
}
rm(t2m_sub)
gc()

t2m_mat=c(t2m_mat, t2m_sub_mean)

}

t2m = matrix(t2m_mat, nrow=length(basins), ncol=length(t2m_sub_mean), byrow=T)
rm(t2m_sub_mean)
gc()


# get centroids
centroids = gCentroid(shp,byid=T)
mylon=centroids@coords[,1]
mylat=centroids@coords[,2]

#=============================================
#	writeNCDF
#=============================================
	if(var=="z"){
	varname=var
	units="?"
	ncfname <- paste0(var,"s.nc")
	title <- "?"
	dlname <- "?"
	data=t2m
	writeNcdfSurf(varname,units,ncfname,title,dlname,data,time_nc,wd)
		}else{
	varname=var
	units="?"
	ncfname <- paste0(var,".nc")
	title <- "?"
	dlname <- "?"
	data=t2m
	writeNcdfSurf(varname,units,ncfname,title,dlname,data,time_nc,wd)
	#tp[tp<0]<-0
	#tp[which(is.infinite(as.vector(tp)))]<-0 # finds and remove Inf's why these exist?! must be TP correct method
	}
}
rm(t2m)
gc()
#=============================================
#	PRESSURE LEVELS VARIABLES
#=============================================
print("compute plevel data per basin")
#read in meteo
nc = nc_open(plevfile)
lev =ncvar_get(nc, 'level')
vars=c('t', 'z','u','v','r')

for ( var in vars){
	print(var)
	t_nc =ncvar_get(nc, var)
	t_mat=c()

		#loop through basins
		for ( b in basins){

		# get grids that exist in basin
		print(b)
		cellIndex=b

		# get row/col of each grid
		myrow = ceiling(cellIndex/ncol(s))
		mycol = cellIndex-((myrow-1)*ncol(s))

		# extract subset of variables corresponding to those grids
		t_sub=t_nc[mycol,myrow,,]
	
			if ( length(dim(t_sub))>2){
				# the mean timeseries for basin
				t_sub_mean = as.vector(t(apply(t_sub, MARGIN=c(3,4), FUN='mean')))

				# get some stats on the distribution
				#t2m_sub_sd = apply(t2m_sub, MARGIN=c(3), FUN='sd')
				#compute p lapse rate here
				
					}else{
				t_sub_mean <- as.vector(t_sub)

			}
			rm(t_sub)
			
		t_mat=c(t_mat, t_sub_mean)
		}

	t = matrix(t_mat, nrow=length(basins), ncol=length(t_sub_mean), byrow=T)
	rm(t_sub_mean)

	# reshape to 3D
	t_array = array(t, dim=c(length(basins), length(time_nc), length(lev)))
	rm(t)


	# get centroids
	centroids = gCentroid(shp,byid=T)
	mylon=centroids@coords[,1]
	mylat=centroids@coords[,2]

	#=============================================
	#	writeNCDF - plevel
	#=============================================
	if(var=="z"){
	varname=var
	units="?"
	ncfname <- paste0(var,"p.nc")
	title <- "?"
	dlname <- "?"
	data=t_array
	writeNcdfPlev(varname,units,ncfname,title,dlname,data,time_nc,wd,lev)
		}else{

	varname=var
	units="?"
	ncfname <- paste0(var,".nc")
	title <- "?"
	dlname <- "?"
	data=t_array
	writeNcdfPlev(varname,units,ncfname,title,dlname,data,time_nc,wd,lev)
}




