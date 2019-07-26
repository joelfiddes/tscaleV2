#DEPENDENCY
require(raster)
require(dplyr)
require(reshape2)
#SOURCE
#source("./rsrc/toposub_src.R")
#====================================================================
# PARAMETERS/ARGS
#====================================================================
args <- 	commandArgs(trailingOnly=TRUE)
wd <- args[1]

writeMeans=FALSE # switch to justr write means of all variables
#wd="/home/joel/sim/adrian_GR_basin/sim/g1/"
landform<-raster(paste0(wd, "/landform.tif")	)

crispSpatialNow = function(resultsVec, landform){
		require(raster)
		l <- length(resultsVec)
		s <- 1:l
		df <- data.frame(s,resultsVec)
		rst <- subs(landform, df,by=1, which=2)
		rst=round(rst,2)
		return(rst)
		}


#make a results cube
### Function to convert a list of dataframes to a 3D array
## All objects in the list will be dataframes with identical column headings.
list2ary = function(input.list){  #input a list of lists
  rows.cols <- dim(input.list[[1]])
  sheets <- length(input.list)
  output.ary <- array(unlist(input.list), dim = c(rows.cols, sheets))
  colnames(output.ary) <- colnames(input.list[[1]])
  row.names(output.ary) <- row.names(input.list[[1]])
  return(output.ary)    # output as a 3-D array
}

# not in correct order!!
sRes.names = list.files(wd, pattern = '^meteo.*csv$', recursive=TRUE, full.name=TRUE) # finds only main csv files

# a load of str parsing to get correct order eg 1-150
mylist = strsplit(sRes.names, '/')
mysplit = unique(rapply(mylist, function(x) tail(x, 1)))
mysub = strsplit(mysplit, "meteoc")
mysplit2 = unique(rapply(mysub, function(x) tail(x, 1)))
mysub2 = strsplit(mysplit2, ".csv")
correctOrder= order(formatC(unlist(mysub2), width=3, flag='0'))


sRes.list = lapply(sRes.names, FUN=read.csv,  sep=',', header=T)
sRes1 <- list2ary(sRes.list)  # convert to array
sRes = sRes1[,,correctOrder]
save(sRes,file = paste0(wd, "/meteocube"))




if (writeMeans==TRUE){
# loop thru all variables (dimension2)
params = names(sRes[1,2:9,1])

# aggregation factors
datetime=read.csv(sRes.names[1])
day_agg = substr(datetime$datetime,1,10)
#generate forcing grids
npars=dim(sRes)[2]
nsamples=dim(sRes)[3]
ndays=length(unique(day_agg))


# init empty array of daily values
sRes_day=array(0,c(ndays, npars,nsamples))
for (par in 1:npars){
for (sample in 1:nsamples){
print(par)
print(sample)
df=sRes[,par,sample]
 out=aggregate(df, list(day_agg), mean)
sRes_day[,par,sample]<- out$x
}}


# now generate grids for every day and parameter (starting at 2)
for (par in 2:npars){

# init rstack with day1
rst =crispSpatialNow(sRes_day[1,par,],landform)
rstack=stack(rst)

for (day in 2:ndays){
print(par)
print(day)
rst =crispSpatialNow(sRes_day[day,par,],landform)

rstack[[day]]<-rst
}

# with dynamic name relating to var
assign(paste0('rstack',par) , rstack)

}


writeRaster(rstack2, paste0(wd,params[1],".tif"),overwrite=TRUE)
writeRaster(rstack3, paste0(wd,params[2],".tif"),overwrite=TRUE)
writeRaster(rstack4, paste0(wd,params[3],".tif"),overwrite=TRUE)
writeRaster(rstack5, paste0(wd,params[4],".tif"),overwrite=TRUE)
writeRaster(rstack6, paste0(wd,params[5],".tif"),overwrite=TRUE)
writeRaster(rstack7, paste0(wd,params[6],".tif"),overwrite=TRUE)
writeRaster(rstack8, paste0(wd,params[7],".tif"),overwrite=TRUE)
writeRaster(rstack9, paste0(wd,params[8],".tif"),overwrite=TRUE)

}

#===============================================================================
#		DIVES CASE BY CASE
#===============================================================================
params = names(sRes[1,,1])
# aggregation factors
datetime=read.csv(sRes.names[1])
day_agg = substr(datetime$datetime,1,10)
npars=dim(sRes)[2]
nsamples=dim(sRes)[3]
ndays=length(unique(day_agg))
day_vec=unique(day_agg)


#=============== DW MEAN ========================================================
param='DW'
par=which(params==param) 

sRes_day=array(0,c(ndays,nsamples)) # init empty array of daily values
	for (sample in 1:nsamples){
	print(sample)
	df=sRes[,par,sample] # subset from main array
	out=aggregate(df, list(day_agg), mean) # aggregate
	sRes_day[,sample]<- out$x
	}
# now generate grids for every day 
# init rstack with day1
rst =crispSpatialNow(sRes_day[1,],landform)
rstack=stack(rst)
	for (day in 2:ndays){
	print(day)
	rst =crispSpatialNow(sRes_day[day,],landform)
	rstack[[day]]<-rst
	}
assign(paste0('rstack',par) , rstack)
writeRaster(rstack, paste0(wd,param,".tif"),overwrite=TRUE)
print(paste0("Written: ", wd,param,".tif"))

#=============== TA MIN =======================================================
param='TA'
par=which(params==param) 

sRes_day=array(0,c(ndays,nsamples)) # init empty array of daily values
	for (sample in 1:nsamples){
	print(sample)
	df=sRes[,par,sample] # subset from main array
	out=aggregate(df, list(day_agg), min) # aggregate
	sRes_day[,sample]<- out$x
	}
# now generate grids for every day 
# init rstack with day1
rst =crispSpatialNow(sRes_day[1,],landform)
rstack=stack(rst)
	for (day in 2:ndays){
	print(day)
	rst =crispSpatialNow(sRes_day[day,],landform)
	rstack[[day]]<-rst
	}
assign(paste0('rstack',par) , rstack)
writeRaster(rstack, paste0(wd,param,"_min.tif"),overwrite=TRUE)
print(paste0("Written: ", wd,param,"_min.tif"))

#=============== TA MAX =======================================================
param='TA'
par=which(params==param) 

sRes_day=array(0,c(ndays,nsamples)) # init empty array of daily values
	for (sample in 1:nsamples){
	print(sample)
	df=sRes[,par,sample] # subset from main array
	out=aggregate(df, list(day_agg), max) # aggregate
	sRes_day[,sample]<- out$x
	}
# now generate grids for every day 
# init rstack with day1
rst =crispSpatialNow(sRes_day[1,],landform)
rstack=stack(rst)
	for (day in 2:ndays){
	print(day)
	rst =crispSpatialNow(sRes_day[day,],landform)
	rstack[[day]]<-rst
	}
assign(paste0('rstack',par) , rstack)
writeRaster(rstack, paste0(wd,param,"_max.tif"),overwrite=TRUE)
print(paste0("Written: ", wd,param,"_max.tif"))

#=============== ISWR =======================================================
param='ISWR'
watts2joules = 0.0036

par=which(params==param) 

sRes_day=array(0,c(ndays,nsamples)) # init empty array of daily values
	for (sample in 1:nsamples){
	print(sample)
	df=sRes[,par,sample] # subset from main array
	out=aggregate(df, list(day_agg), sum) # aggregate
	sRes_day[,sample]<- out$x
	}
# now generate grids for every day 
# init rstack with day1
rst =crispSpatialNow(sRes_day[1,],landform)
rstack=stack(rst)
	for (day in 2:ndays){
	print(day)
	rst =crispSpatialNow(sRes_day[day,],landform)
	rstack[[day]]<-rst
	}
assign(paste0('rstack',par) , rstack)
writeRaster(rstack*watts2joules, paste0(wd,param,"_sum.tif"),overwrite=TRUE)
print(paste0("Written: ", wd,param,"_sum.tif"))

#=============== VW MEAN=======================================================
param='VW'
par=which(params==param) 

sRes_day=array(0,c(ndays,nsamples)) # init empty array of daily values
	for (sample in 1:nsamples){
	print(sample)
	df=sRes[,par,sample] # subset from main array
	out=aggregate(df, list(day_agg), mean) # aggregate
	sRes_day[,sample]<- out$x
	}
# now generate grids for every day 
# init rstack with day1
rst =crispSpatialNow(sRes_day[1,],landform)
rstack=stack(rst)
	for (day in 2:ndays){
	print(day)
	rst =crispSpatialNow(sRes_day[day,],landform)
	rstack[[day]]<-rst
	}
assign(paste0('rstack',par) , rstack)
writeRaster(rstack, paste0(wd,param,"_mean.tif"),overwrite=TRUE)
print(paste0("Written: ", wd,param,"_mean.tif"))

#=============== VW MAX=======================================================
param='VW'
par=which(params==param) 

sRes_day=array(0,c(ndays,nsamples)) # init empty array of daily values
	for (sample in 1:nsamples){
	print(sample)
	df=sRes[,par,sample] # subset from main array
	out=aggregate(df, list(day_agg), max) # aggregate
	sRes_day[,sample]<- out$x
	}
# now generate grids for every day 
# init rstack with day1
rst =crispSpatialNow(sRes_day[1,],landform)
rstack=stack(rst)
	for (day in 2:ndays){
	print(day)
	rst =crispSpatialNow(sRes_day[day,],landform)
	rstack[[day]]<-rst
	}
assign(paste0('rstack',par) , rstack)
writeRaster(rstack, paste0(wd,param,"_max.tif"),overwrite=TRUE)
print(paste0("Written: ", wd,param,"_max.tif"))

#=============== PSUM =======================================================
param='PSUM'
par=which(params==param) 

sRes_day=array(0,c(ndays,nsamples)) # init empty array of daily values
	for (sample in 1:nsamples){
	print(sample)
	df=sRes[,par,sample] # subset from main array
	out=aggregate(df, list(day_agg), sum) # aggregate
	sRes_day[,sample]<- out$x
	}
# now generate grids for every day 
# init rstack with day1
rst =crispSpatialNow(sRes_day[1,],landform)
rstack=stack(rst)
	for (day in 2:ndays){
	print(day)
	rst =crispSpatialNow(sRes_day[day,],landform)
	rstack[[day]]<-rst
	}
assign(paste0('rstack',par) , rstack)
writeRaster(rstack, paste0(wd,param,".tif"),overwrite=TRUE)
print(paste0("Written: ", wd,param,".tif"))

#=============== TA MESAN =======================================================
param='TA'
par=which(params==param) 

sRes_day=array(0,c(ndays,nsamples)) # init empty array of daily values
	for (sample in 1:nsamples){
	print(sample)
	df=sRes[,par,sample] # subset from main array
	out=aggregate(df, list(day_agg), max) # aggregate
	sRes_day[,sample]<- out$x
	}
# now generate grids for every day 
# init rstack with day1
rst =crispSpatialNow(sRes_day[1,],landform)
rstackT=stack(rst)
	for (day in 2:ndays){
	print(day)
	rst =crispSpatialNow(sRes_day[day,],landform)
	rstackT[[day]]<-rst
	}
assign(paste0('rstack',par) , rstack)
writeRaster(rstackT, paste0(wd,param,"_mean.tif"),overwrite=TRUE)
print(paste0("Written: ", wd,param,"_mean.tif"))
#=============== vapour pressure deficit =======================================
param='RH'
par=which(params==param) 

sRes_day=array(0,c(ndays,nsamples)) # init empty array of daily values
	for (sample in 1:nsamples){
	print(sample)
	df=sRes[,par,sample] # subset from main array
	out=aggregate(df, list(day_agg), mean) # aggregate
	sRes_day[,sample]<- out$x
	}
# now generate grids for every day 
# init rstack with day1
rst =crispSpatialNow(sRes_day[1,],landform)
rstackRH=stack(rst)
	for (day in 2:ndays){
	print(day)
	rst =crispSpatialNow(sRes_day[day,],landform)
	rstackRH[[day]]<-rst
	}
assign(paste0('rstack',par) , rstack)
writeRaster(rstackRH, paste0(wd,param,".tif"),overwrite=TRUE)
print(paste0("Written: ", wd,param,".tif"))

# calc saturated vapour pressure (pascals)
# T in deg 
T=rstackT-273.13
svp = 610.7*10^(7.5*T/(237.3+T))

# As VPD is the saturated vapour pressure minus the actual vapour pressure (SVP - VPactual), and VPactual = (RH*SVP)/100
# RH in 0-100
RH=rstackRH*100
VPD = ((100 - RH)/100)*svp

vpd_kpa = VPD/1000 # pa to kpa
writeRaster(VPD,"rstackRH.tif",overwrite=TRUE)
print(paste0("Written: rstackRH.tif"))





dat=sRes_day[,8,]
lf =getValues(landform)
lf2 =rep(lf, ndays)
	for (i in 1:150){
		print(i)
		my1vec=rep(dat[,i], each =length(which(lf==i)))
		lf2[which(lf2==i)]<-my1vec
	}


out = (array(data = lf2, dim = c(ncol(landform), nrow(landform), ndays) ))

new.data <- aperm(out, c(2,1,3))
day=1
par(mfrow=c(1,2))
plot(raster(new.data[,,day]))

plot(crispSpatialNow(sRes_day[day,par,],landform))


r.pts <- rasterToPoints(landform, spatial=TRUE)
lon = sort(unique(r.pts@coords[,1]))
lat = rev(unique(r.pts@coords[,2])))


# create a netCDF file 
# define dimensions
londim <- ncdim_def("lon", "degrees_east", as.double(lon))
latdim <- ncdim_def("lat", "degrees_north", as.double(lat))
tdim <- ncdim_def("time", paste0("days since ",day_vec[1]), as.double(0:(length(day_vec)-1)))
# define variables
varname="T"
units="K"
#dlname <- "test variable -- original"
fillvalue <- -9999
tmp.def <- ncvar_def(varname, units, list(latdim, londim,tdim), fillvalue, prec = "single")

ncfname <- "TA.nc"
ncout <- nc_create(ncfname, list(tmp.def), force_v4 = TRUE)

# put the array
ncvar_put(ncout, tmp.def, new.data)

# put additional attributes into dimension and data variables
ncatt_put(ncout, "lon", "axis", "X")  
ncatt_put(ncout, "lat", "axis", "Y")

# add global attributes
title <- "TSCALE TA"
ncatt_put(ncout, 0, "title", title)

# close the file, writing data to disk
nc_close(ncout)
raster("TA.nc")
