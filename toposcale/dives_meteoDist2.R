#DEPENDENCY
require(raster)
require(dplyr)
require(reshape2)
require(ncdf4)
#SOURCE
#source("./rsrc/toposub_src.R")
#====================================================================
# PARAMETERS/ARGS
#====================================================================
#args <- 	commandArgs(trailingOnly=TRUE)
#wd <- args[1]
wd="/home/joel/sim/adrian_GR_basin/sim/g1/"
wdout="/home/joel/tscale3dout/"
startDate="1979-08-31" # this is one day before first timstamp as ncdf index starts at 1 wjich is 1day since startDate

year_list_prescribed=TRUE
year_list_given=c("2004")

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


# write netcdf function

writeNcdf=function(

	lf2,
	mydatesIndex,
	mydata,
	startDate,
	par,
	nsamp,
	varname,
	longname,
	units,
	filename,
	title
	){



		for (i in 1:nsamp){
			print(i)
			my1vec=rep(mydata[,i], each =length(which(lf==i)))
			lf2[which(lf2==i)]<-my1vec
		}


	out = (array(data = lf2, dim = c(ncol(landform), nrow(landform), ndays_cut) ))
	mydata2 <- aperm(out, c(2,1,3))
	rm(out)
	gc()

	r.pts <- rasterToPoints(landform, spatial=TRUE)
	lon = sort(unique(r.pts@coords[,1]))
	lat = rev(sort(unique(r.pts@coords[,2])))


	# create a netCDF file 
	# define dimensions
	londim <- ncdim_def("lon", "degrees_east", as.double(lon))
	latdim <- ncdim_def("lat", "degrees_north", as.double(lat))
	tdim <- ncdim_def("time", paste0("days since ",startDate), as.double(mydatesIndex))
	# define variables
	varname=varname
	units=units
	#dlname <- "test variable -- original"
	fillvalue <- -9999
	tmp.def <- ncvar_def(name =varname,  units =units, dim =list(latdim, londim,tdim), missval = fillvalue, longname=longname, prec = "single")


	ncfname <- filename
	ncout <- nc_create(ncfname, list(tmp.def), force_v4 = TRUE)

	# put the array
	ncvar_put(ncout, tmp.def, mydata2)

	# put additional attributes into dimension and data variables
	ncatt_put(ncout, "lon", "axis", "X")  
	ncatt_put(ncout, "lat", "axis", "Y")

	# add global attributes
	title <-title
	ncatt_put(ncout, 0, "title", title)

	# close the file, writing data to disk
	nc_close(ncout)


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

rm (sRes.list)
gc()
sRes = sRes1[,,correctOrder]
#save(sRes,file = paste0(wd, "/meteocube"))
rm(sRes1)
gc()
        



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
lf =getValues(landform)
years=substr(day_vec,1,4)
year_list=unique(years)

if(year_list_prescribed==TRUE){
	year_list <- year_list_given
	}
#=============== DW MEAN ========================================================
param='DW'
longname="Wind direction daily mean"
par=which(params==param) 
units="deg"
filename=paste0(param,"_", yr,".nc")

	sRes_day=array(0,c(ndays,nsamples)) # init empty array of daily values
	for (sample in 1:nsamples){
	print(sample)
	df=sRes[,par,sample] # subset from main array
	out=aggregate(df, list(day_agg), mean) # aggregate
	sRes_day[,sample]<- out$x
	}


# build for single years
for (yr in year_list){
	filename=paste0(wdout,"/",param,"_", yr,".nc")
	mydatesIndex=which(years==yr)
	ndays_cut=length(mydatesIndex)
	lf2 =rep(lf, ndays_cut)
	mydata =sRes_day[mydatesIndex,]

	# now generate netcdf of results
	writeNcdf(lf2,mydatesIndex, mydata, startDate, par,nsamp=nsamples,varname=param,longname,units,filename,title=paste0("Tscale ",param))
	print(paste0("Written: ", wd,param,yr,".tif"))
	}


#=============== TA MEAN =======================================================
param='TA'
longname='Air temperature daily mean'
par=which(params==param) 
units='K'


sRes_day=array(0,c(ndays,nsamples)) # init empty array of daily values
	for (sample in 1:nsamples){
	print(sample)
	df=sRes[,par,sample] # subset from main array
	out=aggregate(df, list(day_agg), min) # aggregate
	sRes_day[,sample]<- out$x
	}

# build for single years
for (yr in year_list){
	filename=paste0(wdout,"/",param,"MEAN_", yr,".nc")
	mydatesIndex=which(years==yr)
	ndays_cut=length(mydatesIndex)
	lf2 =rep(lf, ndays_cut)
	mydata =sRes_day[mydatesIndex,]

	# now generate netcdf of results
	writeNcdf(lf2,mydatesIndex, mydata, startDate, par,nsamp=nsamples,varname=param,longname,units,filename,title=paste0("Tscale ",param))
	print(paste0("Written: ", wd,param,yr,".tif"))
	}


#=============== TA MAX =======================================================
param='TA'
longname='Air temperature daily max'
par=which(params==param) 
units='K'


sRes_day=array(0,c(ndays,nsamples)) # init empty array of daily values
	for (sample in 1:nsamples){
	print(sample)
	df=sRes[,par,sample] # subset from main array
	out=aggregate(df, list(day_agg), max) # aggregate
	sRes_day[,sample]<- out$x
	}

# build for single years
for (yr in year_list){
	filename=paste0(wdout,"/",param,"MAX_", yr,".nc")
	mydatesIndex=which(years==yr)
	ndays_cut=length(mydatesIndex)
	lf2 =rep(lf, ndays_cut)
	mydata =sRes_day[mydatesIndex,]

	# now generate netcdf of results
	writeNcdf(lf2,mydatesIndex, mydata, startDate, par,nsamp=nsamples,varname=param,longname,units,filename,title=paste0("Tscale ",param))
	print(paste0("Written: ", wd,param,yr,".tif"))
	}


#=============== ISWR =======================================================
param='ISWR'
longname='ISWR daily sum'
watts2megajoules = 0.0036
par=which(params==param) 
units="MJm**2"


sRes_day=array(0,c(ndays,nsamples)) # init empty array of daily values
	for (sample in 1:nsamples){
	print(sample)
	df=sRes[,par,sample] # subset from main array
	out=aggregate(df, list(day_agg), sum) # aggregate
	sRes_day[,sample]<- out$x
	}
sRes_day=sRes_day*watts2megajoules

# build for single years
for (yr in year_list){
	filename=paste0(wdout,"/",param,"SUM_", yr,".nc")
	mydatesIndex=which(years==yr)
	ndays_cut=length(mydatesIndex)
	lf2 =rep(lf, ndays_cut)
	mydata =sRes_day[mydatesIndex,]

	# now generate netcdf of results
	writeNcdf(lf2,mydatesIndex, mydata, startDate, par,nsamp=nsamples,varname=param,longname,units,filename,title=paste0("Tscale ",param))
	print(paste0("Written: ", wd,param,yr,".tif"))
	}

#=============== VW MEAN=======================================================
param="VW"
longname='Wind velocity daily mean'

par=which(params==param) 
units='ms**1'
filename=paste0(param,"MEAN_", yr,".nc")

sRes_day=array(0,c(ndays,nsamples)) # init empty array of daily values
	for (sample in 1:nsamples){
	print(sample)
	df=sRes[,par,sample] # subset from main array
	out=aggregate(df, list(day_agg), mean) # aggregate
	sRes_day[,sample]<- out$x
	}
# build for single years
for (yr in year_list){
	filename=paste0(wdout,"/",param,"MEAN_", yr,".nc")
	mydatesIndex=which(years==yr)
	ndays_cut=length(mydatesIndex)
	lf2 =rep(lf, ndays_cut)
	mydata =sRes_day[mydatesIndex,]

	# now generate netcdf of results
	writeNcdf(lf2,mydatesIndex, mydata, startDate, par,nsamp=nsamples,varname=param,longname,units,filename,title=paste0("Tscale ",param))
	print(paste0("Written: ", wd,param,yr,".tif"))
	}

#=============== VW MAX=======================================================
param='VW'
longname='Wind velocity daily max'
par=which(params==param) 
units='ms**1'


sRes_day=array(0,c(ndays,nsamples)) # init empty array of daily values
	for (sample in 1:nsamples){
	print(sample)
	df=sRes[,par,sample] # subset from main array
	out=aggregate(df, list(day_agg), max) # aggregate
	sRes_day[,sample]<- out$x
	}
# build for single years
for (yr in year_list){
	filename=paste0(wdout,"/",param,"MAX_", yr,".nc")
	mydatesIndex=which(years==yr)
	ndays_cut=length(mydatesIndex)
	lf2 =rep(lf, ndays_cut)
	mydata =sRes_day[mydatesIndex,]

	# now generate netcdf of results
	writeNcdf(lf2,mydatesIndex, mydata, startDate, par,nsamp=nsamples,varname=param,longname,units,filename,title=paste0("Tscale ",param))
	print(paste0("Written: ", wd,param,yr,".tif"))
	}

#=============== PSUM =======================================================
param='PSUM'
longname='Precipitation daily sum'
par=which(params==param) 
units='mm'



sRes_day=array(0,c(ndays,nsamples)) # init empty array of daily values
	for (sample in 1:nsamples){
	print(sample)
	df=sRes[,par,sample] # subset from main array
	out=aggregate(df, list(day_agg), sum) # aggregate
	sRes_day[,sample]<- out$x
	}

# build for single years
for (yr in year_list){
	filename=paste0(wdout,"/",param,"_", yr,".nc")
	mydatesIndex=which(years==yr)
	ndays_cut=length(mydatesIndex)
	lf2 =rep(lf, ndays_cut)
	mydata =sRes_day[mydatesIndex,]

	# now generate netcdf of results
	writeNcdf(lf2,mydatesIndex, mydata, startDate, par,nsamp=nsamples,varname=param,longname,units,filename,title=paste0("Tscale ",param))
	print(paste0("Written: ", wd,param,yr,".tif"))
	}



#=============== relative humidity =======================================
param='RH'
longname='Relative humidity daily mean'
par=which(params==param) 
units='%'
sRes_day=array(0,c(ndays,nsamples)) # init empty array of daily values
	for (sample in 1:nsamples){
	print(sample)
	df=sRes[,par,sample] # subset from main array
	out=aggregate(df, list(day_agg), mean) # aggregate
	sRes_day[,sample]<- out$x
	}


# build for single years
for (yr in year_list){
	filename=paste0(wdout,"/",param,"MEAN_", yr,".nc")
	mydatesIndex=which(years==yr)
	ndays_cut=length(mydatesIndex)
	lf2 =rep(lf, ndays_cut)
	mydata =sRes_day[mydatesIndex,]

	# now generate netcdf of results
	writeNcdf(lf2,mydatesIndex, mydata, startDate, par,nsamp=nsamples,varname=param,longname,units,filename,title=paste0("Tscale ",param))
	print(paste0("Written: ", wd,param,yr,".tif"))
	}


#=============== vapour pressure deficit =======================================
param='VPD'
longname='Vapour pressure deficit daily mean'
par=which(params==param) 
units='kPa'

for (yr in year_list){

	rh=nc_open(filename=paste0(wdout,"/RH_", yr,".nc"))
	t=nc_open(filename=paste0(wdout,"/TAMEAN_", yr,".nc"))
	rhdat = ncvar_get(rh, "RH")
	tdat = ncvar_get(t, "TA")
	# calc saturated vapour pressure (pascals)
	# T in deg 
	T=tdat-273.13
	svp = 610.7*10^(7.5*T/(237.3+T))
	# As VPD is the saturated vapour pressure minus the actual vapour pressure (SVP - VPactual), and VPactual = (RH*SVP)/100
	# RH in 0-100
	RH=rhdat*100
	VPD = ((100 - RH)/100)*svp

	vpd_kpa = VPD/1000 # pa to kPa

	# now generate netcdf of results vpd_kpa is already in geo format so no need to use main function
	filename=paste0(wdout,"/",param,"MEAN_", yr,".nc")
	mydatesIndex=which(years==yr)
	mydata2 =vpd_kpa
	nsamp=nsamples
	varname=param
	title=paste0("Tscale ",param)

	r.pts <- rasterToPoints(landform, spatial=TRUE)
	lon = sort(unique(r.pts@coords[,1]))
	lat = rev(sort(unique(r.pts@coords[,2])))


	# create a netCDF file 
	# define dimensions
	londim <- ncdim_def("lon", "degrees_east", as.double(lon))
	latdim <- ncdim_def("lat", "degrees_north", as.double(lat))
	tdim <- ncdim_def("time", paste0("days since ",startDate), as.double(mydatesIndex))
	# define variables
	varname=varname
	units=units
	#dlname <- "test variable -- original"
	fillvalue <- -9999
	tmp.def <- ncvar_def(name =varname,  units =units, dim =list(latdim, londim,tdim), missval = fillvalue, longname=longname, prec = "single")


	ncfname <- filename
	ncout <- nc_create(ncfname, list(tmp.def), force_v4 = TRUE)

	# put the array
	ncvar_put(ncout, tmp.def, mydata2)

	# put additional attributes into dimension and data variables
	ncatt_put(ncout, "lon", "axis", "X")  
	ncatt_put(ncout, "lat", "axis", "Y")

	# add global attributes
	title <-title
	ncatt_put(ncout, 0, "title", title)

	# close the file, writing data to disk
	nc_close(ncout)


	print(paste0("Written: ", wdout,param,yr,".tif"))
	}





par(mfrow=c(3,3))
files = list.files(wdout,pattern="*2004.nc", full.names=TRUE)
require(raster)
require(viridis)
for (f in files[2:9]){
r = mean(stack(f))
plot(r, main=f, col=viridis(100))
}

