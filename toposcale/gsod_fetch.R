#gsod_fetch.R
# request gsod data
require(GSODR) #https://docs.ropensci.org/GSODR/reference/get_GSOD.html
YEARS=1979:2018
country="tj"
eleFilter=2000
out="~/"

# latest inventory
inventory <- get_inventory()
gsod <- get_GSOD(years = 2014, country=country) # only works for single years
stat=unique(gsod$NAME[which(gsod$ELEVATION>eleFilter)])
STNID = unique(gsod$STNID[which(gsod$ELEVATION>eleFilter)])
# a nearest station request
# nearest_stations(39,70,200)

snow_list = list()
date_list = list()
plotdim = ceiling(sqrt(length(stat)))
par(mfrow=c(plotdim, plotdim))

for ( i in 1:length(STNID) ){

print(STNID[i])
# get stations of TJ

gsod <- get_GSOD(years = YEARS, station=STNID[i])
#GSOD_SF <- st_as_sf(x = gsod, coords = c("LONGITUDE", "LATITUDE"),  crs = "+proj=longlat +datum=WGS84")
# stations filter here
# stat = unique(gsod$NAME)
#plotdim = ceiling(sqrt(length(stat)))
#par(mfrow=c(plotdim, plotdim))
#for (i in (1:length(stat))){
	print(i)
	ele=unique(gsod$ELEVATION)#[gsod$NAME==stat[i]])
statname=unique(gsod$NAME)
	snow=gsod$SNDP#[gsod$NAME==stat[i]]
	#if(length(is.na(T)[is.na(T)==TRUE])==length(T)){print("skip");next}
	date=gsod$YEARMODA #[gsod$NAME==stat[i]]
	myyears=unique(substr(date,1,4))
# date=gsod$YEARMODA
# snow = gsod$SNDP
# plot(date, snow)

if( sum(snow,na.rm=T)==0  ) {
plot(date,rep(0,length(date)), main=paste(statname, '/',ele,"/"), type='l')
}else{
	plot(date,snow, main=paste(statname, '/',ele,"/"), type='l')
	abline(h=0)
}
snow_list[[i]] <- snow
date_list[[i]] <- date

}
#}












# make dataframe
df=data.frame(name=character(),ele=integer(),lon=double(),lat=double())
for (i in (1:length(stat))){

	print(i)
	lon =unique(gsod$LONGITUDE[gsod$NAME==stat[i]])
	lat = unique(gsod$LATITUDE[gsod$NAME==stat[i]])
	ele=unique(gsod$ELEVATION[gsod$NAME==stat[i]])
	if(length(lon)==0){lon=-999 ; print("Warning! LAT missing!");next}
	if(length(lat)==0){lat=-999; print("Warning! LON missing!");next}
	if(length(ele)==0){ele=-999; print("Warning! ELE missing!")}
	name = stat[i]

	mypos =cbind(name,ele,lon,lat)
	# df$name[i]=stat[i]
	# df$ele[i]=ele
	# df$lon[i]=lon
	# df$lat[i]=lat
	df=rbind(df,mypos)
}
df$pk=1:length(df$name)
outname=paste0(out,"/",country,year,"over", eleFilter)




write.csv(df, paste0(outname,".csv"),row.names=FALSE)
# make shape





# PARAMETERS
longcol=3 #x
latcol=4 # y
outfile=paste0(outname,'.shp')
infile=paste0(outname,'.csv')
skip=0
header=TRUE
sep=","

#EPSG:4326
#WGS 1984 (Google it)
proj="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"


# FUNCTION
makePointShapeGeneriRc=function(lon,lat,data,proj=proj){
library(raster)
library(rgdal)
loc<-data.frame(data[,lon], data[,lat])
spoints<-SpatialPointsDataFrame(loc,as.data.frame(data), proj4string= CRS(proj))
return(spoints)
}

#CODE
df=read.csv(infile, skip=skip, sep=sep, header=header)

# replace values

df$ele[df$name=='DEHAVZ']<-2571
df$lat[df$name=='DEHAVZ']<-39.44810
df$lon[df$name=='DEHAVZ']<-70.194489

df$ele[df$name=='MADRUSHKAT']<-2240
df$lat[df$name=='MADRUSHKAT']<-39.44654
df$lon[df$name=='MADRUSHKAT']<-69.66453

df$ele[df$name=='ISKANDERKUL']<-2206
df$lat[df$name=='ISKANDERKUL']<-39.08476
df$lon[df$name=='ISKANDERKUL']<-68.37087

df$ele[df$name=='ANZQBSKIJ PEREVAL']<-3369
df$lat[df$name=='ANZQBSKIJ PEREVAL']<-39.08317
df$lon[df$name=='ANZQBSKIJ PEREVAL']<-68.86502

df$ele[df$name=='KHOROG']<-2084
df$lat[df$name=='KHOROG']<-37.50138
df$lon[df$name=='KHOROG']<-68.86502

df$ele[df$name=='ISHKASHIM']<-2533
df$lat[df$name=='ISHKASHIM']<-36.72433
df$lon[df$name=='ISHKASHIM']<-71.614622

df$ele[df$name=='KARAKUL']<-3946
df$lat[df$name=='KARAKUL']<-39.00915
df$lon[df$name=='KARAKUL']<-73.56011

df$ele[df$name=='IHRT']<-3278
df$lat[df$name=='IHRT']<-38.16537
df$lon[df$name=='IHRT']<-72.61638

# write out again
write.csv(df, paste0(outname,".csv"),row.names=FALSE)

shp=makePointShapeGeneriRc(lon=3,lat=4,data=df,proj=proj)


#WRITE
shapefile(x=shp,filename=outfile,overwrite=TRUE)


