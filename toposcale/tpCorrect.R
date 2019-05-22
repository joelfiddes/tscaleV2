# takes tpmm.nc and SURF.nc (defined in fetch_era5.py)
# computes difference of monthly totals and apply this as weights to each hourly value of each month
# writes the updated totals back to file
# if run twice by accident then the correction factor is zero and therefor no correction is done (dummy proof)

# the idea is to download 6h tp then correct by monthly means to account for full monthly budget
# monthly means are super fast to download (smaller data and faster server "spinning disks!")
# same idea for cmip5 (monthly cc signal)
# all data? 1 year of hourly distribution, scale by monthly means - run this past seb!!

# documented era5 hgere: https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation#ERA5datadocumentation-Monthlymeans
#values are mean daily for that month, weight by number of days in month eg annualtotal = sum(wfj_m*c(31,28,31,30,31,30,31,31,30,31,30,31))

# THIS IS THE
require(ncdf4)
require(raster)

# pars 
args = commandArgs(trailingOnly=TRUE)
mmfile=args[1] # monthly mean tp file eg "/home/joel/sim/imis/forcing/tpmm.nc"
t6hfile = args[2] # subdaily file likely 6 or 3 h step but can be anything eg 

#mmfile="/home/joel/sim/imis/forcing/tpmm.nc"
#t6hfile ="/home/joel/sim/imis/forcing/SURF.nc"



# make backup of file as replace values in place
system(paste0("cp " , t6hfile," ", t6hfile, "_backup.nc"))


# get monthly dates
nc=nc_open(mmfile)
time = ncvar_get( nc,'time')
z <- time*60*60 #make seconds
origin = unlist(strsplit(nc$dim$time$units, " "))[[3]]
datesMonth<-ISOdatetime(origin,0,0,0,0,0,tz='UTC') + z #dates sequence
tpmm = ncvar_get( nc,'tp')

# number of years in monthly data
nyears=dim(tpmm)[3]/12

# number of days per month vec, ignore leap years pain in ass to implement and insignificant error
monthLengthVec = rep(c(31,28,31,30,31,30,31,31,30,31,30,31),nyears) 

# weighted by days in the month
tpmm2 =tpmm*monthLengthVec




nc1=nc_open(t6hfile, write=T)
time = ncvar_get( nc1,'time')
z <- time*60*60 #make seconds
origin = unlist(strsplit(nc1$dim$time$units, " "))[[3]]
datesHr<-ISOdatetime(origin,0,0,0,0,0,tz='UTC') + z #dates sequence
tphr = ncvar_get( nc1,'tp')
# make monthly totals of subdaily data
mo =substr(datesHr, 1,7)

# calculate timestep
timestep = as.numeric(datesHr[2] -datesHr[1])
dailyCount = 24 / timestep # counts per day

# cut to length of obs by computing startMonth and endMonth of obs relative to 
# full years given by monthly means 
umo = unique(mo)
momm =substr(datesMonth, 1,7)
startMonth = which(momm==umo[1])
endMonth = which(momm==umo[length(umo)])
tpmm3 = tpmm2[,,startMonth:endMonth]


# hourly data aggregated as monthly mean
a =apply(tphr,c(1,2), 'aggregate', list(mo), sum )
plot(tpmm3[1,1,], type='l')
lines(a[[1]]$x, type='l', col='red')

# check sums
sum(tpmm3)
sum(tphr*timestep)
sum(as.numeric(as.vector(unlist(a))),na.rm=T)*timestep

# flatten
tpmm4 = matrix(tpmm3,  dim(tpmm3)[1]*dim(tpmm3)[2], dim(tpmm3)[3])


# generates list of monthly correction factors
cfl=c()
for ( i in 1:length(a) ){

	hrmm = a[[i]]$x*timestep
	mm = tpmm4[i,]

	cfl[[i]]<-(mm/hrmm)

	#plot(hrmm, type='l', ylim=c(0,0.5))
	#lines(mm,col='red')

	}

# flatten six hour data matrix
tph_flat = matrix(tphr, dim(tphr)[1]*dim(tphr)[2], dim(tphr)[3])

corP=c()
for ( i in 1:length(a) ){

	# correction vector same length as hourly data
	corVec = rep(unlist(cfl[i]),times = monthLengthVec[startMonth:endMonth]*dailyCount)

	# in case of unequal vecytors this means one or more leap years is present in dataset
	# implementing these would be costly so just recyle last element now - of course this can lead to shift of up to 10 days in case of full era record
	# needs attention!!

	if (length(corVec)!=length(tph_flat[i,]) ){

		ntimes=length(tph_flat[i,])- length(corVec)

		# paste replicates of last correction factor to extend timeseries
		corVec = c(corVec, rep(corVec[length(corVec)], ntimes))
		}
	# timeseries for grid i of corrected 6 h data adjusted by monthly totals and scaled by timestep
	corP = c(corP,tph_flat[i,]*corVec*timestep )
}


corP_array = matrix(corP , dim(tphr)[1] *dim(tphr)[2], dim(tphr)[3] ,byrow=T)
corP_array2 = array(corP_array , dim = c(dim(tphr)[1] ,dim(tphr)[2], dim(tphr)[3]))

ncvar_put( nc=nc1, varid='tp', vals=corP_array2, start=NA, count=NA, verbose=FALSE )

print(paste("timestep=", timestep))
print(paste("Original domain wide totals multiplied by timestep (rough scaling)=", round(sum(tphr)*timestep,0), "m"))
print(paste("Monthly totals=", round(sum(tpmm3),0), "m"))
print(paste("New totals after correction=", round(sum(corP_array2),0), "m"))
#print(mean(unlist(cfl), na.rm=T))
# check all good
#plot(corP_array[1,1:100], typ='l')
#lines(tph_flat[1,1:100], col='red')

#plot(corP_array2[1,1, 1:100], typ='l')
#lines(tphr[1,1, 1:100], col='red')





