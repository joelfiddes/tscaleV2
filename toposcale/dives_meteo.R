# rscript to convert tscale output to DIVES format
wd='/home/caduff/sim/dives/'
samples = 150
basins=list.files(paste0(wd,'/sim/'))

# be careful!
# psum currently divided by 3 as precip has already been multiplied by 3 in raw - this is a hack that only applies to current DIVES data
timestep =3

for (basinID in basins){

	for ( i in 1:samples){
	print(paste0(basinID,'/',i))
	dat=read.csv(paste0(wd,'/sim/', basinID, '/forcing/meteoc', i, '.csv'))

	#============= DATETIME =============================================================
	params = names(dat)
	# aggregation factors
	datetime=dat$datetime
	day_agg = substr(datetime,1,10)
	#npars=dim(sRes)[2]
	#nsamples=dim(sRes)[3]
	ndays=length(unique(day_agg))
	date=unique(day_agg)


	#=============== Aggregations ========================================================
	watts2Mjoules = (60*60)/1000000 # (sec*min) / 1000000 (j in a megajoule) watts = j/s


	TAMEAN=aggregate(dat$TA, list(day_agg), mean)[,2] # K
	TAMIN=aggregate(dat$TA, list(day_agg), min)[,2] # K
	TAMAX=aggregate(dat$TA, list(day_agg), max)[,2] # K
	ISWRSUM=aggregate(dat$ISWR*watts2Mjoules, list(day_agg), sum)[,2] *timestep  # mj in a day
	VWMEAN=aggregate(dat$VW, list(day_agg), mean)[,2] # m/s
	VWMAX=aggregate(dat$VW, list(day_agg), max)[,2] # m/s
	DWMEAN=aggregate(dat$DW, list(day_agg), mean)[,2] # m/s
	PSUM=aggregate(dat$PSUM, list(day_agg), sum)[,2] /timestep# mm
	RHMEAN=aggregate(dat$RH, list(day_agg), mean)[,2] # 0-1

	# calc saturated vapour pressure (pascals)
	# T in deg 
	T=TAMEAN-273.13
	svp = 610.7*10^(7.5*T/(237.3+T))

	# As VPD is the saturated vapour pressure minus the actual vapour pressure (SVP - VPactual), and VPactual = (RH*SVP)/100
	# RH in 0-100
	RH=RHMEAN*100
	VPD = ((100 - RH)/100)*svp

	VPDMEAN = VPD/1000 # pa to kpa


	df=data.frame(date, TAMIN,TAMAX, ISWRSUM,VWMAX,DWMEAN, VPDMEAN,PSUM)

	write.csv(df, paste0(wd,'/sim/', basinID, '/forcing/',basinID,'_samp', i, '.csv'),quote=FALSE, row.names=FALSE)
	}
}