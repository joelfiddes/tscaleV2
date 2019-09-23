# rscript to convert tscale output to DIVES format
wd='/home/caduff/sim/dives/'

basins=list.files(paste0(wd,'/sim/'))

for (basinID in basins){
samples = list.files(paste0(wd,'/sim/', basinID), pattern='^c00*')


for ( i in 1:length(samples)){
print(i)
dat=read.csv(paste0(wd,'/sim/', basinID, '/forcing/meteoc', i, '.csv'))




#===============================================================================
#		DIVES CASE BY CASE
#===============================================================================
params = names(dat)
# aggregation factors
datetime=dat$datetime
day_agg = substr(datetime,1,10)
#npars=dim(sRes)[2]
#nsamples=dim(sRes)[3]
ndays=length(unique(day_agg))
date=unique(day_agg)


#=============== DW MEAN ========================================================
watts2joules = 0.0036

DWMEAN=aggregate(dat$DW, list(day_agg), mean)[,2] # m/s
TAMEAN=aggregate(dat$TA, list(day_agg), mean)[,2] # K
TAMIN=aggregate(dat$TA, list(day_agg), min)[,2] # K
TAMAX=aggregate(dat$TA, list(day_agg), max)[,2] # K
ISWRSUM=aggregate(dat$ISWR*watts2joules, list(day_agg), sum)[,2] # joules/m2
VWMEAN=aggregate(dat$ISWR, list(day_agg), mean)[,2] # m/s
VWMAX=aggregate(dat$ISWR, list(day_agg), max)[,2] # m/s
PSUM=aggregate(dat$PSUM, list(day_agg), sum)[,2] # mm
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


df=data.frame(date, TAMIN,TAMAX, ISWRSUM,VWMAX,VPDMEAN,PSUM)

write.csv(df, paste0(wd,'/sim/', basinID, '/forcing/',basinID,'_samp', i, '.csv'),quote=FALSE, row.names=FALSE)
}}
