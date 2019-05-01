require(hydroGOF)

lp=read.csv("/home/joel/sim/imis/listpoints.txt")
d1=read.csv("/home/joel/sim/imis/out/meteoc1.csv") # 1d case
d3=read.csv("/home/joel/sim/imis/out/meteoc1.csv") # 3d case

# read in an cut obs
nand=read.table(paste0("~/sim/snowpack/input/WFJ_optimaldataset_v7.smet"),sep='', skip=46) #v7
names(nand) <- c("timestamp", "TA", "RH", "VW", "DW", "ISWR", "OSWR", "ILWR", "OLWR", "PSUM", "HS_M", "HS", "TSG", "TSS", "LYSI" ,"SWE_M")
start = which(nand$timestamp=="2015-09-01T00:00")
end = which(nand$timestamp=="2016-09-01T00:00")
wfj_cut = nand[start:end,]

# aggregation vec obs to 3 h step
agg3h =rep(1:(length(wfj_cut$timestamp)/48), each=48)
wfj_cut = nand[start:(end-1),]
# compute vars
TA = aggregate(wfj_cut$TA,by=list(agg3h), FUN="mean")
TA[TA < 200] <- NA
RH = aggregate(wfj_cut$RH,by=list(agg3h), FUN="mean")*100
RH[RH < 0] <- NA
VW = aggregate(wfj_cut$VW,by=list(agg3h), FUN="mean")
DW = aggregate(wfj_cut$DW,by=list(agg3h), FUN="mean")
ISWR = aggregate(wfj_cut$ISWR,by=list(agg3h), FUN="mean")
ILWR = aggregate(wfj_cut$ILWR,by=list(agg3h), FUN="mean")
PSUM = aggregate(wfj_cut$PSUM,by=list(agg3h), FUN="sum") # 3h TOTAL
PSUM[PSUM<0]<-NA
HS = aggregate(wfj_cut$HS,by=list(agg3h), FUN="mean")
obs = data.frame(TA$x,RH$x, DW$x,VW$x,ISWR$x,ILWR$x,PSUM$x, HS$x)

lwd=1
plot(d1$TA, type='l', col='green',lwd=lwd)
#lines(d3$TA, type='l', col='red',lwd=lwd)
lines(obs$TA, col='black',lwd=lwd)
lwd=3
legend("bottomright", c("1d", "3d", "obs"), col=c("green", "red", "black"), lty=1,lwd=lwd)


orig=read.csv("/home/joel/sim/snowpack/input/WFJ_ERA5_2.smet", sep='\t', skip=12)
names(orig) <- c("timestamp", "PSUM",  "VW", "DW", "RH", "TA", "ISWR", "ILWR") 

 obs$HS.x[obs$HS.x< -5]<-NA


plot(obs$HS.x, col='red', ylim=c(-10,400))
lines(m$hs_mod*100)
abline(h=0)
legend("topright", c("toposcale","obs"), col=c("black", "red"), lty=1,lwd=lwd)


