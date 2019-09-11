# This routine compliments the basin/ chirps tsub approach by bilinear interpolating p timeseries from netcdf (era5) to the centre (long lat) of the sample. Samples now have a centre in space! Due to clustering on Chirps climatology and calc mean lon lat.
require(ncdf4)


surfnc = "/home/joel/sim/tsub_basinSml/forcing/SURF.nc"
Nclust=50
sim="/home/joel/sim/tsub_basinSml/sim/g9"
lp=read.csv(paste0(sim,"/listpoints.txt"))
start="2013-09-01T00:00"
end="2014-09-01T00:00"

for (i in 1:length(lp$lonRst)){

mylon=lp$longRst[i]
mylat=lp$latRst[i]

outfile=paste0(mylon,"_",mylat,".nc")
cmd = paste0("cdo intgridbil,lon=",mylon,"_lat=",mylat," ",surfnc," "  ,outfile,"_tmp.nc")
system(cmd)

cmd=paste0("cdo seldate," ,start,",",end,  " ",outfile,"_tmp.nc", " ",outfile)
system(cmd)

nc = nc_open(outfile)
tp = ncvar_get(nc, 'tp')

tpcut = tp[1:(length(tp)-1)] # drop the last timestamp as our dates always cut to 08-31T06:00 not 09-01T00:00

met= read.csv(paste0("/home/joel/sim/tsub_basinSml/sim/g9/forcing/meteoc",i,".csv"))
met$PINT<-tpcut*6
write.csv(met,paste0("/home/joel/sim/tsub_basinSml/sim/g9/forcing/meteoc",i,".csv"))
print(paste0("replaced meteoc",i,".csv"))
} 

#replace .csv this is then process to csv.gtp by script just need to run "sim" part



