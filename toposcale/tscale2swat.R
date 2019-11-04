# run points spasre with basin centroids
# read output from pints sparse
# convert to swat
require(raster) 
require(getMet)
# this sniippet generates sparse points file for 
require(rgeos)
require(sp)
swat_subbasins <- shapefile("/home/joel/data/zisp/sim/Pete2/Watershed/Shapes/subs1.shp")
crs_ll <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
swat_subbasins_ll <-spTransform(swat_subbasins, crs_ll )
sparse_points <- gCentroid(swat_subbasins_ll ,byid=T)
shapefile(sparse_points, "sparse_points.shp", overwrite=T)

# map basins and points
# swat basins are not ordered sequentially in shp file therefore toposcale output meteoc1.csv 
# does not equt to basin 1 but feature 1 of the sub.shp
# In our example here subbasin 1 is actually feature 3
#R> swat_subbasins$Subbasin                                                                                                  
# [1] "5" "4" "1" "2" "3"   

# so here we log this 'swat order' to use later
swat_station_order = as.numeric(swat_subbasins$Subbasin)  
#Need library called "ncdf4" 

#=======================================================
# cfsr data
#getSWATcfsr(centroids, outDir = getwd())
#========================================================
wd="/home/joel/swatr_demo/sim/"
outDir="/home/joel/swatr_demo/p5.1/Backup/"
#shp= shapefile("/home/joel/centroids_wgs.shp") points file used by toposcale
swat_subbasins=shapefile("/home/joel/data/zisp/sim/Pete2/Watershed/Shapes/subs1.shp")# SWAT generated subbasins polygon
files = list.files(wd)

dat = read.csv(paste0(wd, "/",files[1], "/forcing/meteoc1.csv"))
myday=substr(dat$datetime,1,10)

umyday = unique(myday)
require(lubridate)
x = as.Date(umyday)
jday = formatC(yday(x),width=3, flag='0')

year=year(x)
swatDates=paste0(year,jday)

df=c(swatDates)
eleVec=c()


# init lists
    tmpls <- list()
    pcpls <- list()
    hmdls <- list()
    slrls <- list()
    wndls <- list()

for (myindex in 1:5) {


  #This finds which feature index basin "myindex" actually is
  # in this example basin 1 is in swat_station_order 3
  # R> as.numeric(swat_station_order)                                                                                           
  # [1] 5 4 1 2 3 
  # therefore we use meteoc3.csv for basin 1
  a = which(swat_station_order==myindex) 

    # get elevation
    ele <- read.csv(paste0(wd, "/",files[a], "/listpoints.txt"))$ele
    eleVec=c(eleVec, ele)

    #getMet::getSWATcfsr
    dat <- read.csv(paste0(wd, "/",files[a], "/forcing/meteoc1.csv"))
    date.full <- swatDates

    # TMAX
    TMX=aggregate(dat$TA,list(myday), max)
    tmax.format <- data.frame(date.full, TMX$x)
    colnames(tmax.format) <- c("date", "tmax")
    tmax.format$tmax <- formatC(tmax.format$tmax, digits = 1, width = 5, flag = "0", format = "f")

    # TMIN
    TMN=aggregate(dat$TA,list(myday), min)
    tmin.format <- data.frame(tmax.format$date, TMN$x)
    colnames(tmin.format) <- c("date", "tmin")
    tmin.format$tmin <- formatC(tmin.format$tmin, digits = 1, width = 5, flag = "0", format = "f")
   
    # join dataframes
    temp.format <- data.frame(tmax.format, tmin.format$tmin)
    colnames(temp.format) <- c("date", "tmax", "tmin")


    #RHUM spec: https://swat.tamu.edu/media/69329/ch10_input_hmd.pdf
    #        Daily mean RH [0-1]
    RH=aggregate(dat$RH,list(myday), mean)
    rhum.format <- data.frame(tmax.format$date, RH$x)
    colnames(rhum.format) <- c("date", "rhum")
    rhum.format$rhum <- formatC(rhum.format$rhum, digits = 3, 
        width = 8, flag = "0", format = "f")
    colnames(rhum.format) <- c("date", "rhum")

    # PSUM
    PSUM=aggregate(dat$PSUM,list(myday), sum)
    pcp.format <- data.frame(tmax.format$date, PSUM$x*6)
    colnames(pcp.format) <- c("date", "pcp")
    pcp.format$pcp <- formatC(pcp.format$pcp, digits = 1, 
        width = 5, flag = "0", format = "f")
    colnames(pcp.format) <- c("date", "pcp")

    # wind
    VW=aggregate(dat$VW,list(myday), mean)
    wnd.format <- data.frame(tmax.format$date, VW$x)
    colnames(wnd.format) <- c("date", "wnd")
    wnd.format$wnd <- formatC(wnd.format$wnd, digits = 3, 
        width = 8, flag = "0", format = "f")
    colnames(wnd.format) <- c("date", "wnd")
    
    # ISWR
    # total daily Mj/m2

    # convert watts/m2 to MJ/s / m2
    
    ISWR_mj = dat$ISWR*10^-6
    secsinaday=24*60*60
    ISWR_mjsec=aggregate(ISWR_mj,list(myday), mean)
    ISWR = ISWR_mjsec$x*secsinaday

    slr.format <- data.frame(tmax.format$date, ISWR)
    colnames(slr.format) <- c("date", "slr")
    slr.format$slr <- formatC(slr.format$slr, digits = 3, 
        width = 8, flag = "0", format = "f")

    colnames(slr.format) <- c("date", "slr")
    tmpls[[a]] <- temp.format
    pcpls[[a]] <- pcp.format
    hmdls[[a]] <- rhum.format
    slrls[[a]] <- slr.format
    wndls[[a]] <- wnd.format
    tmpls[[a]][, 1] <- NULL
    tmpls[[a]]$full <- paste(tmpls[[a]][, 1], tmpls[[a]][, 
        2], sep = "")
    tmpls[[a]]$tmax <- NULL
    tmpls[[a]]$tmin <- NULL
    pcpls[[a]][, 1] <- NULL
    hmdls[[a]][, 1] <- NULL
    slrls[[a]][, 1] <- NULL
    wndls[[a]][, 1] <- NULL
}

alltmps <- do.call("cbind", tmpls)
allpcps <- do.call("cbind", pcpls)
allhmds <- do.call("cbind", hmdls)
allslrs <- do.call("cbind", slrls)
allwnds <- do.call("cbind", wndls)
tmpprint <- data.frame(tmax.format$date, alltmps)
pcpprint <- data.frame(tmax.format$date, allpcps)
hmdprint <- data.frame(tmax.format$date, allhmds)
slrprint <- data.frame(tmax.format$date, allslrs)
wndprint <- data.frame(tmax.format$date, allwnds)
Lati <- sparse_points@coords[,2]
Long <- sparse_points@coords[,1]

# create header note:
# only tmp and pcp have lat/lon/ele hdrs
headtmp <- c("Station  tmp degC \tSource era5", 
  paste("Lati", paste(round(Lati,5), collapse = "   "), sep = "   "), 
  paste("Long", paste(round(Long,5),  collapse = "   "), sep = "   "), 
      paste("Elev", paste(eleVec,  collapse = "   "), sep = "   "))
headpcp <- c("Station  pcp mm/day \tSource era5", 
  paste("Lati", paste(round(Lati,5), collapse = "   "), sep = "   "), 
  paste("Long", paste(round(Long,5),  collapse = "   "), sep = "   "), 
  paste("Elev", paste(eleVec,  collapse = "   "), sep = "   "))
headhmd <- "Relative Humidity % \tSource era5"
headwnd <- "Wind Speed m/s \tSource era5"
headslr <- "Solar Radiation MJ/m^2 \tSource era5"

setwd(outDir)
writeLines(headtmp, "tmp1.tmp")
writeLines(headpcp, "pcp1.pcp")
writeLines(headhmd, "hmd.hmd")
writeLines(headwnd, "slr.slr")
writeLines(headslr, "wnd.wnd")
write.table(tmpprint, "tmp1.tmp", append = TRUE, sep = "", 
    row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(pcpprint, "pcp1.pcp", append = TRUE, sep = "", 
    row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(hmdprint, "hmd.hmd", append = TRUE, sep = "", 
    row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(slrprint, "slr.slr", append = TRUE, sep = "", 
    row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(wndprint, "wnd.wnd", append = TRUE, sep = "", 
    row.names = FALSE, col.names = FALSE, quote = FALSE)
























