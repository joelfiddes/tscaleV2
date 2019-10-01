# precip lapse rates

dem=raster("/home/joel/Projects/2018/aral_book/demaralsea.tif")
era5L=raster('/home/joel/data/era5_month_means/era5TP_mm2014_annualsum.nc',varname='tp')
shp=shapefile('/home/joel/data/HMA_SNO/input_hyperion/amu_shp7_FINAL.shp')

demc =crop(dem,shp)
era =crop(era5L,shp)
plot(demc)
plot(era)
plot(demc,add=T)

# a tricky east west fgradient
plot(shp[38,],add=T)

# shows resolution issues nicely
dem25=aggregate(demc,fact=25)
dem10=aggregate(demc,fact=10)



corvec=c()
for (n in 1:150){
	print(n)
	if(n==45 | n==55 | n==104){next}
plot(getValues(mask(era,shp[n,])), getValues(mask(d10,shp[n,])))
cv=cor(getValues(mask(era,shp[n,])), getValues(mask(d10,shp[n,])))
corvec=c(corvec,cv)
}

# lack of trend sign of E-