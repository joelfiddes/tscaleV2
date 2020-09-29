"""main

# this file interpolates forcing from PLEV / SURF to centroid of each sample, thereby accounting for bias in distributuions of samples ie flat desert in North, mountain ridge in south.

Main toposcale run script. Requires existence of:
	- pressurelevel data file: "PLEV.nc"
	- surface data file: "/SURF.nc"
	- listpoints file: "/listpoints.txt"

	- runs on however many pressure levels are downloaded

Object description:

	* pob - pressure level object
	* sob - surface object
	* tob - toposcale object
	* stat - station object (points, cluster centroids or grid centroids)

Example:
	python /home/joel/src/tscaleV2/toposcale/tscale_run.py /home/joel/sim/imis/forcing/ /home/joel/sim/imis/ /home/joel/sim/imis/forcing/ 2000-09-01 2001-09-01 FALSE
	
Args:
	inDir: directory containing input meteo PLEV.nc and SURF.nc 
		"/home/joel/sim/topomapptest/forcing/"
	home: location of listpoints.txt
		"/home/joel/sim/topomapptest/sim/g1m1/"
	outDir: location where output meteo files are written
		"/home/joel/sim/topomapptest/sim/g1m1/forcing/"
	startTime: ISO format #"2016-08-01 18:00:00"
	endTime: ISO format #"2016-08-01 18:00:00"
	windCor: use sebs wind correction str: "TRUE" or "FALSE"
Todo:

"""

import pandas as pd # only required to import listpoints nicely, should be able to remove
import era5 as e5
import tscale as ts
import helper as hp
import numpy as np
import sys
import os
import logging
import netCDF4 as nc
from tqdm import tqdm
import subprocess

#=== ARGS==============================================
inDir=sys.argv[1] # /home/joel/sim/imis/forcing/
home=sys.argv[2] # /home/joel/sim/imis/sim/g1
outDir = sys.argv[3] # /home/joel/sim/imis/sim/g1/forcing
startTime = sys.argv[4]
endTime = sys.argv[5]
windCor=sys.argv[6]
basin=sys.argv[7] 
plapse=sys.argv[8]

basin = int(basin) -1  #python index

# # DEBUG
# cp paste cmd line arg to ipython

# sys.argv = '/home/joel/src//tscaleV2/toposcale/tscale_run_basin_interp.py', '/home/joel/sim/paiku_basin/forcing', '/home/joel/sim/paiku_basin/sim/g14', '/home/joel/sim/paiku_basin/sim/g14/forcing/', '2013-09-01 00:00:00', '2014-09-01 00:00:00', 'FALSE', '14', 'FALSE'

# inDir=sys.argv[1] # /home/joel/sim/imis/forcing/
# home=sys.argv[2] # /home/joel/sim/imis/sim/g1
# outDir = sys.argv[3] # /home/joel/sim/imis/sim/g1/forcing
# startTime = sys.argv[4]
# endTime = sys.argv[5]
# windCor=sys.argv[6]
# basin=sys.argv[7] 
# plapse=sys.argv[8]
# 



lpfile=home + "/listpoints.txt"
# read in lispoints
lp=pd.read_csv(lpfile)

#===============================================================================
#	Logging
#===============================================================================
logging.basicConfig(level=logging.DEBUG, filename=home+"/tscale_logfile", filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")

for i in tqdm(range(lp.id.size)):

	fileout=outDir+"/meteo"+"c"+str(i+1)+".csv"
	if os.path.exists(fileout):
		logging.info(fileout + " exists!")
		continue

	# station attribute structure , tz always =0 for case of ERA5
	stat = hp.Bunch(ele = lp.ele[i], slp = lp.slp[i],asp = lp.asp[i],svf = lp.svf[i],lon = lp.lon[i], lat =lp.lat[i],sro = lp.surfRough[i],tz = lp.tz[i]  )

	#===============================================================================
	#	PLEVEL - mak a POB
	#===============================================================================
	# Pressure level object 
	plev=nc.Dataset(inDir+"/PLEV.nc")
	
	# Datetime structure 
	nctime = plev.variables['time']
	dtime = pd.to_datetime(nc.num2date(nctime[:],nctime.units, calendar="standard"))
	startIndex = int(np.where(dtime==startTime)[0])#"2016-08-01 18:00:00"
	endIndex = int(np.where(dtime==endTime)[0])#"2016-08-01 18:00:00"

	# cut p.dtime to start and end
	dtime=dtime[startIndex:endIndex]



	# bilin interp to coords centroids of sample
	cmd = ["cdo -b F64 -remapbil,lon="+str(stat.lon)+"_lat="+str(stat.lat)+" "+ inDir+"/PLEV.nc "+ outDir+"/PLEVtmp.nc"]
	subprocess.check_output(cmd, shell = "TRUE")
	plevtemp=nc.Dataset(outDir+"/PLEVtmp.nc")

	# extract data to test for neighbouring clls and therefor we best interpolation method (NN or bilin)
	zpdat_test = plevtemp.variables['z'][ startIndex:endIndex,:,:,: ].squeeze().transpose() # dims = levels * time 
	
	# here we test for masked values which is sign of not enough points for bilinear interp, use NN instaed
	if np.ma.is_masked(zpdat_test)==True :
		print 'remapbil fails due to not enough grid points (edge of era5 grid), trying NN interp instead'
			# bilin interp to coords centroids of sample
		cmd = ["cdo -b F64 -remapnn,lon="+str(stat.lon)+"_lat="+str(stat.lat)+" "+ inDir+"/PLEV.nc "+ outDir+"/PLEVtmp.nc"]
		subprocess.check_output(cmd, shell = "TRUE")
		plevtemp=nc.Dataset(outDir+"/PLEVtmp.nc")

		# test data again
		zpdat_testagain = plevtemp.variables['z'][ startIndex:endIndex,:,:,: ].squeeze().transpose() # dims = levels * time 
		if np.ma.is_masked(zpdat_testagain)==True :
			print "NN failed!"
		if np.ma.is_masked(zpdat_testagain)==False :
			print "NN worked!"

	# now we extract all data from bilin or NN interpolated datset
	zpdat = plevtemp.variables['z'][ startIndex:endIndex,:,:,: ].squeeze().transpose() # dims = levels * time
	tdat = plevtemp.variables['t'][ startIndex:endIndex,:,:,: ].squeeze().transpose()  
	udat = plevtemp.variables['u'][ startIndex:endIndex,:,:,: ].squeeze().transpose()  
	vdat = plevtemp.variables['v'][ startIndex:endIndex,:,:,: ].squeeze().transpose()  
	rdat = plevtemp.variables['r'][ startIndex:endIndex,:,:,: ].squeeze().transpose()  

	# get levels
	lev=plev.variables['level']
	# create pob instance/Bunch
	pob = hp.Bunch(z=zpdat.T,t=tdat.T,u=udat.T,v=vdat.T,r=rdat.T, dtime=dtime, levels=lev)

	#===============================================================================
	#	Surface - Make a SOB
	#===============================================================================
	# Pressure level object 
	surf=nc.Dataset(inDir+"/SURF.nc")

	# Datetime structure 
	nctime = surf.variables['time']
	dtime = pd.to_datetime(nc.num2date(nctime[:],nctime.units, calendar="standard"))
	startIndex = int(np.where(dtime==startTime)[0])#"2016-08-01 18:00:00"
	endIndex = int(np.where(dtime==endTime)[0])#"2016-08-01 18:00:00"

	# cut p.dtime to start and end
	dtime=dtime[startIndex:endIndex]



	# test for bilin again for surf data
	if np.ma.is_masked(zpdat_test)==True :
		print 'remapbil fails due to not enough grid points (edge of era5 grid), trying NN interp instead'
			# bilin interp to coords centroids of sample
		cmd = ["cdo -b F64 -remapnn,lon="+str(stat.lon)+"_lat="+str(stat.lat)+" "+ inDir+"/SURF.nc "+ outDir+"/SURFtmp.nc"]
		subprocess.check_output(cmd, shell = "TRUE")
		surftemp=nc.Dataset(outDir+"/SURFtmp.nc")

	if np.ma.is_masked(zpdat_test)==False :
		# bilin interp to coords centroids of sample
		cmd = ["cdo -b F64 -remapbil,lon="+str(stat.lon)+"_lat="+str(stat.lat)+" "+ inDir+"/SURF.nc "+ outDir+"/SURFtmp.nc"]
		subprocess.check_output(cmd, shell = "TRUE")
		surftemp=nc.Dataset(outDir+"/SURFtmp.nc")

	# extract data
	t2mdat = surftemp.variables['t2m'][ startIndex:endIndex,:,:].squeeze().transpose()  
	d2mdat = surftemp.variables['d2m'][ startIndex:endIndex,:,:].squeeze().transpose()  
	zsdat = surftemp.variables['z'][ startIndex:endIndex,:,:].squeeze().transpose()  
	ssrddat = surftemp.variables['ssrd'][ startIndex:endIndex,:,:].squeeze().transpose()  
	strddat = surftemp.variables['strd'][ startIndex:endIndex,:,: ].squeeze().transpose()  
	tpdat = surftemp.variables['tp'][ startIndex:endIndex,:,: ].squeeze().transpose()  
	tisrdat = surftemp.variables['tisr'][ startIndex:endIndex,:,: ].squeeze().transpose()  
	

	""" Convert SWin from accumulated quantities in J/m2 to 
	instantaneous W/m2 see: 
	https://confluence.ecmwf.int/pages/viewpage.action?pageId=104241513

	Args:
		step: timstep in seconds (era5=3600, ensemble=10800)

	Note: both EDA (ensemble 3h) and HRES (1h) are accumulated over the timestep
	and therefore treated here the same ie step = 3600s
	https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation"""
	strddatI = strddat/3600  
	ssrddatI = ssrddat/3600
	tisrdatI = tisrdat/3600 


	""" convert tp from m/timestep (total accumulation over timestep) to rate in mm/h 

		Args:
			step: timstep in seconds (era5=3600, ensemble=10800)

		Note: both EDA (ensemble 3h) and HRES (1h) are accumulated over the timestep
		however in prprocessing TP is scaled to account for missing sums. Therefor tp is sum over actual timestep (1,3,6h)
		https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
	"""
	# compute step in seconds for accumulated surface fields
	a=dtime[2]-dtime[1]
	step = a.seconds
	stepinhr=step/(60*60)
	
	prate = tpdat*1000 # convert m per timestep -> mm/hour [PRATE] use approx scaling method until monthly correction fixed. 
	# Although this is actually pretty accurate. tp is accumulated over hour so prate IS tp at coarser timstep (eg 6h) we introduce uncertainty
	psum = prate* stepinhr # mm/ timestep  [PSUM] scale by stepinhr

	""" compute surface elevation of coarse grid"""
	gridEle = zsdat/9.80665


	# create pob instance/Bunch
	sob = hp.Bunch(t2m=t2mdat,d2m=d2mdat,z=zsdat,ssrd=ssrddatI,strd=strddatI,prate=prate,psum=psum, tisr=tisrdatI,gridEle=gridEle, dtime=dtime)


	#=== toposcale object ====================================================

	"""init object"""
	t = ts.tscale(pob.z)

	""" Downscale pl fields """
	for v in ['t','u','v','r']:
		dat = getattr(pob, v) # allows dynamic call to instance variable
		t.tscale1D(dat,stat)
		t.addVar(v,t.interpVar)
	

	# constrain r to interval 5-100 - do here as required by LWin parameterisation
	t.r[t.r <5]=5
	t.r[t.r>100]=100

	tob = t
	# compute downscaled longwave (LWf) Wm**2
	t.lwin(sob, tob,stat)

	# compute downscaled shortwave (SWf) Wm**2m	
	t.swin(pob, sob,tob, stat,pob.dtime)

	# compute downscaled precipitation in mm/h (PRATE) and m/step (PSUM)
	t.precip(sob,stat, plapse)

	t.wind(tob)
	t.ws[t.ws<0]=0

	if windCor=="TRUE":
		# corrected wind
		blend = 40 # blend height
		statblend=stat # new instance of stat for correction only 
		statblend.ele = statblend.ele + blend # add blend to stat ele
		v='v'
		dat = getattr(p, v) # allows dynamic call to instance variable
		t.tscale1D(dat,statblend)
		t.addVar('vblend',t.interpVar)

		v='u'
		dat = getattr(p, v) # allows dynamic call to instance variable
		t.tscale1D(dat,statblend)
		t.addVar('ublend',t.interpVar)

		# ws from vectors
		t.wsblend = np.sqrt(tob.ublend**2+tob.vblend**2)

		t.windCorRough( tob, pob, sob, stat)

	 # plt.plot(t.SWfdirCor)
	# plt.show()
		# partition prate to rain snow (mm/hr)
	lowthresh=272.15
	highthresh = 274.15
	d = {'prate': t.prate, 'ta': t.t }
	df = pd.DataFrame(data=d)
	snow = df.prate.where(df.ta<lowthresh) 
	rain=df.prate.where(df.ta>=highthresh) 

	mix1S = df.prate.where((df.ta >= lowthresh) & (df.ta<=highthresh), inplace = False)
	mix1T = df.ta.where((df.ta >= lowthresh) & (df.ta<=highthresh), inplace = False)
	mixSno=(highthresh - mix1T) / (highthresh-lowthresh)
	mixRain=1-mixSno
	addSnow=mix1S*mixSno
	addRain=mix1S*mixRain
	
	# nas to 0
	snow[np.isnan(snow)] = 0 	
	rain[np.isnan(rain)] = 0 
	addRain[np.isnan(addRain)] = 0 
	addSnow[np.isnan(addSnow)] = 0 

	# linearly reduce snow to zero in steep slopes
	#if steepSnowReduce=="TRUE": # make this an option if need that in future
	snowSMIN=30.
	snowSMAX=80.
	slope=stat.slp

	k= (snowSMAX-slope)/(snowSMAX - snowSMIN)

	if slope<snowSMIN:
		k=1
	if slope>snowSMAX:
		k=0

	t.snowTot=(snow+addSnow) * k
	t.rainTot=rain + addRain


	df = pd.DataFrame({	"TA":t.t, 
				"RH":t.r*0.01, #meteoio 0-1
				"VW":t.ws,
				"DW":t.wd, 
				"ILWR":t.LWf, 
				"ISWR":t.SWfglob, 
				"PINT":t.prate,
				"PSUM":t.psum,
				"P":t.psf,
				"Rf":np.array(t.rainTot),
				"Sf":np.array(t.snowTot)
				},index=pob.dtime)
	df.index.name="datetime"

	# fill outstanding nan in SW routine with 0 (night)
	df.ISWR = df.ISWR.fillna(0)
	fileout=outDir+"/meteo"+"c"+str(i+1)+".csv"
	df.to_csv(path_or_buf=fileout ,na_rep=-999,float_format='%.3f')
	logging.info(fileout + " complete")


