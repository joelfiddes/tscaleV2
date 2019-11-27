"""main

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

# DEBUG
# inDir= '/home/joel/sim/imis/forcing/'
# home='/home/joel/sim/imis/'
# outDir = inDir 
# startTime = '2000-09-01'
# endTime = '2001-09-01'
# windCor='FALSE'

t2m_file=inDir+"/t2m.nc"
d2m_file=inDir+"/d2m.nc"
zs_file=inDir+"/zs.nc"
ssrd_file=inDir+"/ssrd.nc"
strd_file=inDir+"/strd.nc"
tp_file=inDir+"/tp.nc"
tisr_file=inDir+"/tisr.nc"

zp_file=inDir+"/zp.nc"
t_file=inDir+"/t.nc"
u_file=inDir+"/u.nc"
v_file=inDir+"/v.nc"
r_file=inDir+"/r.nc"

lpfile=home + "/listpoints.txt"
# read in lispoints
lp=pd.read_csv(lpfile)

#===============================================================================
#	Logging
#===============================================================================
logging.basicConfig(level=logging.DEBUG, filename=home+"/tscale_logfile", filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")

for i in range(lp.id.size):

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
	zp = nc.Dataset(zp_file)
	t = nc.Dataset(t_file)
	u = nc.Dataset(u_file)
	v = nc.Dataset(v_file)
	r = nc.Dataset(r_file)
	
	# Datetime structure 
	nctime = t.variables['time']
	dtime = pd.to_datetime(nc.num2date(nctime[:],nctime.units, calendar="standard"))
	startIndex = int(np.where(dtime==startTime)[0])#"2016-08-01 18:00:00"
	endIndex = int(np.where(dtime==endTime)[0])#"2016-08-01 18:00:00"

	# cut p.dtime to start and end
	dtime=dtime[startIndex:endIndex]

	# extract data
	zpdat = zp.variables['z'][:, startIndex:endIndex,basin ] # dims = levels * time 
	tdat = t.variables['t'][ :, startIndex:endIndex,basin ]
	udat = u.variables['u'][ :, startIndex:endIndex,basin ]
	vdat = v.variables['v'][ :,startIndex:endIndex,basin ]
	rdat = r.variables['r'][ :,startIndex:endIndex,basin ]

	# get levels
	lev=zp.variables['levels']
	# create pob instance/Bunch
	pob = hp.Bunch(z=zpdat.T,t=tdat.T,u=udat.T,v=vdat.T,r=rdat.T, dtime=dtime, levels=lev)

#===============================================================================
#	Surface - Make a SOB
#===============================================================================
	# Pressure level object 
	t2m = nc.Dataset(t2m_file)
	d2m = nc.Dataset(d2m_file)
	zs = nc.Dataset(zs_file)
	ssrd = nc.Dataset(ssrd_file)
	strd = nc.Dataset(strd_file)
	tp = nc.Dataset(tp_file)
	tisr = nc.Dataset(tisr_file)

	# Datetime structure 
	nctime = t2m.variables['time']
	dtime = pd.to_datetime(nc.num2date(nctime[:],nctime.units, calendar="standard"))
	startIndex = int(np.where(dtime==startTime)[0])#"2016-08-01 18:00:00"
	endIndex = int(np.where(dtime==endTime)[0])#"2016-08-01 18:00:00"

	# cut p.dtime to start and end
	dtime=dtime[startIndex:endIndex]

	# extract data
	t2mdat = t2m.variables['t2m'][ startIndex:endIndex ,basin] 
	d2mdat = d2m.variables['d2m'][ startIndex:endIndex ,basin]
	zsdat = zs.variables['z'][ startIndex:endIndex ,basin]
	ssrddat = ssrd.variables['ssrd'][ startIndex:endIndex ,basin]
	strddat = strd.variables['strd'][ startIndex:endIndex,basin]
	tpdat = tp.variables['tp'][ startIndex:endIndex,basin]
	tisrdat = tisr.variables['tisr'][ startIndex:endIndex,basin]
	

	""" Convert SWin from accumulated quantities in J/m2 to 
	instantaneous W/m2 see: 
	https://confluence.ecmwf.int/pages/viewpage.action?pageId=104241513

	Args:
		step: timstep in seconds (era5=3600, ensemble=10800)

	Note: both EDA (ensemble 3h) and HRES (1h) are accumulated over the timestep
	and therefore treated here the same ie step = 3600s
	https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
				"""
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
	

	df = pd.DataFrame({	"TA":t.t, 
				"RH":t.r*0.01, #meteoio 0-1
				"VW":t.ws,
				"DW":t.wd, 
				"ILWR":t.LWf, 
				"ISWR":t.SWfglob, 
				"PINT":t.prate,
				"PSUM":t.psum
				},index=pob.dtime)
	df.index.name="datetime"

	# fill outstanding nan in SW routine with 0 (night)
	df.ISWR = df.ISWR.fillna(0)
	fileout=outDir+"/meteo"+"c"+str(i+1)+".csv"
	df.to_csv(path_or_buf=fileout ,na_rep=-999,float_format='%.3f')
	logging.info(fileout + " complete")


