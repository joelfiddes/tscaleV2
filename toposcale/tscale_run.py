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
import logging
#=== ARGS==============================================
inDir=sys.argv[1] # /home/joel/sim/imis/forcing/
home=sys.argv[2] # /home/joel/sim/imis/
outDir = sys.argv[3] # 
startTime = sys.argv[4]
endTime = sys.argv[5]
windCor=sys.argv[6]

# DEBUG
# inDir= '/home/joel/sim/imis/forcing/'
# home='/home/joel/sim/imis/'
# outDir = inDir 
# startTime = '2000-09-01'
# endTime = '2001-09-01'
# windCor='FALSE'

fp=inDir+"/PLEV.nc"
fs=inDir+"/SURF.nc"
lpfile=home + "/listpoints.txt"
# read in lispoints
lp=pd.read_csv(lpfile)

#===============================================================================
#	Logging
#===============================================================================
logging.basicConfig(level=logging.DEBUG, filename=home+"/tscale_logfile", filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")

for i in range(lp.id.size):

	# station attribute structure , tz always =0 for case of ERA5
	stat = hp.Bunch(ele = lp.ele[i], slp = lp.slp[i],asp = lp.asp[i],svf = lp.svf[i],lon = lp.lon[i], lat =lp.lat[i],sro = lp.surfRough[i],tz = lp.tz[i]  )

	#=== Pressure level object =============================================

	""" preprocess pressure level fields and add to structure """
	p=e5.Plev(fp, stat.lat, stat.lon)
	p.getVarNames()

	""" Datetime structure """
	p.addTime()
	startIndex = int(np.where(p.dtime==startTime)[0])#"2016-08-01 18:00:00"
	endIndex = int(np.where(p.dtime==endTime)[0])#"2016-08-01 18:00:00"

	"""cut p.dtime to start and end"""
	p.dtime=p.dtime[startIndex:endIndex]

	for v in p.varnames:

		#v=p.varnames[1]
		#p.getVar(v)
		#name = p.myvar.name

		p.extractCgc(v, startIndex, endIndex) # adds data to self

		p.addVar(v,p.var) # adds data again with correct name - redundancy
		#can now remove p.var it is redundant

	""" dimensions of data """
	p.addShape()



	""" Pressure level values """
	p.plevels()

	pob = p
	#=== Surface object ====================================================

	""" preprocess surface level fields and add to structure """
	s=e5.Surf(fs, stat.lat, stat.lon)
	s.getVarNames()

	""" Datetime structure """
	s.addTime()
	"""cut p.dtime to start and end"""
	s.dtime=s.dtime[startIndex:endIndex]

	# compute step in seconds for accumulated surface fields
	a=s.dtime[2]-s.dtime[1]
	step = a.seconds

	for v in s.varnames:

		#v=p.varnames[1]
		#print v
		s.getVar(v)
		#name = p.myvar.name

		s.extractCgc(v ,startIndex, endIndex) # adds data to self

		s.addVar(v,s.var) # adds data again with correct name - redundancy

	""" rad conversions era5 always accumulated over 1h or 3600s """
	s.instRad(3600)

	""" precip conversions m/timestep to mm/h"""
	s.tp2rate(step)

	""" dimensions of data """
	s.addShape()
	""" Datetime structure """
	#s.addTime()
	""" surface elevation of coarse grid"""
	s.gridEle()
	sob = s

	#=== toposcale object ====================================================

	"""init object"""
	t = ts.tscale(p.z)

	""" Downscale pl fields """
	for v in p.varnames:
		if (v=='z'):
			continue
		#print v
		dat = getattr(p, v) # allows dynamic call to instance variable

		t.tscale1D(dat,stat)
		t.addVar(v,t.interpVar)
	

	# constrain r to interval 5-100 - do here as required by LWin parameterisation
	t.r[t.r <5]=5
	t.r[t.r>100]=100

	tob = t
	# compute downscaled longwave (LWf) Wm**2
	t.lwin(sob, tob, stat)

	# compute downscaled shortwave (SWf) Wm**2
	t.swin(pob, sob,tob, stat,p.dtime)

	# compute downscaled precipitation in mm/h (PRATE) and m/step (PSUM)
	t.precip(sob,stat)

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
				},index=p.dtime)
	df.index.name="datetime"

	# fill outstanding nan in SW routine with 0 (night)
	df.ISWR = df.ISWR.fillna(0)
	fileout=outDir+"/meteo"+"c"+str(i+1)+".csv"
	df.to_csv(path_or_buf=fileout ,na_rep=-999,float_format='%.3f')
	logging.info(fileout + " complete")


