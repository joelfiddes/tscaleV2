"""main

Main toposcale run script. Requires existence of:
	- pressurelevel data file: "PLEV.nc"
	- surface data file: "/SURF.nc"
	- listpoints file: "/listpoints.txt"

Object description:

	* pob - pressure level object
	* sob - surface object
	* tob - toposcale object
	* stat - station object (points, cluster centroids or grid centroids)

Example:

	python  tscale_run.py
				/home/joel/sim/topomapptest/forcing/
				/home/joel/sim/topomapptest/sim/g1m1/
				/home/joel/sim/topomapptest/sim/g1m1/forcing/
Args:
	inDir: directory containing input meteo PLEV.nc and SURF.nc 
		"/home/joel/sim/topomapptest/forcing/"
	listpointsLoc: location of listpoints.txt
		"/home/joel/sim/topomapptest/sim/g1m1/"
	outDir: location where output meteo files are written
		"/home/joel/sim/topomapptest/sim/g1m1/forcing/"
	startTime: ISO format #"2016-08-01 18:00:00"
	endTime: ISO format #"2016-08-01 18:00:00"
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
inDir=sys.argv[1] # /home/joel/sim/topomapptest/forcing/PLEV.nc
listpointsLoc=sys.argv[2]
outDir = sys.argv[3]
startTime = sys.argv[4]
endTime = sys.argv[5]

fp=inDir+"/PLEV.nc"
fs=inDir+"/SURF.nc"
lpfile=listpointsLoc + "/listpoints.txt"
# read in lispoints
lp=pd.read_csv(lpfile)

#===============================================================================
#	Logging
#===============================================================================
logging.basicConfig(level=logging.DEBUG, filename=listpointsLoc+"/tscale_logfile", filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")

for i in range(lp.id.size):

	# station attribute structure
	stat = hp.Bunch(ele = lp.ele[i], slp = lp.slp[i],asp = lp.asp[i],svf = lp.svf[i],lon = lp.lon[i], lat =lp.lat[i],sro = lp.surfRough[i],tz = lp.tz[i]  )

	#=== Pressure level object =============================================

	""" preprocess pressure level fields and add to structure """
	p=e5.Plev(fp, stat.lat, stat.lon)
	p.getVarNames()

	""" Datetime structure """
	p.addTime()
	startIndex = int(np.where(p.dtime==startTime)[0])#"2016-08-01 18:00:00"
	endIndex = int(np.where(p.dtime==endTime)[0])#"2016-08-01 18:00:00"


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

	for v in s.varnames:

		#v=p.varnames[1]
		#print v
		s.getVar(v)
		#name = p.myvar.name

		s.extractCgc(v ,startIndex, endIndex) # adds data to self

		s.addVar(v,s.var) # adds data again with correct name - redundancy

	""" rad conversions """
	s.instRad()
	""" dimensions of data """
	s.addShape()
	""" Datetime structure """
	s.addTime()
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
	tob = t

	# compute downscaled longwave (LWf) Wm**2
	t.lwin(sob, tob)

	# compute downscaled shortwave (SWf) Wm**2
	t.swin(pob, sob,tob, stat,p.dtime)

	# compute downscaled precipitation in ms*1
	t.precip(sob,stat)

	t.wind(tob)


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
				"RH":t.r,
				"WS":t.wscor,
				"WD":t.wd, 
				"LWIN":t.LWf, 
				"SWIN":t.SWfglob, 
				"PRATE":t.TPf
				},index=p.dtime)
	df.index.name="datetime"

	fileout=wd + "/sim/" + simdir + "/forcing/meteo"+"c"+str(i+1)+".csv"
	df.to_csv(path_or_buf=fileout ,na_rep=-999,float_format='%.3f')
	logging.info(fileout + " complete")


