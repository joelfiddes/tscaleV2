"""main

Main toposcale run script.

Object description:

	* pob - pressure level object
	* sob - surface object
	* tob - toposcale object
	* stat - station object (points, cluster centroids or grid centroids)

Example:

	python main.py

Args:
	inDir: directory containing input meteo PLEV.nc and SURF.nc
	listpointsLoc: location of listpoints.txt
	outDir: location where output meteo files are written
	memeber: which ensemble memeber is to be downscaled from era5 ensembles dataset [1-10]
		startTime: ISO format #"2016-08-01 18:00:00"
	endTime: ISO format #"2016-08-01 18:00:00"

Todo:
	merge with tscale_run.py, need to handle optional parameter "memebers" and extractCGC5d fun.
"""

import pandas as pd # only required to import listpoints nicely, should be able to remove
import era5 as e5
import tscale as ts
import helper as hp
import numpy as np
import sys
import logging

#=== ARGS==============================================
inDir=sys.argv[1]
listpointsLoc=sys.argv[2]
outDir = sys.argv[3]
member = int(sys.argv[4])
startTime = sys.argv[5]
endTime = sys.argv[6]

fp=inDir+"/PLEV.nc"
fs=inDir+"/SURF.nc"
lpfile=listpointsLoc + "/listpoints.txt"
# read in lispoints
lp=pd.read_csv(lpfile)

#===============================================================================
#	Logging
#===============================================================================
#logging.basicConfig(level=logging.DEBUG, filename=wd+"/sim/"+ simdir+"/logfile", filemode="a+",
                        #format="%(asctime)-15s %(levelname)-8s %(message)s")


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
	p.dtime = p.dtime[startIndex:endIndex]


	for v in p.varnames:

		#v=p.varnames[1]
		#p.getVar(v)
		#name = p.myvar.name

		p.extractCgc5d(v, member,startIndex, endIndex) # adds data to self
		
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

	# compute step in seconds for accumulated surface fields
	a=s.dtime[2]-s.dtime[1]
	step = a.seconds

	for v in s.varnames:

		#v=p.varnames[1]
		#print v
		s.getVar(v)
		#name = p.myvar.name

		s.extractCgc4d(v,member,startIndex, endIndex) # adds data to self

		s.addVar(v,s.var) # adds data again with correct name - redundancy

	""" rad conversions """
	s.instRad(step)
	""" precip conversions """
	s.tp2rate(step)
	""" dimensions of data """
	s.addShape()

	# cut time vector to required length
	s.dtime = s.dtime[startIndex:endIndex]
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

	fileout=outDir + "/meteo"+"c"+str(i+1)+".csv"
	df.to_csv(path_or_buf=fileout ,na_rep=-999,float_format='%.3f')
	print(fileout + " complete")


