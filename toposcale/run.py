"""main

Main toposcale run script.

Object description:

	* pob - pressure level object
	* sob - surface object
	* tob - toposcale object
	* stat - station object (points, cluster centroids or grid centroids)

Example:

	python main.py

Attributes:

Todo:

"""




import pandas as pd
import era5 as e5
import tscale as ts
import helper as hp
import numpy as np
#=== ARGS==============================================
fp="/home/joel/src/tscaleV2/tests/plevel.nc" #ARG3
fs="/home/joel/src/tscaleV2/tests/surface.nc"# ARG2
lpfile= "/home/joel/src/tscaleV2/tests/listpoints.txt" # ARG1

# read in lispoints
lp=pd.read_csv(lpfile)


for i in range(lp.pk.size):

	# station attribute structure
	stat = hp.Bunch(ele = lp.ele[i], slp = lp.slp[i],asp = lp.asp[i],svf = lp.svf[i],lon = lp.lon[i], lat =lp.lat[i],sro = lp.surfRough[i],tz = lp.tz[i]  )

	#=== Pressure level object =============================================

	""" preprocess pressure level fields and add to structure """
	p=e5.Plev(fp, stat.lat, stat.lon)
	p.getVarNames()

	for v in p.varnames:

		#v=p.varnames[1]
		#p.getVar(v)
		#name = p.myvar.name
		p.extractCgc(v) # adds data to self
		p.addVar(v,p.var) # adds data again with correct name - redundancy
		#can now remove p.var it is redundant

	""" dimensions of data """
	p.addShape()

	""" Datetime structure """
	p.addTime()

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
		s.extractCgc(v) # adds data to self
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
	

	df = pd.DataFrame({	"T":t.t, 
				"RH":t.r,
				"WS":t.wscor,
				"WD":t.wd, 
				"LWIN":t.LWf, 
				"SWIN":t.SWfglob, 
				"TP":t.TPf
				},index=p.dtime)
	df.index.name="datetime"

	fileout="meteo"+str(i)+".csv"
	df.to_csv(path_or_buf=fileout ,na_rep=-999,float_format='%.3f')
	print(fileout + " complete")


