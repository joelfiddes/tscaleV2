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

import sys
import os
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from configobj import ConfigObj
import logging
import era5 as e5
import tscale as ts
import helper as hp
import fetch_era5 as fe
#reload(helper)
#====================================================================
#	Config setup
#====================================================================
#os.system("python writeConfig.py") # update config DONE IN run.sh file
config = ConfigObj(sys.argv[1])
#os.remove("logfile")
logging.basicConfig(level=logging.DEBUG, filename="logfile", filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")

logging.info("Start run")
#====================================================================
#	Define paths
#====================================================================
wd= config["main"]["wd"]
#shpPath =  wd + "/spatial/" + config["main"]["shp"]
eraDir = wd + "/forcing/"
if not os.path.exists(eraDir):
	os.makedirs(eraDir)
#====================================================================
#	Get domain extent from shape
#====================================================================
#shp_extent = hp.get_shp_extent(shpPath)
#lonW = shp_extent[0]
#lonE = shp_extent[1]
#latS = shp_extent[2]
#latN = shp_extent[3]

lonW = config["main"]["lonW"]
lonE = config["main"]["lonE"]
latN = config["main"]["latN"]
latS = config["main"]["latS"]

#====================================================================
#	Fetch forcing
#====================================================================
fe.retrieve_era5_surf(config,eraDir, latN,latS,lonE,lonW)
fe.retrieve_era5_plev(config,eraDir, latN,latS,lonE,lonW)

# concat monthly files
eraCat(eraDir, "SURF")
eraCat(eraDir, "PLEV")



#=== Parameters (for INI) DEFAULT LOCATION==============================================
#fp="/home/joel/src/tscaleV2/tests/plevel.nc"
#fs="/home/joel/src/tscaleV2/tests/surface.nc"

# station attribute structure
stat = hp.Bunch(ele = 4000, slp = 0, asp = 180, svf = 1, long=9, lat=46, tz=1  )

#=== Pressure level object =============================================

""" preprocess pressure level fields and add to structure """
p=e5.Plev(fp, stat.lat, stat.long)
p.getVarNames()

for v in p.varnames:

	#v=p.varnames[1]
	print v
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
s=e5.Surf(fs, stat.lat, stat.long)
s.getVarNames()

for v in s.varnames:

	#v=p.varnames[1]
	print v
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
	print v
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
