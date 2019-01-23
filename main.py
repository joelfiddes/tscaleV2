"""main

Main toposcale run script.


Example:

	python main.py



Attributes:


Todo:




"""

import era5 as e5
import tscale as ts
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
reload(e5)
reload(ts)

#=== Object description
#
#	* pob - pressure level object
#	* sob - surface object
#	* tob - toposcale object
#	* stat - station object

#=== Parameters (for INI) ==============================================
fp="/home/joel/src/tscaleV2/tests/plevel.nc"
fs="/home/joel/src/tscaleV2/tests/surface.nc"

""" that's it!  Now, you can create a Bunch
 whenever you want to group a few variables:
http://code.activestate.com/recipes/52308-the-simple-but-handy-collector-of-a-bunch-of-named/?in=user-97991 """
class Bunch:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

# station attribute structure
stat = Bunch(ele = 4000, slp = 0, asp = 180, svf = 1, long=9, lat=46, tz=1  )

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
