import era5 as e5
import tscale as ts
reload(e5)
reload(ts)

#=== Object description
#
#	* pob - pressure level object
#	* sob - surface object
#	* tob - toposcale object
#	* stat - station object

#=== Parameters (for INI) ==============================================
fp="/home/joel/mnt/myserver/sim/wfj_interim2/eraDat/PLEVEL.nc"
fs="/home/joel/mnt/myserver/sim/wfj_interim2/eraDat/SURF.nc"

# add to stat?
mylon=9
mylat=46

""" that's it!  Now, you can create a Bunch
 whenever you want to group a few variables:
http://code.activestate.com/recipes/52308-the-simple-but-handy-collector-of-a-bunch-of-named/?in=user-97991 """
class Bunch:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

# station attribute structure
stat = Bunch(ele = 2000, slp = 30, asp = 180, svf = 0.9  )

#=== Pressure level object =============================================

""" preprocess pressure level fields and add to structure """
p=e5.Plev(fp, mylat, mylon)
p.getVarNames()

for v in p.varnames:

	#v=p.varnames[1]
	print v
	#p.getVar(v)
	#name = p.myvar.name
	p.extractCgc(v) # adds data to self
	p.addVar(v,p.var) # adds data again with correct name - redundancy
	#can now remove p.var it is redundant

# conversions
#p.convZ()
p.addShape()
p.addTime()
p.plevels()
pob = p
#=== Surface object ====================================================

""" preprocess surface level fields and add to structure """
s=e5.Surf(fs, mylat, mylon)
s.getVarNames()

for v in s.varnames:

	#v=p.varnames[1]
	print v
	s.getVar(v)
	#name = p.myvar.name
	s.extractCgc(v) # adds data to self
	s.addVar(v,s.var) # adds data again with correct name - redundancy

# conversions
s.instRad()
s.addShape()
s.addTime()
sob = s
#=== toposcale object ====================================================

""" Downscale fields """
t = ts.tscale(p.z)

for v in p.varnames:
	if (v=='z'):
		continue
	print v
	dat = getattr(p, v) # allows dynamic call to instance variable

	t.tscale1D(dat,stat)
	t.addVar(v,t.interpVar)
tob = t


# compute downscaled longwave (LWf)
t.lwin(sob, tob)

# compute downscaled shortavae (SWf)
t.swin(pob, sob,tob, stat)
