import era5 as e5
import tscale as ts
reload(e5)
reload(ts)

#=== Parameters (for INI) ==============================================
fp="/home/joel/mnt/myserver/sim/wfj_interim2/eraDat/PLEVEL.nc"
fs="/home/joel/mnt/myserver/sim/wfj_interim2/eraDat/SURF.nc"
mylon=9
mylat=46
statEle = 2000

#=== Pressure level object =============================================

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

#=== Surface object ====================================================

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

#=== toposcale object ====================================================

t = ts.tscale(p.z, statEle)

for v in p.varnames:
	if (v=='z'):
		continue
	print v
	dat = getattr(p, v) # allows dynamic call to instance variable

	t.tscale1D(dat)
	t.addVar(v,t.interpVar)

sob = s
tob= t

# compute downscaled longwave (LWf)
t.lwin(sob, tob)

# compute downscaled shortavae (SWf)
t.swin(sob)
