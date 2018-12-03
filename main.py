
#=== Parameters (for INI) ==============================================================
fp="/home/joel/src/tscaleV2/PLEV.nc"
mylon=9
mylat=46


#=== code ==============================================================
import era5 as e5
p=e5.Plev(fp, mylat, mylon)
p.getVarNames()

for v in p.vars:

	#v=p.vars[1]
	print v
	p.getVar(v)
	#name = p.myvar.name
	p.extractCgc(v) # adds data to self
	p.addVar(v,p.var) # adds data again with correct name - redundancy
	#can now remove p.var it is redundant

fp="/home/joel/src/tscaleV2/SURF.nc"

s=e5.Surf(fp, mylat, mylon)
s.getVarNames()

for v in s.vars:

	#v=p.vars[1]
	print v
	s.getVar(v)
	#name = p.myvar.name
	s.extractCgc(v) # adds data to self
	s.addVar(v,s.var) # adds data again with correct name - redundancy