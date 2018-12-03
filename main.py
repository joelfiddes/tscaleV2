
#=== Parameters (for INI) ==============================================================
fp="/home/joel/src/tscaleV2/PLEV.nc"
mylon=9
mylat=46


#=== code ==============================================================
import era5 as e5
#eraCat("/home/joel/mnt/myserver/sim/wfj_era5/eraDat/", "PLEVEL")

p=e5.Plev(fp)
p.vars()

#for v in myvars:
	v=myvars[1]
	print v
	p.getVar(v)
	name = p.myvar.name
	b = p.extractCgc(v, lon=mylon, lat=mylat)


#testing only

p=Plev(fp)
p.getVarNames()


v=p.vars[1]
print v
p.getVar(v)
#name = p.myvar.name
p.extractCgc(v, lon=mylon, lat=mylat) # adds data to self
p.addVar(v,p.var) # adds data again with correct name - redundancy
