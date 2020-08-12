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
	python /home/caduff/src/tscaleV2/toposcale/tscale_run.py /home/caduff/sim/paiku/forcing/ /home/joel/sim/paiku/sim/g2 /home/caduff/sim/paiku/sim/g2/forcing/ 1979-09-01 2019-09-01 FALSE era5 FALSE

	
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
import tscale as ts
import helper as hp
import numpy as np
import sys
import logging
import netCDF4 as nc
from tqdm import tqdm
#import ipdb
#ipdb.set_trace()
#=== ARGS==============================================
inDir=sys.argv[1] # /home/joel/sim/imis/forcing/
home=sys.argv[2] # /home/joel/sim/imis/
outDir = sys.argv[3] # 
startTime = sys.argv[4]
endTime = sys.argv[5]
windCor=sys.argv[6]
dataset=sys.argv[7]
plapse=sys.argv[8]
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


if dataset=='interim':
	import erai as era

if dataset=='era5':
	import era5 as era
#===============================================================================
#	Logging
#===============================================================================
logging.basicConfig(level=logging.DEBUG, filename=home+"/tscale_logfile", filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")

# add if statement to accomodate single line listpoints with eg id=5
#if lp.id.size  == 1:
#	stat = hp.Bunch(ele = lp.ele[0], slp = lp.slp[0],asp = lp.asp[0],svf = lp.svf[0],lon = lp.lon[0], lat =lp.lat[0],sro = lp.surfRough[0],tz = lp.tz[0]  )

# case of multipoint lispoints

for i in tqdm(range(lp.id.size)):

	# station attribute structure , tz always =0 for case of ERA5
	stat = hp.Bunch(ele = lp.ele[i], slp = lp.slp[i],asp = lp.asp[i],svf = lp.svf[i],lon = lp.lon[i], lat =lp.lat[i],sro = lp.surfRough[i],tz = lp.tz[i]  )

	#=== Init Pressure level object =============================================

	""" preprocess pressure level fields and add to structure """
	pob=era.Plev(fp, stat.lat, stat.lon)
	pob.getVarNames()

	""" Datetime structure """
	pob.addTime()
	startIndexP = int(np.where(pob.dtime==startTime)[0])#"2016-08-01 18:00:00"
	endIndexP = int(np.where(pob.dtime==endTime)[0])#"2016-08-01 18:00:00"
	pob.dtime=pob.dtime[startIndexP:endIndexP]	

	""" Add plevels """
	pob.plevels()
	nlev=pob.levels.shape[0]


	#=== Init Surface level object =============================================
	""" preprocess surface level fields and add to structure """
	s=era.Surf(fs, stat.lat, stat.lon)
	s.getVarNames()

	""" Datetime structure """
	s.addTime()
	startIndexS = int(np.where(s.dtime==startTime)[0])#"2016-08-01 18:00:00"
	endIndexS = int(np.where(s.dtime==endTime)[0])#"2016-08-01 18:00:00"
	s.dtime=s.dtime[startIndexS:endIndexS]

	# check for unequal dtime ie as case for interim
	# extract time vectors to use in interpolation and convert from 
	# hour since to sec since by *3600
	if pob.dtime.shape!=s.dtime.shape:

		f = nc.Dataset(fs)   
		time_out = f.variables['time'][:].astype(np.int64) *3600
		# cut to correct dates
		time_out=time_out[startIndexS:endIndexS]

		f = nc.Dataset(fp)   
		time_in = f.variables['time'][:].astype(np.int64)  *3600
		# cut to correct dates
		time_in=time_in[startIndexP:endIndexP]





	# only required for erai, must have buffered datasets (ie a day before required startDate in raw data) as download will start at 03:00 on given day. 
	# we assume day starts at 00:00 so s.dtime==startTime will not match in this case
	#if p.dtime.shape!=s.dtime.shape:

#	"""cut p.dtime to start and end"""
#	s.dtime=s.dtime[startIndexS:endIndexS]

#	# we always trim p.dtime here
#	startIndexP = int(np.where(p.dtime==startTime)[0])#"2016-08-01 18:00:00"
#	endIndexP = int(np.where(p.dtime==endTime)[0])#"2016-08-01 18:00:00"

#	"""cut p.dtime to start and end"""
#	p.dtime=p.dtime[startIndexP:endIndexP]	

	for v in pob.varnames:

		#v=pob.varnames[1]
		#pob.getVar(v)
		#name = pob.myvar.name

		pob.extractCgc(v, startIndexP, endIndexP) # adds data to self
		
		# do time interpolation of PLEV in case of ERAI
		if pob.dtime.shape!=s.dtime.shape:
			logging.info("PLEV timestep != SURF timestep")
			# init dummy container
			pob.var2 = np.zeros((s.dtime.shape[0],pob.var.shape[1]))
			for lev in range(0,nlev):
				pob.var2[:,lev]= era.series_interpolate(time_out, time_in, pob.var[:,lev])

			pob.var = pob.var2

		pob.addVar(v,pob.var) # adds data again with correct name - redundancy
		#can now remove pob.var it is redundant

	
	""" dimensions of data """
	pob.addShape()

	logging.info("made a POB!")
	#=== Surface object ====================================================


	# compute step in seconds for accumulated surface fields
	a=s.dtime[2]-s.dtime[1]
	step = a.seconds

	for v in s.varnames:

		#v=p.varnames[1]
		#print v
		s.getVar(v)
		#name = p.myvar.name

		s.extractCgc(v ,startIndexS, endIndexS) # adds data to self

		s.addVar(v,s.var) # adds data again with correct name - redundancy

	''' convert accumulated field swin,lwin and P '''

	if dataset=='interim':
		s.ssrd = s.cummulative2total(s.ssrd, s.dtime)
		s.strd = s.cummulative2total(s.strd, s.dtime)
		""" rad interim always accumulated over step """
		s.instRad(step)

	if dataset=='era5':
		""" rad era5 always accumulated over 1h or 3600s """
		s.instRad(3600)


	if dataset=='interim':
		s.tp = s.cummulative2total(s.tp, s.dtime)
		""" precip conversions m/timestep to mm/h (PRATE) and mm/step (PSUM)"""
		s.tp2rate(step)

	if dataset=='era5':
		""" p  era5 always accumulated over 1h or 3600s """
		s.tp2rate(3600)

	""" dimensions of data """
	s.addShape()
	""" Datetime structure """
	#s.addTime()
	""" surface elevation of coarse grid"""
	s.gridEle()
	sob = s
	logging.info("made a SOB!")
	#=== toposcale object ====================================================

	"""init object"""
	t = ts.tscale(pob.z)

	""" Downscale pl fields """
	for v in pob.varnames:
		if (v=='z'):
			continue
		#print v
		dat = getattr(pob, v) # allows dynamic call to instance variable

		t.tscale1D(dat,stat)
		t.addVar(v,t.interpVar)
	

	# constrain r to interval 5-100 - do here as required by LWin parameterisation
	t.r[t.r <5]=5
	t.r[t.r>100]=100

	tob = t
	# compute downscaled longwave (LWf) Wm**2
	t.lwin(sob, tob, stat)

	# compute downscaled shortwave (SWf) Wm**2
	t.swin(pob, sob,tob, stat,s.dtime)

	# compute downscaled precipitation in mm/h (PRATE) and m/step (PSUM) - we dont apply scaling now
	t.precip(sob,stat,plapse)

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
	
	# partition prate to rain snow (mm/hr) for FSM
	lowthresh=272.15
	highthresh = 274.15
	d = {'prate': t.prate, 'ta': t.t }
	df = pd.DataFrame(data=d)
	snow = df.prate.where(df.ta<lowthresh) 
	rain=df.prate.where(df.ta>=highthresh) 

	mix1S = df.prate.where((df.ta >= lowthresh) & (df.ta<=highthresh), inplace = False)
	mix1T = df.ta.where((df.ta >= lowthresh) & (df.ta<=highthresh), inplace = False)
	mixSno=(highthresh - mix1T) / (highthresh-lowthresh)
	mixRain=1-mixSno
	addSnow=mix1S*mixSno
	addRain=mix1S*mixRain
	
	# nas to 0
	snow[np.isnan(snow)] = 0 	
	rain[np.isnan(rain)] = 0 
	addRain[np.isnan(addRain)] = 0 
	addSnow[np.isnan(addSnow)] = 0 

	# linearly reduce snow to zero in steep slopes for FSM
	
	#if steepSnowReduce=="TRUE": # make this an option if need that in future
	snowSMIN=30.
	snowSMAX=80.
	slope=stat.slp

	k= (snowSMAX-slope)/(snowSMAX - snowSMIN)

	if slope<snowSMIN:
		k=1
	if slope>snowSMAX:
		k=0

	t.snowTot=(snow+addSnow) * k
	t.rainTot=rain + addRain


	logging.info("made a TOB!")

	df = pd.DataFrame({	"TA":t.t, 
				"RH":t.r*0.01, #meteoio 0-1
				"VW":t.ws, 
				"DW":t.wd, 
				"ILWR":t.LWf, 
				"ISWR":t.SWfglob, 
				"PINT":t.prate,
				"PSUM":t.psum,
				"P":t.psf,
				"Rf":np.array(t.rainTot),
				"Sf":np.array(t.snowTot)
				},index=s.dtime)
	df.index.name="datetime"

	# fill outstanding nan in SW routine with 0 (night)
	df.ISWR = df.ISWR.fillna(0)
	fileout=outDir+"/meteo"+"c"+str(i+1)+".csv"
	df.to_csv(path_or_buf=fileout ,na_rep=-999,float_format='%.3f')
	logging.info(fileout + " complete")


