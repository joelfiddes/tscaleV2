"""
3D toposcale of surface and pressurelevel fields including surface effects. Can be used to produce grids or point timeseries
at specific locations with lon,lat. Not used by TopoSUB (1D interp).

Non-linear increase in unit time per step. time/step increase with more steps, Why?
y=c(0.55,5.85,13.25,35.01, 87.94)
x=c(216,2160,4368,8784, 17544)
plot(x,y, ylab="T in mins", xlab="timesteps")


calls: tscale2d_era5_src.py and tscale3d_era5_src.py to perform interpolations.
        
Args:
	mode: 'grid' or 'points'

            
Returns:
	Timeseries of points or grids depending on mode
                
Example: 
	python tscale3D.py '/home/joel/sim/cci_perm/cci_test/' 'grid' '2014-09-01' '2014-09-10' "HRES" None

Notes:
- pressure level arrays in EDA (ascending) and inverse to HRES (descending)

example EDA geopotential
In [67]: z_interp[:,1]
Out[67]: 
array([ 1574.43959939,  3683.61224602,  5847.31984393,  8068.19587794,
       10343.13090714, 12663.88264812, 15037.5650956 , 17468.52776509,
       19963.81837095, 22529.31287884, 25170.64426949, 30685.20454384,
       36544.61456106, 42786.12626777, 49468.28618534, 56662.69980236,
       64459.43907668, 72955.92320608, 82299.63734099, 92752.36549403])

example HRES geopotential
In [25]: z_interp[:,1]
Out[25]: 
array([56332.15453   , 49194.88561242, 42571.16299455, 36401.05284972,
       30615.17900161, 25151.02196225, 22524.83301162, 19964.00851129,
       17465.31275158, 15028.14684696, 12649.5725558 , 10332.79436759,
        8070.00635785,  5862.73901075,  3709.39604655,  1604.20734139])

Todo:
	- paralelise by time - just split into 2 year chuncks eg 20 cores for 40 years
	- bring subroutines into single file eg *src.py
"""
import os
import sys
#sys.path.insert(0, "/home/joel/src/tscaleV2/toposcale/")

import pandas as pd
import netCDF4 as nc
import numpy as np
import logging
from datetime import datetime
import helper as hp
import solarGeom as sg
import era5
import tscale2d_era5_src as t2d
import tscale3d_era5_src as t3d
import time
import xarray as xr
#reload(t2d)
#reload(t3d)

#args
#wdir="/home/joel/results"
#mode="points" # "grid"

#test='False'
#start="2014-09-01"
#end="2016-09-01"

# wdir='/home/joel/sim/imis/'
# mode= 'point'
# start= '2004-09-01'
# end='2005-09-01'
# dataset="HRES"
# member=None

def main(wdir, mode, start, end, dataset, member=None):


	#Set up Log
	# make out path for results
	logs = wdir+"/logs/"
	if not os.path.exists(logs):
		os.makedirs(logs)

	logfile=logs+"/logfile"+start
	#if os.path.isfile(logfile) == True:
		#os.remove(logfile)
	logging.basicConfig(level=logging.DEBUG, filename=logfile, filemode="a+",format="%(asctime)-15s %(levelname)-8s %(message)s")
	logging.info("Running member"+ str(member))
	# start timer
	start_time = time.time()


	# make these names standard:
	stationsfile= wdir+'/listpoints.txt'
	demfile = wdir+'/predictors/ele.nc'
	surfile=wdir+'/forcing/SURF.nc'
	plevfile=wdir+ '/forcing/PLEV.nc'
	
	# make out path for results
	out = wdir+"/out/"
	if not os.path.exists(out):
		os.makedirs(out)

	# convert to python indexing and int
	if dataset=='EDA':
		member =int(member)-1

	if mode=="grid":
		dem  = nc.Dataset(demfile)
		dem_ele = dem.variables['Band1'][:]

	
	# this is main dtime used in code
	f = nc.Dataset( plevfile)
	nctime = f.variables['time']
	dtime = pd.to_datetime(nc.num2date(nctime[:],nctime.units, calendar="standard"))
	starti = np.asscalar(np.where(dtime==start+' 00:00:00')[0])
	endi = np.asscalar(np.where(dtime==end+' 00:00:00')[0])
	dtime= dtime[starti:endi,]
	timesteps=len(dtime)

	# this dtime is only to check for different timesteps
	f = nc.Dataset( surfile)
	nctime = f.variables['time']
	dtime2 = pd.to_datetime(nc.num2date(nctime[:],nctime.units, calendar="standard"))
	starti2 = np.asscalar(np.where(dtime2==start+' 00:00:00')[0])
	endi2 = np.asscalar(np.where(dtime2==end+' 00:00:00')[0])
	dtime2= dtime2[starti2:endi2,]
	timesteps2=len(dtime2)

	# catch cases of 3h plev and 1h surf
	if len(dtime)!=len(dtime2):
		print("resampling 1h surf data to 3h")
		df =xr.open_dataset(surfile)
		df2 =df.resample(time='3H').mean()
		df2['tp']=df.tp*3
		df2.to_netcdf(surfile)

	# constants
	g=9.81

	logging.info("Running "+ str(timesteps)+ " timesteps")
	
	#===============================================================================
	# Point stuff
	#
	#
	#
	#
	#===============================================================================
	if mode=="points" or mode=="point":
		if os.path.isfile(stationsfile) == False:
			logging.info("No listpoints.txt found!")

		logging.info("Start points run...")
		
		#open points
		mystations=pd.read_csv(stationsfile)
	#===============================================================================
	# make a pob hack
	#===============================================================================

		import recapp_era5 as rc
		f = nc.Dataset(plevfile)
		lev = f.variables['level'][:]
		 #dataImport.plf_get()    # pressure level air temperature
		var='t'

		# station case
		if dataset=="HRES":
			ds = rc.t3d( pl=plevfile)
		if dataset=="EDA":
			ds = rc.t3d_eda( pl=plevfile)

		#timesteps = ds.pl.variables['time'][:].size
		out_xyz_dem, lats, lons, shape, names= ds.demGrid(stations=mystations)

		# init grid stack
		# init grid stack
		xdim=lev.shape[0]

		ydim=shape[0]
		t_interp_out = np.zeros((xdim,ydim))
		z_interp_out = np.zeros((xdim,ydim))
			

		for timestep in range(starti,endi):
			#logging.info(str(round(float(timestep)/float(timesteps)*100,0))+ "% done")

			gridT,gridZ,gridLat,gridLon=ds.gridValue(var,timestep)
			if dataset=="HRES":
				t_interp, z_interp = ds.inLevelInterp(gridT,gridZ, gridLat,gridLon,out_xyz_dem)
			if dataset=="EDA":
				t_interp, z_interp = ds.inLevelInterp(gridT,gridZ, gridLat,gridLon,out_xyz_dem,member)

			t_interp_out = np.dstack((t_interp_out, t_interp))
			z_interp_out = np.dstack((z_interp_out, z_interp))
		# drop init blank layer
		tinterp =t_interp_out[:,:,1:]
		zinterp =z_interp_out[:,:,1:]

		pob= hp.Bunch(t=tinterp,z=zinterp, levels=lev)

		logging.info("made a pob!")
		#===============================================================================
		# tscale3d (results and input to radiation routines)
		#===============================================================================
		if dataset=="HRES":
			t = t3d.main(wdir, 'point', 't', starti,endi, dataset )
			r = t3d.main(wdir, 'point', 'r', starti,endi,dataset )
			u = t3d.main(wdir, 'point', 'u', starti,endi,dataset )
			v = t3d.main(wdir, 'point', 'v', starti,endi,dataset)
		if dataset=="EDA":
			t = t3d.main(wdir, 'point', 't', starti,endi, dataset, member)
			r = t3d.main(wdir, 'point', 'r', starti,endi,dataset, member)
			u = t3d.main(wdir, 'point', 'u', starti,endi,dataset, member)
			v = t3d.main(wdir, 'point', 'v', starti,endi,dataset, member)
		# compute wind speed and direction
		ws = np.sqrt(u**2+v**2)
		wd =   (180 / np.pi) * np.arctan(u/v) + np.where(v>0,180,np.where(u>0,360,0))

		tob = hp.Bunch(t=t,r=r,u=u, v=v, ws=ws,wd=wd, dtime=dtime)
		logging.info("made a tob!")
		#===============================================================================
		# tscale2d
		#===============================================================================
		if dataset=="HRES":
			t2m = t2d.main(wdir, 'point', 't2m', starti,endi,dataset)
			tp = t2d.main(wdir, 'point', 'tp', starti,endi,dataset)
			ssrd = t2d.main(wdir, 'point', 'ssrd', starti,endi,dataset)
			strd = t2d.main(wdir, 'point', 'strd', starti,endi,dataset)
			tisr = t2d.main(wdir, 'point', 'tisr', starti,endi,dataset)
			d2m = t2d.main(wdir, 'point', 'd2m', starti,endi,dataset)
			z = t2d.main(wdir, 'point', 'z', starti,endi,dataset) # always true as this is time invariant
		if dataset=="EDA":
			t2m = t2d.main(wdir, 'point', 't2m', starti,endi,dataset,member)
			tp = t2d.main(wdir, 'point', 'tp', starti,endi,dataset,member)
			ssrd = t2d.main(wdir, 'point', 'ssrd', starti,endi,dataset,member)
			strd = t2d.main(wdir, 'point', 'strd', starti,endi,dataset,member)
			tisr = t2d.main(wdir, 'point', 'tisr', starti,endi,dataset,member)
			d2m = t2d.main(wdir, 'point', 'd2m', starti,endi,dataset,member)
			z = t2d.main(wdir, 'point', 'z', starti,endi,dataset,member) # always true as this is time invariant

		gridEle=z[0,:]/g

		sob = hp.Bunch(t2m=t2m, tp=tp, ssrd=ssrd, strd=strd, tisr=tisr, d2m=d2m, z=z, gridEle=gridEle, dtime=dtime)
		logging.info("made a sob!")
		# functions

		def tp2rate(tp, step):
			""" convert tp from m/timestep (total accumulation over timestep) to rate in mm/h 

					Args:
						step: timstep in seconds (era5=3600, ensemble=10800)

					Note: both EDA (ensemble 3h) and HRES (1h) are accumulated over the timestep
					and therefore treated here the same.
					https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
			"""
			tp = tp/step*60*60 # convert metres per timestep -> m/hour 
			pmmhr = tp	*1000 # m/hour-> mm/hour
			return pmmhr

		def precipPoint(fineEle, sob):
			'''
			Args:
			fineEle:ele vector from station dataframe
			sob: contains gridEle, tp and dtime
			'''
			# convert TP to mm/hr
			
			lookups = {
				   1:0.35,
				   2:0.35,
				   3:0.35,
				   4:0.3,
				   5:0.25,
				   6:0.2,
				   7:0.2,
				   8:0.2,
				   9:0.2,
				   10:0.25,
				   11:0.3,
				   12:0.35
			}

			# Precipitation lapse rate, varies by month (Liston and Elder, 2006).
			pfis = sob.dtime.month.map(lookups)
			pfis = np.repeat(pfis.values[:,None], fineEle.size, axis=1)

			dz=(fineEle-sob.gridEle)/1e3  # Elevation difference in kilometers between the fine and coarse surface.
			
				   
			lp=(1+pfis.T*dz[:,None])/(1-pfis.T*dz[:,None])# Precipitation correction factor.
			#Pf=sob.pmmhr.T*lp
			prate=sob.pmmhr.T*lp # mm/hour
			psum=sob.tp.T*1000*lp # mm/timestep
			return prate, psum






		#===============================================================================
		# Precip
		#===============================================================================
		a=dtime[2]-dtime[1]
		step = a.seconds
		pmmhr = tp2rate(tp,step)
		sob.pmmhr = pmmhr
		tob.prate , tob.psum = precipPoint(mystations.ele, sob)
		logging.info("made prate!")
		#===============================================================================
		# Longwave
		#===============================================================================
		def instRad(sob, step):
			""" Convert SWin from accumulated quantities in J/m2 to 
			instantaneous W/m2 see: 
			https://confluence.ecmwf.int/pages/viewpage.action?pageId=104241513

			Args:
				step: timstep in seconds (era5=3600, ensemble=10800)

			Note: both EDA (ensemble 3h) and HRES (1h) are accumulated over the timestep
			and therefore treated here the same.
			https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
					"""
			sob.strd = sob.strd/step  
			sob.ssrd = sob.ssrd/step
			sob.tisr = sob.tisr/step 
			
		def lwin(sob,tob):
			"""Convert to RH (should be in a function). Following MG Lawrence 
			DOI 10.1175/BAMS-86-2-225 """
			A1=7.625 
			B1=243.04 
			C1=610.94
			tc=sob.t2m-273.15
			tdc=sob.d2m-273.15
			tf=tob.t-273.15 # fout.T
			c=(A1*tc)/(B1+tc)
			RHc=100*np.exp((tdc*A1-tdc*c-B1*c)/(B1+tdc)) # Inverting eq. 8 in Lawrence.

			""" Calculate saturation vapor pressure at grid and "subgrid" [also
			through function] using the Magnus formula."""
			
			svpf=C1*np.exp(A1*tf/(B1+tf))
			svpc=C1*np.exp(A1*tc/(B1+tc))

			"""Calculate the vapor pressure at grid (c) and subgrid (f)."""
			vpf=tob.r*svpf/1e2 # RHf
			vpc=RHc*svpc/1e2

			"""Use the vapor pressure and temperature to calculate clear sky
			 # emssivity at grid and subgrid. [also function]
			Konzelmann et al. 1994
			Ta in kelvin

			 """
			x1=0.43 
			x2=5.7
			cef=0.23+x1*(vpf/tob.t)**(1/x2) #Pretty sure the use of Kelvin is correct.
			cec=0.23+x1*(vpc/sob.t2m)**(1/x2)

			"""Diagnose the all sky emissivity at grid."""
			sbc=5.67e-8
			aec=sob.strd/(sbc*sob.t2m**4)
			# need to constrain to 1 as original code?

			""" Calculate the "cloud" emissivity at grid, assume this is the same at
		 	subgrid."""
			deltae=aec-cec

			""" Use the former cloud emissivity to compute the all sky emissivity at 
			subgrid. """
			aef=cef+deltae
			LWf=aef*sbc*tob.t**4
			return(LWf)

		instRad(sob,3600)
		tob.lwin = lwin(sob,tob)
		logging.info("made lwin")
		#===============================================================================
		# Shortwave
		#===============================================================================
		#def swin2D(pob,sob,tob, stat, dates): # does not work yet (booleen indexing of 2d arraya fauiles as np.min returns single value when we need 161
		#	
		#	""" toposcale surface pressure using hypsometric equation - move to own 
		#	class """
		#	g=9.81
		#	R=287.05  # Gas constant for dry air.
		#	ztemp = pob.z
		#	Ttemp = pob.t
		#	statz = stat.ele*g
		#	dz=ztemp.transpose()-statz[None,:,None] # transpose not needed but now consistent with CGC surface pressure equations

		#	psf=[]
		#	# loop through timesteps
		#	for i in range(0,dz.shape[0]):
		#		
		#		# 	# find overlying layer
		#		thisp = dz[i,0,:]==np.min(dz[i,0,:][dz[i,0,:]>0])
		#		thispVec = thisp.T.reshape(thisp.size)
		#		TtempVec = Ttemp.reshape(16*161, 20)

		#		# booleen indexing
		#		T1=TtempVec[thispVec,i]


		#		T1=Ttemp[thisp.T,i]
		#		z1=ztemp[thisp.T,i]
		#		p1=pob.levels[thisp]*1e2 #Convert to Pa.
		#		Tbar=np.mean([T1, tob.t[i]],axis=0)
		#		""" Hypsometric equation."""
		#		psf.append(p1*np.exp((z1-statz)*(g/(Tbar*R))))

		#	psf=np.array(psf).squeeze()

		#	## Specific humidity routine.
		#	# mrvf=0.622.*vpf./(psf-vpf); #Mixing ratio for water vapor at subgrid.
		#	#  qf=mrvf./(1+mrvf); # Specific humidity at subgrid [kg/kg].
		#	# fout.q(:,:,n)=qf; 


		#	""" Maybe follow Dubayah's approach (as in Rittger and Girotto) instead
		#	for the shortwave downscaling, other than terrain effects. """

		#	""" Height of the "grid" (coarse scale)"""
		#	Zc=sob.z

		#	""" toa """
		#	SWtoa = sob.tisr 

		#	""" Downwelling shortwave flux of the "grid" using nearest neighbor."""
		#	SWc=sob.ssrd

		#	"""Calculate the clearness index."""
		#	kt=SWc/SWtoa

		#	#kt[is.na(kt)==T]<-0 # make sure 0/0 =0
		#	#kt[is.infinite(kt)==T]<-0 # make sure 0/0 =0
		#	kt[kt<0]=0
		#	kt[kt>1]=0.8 #upper limit of kt
		#	kt=kt


		#	"""
		#	Calculate the diffuse fraction following the regression of Ruiz-Arias 2010 
		#	
		#	"""
		#	kd=0.952-1.041*np.exp(-1*np.exp(2.3-4.702*kt))
		#	kd = kd

		#	""" Use this to calculate the downwelling diffuse and direct shortwave radiation at grid. """
		#	SWcdiff=kd*SWc
		#	SWcdir=(1-kd)*SWc
		#	SWcdiff=SWcdiff
		#	SWcdir=SWcdir

		#	""" Use the above with the sky-view fraction to calculate the 
		#	downwelling diffuse shortwave radiation at subgrid. """
		#	SWfdiff=stat.svf*SWcdiff
		#	SWfdiff.set_fill_value(0)
		#	SWfdiff = SWfdiff.filled()

		#	""" Direct shortwave routine, modified from Joel. 
		#	Get surface pressure at "grid" (coarse scale). Can remove this
		#	part once surface pressure field is downloaded, or just check
		#	for existance. """

		#	ztemp = pob.z
		#	Ttemp = pob.t
		#	dz=ztemp.transpose()-sob.z

		#	psc=[]
		#	for i in range(0,dz.shape[1]):
		#	
		#		#thisp.append(np.argmin(dz[:,i][dz[:,i]>0]))
		#		thisp = dz[:,i]==np.min(dz[:,i][dz[:,i]>0])
		#		z0 = sob.z[i]
		#		T0 = sob.t2m[i]
		#		T1=Ttemp[i,thisp]
		#		z1=ztemp[i,thisp]
		#		p1=pob.levels[thisp]*1e2 #Convert to Pa.
		#		Tbar=np.mean([T0, T1],axis=0)
		#		""" Hypsometric equation."""
		#		psc.append(p1*np.exp((z1-z0)*(g/(Tbar*R))))

		#	psc=np.array(psc).squeeze()

		#	#T1=Ttemp(thisp)
		#	#z1=ztemp(thisp)
		#	#p1=pob.levels(thisp)*1e2 #Convert to Pa.
		#	#Tbar=mean([T0 T1])
		#	
		#	"""compute julian dates"""
		#	jd= sg.to_jd(dates)

		#	"""
		#	Calculates a unit vector in the direction of the sun from the observer 
		#	position.
		#	"""
		#	sunv=sg.sunvector(jd=jd, latitude=stat.lat, longitude=stat.lon, timezone=stat.tz)

		#	"""
		#	Computes azimuth , zenith  and sun elevation 
		#	for each timestamp
		#	"""
		#	sp=sg.sunpos(sunv)
		#	sp=sp

		#	# Cosine of the zenith angle.
		#	sp.zen=sp.zen
		#	#sp.zen=sp.zen*(sp.zen>0) # Sun might be below the horizon.
		#	muz=np.cos(sp.zen) 
		#	muz = muz
		#	# NB! psc must be in Pa (NOT hPA!).
		#	#if np.max(psc<1.5e3): # Obviously not in Pa
		#		#psc=psc*1e2
		#   
		#	
		#	# Calculate the "broadband" absorption coefficient. Elevation correction
		#	# from Kris
		#	ka=(g*muz/(psc))*np.log(SWtoa/SWcdir)	
		#	ka.set_fill_value(0)
		#	ka = ka.filled()
		#	# Note this equation is obtained by inverting Beer's law, i.e. use
		#	#I_0=I_inf x exp[(-ka/mu) int_z0**inf rho dz]
		#	# Along with hydrostatic equation to convert to pressure coordinates then
		#	# solve for ka using p(z=inf)=0.
		#	
		#	
		#	# Now you can (finally) find the direct component at subgrid. 
		#	SWfdir=SWtoa*np.exp(-ka*psf/(g*muz))

		#	""" Then perform the terrain correction. [Corripio 2003 / rpackage insol port]."""

		#	"""compute mean horizon elevation - why negative hor.el possible??? """
		#	horel=(((np.arccos(np.sqrt(stat.svf))*180)/np.pi)*2)-stat.slp
		#	if horel < 0:
		#		horel = 0 
		#	meanhorel = horel

		#	"""
		#	normal vector - Calculates a unit vector normal to a surface defined by 
		#	slope inclination and slope orientation.
		#	"""
		#	nv = sg.normalvector(slope=stat.slp, aspect=stat.asp)

		#	"""
		#	Method 1: Computes the intensity according to the position of the sun (sunv) and 
		#	dotproduct normal vector to slope.
		#	From corripio r package
		#	"""
		#	dotprod=np.dot(sunv ,np.transpose(nv)) 
		#	dprod = dotprod.squeeze()
		#	dprod[dprod<0] = 0 #negative indicates selfshading
		#	dprod = dprod

		#	"""Method 2: Illumination angles. Dozier"""
		#	saz=sp.azi
		#	cosis=muz*np.cos(stat.slp)+np.sin(sp.zen)*np.sin(stat.slp)*np.cos(sp.azi-stat.asp)# cosine of illumination angle at subgrid.
		#	cosic=muz # cosine of illumination angle at grid (slope=0).

		#	"""
		#	SUN ELEVATION below hor.el set to 0 - binary mask
		#	"""
		#	selMask = sp.sel
		#	selMask[selMask<horel]=0
		#	selMask[selMask>0]=1
		#	selMask = selMask

		#	"""
		#	derive incident radiation on slope accounting for self shading and cast 
		#	shadow and solar geometry
		#	BOTH formulations seem to be broken
		#	"""
		#	#SWfdirCor=selMask*(cosis/cosic)*SWfdir
		#	SWfdirCor=selMask*dprod*SWfdir
		#   
		#	SWfglob =  SWfdiff+ SWfdirCor
		#	return SWfglob
		#	""" 
		#	Missing components
		#	- terrain reflection
		#	"""
		def swin1D(pob,sob,tob, stat, dates, index):
				# many arrays transposed
				""" toposcale surface pressure using hypsometric equation - move to own 
				class 

				index: index of station array (numeric)
				"""
				g= 9.81
				R=287.05  # Gas constant for dry air.
				tz=0 # ERA5 is always utc0, used to compute sunvector
				ztemp = pob.z[:,index,:].T
				Ttemp = pob.t[:,index,:].T
				statz = stat.ele[index]*g
				dz=ztemp.T-statz # transpose not needed but now consistent with CGC surface pressure equations

				psf=[]
				# loop through timesteps
				#for i in range(starti,endi):
				for i in range(0,dates.size):
					
					# 	# find overlying layer
					thisp = dz[:,i]==np.min(dz[:,i][dz[:,i]>0])

					# booleen indexing
					T1=Ttemp[i,thisp]
					z1=ztemp[i,thisp]
					p1=pob.levels[thisp]*1e2 #Convert to Pa.
					Tbar=np.mean([T1, tob.t[i, index]],axis=0)
					""" Hypsometric equation."""
					psf.append(p1*np.exp((z1-statz)*(g/(Tbar*R))))

				psf=np.array(psf).squeeze()

				## Specific humidity routine.
				# mrvf=0.622.*vpf./(psf-vpf); #Mixing ratio for water vapor at subgrid.
				#  qf=mrvf./(1+mrvf); # Specific humidity at subgrid [kg/kg].
				# fout.q(:,:,n)=qf; 


				""" Maybe follow Dubayah's approach (as in Rittger and Girotto) instead
				for the shortwave downscaling, other than terrain effects. """

				""" Height of the "grid" (coarse scale)"""
				Zc=sob.z[:,index] # should this be a single value or vector?

				""" toa """
				SWtoa = sob.tisr[:,index]

				""" Downwelling shortwave flux of the "grid" using nearest neighbor."""
				SWc=sob.ssrd[:,index]

				"""Calculate the clearness index."""
				kt=SWc/SWtoa

				#kt[is.na(kt)==T]<-0 # make sure 0/0 =0
				#kt[is.infinite(kt)==T]<-0 # make sure 0/0 =0
				kt[kt<0]=0
				kt[kt>1]=0.8 #upper limit of kt
				kt=kt


				"""
				Calculate the diffuse fraction following the regression of Ruiz-Arias 2010 
				
				"""
				kd=0.952-1.041*np.exp(-1*np.exp(2.3-4.702*kt))
				kd = kd

				""" Use this to calculate the downwelling diffuse and direct shortwave radiation at grid. """
				SWcdiff=kd*SWc
				SWcdir=(1-kd)*SWc
				SWcdiff=SWcdiff
				SWcdir=SWcdir

				""" Use the above with the sky-view fraction to calculate the 
				downwelling diffuse shortwave radiation at subgrid. """
				SWfdiff=stat.svf[index]*SWcdiff
				SWfdiff = np.nan_to_num(SWfdiff) # convert nans (night) to 0


				""" Direct shortwave routine, modified from Joel. 
				Get surface pressure at "grid" (coarse scale). Can remove this
				part once surface pressure field is downloaded, or just check
				for existance. """

				ztemp = pob.z[:,index,:].T
				Ttemp = pob.t[:,index,:].T
				dz=ztemp.transpose()-sob.z[:,index]

				psc=[]
				for i in range(0,dz.shape[1]):
				
					#thisp.append(np.argmin(dz[:,i][dz[:,i]>0]))
					thisp = dz[:,i]==np.min(dz[:,i][dz[:,i]>0])
					z0 = sob.z[i,index]
					T0 = sob.t2m[i,index]
					T1=Ttemp[i,thisp]
					z1=ztemp[i,thisp]
					p1=pob.levels[thisp]*1e2 #Convert to Pa.
					Tbar=np.mean([T0, T1],axis=0)
					""" Hypsometric equation."""
					psc.append(p1*np.exp((z1-z0)*(g/(Tbar*R))))

				psc=np.array(psc).squeeze()

				#T1=Ttemp(thisp)
				#z1=ztemp(thisp)
				#p1=pob.levels(thisp)*1e2 #Convert to Pa.
				#Tbar=mean([T0 T1])
				
				"""compute julian dates"""
				jd= sg.to_jd(dates)

				"""
				Calculates a unit vector in the direction of the sun from the observer 
				position.
				"""
				sunv=sg.sunvector(jd=jd, latitude=stat.lat[index], longitude=stat.lon[index], timezone=tz)

				"""
				Computes azimuth , zenith  and sun elevation 
				for each timestamp
				"""
				sp=sg.sunpos(sunv)
				sp=sp

				# Cosine of the zenith angle.
				sp.zen=sp.zen
				#sp.zen=sp.zen*(sp.zen>0) # Sun might be below the horizon.
				muz=np.cos(sp.zen) 
				muz = muz
				# NB! psc must be in Pa (NOT hPA!).
				#if np.max(psc<1.5e3): # Obviously not in Pa
					#psc=psc*1e2
			   
				
				# Calculate the "broadband" absorption coefficient. Elevation correction
				# from Kris
				ka=(g*muz/(psc))*np.log(SWtoa/SWcdir)	
				#ka.set_fill_value(0)
				#ka = ka.filled()
				# set inf (from SWtoa/SWcdir at nigh, zero division) to 0 (night)
				ka[ka == -np.inf] = 0
				ka[ka == np.inf] = 0

				# Note this equation is obtained by inverting Beer's law, i.e. use
				#I_0=I_inf x exp[(-ka/mu) int_z0**inf rho dz]
				# Along with hydrostatic equation to convert to pressure coordinates then
				# solve for ka using p(z=inf)=0.
				
				
				# Now you can (finally) find the direct component at subgrid. 
				SWfdir=SWtoa*np.exp(-ka*psf/(g*muz))

				""" Then perform the terrain correction. [Corripio 2003 / rpackage insol port]."""

				"""compute mean horizon elevation - why negative hor.el possible??? """
				horel=(((np.arccos(np.sqrt(stat.svf[index]))*180)/np.pi)*2)-stat.slp[index]
				if horel < 0:
					horel = 0 
				meanhorel = horel

				"""
				normal vector - Calculates a unit vector normal to a surface defined by 
				slope inclination and slope orientation.
				"""
				nv = sg.normalvector(slope=stat.slp[index], aspect=stat.asp[index])

				"""
				Method 1: Computes the intensity according to the position of the sun (sunv) and 
				dotproduct normal vector to slope.
				From corripio r package
				"""
				dotprod=np.dot(sunv ,np.transpose(nv)) 
				dprod = dotprod.squeeze()
				dprod[dprod<0] = 0 #negative indicates selfshading
				dprod = dprod

				"""Method 2: Illumination angles. Dozier"""
				saz=sp.azi
				cosis=muz*np.cos(stat.slp[index])+np.sin(sp.zen)*np.sin(stat.slp[index])*np.cos(sp.azi-stat.asp[index])# cosine of illumination angle at subgrid.
				cosic=muz # cosine of illumination angle at grid (slope=0).

				"""
				SUN ELEVATION below hor.el set to 0 - binary mask
				"""
				selMask = sp.sel
				selMask[selMask<horel]=0
				selMask[selMask>0]=1
				selMask = selMask

				"""
				derive incident radiation on slope accounting for self shading and cast 
				shadow and solar geometry
				BOTH formulations seem to be broken
				"""
				#SWfdirCor=selMask*(cosis/cosic)*SWfdir
				SWfdirCor=selMask*dprod*SWfdir
			   
				SWfglob =  SWfdiff+ SWfdirCor
				return SWfglob
				""" 
				Missing components
				- terrain reflection
				"""
				# init grid stack
		#init first row
		ntimestamps=dtime.shape[0]
		ts_swin = np.zeros((ntimestamps))
		for stationi in range(0, mystations.shape[0]):

			'''here we test for array full of NaNs due to points being 
			outside of grid. The NaN arrays are created in the 
			interpolation but only break the code here. If array is
 			all NaN then we just fill that stations slot 
			with NaN'''
			testNans = np.count_nonzero(~np.isnan(pob.z[:,stationi,:]))
			if testNans != 0:
				ts= swin1D(pob=pob,sob=sob,tob=tob, stat=mystations, dates=dtime, index=stationi)
				ts_swin=np.column_stack((ts_swin,ts))
			if testNans == 0:
				nan_vec = np.empty(ntimestamps) * np.nan
				ts_swin=np.column_stack((ts_swin,nan_vec))
			

		# drop init row
		tob.swin =ts_swin[:,1:]
		
		#===============================================================================
		# make dataframe (write individual files plus netcdf)
		#==============================================================================

		logging.info("Writing toposcale files...")
		for i in range(0,tob.t.shape[1]):
			df = pd.DataFrame({	"TA":tob.t[:,i], 
						"RH":tob.r[:,i]*0.01, #meteoio 0-1
						"VW":tob.ws[:,i],
						"DW":tob.wd[:,i],
						"ILWR":tob.lwin[:,i], 
						"ISWR":tob.swin[:,i], 
						"PINT":tob.prate[i,:],
						"PSUM":tob.psum[i,:]
						},index=tob.dtime)
			df.index.name="datetime"

			# fill outstanding nan in SW routine with 0 (night)
			df.ISWR = df.ISWR.fillna(0)

			if dataset=='EDA':
				fileout=wdir+"/out/meteo"+str(i)+"_"+start+"_"+str(member+1)+"_.csv" # convert member index back to 1-10

			if dataset=='HRES':
				fileout=wdir+"/out/meteo"+"c"+str(i+1)+".csv" # convert member index back to 1-10
			column_order = ['TA', 'RH', 'VW', 'DW', 'ILWR', 'ISWR', 'PINT', 'PSUM']
			df[column_order].to_csv(path_or_buf=fileout ,na_rep=-999,float_format='%.3f')
		#logging.info(fileout + " complete")

		#===============================================================================
		# Grid stuff
		#
		#===============================================================================
	if mode=="grid":
		logging.info("Running TopoSCALE3D grid")
		#===============================================================================
		# tscale3d
		#===============================================================================

		
		t = t3d.main( wdir, 'grid', 't', starti,endi,dataset)
		r = t3d.main( wdir, 'grid', 'r', starti,endi,dataset)
		gtob = hp.Bunch(t=t,r=r, dtime=dtime)

		#===============================================================================
		# tscale2d
		#===============================================================================
		t2m = t2d.main( wdir, 'grid', 't2m', starti,endi,dataset)
		tp = t2d.main( wdir, 'grid', 'tp', starti,endi,dataset)
		ssrd = t2d.main( wdir, 'grid', 'ssrd', starti,endi,dataset)
		strd = t2d.main( wdir, 'grid', 'strd', starti,endi,dataset)
		tisr = t2d.main( wdir, 'grid', 'tisr', starti,endi,dataset)
		d2m = t2d.main( wdir, 'grid', 'd2m', starti,endi,dataset)
		z = t2d.main( wdir, 'grid', 'z', starti,endi,dataset) # always true as this is time invariant
		gridEle=z[:,:,0]/g
		gsob = hp.Bunch(t2m=t2m, tp=tp, ssrd=ssrd, strd=strd, tisr=tisr, d2m=d2m, z=z, gridEle=gridEle, dtime=dtime)


		def precipGrid(fineEle, gsob):
			'''
			Args:
				fineEle
				gridEle:
			'''
			# convert TP to mm/hr
			
			lookups = {
				   1:0.3,
				   2:0.3,
				   3:0.3,
				   4:0.3,
				   5:0.25,
				   6:0.2,
				   7:0.2,
				   8:0.2,
				   9:0.2,
				   10:0.25,
				   11:0.3,
				   12:0.3
			}

			# Precipitation lapse rate, varies by month (Liston and Elder, 2006).
			pfis = gsob.dtime.month.map(lookups)
			

			dz=(fineEle-gsob.gridEle)/1e3  # Elevation difference in kilometers between the fine and coarse surface.
			dz2 =dz.reshape(dz.size) #make grid a vector
			pfis2 = np.repeat(pfis.values[:,None], dz2.size, axis=1)
				   
			lp=(1+pfis2.T*dz2[:,None])/(1-pfis2.T*dz2[:,None])# Precipitation correction factor.
			lp2 = lp.reshape(gsob.pmmhr.shape)
			Pf=gsob.pmmhr*lp2
			
			return Pf

		def tp2rate(tp, step):
			""" convert tp from m/timestep (total accumulation over timestep) to rate in mm/h 

					Args:
						step: timstep in seconds (era5=3600, ensemble=10800)

					Note: both EDA (ensemble 3h) and HRES (1h) are accumulated over the timestep
					and therefore treated here the same.
					https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
			"""
			tp = tp/step*60*60 # convert metres per timestep (in secs) -> m/hour 
			pmmhr = tp	*1000 # m/hour-> mm/hour
			return pmmhr

		a=dtime[2]-dtime[1]
		step = a.seconds
		gsob.pmmhr = tp2rate(tp,step)
		grid_prate = precipGrid(dem_ele,gsob)
		
		dem  = nc.Dataset(demfile)
		lon = dem.variables['lon'][:]
		lat = dem.variables['lat'][:]
		# These packages problem on cluster
		from osgeo import gdal
		from osgeo import gdal_array
		from osgeo import osr

		for i in range(0, grid_prate.shape[2]):

			myname=wdir+'/out/prate'+str(i)+'.tif'
			array = grid_prate[::-1,:,i]
			#lat = out_xyz_dem[:,0].reshape(l.shape)
			#lon = out_xyz_dem[:,1].reshape(l.shape)

			xmin,ymin,xmax,ymax = [lon.min(),lat.min(),lon.max(),lat.max()]
			nrows,ncols = np.shape(array)
			xres = (xmax-xmin)/float(ncols)
			yres = (ymax-ymin)/float(nrows)
			geotransform=(xmin,xres,0,ymax,0, -yres)   

			output_raster = gdal.GetDriverByName('GTiff').Create(myname,ncols, nrows, 1 ,gdal.GDT_Float32)# Open the file
			output_raster.GetRasterBand(1).WriteArray( array )  # Writes my array to the raster
			output_raster.SetGeoTransform(geotransform)# Specify its coordinates
			srs = osr.SpatialReference()# Establish its coordinate encoding
			srs.ImportFromEPSG(4326)   # This one specifies WGS84 lat long.
			output_raster.SetProjection(srs.ExportToWkt())# Exports the coordinate system 
			output_raster = None

		# for i in range(0, t.shape[2]):

		# 	myname=wdir+'/out/t'+str(i)+'.tif'
		# 	array = t[:,:,i]
		# 	#lat = out_xyz_dem[:,0].reshape(l.shape)
		# 	#lon = out_xyz_dem[:,1].reshape(l.shape)

		# 	xmin,ymin,xmax,ymax = [lon.min(),lat.min(),lon.max(),lat.max()]
		# 	nrows,ncols = np.shape(array)
		# 	xres = (xmax-xmin)/float(ncols)
		# 	yres = (ymax-ymin)/float(nrows)
		# 	geotransform=(xmin,xres,0,ymax,0, -yres)   

		# 	output_raster = gdal.GetDriverByName('GTiff').Create(myname,ncols, nrows, 1 ,gdal.GDT_Float32)# Open the file
		# 	output_raster.GetRasterBand(1).WriteArray( array )  # Writes my array to the raster
		# 	output_raster.SetGeoTransform(geotransform)# Specify its coordinates
		# 	srs = osr.SpatialReference()# Establish its coordinate encoding
		# 	srs.ImportFromEPSG(4326)   # This one specifies WGS84 lat long.
		# 	output_raster.SetProjection(srs.ExportToWkt())# Exports the coordinate system 
		# 	output_raster = None
		#===============================================================================
		# Longwave
		#===============================================================================
		def instRad(sob, step):
			""" Convert SWin from accumulated quantities in J/m2 to 
			instantaneous W/m2 see: 
			https://confluence.ecmwf.int/pages/viewpage.action?pageId=104241513

			Args:
				step: timstep in seconds (era5=3600, ensemble=10800)

			Note: both EDA (ensemble 3h) and HRES (1h) are accumulated over the timestep
			and therefore treated here the same.
			https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
					"""
			sob.strd = sob.strd/step  
			sob.ssrd = sob.ssrd/step
			sob.tisr = sob.tisr/step 
		
		def lwin(sob,tob):
			"""Convert to RH (should be in a function). Following MG Lawrence 
			DOI 10.1175/BAMS-86-2-225 """
			A1=7.625 
			B1=243.04 
			C1=610.94
			tc=sob.t2m-273.15
			tdc=sob.d2m-273.15
			tf=tob.t-273.15 # fout.T
			c=(A1*tc)/(B1+tc)
			RHc=100*np.exp((tdc*A1-tdc*c-B1*c)/(B1+tdc)) # Inverting eq. 8 in Lawrence.

			""" Calculate saturation vapor pressure at grid and "subgrid" [also
			through function] using the Magnus formula."""
			
			svpf=C1*np.exp(A1*tf/(B1+tf))
			svpc=C1*np.exp(A1*tc/(B1+tc))

			"""Calculate the vapor pressure at grid (c) and subgrid (f)."""
			vpf=tob.r*svpf/1e2 # RHf
			vpc=RHc*svpc/1e2

			"""Use the vapor pressure and temperature to calculate clear sky
			 # emssivity at grid and subgrid. [also function]
			Konzelmann et al. 1994
			Ta in kelvin

			 """
			x1=0.43 
			x2=5.7
			cef=0.23+x1*(vpf/tob.t)**(1/x2) #Pretty sure the use of Kelvin is correct.
			cec=0.23+x1*(vpc/sob.t2m)**(1/x2)

			"""Diagnose the all sky emissivity at grid."""
			sbc=5.67e-8
			aec=sob.strd/(sbc*sob.t2m**4)
			# need to constrain to 1 as original code?

			""" Calculate the "cloud" emissivity at grid, assume this is the same at
		 	subgrid."""
			deltae=aec-cec

			""" Use the former cloud emissivity to compute the all sky emissivity at 
			subgrid. """
			aef=cef+deltae
			LWf=aef*sbc*tob.t**4
			return(LWf)

		instRad(gsob,3600)
		ts_lwin = lwin(gsob,gtob)

		dem  = nc.Dataset(demfile)
		lon = dem.variables['lon'][:]
		lat = dem.variables['lat'][:]
		# These packages problem on cluster
		from osgeo import gdal
		from osgeo import gdal_array
		from osgeo import osr

		for i in range(0, t.shape[2]):

			myname=wdir+'/out/lwin'+str(i)+'.tif'
			array = t[:,:,i]
			#lat = out_xyz_dem[:,0].reshape(l.shape)
			#lon = out_xyz_dem[:,1].reshape(l.shape)

			xmin,ymin,xmax,ymax = [lon.min(),lat.min(),lon.max(),lat.max()]
			nrows,ncols = np.shape(array)
			xres = (xmax-xmin)/float(ncols)
			yres = (ymax-ymin)/float(nrows)
			geotransform=(xmin,xres,0,ymax,0, -yres)   

			output_raster = gdal.GetDriverByName('GTiff').Create(myname,ncols, nrows, 1 ,gdal.GDT_Float32)# Open the file
			output_raster.GetRasterBand(1).WriteArray( array )  # Writes my array to the raster
			output_raster.SetGeoTransform(geotransform)# Specify its coordinates
			srs = osr.SpatialReference()# Establish its coordinate encoding
			srs.ImportFromEPSG(4326)   # This one specifies WGS84 lat long.
			output_raster.SetProjection(srs.ExportToWkt())# Exports the coordinate system 
			output_raster = None
	logging.info("Toposcale complete!")
	logging.info("%f minutes" % round((time.time()/60 - start_time/60),2) )

#===============================================================================
#	Calling Main
#===============================================================================
if __name__ == '__main__':
	wdir = sys.argv[1]
	mode = sys.argv[2]
	start = sys.argv[3]
	end=sys.argv[4]
	dataset=sys.argv[5]
	member=sys.argv[6]
	main(wdir, mode, start, end, dataset, member)
