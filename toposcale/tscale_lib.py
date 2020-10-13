# tscale_lib

"""
python2

Core tscale methods

Example:


Vars:


Details:


"""

import tscale as ts
import helper as hp
import subprocess
import glob
import pandas as pd
#import tscale2d_era5_src as t2d
#import tscale3d_era5_src as t3d
import recapp_era5 as rc
import solarGeom as sg
from tqdm import tqdm
import netCDF4 as nc
import numpy as np
import os
from joblib import Parallel, delayed
import xarray as xr
from configobj import ConfigObj
import tscale_lib as tlib


def tscale_1d(tscaleSrc, wdir, startDate, endDate, dataset, windCor, plapse):
	'''
	Make listpoints

	Runs: makeListpoints2_points

	Args:
		tmappSrc:
		wdir:
	'''

	cmd = [
			"python",
			tscaleSrc + "/toposcale/tscale_run.py",
			wdir,
			home,
			startDate,
			endDate,
			dataset,
			windCor,
			plapse
		]
	subprocess.check_output(cmd)

def tscale_3d(tscaleSrc, wdir, runmode, startDate, endDate):
	'''
	Make listpoints

	Runs: makeListpoints2_points

	Args:
		tmappSrc:
		wdir:
	'''
	cmd = [
			"python",
			tscaleSrc + "/toposcale/tscale3D.py",
			wdir,
			runmode,
			startDate,
			endDate,
			"HRES",
			'1'

		]

	subprocess.check_output(cmd)



def download_sparse_dem(tmappSrc, wdir, pointsShp, demRes, pointsBuffer):

	'''
	Downloads dem grids corresponding to points ie a sparse network of grids - 
	saves time in case points are very sparse

	Runs: getDEM_points.R

	Args:
		tmappSrc:
		wd:
		demdir:
		pointsShp:
		demRes:
		pointsBuffer
	'''
	cmd = ["Rscript", tmappSrc+"/rsrc/getDEM_points.R" ,wdir , pointsShp,
		str(demRes), str(pointsBuffer)]
	subprocess.check_output(cmd)

def setup_sim_dirs(tmappSrc, wdir):
	'''
	Sets up sim directories.

	Runs: prepClust_pointsSparse.R

	Args:
		tmappSrc:
		wdir:
	'''
	cmd = ["Rscript", tmappSrc+"/rsrc/prepClust_pointsSparse.R", wdir]
	subprocess.check_output(cmd)


def compute_terrain(tmappSrc, home, svfSectors, svfMaxDist):
	'''
	Compute asp,slp,svf from dem.

	Runs: prepClust_pointsSparse.R

	Args:
		tmappSrc:
		wdir:
	'''
	cmd = ["Rscript", tmappSrc+"/rsrc/computeTopo_SVF_points.R", home, svfSectors, svfMaxDist]
	subprocess.check_output(cmd)

def make_listpoints(tmappSrc, home, pointsShp):
	'''
	Make listpoints

	Runs: makeListpoints2_points

	Args:
		tmappSrc:
		wdir:
	'''
	cmd = ["Rscript", tmappSrc+"/rsrc/makeListpoints2_points.R", home, pointsShp]
	subprocess.check_output(cmd)




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

	instRad(sob,3600) # 1h is used as even tho we use 3h or 6h data the value is accumulated over 1h - we do not lose budget in same way as P, just resolution.
	tob.lwin = lwin(sob,tob)
	print("made lwin")
#===============================================================================
# Shortwave -vectorised
#===============================================================================
def swin2D(pob,sob,tob, lp, dtime): # does not work yet (booleen indexing of 2d arraya fauiles as np.min returns single value when we need 161


	timesize=len(dtime)
	statsize=len(lp.ele)
	""" toposcale surface pressure using hypsometric equation - move to own 
	class """
	g=9.81
	R=287.05  # Gas constant for dry air.
	#ztemp = pob.z # geopotential height b = np.transpose(a, (2, 0, 1))
	ztemp = np.transpose(pob.z, (2, 0, 1))
	#Ttemp = pob.t
	Ttemp = np.transpose(pob.t, (2, 0, 1))
	statz = np.array(lp.ele)*g
	#dz=ztemp.transpose()-statz[None,:,None] # transpose not needed but now consistent with CGC surface pressure equations
	dz=ztemp-statz # dimensions of dz : time, levels, stations

	# set all levels below surface to very big number so they canot be found by min
	newdz=dz
	newdz[newdz<0]=999999


	psf =np.zeros( (dtime.size, statz.shape[0]) )

	# loop through timesteps
	for i in range(0,dtime.size):
		
		# find overlying layer
		thisp = dz[i,:,:]==np.min(newdz[i,:,:],axis=0) # thisp is a booleen matrix of levels x stations with true indicating overlying plevel over station surface ele

		# flatten to 1 dimension order='Fortran' or row major
		thispVec =thisp.reshape(thisp.size,order='F')
		TtempVec = Ttemp.reshape(Ttemp.shape[0], Ttemp.shape[1]*Ttemp.shape[2], order='F')
		ztempVec = ztemp.reshape(ztemp.shape[0], ztemp.shape[1]*ztemp.shape[2], order='F')
		
		# booleen indexing to find temp and geopotential that correspond to lowesest overlying layer
		T1=TtempVec[i,thispVec]
		z1=ztempVec[i,thispVec]


		p1=np.tile(pob.levels[::-1],statz.shape[0])[thispVec]*1e2 #Convert to Pa. Reverse levels to ensure low ele (hig pressure) to high elel (low pressure)
		Tbar=np.mean([T1, tob.t[i, :]],axis=0) # temperature midway between surface (toposcale T) and loweset overlying level (T1)
		""" Hypsometric equation.""" #P1 is above surface is this correct? Yes!
		psf[i,:]=(p1*np.exp((z1-statz)*(g/(Tbar*R)))) # exponent is positive ie increases pressure as surface is lower than pressure level

	#psf=np.array(psf).squeeze()




	## Specific humidity routine.
	# mrvf=0.622.*vpf./(psf-vpf); #Mixing ratio for water vapor at subgrid.
	#  qf=mrvf./(1+mrvf); # Specific humidity at subgrid [kg/kg].
	# fout.q(:,:,n)=qf; 


	""" Maybe follow Dubayah's approach (as in Rittger and Girotto) instead
	for the shortwave downscaling, other than terrain effects. """

	""" Height of the "grid" (coarse scale)"""
	Zc=sob.z

	""" toa """
	SWtoa = sob.tisr 

	""" Downwelling shortwave flux of the "grid" using nearest neighbor."""
	SWc=sob.ssrd

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
	SWfdiff=np.array(lp.svf)*SWcdiff
	SWfdiff = np.nan_to_num(SWfdiff) # convert nans (night) to 0

	""" Direct shortwave routine, modified from Joel. 
	Get surface pressure at "grid" (coarse scale). Can remove this
	part once surface pressure field is downloaded, or just check
	for existance. """

	#ztemp = pob.z
	#Ttemp = pob.t
	#dz=ztemp.transpose()-sob.z

	ztemp = np.transpose(pob.z, (0, 2, 1))
	Ttemp = np.transpose(pob.t, (0, 2, 1))
	dz=ztemp-sob.z # dimensions of dz : levels, time, stations

	# set all levels below surface to very big number so they canot be found by min
	newdz=dz
	newdz[newdz<0]=999999

	psc =np.zeros( (dtime.size, statz.shape[0]) )
	for i in range(0,dtime.size):

		#thisp.append(np.argmin(dz[:,i][dz[:,i]>0]))
		# find overlying layer
		thisp = dz[:,i,:]==np.min(newdz[:,i,:],axis=0) # thisp is a booleen matrix of levels x stations with true indicating overlying plevel over station surface ele !! time index in middle this time!!!
		z0 = sob.z[i,:]
		T0 = sob.t2m[i,:]

		# flatten to 1 dimension order='Fortran' or row major
		thispVec =thisp.reshape(thisp.size,order='F')
		TtempVec = Ttemp.reshape(Ttemp.shape[1], Ttemp.shape[0]*Ttemp.shape[2], order='F') # !! order of permutations is different from pressure at finegrid routine (time is middle dimension)
		ztempVec = ztemp.reshape(ztemp.shape[1], ztemp.shape[0]*ztemp.shape[2], order='F')# !! order of permutations is different from pressure at finegrid routine (time is middle dimension)
		
		# booleen indexing to find temp and geopotential that correspond to lowesest overlying layer
		T1=TtempVec[i,thispVec]
		z1=ztempVec[i,thispVec]

		p1=np.tile(pob.levels[::-1],statz.shape[0])[thispVec]*1e2 #Convert to Pa.
		Tbar=np.mean([T0, T1],axis=0)
		""" Hypsometric equation."""
		psc[i,:] = (p1*np.exp((z1-z0)*(g/(Tbar*R))))



	"""compute julian dtime"""
	jd= sg.to_jd(dtime)

	"""
	Calculates a unit vector in the direction of the sun from the observer 
	position.
	"""
	#sunv=sg.sunvector(jd=jd, latitude=stat.lat, longitude=stat.lon, timezone=stat.tz)
	svx,svy,svz =	sg.sunvectorMD(jd, lp.lat, lp.lon, lp.tz,statsize,timesize)
	"""
	Computes azimuth , zenith  and sun elevation 
	for each timestamp !!! NOW THIS NEEDS ACCEPT SEPARATE VECTORS!!
	"""
	sp=sg.sunposMD(svx,svy,svz)


	# Cosine of the zenith angle.
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
	ka = np.nan_to_num(ka) 


	# Note this equation is obtained by inverting Beer's law, i.e. use
	#I_0=I_inf x exp[(-ka/mu) int_z0**inf rho dz]
	# Along with hydrostatic equation to convert to pressure coordinates then
	# solve for ka using p(z=inf)=0.


	# Now you can (finally) find the direct component at subgrid. 
	SWfdir=SWtoa*np.exp(-ka*psf/(g*muz))

	""" Then perform the terrain correction. [Corripio 2003 / rpackage insol port]."""

	"""compute mean horizon elevation - why negative hor.el possible??? """
	horel=(((np.arccos(np.sqrt(lp.svf))*180)/np.pi)*2)-lp.slp
	horel[horel<0]=0


	"""
	normal vector - Calculates a unit vector normal to a surface defined by 
	slope inclination and slope orientation.
	"""
	nv = sg.normalvector(slope=lp.slp, aspect=lp.asp)

	"""
	Method 1: Computes the intensity according to the position of the sun (sunv) and 
	dotproduct normal vector to slope.
	From corripio r package
	"""

	""" need to reconstruct sunv matrix for multipoint case here"""
	sunv=np.array((svx,svy,svz))

	#dotprod=np.dot(sunv ,np.transpose(nv))\
	dotprod=np.tensordot(sunv  ,nv , axes=([0],[1]))

	dprod=dotprod[:,:,0]
	dprod[dprod<0] = 0 #negative indicates selfshading
	#dprod = dprod

	"""Method 2: Illumination angles. Dozier""" #SEE tscale.py
	#saz=sp.azi

	# a=np.array(stat.asp)
	# b =np.tile(a,timesize)
	# asp=b.reshape(timesize,statsize)
	# a=np.array(stat.slp)
	# b =np.tile(a,timesize)
	# slp=b.reshape(timesize,statsize)
	# cosis=muz*np.cos(slp)+np.sin(sp.zen)*np.sin(slp)*np.cos(sp.azi-asp)# cosine of illumination angle at subgrid.
	# cosic=muz # cosine of illumination angle at grid (slope=0).

	"""
	SUN ELEVATION below hor.el set to 0 - binary mask
	"""
	a=np.array(horel)
	b =np.tile(a,timesize)
	horel2=b.reshape(timesize,statsize)

	selMask = sp.sel
	selMask[selMask<horel2]=0
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

	return SWfglob, psf
	""" 
	Missing components
	- terrain reflection
	"""


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
		for i in range(0,dates.size):
		
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
		#ka[ka == -np.inf] = 0
		#ka[ka == np.inf] = 0
		ka = np.nan_to_num(ka) 

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
		return SWfglob, psf
		""" 
		Missing components
		- terrain reflection
		"""
		# init grid stack

def snowPartition(TA,PINT): 

	# partition prate to rain snow (mm/hr)	
	lowthresh=272.15
	highthresh = 274.15
	d = {'prate': PINT, 'ta': TA}
	df = pd.DataFrame(data=d)
	snow = df.prate.where(df.ta<lowthresh) 
	rain=df.prate.where(df.ta>=highthresh) 

	mix1S = df.prate.where((df.ta >= lowthresh) & (df.ta<=highthresh), inplace = False)
	mix1T = df.ta.where((df.ta >= lowthresh) & (df.ta<=highthresh), inplace = False)
	mixSno=(highthresh - mix1T) / (highthresh-lowthresh)
	mixRain=1.-mixSno
	addSnow=mix1S*mixSno
	addRain=mix1S*mixRain

	# nas to 0
	snow[np.isnan(snow)] = 0 	
	rain[np.isnan(rain)] = 0 
	addRain[np.isnan(addRain)] = 0 
	addSnow[np.isnan(addSnow)] = 0 

	snowTot=(snow+addSnow)
	rainTot=rain + addRain

	return snowTot,rainTot

def tscale3dmain(mymonth):
	'''
	main function that runs tscale3D for a single month of era5 input. This
	avoids massive i/o overhead of eg. 72gb PLEV.nc files.
	
	Args: mymonth (str)
	
	Example tscale3dmain("_197909.nc")
	'''
	surfile=wdir+'/forcing/SURF_'+mymonth
	plevfile=wdir+ '/forcing/PLEV_'+mymonth

	#===========================================================================
	# Init t3d object
	#===========================================================================
	f = nc.Dataset(plevfile)

	lev = f.variables['level'][:]
	var='t'

	t3d = rc.t3d( pl=plevfile, sa=surfile)
	t3d.addTime()
	a=t3d.dtime[2]-t3d.dtime[1]
	t3d.step = a.seconds
	t3d.lp=lp
	#timesteps = t3d.pl.variables['time'][:].size
	t3d.out_xyz_dem, t3d.lats, t3d.lons, t3d.shape, t3d.names= t3d.demGrid(stations=lp)

	#===========================================================================
	# pob
	#===========================================================================
	xdim=lev.shape[0]
	ydim=t3d.shape[0]

	t_interp_out = np.zeros((xdim,ydim))
	z_interp_out = np.zeros((xdim,ydim))


	for timestep in range(len(t3d.dtime )):
		#print(str(round(float(timestep)/float(timesteps)*100,0))+ "% done")

		gridT,gridZ,gridLat,gridLon=t3d.gridValue(var,timestep)

		t_interp, z_interp = t3d.inLevelInterp(gridT,gridZ, gridLat,gridLon,t3d.out_xyz_dem)

		t_interp_out = np.dstack((t_interp_out, t_interp))
		z_interp_out = np.dstack((z_interp_out, z_interp))

	# drop init blank layer
	tinterp =t_interp_out[:,:,1:]
	zinterp =z_interp_out[:,:,1:]
	pob= hp.Bunch(t=tinterp,z=zinterp, levels=lev)

	print("made a pob!")

	#===========================================================================
	# tob
	#===========================================================================

	starti=0
	endi=len(t3d.dtime)
	t =t3d.tscale3d('t',starti, endi)
	r =t3d.tscale3d('r',starti, endi)
	u =t3d.tscale3d('u',starti, endi)
	v =t3d.tscale3d('v',starti, endi)

	# compute wind speed and direction
	ws = np.sqrt(u**2+v**2)
	wd =   (180 / np.pi) * np.arctan(u/v) + np.where(v>0,180,np.where(u>0,360,0))

	tob = hp.Bunch(t=t,r=r,u=u, v=v, ws=ws,wd=wd, dtime=t3d.dtime)

	# Physical filters

	# constrain RH to 0-100 interval
	# constrain r to interval 5-100 - do here as required by LWin parameterisation
	tob.r[tob.r <5]=5
	tob.r[tob.r>100]=100
	tob.ws[tob.ws<0]=0

	print("made a tob!")

	#===========================================================================
	# sob
	#===========================================================================

	t2m = t3d.tscale2d('t2m', starti,endi)
	tp = t3d.tscale2d('tp', starti,endi)
	ssrd = t3d.tscale2d('ssrd', starti,endi)
	strd = t3d.tscale2d('strd', starti,endi)
	tisr = t3d.tscale2d('tisr', starti,endi)
	d2m = t3d.tscale2d('d2m', starti,endi)
	z = t3d.tscale2d('z', starti,endi)

	gridEle=z[0,:]/g

	sob = hp.Bunch(t2m=t2m, tp=tp, ssrd=ssrd, strd=strd, tisr=tisr, d2m=d2m, z=z, gridEle=gridEle, dtime=t3d.dtime)
	print("made a sob!")

	#===============================================================================
	# Precip
	#===============================================================================

	pmmhr = tp2rate(tp,t3d.step)
	sob.pmmhr = pmmhr
	tob.prate , tob.psum = precipPoint(lp.ele, sob)
	print("made prate!")

	instRad(sob,3600) # 1h is used as even tho we use 3h or 6h data the value is accumulated over 1h - we do not lose budget in same way as P, just resolution.
	tob.lwin = lwin(sob,tob)
	print("made lwin")

	# vector method
	tob.swin, tob.psf = swin2D(pob,sob,tob, lp, dtime=t3d.dtime)

	# loop method
	# init first row
	# ntimestamps=ds.dtime.shape[0]
	# ts_swin = np.zeros((ntimestamps))
	# for stationi in range(0, lp.shape[0]):
	# 	print stationi 
	# 	'''here we test for array full of NaNs due to points being 
	# 	outside of grid. The NaN arrays are created in the 
	# 	interpolation but only break the code here. If array is
	# 		all NaN then we just fill that stations slot 
	# 	with NaN'''
	# 	testNans = np.count_nonzero(~np.isnan(pob.z[:,stationi,:]))
	# 	if testNans != 0:
	# 		ts= swin1D(pob=pob,sob=sob,tob=tob, stat=lp, dates=ds.dtime, index=stationi)
	# 		ts_swin=np.column_stack((ts_swin,ts))
	# 	if testNans == 0:
	# 		nan_vec = np.empty(ntimestamps) * np.nan
	# 		ts_swin=np.column_stack((ts_swin,nan_vec))
		

	# # # drop init row
	# tob.swin =ts_swin[:,1:]
	print("Made Swin")




	#===========================================================================
	# make dataframe (write individual files plus netcdf)
	#===========================================================================
	start=1
	print("Writing toposcale files...")
	for i in range(0,tob.t.shape[1]):

		# partition
		Sf,Rf = snowPartition(tob.t[:,i], tob.prate[i,:])

		df = pd.DataFrame({	"TA":tob.t[:,i], 
					"RH":tob.r[:,i],
					"P":tob.psf[:,i],
					"VW":tob.ws[:,i],
					"DW":tob.wd[:,i],
					"ILWR":tob.lwin[:,i], 
					"ISWR":tob.swin[:,i], 
					"PINT":tob.prate[i,:],
					"PSUM":tob.psum[i,:],
					"Snowf":np.array(Sf),
					"Rainf":np.array(Rf)
					},index=tob.dtime)
		df.index.name="datetime"

		# fill outstanding nan in SW routine with 0 (night)
		df.ISWR = df.ISWR.fillna(0)




		fileout=wdir+"/out/meteo"+"c"+str(i+1)+"_"+mymonth+".csv" # convert python index back to g* sim dir index
		column_order = ['TA', 'RH', 'P', 'VW', 'DW', 'ILWR', 'ISWR', 'PINT', 'PSUM', 'Snowf', 'Rainf']
		df[column_order].to_csv(path_or_buf=fileout ,na_rep=-999,float_format='%.3f')
		print("written "+fileout)




def split_netcdf(file2split):
	''' splits concatenated file (SURF.nc, PLEV.ncin to monthly files 
	as originally downloaded from CDS, for use in tscale_fast'''
	basename = file2split.split(".nc")[0]
	ds = xr.open_dataset(file2split) 
	dates, datasets = zip(*ds.resample(time='1M').mean('time').groupby('time'))
	filenames = [basename+"_"+pd.to_datetime(date).strftime('%Y%m') + '.nc' for date in dates]
	xr.save_mfdataset(datasets, filenames)