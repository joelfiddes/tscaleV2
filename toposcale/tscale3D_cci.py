"""
Modified version of tscale3D.py specifically for CCI permafrost. 
 - Accepts single timestep as argument.
 - Removes dependency on slp,asp,svf that are not avvailable in CCI. This effects LW (svf) and SWin (Horizon, selfshading).

Args:
	wdir = toplevel directory of simulation
	start = timstep to compute in ISO
	tz = timezone of location (not era5, which is UTC+0) used for solar calcs
          
Returns:
	single nc file of each variable at timestep
          
Example: 
	python tscale3D_cci.py /home/joel/sim/tscale3Dbig/ '2013-09-02 06:00:00'


"""
import os
import sys
import pandas as pd
import netCDF4 as nc
import numpy as np
import logging
from datetime import datetime
import helper as hp
import solarGeom as sg
#import era5
import tscale2d_era5_src as t2d
import tscale3d_era5_src as t3d
import recapp_era5 as rc
import time


# wdir="/home/joel/sim/tscale3Dbig/"
# start="2013-09-01 06:00:00"

def main(wdir, start, tz):

	# fixed stuff
	dataset="HRES"
	tz=int(tz)
	g=9.81

	#Set up Log
	logs = wdir+"/logs/"
	if not os.path.exists(logs):
		os.makedirs(logs)
	logfile=logs+"/logfile"+start
	logging.basicConfig(level=logging.DEBUG, filename=logfile, filemode="a+",format="%(asctime)-15s %(levelname)-8s %(message)s")
	
	#start timer
	start_time = time.time()

	# standard inputs
	demfile = wdir+'/predictors/ele.nc'
	surfile=wdir+'/forcing/SURF.nc'
	plevfile=wdir+ '/forcing/PLEV.nc'
	
	# make out path for results
	out = wdir+"/out/"
	if not os.path.exists(out):
		os.makedirs(out)

	# open DEM
	dem  = nck.Dataset(demfile)
	dem_ele = dem.variables['ele'][:]

	# time stuff
	f = nc.Dataset( plevfile)
	nctime = f.variables['time']
	dtime = pd.to_datetime(nc.num2date(nctime[:],nctime.units, calendar="standard"))
	starti = np.asscalar(np.where(dtime==start)[0]) # so can run on a single timestep
	endi = starti+1

	# compute timestep before we cut timeseries
	a=dtime[2]-dtime[1]
	step = a.seconds

	# extract timestep
	dtime= dtime[starti:endi,] 

	print("Running timestep "+ start)

	#===============================================================================
	# tscale3d - 3D interpolation of pressure level fields
	#===============================================================================

	t = t3d.main( wdir, 'grid', 't', starti,endi,dataset)
	r = t3d.main( wdir, 'grid', 'r', starti,endi,dataset)
	u = t3d.main( wdir, 'grid', 'u', starti,endi,dataset)
	v = t3d.main( wdir, 'grid', 'v', starti,endi,dataset)
	gtob = hp.Bunch(t=t,r=r,u=u,v=v, dtime=dtime)

	#===============================================================================
	# tscale2d - Generates 2D interpolations from coarse (ERA5) to fine (1km) grid
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

	#===============================================================================
	# Precip
	#===============================================================================

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

	def precipGrid(fineEle, gsob):
		'''

		'''
		
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


	gsob.pmmhr = tp2rate(tp,step)
	grid_prate = precipGrid(dem_ele,gsob)
	

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
		and therefore treated here the same ie step=3600s (1h)
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


#===============================================================================
# make a pob - required as input to swin routine. This object is data on each 
# pressure level interpolated to the fine grid ie has dimensions xy (finegrid) x plev
#===============================================================================

	f = nc.Dataset(plevfile)
	lev = f.variables['level'][:]
	var='t'
	ds = rc.t3d( pl=plevfile, dem =demfile)
	out_xyz_dem, lats, lons, shape= ds.demGrid()
	xdim=lev.shape[0]

	ydim=shape[0]*shape[1]
	t_interp_out = np.zeros((xdim,ydim))
	z_interp_out = np.zeros((xdim,ydim))
		
	for timestep in range(starti,endi):
		gridT,gridZ,gridLat,gridLon=ds.gridValue(var,timestep)
		t_interp, z_interp = ds.inLevelInterp(gridT,gridZ, gridLat,gridLon,out_xyz_dem)
		t_interp_out = np.dstack((t_interp_out, t_interp))
		z_interp_out = np.dstack((z_interp_out, z_interp))

	# drop init blank layer
	tinterp =t_interp_out[:,:,1:]
	zinterp =z_interp_out[:,:,1:]

	gpob= hp.Bunch(t=tinterp,z=zinterp, levels=lev)

	logging.info("made a pob!")

	dem  = nc.Dataset(demfile)
	lon = dem.variables['longitude'][:]
	lat = dem.variables['latitude'][:]
	demv = dem.variables['ele'][:]
	demlon1 = dem.variables['longitude'][:]
	demlat1 = dem.variables['latitude'][:]
	demlon = np.tile(demlon1,demlat1.size)
	demlat = np.repeat(demlat1,demlon1.size)
	demv=np.reshape(demv,demv.size)

	# why are these masked values generated?
	demv =np.ma.filled(demv, fill_value=1)
	tz = np.repeat(tz,demv.size)

	stat = pd.DataFrame({	"ele":demv, 
					"lon":demlon,
					"lat":demlat,
					"tz":tz					
					})
#===============================================================================
# Compute Shortwave
#===============================================================================	

	def swin2D(pob,sob,tob, stat, dates): 
		'''
		main edit over standard function for points:
		- 3d tob,sob reduce to 2D (reshape)
		'''
		
		timesize=len(dates)
		statsize=len(stat.ele)
		""" toposcale surface pressure using hypsometric equation - move to own 
		class """
		g=9.81
		R=287.05  # Gas constant for dry air.
		#ztemp = pob.z # geopotential height b = np.transpose(a, (2, 0, 1))
		ztemp = np.transpose(pob.z, (2, 0, 1)) # pob is originally ordered levels,stations/cells,time
		#Ttemp = pob.t
		Ttemp = np.transpose(pob.t, (2, 0, 1))# pob is originally ordered levels,stations/cells,time
		statz = np.array(stat.ele)*g
		#dz=ztemp.transpose()-statz[None,:,None] # transpose not needed but now consistent with CGC surface pressure equations
		dz=ztemp-statz # dimensions of dz : time, levels, stations
		
		# set all levels below surface to very big number so they canot be found by min
		newdz=dz
		newdz[newdz<0]=999999
		
		psf =np.zeros( (dates.size, statz.shape[0]) )

		# reshape tob.t here
		tob.tr=tob.t.reshape(tob.t.shape[0]*tob.t.shape[1], tob.t.shape[2], order='F')
		tob.trT=tob.tr.T # transpose to get right order


		# loop through timesteps
		for i in range(0,dates.size):
			
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
			Tbar=np.mean([T1, tob.trT[i, :]],axis=0) # temperature midway between surface (toposcale T) and loweset overlying level (T1)
			""" Hypsometric equation.""" #P1 is above surface is this correct? Yes!
			psf[i,:]=(p1*np.exp((z1-statz)*(g/(Tbar*R)))) # exponent is positive ie increases pressure as surface is lower than pressure level


		""" Maybe follow Dubayah's approach (as in Rittger and Girotto) instead
		for the shortwave downscaling, other than terrain effects. """

		""" Height of the "grid" (coarse scale)"""
		Zc=sob.z.reshape(sob.z.shape[0]*sob.z.shape[1], sob.z.shape[2]).T # reshape and transpose to remove dimension and make broadcastable

		""" toa """
		SWtoa = sob.tisr.reshape(sob.tisr.shape[0]*sob.tisr.shape[1], sob.tisr.shape[2]).T  

		""" Downwelling shortwave flux of the "grid" using nearest neighbor."""
		SWc=sob.ssrd.reshape(sob.ssrd.shape[0]*sob.ssrd.shape[1], sob.ssrd.shape[2]).T 

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
		SWfdiff=SWcdiff
		SWfdiff = np.nan_to_num(SWfdiff) # convert nans (night) to 0

		""" Direct shortwave routine, modified from Joel. 
		Get surface pressure at "grid" (coarse scale). Can remove this
		part once surface pressure field is downloaded, or just check
		for existance. """

		ztemp = np.transpose(pob.z, (0, 2, 1))
		Ttemp = np.transpose(pob.t, (0, 2, 1))
		dz=ztemp-Zc # dimensions of dz : levels, time, stations

		# set all levels below surface to very big number so they canot be found by min
		newdz=dz
		newdz[newdz<0]=999999

		psc =np.zeros( (dates.size, statz.shape[0]) )
		for i in range(0,dates.size):
		
			#thisp.append(np.argmin(dz[:,i][dz[:,i]>0]))
			# find overlying layer
			thisp = dz[:,i,:]==np.min(newdz[:,i,:],axis=0) # thisp is a booleen matrix of levels x stations with true indicating overlying plevel over station surface ele !! time index in middle this time!!!
			z0 = Zc[i,:]
			T0 = sob.t2m.reshape(sob.t2m.shape[0]*sob.t2m.shape[1], sob.t2m.shape[2]).T[i,:]

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


		
		"""compute julian dates"""
		jd= sg.to_jd(dates)

		"""
		Calculates a unit vector in the direction of the sun from the observer 
		position.
		"""
		svx,svy,svz =	sg.sunvectorMD(jd, stat.lat, stat.lon, stat.tz,statsize,timesize)
		sp=sg.sunposMD(svx,svy,svz)


		# Cosine of the zenith angle.
		#sp.zen=sp.zen*(sp.zen>0) # Sun might be below the horizon.
		muz=np.cos(sp.zen) 
		muz = muz
		# NB! psc must be in Pa (NOT hPA!).
	  
		# Calculate the "broadband" absorption coefficient. Elevation correction
		ka=(g*muz/(psc))*np.log(SWtoa/SWcdir)	
		ka = np.nan_to_num(ka) 

		# Now you can (finally) find the direct component at subgrid. 
		SWfdir=SWtoa*np.exp(-ka*psf/(g*muz))
		SWfdirCor=SWfdir #*dprod
	  
		SWfglob =  SWfdiff+ SWfdirCor
		print(" %f minutes for VECTORISED interpolation %s" % (round((time.time()/60 - start_time/60),2),"swin") )
		return SWfglob

	gtob.swin = swin2D(gpob,gsob,gtob, stat, dtime)
	gtob.swin =gtob.swin.reshape(gtob.swin.shape[0], gsob.ssrd.shape[0], gsob.ssrd.shape[1])

	ntime = len(dtime)
	stephr =a.seconds/60/60
	rtime=np.array(range(len(dtime)))*stephr


#===============================================================================
# write results
#===============================================================================	

	
	## Longwave

	#open
	f = nc.Dataset('lwin.nc','w', format='NETCDF4')

	#make dimensions
	f.createDimension('lon', len(lon))
	f.createDimension('lat', len(lat))
	f.createDimension('time', ntime)

	#make dimension variables
	longitude = f.createVariable('lon',    'f4',('lon',))
	latitude  = f.createVariable('lat',    'f4',('lat',))
	mytime = f.createVariable('time', 'i', ('time',))
	lwin = f.createVariable('lwin',    'f4',('lon','lat','time'))

	#assign dimensions
	longitude[:] = lon
	latitude[:]  = lat
	mytime[:] = rtime

	varT = np.transpose(ts_lwin ,(1, 0, 2))
	lwin[:] = np.flip(varT, 1)
	
	#metadata
	f.history = 'Created by toposcale on '+time.ctime()
	mytime.units = 'hours since '+str(dtime[0])

	f.close()

	## Shortwave

	#open
	f = nc.Dataset('swin.nc','w', format='NETCDF4')

	#make dimensions
	f.createDimension('lon', len(lon))
	f.createDimension('lat', len(lat))
	f.createDimension('time', ntime)

	#make dimension variables
	longitude = f.createVariable('lon',    'f4',('lon',))
	latitude  = f.createVariable('lat',    'f4',('lat',))
	mytime = f.createVariable('time', 'i', ('time',))
	lwin = f.createVariable('swin',    'f4',('lon','lat','time'))

	#assign dimensions
	longitude[:] = lon
	latitude[:]  = lat
	mytime[:] = rtime

	varT = np.transpose(gtob.swin ,(2, 1, 0)) #t,lat,lon -> lon,lat,t
	lwin[:] = varT
	
	#metadata
	f.history = 'Created by toposcale on '+time.ctime()
	mytime.units = 'hours since '+str(dtime[0])

	f.close()



	## ta
	
	#open
	f = nc.Dataset('ta.nc','w', format='NETCDF4')

	#make dimensions
	f.createDimension('lon', len(lon))
	f.createDimension('lat', len(lat))
	f.createDimension('time', ntime)

	#make dimension variables
	longitude = f.createVariable('lon',    'f4',('lon',))
	latitude  = f.createVariable('lat',    'f4',('lat',))
	mytime = f.createVariable('time', 'i', ('time',))
	lwin = f.createVariable('ta',    'f4',('lon','lat','time'))

	#assign dimensions
	longitude[:] = lon
	latitude[:]  = lat
	mytime[:] = rtime

	varT = np.transpose(gtob.t ,(1, 0, 2))
	lwin[:] = np.flip(varT, 1)
	
	#metadata
	f.history = 'Created by toposcale on '+time.ctime()
	mytime.units = 'hours since '+str(dtime[0])

	f.close()

	# rh
	
	#open
	f = nc.Dataset('rh.nc','w', format='NETCDF4')

	#make dimensions
	f.createDimension('lon', len(lon))
	f.createDimension('lat', len(lat))
	f.createDimension('time', ntime)

	#make dimension variables
	longitude = f.createVariable('lon',    'f4',('lon',))
	latitude  = f.createVariable('lat',    'f4',('lat',))
	mytime = f.createVariable('time', 'i', ('time',))
	lwin = f.createVariable('rh',    'f4',('lon','lat','time'))

	#assign dimensions
	longitude[:] = lon
	latitude[:]  = lat
	mytime[:] = rtime

	varT = np.transpose(gtob.r ,(1, 0, 2))
	lwin[:] = np.flip(varT, 1)
	
	#metadata
	f.history = 'Created by toposcale on '+time.ctime()
	mytime.units = 'hours since '+str(dtime[0])

	f.close()

	# u
	
	#open
	f = nc.Dataset('u.nc','w', format='NETCDF4')

	#make dimensions
	f.createDimension('lon', len(lon))
	f.createDimension('lat', len(lat))
	f.createDimension('time', ntime)

	#make dimension variables
	longitude = f.createVariable('lon',    'f4',('lon',))
	latitude  = f.createVariable('lat',    'f4',('lat',))
	mytime = f.createVariable('time', 'i', ('time',))
	lwin = f.createVariable('u',    'f4',('lon','lat','time'))

	#assign dimensions
	longitude[:] = lon
	latitude[:]  = lat
	mytime[:] = rtime

	varT = np.transpose(gtob.u ,(1, 0, 2))
	lwin[:] = np.flip(varT, 1)
	
	#metadata
	f.history = 'Created by toposcale on '+time.ctime()
	mytime.units = 'hours since '+str(dtime[0])

	f.close()

	# v
	
	#open
	f = nc.Dataset('v.nc','w', format='NETCDF4')

	#make dimensions
	f.createDimension('lon', len(lon))
	f.createDimension('lat', len(lat))
	f.createDimension('time', ntime)

	#make dimension variables
	longitude = f.createVariable('lon',    'f4',('lon',))
	latitude  = f.createVariable('lat',    'f4',('lat',))
	mytime = f.createVariable('time', 'i', ('time',))
	lwin = f.createVariable('v',    'f4',('lon','lat','time'))

	#assign dimensions
	longitude[:] = lon
	latitude[:]  = lat
	mytime[:] = rtime

	varT = np.transpose(gtob.v ,(1, 0, 2))
	lwin[:] = np.flip(varT, 1)
	
	#metadata
	f.history = 'Created by toposcale on '+time.ctime()
	mytime.units = 'hours since '+str(dtime[0])

	f.close()

	#open
	f = nc.Dataset('pint.nc','w', format='NETCDF4')

	#make dimensions
	f.createDimension('lon', len(lon))
	f.createDimension('lat', len(lat))
	f.createDimension('time', ntime)

	#make dimension variables
	longitude = f.createVariable('lon',    'f4',('lon',))
	latitude  = f.createVariable('lat',    'f4',('lat',))
	mytime = f.createVariable('time', 'i', ('time',))
	lwin = f.createVariable('pint',    'f4',('lon','lat','time'))

	#assign dimensions
	longitude[:] = lon
	latitude[:]  = lat
	mytime[:] = rtime

	varT = np.transpose(grid_prate ,(1, 0, 2))
	lwin[:] = np.flip(varT, 1)
	
	#metadata
	f.history = 'Created by toposcale on '+time.ctime()
	mytime.units = 'hours since '+str(dtime[0])

	f.close()
	

	logging.info("Toposcale complete!")
	#print("%f minutes" % round((time.time()/60 - start_time/60),2) )
	print("%f seconds for run" % round((time.time() - start_time),2) )




#===============================================================================
#	Calling Main
#===============================================================================
if __name__ == '__main__':
	wdir = sys.argv[1]
	start = sys.argv[2]
	tz = sys.argv[3]
	main(wdir, start,tz )
