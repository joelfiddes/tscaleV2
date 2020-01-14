"""
Modified version of tscale3D.py specifically for CCI permafrost. 
- Accepts single timestep as argument.
- Removes dependency on slp,asp,svf that are not avvailable in CCI. This effects LW (svf) and SWin (Horizon, selfshading).

update Dec 2019:
- accepts assci list
- 9 day chunks
- writes netcdf
- includes sublimation
- reads global single parameter files
- writes 9 days netcdf

Args:
	 coords = '/home/joel/sim/cci_perm_final/coordinates.dat'
	 eraDir = '/home/joel/sim/cci_perm_final/era5'
	 outDir = '/home/joel/sim/cci_perm_final/era5' + '/out'
	 start='1980-01-01 00:00:00'
	 end='1980-01-10 00:00:00'

Example:
	tscale_cci.py  "/home/joel/sim/cci_perm_final/coordinates.dat" "/home/joel/sim/cci_perm_final/era5/"  "/home/joel/sim/cci_perm_final/era5/out/" '1980-01-01 00:00:00' '1980-01-10 00:00:00'

"""

import sys
import os
import pandas as pd
import netCDF4 as nc
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from bisect import bisect_left
from bisect import bisect_right
import helper as hp
import solarGeom as sg


def main(coords,eraDir, outDir,start, end, startIndex):
	print(start)
	print(end)
	g=9.81 # geopotential constant
	tz=0 # timezone always utc0 for era5 data

	zp_file=eraDir+"/PLEV_geopotential_1980.nc"

	plevDict={
	eraDir+"/PLEV_temperature_1980.nc" : "t",
	eraDir+"/PLEV_u_component_of_wind_1980.nc": "u",
	eraDir+"/PLEV_v_component_of_wind_1980.nc": "v",
	eraDir+"/PLEV_relative_humidity_1980.nc" : "r"
	}

	surfDict={
	eraDir+"/SURF_2m_temperature_1980.nc" : "t2m",
	eraDir+"/SURF_2m_dewpoint_temperature_1980.nc": "d2m",
	eraDir+"/SURF_geopotential_1980.nc": "z",
	eraDir+"/SURF_surface_solar_radiation_downwards_1980.nc" : "ssrd",
	eraDir+"/SURF_surface_thermal_radiation_downwards_1980.nc" :"strd",
	eraDir+"/SURF_Total_precipitation_1980.nc" : "tp",
	eraDir+"/SURF_TOA_incident_solar_radiation_1980.nc": "tisr"
	}

	# read in lispoints
	lp=pd.read_csv(coords)

	# make out path for results
	out = outDir
	if not os.path.exists(out):
		os.makedirs(out)

	# time stuff
	f = nc.Dataset(zp_file)
	nctime = f.variables['time']
	dtime = pd.to_datetime(nc.num2date(nctime[:],nctime.units, calendar="standard"))
	starti = np.asscalar(np.where(dtime==start)[0]) # so can run on a single timestep
	print(end)
	endi = np.asscalar(np.where(dtime==end)[0]) 
	year = dtime.year[0] 

	# compute timestep before we cut timeseries
	a=dtime[2]-dtime[1]
	step = a.seconds
	stephr=step/(60*60)
	# extract timestep
	dtime= dtime[starti:endi,] 

	print("Running timestep "+ start)

	#===============================================================================
	# tscale3d - 3D interpolation of pressure level fields
	#===============================================================================
	ele= lp.iloc[:,2]#[:,3]
	lats = lp.iloc[:,0]#[:,2] #[s['lat'] for s in stations]
	lons = lp.iloc[:,1] #[:,1] #[s['lon'] for s in stations]

	lp = hp.Bunch(ele=ele, lat=lats, lon=lons)
	#idstat = stations.id#[:,0]#[s['id'] for s in stations]
	out_xyz_dem = np.asarray([lats,lons,ele*g],order="F" ).transpose()

	# init gtob object
	gtob = hp.Bunch(time=dtime)

	for plev in plevDict:

		t= nc.Dataset(plev) # key is filename
		varname=plevDict[plev] # value of key is par shortname
		z=nc.Dataset(zp_file) # could be outside loop

		# init grid stack
		xdim=out_xyz_dem.shape[0]
		sa_vec = np.zeros(xdim)
		#names=1:n

		for timestep in range(starti, endi):

			"""
	        Return original grid temperatures and geopotential of differnet
	        pressure levels. The function are called by inLevelInterp() to
	        get the input ERA-Interim values.
	        
	        Args: 
	            variable: Given interpolated climate variable
	            timestep: Time need to be interpolated. Time is in interger (e.g.
	            0, 1, 2)
	            
	        Returns:
	            gridT: Grid temperatures of different pressure levels. Retruned 
	            temperature are formated in [level, lat, lon]
	            gridZ: Grid geopotential of different pressure levels. Retruned 
	            temperature are formated in [level, lat, lon]
	            gridLon: Grid longitude of pressure level variables
	            gridLat: Grid latitude of pressure level variables
	        
	        Example:
	            gridT,gridZ,gridLat,gridLon=downscaling.gridValue('Temperature',0)
	            
			"""

			gridT = t.variables[varname][timestep,:,:,:]
			gridZ = z.variables['z'][timestep,:,:,:]
			gridLat = t['latitude'][:]
			gridLat = gridLat[::-1] # reverse to deal with ERA5 order
			gridLon = t['longitude'][:]
			gridLev = t.variables['level'][::-1] 
	        #return gridT,gridZ,gridLat,gridLon

			"""
			This is a 2D interpolatation, and returns interpolated temperatures
			of different pressure levels.

			Interpolated domain is smaller than original domain - original (ERA5) domain
			should be one cell larger than expected point or grid domain.

			Args:
			    gridT: Grid temperatures of different pressure levels. Retruned 
			        temperature are formated in [level, lat, lon]
			    gridZ: Grid geopotential of different pressure levels. Retruned 
			        temperature are formated in [level, lat, lon]
			    gridLat: Grid longitude of pressure level variables
			    gridLon: Grid latitude of pressure level variables
			    out_xyz: Given sites, which will be interpolated.
			    
			Returns:
			    t_interp: Interpolated temperatre of different pressure levels. 
			        The returned values are fomrated in [level, lat, lon]
			    z_interp: Interpolated geopotential of different pressure levels. 
			        The returned values are fomrated in [level, lat, lon]

			Examples:
			    downscaling = downscaling(dem, geop, sa, pl)

			    out_xyz_dem, lats, lons, shape = downscaling.demGrid()
			    out_xyz_sur = downscaling.surGrid(lats, lons, None)

			    #interpolate 2-meter temperature
			    surTa = downscaling.surTa(0, out_xyz_sur)
			    #original ERA-I values
			    gridT,gridZ,gridLat,gridLon = downscaling.gridValue(variable,0)
			    #interpolate temperatures and geopotential of different 
			    pressure levels.

			    t_interp, z_interp = downscaling.inLevelInterp(gridT,gridZ,
			                                                   gridLat,gridLon,
			                                                   out_xyz_dem)
			"""



			shape = gridT.shape
			#create array to hold interpolation resultes
			t_interp = np.zeros([shape[0], len(out_xyz_dem)])
			z_interp = np.zeros([shape[0], len(out_xyz_dem)]) # HOW MANY TIMES DO WE REALLY NEED TO COMPUTE THIS?

			#temperatue and elevation interpolation 2d
			for i in range(shape[0]):
				ft = RegularGridInterpolator((gridLat,gridLon), 
			                                  gridT[i,:,:], 'linear', bounds_error=False)
				fz = RegularGridInterpolator((gridLat,gridLon), 
			                                  gridZ[i,:,:], 'linear', bounds_error=False)
				t_interp[i,:] = ft(out_xyz_dem[:,:2] )#temperature

				z_interp[i,:] = fz(out_xyz_dem[:,:2])#elevation

				# invert pressure levels
				#t_interp = t_interp[::-1,:]
				#z_interp = z_interp[::-1,:]

				"""This is a 1D interpoation. The function return interpolated 
				upper air temperature at the given sites by 
				interpolation between different pressure levels.
				"""
			ele = out_xyz_dem[:,2]
			size = np.arange(out_xyz_dem.shape[0])
			n = [bisect_left(z_interp[:,i], ele[i]) for i in size]
			n = [x+1 if x == 0 else x for x in n]

			lowN = [l-1 for l in n]

			upperT = t_interp[n,size]
			upperZ = z_interp[n,size]
			dG  = upperT-t_interp[lowN,size]#<0
			dG /= upperZ-z_interp[lowN,size]#<0
			dG *= out_xyz_dem[:,2] - upperZ#>0
			dG += upperT
			     
			pl_obs=dG


			sa_vec=np.column_stack((sa_vec,pl_obs))
		
		# drop init row
		sa_vec =sa_vec[:,1:]
		# rename to variable
		setattr( gtob, varname, sa_vec)
	print("t,r,u,v done")
	#===============================================================================
	# tscale2d - Generates 2D interpolations from coarse (ERA5) to fine (1km) grid
	#===============================================================================
	gsob = hp.Bunch(dtime=dtime)
	# init grid stack
	xdim=shape[0]
	sa_vec = np.zeros((xdim))

	for surf in surfDict:
		
		t= nc.Dataset(surf) # key is filename
		varname=surfDict[surf] # value of key is par shortname
		z=nc.Dataset(zp_file) # could be outside loop

		# init grid stack
		xdim=out_xyz_dem.shape[0]
		sa_vec = np.zeros(xdim)
		#names=1:n

		for timestep in range(starti, endi):

			"""
			2D interpolated of surface firelds.
				Args:
					timestep: Timestep of interpolation as an interger (index)
					stations: pandas dataframe of input station csv file (id,lon,lat,ele)
					var: surface variarble eg "ssrd"
				Returns:
					t_sp: 2D interpolation of ERA surface field to stations points

				Example:
					surTa = ds.surTaPoint(0, mystations, 't2m')
			"""
			# read in data from variable 'varname'
			in_v = t[varname][timestep,:,:]#geopotential
			lat = t.variables['latitude'][:]
			lat=lat[::-1]
			lon = t.variables['longitude'][:]

			# 2d interpolation
			f_sa = RegularGridInterpolator((lat,lon), in_v, 'linear', bounds_error=False)
			out_xy = np.asarray([lons,lats]).T
			sa_t = f_sa(out_xy) 

			# stack timepoint to existing
			sa_vec=np.column_stack((sa_vec,sa_t))

		# drop init row
		sa_vec =sa_vec[:,1:]

		# Add to gsob
		setattr( gsob, varname, sa_vec)



	print("Made a sob")
	#===============================================================================
	# Conversions
	#===============================================================================


	""" convert tp from m/timestep (total accumulation over timestep) to rate in mm/h 

				Args:
					step: timstep in seconds (era5=3600, ensemble=10800)

				Note: both EDA (ensemble 3h) and HRES (1h) are accumulated over the timestep
				and therefore treated here the same.
			https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
	"""
	tphrm = gsob.tp #/step*60*60 # convert metres per timestep (in secs) -> m/hour 
	gsob.pmmhr = tphrm	*1000 # m/hour-> mm/hour

	""" Convert SWin from accumulated quantities in J/m2 to 
	instantaneous W/m2 see: 
	https://confluence.ecmwf.int/pages/viewpage.action?pageId=104241513

	Args:
		step: timstep in seconds (era5=3600, ensemble=10800)

	Note: both EDA (ensemble 3h) and HRES (1h) are accumulated over the timestep
	and therefore treated here the same ie step=3600s (1h)
	https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
	"""
	gsob.strd = gsob.strd/3600  
	gsob.ssrd = gsob.ssrd/3600
	gsob.tisr = gsob.tisr/3600 

	gsob.gridEle=gsob.z[:,0]/g
	gtob.prate=gsob.pmmhr #*pcf # mm/hour
	gtob.psum=gsob.tp*1000*stephr #pcf mm/timestep
	print("conversions done")
	#===============================================================================
	# TopoSCALE
	#=============================================================================


	#===============================================================================
	# Precip
	#===============================================================================

	'''
	Args:
	fineEle:ele vector from station dataframe
	sob: contains gridEle, tp and dtime
	'''
	# convert TP to mm/hr

	# lookups = {
	# 	   1:0.35,
	# 	   2:0.35,
	# 	   3:0.35,
	# 	   4:0.3,
	# 	   5:0.25,
	# 	   6:0.2,
	# 	   7:0.2,
	# 	   8:0.2,
	# 	   9:0.2,
	# 	   10:0.25,
	# 	   11:0.3,
	# 	   12:0.35
	# }

	# # Precipitation lapse rate, varies by month (Liston and Elder, 2006).
	# pfis = gsob.dtime.month.map(lookups)
	# pfis = np.repeat(pfis.values[:,None], lp.ele.size, axis=1)

	# dz=(lp.ele-gsob.gridEle)/1e3  # Elevation difference in kilometers between the fine and coarse surface.

		   
	# pcf=(1+pfis.T*dz[:,None])/(1-pfis.T*dz[:,None])# Precipitation correction factor.
	#Pf=sob.pmmhr.T*lp



	#===============================================================================
	# Longwave
	#===============================================================================




	"""Convert to RH (should be in a function). Following MG Lawrence 
	DOI 10.1175/BAMS-86-2-225 """
	A1=7.625 
	B1=243.04 
	C1=610.94
	tc=gsob.t2m-273.15
	tdc=gsob.d2m-273.15
	tf=gtob.t-273.15 # fout.T
	c=(A1*tc)/(B1+tc)
	RHc=100*np.exp((tdc*A1-tdc*c-B1*c)/(B1+tdc)) # Inverting eq. 8 in Lawrence.

	""" Calculate saturation vapor pressure at grid and "subgrid" [also
	through function] using the Magnus formula."""

	svpf=C1*np.exp(A1*tf/(B1+tf))
	svpc=C1*np.exp(A1*tc/(B1+tc))

	"""
	Calculate the vapor pressure at grid (c) and subgrid (f).
	"""
	vpf=gtob.r*svpf/1e2 # RHf
	vpc=RHc*svpc/1e2

	"""
	Use the vapor pressure and temperature to calculate clear sky
	# emssivity at grid and subgrid. [also function]
	Konzelmann et al. 1994
	Ta in kelvin

	"""
	x1=0.43 
	x2=5.7
	cef=0.23+x1*(vpf/gtob.t)**(1/x2) #Pretty sure the use of Kelvin is correct.
	cec=0.23+x1*(vpc/gsob.t2m)**(1/x2)

	"""Diagnose the all sky emissivity at grid."""
	sbc=5.67e-8
	aec=gsob.strd/(sbc*gsob.t2m**4)
	# need to constrain to 1 as original code?

	""" 
	Calculate the "cloud" emissivity at grid, assume this is the same at
		subgrid.
	"""
	deltae=aec-cec

	""" 
	Use the former cloud emissivity to compute the all sky emissivity at 
	subgrid. 
	"""
	aef=cef+deltae
	gtob.lwin=aef*sbc*gtob.t**4


	print("Lwin done")
	#===============================================================================
	# make a pob - required as input to swin routine. This object is data on each 
	# pressure level interpolated to the fine grid ie has dimensions xy (finegrid) x plev
	#===============================================================================

	f = nc.Dataset(zp_file)
	lev = f.variables['level'][:]
	var='t'
	#	ds = rc.t3d( pl=plevfile, dem =demfile)
	#	out_xyz_dem, lats, lons, shape= ds.demGrid()
	xdim=lev.shape[0]
	ydim=out_xyz_dem.shape[0]  
	t_interp_out = np.zeros((xdim,ydim))
	z_interp_out = np.zeros((xdim,ydim))
		
	# for timestep in range(starti,endi):
	# 	gridT,gridZ,gridLat,gridLon=ds.gridValue(var,timestep)
	# 	t_interp, z_interp = ds.inLevelInterp(gridT,gridZ, gridLat,gridLon,out_xyz_dem)
	# 	t_interp_out = np.dstack((t_interp_out, t_interp))
	# 	z_interp_out = np.dstack((z_interp_out, z_interp))


	#======================= all this can be done in first instance need to gen t/z_interp_out additionally
	tfile = plevDict.keys()[plevDict.values().index('t')] 
	t= nc.Dataset(tfile) # key is filename
	z=nc.Dataset(zp_file) # could be outside loop

	for timestep in range(starti, endi):

		gridT = t.variables['t'][timestep,:,:,:]
		gridZ = z.variables['z'][timestep,:,:,:]
		gridLat = t['latitude'][:]
		gridLat = gridLat[::-1] # reverse to deal with ERA5 order
		gridLon = t['longitude'][:]
		gridLev = t.variables['level'][::-1] 
		#return gridT,gridZ,gridLat,gridLon


		shape = gridT.shape
		#create array to hold interpolation resultes
		t_interp = np.zeros([shape[0], len(out_xyz_dem)])
		z_interp = np.zeros([shape[0], len(out_xyz_dem)]) # HOW MANY TIMES DO WE REALLY NEED TO COMPUTE THIS?

		#temperatue and elevation interpolation 2d
		for i in range(shape[0]):
			ft = RegularGridInterpolator((gridLat,gridLon), 
		                                  gridT[i,:,:], 'linear', bounds_error=False)
			fz = RegularGridInterpolator((gridLat,gridLon), 
		                                  gridZ[i,:,:], 'linear', bounds_error=False)
			t_interp[i,:] = ft(out_xyz_dem[:,:2] )#temperature
			z_interp[i,:] = fz(out_xyz_dem[:,:2])#elevation

		t_interp_out = np.dstack((t_interp_out, t_interp))
		z_interp_out = np.dstack((z_interp_out, z_interp))

		# invert pressure levels
		#t_interp = t_interp[::-1,:]
		#z_interp = z_interp[::-1,:]


	# drop init blank layer
	tinterp =t_interp_out[:,:,1:]
	zinterp =z_interp_out[:,:,1:]

	gpob= hp.Bunch(t=tinterp,z=zinterp, levels=lev)

	# dem  = nc.Dataset(demfile)
	# lon = dem.variables['longitude'][:]
	# lat = dem.variables['latitude'][:]
	# demv = dem.variables['ele'][:]
	# demlon1 = dem.variables['longitude'][:]
	# demlat1 = dem.variables['latitude'][:]
	# demlon = np.tile(demlon1,demlat1.size)
	# demlat = np.repeat(demlat1,demlon1.size)
	# demv=np.reshape(demv,demv.size)

	# # why are these masked values generated?
	# demv =np.ma.filled(demv, fill_value=1)
	# tz = np.repeat(tz,demv.size)

	# stat = pd.DataFrame({	"ele":demv, 
	# 				"lon":demlon,
	# 				"lat":demlat,
	# 				"tz":tz					
	# 				})
	#===============================================================================
	# Compute Shortwave
	#===============================================================================	

	#def swin2D(pob,sob,tob, stat, dates): 
	'''
	main edit over standard function for points:
	- 3d tob,sob reduce to 2D (reshape)
	'''


	""" toposcale surface pressure using hypsometric equation - move to own 
	class """
	g=9.81
	R=287.05  # Gas constant for dry air.
	#ztemp = pob.z # geopotential height b = np.transpose(a, (2, 0, 1))
	ztemp = np.transpose(gpob.z, (2, 0, 1)) # pob is originally ordered levels,stations/cells,time
	#Ttemp = pob.t
	Ttemp = np.transpose(gpob.t, (2, 0, 1))# pob is originally ordered levels,stations/cells,time
	statz = np.array(lp.ele)*g
	#dz=ztemp.transpose()-statz[None,:,None] # transpose not needed but now consistent with CGC surface pressure equations
	dz=ztemp-statz # dimensions of dz : time, levels, stations

	# set all levels below surface to very big number so they canot be found by min
	newdz=dz
	newdz[newdz<0]=999999

	psf =np.zeros( (gsob.dtime.size , statz.shape[0]) )

	# reshape tob.t here
	#gtob.tr=gtob.t.reshape(gtob.t.shape[0]*gtob.t.shape[1], gtob.t.shape[2], order='F')
	#tob.trT=tob.tr.T # transpose to get right order


	# loop through timesteps
	for i in range(0,gsob.dtime.size ):
		
		# find overlying layer
		thisp = dz[i,:,:]==np.min(newdz[i,:,:],axis=0) # thisp is a booleen matrix of levels x stations with true indicating overlying plevel over station surface ele

		# flatten to 1 dimension order='Fortran' or row major
		thispVec =thisp.reshape(thisp.size,order='F')
		TtempVec = Ttemp.reshape(Ttemp.shape[0], Ttemp.shape[1]*Ttemp.shape[2], order='F')
		ztempVec = ztemp.reshape(ztemp.shape[0], ztemp.shape[1]*ztemp.shape[2], order='F')
		
		# booleen indexing to find temp and geopotential that correspond to lowesest overlying layer
		T1=TtempVec[i,thispVec]
		z1=ztempVec[i,thispVec]


		p1=np.tile(gpob.levels[::-1],statz.shape[0])[thispVec]*1e2 #Convert to Pa. Reverse levels to ensure low ele (hig pressure) to high elel (low pressure)
		Tbar=np.mean([T1, gtob.t[:, i]],axis=0) # temperature midway between surface (toposcale T) and loweset overlying level (T1)
		""" Hypsometric equation.""" #P1 is above surface is this correct? Yes!
		psf[i,:]=(p1*np.exp((z1-statz)*(g/(Tbar*R)))) # exponent is positive ie increases pressure as surface is lower than pressure level


	""" Maybe follow Dubayah's approach (as in Rittger and Girotto) instead
	for the shortwave downscaling, other than terrain effects. """

	""" Height of the "grid" (coarse scale)"""
	Zc=gsob.z.T #.reshape(gsob.z.shape[0]*gsob.z.shape[1], gsob.z.shape[2]).T # reshape and transpose to remove dimension and make broadcastable

	""" toa """
	SWtoa = gsob.tisr#.reshape(gsob.tisr.shape[0]*gsob.tisr.shape[1], gsob.tisr.shape[2]).T  

	""" Downwelling shortwave flux of the "grid" using nearest neighbor."""
	SWc=gsob.ssrd#.reshape(sob.ssrd.shape[0]*sob.ssrd.shape[1], sob.ssrd.shape[2]).T 

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

	ztemp = np.transpose(gpob.z, (0, 2, 1))
	Ttemp = np.transpose(gpob.t, (0, 2, 1))
	dz=ztemp-Zc # dimensions of dz : levels, time, stations

	# set all levels below surface to very big number so they canot be found by min
	newdz=dz
	newdz[newdz<0]=999999

	psc =np.zeros( (gsob.dtime.size, statz.shape[0]) )
	for i in range(0,gsob.dtime.size):

		#thisp.append(np.argmin(dz[:,i][dz[:,i]>0]))
		# find overlying layer
		thisp = dz[:,i,:]==np.min(newdz[:,i,:],axis=0) # thisp is a booleen matrix of levels x stations with true indicating overlying plevel over station surface ele !! time index in middle this time!!!
		z0 = Zc[i,:]
		T0 = gsob.t2m.T[i,:]

		# flatten to 1 dimension order='Fortran' or row major
		thispVec =thisp.reshape(thisp.size,order='F')
		TtempVec = Ttemp.reshape(Ttemp.shape[1], Ttemp.shape[0]*Ttemp.shape[2], order='F') # !! order of permutations is different from pressure at finegrid routine (time is middle dimension)
		ztempVec = ztemp.reshape(ztemp.shape[1], ztemp.shape[0]*ztemp.shape[2], order='F')# !! order of permutations is different from pressure at finegrid routine (time is middle dimension)
		
		# booleen indexing to find temp and geopotential that correspond to lowesest overlying layer
		T1=TtempVec[i,thispVec]
		z1=ztempVec[i,thispVec]

		p1=np.tile(gpob.levels[::-1],statz.shape[0])[thispVec]*1e2 #Convert to Pa.
		Tbar=np.mean([T0, T1],axis=0)
		""" Hypsometric equation."""
		psc[i,:] = (p1*np.exp((z1-z0)*(g/(Tbar*R))))



	"""compute julian dates"""
	jd= sg.to_jd(gsob.dtime)

	"""
	Calculates a unit vector in the direction of the sun from the observer 
	position.
	"""
	svx,svy,svz =	sg.sunvectorMD(jd, lp.lat, lp.lon, tz,lp.ele.size,gsob.dtime.size )
	sp=sg.sunposMD(svx,svy,svz)


	# Cosine of the zenith angle.
	#sp.zen=sp.zen*(sp.zen>0) # Sun might be below the horizon.
	muz=np.cos(sp.zen) 
	muz = muz
	# NB! psc must be in Pa (NOT hPA!).

	# Calculate the "broadband" absorption coefficient. Elevation correction
	ka=(g*muz/(psc))*np.log(SWtoa.T/SWcdir.T)	
	ka = np.nan_to_num(ka) 

	# Now you can (finally) find the direct component at subgrid. 
	SWfdir=SWtoa.T*np.exp(-ka*psf/(g*muz))
	SWfdirCor=SWfdir #*dprod

	gtob.swin  =  SWfdiff+ SWfdirCor.T

	print("Swin done")


	#gtob.swin = swin2D(gpob,gsob,gtob, stat, dtime)
	#gtob.swin =gtob.swin.reshape(gtob.swin.shape[0], gsob.ssrd.shape[0], gsob.ssrd.shape[1])

	ntime = len(gsob.dtime)
	stephr =a.seconds/60/60
	rtime=np.array(range(len(gsob.dtime)))*stephr


	#===============================================================================
	# write results
	#===============================================================================	
	## conversuions
	T = np.single(gtob.t -273.15)
	gtob.r[gtob.r>100]=100
	RH = np.single(gtob.r)
	Abs_Hum = np.single( (6.112 * np.exp( (17.67 * T)/(T+243.5) ) * RH * 2.1674 )/(273.15+T) )
	ws = np.single( np.sqrt(gtob.u**2+gtob.v**2) )
	prate = np.single(gtob.prate)

	# Dictionary to loop over
	varDict={
	"t":T,
	"ws"   : ws,
	"shum" : Abs_Hum,  
	"swin" : gtob.swin,
	"lwin" : gtob.lwin, 
	"prate": prate 
	}

	for var in varDict:
		#open
		f = nc.Dataset(outDir+"/"+var+"_"+str(startIndex+1)+"_"+str(year)+".nc",'w', format='NETCDF4')

		#make dimensions
		f.createDimension('lon', len(lp.lon))
		f.createDimension('lat', len(lp.lat))
		f.createDimension('time', ntime)

		#make dimension variables
		longitude = f.createVariable('lon',    'f4',('lon',))
		latitude  = f.createVariable('lat',    'f4',('lat',))
		mytime = f.createVariable('time', 'i', ('time',))
		myvar = f.createVariable(var,    'f4',('lon','time'))

		#assign dimensions
		longitude = lp.lon
		latitude  = lp.lat
		mytime[:] = rtime
		myvar[:] = varDict[var]

		#metadata
		f.history = 'Created by toposcale on '
		mytime.units = 'hours since '+str(gsob.dtime[0])
		f.close()




	print("Toposcale complete!")





#===============================================================================
#	Calling Main
#===============================================================================
if __name__ == '__main__':
	coords = sys.argv[1]
	eraDir = sys.argv[2]
	outDir = sys.argv[3]
	start = sys.argv[4]
	end = sys.argv[5]

	main(coords, eraDir, outDir, start, end )
