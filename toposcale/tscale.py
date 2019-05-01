"""tscale.py

Downscale forcing from atmospheric grid.

This module downscales all required fields, TA, RH, WS, WD, LW, SW based on
properties of a stat object that represents a station (Mode=POINT) or a cluster 
centroid (mode=TSUB) or a pixel (mode=GRID).

Example:

	Initialise a new tscale instance with::
        	$ p.plevels()

Attributes:

Todo:
    * Precipitation look up table hardcoded in method



"""


import numpy as np
import netCDF4 as nc
import solarGeom as sg
reload(sg)
class tscale(object):
	"""
	Interpolate pressure level vector from single cgc in the TopoSUB use case

	Args:

	Example:

	"""	

	def __init__(self, z):
		"""Example function with types documented in the docstring.

	    `PEP 484`_ type annotations are supported. If attribute, parameter, and
	    return types are annotated according to `PEP 484`_, they do not need to be
	    included in the docstring:

	    Args:
	        param1 (int): The first parameter.
	        param2 (str): The second parameter.

	    Returns:
	        bool: The return value. True for success, False otherwise.

	    .. _PEP 484:
	        https://www.python.org/dev/peps/pep-0484/

	    """

		self.g= 9.80665 #m s-2
		#self.statEle = statEle
		self.z=z
		self.R=287.05  # Gas constant for dry air.

		

	def tscale1D(self,dat, stat):
		"""Example function with types documented in the docstring.

	    `PEP 484`_ type annotations are supported. If attribute, parameter, and
	    return types are annotated according to `PEP 484`_, they do not need to be
	    included in the docstring:

	   	Interpolates vector of pressure level data to station elevation at 
		each timestep to create timeseries of downscaled values at station 
		elevation

	    Args:
	        param1 (int): The first parameter.
	        param2 (str): The second parameter.

	    Returns:
	        bool: The return value. True for success, False otherwise.

	    .. _PEP 484:
	        https://www.python.org/dev/peps/pep-0484/

	    """
		self.dat =dat


		self.ele=self.z/self.g
		self.interpVar=[]

		# this is a timestep loop
		for i in range(0,(self.z.shape)[0]): 
			
			y = self.dat[i,]
			x = self.ele[i,]

			# test if x (ele) is ascending (EDA), if not (HRES) reverse vectors
			# required for interp method https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.interp.html
			if not np.all(np.diff(x)>0):
				x =x[::-1]
				y =y[::-1]
			
			self.interpVar.append(np.interp(stat.ele,x, y ))
		self.interpVar = np.array(self.interpVar) # convert list to np array

	def addVar(self,varname,datInterp):
		""" Assign correct attribute name """
		setattr(self, varname, datInterp)
		
	# class interp3D(object): OR FUNCYTION HERE?
# 	"""
# 	Return interpolated cgc to point of interest

# 	Args:
# 		lon: POI lon
# 		lat: POI lat
# 		nc: data netcdf

# 	Example:
# 		lon=7.5
# 		lat=48.5
# 		nc= /home/joel/PLEV.nc

# 		cgcExtract(lon, lat, nc)
		
# 	""" 

# class lwin(object):
# 	""" Methods for downscaling longwave radiation """
	
# 	def __init__(self, statEle):
# 		self.statEle = statEle
	def wind (self, tob):
		""" methods for working with wind"""
		self.ws = np.sqrt(tob.u**2+tob.v**2)
		self.wd =   (180 / np.pi) * np.arctan(tob.u/tob.v) + np.where(tob.v>0,180,np.where(tob.u>0,360,0))

	def windCorRough(self, tob, pob, sob, stat ):
		"""corrects wind speed according to true rougness values"""
		target_height = 2
		blending_height = 40
		z0 = 1e-3#rough
		#av_wind = 0
		#wind_over_10 = 0
		#timesteps = 365*4
		#stationEle=2540
		#lon= 4#box index in dowwloaded era
		#lat= 2


		# surface
		#nc = nc_open(surface_file)
		#u_starvec = ncvar_get(nc, 'zust')
		u_starvec = sob.zust
		#nc = nc_open(surface_file)
		#Qhvec = ncvar_get(nc, 'ishf')
		Qhvec = sob.ishf
		#nc = nc_open(surface_file)
		#Qe = ncvar_get(nc, 'ie')
		Qe = sob.ie
		Qevec = (Qe*2.264e6)
		Tvec = tob.t		# T is air temperature (degC) downscaled to blending_height above station
		windvec= np.sqrt(tob.u**2 +tob.v**2) # winddownscaled to blending_height above station
		



		def	get_L_star (u_star, Qh, Qe, T):
			rho=1.293
			cp=1005
			L=2.26e6
			kappa=0.4
			g=9.81
			T = T #+273.15
			one_over_L_star = 1/(rho*cp*T/kappa/g*u_star**3/(Qh + 0.61*cp/L*T*Qe))
			return one_over_L_star

		def	psi_M(zeta1, zeta2):
			a=1
			b=2/3
			c=5
			d=0.35
			zeta1 = zeta1  + (zeta1>5)* (5-zeta1)
			zeta2 = zeta2 + (zeta2>5) * (5-zeta2)
			

			if (zeta1 <= 0):
				res=-2*np.arctan((1 - 16*zeta1)**(1/4)) + np.log( (1 + (1 - 16*zeta1)**(1/4))**2 * (1 + (1 - 16*zeta1)**0.5)/8) - (-2*np.arctan((1 - 16*zeta2)**(1/4)) + np.log( (1 + (1 - 16*zeta2)**(1/4))**2 * (1 + (1 - 16*zeta2)**0.5)/8))
			else:
				res = -b*(zeta1 - c/d) *np.exp(-d*zeta1) - a* zeta1 - (-b*(zeta2 -c/d) *np.exp(-d*zeta2) - a* zeta2)
			return res

		wind_downscaled=[]	
		for i in range(0,len(Tvec)):
			u_star = u_starvec[i]
			Qh = Qhvec[i]
			Qe = Qevec[i]

			T=Tvec[i]
			wind=windvec[i]
			one_over_L_star = get_L_star(u_star, Qh, Qe, T)

			valid_flag = ( (np.log(target_height / z0) - psi_M(target_height 
			* one_over_L_star, z0 * one_over_L_star)) 
			/(np.log(blending_height / z0) 
			- psi_M(blending_height * one_over_L_star, z0* one_over_L_star)))

			if (np.isnan(valid_flag)):
				valid_flag = 1
			if(valid_flag<1e-3):
				valid_flag = 1
			if(valid_flag>1):
				valid_flag = 1

			wind_downscaled.append(wind * valid_flag)

		self.ws = wind_downscaled	
			
			





	def lwin(self, sob,tob):
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
		aec[aec>1]<-1	
		# need to constrain to 1 as original code? I think so

		""" Calculate the "cloud" emissivity at grid, assume this is the same at
	 	subgrid."""
		deltae=aec-cec

		""" Use the former cloud emissivity to compute the all sky emissivity at 
		subgrid. """
		aef=cef+deltae
		self.LWf=aef*sbc*tob.t**4
		#fout.LW(:,:,n)=LWf


	def Rh2Wvp(tair,RH):
		''''compute water vapour pressure

		Args:
			RH=relative humidity (%: 0-100)
			tair=air temperature (kelvin)

		'''
		#constants
		es0=6.11 #reference saturation vapor pressure (es at a certain temp, usually 0 deg C)
		T0=273.15 #(273.15 Kelvin,  Kelvin = degree C +273.15)
		lv=2.5E6 #latent heat of vaporization of water (2.5 * 10^6 joules per kilogram)
		Rv=461.5 #gas constant for water vapor (461.5 joules* Kelvin / kilogram)

		
		es = es0 * np.exp( lv/Rv * (1/T0 - 1/tair))
		e=(RH*es)/100
		return(e)
	

	def relHumCalc(tair,tdew):
		'''convert ERA interim surface dewpoint temp to relative humidity based on d2m and t2m
			Args:
				tdew=dewpoint temp (kelvin)
				tair=air temperature (kelvin)
		'''
		#constants
		es0=6.11 #reference saturation vapor pressure (es at a certain temp, usually 0 deg C)
		T0=273.15 #(273.15 Kelvin,  Kelvin = degree C +273.15)
		lv=2.5E6 #latent heat of vaporization of water (2.5 * 10^6 joules per kilogram)
		Rv=461.5 #gas constant for water vapor (461.5 joules* Kelvin / kilogram)



		e  = es0 * np.exp( lv/Rv * (1/T0 - 1/tdew))
		es = es0 * np.exp( lv/Rv * (1/T0 - 1/tair))

		RH = e/es * (100) 

		return(RH)


	def lwin_joel(self, sob,tob):
		#All T in kelvin

		T_sub=tob.t # K
		RH=tob.r # 0-100
		lwGrid= sob.strd#Wm**-2
		tGrid= sob.t2m #T
		rhGrid= relHumCalc(sob.t2m, sob.d2m)# 0-100
		#calculate water vapour pressure at sub
		pv_sub=Rh2Wvp(tair=T_sub, RH=RH)

		#parameters
		x1=0.484 # 0.484 or 0.43 Gubler 2012
		x2=8 # 8 or 5.7 Gubler 2012
		sb=5.67*10**-8

		#calculate clear-sky emissivity subgrid according Konzelmann
		Ecl_sub=0.23+ x1*(pv_sub*100/T_sub)**(1/x2)

		#calculate clear-sky LWin
		lwcl=Ecl_sub*sb*T_sub**4

		#T grid in Kelvin
		T_grid=tGrid
		RH=rhGrid
		#calculate all-sky emissivity at grid from SB equation
		Eas_grid=lwGrid/(sb*T_grid**4)	
		#constrain Eas_grid to max. 1
		Eas_grid[Eas_grid>1]<-1	

		#calculate water vapour pressure at grid
		pv_grid=Rh2Wvp(tair=T_grid,RH=RH )
		#calculate clear-sky emissivity grid according Konzelmann
		Ecl_grid=0.23+ x1*((pv_grid*100)/T_grid)**(1/x2)	#Ecl_grid

		#calculate cloud emissivivity at grid
		deltaE=Eas_grid-Ecl_grid

		#calculate all-sky emissivity subgrid 
		Eas_sub=Ecl_sub+deltaE				#Eas_sub
		Eas_sub[Eas_sub>1]<-1	

		#calculate lwin subgrid
		lwsub=Eas_sub*sb*T_sub**4

		return(lwsub)
		

	def swin(self,pob,sob,tob, stat, dates):
		
		""" toposcale surface pressure using hypsometric equation - move to own 
		class """

		ztemp = pob.z
		Ttemp = pob.t
		statz = stat.ele*self.g
		dz=ztemp.transpose()-statz # transpose not needed but now consistent with CGC surface pressure equations

		self.psf=[]
		# loop through timesteps
		for i in range(0,dates.size):
			
			# 	# find overlying layer
			thisp = dz[:,i]==np.min(dz[:,i][dz[:,i]>0])

			# booleen indexing
			T1=Ttemp[i,thisp]
			z1=ztemp[i,thisp]
			p1=pob.levels[thisp]*1e2 #Convert to Pa.
			Tbar=np.mean([T1, tob.t[i]],axis=0)
			""" Hypsometric equation."""
			self.psf.append(p1*np.exp((z1-statz)*(self.g/(Tbar*self.R))))

		self.psf=np.array(self.psf).squeeze()

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
		self.kt=kt


		"""
		Calculate the diffuse fraction following the regression of Ruiz-Arias 2010 
		
		"""
		kd=0.952-1.041*np.exp(-1*np.exp(2.3-4.702*kt))
		self.kd = kd

		""" Use this to calculate the downwelling diffuse and direct shortwave radiation at grid. """
		SWcdiff=kd*SWc
		SWcdir=(1-kd)*SWc
		self.SWcdiff=SWcdiff
		self.SWcdir=SWcdir

		""" Use the above with the sky-view fraction to calculate the 
		downwelling diffuse shortwave radiation at subgrid. """
		self.SWfdiff=stat.svf*SWcdiff
		self.SWfdiff.set_fill_value(0)
		self.SWfdiff = self.SWfdiff.filled()

		""" Direct shortwave routine, modified from Joel. 
		Get surface pressure at "grid" (coarse scale). Can remove this
		part once surface pressure field is downloaded, or just check
		for existance. """

		ztemp = pob.z
		Ttemp = pob.t
		dz=ztemp.transpose()-sob.z

		self.psc=[]
		for i in range(0,dz.shape[1]):
		
			#thisp.append(np.argmin(dz[:,i][dz[:,i]>0]))
			thisp = dz[:,i]==np.min(dz[:,i][dz[:,i]>0])
			z0 = sob.z[i]
			T0 = sob.t2m[i]
			T1=Ttemp[i,thisp]
			z1=ztemp[i,thisp]
			p1=pob.levels[thisp]*1e2 #Convert to Pa.
			Tbar=np.mean([T0, T1],axis=0)
			""" Hypsometric equation."""
			self.psc.append(p1*np.exp((z1-z0)*(self.g/(Tbar*self.R))))

		self.psc=np.array(self.psc).squeeze()

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
		sunv=sg.sunvector(jd=jd, latitude=stat.lat, longitude=stat.lon, timezone=stat.tz)

		"""
		Computes azimuth , zenith  and sun elevation 
		for each timestamp
		"""
		sp=sg.sunpos(sunv)
		self.sp=sp

		# Cosine of the zenith angle.
		sp.zen=sp.zen
		#sp.zen=sp.zen*(sp.zen>0) # Sun might be below the horizon.
		muz=np.cos(sp.zen) 
		self.muz = muz
		# NB! psc must be in Pa (NOT hPA!).
		#if np.max(psc<1.5e3): # Obviously not in Pa
			#psc=psc*1e2
	   
		
		# Calculate the "broadband" absorption coefficient. Elevation correction
		# from Kris
		ka=(self.g*muz/(self.psc))*np.log(SWtoa/SWcdir)	
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
		self.SWfdir=SWtoa*np.exp(-ka*self.psf/(self.g*muz))

		""" Then perform the terrain correction. [Corripio 2003 / rpackage insol port]."""

		"""compute mean horizon elevation - why negative hor.el possible??? """
		horel=(((np.arccos(np.sqrt(stat.svf))*180)/np.pi)*2)-stat.slp
		if horel < 0:
			horel = 0 
		self.meanhorel = horel

		"""
		normal vector - Calculates a unit vector normal to a surface defined by 
		slope inclination and slope orientation.
		"""
		nv = sg.normalvector(slope=stat.slp, aspect=stat.asp)

		"""
		Method 1: Computes the intensity according to the position of the sun (sunv) and 
		dotproduct normal vector to slope.
		From corripio r package
		"""
		dotprod=np.dot(sunv ,np.transpose(nv)) 
		dprod = dotprod.squeeze()
		dprod[dprod<0] = 0 #negative indicates selfshading
		self.dprod = dprod

		"""Method 2: Illumination angles. Dozier"""
		saz=sp.azi
		cosis=muz*np.cos(stat.slp)+np.sin(sp.zen)*np.sin(stat.slp)*np.cos(sp.azi-stat.asp)# cosine of illumination angle at subgrid.
		cosic=muz # cosine of illumination angle at grid (slope=0).

		"""
		SUN ELEVATION below hor.el set to 0 - binary mask
		"""
		selMask = sp.sel
		selMask[selMask<horel]=0
		selMask[selMask>0]=1
		self.selMask = selMask

		"""
		derive incident radiation on slope accounting for self shading and cast 
		shadow and solar geometry
		BOTH formulations seem to be broken
		"""
		#self.SWfdirCor=selMask*(cosis/cosic)*self.SWfdir
		self.SWfdirCor=selMask*dprod*self.SWfdir
	   
		self.SWfglob = self. SWfdiff+ self.SWfdirCor

		""" 
		Missing components
		- terrain reflection
		"""

	def precip(self, sob, stat):

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
		dz=(stat.ele-sob.gridEle)/1e3     # Elevation difference in kilometers between the fine and coarse surface.
		lp=(1+pfis*dz)/(1-pfis*dz)# Precipitation correction factor.
		
		self.prate=sob.pmmhr*lp # mm/hour
		self.psum=sob.tp*1000*lp # mm/timestep
