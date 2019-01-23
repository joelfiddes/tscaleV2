"""Example Google style docstrings.

This module demonstrates documentation as specified by the `Google Python
Style Guide`_. Docstrings may extend over multiple lines. Sections are created
with a section header and a colon followed by a block of indented text.

Example:
    Examples can be given using either the ``Example`` or ``Examples``
    sections. Sections support any reStructuredText formatting, including
    literal blocks::

        $ python example_google.py

Section breaks are created by resuming unindented text. Section breaks
are also implicitly created anytime a new section starts.

Attributes:
    module_level_variable1 (int): Module level variables may be documented in
        either the ``Attributes`` section of the module docstring, or in an
        inline docstring immediately following the variable.

        Either form is acceptable, but the two should not be mixed. Choose
        one convention to document module level variables and be consistent
        with it.

Todo:
    * For module TODOs
    * You have to also use ``sphinx.ext.todo`` extension

.. _Google Python Style Guide:
   http://google.github.io/styleguide/pyguide.html

"""

# TopoSCALE: Topograhic based downscaling of atmospheric datasets
#
# === DESCRIPTION ==============================================================

# === COPYRIGHT AND LICENCE ====================================================
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# === CONTRIBUTIONS ============================================================
#
# === NOTES ====================================================================

# === REQUIREMENTS =============================================================
#	* PLEV.nc - all pressure level variables for entire domain (single/multiple 
#		CGCs) [time X grid cell X  level X variable]
#	* SURF.nc - all surface variables for entire domain (single/multiple 
#		CGCs) [time X grid cell X variable]

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
		for i in range(0,(self.z.shape)[0]): 

			# first reverse vectors as elevation has to be ascending in order for interp to work
			y = self.dat[i,]
			y =y[::-1]
			x = self.ele[i,]
			x =x[::-1]

			# check that elevation is ascending, required for interp method https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.interp.html
			if np.all(np.diff(x) > 0) == False:
				print "ERROR!! elevation vector not ascending, interpolation method cannot run"
				break

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

		self.wscor = wind_downscaled	
			
			





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
		 # emssivity at grid and subgrid. [also function]"""
		x1=0.43 
		x2=5.7
		cef=0.23+x1*(vpf/tf)**(1/x2) #Pretty sure the use of Kelvin is correct.
		cec=0.23+x1*(vpc/sob.t2m)**(1/x2)

		"""Diagnose the all sky emissivity at grid."""
		sbc=5.67e-8
		aec=sob.strd/(sbc*sob.t2m**4)

		""" Calculate the "cloud" emissivity at grid, assume this is the same at
	 	subgrid."""
		deltae=aec-cec

		""" Use the former cloud emissivity to compute the all sky emissivity at 
		subgrid. """
		aef=cef+deltae
		self.LWf=aef*sbc*tob.t**4
		#fout.LW(:,:,n)=LWf



	def swin(self,pob,sob,tob, stat, dates):
		
		""" toposcale surface pressure using hypsometric equation - move to own 
		class """

		ztemp = pob.z
		Ttemp = pob.t
		statz = stat.ele*self.g
		dz=ztemp.transpose()-statz # transpose not needed but now consistent with CGC surface pressure equations

		self.psf=[]
		# loop through timesteps
		for i in range(0,dz.shape[1]):
			print i
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
		self. SWfdiff=stat.svf*SWcdiff

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
		sunv=sg.sunvector(jd=jd, latitude=stat.lat, longitude=stat.long, timezone=stat.tz)

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
		
		# Note this equation is obtained by inverting Beer's law, i.e. use
		#I_0=I_inf x exp[(-ka/mu) int_z0**inf rho dz]
		# Along with hydrostatic equation to convert to pressure coordinates then
		# solve for ka using p(z=inf)=0.
		
		
		# Now you can (finally) find the direct component at subgrid. 
		SWfdir=SWtoa*np.exp(-ka*self.psf/(self.g*muz))
		self.SWfdir = SWfdir

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
		Pf=sob.tp*lp
		self.TPf=Pf
