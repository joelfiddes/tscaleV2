
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
		self.g= 9.80665 #m s-2
		#self.statEle = statEle
		self.z=z
		self.R=287.05  # Gas constant for dry air.

		

	def tscale1D(self,dat, stat):
		""" Interpolates vector of pressure level data to station elevation at 
		each timestep to create timeseries of downscaled values at station 
		elevation"""
		self.dat =dat
		self.ele=self.z/self.g
		self.interpVar=[]
		for i in range(0,(self.z.shape)[0]): 
			self.interpVar.append(np.interp(stat.ele,self.ele[i,], self.dat[i,] ))
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

#class wind(object)
	""" methods for working with wind"""

	def swin(self,pob,sob,tob, stat, dates):
		
		""" toposcale surface pressure using hypsometric equation - move to own 
		class """

		ztemp = pob.z
		Ttemp = pob.t
		statz = stat.ele*self.g
		dz=ztemp.transpose()-statz # transpose not needed but now consistent with CGC surface pressure equations

		self.psf=[]
		for i in range(0,dz.shape[1]):
		
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
