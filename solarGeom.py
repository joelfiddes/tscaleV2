#solarGeom.py

# python implementation of corripio 2003 / R package 
import datetime
import numpy as np
import pandas as pd



def to_jd(dt):
	"""
	dt: python dattime object
	Converts a given datetime object (dt) to Julian date.
	Algorithm is copied from https://en.wikipedia.org/wiki/Julian_day
	All variable names are consistentdatetime.datetime.now() with the notation 
	on the wiki page.
	
	cite:https://github.com/dannyzed/julian/blob/master/julian/julian.py

	test : dt = datetime.datetime.now()

	Parameters
	----------
	fmt
	dt: datetime
		Datetime object to convert to MJD
	Returnsstokkang
	
	-------
	jd: float
	"""
	a = np.floor((14-dt.month)/12)
	y = dt.year + 4800 - a
	m = dt.month + 12*a - 3

	jdn = dt.day + np.floor((153*m + 2)/5) + 365*y + np.floor(y/4) - np.floor(y/100) + np.floor(y/400) - 32045
	jd = jdn + (dt.hour - 12.) / 24. + dt.minute / 1440. + dt.second / 86400. + dt.microsecond / 86400000000.
	return (jd)



def eqtime (jd): 
	""" 
	Computes the equation of time for a given Julian Day.

	jd: Julian Day and decimal fraction.
	cat("USAGE: eqtime(jd)\n") 
	"""
	jdc = (jd - 2451545.)/36525.
	sec = 21.448 - jdc * (46.815 + jdc * (0.00059 - jdc * (0.001813)))
	e0 = 23. + (26. + (sec/60.))/60.
	ecc = 0.016708634 - jdc * (4.2037e-05 + 1.267e-07 * jdc)
	oblcorr = e0 + 0.00256 * np.cos(np.radians(125.04 - 1934.136 * jdc))
	y = (np.tan(np.radians(oblcorr)/2.))**2.
	l0 = 280.46646 + jdc * (36000.76983 + jdc * (0.0003032))
	l0 = (l0 - 360. * (l0//360.))%360.
	rl0 = np.radians(l0)
	gmas = 357.52911 + jdc * (35999.05029 - 0.0001537 * jdc)
	gmas = np.radians(gmas)
	EqTime = y * np.sin(2. * rl0) - 2. * ecc * np.sin(gmas) + 4. * ecc * y * np.sin(gmas) * np.cos(2. * rl0) - 0.5 * y**2. * np.sin(4. * rl0) - 1.25 * ecc**2. * np.sin(2. * gmas)
	return(np.degrees(EqTime) * 4.)

def hourangle(jd, longitude, timezone): 
	""" 
	Hour angle, internal function for solar position.

	jd: Julian Day and decimal fraction.
	latitude: Latitude of observer in degrees and decimal fraction.
	longitude: Longitude of observer in degrees and decimal fraction.

	cat("USAGE: hourangle(jd,longitude,timezone)\n julian day, degrees, hours. 
	Return radians \n") 
	"""
	hour = ((jd - np.floor(jd)) * 24. + 12.)%24.
	myeqtime = eqtime(jd)
	stndmeridian = timezone * 15.
	deltalontime = longitude - stndmeridian
	deltalontime = deltalontime * 24./360.
	omegar = np.pi * (((hour + deltalontime + myeqtime/60.)/12.) - 1.)
	return(omegar)

def declination(jd): 
	""" 
	Computes the declination of the Sun for a given Julian Day.

	jd: Julian Day and decimal fraction.
	cat("USAGE: declination(jd) \n")
	"""
	jdc = (jd - 2451545)/36525
	sec = 21.448 - jdc * (46.815 + jdc * (0.00059 - jdc * (0.001813)))
	e0 = 23. + (26. + (sec/60.))/60.
	oblcorr = e0 + 0.00256 * np.cos(np.radians(125.04 - 1934.136 * jdc))
	l0 = 280.46646 + jdc * (36000.76983 + jdc * (0.0003032))
	l0 = (l0 - 360. * (l0//360.))%360.
	gmas = 357.52911 + jdc * (35999.05029 - 0.0001537 * jdc)
	gmas = np.radians(gmas)
	seqcent = np.sin(gmas) * (1.914602 - jdc * (0.004817 + 1.4e-05 * jdc)) + np.sin(2 * gmas) * (0.019993 - 0.000101 * jdc) + np.sin(3 * gmas) * 0.000289
	suntl = l0 + seqcent
	sal = suntl - 0.00569 - 0.00478 * np.sin(np.radians(125.04 - 1934.136 * jdc))
	delta = np.arcsin(np.sin(np.radians(oblcorr)) * np.sin(np.radians(sal)))
	return(np.degrees(delta))

def sunvector (jd, latitude, longitude, timezone): 
	"""
	Calculates a unit vector in the direction of the sun from the observer position

	jd: Julian Day and decimal fraction.
	latitude: Latitude of observer in degrees and decimal fraction.
	longitude: Longitude of observer in degrees and decimal fraction.
	timezone: Time zone in hours, west is negative.
	# cat("USAGE: sunvector(jd,latitude,longitude,timezone)\n values in jd, degrees, hours\n")
	"""
	omegar = hourangle(jd, longitude, timezone)
	deltar = np.radians(declination(jd))
	lambdar = np.radians(latitude)
	svx = -np.sin(omegar) * np.cos(deltar)
	svy = np.sin(lambdar) * np.cos(omegar) * np.cos(deltar) - np.cos(lambdar) * np.sin(deltar)
	svz = np.cos(lambdar) * np.cos(omegar) * np.cos(deltar) + np.sin(lambdar) * np.sin(deltar)
	# class Bunch:
	# 	def __init__(self, **kwds):
	# 		self.__dict__.update(kwds)
	# sv= Bunch(x=svx,y=svy,z=svz)   
	sv = np.c_[svx,svy,svz]
	return(sv)

def normalvector (slope, aspect): 
	"""
	Calculates a unit vector normal to a surface defined by slope inclination 
	and slope orientation.

	slope: slope of position in degrees
	aspect: aspect of position in degrees
	print("USAGE: normalvector(slope,aspect) \n")
	"""
	sloper = np.radians(slope)
	aspectr = np.radians(aspect)
	nvx = np.sin(aspectr) * np.sin(sloper)
	nvy = -np.cos(aspectr) * np.sin(sloper)
	nvz = np.cos(sloper)
	# class Bunch:
	# 	def __init__(self, **kwds):
	# 		self.__dict__.update(kwds)
	# nv= Bunch(x=nvx,y=nvy,z=nvz) 
	nv = np.c_[nvx,nvy,nvz]
	return(nv)

def sunpos(sunv): 
	"""
	Returns a matrix of azimuth and zenith angles of the sun given the unit 
	vectors from the observer to
	the direction of the sun. Plus sun elevation.

	sunv: sunvector
	#print("USAGE: sunpos(sunvector) 3D vector")
	"""
	azimuth = np.degrees(np.pi - np.arctan2(sunv[:,0], sunv[:,1]))
	zenith = np.degrees(np.arccos(sunv[:,2]))
	sunel = 90-zenith # sun elevation
	class Bunch:
		def __init__(self, **kwds):
			self.__dict__.update(kwds)
	sp=Bunch(azi=azimuth, zen=zenith, sel=sunel)
	return(sp)

