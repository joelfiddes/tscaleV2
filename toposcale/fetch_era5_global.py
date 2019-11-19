''''
Retrieve ecmwf data from new CDS api for global scale:

- retrieve 1 variable per year per file
- correct 6 hr precip by monthly totals
- all fields required to drive model
- TA, SWIN, LWIN, RH, WD, WS, TP, P (stability stuff)



'''
#!/usr/bin/env python
from datetime import datetime, timedelta
from collections import OrderedDict
import calendar
import sys
#from ecmwfapi import ECMWFDataServer
import cdsapi
#from dateutil.relativedelta import *
#from retrying import retry
import logging
import glob	
from joblib import Parallel, delayed 
import subprocess
import os
#import multiprocessing 
step = 6
plevels =["1000", "700", "500", "200"]
    # ===============================================================================
    # Log
    # ===============================================================================
    # logfile=wd+"/sim/"+ simdir+"/logfile"

def retrieve_era5_surf(startYear,endYear,eraDir, step):
	""" Sets up era5 surface retrieval.
	* Creates list of year/month pairs to iterate through. 
	* MARS retrievals are most efficient when subset by time. 
	* Identifies preexisting downloads if restarted. 
	* Calls api using parallel function.

	Args:
		product: "reanalysis" (HRES) or "ensemble_members" (EDA)
		startDate:
		endDate:
		eraDir: directory to write output

	Returns:
		Monthly era surface files.		 

	"""
	# product  = config["forcing"]["product"]
	# print(product)
	# startDate = config["main"]["startDate"]
	# endDate = config["main"]["endDate"]
	#grd =   config["era-interim"]["grid"]
	#dataset = config["forcing"]["dataset"]
	#grid=str(grd) + "/" + str(grd)
	logfile = eraDir + "/surf_logfile"
	if os.path.isfile(logfile):
		os.remove(logfile)
	logging.basicConfig(level=logging.DEBUG, filename=logfile, filemode="a+",format="%(asctime)-15s %(levelname)-8s %(message)s")

	num_cores = 4 #config['main']['num_cores']
	myvars=['geopotential', '2m_dewpoint_temperature', 'surface_thermal_radiation_downwards', 'surface_solar_radiation_downwards',
		'Total precipitation','2m_temperature', 'TOA incident solar radiation',
			'friction_velocity','instantaneous_moisture_flux','instantaneous_surface_sensible_heat_flux']
	for var in myvars:

	
		if (step == 1):
			time = [		'00:00','01:00','02:00',\
					'03:00','04:00','05:00',\
					'06:00','07:00','08:00',\
					'09:00','10:00','11:00',\
					'12:00','13:00','14:00',\
					'15:00','16:00','17:00',\
					'18:00','19:00','20:00',\
					'21:00','22:00','23:00']

		if (step == 3):
			time = ['00:00','03:00','06:00','09:00','12:00','15:00','18:00','21:00']

		if (step == 6):
			time = ['00:00','06:00','12:00','18:00']

		# download buffer of +/- 1 month to ensure all necessary timestamps are there for interpolations and consistency between plevel and surf
		dates = [str(startYear), str(endYear)]
		start = datetime.strptime(dates[0], "%Y")
		end = datetime.strptime(dates[1], "%Y")
		dateList = OrderedDict(((start + timedelta(_)).strftime(r"%Y"), None) for _ in xrange((end - start).days)).keys()


		requestDatesVec = []
		targetVec=[]
		yearVec=[]
		for date in dateList:	
			strsplit = date.split('-' )
			year =  int(date)
			
			if eraDir.endswith('/'):
				s = eraDir[:-1]

			target = s + "/SURF_"+var+"_%04d.nc" % (year)
			targetVec.append(target) 
			yearVec.append(year)

		# find files that already downloaded if any with exact matches (in case of restarts)
		# NEW pattern matching for global parameter/year files eg. SURF_2m_dewpoint_temperature_2008.nc
		dataExists = glob.glob(eraDir +"/SURF_"+var+"_????.nc")

		# list only files that dont exist
		targetVecNew = [x for x in targetVec if x not in dataExists]
		logging.info("ECWMF SURF data found:" )
		logging.info(dataExists)
		logging.info("Downloading SURF from ECWMF:")
		logging.info(targetVecNew)

		# Amend requestDatesVec
		index = [targetVec.index(x) for x in targetVecNew]
		yearVecNew  = [yearVec[i] for i in index]

		# https://zacharyst.com/2016/03/31/parallelize-a-multifunction-argument-in-python/	
		Parallel(n_jobs=int(num_cores))(delayed( era5_request_surf)(var , int(yearVecNew[i]), targetVecNew[i],time) for i in range(0,len(yearVecNew)))
		logging.info(var + " "+ str(startYear) + " to "+str(endYear) +" complete!")






def era5_request_surf(var, year, target, time):
	"""CDS surface api call"""
	c = cdsapi.Client()

	c.retrieve(
		'reanalysis-era5-single-levels',
		{
		'variable':[var
		],
		'product_type': 'reanalysis',
		'year':int(year),
		'month':[
			'01','02','03',
			'04','05','06',
			'07','08','09',
			'10','11','12'
		],
		'day':[
			'01','02','03',
			'04','05','06',
			'07','08','09',
			'10','11','12',
			'13','14','15',
			'16','17','18',
			'19','20','21',
			'22','23','24',
			'25','26','27',
			'28','29','30',
			'31'
		],
		'time':time,
		'format':'netcdf'

		},
		target)
	print(target+ " complete")

def retrieve_era5_plev(startYear,endYear,eraDir, step, plevels):
	""" Sets up era5 pressure level retrieval.
	* Creates list of year/month pairs to iterate through. 
	* MARS retrievals are most efficient when subset by time. 
	* Identifies preexisting downloads if restarted. 
	* Calls api using parallel function.

	Args:
		config: config object defining INI
		eraDir: directory to write output
	latN: north latitude of bbox
		latS: south latitude of bbox
		lonE: easterly lon of bbox
		lonW: westerly lon of bbox

	Returns:
		Monthly era pressure level files.		 

	"""
	#product  = config["forcing"]["product"]
	#startDate = config["main"]["startDate"]
	#endDate = config["main"]["endDate"]
	#grd =   config["era-interim"]["grid"]
	#dataset = config["forcing"]["dataset"]
	#grid=str(grd) + "/" + str(grd)
	logfile = eraDir + "/plev_logfile"
	if os.path.isfile(logfile):
		os.remove(logfile)
	logging.basicConfig(level=logging.DEBUG, filename=logfile, filemode="a+",format="%(asctime)-15s %(levelname)-8s %(message)s")

	num_cores = 4 #config['main']['num_cores']
	myvars=['geopotential','temperature','u_component_of_wind',
			'v_component_of_wind', 'relative_humidity']
	for var in myvars:

	
		if (step == 1):
			time = [		'00:00','01:00','02:00',\
					'03:00','04:00','05:00',\
					'06:00','07:00','08:00',\
					'09:00','10:00','11:00',\
					'12:00','13:00','14:00',\
					'15:00','16:00','17:00',\
					'18:00','19:00','20:00',\
					'21:00','22:00','23:00']

		if (step == 3):
			time = ['00:00','03:00','06:00','09:00','12:00','15:00','18:00','21:00']

		if (step == 6):
			time = ['00:00','06:00','12:00','18:00']

		# download buffer of +/- 1 month to ensure all necessary timestamps are there for interpolations and consistency between plevel and surf
		dates = [str(startYear), str(endYear)]
		start = datetime.strptime(dates[0], "%Y")
		end = datetime.strptime(dates[1], "%Y")
		dateList = OrderedDict(((start + timedelta(_)).strftime(r"%Y"), None) for _ in xrange((end - start).days)).keys()

		requestDatesVec = []
		targetVec=[]
		yearVec=[]
		for date in dateList:	
			strsplit = date.split('-' )
			year =  int(date)
			if eraDir.endswith('/'):
				s = eraDir[:-1]
			target = s + "/PLEV_"+var+"_%04d.nc" % (year)
			targetVec.append(target) 
			yearVec.append(year)

		# find files that already downloaded if any with exact matches (in case of restarts)
		# NEW pattern matching for global parameter/year files eg. SURF_2m_dewpoint_temperature_2008.nc
		dataExists = glob.glob(eraDir +"/PLEV_"+var+"_????.nc")

		# list only files that dont exist
		targetVecNew = [x for x in targetVec if x not in dataExists]
		logging.info("ECWMF PLEV data found:" )
		logging.info(dataExists)
		logging.info("Downloading PLEV from ECWMF:")
		logging.info(targetVecNew)

		# Amend requestDatesVec
		index = [targetVec.index(x) for x in targetVecNew]
		yearVecNew  = [yearVec[i] for i in index]

	# https://zacharyst.com/2016/03/31/parallelize-a-multifunction-argument-in-python/	
	Parallel(n_jobs=int(num_cores))(delayed( era5_request_plev)(var, int(yearVecNew[i]), targetVecNew[i], time,plevels) for i in range(0,len(yearVecNew)))
	logging.info(var + " "+str(startYear) + " to "+str(endYear) +" complete!")

#@retry(wait_random_min=10000, wait_random_max=20000)




def era5_request_plev(var, year, target, time,plevels):
	"""CDS plevel api call"""
	c = cdsapi.Client()

	c.retrieve(
		'reanalysis-era5-pressure-levels',
		{
		'product_type': "reanalysis",
		'format':'netcdf',
	
		'variable':[
			var
		],
		'pressure_level':plevels,
		'year':int(year),
		'month':[
			'01','02','03',
			'04','05','06',
			'07','08','09',
			'10','11','12'
		],
		'day':[
			'01','02','03',
			'04','05','06',
			'07','08','09',
			'10','11','12',
			'13','14','15',
			'16','17','18',
			'19','20','21',
			'22','23','24',
			'25','26','27',
			'28','29','30',
			'31'
		],
		'time':time,
		'format':'netcdf'
		},
		target)
	print(target+ " complete")


#if __name__ == "__main__":

#	config	= sys.argv[1]
#	eraDir	 = sys.argv[2]
#	latNorth	= str(float(sys.argv[3]))
#	latSouth	=  str(float(sys.argv[4]))
#	lonEast	 = str(float(sys.argv[5]))
#	lonWest	 = str(float(sys.argv[6]))

#	retrieve_interim(config,eraDir, latNorth,latSouth,lonEast,lonWest)



	
def retrieve_era5_tpmm( startYear,endYear,eraDir):
	"""
	retrive monthly tp means to correct 6 or 3h TP retrieval
	documented era5 here: https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation#ERA5datadocumentation-Monthlymeans
	values are mean daily for that month, weight by number of days in month eg annualtotal = sum(wfj_m*c(31,28,31,30,31,30,31,31,30,31,30,31))
	"""
	# eraDir='/home/joel/sim/imis/forcing/'	
	# latN=47.9
	# latS=45.75
	# lonW=5.7
	# lonE=10.6

	# startDate = '1979-09-01'                        # adjust to your requirements - as of 2017-07 only 2010-01-01 onwards is available
	# endDate = '2018-09-01'   
	param='total_precipitation'                     # adjust to your requirements
	months = [1,2,3,4,5,6,7,8,9,10,11,12] 
	
	dates = [str(startYear), str(endYear)]
	start = datetime.strptime(dates[0], "%Y")
	end = datetime.strptime(dates[1], "%Y")
	dateList = OrderedDict(((start + timedelta(_)).strftime(r"%Y"), None) for _ in xrange((end - start).days)).keys()

	target = eraDir + "tpmm.nc"

	
#	yearVec=[]
#	for date in dateList:	
#		
#		year =  int(date)
#		yearVec.append(year)
#	

#	myyears = list(sorted(set(yearVec)))

	c = cdsapi.Client()

	c.retrieve(
	    'reanalysis-era5-single-levels-monthly-means',
	    {
		'product_type':'reanalysis-monthly-means-of-daily-means',
		'variable': param,
		'year':	myyears,
		'month':[
		    '01','02','03',
		    '04','05','06',
		    '07','08','09',
		    '10','11','12'
		],
		'time':'00:00',
		'format':'netcdf'
	    },
	    target)
