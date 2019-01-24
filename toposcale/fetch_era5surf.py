#!/usr/bin/env python
from datetime import datetime, timedelta
from collections import OrderedDict
import calendar
import sys
#from ecmwfapi import ECMWFDataServer
import cdsapi
from dateutil.relativedelta import *
from retrying import retry
import logging
import glob

def retrieve_interim(config,eraDir, latNorth,latSouth,lonEast,lonWest):
    """      
       A function to demonstrate how to iterate efficiently over several years and months etc    
       for a particular interim_request.     
       Change the variables below to adapt the iteration to your needs.
       You can use the variable "target" to organise the requested data in files as you wish.
       In the example below the data are organised in files per month. (eg "interim_daily_201510.grb")

       useful sources:

       Step and time
       https://software.ecmwf.int/wiki/pages/viewpage.action?pageId=56658233
    """
    startDate = config["main"]["startDate"]
    endDate = config["main"]["endDate"]
    grd =   config["era-interim"]["grid"]
    dataset = config["era-interim"]["dataset"]
    num_cores = config['geotop']['num_cores']
    grid=str(grd) + "/" + str(grd)
    bbox=(str(latNorth) + "/" + str(lonWest) + "/" + str(latSouth) + "/" + str(lonEast)) 


    # download buffer of +/- 1 month to ensure all necessary timestamps are there for interpolations and consistency between plevel and surf
    dates = [str(startDate), str(endDate)]
    start = datetime.strptime(dates[0], "%Y-%m-%d")
    end = datetime.strptime(dates[1], "%Y-%m-%d")
    start = start+relativedelta(months=-1)
    end = end+relativedelta(months=+1)
    dateList = OrderedDict(((start + timedelta(_)).strftime(r"%Y-%m"), None) for _ in xrange((end - start).days)).keys()

    print("Retrieving"+dataset+" dataset")
    print("Bbox = " + bbox)
    print("Grid = " + grd)
    print("Start date = " , dateList[0])
    print("End date = " , dateList[len(dateList)-1])
    print("cores used = " + num_cores)

    requestDatesVec = []
    targetVec=[]
    for date in dateList:    
        strsplit = date.split('-' )
        year =  int(strsplit[0])
        month = int(strsplit[1])   
        firstDate = "%04d%02d%02d" % (year, month, 1)
        numberOfDays = calendar.monthrange(year, month)[1]
        lastDate = "%04d%02d%02d" % (year, month, numberOfDays)
        target = eraDir + "/SURF_%04d%02d.nc" % (year, month)
        requestDates = (firstDate + "/TO/" + lastDate)
        
        requestDatesVec.append(requestDates)
        targetVec.append(target)  

    # find files that already downloaded if any with exact matches (in case of restarts)
    dataExists = glob.glob(eraDir +"/SURF_??????.nc")

    # list only files that dont exist
    targetVecNew = [x for x in targetVec if x not in dataExists]
    logging.info("ECWMF SURF data found:" )
    logging.info(dataExists)
    logging.info("Downloading SURF from ECWMF:")
    logging.info(targetVecNew)

    # Amend requestDatesVec
    index = [targetVec.index(x) for x in targetVecNew]
    requestDatesVecNew  = [requestDatesVec[i] for i in index]

    # https://zacharyst.com/2016/03/31/parallelize-a-multifunction-argument-in-python/
    from joblib import Parallel, delayed 
    import multiprocessing 

    Parallel(n_jobs=int(num_cores))(delayed(interim_request)(requestDatesVecNew[i], targetVecNew[i] , grid, bbox, dataset,timeVec, step, eraClass) for i in range(0,len(requestDatesVecNew)))


@retry(wait_random_min=10000, wait_random_max=20000)












def interim_request(requestDates, target, grid, bbox, dataset, time, step, eraClass):

	c = cdsapi.Client()

	c.retrieve(
	    'reanalysis-era5-single-levels',
	    {
		'variable':['geopotential', '2m_dewpoint_temperature', 'surface_thermal_radiation_downwards', 'surface_solar_radiation_downwards',
		'Total precipitation','2m_temperature', 'TOA incident solar radiation',
		    'friction_velocity','instantaneous_moisture_flux','instantaneous_surface_sensible_heat_flux'
		],
		'product_type':'reanalysis',
			"area": "47/9/46/10",
		'year':'2015',
		'month':[
		    '01'
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
		'time':[
		    '00:00','01:00','02:00',
		    '03:00','04:00','05:00',
		    '06:00','07:00','08:00',
		    '09:00','10:00','11:00',
		    '12:00','13:00','14:00',
		    '15:00','16:00','17:00',
		    '18:00','19:00','20:00',
		    '21:00','22:00','23:00'
		],
		'format':'netcdf'

	    },
	    'surface.nc')

if __name__ == "__main__":

    config    = sys.argv[1]
    eraDir     = sys.argv[2]
    latNorth    = str(float(sys.argv[3]))
    latSouth    =  str(float(sys.argv[4]))
    lonEast     = str(float(sys.argv[5]))
    lonWest     = str(float(sys.argv[6]))

    retrieve_interim(config,eraDir, latNorth,latSouth,lonEast,lonWest)
