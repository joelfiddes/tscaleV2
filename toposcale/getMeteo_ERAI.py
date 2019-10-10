"""main

Get meteo forcing script.

This module controls download and preprocessing of meteo files. Currently 
focussed on ERA5 but could be expanded to include others

Example:

	python main.py

Dependency:
	- CDO
	- NCO
	- CDS API
	- see note below

Todo:
	eracat and eracat5d could be combined as robust NCO based approach that accepts 
	all dimensions - reduce depencdy by -1 (remove CDO)

	- add time step args (config)
	- add plevel args (config)
"""

import sys
import os
from configobj import ConfigObj
import logging
import eraiRetrievePLEVEL as plev
import eraiRetrieveSURFACE as surf
import subprocess

#====================================================================
#	Get domain extent and parameters
#====================================================================

latNorth=52
latSouth=43
lonWest=88
lonEast=120
startDate = '1979-01-01'
endDate = '2019-01-01'
eraDir='/home/joel/sim/mongolia/ERAI_data'
#====================================================================
#	Config setup
#====================================================================
#os.system("python writeConfig.py") # update config DONE IN run.sh file

#os.remove("logfile")
logging.basicConfig(level=logging.DEBUG, filename="getmeteolog", filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")

logging.info("Start run")
#====================================================================
#	Define paths
#====================================================================
#wd= config["main"]["wd"]
#eraDir = wd + "/forcing/"
#if not os.path.exists(eraDir):
#	os.makedirs(eraDir)



#====================================================================
#	Fetch forcing
#====================================================================
surf.retrieve_interim(  startDate,endDate,eraDir, latNorth,latSouth,lonEast,lonWest
	)



plev.retrieve_interim(  startDate,endDate,eraDir, latNorth,latSouth,lonEast,lonWest
	)


#fe.eraCat(eraDir, "SURF")
#fe.eraCat(eraDir, "PLEV")






