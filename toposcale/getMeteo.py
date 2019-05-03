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
import fetch_era5 as fe

#====================================================================
#	Config setup
#====================================================================
#os.system("python writeConfig.py") # update config DONE IN run.sh file
config = ConfigObj(sys.argv[1])
#os.remove("logfile")
logging.basicConfig(level=logging.DEBUG, filename="getmeteolog", filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")

logging.info("Start run")
#====================================================================
#	Define paths
#====================================================================
wd= config["main"]["wd"]
eraDir = wd + "/forcing/"
if not os.path.exists(eraDir):
	os.makedirs(eraDir)
#====================================================================
#	Get domain extent and parameters
#====================================================================
lonW = config["main"]["lonW"]
lonE = config["main"]["lonE"]
latN = config["main"]["latN"]
latS = config["main"]["latS"]

# timsteps to retrieve (1,3,6)
step = config["forcing"]["step"]

# pressure levels to retriev
plevels = config["forcing"]["plevels"]
#====================================================================
#	Fetch forcing
#====================================================================
fe.retrieve_era5_surf(   
	config["forcing"]["product"] ,
	config["main"]["startDate"], 
	config["main"]["endDate"],
	eraDir, 
	latN,latS,lonE,lonW,step
	)

fe.retrieve_era5_tpmm(   
	config["main"]["startDate"], 
	config["main"]["endDate"],
	eraDir, 
	latN,latS,lonE,lonW
	)

fe.retrieve_era5_plev(   
	config["forcing"]["product"] ,
	config["main"]["startDate"], 
	config["main"]["endDate"],
	eraDir, 
	latN,latS,lonE,lonW,step,plevels
	)

''' concat monthly files using CDO'''
if (config["forcing"]["product"] == "reanlysis"):
	fe.eraCat(eraDir, "SURF")
	fe. eraCat(eraDir, "PLEV")

'''5d ensemble product requires NCO operators'''
if (config["forcing"]["product"] == "ensemble_members"):
	fe.eraCat5d(eraDir, "SURF")
	fe. eraCat5d(eraDir, "PLEV")


cmd = ["Rscript",  config['main']['tscale_root']+"tpCorrect.R",eraDir+ "tpmm.nc" , eraDir, "SURF.nc"]
subprocess.check_output(cmd)


