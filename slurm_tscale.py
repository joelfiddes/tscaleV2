#!/usr/bin/env python3
import os
import sys
sys.path.append(os.getcwd())

from joblib import Parallel, delayed
import numpy as np
import time
import sys
import tscale3D

num_cores = int(sys.argv[1]) # Must input number of cores
wd=sys.argv[2]
startDate=sys.argv[3]
endDate=sys.argv[4]
# args to set
wd="/home/caduff/tscale/"
startDate ="1979-09-01"
endDate = "2018-09-01"

# get years
sy=startDate.split("-")[0]
ey = endDate.split("-")[0]

# get month
sm = startDate.split("-")[1]
em = endDate.split("-")[1]

# get day
sd = startDate.split("-")[2]
ed = endDate.split("-")[2]


a = range(int(sy),int(ey)+1)
datelist = [str(s) + "-09-01"  for s in a]

njobs = len(a)-1





print("running tscale jobs: "+str(njobs))


# run all memeber = 1 first
Parallel(n_jobs=int(num_cores))(delayed(tscale3D.main)(wdir=wd, mode='points', start=datelist[i], end=datelist[i+1]) for i in range(0,njobs))

print("All cluster jobs complete!")

