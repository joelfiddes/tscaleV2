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
print(num_cores)
wd=sys.argv[2]
startDate=sys.argv[3]
endDate=sys.argv[4]
dataset=sys.argv[5]
mymember=sys.argv[6]

# args to set


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




print(num_cores)
print("running tscale jobs: "+str(njobs))


# run all memeber = 1 first
if dataset=="EDA":
        Parallel(n_jobs=int(njobs))(delayed(tscale3D.main)(wdir=wd, mode='points', start=datelist[i], end=datelist[i+1], dataset=dataset, member=mymember) for i in range(0,num_cores))
if dataset=="HRES":
        Parallel(n_jobs=int(njobs))(delayed(tscale3D.main)(wdir=wd, mode='points', start=datelist[i], end=datelist[i+1], dataset=dataset) for i in range(0,njobs))

print("All cluster jobs complete!")


