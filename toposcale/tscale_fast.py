# tscale3D

"""
python2

This module starts with list of lon and lat and generates forcing fields.
requires that era5 is downloaded.

Example:


Vars:


Details:


"""

import tscale as ts
import helper as hp
import subprocess
import glob
import pandas as pd
#import tscale2d_era5_src as t2d
#import tscale3d_era5_src as t3d
import recapp_era5 as rc
import solarGeom as sg
from tqdm import tqdm
import netCDF4 as nc
import numpy as np
import os
from joblib import Parallel, delayed
import xarray as xr
from configobj import ConfigObj
import tscale_lib as tlib
# all parameters from tscale config ob


wd="/home/joel/sim/qmap/ch_tmapp"
config = ConfigObj(wd + "/config.ini")


wdir="/home/joel/mnt/myserver/sim/tamara"
pointsShp=wdir+"/coords.shp"

# now we just do everything that has been dowbnloaded
#startDate = "1979-02-01"
#endDate = "2019-12-31"

# not implemented yet
#windCor="FALSE"
#plapse="FALSE"

tmappSrc="/home/joel/src/tmapp2/"
tscaleSrc="/home/joel/src/tscaleV2/"
demRes=3
pointsBuffer = 0.08
svfSectors = 8
svfMaxDist = 5000
#dataset='era5'
g=9.81
num_cores=4



#===============================================================================
# KODE
#===============================================================================
outdir = wdir+"/out/"
if not os.path.exists(outdir):
	os.makedirs(outdir)

#print("Downloading DEM")
tlib.download_sparse_dem(tmappSrc, wdir, pointsShp, demRes, pointsBuffer)
tlib.setup_sim_dirs(tmappSrc,wdir)

# now loop over pointd
homes = tqdm(sorted(glob.glob(wdir+"/sim/*")))

# for home in tqdm(homes):
print("computing terrain and listpoints")
# 	compute_terrain(tmappSrc, home, str(svfSectors), str(svfMaxDist))
# 	make_listpoints(tmappSrc, home, pointsShp)
# 	#lp = pd.read_csv(home + "/listpoints.txt") # this shouple be all points
Parallel(n_jobs=int(num_cores))(delayed(tlib.compute_terrain)(tmappSrc,home, str(svfSectors), str(svfMaxDist)) for home in homes)
Parallel(n_jobs=int(num_cores))(delayed(tlib.make_listpoints)(tmappSrc,home, pointsShp) for home in homes)

# concat all listpoint files
filenames = sorted(glob.glob(wdir + "/*/*/listpoints.txt"))

dfs = []
for filename in filenames:
	dfs.append(pd.read_csv(filename))

# Concatenate all lp data into one DataFrame
lp = pd.concat(dfs, ignore_index=True)

# running tscale
print("Running tscale")
mylist = glob.glob(wdir+'/forcing/SURF_*')
mymonths =tqdm(sorted([i.split('_', 1)[1] for i in mylist]))
Parallel(n_jobs=int(num_cores))(delayed(tlib.tscale3dmain)(mymonth) for mymonth in mymonths)

# agregate results
print("Aggregate results")
for id in tqdm(lp.id):

	filenames =	sorted(glob.glob(wdir+"/out/meteo"+"c"+str(id)+"*"))


	dfs = []
	for filename in filenames:
		dfs.append(pd.read_csv(filename))
	# Concatenate all data into one DataFrame
	met_agg = pd.concat(dfs, ignore_index=True)
	fileout=wdir+"/out/allmeteo"+"c"+str(id)+".csv" # convert member index back to 1-10

	met_agg.to_csv(path_or_buf=fileout ,na_rep=-999,float_format='%.3f')

	# clean up
	#for file in filenames:
		#os.remove(file)



# validate /plot


# generate listpoint files
print("Generating report")
filenames = sorted(glob.glob(wdir + "/*/*/listpoints.txt"))

dfs = []
for filename in filenames:
	dfs.append(pd.read_csv(filename))

# Concatenate all lp data into one DataFrame
lp = pd.concat(dfs, ignore_index=True)

# round values for pretty printing
lp=lp.iloc[:,0:9]
decimals = pd.Series([0, 0,0,0,2,4,4,2,0], index=lp.columns )
lp =lp.round(decimals)
mydtype = pd.Series(['int32', 'int32','int32','int32','float','float','float','float','int32'], index=lp.columns )
lp.astype(mydtype).dtypes
decimals = pd.Series([0, 0,0,0,2,4,4,2,0], index=lp.columns )
lp =lp.round(decimals)

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
with PdfPages(wdir+'/tscale_report.pdf') as pdf:

	table = lp
	header = table.columns
	table = np.asarray(table)
	fig = plt.figure(figsize=(8, 6))
	ax = plt.Axes(fig, [0., 0., 1., 1.])
	ax.set_axis_off()
	fig.add_axes(ax)
	tab = plt.table(cellText=table, colLabels=header, cellLoc='center', loc='center')
	tab.auto_set_font_size(True)
	#tab.set_fontsize(8)
	#tab.scale(0.7, 2.5)
	pdf.savefig(fig)
	plt.close()




	filenames =	sorted(glob.glob(wdir+"/out/allmeteoc*"))


	meanTa=[]
	for file in filenames:
		df = pd.read_csv(file) 
		meanTa.append(df.TA.mean() )


	plt.scatter(lp.ele , meanTa)

	for i, txt in enumerate(lp.id):
		plt.annotate(txt, (lp.ele[i], meanTa[i]))
	plt.title("station ID plot")
	plt.xlabel('Elevation (m asl)')
	plt.ylabel('MAAT (K)')
	plt.hlines(273.15,lp.ele.min() ,lp.ele.max() )
	#plt.show() 
	pdf.savefig()  # saves the current figure into a pdf page
	plt.close()


	for file in filenames:

		df = pd.read_csv(file) 
		df.drop(df.columns[[0]], axis=1, inplace=True ) 
		df.set_index(df.iloc[:,0], inplace=True)  
					
		stationid = file.split("allmeteoc")[1].split(".csv")[0]  
		df.plot( subplots=True, title="Station ID:"+ stationid  )
		#plt.xlim(pd.Timestamp('1979-01-01'), pd.Timestamp('1980-01-01'))
		#plt.show() 
		pdf.savefig()  # saves the current figure into a pdf page
		plt.close()



print("tscale3D_fast complete!")
print("Total runtime=")