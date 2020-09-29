# script to generate cryogrid data drom native tscale meteo format
import glob
import os
import shutil
import metpy.calc
import pandas as pd
from metpy.units import units

# Find and copy everything out of sim directory
wd="/home/caduff/sim/paiku"
outpath =wd+"/cryogrid_forcing"

if not os.path.exists(outpath):
	os.makedirs(outpath)

metfiles =glob.glob(wd+"/sim/g*/forcing/meteoc*.csv")
metfiles_sorted = sorted(metfiles)

mydirs =  glob.glob(wd+"/sim/g*/")

for mydir in mydirs:
	print (mydir)

	outdir = mydir.split('/')[-2]
	myout = outpath+"/"+outdir
	inpath= wd+"/sim/"+outdir

	if not os.path.exists(myout):
		os.makedirs(myout)

	shutil.copy2(inpath+"/landform.tif", myout)
	shutil.copy2(inpath+"/listpoints.txt", myout)

	metfiles =glob.glob(mydir+"/forcing/meteoc*.csv")
	for metfile in metfiles:
		print(metfile)
		shutil.copy2(metfile, myout)



# process and write new files
metfiles =glob.glob(outpath+"/g*//meteoc*.csv")

for file in metfiles:
	file="/home/joel/sim/paiku/cryogrid_forcing/g1/meteoc1.csv"
	outname = file.split(".csv")[0]+"_cryo.csv"
	df= pd.read_csv(file, index_col=0, parse_dates=True)


	# specific humidity calc
	temp = np.array(df.TA) *units('K')
	rh=np.array(df.RH) *units('percent')
	P=np.array(df.P) *units('Pa')

	dp = metpy.calc.dewpoint_from_relative_humidity(temp, rh*100)
	df["SH"] = metpy.calc.specific_humidity_from_dewpoint(dp, P)

	# wind spedd threshold for stability
	df.VW[df.VW<0.5] = 0.5
	df.to_csv(path_or_buf=outname ,na_rep=-999,float_format='%.4f', header=True, sep=',')

