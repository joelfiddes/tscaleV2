# ## TopoPAR.m
# # Calculates topographic parameters, including solar geometry, for a
# # given region and time period as specified by an input time array and
# # an input digital elevation model (DEM) structure.
# # Inputs: 
# #		   t = [Nt x 1] time array. Time is in UTC (Zulu time) as the
# #		   number of days since 01-Jan-0000 + 1. This is the standard
# #		   MATLAB datenum format.
# #		   DEM = A DEM structure. This must contain a 'utmz' field
# #		   specifiying the UTM zone, a 'domain' field that contains the UTM
# #		   coordinate arrays  (x,y) and elevation (z) for the raster
# #		   domain that the topographic parameters will be calculated over,
# #		   and a 'pad' field that comtains the same arrays as the 'domain'
# #		   field but for a larger padded domain that encompasses the true
# #		   domain.
# # Outputs:
# #		   tp = A topographic parameter strucutre that contains the
# #		   aspect, slope, elevation, solar zenith angle (time series),
# #		   solar azimuth angle (time series), horizon angle and sky view
# #		   factor for the raster domain.
# # Requires:
# #		   pvl_ephemeris.m : A function to calculate solar geometry based
# #		   on time and position. Available through
# #		   https://pvpmc.sandia.gov/applications/pv_lib-toolbox/
# # References:
# # Zaksek et al. (2011), doi:10.3390/rs3020398 
# # Dozier and Frew (1990), doi:10.1109/36.58986
# #
# # Adaptated from matlab Code by Kristoffer Aalstad (March 7, 2018).

import gdal as gd

# read tif
ds = gd.Open('ele_utm_sml.tif')

# convert to array
a = ds.ReadAsArray()

# Extract DEM data from the DEM array.
width = ds.RasterXSize
height = ds.RasterYSize
gt = ds.GetGeoTransform()
xmin = gt[0]
ymin = gt[3] + width*gt[4] + height*gt[5] 
xmax = gt[0] + width*gt[1] + height*gt[2]
ymax = gt[3]
xres= gt[1]
yres=-gt[5]

x0 = np.arange(xmin, xmax, xres)
y0 = np.arange(ymin,ymax, yres)
z0=a

x = np.arange(xmin, xmax, xres)
y = np.arange(ymin,ymax, yres)
z=a
zall=z

#z=z[2:-2,2:-2]															#z,x,y, subsets cause problems by removing required pixels for dz calc in main loop
#x=x[2:-2]
#y=y[2:-2]
X,Y=np.meshgrid(x,y)


del1=x[2]-x[1]

thesex=np.array([(x>=xmin) & (x<=xmax)])
thesey=np.array([(y>=ymin) & (y<=ymax)])


# Elevation derivatives. SUB zres/yres here instead of del1 as these could be different

dzdy=-1*(zall[3:-1,2:-2]-zall[1:-3,2:-2])/(2*yres) # N.B. y decreases with i.
dzdx=(zall[2:-2,3:-1]-zall[2:-2,1:-3])/(2*xres)

# Slope and aspect. 
# Aspect definition: N,E,S,W = pi, pi/2, 0, -pi/2
asp1=np.arctan2(-1*dzdx,dzdy)
slp1=np.arctan(np.sqrt(dzdx**2+dzdy**2))
#asp=asp1[thesey[0:-1],thesex] # why is thesey 1 too long?							#!!!VECTOR LENGTH WRONG!!!
#slp=slp1[thesey[0:-1],thesex]


# Horizon angles.
# Quite slow due to detailed computation, but you only have to do this
# once per site.
nbins=360
#bins=0:1:359
bins = np.arange(0,360,1)
ncbins=36
#cbins=0:10:350
cbins = np.arange(0,351,10)

ii=len(y0)
jj=len(x0)
h=np.zeros((ii,jj,ncbins), dtype=np.uint8)
ea=np.zeros((nbins,1))


# this could be interesting https://ch.mathworks.com/matlabcentral/fileexchange/59421-dem-based-topography-horizon-model
for j in range(0:jj):
	d=str(round(1e2*j/jj))
	print ('Terrain calculation is '+ d + ' percent complete')

	for i in range(0:ii):
		xis=x0[j]
		yis=y0[i]
		d=np.sqrt((X-xis)**2+(Y-yis)**2)[:,0:-1]									# !!!length fudge here!!!
		dX=X-xis
		dY=Y-yis
		bearing=np.arctan2(dX,dY)
		bearing=np.degrees(bearing+2*np.pi*(bearing<0))
		
		dz=z-z[d==0] # Change in elevation.										# !!!VECTOR LENGTH WRONG d!=z!!!
		eang=np.degrees(np.arctan(dz/d) ) # Elevation angle.					NAN produced at arctan(0)
		
		for n in range(0:nbins):
			if n!=nbins:
				blim=[bins[n], bins[n+1]]
			else:
				blim=[bins[n] ,360]
			
			# Define a 1 km long central vector for this bin,
			# to also check the closest grid cells whos centroids
			# may fall just outside the bearing bin, even if
			# the cell is inside the bin. 
			dxb=np.mean(np.sin(np.radians(blim))*1e2)
			dyb=np.mean(np.cos(np.radians(blim))*1e2)
			#xb=xis:dxb:(xis+10*dxb)
			xb = np.arange(xis,(xis+10*dxb),dxb)
			#yb=yis:dyb:(yis+10*dyb)
			yb = np.arange(yis,(yis+10*dyb),dyb)


			for k in range(0:len(xb)):
				db=np.sqrt((X-xb[k])**2+(Y-yb[k])**2)
				inbin=[db<=np.sqrt(2*5e1**2)] 								# booleen away same dims as grid right? Not a subset....
			
			# Find the closest cell with the maximum elevation angle within this bin.
			#inbin=inbin [ bearing>=blim[0] & bearing<blim[1]]
			inbin2=np.array([np.logical_and(bearing>=blim[0] ,  bearing<blim[1])]).squeeze()
			inbin3 = inbin2[:,0:-1]													# size fudge
			thesea=eang[inbin3]
			maxa=np.array([np.logical_and(inbin3==np.nanmax(thesea), eang==np.nanmax(thesea))]).squeeze()
			maxa = np.array(maxa)
			thesed=d[maxa]
			#here=maxa&d==np.nanmin(thesed)
			here = [np.logical_and(maxa==np.nanmin(thesed), d==np.nanmin(thesed))]
			ea[n]=eang[here]*(eang[here]>0) # 0^deg as lower bound on elevation angle.
		
		
		# Take average of high res bins in coarser 10^deg bins.
		 for n in range(1:ncbins):
			 if n!=ncbins:
				blim=[cbins(n) cbins(n+1)]  
			else:
				blim=[cbins(n) 360]
			
			these=bins>=blim(1)&bins<blim(2)
			avh=mean(ea(these))
			h(i,j,n)=np.uint8(avh)
		 end
	end
end

# Save the center of the coarse bins. Convert these to the same coordinate
# system as the aspect. 
hbins=cbins+5
hbins=deg2rad(hbins)
hbins=(5*pi/2)-hbins
hbins=hbins-2.*pi.*(hbins>2.*pi)
hbins=hbins+pi/2
hbins=hbins-2.*pi.*(hbins>pi)

