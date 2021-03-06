'''

Set of methods to interpolate ERA-Interim data based on station attributes or high res 
dem.

Methods adapted from REDCAPP available at: https://github.com/geocryology/REDCAPP

'''
import numpy as np

import netCDF4 as nc
#import pygrib  as pg
import csv
import re

from scipy.interpolate import RegularGridInterpolator
#from scipy.ndimage import gaussian_filter,generic_filter,convolve,minimum_filter,maximum_filter
from math import radians, exp, floor

from os import path, remove
import glob as gl


class rawData(object):
    """
    Args:
        dir_data: directory containing all raw data and output data
        
    Example:
        dir_data = 'C:/users/bincao/Desktop/data'
        dataImport = rawData(dir_data)
        sa = dataImport.saf_get()#get sa file in the given directory
        pl = dataImport.plf_get()#get pl file in the given directory
    """
    
    def __init__(self, dir_data):
        self.dir = dir_data
        

            
    def asciiDemHeader(self, demAsiccf):
        """Returns header information of input DEM in ASCIIGRID format"""
        header = []
        with open(demAsiccf) as f:
            reader = csv.reader(f)
            i=0
            for row in reader:
                i = i+1
                if (i <= 5):
                    header.append(row[0])
            return header
            
    def asciiDemEle(self, demAsciif):
        """
        
        Args:
           demAsciif: DEM in ASCIIGRID format and lat/lon WGS84 grid
           
         Returns:
             headerInfo: header information of input ascii file
             ele: array-like elevation of input ascii file
        """

        #read elevation
        return np.loadtxt(demAsciif, delimiter=' ', skiprows = 5)
        
    def ascii2ncdf(self, demAsciif, dem_out):
        """convert input DEM in ASCIIGRID to netcdf file
        
        Args:
            demAsciif: DEM in ASCIIGRID format and lat/lon WGS84 grid. 
                the header of ascii file should be same as the example data,
                xllcorner is the lontitude in the left low point, while
                yllcorner is the latitude in the left low point.
            dem_out: output DEM in netcdf format
            
        Returns a ASCIIGRID-based DEM in netcdf
        
        Example:
            dir_data = 'C:/users/bincao/Desktop/data'
            demAsciif = 'DEM_testArea.asc'
            dem_out  = 'C:/users/bincao/Desktop/data/DEM_fine-scale.nc'
            
            dataImport =  dataImport = rawData(dir_data)
            dataImport.ascii2ncdf(dem_file, dem_out)
            
        """
        
        #meta information
        header = self.asciiDemHeader(demAsciif)
        ncol   = int(re.findall(r'\d+', header[0])[0])
        nrow   = int(re.findall(r'\d+', header[1])[0])
        xllcorner = float(re.findall(r'\d+\.\d+', header[2])[0])
        yllcorner = float(re.findall(r'\d+\.\d+', header[3])[0])
        cellsize  = float(re.findall(r'\d+\.\d+', header[4])[0])
        
        #get variables
        ele = self.asciiDemEle(demAsciif)# elevation
        lats = np.linspace(yllcorner, yllcorner+cellsize*nrow, 
                           nrow, endpoint= True)# latitude
        lons = np.linspace(xllcorner, xllcorner+cellsize*ncol, 
                           ncol, endpoint= True)# lontitude
        
        #create nc file
        nc_root = nc.Dataset(dem_out ,'w', format = 'NETCDF4_CLASSIC')
        
        #create dimensions
        nc_root.createDimension('lat', nrow)
        nc_root.createDimension('lon', ncol)
        
        #create variables
        longitudes = nc_root.createVariable('lon', 'f4', ('lon'))
        latitudes  = nc_root.createVariable('lat', 'f4', ('lat'))
        elevation  = nc_root.createVariable('elevation', 
                                           'f4', ('lat', 'lon'), zlib = True)
        
        #assign variables
        longitudes[:] = lons
        latitudes[:]  = lats[::-1]
        elevation[:]  = ele
        
        #attribute
        nc_root.description = "high-resolution topography file"
        #resolution = cellsize
        longitudes.units = 'degree_east (decimal)'
        latitudes.units  = 'degree_north (decimal)'
        elevation.units  = 'm'
        
        nc_root.close()



class tscale3dPl(object):
    """
    Return object for downscaling that has methods for interpolationg
    upper-air temperature and surface influences at surface level
    based on disaggregating coarse-grid reanalysis and dem.
    
    Args:
        dem: A required fine-scale dem in netcdf format
    
    Example:
        dem  = 'example_alps.nc'
        geop = 'alps_geop.nc'
        sa   = 'alps_sa_79_15.nc'
        pl   = 'alps_pl_79_15.nc'
        
        downscaling = downscaling(dem, geop, sa, pl)
        
    """
    
    def __init__(self, pl, dem = None):
        self.g    = 9.80665 #Gravitational acceleration [m/s2]
     
        self.pl   = nc.Dataset(pl)
        if not (dem is None):
            self.dem  = nc.Dataset(dem)
        
        
    def demGrid(self, stations = None):     
        """Return metadata of given stations or dem. 
        Format of stations desired [lat, lon, geop].
        
        Args: 
            stations: A list of dictionaries describing stations. If not given,
                metadata derived from given dem.
            
        Returns:
            out_xyz_dem: Metadata [lat, lon, geop] of input dem or sites
            lons: Longitude of input sites
            lats: Latitude of input sites
            shape: Shape of input dem or sites
            names: Name of input stations. Names will only be avaiable when 
                stations are inputted
        """
        
        if not (stations is None):
            lats = [s['lat'] for s in stations]
            lons = [s['lon'] for s in stations]
            names = [s['name'] for s in stations]
            siteLocation = np.asarray([[s['lat'],s['lon'],
                                        s['ele']*self.g] for s in stations])
            shape = siteLocation.shape
            return siteLocation, lats, lons, shape, names        
        else:
            #out_xyz based on dem
            lons = self.dem.variables['lon'][:]
            lats = self.dem.variables['lat'][:]
            geop = self.dem.variables['elevation'][:]*self.g
            shape = geop.shape
        
            lons, lats = np.meshgrid(lons, lats)
            lons = lons.reshape(lons.size)
            lats = lats.reshape(lats.size)
            geop = geop.reshape(geop.size)
            out_xyz_dem = np.array([lats, lons, geop]).T
            return out_xyz_dem, lats, lons, shape
        
       

    def gridValue(self, variable, ind_time):
        """
        Return original grid temperatures and geopotential of differnet
        pressure levels. The function are called by inLevelInterp() to
        get the input ERA-Interim values.
        
        Args: 
            variable: Given interpolated climate variable
            ind_time: Time need to be interpolated. Time is in interger (e.g.
            0, 1, 2)
            
        Returns:
            gridT: Grid temperatures of different pressure levels. Retruned 
            temperature are formated in [level, lat, lon]
            gridZ: Grid geopotential of different pressure levels. Retruned 
            temperature are formated in [level, lat, lon]
            gridLon: Grid longitude of pressure level variables
            gridLat: Grid latitude of pressure level variables
        
        Example:
            gridT,gridZ,gridLat,gridLon=downscaling.gridValue('Temperature',0)
            
        """
        
        gridT = self.pl.variables[variable][ind_time,:,:,:]
        gridZ = self.pl.variables['Geopotential'][ind_time,:,:,:]
        #x and y
        
        gridLat = self.pl['lat'][:]
        gridLon = self.pl['lon'][:]
        

        return gridT,gridZ,gridLat,gridLon


    def inLevelInterp(self,gridT, gridZ, gridLat, gridLon, out_xyz):
        """
        This is a 2D interpolatation, and returns interpolated temperatures
        of different pressure levels.
        
        Args:
            gridT: Grid temperatures of different pressure levels. Retruned 
                temperature are formated in [level, lat, lon]
            gridZ: Grid geopotential of different pressure levels. Retruned 
                temperature are formated in [level, lat, lon]
            gridLat: Grid longitude of pressure level variables
            gridLon: Grid latitude of pressure level variables
            out_xyz: Given sites, which will be interpolated.
            
        Returns:
            t_interp: Interpolated temperatre of different pressure levels. 
                The returned values are fomrated in [level, lat, lon]
            z_interp: Interpolated geopotential of different pressure levels. 
                The returned values are fomrated in [level, lat, lon]
        
        Examples:
            downscaling = downscaling(dem, geop, sa, pl)

            out_xyz_dem, lats, lons, shape = downscaling.demGrid()
            out_xyz_sur = downscaling.surGrid(lats, lons, None)

            #interpolate 2-meter temperature
            surTa = downscaling.surTa(0, out_xyz_sur)
            #original ERA-I values
            gridT,gridZ,gridLat,gridLon = downscaling.gridValue(variable,0)
            #interpolate temperatures and geopotential of different 
            pressure levels.

            t_interp, z_interp = downscaling.inLevelInterp(gridT,gridZ,
                                                           gridLat,gridLon,
                                                           out_xyz_dem)
        """
        
        shape = gridT.shape
        #create array to hold interpolation resultes
        t_interp = np.zeros([shape[0], len(out_xyz)])
        z_interp = np.zeros([shape[0], len(out_xyz)])

        #temperatue and elevation interpolation 2d
        for i in range(shape[0]):
            ft = RegularGridInterpolator((gridLat,gridLon), 
                                          gridT[i,:,:], 'linear')
            fz = RegularGridInterpolator((gridLat,gridLon), 
                                          gridZ[i,:,:], 'linear')
            t_interp[i,:] = ft(out_xyz[:,:2])#temperature
            z_interp[i,:] = fz(out_xyz[:,:2])#elevation

        t_interp -= 273.15

        return t_interp[::-1,:], z_interp[::-1,:]
        
    

        
    def fast1d(self, t_interp, z_interp, out_xyz):
        """This is a 1D interpoation. The function return interpolated 
        upper air temperature at the given sites by 
        interpolation between different pressure levels.
        
        Args:
            t_interp: Interpolated temperatre of different pressure levels. 
                The returned values are fomrated in [level, lat, lon]
            z_interp: Interpolated geopotential of different pressure levels. 
                The returned values are fomrated in [level, lat, lon]
            out_xyz: Given sites with elevation, which will be interpolated.
            
        Returns:
            dG:upper-air temperature at given sites
                
        Example: 
            downscaling = downscaling(dem, geop, sa, pl)

            out_xyz_dem, lats, lons, shape = downscaling.demGrid()
            out_xyz_sur = downscaling.surGrid(lats, lons, None)

            
            surTa = downscaling.surTa(0, out_xyz_sur)
            #original ERA-I values
            gridT,gridZ,gridLat,gridLon = downscaling.gridValue(variable,0)
            #interpolate temperatures and geopotential of different 
            pressure levels.

            t_interp, z_interp = downscaling.inLevelInterp(gridT,gridZ,
                                                           gridLat,gridLon,
                                                           out_xyz_dem)
            
            #upper air temperature at the coarse and fine scale of elevation
            pl_sa = fast1d(t_interp, z_interp, out_xyz_sur)
            pl_obs = fast1d(t_interp, z_interp, out_xyz_dem)
        """
    
        ele = out_xyz[:,2]
        size = np.arange(out_xyz.shape[0])
        n = [bisect_left(z_interp[:,i], ele[i]) for i in size]
        n = [x+1 if x == 0 else x for x in n]
        
        lowN = [l-1 for l in n]
        
        upperT = t_interp[n,size]
        upperZ = z_interp[n,size]
        dG  = upperT-t_interp[lowN,size]#<0
        dG /= upperZ-z_interp[lowN,size]#<0
        dG *= out_xyz[:,2] - upperZ#>0
        dG += upperT
             
        return dG
    
          
        


