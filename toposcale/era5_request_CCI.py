
"""
This module downloads global ERA5 pressure level and surface data using the CDS API. Options are set as variables directly in this script.

Example:
        $ python era5_request_CCI.py

Vars:
     	startYear: start year is download  (int)
        endYear: end year is not downloaded (int)
        era5dir:  directory to write    (string)
        timestep: timestep of data in h (string)
        plevels: list of pressure levels to retrieve (list)


"""
startYear = 2019
endYear = 2020
era5dir= "/cluster/projects/nn9606k/era5/"
timestep = "6"
plevels = ["1000", "700", "500", "200"]

import fetch_era5_global as e5
e5.retrieve_era5_surf(startYear,endYear, era5dir , timestep)
e5.retrieve_era5_plev(startYear,endYear, era5dir, timestep, plevels)
