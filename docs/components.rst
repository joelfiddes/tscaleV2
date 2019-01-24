TopoSCALE components
--------------------
TopoSCALE has 5 components:

INI file
^^^^^^^^
The INI file defines all required options for a run.

Input plugins
^^^^^^^^^^^^^
	* preprocesses product (resampling, conversions)
	* converts to generic python class structure

:py:mod:`era5`

Core engine
^^^^^^^^^^^
	* accepts generic data structures
	* algorithms
	* outputs generic  python class structure


Output plugins
^^^^^^^^^^^^^^
	* writes specific output formats (CSV, NetCDF)
	* writes model specific outputs (SMET, GEotop, Cryogrid) 
	* writes grids (NetCDF, tiff)

Optional modules
^^^^^^^^^^^^^^^^
	* retrieve and preprocess data products (various api dependencies here)

