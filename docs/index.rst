.. tscale documentation master file, created by
   sphinx-quickstart on Tue Jan 22 22:33:03 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TopoSCALE
=========

TopoSCALE is a downscaling tool that uses the well-resolved description of the atmospheric column provided by atmospheric models, together with high-resolution digital elevation models (DEMs), to downscale coarse-grid climate variables to a fine-scale subgrid. The main aim of this approach is to provide high-resolution driving data for a land-surface model (LSM). Full desciption available in this publication: http://www.geosci-model-dev.net/7/387/2014/

Dependencies
------------
Python specific dependencies are bundled in the virtual env distributed with the code repository and defined in the 'requirements.txt file <https://github.com/joelfiddes/tscaleV2/blob/master/requirements.txt>'. Additional system requirments:

- linux (not tested on other platforms)
- python 2.7 (to be updated to py3)
- pip::

	apt-get install pip

- Climate Data Operators (CDO) used to concat monthly netcdfs downloaded from ECMWF::

	apt-get install cdo

- ECMWF CDS API set up::

	https://confluence.ecmwf.int/display/CKB/How+to+migrate+from+ECMWF+Web+API+to+CDS+API


Quickstart
----------
Get the code using git::

	git clone https://github.com/joelfiddes/tscaleV2.git

or direct download:: 

	https://github.com/joelfiddes/tscaleV2/archive/master.zip

Activate virtual environment::

	cd ./tscale
	source env/bin/activate

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


Documentation
-------------

.. toctree::
   :maxdepth: 2


   main.rst
   era5.rst 
   tscale.rst
   solarGeom.rst
   writeDocs.rst


Search code
-----------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`




