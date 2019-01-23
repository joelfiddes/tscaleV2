.. tscale documentation master file, created by
   sphinx-quickstart on Tue Jan 22 22:33:03 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TopoSCALE
=========

TopoSCALE is a downscaling tool that uses the well-resolved description of the atmospheric column provided by atmospheric models, together with high-resolution digital elevation models (DEMs), to downscale coarse-grid climate variables to a fine-scale subgrid. The main aim of this approach is to provide high-resolution driving data for a land-surface model (LSM). Full desciption available in this publication: http://www.geosci-model-dev.net/7/387/2014/

Dependencies
------------
Python specific dependencies are bundled in the virtual env distributed with the code repository and defined in the `requirements.txt <https://github.com/joelfiddes/tscaleV2/blob/master/requirements.txt>`_ file . Additional system requirments:

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

Modes of operation
------------------
POINT: 
	Generates downscaled timeseries for a specific point defined by long/lat
TSUB: 
	Generates downsclaed timeseried for TopoSUB cluster centroids
GRID:
	Generates a 2D grid of surface fields corresponding to input DEM.

Tutorial
--------
The github repo contains a full example to set up a first prototype and 
understand how the code runs.

Documentation
-------------

.. toctree::
   :maxdepth: 2

   components.rst
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




