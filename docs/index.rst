.. tscale documentation master file, created by
   sphinx-quickstart on Tue Jan 22 22:33:03 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TopoSCALE
=========

TopoSCALE is a downscaling tool that uses the well-resolved description of the atmospheric column provided by atmospheric models, together with high-resolution digital elevation models (DEMs), to downscale coarse-grid climate variables to a fine-scale subgrid. The main aim of this approach is to provide high-resolution driving data for a land-surface model (LSM). Full desciption available in this publication: http://www.geosci-model-dev.net/7/387/2014/

Dependencies
------------
Python spefic dependencies are bundled in the virtual env distributed with the code repository. Additional system requirments:

- linux (not tested on other platforms)
- python 2.7 (to be updated to py3)
- pip::

	apt-get install pip

- cdo::

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
Documentation
-------------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

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




