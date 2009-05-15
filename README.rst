:Date: 15 May 2009
:Author: Anand Patil
:Contact: anand.prabhakar.patil@gmail.com
:Web site: github.com/malaria-atlas-project/generic-mbg
:Copyright: Anand Patil, 2009.
:License: Creative Commons BY-NC-SA, see LICENSE


The generic MBG package allows us to write a PyMC probability model for a project's 
spatial or spatiotemporal count data (to a certain specification), then easily use 
it to fit a dataset & predict using the following three shell commands:

* ``mbg-infer`` runs the MCMC algorithm using the given model & an input dataset,
  stored in a csv file, and stores the traces in an HDF5 archive.

* ``mbg-map`` takes the HDF5 archive produced by mbg-infer, and an ASCII file with
  a MISSING entry in its header. Produces a set of bespoke summary maps on the grid
  expressed by the ASCII header. The missing pixels are missing in the output also.
  
* ``mbg-validate`` takes the HDF5 archive produced by mbg-infer and a 'holdout'
  dataset, stored in a csv file, and returns a set of predictive samples at the
  holdout locations.
  
If the project's members are interested in changing the model or specifying a
subjective prior, there are two additional shell commands available:

* ``mbg-scalar-priors`` draws samples from the prior for all scalar parameters
  (including deterministics) and plots histograms for inspection.
  
* ``mbg-realize-prior`` draws all scalar parameters from the prior, and realizes
  and plots the random field on grids matching a number of input ASCIIs.

All shell commands can be run with only the ``-h`` option to print some help to the
screen. However, if you're reading this document, you don't really need to do that.

***************************
Detailed usage instructions
***************************

If you want to use the shell commands, this section is for you.

``mbg-infer``
=============
::

    mbg-infer module database-file input [options]
    
Produces the requested database file. Also produces plots of the dynamic traces of all
scalar parameters as PDF's, and saves them in the folder ``name-plots``, where ``name``
is the name of the database file. You will need to inspect these plots to determine how
many 'burnin' iterations should be discarded when making maps.

In the future, if the database file already exists, ``mbg-infer`` will give you the option
to append samples to it. For now, it will only give you the option to overwrite it.

Required arguments
------------------

1. The name of the module containing the model specification.

2. The name of the database file to be produced. If you do not want it to go in the current
   directory, specify a path, eg ``/home/anand/traces/run-01-04-2009``.

3. The name of a csv file containing the input data. If it is a different directory, specify
   the path to it, eg ``/home/anand/data/query-01-04-2009.csv``. This csv file must have the
   following columns:
     
     * ``lon``, ``lat`` : The coordinates of the observation in decimal degrees
     
     * ``pos``, ``neg`` : The number of 'positive' and 'negative' observations.
     
     * ``t`` : Time in decimal years since 2009. This is only required for 
       spatiotemporal models.

   All other columns are interpreted as covariates, eg ``ndvi`` etc.
   

Options
-------

* ``-t`` or ``--thin`` : If thin is 10, every 10th MCMC iteration will be stored in the 
  database. Small values are good but slow. 1 is best.

* ``-i`` or ``--iter`` : The number of MCMC iterations to perform. Large values are good
  but slow.

* ``-n`` or ``-ncpus`` : The maximum number of CPU cores to make available to the MCMC 
  algorithm. Should be less than or equal to the number of cores in your computer. The 
  All the cores you make available may not be utilized. Use top or the Activity Monitor
  to monitor your actual CPU usage. Large values are good but tie up more of your computer.



``mbg-map``
===========
::

    mbg-map module database-file burn mask [options]

Produces a folder called ``name-maps`` where ``name`` is the name of the database file.
Puts the requested maps in the folder in ascii format. Also produces PDF images of all
the requested maps for quick viewing.

Required arguments
------------------

1. The name of the module containing the model specification.

2. The name of the database file (produced by mbg-infer) to be used to generate the 
   maps. If you do not want it to go in the current directory, specify a path.
   
3. The number of burnin iterations to discard from the trace before making the maps.
   You will need to figure this out by inspecting the traces produced by ``mbg-infer``.
   
4. The name of an ASCII file. The maps will be produced in ASCII files with identical
   headers and identical MISSING pixels. If the file is in a different directory, specify
   the path to it.

Options
-------

* ``-n`` or ``--n-bins`` : The number of bins to use in the histogram from which quantiles
  are computed. Large values are good, but use up more system memory. Decrease this if you
  see memory errors.

* ``-b`` or ``--bufsize`` : The number of buffer pixels to render around the edges of the
  continents. Set to zero unless the ``raster-thin`` option is greater than 1. The buffer
  will not be very good. In general, if you want a buffer you're better off making your 
  own in ArcView rather than using this option.

* ``-q`` or ``--quantiles`` : A string containing the quantiles you want. For example,
  ``'0.25 0.5 0.75'`` would map the lower and upper quartiles and the medial. Default is 
  ``'0.05 0.25 0.5 0.75 0.95'``.

* ``-r`` or ``--raster-thin`` : If you think your map is going to be smooth (eg because you
  aren't using any covariates), you can use this option to render the maps on a degraded grid,
  then interpolate back to the original grid using splines. For instance, if your input ASCII
  is on a 5km grid, and you use ``-r 5``, the maps will be rendered on a 25km grid, then
  interpolated back to a 5km grid when it is time to produce the output ASCIIs. Small values
  are good but slow. 1 is best.

* ``-t`` or ``--thin`` : The factor by which to thin the MCMC trace stored in the database.
  If you use ``-t 10``, only every 10th stored MCMC iteration will be used to produce the maps.
  Small values are good but slow. 1 is best.

* ``-i`` or ``--iter`` : The total number of predictive samples to use in generating the maps.
  Large values are good but slow. Defaults to 20000.

* ``-c`` or ``--covariates`` : A list of names of ASCII files containing the covariate rasters.
  These files' headers must match those of the input raster, and their missing pixels must match
  those of the input raster also. There must be a file corresponding to every covariate column
  in input 3 of mbg-infer. For example, if you used ``rain`` and ``ndvi`` as your column headers,
  you should use ``-c 'rain.asc ndvi.asc'``. If the rasters are in another folder, specify the path,
  ie ``-c '/home/noor/rain.asc /home/noor/ndvi.asc'``

* ``-y`` or ``--year`` : If your model is spatiotemporal, you must provide the decimal year since
  2009 at which you want your map produced. For example, Jan 1 2008 would be ``-y -1.0``.


``mbg-validate``
================
::

    mbg-validate module database-file burn pred-pts [options]
    
The output format is likely to change in the future. For the time being, the output is 
simply a csv file containing a number of posterior predictive samples at the locations
specified in the ``pred-pts`` file.

Required arguments
------------------

1. The name of the module containing the model specification.

2. The name of the database file (produced by mbg-infer) to be used to generate the 
   maps. If you do not want it to go in the current directory, specify a path.
   
3. The number of burnin iterations to discard from the trace before making the maps.
   You will need to figure this out by inspecting the traces produced by ``mbg-infer``.
   
4. A csv file containing the 'holdout' dataset. It should be in exactly the same format
   as the third required input to ``mbg-infer``.

Options
-------

* ``-t`` or ``--thin`` : The factor by which to thin the MCMC trace stored in the database.
  Small values are good but slow. 1 is best.

* ``-i`` or ``--iter`` : The total number of predictive samples you want to generate. Large
  values are good but slow. Defaults to 20000.


``mbg-scalar-priors``
=====================
::

    mbg-scalar-priors module [options]

Required arguments
------------------

1. The name of the module containing the model specification.

Options
-------

* ``-i`` or ``--iter`` : The total number of predictive samples you want to generate. Large
  values are good but slow. Defaults to 20000.

``mbg-realize-prior``
=====================

This one is not implemented yet.




*******************
Module requirements
*******************

This section tells you how to write new modules that will work with the shell commands.


*******************************************
Questions asked with high prior probability
*******************************************

* **Q**: Can you make this work on Windows?

  **A**: No.