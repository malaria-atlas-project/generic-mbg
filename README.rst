:Date: 15 May 2009
:Author: Anand Patil
:Contact: anand.prabhakar.patil@gmail.com
:Web site: github.com/malaria-atlas-project/generic-mbg
:Copyright: Anand Patil, 2009.
:License: GPL, see GPL in this directory.

The generic MBG package allows us to write PyMC probability models for each 
project that works with some kind of spatial GLM, then turn the model over 
to the project team for testing, fitting, mapping and experimentation using 
a few easy shell commands:

* ``mbg-infer`` runs the MCMC algorithm using the given model & an input dataset,
  stored in a csv file, and stores the traces in an HDF5 archive.

* ``mbg-map`` takes the HDF5 archive produced by mbg-infer, and a raster with
  some pixels missing. Produces a set of bespoke summary maps of each predicted
  quantity that match the raster in terms of grid, missingness pattern and file 
  format. 
  
* ``mbg-3dmap`` takes the HDF5 archive produced by mbg-infer, and a raster with
  some pixels missing. Outputs the full probability density function of each
  predicted quantity over a thinned version of the raster, and stores it under 
  compression in HDF5 format. This file can be opened and examined graphically 
  using MayaVI, or converted into two-dimensional maps later.
  
* ``mbg-areal-predict`` takes the HDF5 archive produced by mbg-infer, a raster
  with some pixels missing and a text file containing one or more multipolygons
  in geojson format, which must be tagget with time intervals for space-time models. 
  Produces samples from the predictive distribution of certain integrals over the 
  given multipolygons. See the PDF documentation for more detail.
  
* ``mbg-validate`` takes the HDF5 archive produced by mbg-infer and a 'holdout'
  dataset, stored in a csv file, and creates a set of predictive samples at the
  holdout locations and some validation plots.
  
* ``mbg-decluster`` partitions a CSV datafile into 'kept' and 'holdout' portions.

* ``mbg-describe-tracefile`` examines an HDF5 archive produced by mbg-infer, and
  tells you when it was produced, which versions of the code were used, how many
  iterations it contains and what the input data were.
  
If the project's members are interested in changing the model or specifying a
subjective prior, there are two additional shell commands available to help:

* ``mbg-scalar-priors`` draws samples from the prior for all scalar parameters
  (including deterministics) and plots histograms for inspection.
  
* ``mbg-realize-prior`` draws all scalar parameters from the prior, and realizes
  and plots the random field on grids matching a number of input ASCIIs.
  
* ``mbg-describe-tracefile`` provides information about the circumstances under which
  traces were produced.

All shell commands can be run with only the ``-h`` option to print some help to the
screen. However, if you're reading this document, you don't really need to do that.

For documentation on how to install and administer the package \& write specializing
modules, run ``builddocs.py`` in the ``docs`` directory and then see ``manual.pdf``.


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

If you determine that more MCMC samples are needed, simply run mbg-infer with the same 
database file argument to pick up where you left off and keep sampling.

Required arguments
------------------

1. The name of the module containing the model specification.

2. The name of the database file to be produced. If you do not want it to go in the current
   directory, specify a path, eg ``/home/anand/traces/run-01-04-2009``. If the database file
   already exists, you will be prompted about whether you want to continue sampling into it
   or remove it.

3. The name of a csv file containing the input data. If it is a different directory, specify
   the path to it, eg ``/home/anand/data/query-01-04-2009.csv``. This csv file must have the
   following columns:
     
   * ``lon``, ``lat`` : The coordinates of the observation in decimal degrees
     
   * ``t`` : Time in decimal years. This is only required for spatiotemporal models.

   All other columns are interpreted as covariates, eg ``ndvi`` etc., UNLESS the module 
   implements the ``non_cov_columns`` attribute. See the PDF documentation.
   

Options
-------

* ``-t`` or ``--thin`` : If thin is 10, every 10th MCMC iteration will be stored in the 
  database. Small values are good but slow. 1 is best.

* ``-i`` or ``--iter`` : The number of MCMC iterations to perform. Large values are good
  but slow. Defaults to 100000.

* ``-n`` or ``-ncpus`` : The maximum number of CPU cores to make available to the MCMC 
  algorithm. Should be less than or equal to the number of cores in your computer. The 
  All the cores you make available may not be utilized. Use top or the Activity Monitor
  to monitor your actual CPU usage. Large values are good but tie up more of your computer.
  Defaults to the value of the environment variable ``OMP_NUM_THREADS``.

``mbg-describe-tracefile``
==========================
::

    mbg-describe-tracefile path

If path is a database file, inspects the database file. Prints out the version of the 
generic package, the module that produced the file and the date the run was started. 
Writes the input data to csv with filename ``database-file-input-csv``, substituting 
the actual filename.

If the path is a directory, walks the filesystem starting from the directory, inspecting
every database file it finds. Does not produce any csv's.

Required arguments
------------------

1. The name of the database file or path to be inspected.


``mbg-covariate-traces``
========================
::

    mbg-covariate-traces module database-file [options]

Postprocesses the given database file to produce MCMC traces for the covariate 
coefficients. Produces a directory called database-file-covariate-traces, and populates 
it with pdf images of the covariate coefficient traces and  


Required arguments
------------------

1. The name of the module containing the model specification.

2. The name of the database file containing the MCMC trace.


Options
-------

* ``-t`` or ``--thin`` : If thin is 10, samples of the covariate coefficients will be
  produced for every 10th MCMC sample. Defaults to 1, meaning no thinning.

* ``-b`` or ``--burn`` : Samples of the covariate coefficients will begin after this
  many 'burnin' iterations are discarded. Defaults to 0, meaning no burnin.



``mbg-decluster``
========================
::

    mbg-decluster input prop [options]

A wrapper for the R function getdeclusteredsample that results in two new tables with 
suffix HOLDOUT and THINNED outut to same directory as tablepath  


Required arguments
------------------

1. (string) path to input table. must include columns 'lon' and 'lat'. If
   also 't' will treat as space-time. If only filename given (no path) assumes file
   in current working directory.

2. (float) what proportion of the full data set will be used for hold-out set.


Options
-------

* ``-m`` or ``--minsample`` : (int) optional minimum sample size (supercedes prop.
  if larger)

* ``-d`` or ``--decluster`` : (logical) do we want to draw spatially declustered
  sample (default) or just simple random.

* ``-p`` or ``--makeplot`` : (logical) do we want to export a pdf map showing
  location of data and selected points. This is exported to same directory as
  tablepathoptional minimum sample size (supercedes prop if larger).


``mbg-map``
===========
::

    mbg-map module database-file burn mask [options]

Produces a folder called ``name-maps`` where ``name`` is the name of the database file.
Puts the requested maps in the folder in format matching the mask. Also produces PDF 
images of all the requested maps for quick viewing.

Required arguments
------------------

1. The name of the module containing the model specification.

2. The name of the database file (produced by mbg-infer) to be used to generate the 
   maps. If you do not want it to go in the current directory, specify a path.
   
3. The number of burnin iterations to discard from the trace before making the maps.
   You will need to figure this out by inspecting the traces produced by ``mbg-infer``.
   
4. The name of a raster, without extension. The maps will be produced in raster files
   in the same format, on identical grids, with identical missing pixels. If the file 
   is in a different directory, specify the path to it.

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

* ``-t`` or ``--thin`` : The factor by which to thin the MCMC trace stored in the database.
  If you use ``-t 10``, only every 10th stored MCMC iteration will be used to produce the maps.
  Small values are good but slow. 1 is best. Defaults to 50.

* ``-i`` or ``--iter`` : The total number of predictive samples to use in generating the maps.
  Large values are good but slow. Defaults to 50000.

* ``-p`` or ``--raster-path`` : The path to the files containing the covariate rasters. These 
  files' headers must match those of the input raster, and their missing pixels must match
  those of the input raster also. There must be a file corresponding to every covariate column
  in input 3 of mbg-infer. For example, if you used ``rain`` and ``ndvi`` as your column headers,
  files ``rain.asc`` and ``ndvi.flt`` and ``temp.hdf5`` should be present in the raster path. 
  Defaults to the current working directory.

* ``-y`` or ``--year`` : If your model is spatiotemporal, you must provide the decimal year at 
  which you want your map produced. For example, Jan 1 2008 would be ``-y 2008``.
  
* ``-d`` or ``--ignore-npd`` : If ``1``, MCMC iterations whose covariance functions are non-
  positive-definite on the data locations plus the prediction locations will be ignored. If
  ``0``, any such iterations will result in errors. Defaults to 0.
  
* ``-u`` or ``--quantile-uplim`` : The upper limit of the mapped quantiles. Defaults to 1.

* ``-l`` or ``--quantile-lolim`` : The lower limit of the mapped quantiles. Defaults to 0.


``mbg-3dmap``
=============
::

    mbg-3dmap module database-file burn mask [options]

Produces a folder called ``name-3dmaps`` where ``name`` is the name of the database file.
Puts a HDF5 file- containing the probability density field of the output of each function
in the specializing module's ``map_postproc`` list in the folder. This data can be examined 
interactively using MayaVI. File ``display_3dmap.py``, included with the package, provides 
a template for scene generation.

Required arguments
------------------

1. The name of the module containing the model specification.

2. The name of the database file (produced by mbg-infer) to be used to generate the 
   maps. If you do not want it to go in the current directory, specify a path.

3. The number of burnin iterations to discard from the trace before making the maps.
   You will need to figure this out by inspecting the traces produced by ``mbg-infer``.

4. The name of a raster, without extension. The maps will be produced in raster files
   in the same format, on identical grids, with identical missing pixels. If the file 
   is in a different directory, specify the path to it.

Options
-------

* ``-n`` or ``--n-bins`` : The number of bins to use in the histogram from which quantiles
  are computed. Large values are good, but use up more system memory. Decrease this if you
  see memory errors. Defaults to 100.

* ``-b`` or ``--bufsize`` : The number of buffer pixels to render around the edges of the
  continents. Set to zero unless the ``raster-thin`` option is greater than 1. The buffer
  will not be very good. In general, if you want a buffer you're better off making your 
  own in ArcView rather than using this option. Defaults to 0.

* ``-q`` or ``--quantiles`` : A string containing the quantiles you want. For example,
  ``'0.25 0.5 0.75'`` would map the lower and upper quartiles and the medial. Default is 
  ``'0.05 0.25 0.5 0.75 0.95'``.

* ``-t`` or ``--thin`` : The factor by which to thin the MCMC trace stored in the database.
  If you use ``-t 10``, only every 10th stored MCMC iteration will be used to produce the maps.
  Small values are good but slow. 1 is best. Defaults to 50.
  
* ``-r`` or ``--raster-thin``: The 3d data cube takes up much more disk space and memory than
  the scalar maps. You might need to degrade the input raster to lower resolution. A value of
  10 means that the 3d maps will have 1/10 the spatial resolution of the input raster. Defaults
  to 1.

* ``-i`` or ``--iter`` : The total number of predictive samples to use in generating the maps.
  Large values are good but slow. Defaults to 50000.

* ``-p`` or ``--raster-path`` : The path to the files containing the covariate rasters. These 
  files' headers must match those of the input raster, and their missing pixels must match
  those of the input raster also. There must be a file corresponding to every covariate column
  in input 3 of mbg-infer. For example, if you used ``rain`` and ``ndvi`` as your column headers,
  files ``rain.asc`` and ``ndvi.flt`` and ``temp.hdf5`` should be present in the raster path. 
  Defaults to the current working directory.

* ``-y`` or ``--year`` : If your model is spatiotemporal, you must provide the decimal year at 
  which you want your map produced. For example, Jan 1 2008 would be ``-y 2008``.
  
* ``-d`` or ``--ignore-npd`` : If ``1``, MCMC iterations whose covariance functions are non-
  positive-definite on the data locations plus the prediction locations will be ignored. If
  ``0``, any such iterations will result in errors. Defaults to 0.


``mbg-areal-predict``
=====================
::

  mbg-areal-predict module database-file burn polyfile [options]

Produces a folder called ``name-areal-samples`` where ``name`` is the name of the 
database file. Populates this folder with files called ``fname-aname-samples`` and
``fname-aname-estimates``. ``fname-aname-samples`` is a CSV file with no column 
headers whose columns correspond to reps, and whose rows correspond to Monte Carlo
trials. ``fname-aname-estimates`` is a CSV file whose columns have titles 
``sname-estimate`` and ``sname-mcse``. ``aname`` iterates over the module's 
``areal-postproc`` functions and ``sname`` iterates over the summaries generated for
``mbg-map``, eg ``mean``, ``quantile-0.5``, etc.

Required arguments
------------------

1. The name of the module containing the model specification.

2. The name of the database file (produced by mbg-infer) to be used to generate the 
  maps. If you do not want it to go in the current directory, specify a path.

3. The number of burnin iterations to discard from the trace before making the maps.
  You will need to figure this out by inspecting the traces produced by ``mbg-infer``.
 
4. The name of a text file containing one or more multipolygons in geojson format. The
  ``properties`` of each multipolygon must contain a unique ``name`` key. For
  spatiotemporal models, they must also contain ``tmin`` and ``tmax`` keys.
  

Options
-------

* ``-r`` or ``--reps`` : The number of repetitions to do, for purposes of estimating Monte
  Carlo standard error. Defaults to 10.
  
* ``x`` or ``--points`` : The number of points in space or space-time to use to estimate the
  integrals at each repetition. Defaults to 100.

* ``-q`` or ``--quantiles`` : A string containing the quantiles you want. For example,
  ``'0.25 0.5 0.75'`` would produce the lower and upper quartiles and the median of the
  posterior for each areal summary. Default is ``'0.05 0.25 0.5 0.75 0.95'``.

* ``-t`` or ``--thin`` : The factor by which to thin the MCMC trace stored in the database.
  If you use ``-t 10``, only every 10th stored MCMC iteration will be used to produce the 
  estimates. Small values are good but slow. 1 is best. Defaults to 10.
  
* ``-i`` or ``--iter`` : The total number of predictive samples to use in generating the 
  estimates. Large values are good but slow. Defaults to 1000.

* ``-p`` or ``--raster-path`` : The path to the files containing the covariate rasters, if any. 
  These files' headers must match one another, and their missing pixels must match also. There 
  must be a file corresponding to every covariate column in input 3 of mbg-infer. If any of the
  multipolygons extend outside the given rasters or contain missing pixels, an error will result.
  
* ``-w`` or ``--weight-raster`` : The name of a raster file in the ``raster-path``, with no
  extension. See the PDF documentation.

* ``-c`` or ``--coordinate-time`` : If ``1``, sampling points line up in the temporal dimension.
  See the PDF documentation.

* ``-d`` or ``--ignore-npd`` : If ``1``, MCMC iterations whose covariance functions are non-
  positive-definite on the data locations plus the prediction locations will be ignored. If
  ``0``, any such iterations will result in errors.

``mbg-validate``
================
::

    mbg-validate module database-file burn pred-pts [options]
    
mbg-validate produces a folder called ``name-validation``, ``name`` being the name of the database file.
It populates this folder with two csv files called ``p-samps`` and ``n-samps`` containing posterior
predictive samples of the probability of positivity and the number of individuals positive at each 
prediction location.

It also writes three of the four MBG world validation panels into the folder as PDF's.

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
  Small values are good but slow. 1 is best. Defaults to 50.

* ``-i`` or ``--iter`` : The total number of predictive samples you want to generate. Large
  values are good but slow. Defaults to 50000.
  
* ``-d`` or ``--ignore-npd`` : If ``1``, MCMC iterations whose covariance functions are non-
  positive-definite on the data locations plus the prediction locations will be ignored. If
  ``0``, any such iterations will result in errors. Defaults to 0.



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
::

    mbg-realize-prior module ascii0.asc ascii1.asc ... [options]
    
mbg-realize-prior produces a number of prior realizations of the target surface (eg parasite
rate, gene frequency, etc). on several different asciis. Joint or 'conditional' simulations
of surfaces are very expensive, so you can only afford to evaluate them on a few thousand
pixels. 

The multiple asciis are meant to be at multiple resolutions: you can make a coarse one over 
your entire area of interest, a medium-resolution one on a zoomed-in subset, and a few fine 
ones over small areas scattered around. That way you can see the large- and small-scale
properties of the surface allowed by your prior without having to render the entire surface
at full resolution.

Outputs a number of surfaces, evaluated onto the masks indicated by the input asciis. Each set
of realizations is coherent across the input asciis; that is, the 'same' surface is evaluated
on each ascii. That means you can meaningfully overlay the output asciis at different
resolutions.

NOTE: All the parameters of the model will be drawn from the prior before generating each
realization. If you want to fix a variable, you must set its ``observed`` flag.

Required arguments
------------------

1. The name of the module containing the model specification.

2. Several ascii files. Realizations will be evaluated on the union of the unmasked regions
   of these files.
   
Options
-------

* ``-n`` or ``--n-realizations`` : The number of realizations to generate. Defaults to 5.

* ``-m`` or ``--mean`` : The value of the global mean to use. Defaults to 0.

* ``-y`` or ``-year`` : If your model is spatiotemporal, you must provide the decimal year at 
  which you want your realizations produced. For example, Jan 1 2008 would be ``-y 2008``.