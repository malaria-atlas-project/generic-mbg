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

===========================
Detailed usage instructions
===========================


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
     
     * ``t`` : Time in years since 2009. This is only required for spatiotemporal models.

   All other columns are interpreted as covariates, eg ``ndvi`` etc.
   

Options
-------

* ``-t`` or ``--thin`` :

* ``-i`` or ``--iter`` :

* ``-n`` or ``-ncpus`` :

* ``-d`` or ``--delay`` :



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

* ``-n`` or ``--n-bins`` :

* ``-b`` or ``--bufsize`` : 

* ``-q`` or ``--quantiles`` : 

* ``-r`` or ``--raster-thin`` :

* ``-t`` or ``--thin`` :

* ``-i`` or ``--iter`` :

* ``-c`` or ``--covariates`` :

* ``-y`` or ``--year`` :


``mbg-validate``
================
::

    mbg-validate module database-file burn pred-pts [options]

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

* ``-t`` or ``--thin`` :

* ``-i`` or ``--iter`` :


``mbg-realize-prior``
=====================

This one is not implemented yet.


``mbg-scalar-priors``
=====================
::

    mbg-scalar-priors module [options]

Required arguments
------------------

1. The name of the module containing the model specification.

Options
-------

* ``-i`` or ``--iter`` :


===================
Module requirements
===================


===========================================
Questions asked with high prior probability
===========================================

* Q: Can you make this work on Windows?
  A: No.