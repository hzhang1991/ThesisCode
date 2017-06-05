DoGPack plotting routines based on Python's MATPLOTLIB
======================================================

This library contains plotting routines for [DoGPack](../../README.md) that
are written in [Python](https://www.python.org/).

See also: ['$DOGPACK/viz/matlab'](../matlab/README.md) for the 
MATLAB plotting routines.

Basic overview
--------------

To execute a plotting routine, type

    $ python $DOGPACK/viz/python/plotdog[n].py

from an application directory, after running `dog.exe', where [n] = dimension of
application.  

For a full list of options, type

        $ python $DOGPACK/viz/python/plotdog[n].py -h

For example, running


        $ python $DOGPACK/viz/python/plotdog1.py -p 4 -o output3000

executes the 1D plotting script and request four points per element and will
read data from the output directory `output3000'.

User supplied script
--------------------

Each application has a single script, that is called once per each frame that
is plotted.  For example, these scripts have the following format:

    plotq1 - 1D script
    plotq2\_cart - 2D, Cartesian script
    plotq2\_unst - 2D, Unstructured script

These routines allow the user to supply items such as line colors, axes, grid
labels, color bars, etc., everything that is application specific to the given
problem.

In the case of a user not supplying one of these function, a single generic
routine is supplied here in the library.
