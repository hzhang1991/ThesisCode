DoGPack plotting routines based on the MATrix LABoratory (MATLAB)
=================================================================

This library contains plotting routines for [DoGPack](../../README.md) that
are written in [MATLAB](http://www.mathworks.com/).  Many of these work
may work out of the box with the open-source alternative to MATLAB called 
[Octave](http://www.gnu.org/software/octave/), but the developers do not keep
up to date with this option.

See also: ['$DOGPACK/viz/python'](../python/README.md) for a complete set of 
[Python](http://www.python.org/)
plotting routines based on [matplotlib](http://matplotlib.org/).

Basic overview
--------------

The basic calling sequence for any applicaiton is based on calling one of the
following functions from the application folder:

plotdog1(), plotdog2(), plotdog3()

for 1, 2 and 3D applications.  The basic format is the following:

    plotdog[n]( points\_per\_dir, outputdir\_in, point\_type),

where points\_per\_dir = number of points in each direction, output\_dir is the
location of the output directory, and point\_type = 1, or 2 (Uniform or
Gaussian).  Without any arguments, this call defaults to 1 point per
direction, output\_dir\_in = 'output', and point\_type = 1.

User supplied script
--------------------

Each application has a single script, that is called once per each frame that
is plotted.  For example, these scripts have the following format:

    plotq1 - 1D script
    plotq2\_cart - 2D, Cartesian script
    plotq2\_unst - 2D, Unstructured script
    plotq3\_cart - 3D, Cartesian script
    plotq3\_cart\_slices - 3D, Cartesian script for plotting slices

These routines allow the user to supply items such as line colors, axes, grid
labels, color bars, etc., everything that is application specific to the given
problem.

In the case of a user not supplying one of these function, a single generic
routine is supplied here in the library.
