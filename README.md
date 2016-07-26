# PLS-to-Python (pls2py)

A collection of scripts to facillitate moving PLS results out of MATLAB and into Python.

As of now, these tools can be used to:

1. Create temporal brainscore plots from event-related PLS analyses
2. Create seed-based functional connectivity matrices for the output of multiple-voxel extraction in an event-related PLS analysis.
3. Project brainscores from significant latent variables in event-related or blocked PLS analyses. 
4. Read in PLS batch txt files and manipulate timing and/or extract onsets for AFNI (or other programs of interest).

All scripts are designed to work with PLSGUI v5.1206281 results derived in MATLAB R2012b. Some have scripts have been tested for compatibility with PLS v6.1311050, but mileage may vary!

These scripts will hopefully become obsolete once PLS transitions to Python, expected later this year (2016). 

If you have any questions, please open a pull request! 