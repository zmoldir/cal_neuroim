This implementation is intended for the deconvolution of calcium neuro-imaging data.
It can either be used as a runnable script as-is, or used in parts by importing from the API,
which it is originally intended for.

_______________________________________________________________________________________________________
See the section "running the script" of the PDF contained in this git for a more elaborate explanation!
_______________________________________________________________________________________________________

The distribution consists of several files:

main.py - a short script which deconvolves input traces and returns spike trains. 
  The required syntax is printed if called without it
  
cal_neuroIm.py - contains the methods used in the script, which it is imported by. 
  Can be used as a library for new implementations.
  
MAstuff.ipynb - ipython notebook with multiple tests, data simulation and and read, 
  as well as visualization routines.
  
REQUIREMENTS: Pandas, SciPy, NumPy, Python 2.7 or above

I strongly recommend taking a look at the PDF if you intend to use this script.
