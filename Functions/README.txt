1) Prior to running any MATLAB code, you need to build the MEX function “msfm2d” by running the command “mex msfm2d.c” —- NOTE: this only works on MATLAB R2016a or earlier. I have included the compiled MEX function for Mac OS X, but for other operating systems you have to run “mex msfm2d.c” —- Once again, this “mex msfm2d.c” compilation step may not work in the newer versions of MATLAB.

2) To simulate the travel times in this ring, run the “ringct.m” script. This will save to a MAT file called “ringct.mat”

3) Run “ctReconBayesianCG.m” to see the reconstruction of the sound speed in the medium based on the travel times.


