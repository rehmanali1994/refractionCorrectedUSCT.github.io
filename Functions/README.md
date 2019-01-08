Prior to running any MATLAB code, you need to build the MEX function “msfm2d” by running the command “mex msfm2d.c” —- NOTE: this only seems to work on MATLAB R2016a or earlier. I have included the compiled MEX function for Mac OS X and Linux, but for other operating systems you have to run “mex msfm2d.c” —- Once again, this “mex msfm2d.c” compilation step may not work in the newer versions of MATLAB.

The MEX code for the fast marching method came from the Accurate Fast Marching Code by Dirk-Jan Kroon: https://www.mathworks.com/matlabcentral/fileexchange/24531-accurate-fast-marching
