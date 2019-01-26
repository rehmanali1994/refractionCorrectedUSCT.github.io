Go to Function folder and run “mex -compatibleArrayDims msfm2d.c” for Fast-Marching Method code

Run travelTimeSimulation.m to simulate travel time dataset.
Dataset is saved to a file called travelTimes.mat

All reconstruction scripts use a Gauss-Newton method to reconstruct the speed of sound in the medium from travel times between elements on the ring. There are three reconstruction scripts for each of three variations on the basic approach: 1) Laplacian regularization (ctReconLaplacianRegularized.m), 2) the Bayesian approach (ctReconBayesian.m), and 3) the resolution-filling technique (ctReconResolutionFilling.m).

ENJOY! 

