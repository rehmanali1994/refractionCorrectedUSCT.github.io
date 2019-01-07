clear
clc

% Add Functions to Path
addpath(genpath('Functions'));

% Sound Speed Map 
[x, z, c, c_bkgnd] = soundSpeedPhantom();
[X, Z] = meshgrid(x,z); 
dxi = 0.12e-3; xmax = 60e-3;
xi = -xmax:dxi:xmax; zi = xi;
Nxi = numel(xi); Nzi = numel(zi);
[Xi, Zi] = meshgrid(xi, zi);
R = sqrt(Xi.^2 + Zi.^2); 
rotAngle = 2.85*pi; T = atan2(Zi, Xi) - rotAngle; 
C = interp2(X, Z, c, R.*cos(T), R.*sin(T), 'linear', c_bkgnd);

% Create Transducer Ring
circle_radius = 55e-3; numElements = 128;
circle_rad_pixels = floor(circle_radius/dxi);
theta = -pi:2*pi/numElements:pi-2*pi/numElements;
x_circ = circle_radius*cos(theta); 
z_circ = circle_radius*sin(theta); 
[x_idx, z_idx, ind] = sampled_circle(Nxi, Nzi, circle_rad_pixels, theta);
msk = zeros(Nzi, Nxi); msk(ind) = 1;

% Show Sound Speed Map with Ring Transducers
figure; imagesc(xi, zi, C+50*msk); 
axis image; colormap gray; colorbar; 
hold on; plot(x_circ, z_circ, 'w.'); 
xlabel('X Coordinate [m]'); 
ylabel('Z Coordinate [m]'); 

% Build Sparse System Matrix for Straight-Path Forward Projection
times = eikTimes(xi, zi, C, ind);

% Save Results to File
save('travelTimes.mat');
