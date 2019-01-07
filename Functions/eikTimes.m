function times = eikTimes(xi, zi, C, ind)
%EIKTIMES Calculates Travel Times Between Elements Using Eikonal Equation
% T = eikProjMat(xi, zi, C, ind)
% INPUTS:
%   xi = x computational grid [m]
%   xi = z computational grid [m]
%   C = sound speed [m/s] map over computational grid
%   ind = element locations as linear indices over grid locations
% OUTPUTS:
%   times = numel(ind) x numel(ind) array of travel times

% Create Computational Grid for Tomography
numElements = numel(ind);
dxi = mean(diff(xi)); 

% Build Sparse System Matrix for Straight-Path Forward Projection
times = zeros(numel(ind), numel(ind));
for elmt = 1:numElements
    % Calculate Eikonal Solution and Ray Tracing
    T = eikonalField(xi, zi, C, ind(elmt)); times(:,elmt) = T(ind);
    % Accumulate Sparse Matrix for these Projections
    disp(['Travel Times From Element ', num2str(elmt), ' Completed']);
end

end

