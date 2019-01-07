function T = eikonalField(x, z, c, elmts)
%EIKONALFIELD Eikonal Solution for Each Ring Element
% T = eikonalField(x, z, c, elmts)
%   x, z -- Lateral and Axial Spatial Grids [m]
%   c -- Sound Speed on Spatial Grid [m/s]
%   elmts -- Locations of Transducer Elements on Grid [Linear Indexing]
%   T -- Eikonal Surfaces for Each Element

% Spatial Grid
dx = mean(diff(x)); 
dz = mean(diff(z));
dr = mean([dx, dz]);
Nx = numel(x); x_idx = 1:Nx;
Nz = numel(z); z_idx = 1:Nz;
[X_IDX, Z_IDX] = meshgrid(x_idx, z_idx);

% Transmission Sources 
tx_srcs = [Z_IDX(elmts)'; X_IDX(elmts)'];

% Arrival Times (Number of Sources x Number of Receivers Array)
T = zeros(numel(z), numel(x), numel(elmts));

% Get Arrival Times After Iterating Over Transmitting Elements
for elmt_idx = 1:numel(elmts)
    % Get Arrival Times Using The Eikonal Equation
    T(:, :, elmt_idx) = dr*msfm2d(c, tx_srcs(:, elmt_idx), true, true);
end

end

