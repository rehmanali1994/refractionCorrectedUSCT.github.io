function H = eikProjMat(xi, zi, C, ind)
%EIKPROJMAT Creates Sparse Matrix for Projections Over Refracted Paths
% H = eikProjMat(xi, zi, C, ind)
% INPUTS:
%   xi = x computational grid [m]
%   xi = z computational grid [m]
%   C = sound speed [m/s] map over computational grid
%   ind = element locations as linear indices over grid locations
% OUTPUTS:
%   H = sparse matrix used to calculate travel times from slowness vector

% Create Computational Grid for Tomography
[Xi, Zi] = meshgrid(xi, zi); 
Nxi = numel(xi); Nzi = numel(zi);
numElements = numel(ind);
dxi = mean(diff(xi)); 

% Build Sparse System Matrix for Straight-Path Forward Projection
H = sparse(numElements*numElements, Nzi*Nxi);
ds = dxi/10; % Step Size for Ray Tracing
for elmt = 1:numElements
    % Calculate Eikonal Solution and Ray Tracing
    T = eikonalField(xi, zi, C, ind(elmt));
    tic; [path_x, path_z] = rayTracEikonDescnd(xi, zi, T, ...
        Xi(ind), Zi(ind), Xi(ind(elmt)), Zi(ind(elmt)), ds); toc;
    % Accumulate Sparse Matrix for these Projections
    tic; path_x_idx = round(min(Nxi, max(1, 1+(path_x-xi(1))/dxi)));
    path_z_idx = round(min(Nzi, max(1, 1+(path_z-zi(1))/dxi)));
    path_x_idx(isnan(path_x)) = NaN;
    path_z_idx(isnan(path_z)) = NaN;
    path_ind = sub2ind([Nzi, Nxi], path_z_idx, path_x_idx);
    elmt_ind = ones(size(path_ind,1),1)*(1:numElements)+(elmt-1)*numElements;
    elmt_ind_flat = elmt_ind(~isnan(path_ind));
    path_ind_flat = path_ind(~isnan(path_ind));
    H = H + ds*sparse(elmt_ind_flat, path_ind_flat, ...
        ones(numel(elmt_ind_flat),1), ...
        numElements*numElements, Nzi*Nxi); toc;
    disp(['Sparse Matrix Tx Element ', num2str(elmt), ' Completed']);
end

end

