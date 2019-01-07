function [path_x, path_y] = rayTracEikonDescnd(x, y, t, xdst, ydst, xsrc, ysrc, ds)
%RAYTRACEIKONDESCND Ray-Tracing Descent on Eikonal Calculated Arrival Times
% paths = rayTracEikonDescnd(x, y, t, xend, yend, xs, ys)
%   x, y -- (x, y) grid on which arrival times are given [m]
%   t -- arrival times [s]
%   xdst, ydst -- (x, y) start positions (vectors) 
%   xsrc, ysrc -- (x, y) src position
%   ds -- arc-length step size [m]

% Arrival Time Surface Grid: Assumes Uniform Sampling
x_start = x(1); y_start = y(1);
dx = x(2)-x(1); dy = y(2)-y(1);
Nx = numel(x); Ny = numel(y);

% Ray-Length Pixel Intersection Matrix
numElements = numel(xdst);

% Arrival Time Gradient
[delT_x, delT_y] = gradient(t);

% Maximum Number of Points in Each Ray Tracing Path
maxNumPoints = ceil(((max(x)-min(x))+(max(y)-min(y))^2)/ds);

% Initialize Paths
path_x = NaN*ones(numElements, maxNumPoints);
path_y = NaN*ones(numElements, maxNumPoints);
i = 1; path_x(:,i) = xdst; path_y(:,i) = ydst; 

while any(sqrt((path_x(:,i)-xsrc).^2 + (path_y(:,i)-ysrc).^2)>ds) && i < maxNumPoints
    % Indices for Paths that Have Not Terminated Yet
    pidx = (sqrt((path_x(:,i)-xsrc).^2 + (path_y(:,i)-ysrc).^2)>ds);
    
    % Convert Path Locations to Grid Index Coordinates
    path_x_idx = min(Nx, max(1, 1+(path_x(pidx,i)-x_start)/dx));
    path_y_idx = min(Ny, max(1, 1+(path_y(pidx,i)-y_start)/dy));
    
    % Bilinear Interpolation of Gradient Using Above Grid Index Coordinates
    dx1 = abs(path_x_idx - floor(path_x_idx));
    dx2 = abs(floor(path_x_idx) - path_x_idx + 1);
    dy1 = abs(path_y_idx - floor(path_y_idx));
    dy2 = abs(floor(path_y_idx) - path_y_idx + 1);
    w11 = dx2.*dy2; w12 = dx2.*dy1;
    w21 = dx1.*dy2; w22 = dx1.*dy1;
    idx11 = sub2ind([Ny, Nx], floor(path_y_idx), floor(path_x_idx));
    idx12 = sub2ind([Ny, Nx], ceil(path_y_idx), floor(path_x_idx));
    idx21 = sub2ind([Ny, Nx], floor(path_y_idx), ceil(path_x_idx));
    idx22 = sub2ind([Ny, Nx], ceil(path_y_idx), ceil(path_x_idx));
    ds_x = w11.*delT_x(idx11) + w12.*delT_x(idx12) + w21.*delT_x(idx21) + w22.*delT_x(idx22);
    ds_y = w11.*delT_y(idx11) + w12.*delT_y(idx12) + w21.*delT_y(idx21) + w22.*delT_y(idx22);
        
    % Get Descent Direction
    dr = [ds_y, ds_x]; dl = sqrt(sum(dr.^2, 2));
    dr(:,1) = dr(:,1)./dl; dr(:,2) = dr(:,2)./dl;
    
    % Descend
    path_y(pidx,i+1) = path_y(pidx,i) - ds*dr(:,1);
    path_x(pidx,i+1) = path_x(pidx,i) - ds*dr(:,2); 
    i = i+1;
end

% Transpose for Output
path_x = path_x';
path_y = path_y';

end

