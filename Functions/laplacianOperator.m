function L = laplacianOperator(Nx, Nz)
%LAPLACIANOPERATOR Del-Squared Operator Over Grid
% L = laplacianOperator(Nx, Nz)
%   Nx = Number of X Positions on Grid
%   Nz = Number of Z Positions on Grid
%   L = Laplacian Operator (as Sparse Matrix) Over Grid

% Linear Indexing Function
lindex = @(r,c) r+(c-1)*Nz;

% Reflecting Constraint in Indices
refl = @(mindex,maxdex,v) ...
    maxdex-abs((abs(v-mindex)+mindex)-maxdex);

% Row and Column Coordinates
r = 1:Nz; c = 1:Nx;
[C, R] = meshgrid(c, r);
C = C(:); R = R(:);

% Construct Laplacian as Sparse Matrix
L = sparse(lindex(R,C), lindex(R,C), ...
    -4*ones(Nx*Nz,1), Nx*Nz, Nx*Nz);
L = L + sparse(lindex(R,C), ...
    lindex(R,refl(1,Nx,C-1)), ...
    ones(Nx*Nz,1), Nx*Nz, Nx*Nz);
L = L + sparse(lindex(R,C), ...
    lindex(R,refl(1,Nx,C+1)), ...
    ones(Nx*Nz,1), Nx*Nz, Nx*Nz);
L = L + sparse(lindex(R,C), ...
    lindex(refl(1,Nz,R-1),C), ...
    ones(Nx*Nz,1), Nx*Nz, Nx*Nz);
L = L + sparse(lindex(R,C), ...
    lindex(refl(1,Nz,R+1),C), ...
    ones(Nx*Nz,1), Nx*Nz, Nx*Nz);

end

