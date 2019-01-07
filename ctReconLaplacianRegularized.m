clear
clc

% Add Functions to Path
addpath(genpath('Functions'));

% Load Ring CT Data
load('travelTimes.mat');

% Show Sound Speed Map with Ring Transducers
figure; imagesc(xi, zi, C); 
axis image; colormap gray; colorbar;
xlabel('X Coordinate [m]'); 
ylabel('Z Coordinate [m]'); 
title('Sound Speed Map');
hold on; plot(x_circ, z_circ, 'w.');

% Simulated Travel Times Using the Eikonal Equation
times = times(:); % Add Random Measurement Errors
times = times + (1e-8)*randn(size(times)); 

%% Solve Using Laplacian Regularization Approach

% Reconstruction Parameters
c_guess = 1540; % Uniform Sound Speed Guess
m = (1/c_guess)*ones(Nzi*Nxi,1); % Running Reconstruction
M = numElements*numElements; % Number of Observations
d = [times; zeros(Nxi*Nzi,1)]; % Observations Vector

% Laplacian Matrix for Regularization
mu_reg = 1e-2; 
L = mu_reg*laplacianOperator(Nxi, Nzi);

% Gauss-Newton Method
numGaussNewtonIterations = 4;
for iter_gn = 1:numGaussNewtonIterations
    % Matrix Encoding Projections Along Integration Paths
    C_current = reshape(1./m, [Nzi, Nxi]);
    H = eikProjMat(xi, zi, C_current, ind);
    % Forward Model
    G = @(m) [H*m; L*m];
    GT = @(s) H'*s(1:M) + L'*s(M+1:end);
    % Initial Conditions
    p = zeros(size(m)); % Conjugate Direction Vector
    beta = 0; % Variable for Updating Conjugate Direction
    s = G(m)-d; % Current Residual Vector
    r = GT(s); % Current Gradient Direction
    % Iterative Reconstruction by Conjugate-Gradient Least-Squares
    figure; numIter = 400; % Total Number of Iterations
    for iter = 1:numIter
        % CGLS Updates
        p = -r + beta*p;
        Gp = G(p); r_norm_sq_last = r'*r;
        alpha = r_norm_sq_last/(Gp'*Gp);
        m = m + alpha*p;
        s = s + alpha*Gp;
        r = GT(s);
        beta = (r'*r)/r_norm_sq_last;  
        % Plot Result of Each Iteration
        imagesc(xi, zi, reshape(1./m, [Nzi, Nxi])); 
        axis image; colormap gray; colorbar(); caxis([min(C(:)), max(C(:))]);
        xlabel('X coordinate'); ylabel('Y coordinate'); 
        title(['Conjugate Gradient Reconstruction Iteration: ', num2str(iter+1)]);
        getframe;
    end
    % Show Final Result
    imagesc(xi, zi, reshape(1./m, [Nzi, Nxi]));
    xlabel('X Coordinate [m]'); ylabel('Z Coordinate [m]'); 
    axis image; colormap gray; colorbar; caxis([min(C(:)), max(C(:))]);
    hold on; plot(x_circ, z_circ, 'w.');
    title('Sound Speed Reconstruction by Conjugate Gradient');
end