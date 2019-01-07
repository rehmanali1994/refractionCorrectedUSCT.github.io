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

%% Solve Using Bayesian Approach

% Dimensions of H
M = numElements*numElements; 
N = Nzi*Nxi;

% Prior Covariance Matrix Q as a Linear Operator
dr = 2*pi*circle_radius/numElements;
blurKernelLen = 0.5*dr/dxi; % Blurring Kernel Length
Q = @(x) reshape(imgaussfilt(reshape(x,[Nzi,Nxi]),blurKernelLen),[N,1]);

% Observation Covariance Matrix R as Sparse Matrix
R = (1e-7)*speye(M);

% Prior Mean Vector
c_guess = 1540; % Uniform Sound Speed Guess
mu = (1/c_guess)*ones(Nzi*Nxi,1); % Running Reconstruction

% Gauss-Newton Method
numGaussNewtonIterations = 4;
for iter_gn = 1:numGaussNewtonIterations
    % Matrix Encoding Projections Along Integration Paths
    C_current = reshape(1./mu, [Nzi, Nxi]);
    H = eikProjMat(xi, zi, C_current, ind);
    % Define A (Symmetric) Matrix as Linear Operator
    A = @(x) H*(Q(H'*x))+R*x;
    % Calculate Vector b
    y = times;
    b = y - H*mu;
    % Initial Conditions for Conjugate Gradient Algorithm
    x = zeros(M, 1); % Initial Solution
    beta = 0; % Variable for Updating Conjugate Direction
    p = zeros(size(x)); % Conjugate Direction Vector
    r = A(x)-b; % Current Gradient Direction
    % Iterative Reconstruction: Use Conjugate Gradient to Solve Ax = b
    figure; numIter = 100;
    for iter = 1:numIter
        % Conjugate Gradient Updates
        p = -r + beta*p;
        r_norm_sq_last = r'*r; Ap = A(p);
        alpha = r_norm_sq_last/(p'*Ap);
        x = x + alpha*p;
        r = r + alpha*Ap;
        beta = (r'*r)/r_norm_sq_last;  
        % Calculate Posterior Mean (Reconstructed Image
        mupost = mu + Q(H'*x);
        % Plot Result of Each Iteration
        imagesc(xi, zi, reshape(1./mupost, [Nzi, Nxi])); 
        axis image; colormap gray; colorbar(); caxis([min(C(:)), max(C(:))]);
        xlabel('X coordinate'); ylabel('Y coordinate'); 
        title(['Conjugate Gradient Reconstruction Iteration: ', num2str(iter+1)]);
        getframe;
    end
    % Show Final Result for this Gauss-Newton Iteration
    imagesc(xi, zi, reshape(1./mupost, [Nzi, Nxi]));
    xlabel('X Coordinate [m]'); ylabel('Z Coordinate [m]'); 
    axis image; colormap gray; colorbar; caxis([min(C(:)), max(C(:))]);
    hold on; plot(x_circ, z_circ, 'w.');
    title('Sound Speed Reconstruction by Conjugate Gradient');
    % Set Prior Mean for Next Iteration to Posterior from Current Iteration
    mu = mupost;
end
