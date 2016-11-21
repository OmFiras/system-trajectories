% p = xy = [x0; y0] column vector
% gamma = column vector

function[KE] = kenergy(p,N,gamma)

% vortex position vectors
x_p = p(1:N);           % column vector of x positions
y_p = p(N+1:2*N);       % column vector of y positions

% Create grid of vector coordinates
[xpmesh, ypmesh] = meshgrid(x_p, y_p);

% Subtraction of vortex coordinates, removing duplicates
xdist = triu(xpmesh' - xpmesh,1);
ydist = triu(ypmesh - ypmesh',1);

% Distances squared between vortexes
r = xdist.^2 + ydist.^2;

% matrix of circulations and ignore duplicate products with triu()
gamma_matrix = triu(gamma*gamma', 1);

% compute energy
KE =  (-1/(4*pi))*sum(nonzeros(gamma_matrix).*log(nonzeros(r)));