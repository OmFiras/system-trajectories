% Computes linear momentum of vortexes

function XLM = xlinearmomentum(p, N, gamma)
% p = [xtraj(i,:)'; ytraj(i,:)'] = column vector

% column vector of x positions
xpp = p(1:N);

% sum of linear momentum in x direction
XLM = sum(gamma'.*(xpp'));