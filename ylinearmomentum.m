% Computes linear momentum of vortexes

function YLM = ylinearmomentum(p, N, gamma)
% p = [xtraj(i,:)'; ytraj(i,:)'] = column vector

% column vector of y positions
ypp = p(1:N);

% sum of linear momentum in x direction
YLM = sum(gamma'.*(ypp'));