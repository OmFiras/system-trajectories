% Computes angular momentum of vortexes

function AM = angularmomentum(p, N, gamma)
% p = [xtraj(i,:)'; ytraj(i,:)']     p is a column vector

% vortex position vectors
xp = p(1:N);            % column vector of x positions
yp = p(N+1:2*N);        % column vector of y positions

AM = sum(gamma'.*(xp'.^2 + yp'.^2));