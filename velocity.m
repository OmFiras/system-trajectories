% Credit for the general method used in this approach belongs to Dr. Tom Ashbee.

function xyvel = velocity(t,p,N,gamma)

% Create grid of vector coordinates
[x0mesh, y0mesh] = meshgrid(p(1:N,1), p(N+1:2*N,1));

% Subtraction of vortex coordinates respectively
xdis = triu(x0mesh' - x0mesh,0);
ydis = triu(y0mesh - y0mesh',0);

% Distances squared between vortexes,
r = triu(xdis.^2 + ydis.^2,-1)' + (xdis.^2 + ydis.^2);

% Redefine gamma vector to proper matrix dimensions
gamma2 = repmat(gamma',N,1);

% Compute some terms early, eliminate NaN
xvelpre = zeros(size(r));
yvelpre = zeros(size(r));
n = r~=0;
tempy = (-triu(ydis,-1)'+ydis);
tempx = (-triu(xdis,-1)'+xdis);
xvelpre(n) = gamma2(n).*tempy(n)./r(n);
yvelpre(n) = gamma2(n).*tempx(n)./r(n);

% Velocities of vortexes
xvel = -(1/(2*pi))*sum(xvelpre,2);
yvel = (1/(2*pi))*sum(yvelpre,2);

% Vortex velocities in a single column
xyvel = [xvel; yvel];