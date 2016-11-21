% Credit for the general method used in this approach belongs to Dr. Tom Ashbee.

tic
% Initial vortex values
N = 4;                      % number of vortices
gamma = [-1;2;-3;4];        % column vector of circulations of each vortex 
% gamma = 10*rand(N,1);     % for randomized circulations


% Initial vortex coordinates on a unit circle
angle = [0;pi/2;pi;3*pi/2];        % column vector of angles between 0 and 2pi radians
x0 = cos(angle);            % column vector of x-coordinates of vortex
y0 = sin(angle);            % column vector of y-coordinates of vortex
xy = [x0; y0];              % input for "p" (position) in velocity function
% angle = 2*pi*rand(N,1);   % for randomized angles


% Generalized timestep parameters
tmax = 100;             % total time
dt = 0.001;              % interval step, 10^-5 crashes
tspan = 0:dt:tmax;      % row time vector


% Runge-Kutta integration
[t,p] = ode45(@velocity, tspan, xy,[],N,gamma);
% First N columns are x-coordinates, next N columns are y-coordinates for
% each vortex respectively


% Trajectories given by x and y positions
xtraj = p(:,1:N);
ytraj = p(:,N+1:2*N);


% Initializing momentum and energy functions
init_angularmomentum = angularmomentum(xy, N, gamma);
init_xlinearmomentum = xlinearmomentum(xy, N, gamma);
init_ylinearmomentum = ylinearmomentum(xy, N, gamma);
init_kenergy = kenergy(xy, N, gamma);


% Vectors demonstrating conservation
cons_angularmomentum = zeros(numel(tspan),1);
cons_xlinearmomentum = zeros(numel(tspan),1);
cons_ylinearmomentum = zeros(numel(tspan),1);
cons_kenergy = zeros(numel(tspan),1);

for i = 1:numel(tspan)
    cons_angularmomentum(i) = abs(init_angularmomentum - angularmomentum([xtraj(i,:)'; ytraj(i,:)'], N, gamma));
    cons_xlinearmomentum(i) = abs(init_xlinearmomentum - xlinearmomentum([xtraj(i,:)'; ytraj(i,:)'], N, gamma));
    cons_ylinearmomentum(i) = abs(init_ylinearmomentum - ylinearmomentum([xtraj(i,:)'; ytraj(i,:)'], N, gamma));
    cons_kenergy(i) = abs(init_kenergy - kenergy([xtraj(i,:)'; ytraj(i,:)'], N, gamma));
end

% Plots of trajectories
figure
plot(x0,y0,'ko')            % circles indicate initial coordinates
hold on
plot(xtraj,ytraj)
hold on
xlabel('x(t)');
ylabel('y(t)');
title('Interaction of Vortex Trajectories');
axis square

% Plots of absolute error
figure
subplot(2,2,1)
semilogy(tspan, cons_angularmomentum, 'g')
xlabel('time')
title('Absolute error consv angular momentum')
subplot(2,2,2)
semilogy(tspan, cons_kenergy, 'k')
xlabel('time')
title('Absolute error consv energy')
subplot(2,2,3)
semilogy(tspan, cons_xlinearmomentum, 'r')
xlabel('time')
title('Absolute error consv x-linear momentum')
subplot(2,2,4)
semilogy(tspan, cons_ylinearmomentum, 'b')
xlabel('time')
title('Absolute error consv y-linear momentum')
toc
