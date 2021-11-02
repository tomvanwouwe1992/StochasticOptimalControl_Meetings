clear all; close all; clc;

% 1D pointmass (mass = m) perturbed by a random (Gaussian) force (F)
% dx/dt = xdot
% dxdot/dt = F/m

% We approximate the evolution of the full stochastic system over 3s by performing
% Monte Carlo simulations (nr_sim simulations). We use a simple forward Euler
% integration scheme with a timestep of 10ms.
T = 10; dt = 0.01; N = ceil(T/dt); t = 0:dt:T;
m = 1; nr_sim = 10000; plotted_index = ceil(nr_sim/20);
x_sim = NaN(nr_sim,N+1); x_sim(:,1) = 0;
xdot_sim = NaN(nr_sim,N+1); xdot_sim(:,1) = 0;
F_var_spectral = 1; % 1 N².t

% Perturbation force drawn from a random distribution. Discretized form the
% (continuous) spectral density.
F = sqrt(F_var_spectral/dt)*randn(nr_sim,N);

% Perform simulation of all samples (Euler)
for i = 1:N
    x_sim(:,i+1) = x_sim(:,i) + xdot_sim(:,i)*dt;
    xdot_sim(:,i+1) = xdot_sim(:,i) + F(:,i)/m*dt;
end

figure()
subplot(3,3,1)
plot(t,x_sim(1:plotted_index:end,:)');
ylabel('[m]');
title('position');
subplot(3,3,2)
plot(t,mean(x_sim));
ylabel('[m]');
title('mean position');
subplot(3,3,3)
plot(t,std(x_sim));
ylabel('[m]');
title('std position');
subplot(3,3,4)
plot(t,xdot_sim(1:plotted_index:end,:)');
ylabel('[m/s]');
title('velocity');
subplot(3,3,5)
plot(t,mean(xdot_sim));
ylabel('[m/s]');
title('mean velocity');
subplot(3,3,6)
plot(t,std(xdot_sim));
ylabel('[m/s]');
title('std velocity');
subplot(3,3,7)
plot(t(1:end-1),F(1:plotted_index:end,:)/m);
ylabel('[m/s²]');
xlabel('time [s]')
title('acceleration');
subplot(3,3,8)
plot(t(1:end-1),mean(F/m));
ylabel('[m/s²]');
xlabel('time [s]')
title('mean acceleration');
subplot(3,3,9)
plot(t(1:end-1),std(F/m));
ylabel('[m/s²]');
xlabel('time [s]')
title('std acceleration');



