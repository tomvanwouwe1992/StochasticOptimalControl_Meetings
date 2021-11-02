clear all; close all; clc;

% 1D pointmass (mass = m) perturbed by a random (Gaussian) force (F)
% dx/dt = xdot
% dxdot/dt = F/m

% We approximate the evolution of the full stochastic system over 3s by performing
% Monte Carlo simulations (nr_sim simulations). We use a simple forward Euler
% integration scheme with a timestep of 10ms.
T = 10; dt = 0.01; N = ceil(T/dt); t = 0:dt:T;
m = 1; nr_sim = 100; plotted_index = ceil(nr_sim/20);
x_sim = NaN(nr_sim,N+1); x_sim(:,1) = 0;
xdot_sim = NaN(nr_sim,N+1); xdot_sim(:,1) = 0;
acc_i = NaN(nr_sim,N); 
F_var_spectral = 7.2; % 1 N².t

options = optimoptions('fmincon');
x = fsolve(@computePdot_KD,[0;0],options,1,[0.2651 0;0 1.6514],F_var_spectral);


K = x(1);
D = x(2);

% Computing mean and covariance propagation based on Lyapunov equations
sigma_W = F_var_spectral;
A = [0 1; K/m D/m];
B = [0; 0];
C = [0; 1/m];
P_0 = zeros(2,2);
P = NaN(2,2,N+1); P(:,:,1) = P_0;
x_mean = NaN(1,N+1); x_mean(1,1) = 0;
xdot_mean = NaN(1,N+1); xdot_mean(1,1) = 0;
for i = 1:N
    state_next = (1+A*dt)*[x_mean(1,i);xdot_mean(1,i)];
    x_mean(1,i+1) = state_next(1); xdot_mean(1,i+1) = state_next(2);
    
    Pdot = A*P(:,:,i) + P(:,:,i)*A' + C*sigma_W*C';
    P(:,:,i+1) = P(:,:,i) + dt*Pdot;
end


figure()
subplot(3,3,1)
ylabel('[m]');
title('position');
subplot(3,3,2)
plot(t,mean(x_sim)); hold on;
plot(t,x_mean,'k');
ylabel('[m]');
title('mean position');
subplot(3,3,3)
plot(t, sqrt(squeeze(P(1,1,:))), 'k');
ylabel('[m]');
title('std position');
subplot(3,3,4)
ylabel('[m/s]');
title('velocity');
subplot(3,3,5)
plot(t,mean(xdot_sim)); hold on;
plot(t,xdot_mean,'k');
ylabel('[m/s]');
title('mean velocity');
subplot(3,3,6)
ylabel('[m/s]');
title('std velocity'); hold on;
plot(t, sqrt(squeeze(P(2,2,:))), 'k');
figure()
plot(t, squeeze(P(1,2,:))./(sqrt(squeeze(P(1,1,:))).*sqrt(squeeze(P(2,2,:)))), 'k');
