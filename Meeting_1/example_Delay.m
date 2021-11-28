clear all; close all; clc;

% 1D pointmass (mass = m) perturbed by a random (Gaussian) force (F)
% dx/dt = xdot
% dxdot/dt = F/m

% We approximate the evolution of the full stochastic system over 3s by performing
% Monte Carlo simulations (nr_sim simulations). We use a simple forward Euler
% integration scheme with a timestep of 10ms.
T = 1; dt = 0.01; N = ceil(T/dt); t = 0:dt:T;
m = 1; nr_sim = 100000; plotted_index = ceil(nr_sim/20);
x_sim = NaN(nr_sim,N+1); x_sim(:,1) = 0;
x_sim_delayed = NaN(nr_sim,N+1); x_sim_delayed(:,1) = 0;
F_var_spectral = 1;
K = -10;
% Perturbation force drawn from a random distribution. Discretized form the
% (continuous) spectral density.
F = sqrt(F_var_spectral/dt)*randn(nr_sim,N);
delay = 5;
for i = 1:N
    vel_i(:,i) = F(:,i)+K*x_sim_delayed(:,i);
    x_sim(:,i+1) = x_sim(:,i) + vel_i(:,i)*dt;
    if i < delay
        x_sim_delayed(:,i+1) = x_sim_delayed(:,i);
    else
        x_sim_delayed(:,i+1) = x_sim(:,i-delay+1);
    end
end

% Computing mean and covariance propagation based on Lyapunov equations
sigma_W = F_var_spectral;
A = [0  K/m ;
     0 0 ];
%  Areg = [ 0 1; 
%           K/m D/m];
B = [0; 0];
C = [1/m; 0];

Cov_x_xd_0 = zeros(1,1);
Cov_x_xd = NaN(1,1,N+1); Cov_x_xd(:,:,1) = Cov_x_xd_0;

P_0 = zeros(1,1);
P = NaN(1,1,N+1); P(:,:,1) = P_0;

Pd_0 = zeros(1,1);
Pd = NaN(1,1,N+1); Pd(:,:,1) = Pd_0;

Pfull_0 = zeros(2,2);
Pfull = NaN(2,2,N+1); Pfull(:,:,1) = Pfull_0;

x_mean = NaN(1,N+1); x_mean(1,1) = 0;
xdot_mean = NaN(1,N+1); xdot_mean(1,1) = 0;
x_delayed_mean = NaN(1,N+1); x_delayed_mean(1,1) = 0;
xdot_delayed_mean = NaN(1,N+1); xdot_delayed_mean(1,1) = 0;

for i = 1:N
%     state_next = (1+Areg*dt)*[x_mean(1,i);xdot_mean(1,i)];
%     x_mean(1,i+1) = state_next(1); xdot_mean(1,i+1) = state_next(2);
    
    Pdot = A*Pfull(:,:,i) + Pfull(:,:,i)*A' + C*sigma_W*C'  ;
    Pfull(:,:,i+1) = Pfull(:,:,i)+ dt*Pdot;

    if i < delay + 1
        x_delayed_mean(:,i+1) = x_delayed_mean(:,i);
        xdot_delayed_mean(:,i+1) = xdot_delayed_mean(:,i);
        Pd(:,:,i+1) = Pd(:,:,i);
        Pfull(2,2,i+1) = Pd(:,:,i+1);
    else
        x_delayed_mean(:,i+1) = x_mean(:,i-delay+1);
        xdot_delayed_mean(:,i+1) = xdot_mean(:,i-delay+1);
        Pd(:,:,i+1) = Pfull(1,1,i-delay+1);
        Pfull(2,2,i+1) = Pd(:,:,i+1);
    end
    
end


figure()
subplot(3,3,1)
plot(t,x_sim(1:plotted_index:end,:)');
ylabel('[m]');
title('position');
subplot(3,3,2)
plot(t,mean(x_sim)); hold on;
% plot(t,x_mean,'k');
ylabel('[m]');
title('mean position');
subplot(3,3,3)
plot(t,std(x_sim)); hold on;
plot(t, sqrt(squeeze(Pfull(1,1,:))), 'k');
ylabel('[m]');
title('std position');



for i = 1:N+1
P_sim(:,:,i) = cov(x_sim(:,i),xdot_sim(:,i));
end
figure()
plot(t, squeeze(P_sim(1,2,:))./(sqrt(squeeze(P_sim(1,1,:))).*sqrt(squeeze(P_sim(2,2,:))))); hold on;
plot(t, squeeze(P(1,2,:))./(sqrt(squeeze(P(1,1,:))).*sqrt(squeeze(P(2,2,:)))), 'k');
