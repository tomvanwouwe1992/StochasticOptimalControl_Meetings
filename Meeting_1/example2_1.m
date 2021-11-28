clear all; close all; clc;

% 1D pointmass (mass = m) perturbed by a random (Gaussian) force (F)
% dx/dt = xdot
% dxdot/dt = F/m

% We approximate the evolution of the full stochastic system over 3s by performing
% Monte Carlo simulations (nr_sim simulations). We use a simple forward Euler
% integration scheme with a timestep of 10ms.
T = 10; dt = 0.01; N = ceil(T/dt); t = 0:dt:T;
% m = 1; nr_sim = 10000; plotted_index = ceil(nr_sim/20);
% x_sim = NaN(nr_sim,N+1); x_sim(:,1) = 0;
% xdot_sim = NaN(nr_sim,N+1); xdot_sim(:,1) = 0;
% acc_i = NaN(nr_sim,N); 
% K = -10;
% D = -2*sqrt(abs(K));
% % Perturbation force drawn from a random distribution. Discretized form the
% % (continuous) spectral density.
% F = sqrt(F_var_spectral/dt)*randn(nr_sim,N);
% for i = 1:N
%     x_sim(:,i+1) = x_sim(:,i) + xdot_sim(:,i)*dt;
%     acc_i(:,i) = (F(:,i)+K*x_sim(:,i)+D*xdot_sim(:,i))/m;
%     xdot_sim(:,i+1) = xdot_sim(:,i) + (acc_i(:,i))*dt;
% end

% Computing mean and covariance propagation based on Lyapunov equations
F_var_spectral = 1; % 1 N².t
sigma_W = F_var_spectral;
g = 9.81; l = 1; m = 70;
C = [0; 1/(m*l^2)];
P_0 = zeros(2,2);
P = NaN(2,2,N+1); P(:,:,1) = P_0;
x_mean = NaN(1,N+1); x_mean(1,1) = pi-pi/8;
xdot_mean = NaN(1,N+1); xdot_mean(1,1) = 0;
for i = 1:N
    state_next = (1+A*dt)*[x_mean(1,i);xdot_mean(1,i)];
    x_mean(1,i+1) = state_next(1); xdot_mean(1,i+1) = state_next(2);
    A = [0 1; 
        cos(x_mean(1,i))*g/l 0];
    Pdot = A*P(:,:,i) + P(:,:,i)*A' + C*sigma_W*C';
    P(:,:,i+1) = P(:,:,i) + dt*Pdot;
end


figure()
subplot(3,3,1)
plot(t,x_sim(1:plotted_index:end,:)');
ylabel('[m]');
title('position');
subplot(3,3,2)
plot(t,mean(x_sim)); hold on;
plot(t,x_mean,'k');
ylabel('[m]');
title('mean position');
subplot(3,3,3)
plot(t,std(x_sim)); hold on;
plot(t, sqrt(squeeze(P(1,1,:))), 'k');
ylabel('[m]');
title('std position');
subplot(3,3,4)
plot(t,xdot_sim(1:plotted_index:end,:)');
ylabel('[m/s]');
title('velocity');
subplot(3,3,5)
plot(t,mean(xdot_sim)); hold on;
plot(t,xdot_mean,'k');
ylabel('[m/s]');
title('mean velocity');
subplot(3,3,6)
plot(t,std(xdot_sim));
ylabel('[m/s]');
title('std velocity'); hold on;
plot(t, sqrt(squeeze(P(2,2,:))), 'k');
subplot(3,3,7)
plot(t(1:end-1),acc_i(1:plotted_index:end,:)');
ylabel('[m/s²]');
xlabel('time [s]')
title('acceleration');
subplot(3,3,8)
plot(t(1:end-1),mean(acc_i));
ylabel('[m/s²]');
xlabel('time [s]')
title('mean acceleration');
subplot(3,3,9)
plot(t(1:end-1),std(acc_i));
ylabel('[m/s²]');
xlabel('time [s]')
title('std acceleration');

for i = 1:N+1
P_sim(:,:,i) = cov(x_sim(:,i),xdot_sim(:,i));
end
figure()
plot(t, squeeze(P_sim(1,2,:))./(sqrt(squeeze(P_sim(1,1,:))).*sqrt(squeeze(P_sim(2,2,:))))); hold on;
plot(t, squeeze(P(1,2,:))./(sqrt(squeeze(P(1,1,:))).*sqrt(squeeze(P(2,2,:)))), 'k');
