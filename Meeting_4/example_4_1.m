clear all; clc; close all;

import casadi.*
dt = 0.167; T = 10.02; t = 0:dt:T;
N = round(T/dt);



%- Stochastic Dynamics of the controlled system
% Variables initialisaiton
x_mean = MX.sym('x_mean',2);
x_mean_dot = MX.sym('x_mean_dot',2);
w =  MX.sym('w',1); % disturbing force
F =  MX.sym('F',1);
K = MX.sym('X',1,1);
B = MX.sym('X',1,1);

A = MX.sym('A',2,2);
C = MX.sym('A',2,1);

% Dynamic System parameters
g = 9.81; l = 1.2; m = 5;

% Imlicit dynamics of stochastic problem ( xdot = f(x,u,...) )
implicitDyn = x_mean_dot - [ x_mean(2); (-m*g*l*sin(x_mean(1)) + F + w + K*(x_mean(1)-pi) + B*x_mean(2))/(m*l^2)];
fcn_dyn = Function('fcn_dyn',{x_mean_dot,x_mean,F,w,K,B},{implicitDyn});

% Linearization of dynamics around mean state - IMPLICIT CONSTRAINTS TO
% COMPUTE A AND C FROM
A_implicit_Error = jacobian(implicitDyn,x_mean) + jacobian(implicitDyn,x_mean_dot)*A;
fcn_A_implicit = Function('fcn_A_implicit',{x_mean_dot,x_mean,F,w,K,B,A},{A_implicit_Error});

C_implicit_Error = jacobian(implicitDyn,w) + jacobian(implicitDyn,x_mean_dot)*C;
fcn_C_implicit = Function('fcn_C_implicit',{x_mean,F,w,K,B},{C_implicit_Error});

clear x_mean F K B w

% Plot variables
figure();
cs = linspecer(8,'blue');
cs = cs(2:end,:);

% OCP for seven different noise levels
for noise_factor = 1:7
    
    sigma_w = 10^(noise_factor-3);

    opti = casadi.Opti(); % Create opti instance
    
    F = opti.variable(1,N+1);      % force/torque
    x_mean = opti.variable(2,N+1); % mean trajectory
    Pvar = opti.variable(4,N+1);   % covariance matrix
    K = opti.variable(1,1);        % stiffness
    D = opti.variable(1,1);        % damping
    opti.subject_to(K == 0);
    opti.subject_to(D == 0);

    %%Define variables A and C
    
    P_k = 0.0001*eye(2);        
    opti.subject_to(Pvar(:,1) == P_k(:));

    % expected effort
    obj = (F(:,1) + K*(x_mean(1,1)-pi) + D*x_mean(2,1))^2 + [K D]*P_k*[K;D];
    for i = 1:N
        % integration of mean state trajectory
        x_mean_dot = fcn_dyn(x_mean(:,i),F(i),0,K,D);
        x_mean_dot_next = fcn_dyn(x_mean(:,i+1),F(i+1),0,K,D);
        opti.subject_to(G_Trapezoidal(x_mean(:,i),x_mean(:,i+1),x_mean_dot,x_mean_dot_next,dt) == 0);
        
        % Add constraints to determine A_k and C_k!
        % .....
        % .....
        
        
        % integration of covariance matrix
        % - positive definiteness preserving Lyapunov discretisation (ignore for now)
        DG_DW = DG_DW_Trapezoidal(C_k,dt);
        DG_DX = DG_DX_Trapezoidal(A_k,dt);
        DG_DZ = DG_DZ_Trapezoidal(A_k_next,dt);
        M_k = (DG_DZ)^(-1);
        % - integration step
        P_k = M_k*(DG_DX*P_k*DG_DX' + DG_DW*sigma_w*DG_DW')*M_k';
        opti.subject_to(Pvar(:,i+1) == P_k(:));
        
        % expected effort
        obj = obj + (F(:,i+1) + K*(x_mean(1,i+1)-pi) + D*x_mean(2,i+1))^2 + [K D]*P_k*[K;D]; % expected effort
    end
    
    % boundary conditions mean state
    opti.subject_to(x_mean(:,end) == [pi; 0]);
    opti.subject_to(x_mean(:,1) == 0);
    
    % cost function
    opti.minimize(obj + 1e4*(P_k(1,1)));  % Add minimization of variance, otherwise feedback gains will be zero and feedforward only is optimal solution
    
    optionssol.ipopt.linear_solver = 'mumps';
    optionssol.ipopt.tol = 1e-5; 
    optionssol.ipopt.dual_inf_tol = 1e-5;
    optionssol.ipopt.constr_viol_tol = 1e-7;
    
    opti.solver('ipopt',optionssol);
    
    sol = opti.solve();
    x_mean_sol = sol.value(x_mean);
    F_sol = sol.value(F);
    P_sol = sol.value(Pvar);
    K_sol = sol.value(K);
    D_sol = sol.value(D);
    F_sol = F_sol + K_sol*(x_mean_sol(1,:)-pi) + D_sol*x_mean_sol(2,:);
    
    subplot(3,2,1)
    plot(t,x_mean_sol(1,:),'Color',cs(noise_factor,:));
    title('position');
    hold on;
    subplot(3,2,2)
    plot(t,180/pi*sqrt(P_sol(1,:)),'Color',cs(noise_factor,:));
    title('position SD');
    hold on;
    set(gca, 'YScale', 'log')
    
    subplot(3,2,3)
    plot(t,x_mean_sol(2,:),'Color',cs(noise_factor,:));
    title('velocity');
    hold on;
    
    subplot(3,2,4)
    plot(t,180/pi*sqrt(P_sol(4,:)),'Color',cs(noise_factor,:));
    title('velocity SD');
    hold on;
    set(gca, 'YScale', 'log')
    
    subplot(3,2,5)
    plot(t,F_sol(1,:),'Color',cs(noise_factor,:));
    title('torque');
    hold on;
    
end