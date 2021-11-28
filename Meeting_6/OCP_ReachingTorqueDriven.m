function [result,succes] = OCP_ReachingTorqueDriven(target,forceField,wM_SD,wEE_SD,wEEdot_SD,guessName)
%% SECTION TITLE
% DESCRIPTIVE TEXT
if strcmp(target,'CIRCLE')
    targetNR = 1;
    %     guessName = 'result_CIRCLE.mat';
elseif strcmp(target,'BAR')
    targetNR = 2;
    %     guessName = 'result_BAR.mat';
elseif strcmp(target,'OBSTACLE')
    targetNR = 3;
    %     guessName = 'result_OBSTACLE.mat';
else
    error('Unknown target specified')
end

import casadi.*

% Set-up structure with data that specifies model
auxdata = initializeModelParameters();

% Additional simulation settings
T = auxdata.T;
dt = 0.01; auxdata.dt = dt; % time step
N = round(T/dt); auxdata.N = N;
time = 0:dt:T; auxdata.time = time;
nStates = 6; auxdata.nStates = nStates; % #states (6 muscle states + 2 positions + 2 velocities)
wM = (wM_SD*ones(2,1)).^2/dt; auxdata.wM = wM; % Motor noise: go from SD of continuous noise source to variance of discrete sample
wEE = (wEE_SD*ones(2,1)).^2/dt; auxdata.wEE = wEE; % Sensor position noise: go from SD of continuous noise source to variance of discrete sample
wEEdot = (wEEdot_SD*ones(2,1)).^2/dt; auxdata.wEEdot = wEEdot; % Sensor velocity noise: go from SD of continuous noise source to variance of discrete sample

sigma_w = [wM; wEE; wEEdot].*eye(6); auxdata.sigma_w = sigma_w; % Collect noise in covariance matrix

auxdata.forceField = forceField;

%%%% Define CasADi functions - for reasons of efficiency and to compute sensitivities (jacobians) of functions
functions = generateFunctions_OCP_ReachingTorqueDriven(auxdata);


opti = casadi.Opti(); % Create opti instance

% Initial and final position of the reaching movement
fsolve_options = optimoptions('fsolve','FiniteDifferenceType','central','StepTolerance',1e-10,'OptimalityTolerance',1e-10);
shoulder_pos_init = 20*pi/180;
shoulder_pos_final = 55*pi/180;
f = @(x)get_px(x,auxdata,shoulder_pos_init);
initial_pos = fsolve(f,ones,fsolve_options);
initial_pos = [shoulder_pos_init; initial_pos];

f = @(x)get_px(x,auxdata,shoulder_pos_final);
final_pos = fsolve(f,ones,fsolve_options);
final_pos = [shoulder_pos_final; final_pos];
EE_final = EndEffectorPos(final_pos,auxdata);

pos_guess = interp1( [0 auxdata.T] , [initial_pos final_pos]' , time );
X_init_guess = [0; 0; initial_pos; 0; 0];
e_ff_guess = -0.01*ones(2,N+1);
Pmat_init = [1e-6;1e-6;1e-4;1e-4;1e-7;1e-7;].*eye(6);
K_guess = zeros(8,N+1);
[X_guess, M_guess, EE_ref_guess, Pmat_guess] = approximateForwardSim(X_init_guess,Pmat_init,e_ff_guess,K_guess,auxdata,functions);

% Generate variables + initial guess
X_init = opti.variable(nStates,1); opti.set_initial(X_init, X_init_guess);
X = opti.variable(nStates,N+1); opti.set_initial(X, X_guess);
e_ff = opti.variable(2,N+1); opti.set_initial(e_ff, e_ff_guess);
K = opti.variable(2*4,N+1); opti.set_initial(K, K_guess);
EE_ref = opti.variable(4,N+1); opti.set_initial(EE_ref, EE_ref_guess);
M = opti.variable(nStates,nStates*N); opti.set_initial(M, M_guess);

% % Initial guess controls
% guessName = 'result_time_0.8_BAR_forceField_0_0.05_0.0003_0.0024_.mat';
% load(guessName);
% 
% % Interpolate guess to new mesh or new timing (if reaching movement timing
% % has changed)
% timeStepsGuess_new = size(result.X,2)-1;
% timeVec_new = 0:T/timeStepsGuess_new:T;
% e_ff_guess = interp1(timeVec_new,result.e_ff,time);
% K_guess = interp1(timeVec_new,result.K',time)';
% 
% X_init_guess = result.X(:,1);
% if abs(EE_ref_guess(1,end)) > 0.1 % control of initial guess is not very stable/precise, better to violate the dynamics and provide trajectories that satisfy the constraints
%     X_init = opti.variable(nStates,1); opti.set_initial(X_init, X_init_guess);
%     X = opti.variable(nStates,N+1); opti.set_initial(X, interp1(timeVec_new,result.X',time)');
%     e_ff = opti.variable(6,N+1); opti.set_initial(e_ff, e_ff_guess');
%     K = opti.variable(6*4,N+1); opti.set_initial(K, K_guess);
%     EE_ref = opti.variable(4,N+1); opti.set_initial(EE_ref, result.EE_ref);
%     M = opti.variable(nStates,nStates*N); opti.set_initial(M, result.M);
% else
%     X_init = opti.variable(nStates,1); opti.set_initial(X_init, X_init_guess);
%     X = opti.variable(nStates,N+1); opti.set_initial(X, X_guess);
%     e_ff = opti.variable(6,N+1); opti.set_initial(e_ff, e_ff_guess');
%     K = opti.variable(6*4,N+1); opti.set_initial(K, K_guess);
%     EE_ref = opti.variable(4,N+1); opti.set_initial(EE_ref, EE_ref_guess);
%     M = opti.variable(nStates,nStates*N); opti.set_initial(M, M_guess);
% end

Pmat_i = Pmat_init;
J_fb = 0;
for i = 1:N
    X_i = X(:,i);
    X_i_plus = X(:,i+1);
    e_ff_i = e_ff(:,i);
    e_ff_i_plus = e_ff(:,i+1);
    
    dX_i = functions.f_forwardMusculoskeletalDynamics(X_i,e_ff_i,0,0,0*wM,0*wEE,0*wEEdot);
    dX_i_plus = functions.f_forwardMusculoskeletalDynamics(X_i_plus,e_ff_i_plus,0,0,0*wM,0*wEE,0*wEEdot);
    opti.subject_to(functions.f_G_Trapezoidal(X_i,X_i_plus,dX_i,dX_i_plus)*1e3 == 0);
    
    M_i = M(:,(i-1)*nStates + 1:i*nStates);
    EE_ref_i = EE_ref(:,i);
    EE_ref_i_plus = EE_ref(:,i+1);
    
    K_i = reshape(K(:,i),2,4);
    K_i_plus = reshape(K(:,i+1),2,4);
    
    DdX_DX_i = functions.f_DdX_DX(X_i,e_ff_i,K_i,EE_ref_i,wM,wEE,wEEdot);
    DdZ_DX_i = functions.f_DdX_DX(X_i_plus,e_ff_i_plus,K_i_plus,EE_ref_i_plus,wM,wEE,wEEdot);
    DdX_Dw_i = functions.f_DdX_Dw(X_i,e_ff_i,K_i,EE_ref_i,wM,wEE,wEEdot);
    
    DG_DX_i = functions.f_DG_DX(DdX_DX_i);
    DG_DZ_i = functions.f_DG_DZ(DdZ_DX_i);
    DG_DW_i = functions.f_DG_DW(DdX_Dw_i);
    
    opti.subject_to(M_i*DG_DZ_i - eye(nStates) == 0);
    J_fb = J_fb + (functions.f_expectedEffort_fb(X_i,Pmat_i,K_i,EE_ref_i,wEE,wEEdot) + trace(Pmat_i(1:2,1:2)))/2; %expectedEffort_fb(i);
    
    % Obstacle
    if targetNR == 3
        if i*dt > T*5/8
            P_q_i = Pmat_i(3:4,3:4);
            P_EEPos_i = functions.f_P_EEPos(X(3:4,i),P_q_i);
            opti.subject_to(P_EEPos_i(1,1) < 0.004^2);
            opti.subject_to(-1e-4 < EE_ref_i(1) < 1e-4)
        end
    end
    
    Pmat_i = M_i*(DG_DX_i*Pmat_i*DG_DX_i' + DG_DW_i*sigma_w*DG_DW_i')*M_i'; % + dGdW*sigmaW*dGdW'
    
end
J_fb = J_fb + (functions.f_expectedEffort_fb(X_i_plus,Pmat_i,K_i_plus,EE_ref_i_plus,wEE,wEEdot) + trace(Pmat_i(1:2,1:2)))/2;

%% Boundary conditions
 
% % Initial conditions
EE = [EndEffectorPos(X(3:4,:),auxdata); EndEffectorVel(X(3:4,:),X(5:6,:),auxdata)];
opti.subject_to(EE - EE_ref == 0);
opti.subject_to(X(:,1) == X_init); % Set X_init to be equal to the first state
opti.subject_to(X_init(3:6) - [initial_pos;0;0] == 0); % Initial position and velocity
dX_init = functions.f_forwardMusculoskeletalDynamics(X_init,e_ff(:,1),0,0,0*wM,0*wEE,0*wEEdot); % Initial state derivative
opti.subject_to(dX_init([5:6]) == 0); % Initial acceleration equals zero (activations need to fullfill requirement)

% Reaching motion must end in the final reach position with zero angular joint velocity
opti.subject_to(functions.f_EEPos(X(3:4,end)) - EE_final == 0);
opti.subject_to(X(5:6,end) == [0;0]);

% Final acceleration equals zero (activations balanced)
dX_end = functions.f_forwardMusculoskeletalDynamics(X(:,end),e_ff(:,end),0,0,0*wM,0*wEE,0*wEEdot); % Initial state derivative
opti.subject_to(dX_end(5:6) == 0); % Initial acceleration equals zero (activations need to fullfill requirement)

   
% End effector endpoint accuracy
% Constrain the end point position and velocity standard deviation in the x and y
% directions to be below 0.4cm and 2cm/s respectively (depending on the
% target shape)
P_q_final = Pmat_i(3:4,3:4);
P_EEPos_final = functions.f_P_EEPos(X(3:4,end),P_q_final);
P_q_qdot_final = Pmat_i(3:6,3:6);
P_EEVel_final = functions.f_P_EEVel(X(3:4,end),X(5:6,end),P_q_qdot_final);
if targetNR == 1 || targetNR == 3
    opti.subject_to(P_EEPos_final(1,1) < 0.004^2);

else
%     opti.subject_to(P_EEPos_final(1,1) < 0.025^2);
end
opti.subject_to(P_EEPos_final(2,2) < 0.004^2);


% Limit variance on activations in endpoint
% opti.subject_to((Pmat_i(1,1) - 0.01^2) < 0); opti.subject_to((Pmat_i(2,2) - 0.01^2) < 0); opti.subject_to((Pmat_i(3,3) - 0.01^2) < 0);
% opti.subject_to((Pmat_i(4,4) - 0.01^2) < 0); opti.subject_to((Pmat_i(5,5) - 0.01^2) < 0); opti.subject_to((Pmat_i(6,6) - 0.01^2) < 0);

% Bounds on the feedforward excitations, activations and joint angles
opti.subject_to(-1 < e_ff(:) < 1);
opti.subject_to(-1 < X(1:2,:) < 1);
opti.subject_to(0 < X(3,:)*180/pi < 180);
opti.subject_to(0 < X(4,:)*180/pi < 180);

%% Cost function
opti.minimize(1e3*((sumsqr(e_ff)+sumsqr(X(1:2,:)))/2+J_fb)*dt);

%% Setup solver
% optionssol.ipopt.nlp_scaling_method = 'gradient-based';
optionssol.ipopt.linear_solver = 'ma57';
optionssol.ipopt.tol = 1e-3; % output.setup.nlp.ipoptoptions.tolerance;
optionssol.ipopt.dual_inf_tol = 3e-4;
optionssol.ipopt.constr_viol_tol = 1e-7;
optionssol.ipopt.max_iter = 10000;

optionssol.ipopt.hessian_approximation = 'limited-memory';
opti.solver('ipopt',optionssol);

% try
    tic;
%     result = solve_NLPSOL(opti,optionssol);
    sol = opti.solve();
    toc
    e_ff_sol = sol.value(e_ff);
    K_sol = sol.value(K);
    X_init_sol = sol.value(X_init);
    [X_sol, ~, EE_ref_sol, Pmat_sol] = approximateForwardSim(X_init_sol,Pmat_init,e_ff_sol,K_sol,auxdata,functions);
    
    for i = 1:N+1
        Pmat_sol_i = Pmat_sol(:,:,i);
        P_q_sol_i = Pmat_sol_i(3:4,3:4);
        P_EEPos_mat_i = functions.f_P_EEPos(X_sol(3:4,i),P_q_sol_i);
        P_EEPos_sol(:,i) = full([P_EEPos_mat_i(1,1); P_EEPos_mat_i(2,1); P_EEPos_mat_i(2,2)]);
        P_qdot_sol_i = Pmat_sol_i(3:6,3:6);
        P_EEVel_mat_i = functions.f_P_EEVel(X_sol(3:4,i),X_sol(5:6,i),P_qdot_sol_i);
        P_EEVel_sol(:,i) = full([P_EEVel_mat_i(1,1); P_EEVel_mat_i(2,1); P_EEVel_mat_i(2,2)]);
    end
    
    EEPos_sol = EndEffectorPos(X_sol(3:4,:),auxdata)';
    EEVel_sol = EndEffectorVel(X_sol(3:4,:),X_sol(5:6,:),auxdata)';
    
    clear result;
    result.e_ff = e_ff_sol';
    result.X = X_sol;
    result.a = X_sol(1:2,:)';
    result.q = X_sol(3:4,:)';
    result.qdot = X_sol(5:6,:)';
    result.K = K_sol;
    result.M = sol.value(M);
    result.Pmat = Pmat_sol;
    result.time = 0:dt:T;
    result.auxdata = auxdata;
    result.EEPos = EEPos_sol;
    result.EEVel = EEVel_sol;
    result.P_EEPos = P_EEPos_sol;
    result.P_EEVel = P_EEVel_sol;
    result.EE_ref = EE_ref_sol;
    
    % Cost function
    %-feedback
    J_fb = NaN(auxdata.N+1,1);
    J_fb_corrective = NaN(auxdata.N+1,1);
    J_fb_sensoryNoise = NaN(auxdata.N+1,1);
    J_var_act = NaN(auxdata.N+1,1);
    e_fb_corrective = NaN(2,auxdata.N+1); 
    
    for i = 1:auxdata.N+1
        J_fb(i) = full(functions.f_expectedEffort_fb(X_sol(:,i),Pmat_sol(:,:,i),reshape(K_sol(:,i),2,4),EE_ref_sol(:,i),wEE,wEEdot));
        J_fb_corrective(i) = full(functions.f_expectedEffort_fb_corrective(X_sol(:,i),Pmat_sol(:,:,i),reshape(K_sol(:,i),2,4),EE_ref_sol(:,i)));
        J_fb_sensoryNoise(i) = full(functions.f_expectedEffort_fb_sensoryNoise(reshape(K_sol(:,i),2,4),wEE,wEEdot));
        J_var_act(i) = trace(Pmat_sol(1:2,1:2,i));
        e_fb_corrective(:,i) = full(functions.f_fb_corrective(X_sol(:,i),Pmat_sol(:,:,i),reshape(K_sol(:,i),2,4),EE_ref_sol(:,i)));
    end
    J_fb_total = (sum(J_fb)+sum(J_var_act))/2;
    J_ff_total = (sumsqr(e_ff_sol)+sumsqr(X_sol(1:2,:)))/2;
    result.J_fb = J_fb;
    result.J_var_act = J_var_act;

    result.J_fb_corrective = J_fb_corrective;
    result.J_fb_sensoryNoise = J_fb_sensoryNoise;
    result.J_fb_total = J_fb_total;
    result.J_ff_total = J_ff_total;
    
    succes = true;
    
% catch
%     succes = false;
end



