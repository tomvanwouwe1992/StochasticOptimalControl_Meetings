function functions = generateFunctions_OCP_ReachingTorqueDriven(auxdata)
import casadi.*

% FORWARD DYNAMICS
e_ff_MX = MX.sym('u_MX',2);       % Controls (torque excitations - feedforward)
X_MX = MX.sym('X_MX',auxdata.nStates);    % States (torque activations (2) and joint kinematics)
wM_MX = MX.sym('wM_MX',2);        % Motor noise
wEE_MX = MX.sym('wEE_MX',2);      % Sensor position noise
wEEdot_MX = MX.sym('wEEdot_MX',2);% Sensor velocity noise
EE_ref_MX = MX.sym('EE_ref_MX', 4); % Reference end-effector kinematics
EE_MX = [EndEffectorPos(X_MX(3:4),auxdata); EndEffectorVel(X_MX(3:4),X_MX(5:6),auxdata)]; % End-effector kinematics computed from joint kinematics 
K_MX = MX.sym('K_MX',2,4);        % Feedback gains
e_fb_MX = K_MX*((EE_MX - EE_ref_MX) + [wEE_MX;wEEdot_MX]); % Feedback excitations composed from the feedback of the end-effector kinematics error, corrupted by sensor noise
u_MX = e_ff_MX + e_fb_MX;         % Total muscle excitations (ff + fb)

% unperturbed stochastic forward dynamics

F_forceField = auxdata.forceField*(auxdata.l1*cos(X_MX(3,:)) + auxdata.l2*cos(X_MX(3,:)+X_MX(4,:)));
T_forceField = -F_forceField*[auxdata.l2*sin(X_MX(3,:)+X_MX(4,:))+auxdata.l1*sin(X_MX(3,:));auxdata.l2*sin(X_MX(3,:)+X_MX(4,:))];

dX_MX=   [  (u_MX - X_MX(1:2))/auxdata.tau; % activation dyn
             X_MX(5:6);                     % velocity
             armForwardDynamics([auxdata.T_shoulder_max;auxdata.T_elbow_max].*X_MX(1:2),X_MX(4),X_MX(5:6),wM_MX,T_forceField,auxdata)]; % acceleration

f_forwardMusculoskeletalDynamics = Function('f_forwardMusculoskeletalDynamics',{X_MX,e_ff_MX,K_MX,EE_ref_MX,wM_MX,wEE_MX,wEEdot_MX},{dX_MX});  
functions.f_forwardMusculoskeletalDynamics = f_forwardMusculoskeletalDynamics;

% Sensitivity of forward dynamics to states
DdX_DX_MX = jacobian(dX_MX, X_MX');
f_DdX_DX = Function('f_DdX_DX',{X_MX,e_ff_MX,K_MX,EE_ref_MX,wM_MX,wEE_MX,wEEdot_MX},{DdX_DX_MX});
functions.f_DdX_DX = f_DdX_DX;

% Sensitivity of forward dynamics to motor noise
DdX_DwM_MX = jacobian(dX_MX, [wM_MX' wEE_MX' wEEdot_MX']);
f_DdX_Dw = Function('f_DdX_Dw',{X_MX,e_ff_MX,K_MX,EE_ref_MX,wM_MX,wEE_MX,wEEdot_MX},{DdX_DwM_MX});
functions.f_DdX_Dw = f_DdX_Dw;

% Trapezoidal integration scheme (implicit)
X_plus_MX = MX.sym('X_plus_MX',auxdata.nStates); % States at mesh end
dX_MX = MX.sym('dX_MX',auxdata.nStates); % State derivative
dX_plus_MX = MX.sym('dX_plus_MX',auxdata.nStates); % State derivative at mesh end
f_G_Trapezoidal = Function('f_G_Trapezoidal',{X_MX,X_plus_MX,dX_MX,dX_plus_MX},{G_Trapezoidal(X_MX,X_plus_MX,dX_MX,dX_plus_MX,auxdata.dt)});
functions.f_G_Trapezoidal = f_G_Trapezoidal;

% Sensitivity of trapezoidal integration scheme to changes in initial state
DdX_DX_MX = MX.sym('DdX_DX_MX',auxdata.nStates,auxdata.nStates); 
DG_DX_MX = DG_DX_Trapezoidal(DdX_DX_MX,auxdata.dt);
f_DG_DX = Function('f_DG_DX',{DdX_DX_MX},{DG_DX_MX});
functions.f_DG_DX = f_DG_DX;

% Sensitivity of trapezoidal integration scheme to changes in final state
DG_DZ_MX = DG_DZ_Trapzoidal(DdX_DX_MX,auxdata.dt);
f_DG_DZ = Function('f_DG_DZ',{DdX_DX_MX},{DG_DZ_MX});
functions.f_DG_DZ = f_DG_DZ;

% Sensitivity of trapezoidal integration scheme to changes in motor noise
DdX_Dw_MX = MX.sym('DdX_Dw_MX',auxdata.nStates,6); 
DG_DW_MX = DG_DW_Trapezoidal(DdX_Dw_MX,auxdata.dt);
f_DG_DW = Function('f_DG_DW',{DdX_Dw_MX},{DG_DW_MX});
functions.f_DG_DW = f_DG_DW;

% Some other useful function definitions
% End effector position and variability
q_MX = MX.sym('q_MX',2); P_q_MX = MX.sym('P_q_MX',2,2);
EEPos_MX = EndEffectorPos(q_MX,auxdata); f_EEPos = Function('f_EEPos',{q_MX},{EEPos_MX}); % End effector position
functions.f_EEPos = f_EEPos;
P_EEPos_MX = jacobian(EEPos_MX,q_MX)*P_q_MX*jacobian(EEPos_MX,q_MX)'; f_P_EEPos = Function('f_P_EEPos',{q_MX,P_q_MX},{P_EEPos_MX}); % Covariance of end effector position
functions.f_P_EEPos = f_P_EEPos;
% End effector velocity and variability
q_MX = MX.sym('q_MX',2); qdot_MX = MX.sym('q_MX',2); P_qdot_MX = MX.sym('P_qdot_MX',4,4);
EEVel_MX = EndEffectorVel(q_MX,qdot_MX,auxdata); f_EEVel = Function('f_EEVel',{q_MX,qdot_MX},{EEVel_MX}); % End effector position
functions.f_EEVel = f_EEVel;
P_EEVel_MX = jacobian(EEVel_MX,[q_MX qdot_MX])*P_qdot_MX*jacobian(EEVel_MX,[q_MX qdot_MX])'; f_P_EEVel = Function('f_P_EEVel',{q_MX,qdot_MX,P_qdot_MX},{P_EEVel_MX}); % Covariance of end effector position
functions.f_P_EEVel = f_P_EEVel;

%
P_MX = MX.sym('P_MX',auxdata.nStates,auxdata.nStates);
sensoryNoise_MX = [wEE_MX;wEEdot_MX].*eye(4);
expectedEffort_fb_MX = trace(jacobian(e_fb_MX,X_MX)*P_MX*jacobian(e_fb_MX,X_MX)') + trace(K_MX*sensoryNoise_MX*K_MX');
f_expectedEffort_fb = Function('f_expectedEffort_fb',{X_MX,P_MX,K_MX,EE_ref_MX,wEE_MX,wEEdot_MX},{expectedEffort_fb_MX});
functions.f_expectedEffort_fb = f_expectedEffort_fb;

f_fb_corrective = Function('f_fb_corrective',{X_MX,P_MX,K_MX,EE_ref_MX},{diag(jacobian(e_fb_MX,X_MX)*P_MX*jacobian(e_fb_MX,X_MX)')});
functions.f_fb_corrective = f_fb_corrective;

f_expectedEffort_fb_corrective = Function('f_expectedEffort_fb_corrective',{X_MX,P_MX,K_MX,EE_ref_MX},{trace(jacobian(e_fb_MX,X_MX)*P_MX*jacobian(e_fb_MX,X_MX)')});
functions.f_expectedEffort_fb_corrective = f_expectedEffort_fb_corrective;

f_expectedEffort_fb_sensoryNoise = Function('f_expectedEffort_fb_sensoryNoise',{K_MX,wEE_MX,wEEdot_MX},{trace(K_MX*sensoryNoise_MX*K_MX')});
functions.f_expectedEffort_fb_sensoryNoise = f_expectedEffort_fb_sensoryNoise;
end