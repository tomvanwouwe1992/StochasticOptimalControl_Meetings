clear all; close all;
import casadi.*
load('resultAdapted_dt0.05_fast1.mat')
N_robust = length(result.q1_sol)-1;
dt_robust = 0.03125;
time = 0:dt_robust:0.5;
indicesP = NaN(10,10);
vecindicesP = NaN(55,1);

q1_sol = result.q1_sol; q2_sol = result.q2_sol; q3_sol = result.q3_sol; q4_sol = result.q4_sol; q5_sol = result.q5_sol;
dq1_sol = result.dq1_sol; dq2_sol = result.dq2_sol; dq3_sol = result.dq3_sol; dq4_sol = result.dq4_sol; dq5_sol = result.dq5_sol;
T1_sol = result.T1_sol; T2_sol = result.T2_sol; T3_sol = result.T3_sol; T4_sol = result.T4_sol; T5_sol = result.T5_sol;
refTraj_sol = result.refTraj_sol;
K_sol = result.K_sol;
ct = 1;
for i = 1:10
    for j = i:10
        vecindicesP(ct,1) = (i-1)*10+j;
        indicesP(i,j) = ct;
        ct = ct + 1;
    end
end
for i = 1:9
    for j = i+1:10
        indicesP(j,i) = indicesP(i,j);
    end
end

varianceIndices = diag(indicesP);
covarianceIndices = 1:55;
covarianceIndices(varianceIndices) = [];

K_sol = result.K_sol;
feedbackTransform = [1  0 0 0 0 0 0 0 0 0;...
                     1 -1 0 0 0 0 0 0 0 0;...
                     0 1 -1 0 0 0 0 0 0 0;...
                     0 0 -1 1 0 0 0 0 0 0;...
                     0 0 0 -1 1 0 0 0 0 0;...
                     0 0 0 0 0 1  0 0 0 0;...
                     0 0 0 0 0 1 -1 0 0 0;...
                     0 0 0 0 0 0 1 -1 0 0;...
                     0 0 0 0 0 0 0 -1 1 0;...
                     0 0 0 0 0 0 0 0 -1 1];
                 
   
X_mean_MX = MX.sym('X_mean_MX',10,1);
X_conf1_MX = MX.sym('X_conf1_MX',1,1);
X_conf2_MX = MX.sym('X_conf2_MX',1,1);
X_conf3_MX = MX.sym('X_conf3_MX',1,1);
X_conf4_MX = MX.sym('X_conf4_MX',1,1);
X_conf5_MX = MX.sym('X_conf5_MX',1,1);
X_conf6_MX = MX.sym('X_conf6_MX',1,1);
X_conf7_MX = MX.sym('X_conf7_MX',1,1);
X_conf8_MX = MX.sym('X_conf8_MX',1,1);
X_conf9_MX = MX.sym('X_conf9_MX',1,1);
X_conf10_MX = MX.sym('X_conf10_MX',1,1);
X_conf_MX = [X_conf1_MX;X_conf2_MX;X_conf3_MX;X_conf4_MX;X_conf5_MX;X_conf6_MX;X_conf7_MX;X_conf8_MX;X_conf9_MX;X_conf10_MX]';




% 
%    P_sol = result.P_sol;              
% for k = 1:N_robust+1
%     Pvec_sol = P_sol(:,k);
%     Pmat_sol = Pvec_sol(indicesP);
%     Kk = K_sol(:,(k-1)*10+1:10*k);
%     cost_expectedEffort = cost_expectedEffort + sum(diag((Kk*feedbackTransform)*Pmat_sol*(Kk*feedbackTransform)'))*0.1;
%     T_std(:,k) = sqrt(diag((Kk*feedbackTransform)*Pmat_sol*(Kk*feedbackTransform)'));
%     q_std(:,k) = sqrt(diag(Pmat_sol(1:5,1:5)));
%     dq_std(:,k) = sqrt(diag(Pmat_sol(6:10,6:10)));
%     [~,flag] = chol(Pmat_sol);
%     flag
% end
L_sol = result.L_sol;
for k = 1:N_robust+1
    Lvec_sol = L_sol(:,k);
    Lmat_sol = tril(ones(10)).*Lvec_sol(indicesP);
    Pmat_sol = Lmat_sol*Lmat_sol';
    Kk = K_sol(:,(k-1)*10+1:10*k);
    cost_expectedEffort(k) = sum(diag((Kk*feedbackTransform)*Pmat_sol*(Kk*feedbackTransform)'))*0.1;
    T_std(:,k) = sqrt(diag(Kk*feedbackTransform*Pmat_sol*(Kk*feedbackTransform)'));
    q_std(:,k) = sqrt(diag(Pmat_sol(1:5,1:5)));
        dq_std(:,k) = sqrt(diag(Pmat_sol(6:10,6:10)));
        Pmat_rel = feedbackTransform*Pmat_sol*feedbackTransform';
    qrel_std(:,k) = sqrt(diag(Pmat_rel(1:5,1:5)));
    dqrel_std(:,k) = sqrt(diag(Pmat_rel(6:10,6:10)));

    [~,flag] = chol(Pmat_sol);
    flag
end
cost_Torque = 0.1*(sumsqr(T1_sol) + sumsqr(T2_sol) + sumsqr(T3_sol) + sumsqr(T4_sol) + sumsqr(T5_sol)) ;
figure(1)
subplot(3,2,1)
plot(time,T1_sol,'LineWidth',2,'Color','b');hold on;
plot(time,T1_sol+3*T_std(1,:),'--','LineWidth',1,'Color','b');hold on;
plot(time,T1_sol-3*T_std(1,:),'--','LineWidth',1,'Color','b');hold on;

subplot(3,2,2)
plot(time,T2_sol,'LineWidth',2);hold on;
plot(time,T2_sol+3*T_std(2,:),'--','LineWidth',1,'Color','b');hold on;
plot(time,T2_sol-3*T_std(2,:),'--','LineWidth',1,'Color','b');hold on;

subplot(3,2,3)
plot(time,T3_sol,'LineWidth',2);hold on;
plot(time,T3_sol+3*T_std(3,:),'--','LineWidth',1,'Color','b');hold on;
plot(time,T3_sol-3*T_std(3,:),'--','LineWidth',1,'Color','b');hold on;

subplot(3,2,4)
plot(time,T4_sol,'LineWidth',2);hold on;
plot(time,T4_sol+3*T_std(4,:),'--','LineWidth',1,'Color','b');hold on;
plot(time,T4_sol-3*T_std(4,:),'--','LineWidth',1,'Color','b');hold on;

subplot(3,2,5)
plot(time,T5_sol,'LineWidth',2); hold on;
plot(time,T5_sol+3*T_std(5,:),'--','LineWidth',1,'Color','b');hold on;
plot(time,T5_sol-3*T_std(5,:),'--','LineWidth',1,'Color','b');hold on;


m1 = 5; m5 = 5;
m2 = 9.3; m4 = 9.3;
m3 = 34;
I1 = 0.1; I5 = 0.1;
I2 = 0.14; I4 = 0.14;
I3 = 1.43
l1 = 0.45; l5 = 0.45;
l2 = 0.43; l4 = 0.43;
l3 = 0.625;
lc1 = 0.2; lc5 = 0.25;
lc2 = 0.22; lc4 = 0.21;
lc3 = 0.32;
g = 9.81;
nStates = 10;
% Forward integration with solution

q1_MX = MX.sym('q1_MX',1); q2_MX = MX.sym('q2_MX',1); q3_MX = MX.sym('q3_MX',1); q4_MX = MX.sym('q4_MX',1); q5_MX = MX.sym('q5_MX',1);
dq1_MX = MX.sym('dq1_MX',1); dq2_MX = MX.sym('dq2_MX',1); dq3_MX = MX.sym('dq3_MX',1); dq4_MX = MX.sym('dq4_MX',1); dq5_MX = MX.sym('dq5_MX',1);
q1_min_MX = MX.sym('q1_min_MX',1); q2_min_MX = MX.sym('q2_min_MX',1); q3_min_MX = MX.sym('q3_min_MX',1); q4_min_MX = MX.sym('q4_min_MX',1); q5_min_MX = MX.sym('q5_min_MX',1);
q1_plus_MX = MX.sym('q1_plus_MX',1); q2_plus_MX = MX.sym('q2_plus_MX',1); q3_plus_MX = MX.sym('q3_plus_MX',1); q4_plus_MX = MX.sym('q4_plus_MX',1); q5_plus_MX = MX.sym('q5_plus_MX',1);

dq1_min_MX = MX.sym('dq1_min_MX',1); dq2_min_MX = MX.sym('dq2_min_MX',1); dq3_min_MX = MX.sym('dq3_min_MX',1); dq4_min_MX = MX.sym('dq4_min_MX',1); dq5_min_MX = MX.sym('dq5_min_MX',1);
dq1_plus_MX = MX.sym('dq1_plus_MX',1); dq2_plus_MX = MX.sym('dq2_plus_MX',1); dq3_plus_MX = MX.sym('dq3_plus_MX',1); dq4_plus_MX = MX.sym('dq4_plus_MX',1); dq5_plus_MX = MX.sym('dq5_plus_MX',1);
state_plus_sol = SX.sym('state_plus_sol',10);


ddq1_MX = MX.sym('ddq1_MX',1); ddq2_MX = MX.sym('ddq2_MX',1); ddq3_MX = MX.sym('ddq3_MX',1); ddq4_MX = MX.sym('ddq4_MX',1); ddq5_MX = MX.sym('ddq5_MX',1);
T1_MX = MX.sym('T1_MX',1); T2_MX = MX.sym('T2_MX',1); T3_MX = MX.sym('T3_MX',1); T4_MX = MX.sym('T4_MX',1); T5_MX = MX.sym('T5_MX',1);
refTraj_MX = MX.sym('refTraj_MX',10,1);
K_MX = MX.sym('K_MX',5,10);
relativeJointPos_MX = relativeJointPos(q1_MX,q2_MX,q3_MX,q4_MX,q5_MX);
relativeJointVel_MX = relativeJointVel(dq1_MX,dq2_MX,dq3_MX,dq4_MX,dq5_MX);
kinFB_MX = [relativeJointPos_MX; relativeJointVel_MX] - refTraj_MX;
T_FB_MX = K_MX*kinFB_MX;
T1_FF_FB_MX = T_FB_MX(1) + T1_MX;
T2_FF_FB_MX = T_FB_MX(2) + T2_MX;
T3_FF_FB_MX = T_FB_MX(3) + T3_MX;
T4_FF_FB_MX = T_FB_MX(4) + T4_MX;
T5_FF_FB_MX = T_FB_MX(5) + T5_MX;

T_full_out = [T1_FF_FB_MX; T2_FF_FB_MX; T3_FF_FB_MX; T4_FF_FB_MX; T5_FF_FB_MX];

f_T_Full = Function('f_T_Full',{T1_MX,T2_MX,T3_MX,T4_MX,T5_MX,refTraj_MX,q1_MX,q2_MX,q3_MX,q4_MX,q5_MX,dq1_MX,dq2_MX,dq3_MX,dq4_MX,dq5_MX,K_MX},{T_full_out});

f_feedbackSignal_Full = Function('f_feedbackSignal_Full',{refTraj_MX,q1_MX,q2_MX,q3_MX,q4_MX,q5_MX,dq1_MX,dq2_MX,dq3_MX,dq4_MX,dq5_MX},{kinFB_MX});

eq_SysDyn_Error = eq_SysDyn(I1,I2,I3,I4,I5,T1_FF_FB_MX,T2_FF_FB_MX,T3_FF_FB_MX,T4_FF_FB_MX,T5_FF_FB_MX,ddq1_MX,ddq2_MX,ddq3_MX,ddq4_MX,ddq5_MX,dq1_MX,dq2_MX,dq3_MX,dq4_MX,dq5_MX,g,l1,l2,l4,lc1,lc2,lc3,lc4,lc5,m1,m2,m3,m4,m5,q1_MX,q2_MX,q3_MX,q4_MX,q5_MX)
f_eq_SysDyn_Error = Function('f_eq_SysDyn_Error',{T1_MX,T2_MX,T3_MX,T4_MX,T5_MX,ddq1_MX,ddq2_MX,ddq3_MX,ddq4_MX,ddq5_MX,dq1_MX,dq2_MX,dq3_MX,dq4_MX,dq5_MX,q1_MX,q2_MX,q3_MX,q4_MX,q5_MX,refTraj_MX,K_MX},{eq_SysDyn_Error});

heelStrike_error = eq_HeelStrike(I1,I2,I3,I4,I5,dq1_min_MX,dq2_min_MX,dq3_min_MX,dq4_min_MX,dq5_min_MX,dq1_plus_MX,dq2_plus_MX,dq3_plus_MX,dq4_plus_MX,dq5_plus_MX,l1,l2,l4,l5,lc1,lc2,lc3,lc4,lc5,m1,m2,m3,m4,m5,q1_min_MX,q2_min_MX,q3_min_MX,q4_min_MX,q5_min_MX,q1_plus_MX,q2_plus_MX,q3_plus_MX,q4_plus_MX,q5_plus_MX)
f_heelStrike_error = Function('f_heelStrike_error',{dq1_min_MX,dq2_min_MX,dq3_min_MX,dq4_min_MX,dq5_min_MX,dq1_plus_MX,dq2_plus_MX,dq3_plus_MX,dq4_plus_MX,dq5_plus_MX,q1_min_MX,q2_min_MX,q3_min_MX,q4_min_MX,q5_min_MX,q1_plus_MX,q2_plus_MX,q3_plus_MX,q4_plus_MX,q5_plus_MX},{heelStrike_error});




figure(1)
K_sol_11 = K_sol(2,1:10:end);
K_sol_12 = K_sol(2,2:10:end);
K_sol_13 = K_sol(2,3:10:end);
K_sol_14 = K_sol(2,4:10:end);
K_sol_15 = K_sol(2,5:10:end);
K_sol_16 = K_sol(2,6:10:end);
K_sol_17 = K_sol(2,7:10:end);
K_sol_18 = K_sol(2,8:10:end);
K_sol_19 = K_sol(2,9:10:end);
K_sol_110 = K_sol(2,10:10:end);
colors = linspecer(5);
K_sol_ank = [K_sol_11' K_sol_12' K_sol_13' K_sol_14' K_sol_15' K_sol_16' K_sol_17' K_sol_18' K_sol_19' K_sol_110']


q_sol = [q1_sol;q2_sol;q3_sol;q4_sol;q5_sol];

dq_sol = [dq1_sol;dq2_sol;dq3_sol;dq4_sol;dq5_sol];
for i = 1:5
% Reference simulation

    Lvec_sol = L_sol(:,1);
    Lmat_sol = tril(ones(10)).*Lvec_sol(indicesP);
    Pmat_sol = Lmat_sol*Lmat_sol';
    Pvec_sol = Pmat_sol(vecindicesP);
    initPert = mvnrnd(zeros(10,1),Pmat_sol)
q_fw = NaN(5,N_robust+1); q_fw(:,1) = [q1_sol(1);q2_sol(1);q3_sol(1);q4_sol(1);q5_sol(1)] + 0*initPert(1:5)'; % Deviation on initial positions of 0.6° std 
dq_fw = NaN(5,N_robust+1); dq_fw(:,1) = [dq1_sol(1);dq2_sol(1);dq3_sol(1);dq4_sol(1);dq5_sol(1)] + 0*initPert(6:10)'; % Deviation on initial positions of 0.6°/s std
T_noise = 0*0.2/0.03125*randn(5,N_robust+1); % No motor noise (for now)
T_fw = NaN(5,N_robust+1);
feedbackSignal_fw = NaN(10,N_robust+1);
for k = 1:N_robust
    COM_segments = COM_segments_fcn(l1,l2,l4,lc1,lc2,lc3,lc4,lc5,q_fw(1,k),q_fw(2,k),q_fw(3,k),q_fw(4,k),q_fw(5,k));
COM_segments_y = COM_segments(2:3:end);
if k == 1
    T_pert = -2*[m1 m2 m3 m4 m5;...
          0  m2 m3 m4 m5;...
          0  0  m3 0  0;...
          0  0  0  m4 m5;...
          0  0  0  0  m5]*COM_segments_y;
      T_noise(:,k) = T_noise(:,k) + T_pert;
      T_noise(:,k+1) = T_noise(:,k+1) + T_pert;
end
if k == 3
    T_pert = -2*[m1 m2 m3 m4 m5;...
          0  m2 m3 m4 m5;...
          0  0  m3 0  0;...
          0  0  0  m4 m5;...
          0  0  0  0  m5]*COM_segments_y;
      T_noise(:,k) = T_noise(:,k) + T_pert;
      T_noise(:,k+1) = T_noise(:,k+1) + T_pert;
end
    Urf = SX.sym('Urf',30);
    dX =  Urf(1:10);
    dX_plus = Urf(11:20);
    X_plus = Urf(21:30);
    Kk = K_sol(:,(k-1)*10+1:10*k);
    Kk_plus = K_sol(:,(k)*10+1:10*(k+1));
    rf = rootfinder('rf','newton',struct('x',Urf,'g',[f_eq_SysDyn_Error(T1_sol(k) + T_noise(1,k),T2_sol(k) + T_noise(2,k),T3_sol(k) + T_noise(3,k),T4_sol(k) + T_noise(4,k),T5_sol(k) + T_noise(5,k), ...
                                         dX(6),dX(7),dX(8),dX(9),dX(10),...
                                         dX(1),dX(2),dX(3),dX(4),dX(5),...
                                         q_fw(1,k),q_fw(2,k),q_fw(3,k),q_fw(4,k),q_fw(5,k),refTraj_sol(:,k),Kk); ...
                                         f_eq_SysDyn_Error(T1_sol(k+1) + T_noise(1,k+1),T2_sol(k+1) + T_noise(2,k+1),T3_sol(k+1) + T_noise(3,k+1),T4_sol(k+1) + T_noise(4,k+1),T5_sol(k+1) + T_noise(5,k+1), ...
                                         dX_plus(6),dX_plus(7),dX_plus(8),dX_plus(9),dX_plus(10),...
                                         dX_plus(1),dX_plus(2),dX_plus(3),dX_plus(4),dX_plus(5),...
                                         X_plus(1),X_plus(2),X_plus(3),X_plus(4),X_plus(5),refTraj_sol(:,k+1),Kk_plus); ...
                                         dX(1:5) - [dq_fw(1,k);dq_fw(2,k);dq_fw(3,k);dq_fw(4,k);dq_fw(5,k)];...
                                         X_plus - ([q_fw(1,k);q_fw(2,k);q_fw(3,k);q_fw(4,k);q_fw(5,k);dq_fw(1,k);dq_fw(2,k);dq_fw(3,k);dq_fw(4,k);dq_fw(5,k)] + 0.5*dt_robust*(dX_plus + dX));...
                                         dX_plus(1:5) - X_plus(6:10);...
                                         ]),struct('abstol',1e-12));
    solution_u = rf(zeros(30,1),[]);
    solution_u = full(solution_u);
    dX_sol = solution_u(1:10);
    dX_plus_sol = solution_u(11:20);
    X_plus_sol = solution_u(21:30);
    q_fw(:,k+1) = X_plus_sol(1:5);
    dq_fw(:,k+1) = X_plus_sol(6:10);
    feedbackSignal_fw(:,k) = full(f_feedbackSignal_Full(refTraj_sol(:,k),q_fw(1,k),q_fw(2,k),q_fw(3,k),q_fw(4,k),q_fw(5,k),dq_fw(1,k),dq_fw(2,k),dq_fw(3,k),dq_fw(4,k),dq_fw(5,k)));
    T_fw(:,k) = full(f_T_Full(T1_sol(k) ,T2_sol(k) ,T3_sol(k) ,T4_sol(k),T5_sol(k),...
                         refTraj_sol(:,k),q_fw(1,k),q_fw(2,k),q_fw(3,k),q_fw(4,k),q_fw(5,k),dq_fw(1,k),dq_fw(2,k),dq_fw(3,k),dq_fw(4,k),dq_fw(5,k),Kk));
P_J = JointPos(l1,l2,l3,l4,l5,q_fw(1,k+1),q_fw(2,k+1),q_fw(3,k+1),q_fw(4,k+1),q_fw(5,k+1));
P_J(10)
end
k = k+1;
T_fw(:,k) = full(f_T_Full(T1_sol(k) ,T2_sol(k) ,T3_sol(k) ,T4_sol(k) ,T5_sol(k),...
                         refTraj_sol(:,k),q_fw(1,k),q_fw(2,k),q_fw(3,k),q_fw(4,k),q_fw(5,k),dq_fw(1,k),dq_fw(2,k),dq_fw(3,k),dq_fw(4,k),dq_fw(5,k),Kk_plus));
feedbackSignal_fw(:,k) = full(f_feedbackSignal_Full(refTraj_sol(:,k),q_fw(1,k),q_fw(2,k),q_fw(3,k),q_fw(4,k),q_fw(5,k),dq_fw(1,k),dq_fw(2,k),dq_fw(3,k),dq_fw(4,k),dq_fw(5,k)));
refSim.q_fw = q_fw;
refSim.dq_fw = dq_fw;
T_sol = [T1_sol; T2_sol; T3_sol; T4_sol; T5_sol];
P_J = JointPos(l1,l2,l3,l4,l5,q_fw(1,end),q_fw(2,end),q_fw(3,end),q_fw(4,end),q_fw(5,end));
P_J(10)
figure(1)

for i = 1:5
subplot(3,2,i)
plot(time,T_fw(i,:)'); hold on;
% plot(time,T_sol(i,:)','--'); hold on;
end
figure(2)
for i = 1:5
subplot(3,2,i)
plot(time,q_fw(i,:)'); hold on;
% plot(time,q_sol(i,:)','--'); hold on;
end

figure(6)
qfw_rel = relativeJointPos(q_fw(1,:),q_fw(2,:),q_fw(3,:),q_fw(4,:),q_fw(5,:));

for i = 1:5
subplot(3,2,i)
plot(time,qfw_rel(i,:)'); hold on;
% plot(time,q_sol(i,:)','--'); hold on;
end

figure(7)
dqfw_rel = relativeJointPos(dq_fw(1,:),dq_fw(2,:),dq_fw(3,:),dq_fw(4,:),dq_fw(5,:));

for i = 1:5
subplot(3,2,i)
plot(time,dqfw_rel(i,:)'); hold on;
% plot(time,q_sol(i,:)','--'); hold on;
end

figure(3)
for i = 1:5
subplot(3,2,i)
plot(time,dq_fw(i,:)'); hold on;
% plot(time,dq_sol(i,:)','--'); hold on;
end

end

figure(2)
for i = 1:5
subplot(3,2,i)
plot(time,q_sol(i,:)','LineWidth',2); hold on;
plot(time,q_sol(i,:)' + 3*q_std(i,:)','--'); hold on;
plot(time,q_sol(i,:)' - 3*q_std(i,:)','--'); hold on;
end

figure(6)
qsol_rel = relativeJointPos(q_sol(1,:),q_sol(2,:),q_sol(3,:),q_sol(4,:),q_sol(5,:));
for i = 1:5
subplot(3,2,i)
plot(time,qsol_rel(i,:)'); hold on;
plot(time,qsol_rel(i,:)' + 3*qrel_std(i,:)','--'); hold on;
plot(time,qsol_rel(i,:)' - 3*qrel_std(i,:)','--'); hold on;
% plot(time,q_sol(i,:)','--'); hold on;
end

figure(7)
dqsol_rel = relativeJointPos(dq_sol(1,:),dq_sol(2,:),dq_sol(3,:),dq_sol(4,:),dq_sol(5,:));
for i = 1:5
subplot(3,2,i)
plot(time,dqsol_rel(i,:)'); hold on;
plot(time,dqsol_rel(i,:)' + 3*dqrel_std(i,:)','--'); hold on;
plot(time,dqsol_rel(i,:)' - 3*dqrel_std(i,:)','--'); hold on;
% plot(time,q_sol(i,:)','--'); hold on;
end

figure(3)
for i = 1:5
subplot(3,2,i)
plot(time,dq_sol(i,:)','LineWidth',2); hold on;
plot(time,dq_sol(i,:)' + 3*dq_std(i,:)','--'); hold on;
plot(time,dq_sol(i,:)' - 3*dq_std(i,:)','--'); hold on;
end

figure(4)
subplot(1,2,1)
plot(time,K_sol(4,1:10:end)); hold on;
plot(time,K_sol(4,2:10:end));hold on;
plot(time,K_sol(4,3:10:end));hold on;
plot(time,K_sol(4,4:10:end));hold on;
plot(time,K_sol(4,5:10:end));hold on;

subplot(1,2,2)
plot(time,K_sol(4,6:10:end));hold on;
plot(time,K_sol(4,7:10:end));hold on;
plot(time,K_sol(4,8:10:end));hold on;
plot(time,K_sol(4,9:10:end));hold on;
plot(time,K_sol(4,10:10:end));

figure(5)
subplot(1,2,1)
plot(time,feedbackSignal_fw(1:5,:))'; hold on;
subplot(1,2,2)
plot(time,feedbackSignal_fw(6:10,:))';

figure(20)
COM = COM_fcn(l1,l2,l4,lc1,lc2,lc3,lc4,lc5,m1,m2,m3,m4,m5,q_sol(1,:),q_sol(2,:),q_sol(3,:),q_sol(4,:),q_sol(5,:));


P_J = JointPos(l1,l2,l3,l4,l5,q1_sol,q2_sol,q3_sol,q4_sol,q5_sol);
P_J = spline(time,P_J,0:0.005:0.8)
animation(P_J',0.005)


k_full = 1;
nr_sim = 1;
for i_sim = 1:nr_sim
    q_fw = NaN(5,1*(N_robust+1)); q_fw(:,1) = [q1_sol(1);q2_sol(1);q3_sol(1);q4_sol(1);q5_sol(1)] + 0.01*randn(5,1); % Deviation on initial positions of 0.6° std 
    dq_fw = NaN(5,1*(N_robust+1)); dq_fw(:,1) = [dq1_sol(1);dq2_sol(1);dq3_sol(1);dq4_sol(1);dq5_sol(1)] + 0.01*randn(5,1); % Deviation on initial positions of 0.6°/s std
    T_FB_fw = NaN(5,N_robust);
    contactInfo = NaN(10,3);
    for gaitcylce = 1:30
        for k = 1:N_robust+50
            if k < N_robust
                Kk = K_sol(:,(k-1)*10+1:10*k);
                refTraj_solk = refTraj_sol(:,k);
                T1_solk = T1_sol(k);
                T2_solk = T2_sol(k);
                T3_solk = T3_sol(k);
                T4_solk = T4_sol(k);
                T5_solk = T5_sol(k);
            else
                Kk = K_sol(:,(N_robust-1)*10+1:10*N_robust);
                refTraj_solk = refTraj_sol(:,N_robust);
                T1_solk = T1_sol(N_robust);
                T2_solk = T2_sol(N_robust);
                T3_solk = T3_sol(N_robust);
                T4_solk = T4_sol(N_robust);
                T5_solk = T5_sol(N_robust);
            end
            Kk = Kk;
            T_noise = 0.5*randn(5,1); % No motor noise (for now)
            relativeJointPos_k = relativeJointPos(q_fw(1,k_full),q_fw(2,k_full),q_fw(3,k_full),q_fw(4,k_full),q_fw(5,k_full));
            relativeJointVel_k = relativeJointVel(dq_fw(1,k_full),dq_fw(2,k_full),dq_fw(3,k_full),dq_fw(4,k_full),dq_fw(5,k_full));
            kinFB_k = [relativeJointPos_k; relativeJointVel_k] - refTraj_solk;
            T_FB_k =Kk*kinFB_k;
            rf = rootfinder('rf','newton',struct('x',ddq_sol,'g',f_eq_SysDyn_Error(T1_solk + T_noise(1),T2_solk + T_noise(2),T3_solk + T_noise(3),T4_solk + T_noise(4),T5_solk + T_noise(5),ddq_sol(1),ddq_sol(2),ddq_sol(3),ddq_sol(4),ddq_sol(5),dq_fw(1,k_full),dq_fw(2,k_full),dq_fw(3,k_full),dq_fw(4,k_full),dq_fw(5,k_full),q_fw(1,k_full),q_fw(2,k_full),q_fw(3,k_full),q_fw(4,k_full),q_fw(5,k_full),refTraj_solk,Kk)),struct('abstol',1e-16));
            solution_u = rf(zeros(5,1),[]);
            solution_u = full(solution_u);
%             solution_u = [result.ddq1_sol(k); result.ddq2_sol(k); result.ddq3_sol(k); result.ddq4_sol(k); result.ddq5_sol(k)];
            dq_fw(:,k_full+1) = dq_fw(:,k_full) + dt_robust*solution_u;
            q_fw(:,k_full+1) = q_fw(:,k_full) + dt_robust*dq_fw(:,k_full);
            T_FB_fw(:,k_full) = T_FB_k;
            P_J = JointPos(l1,l2,l3,l4,l5,q_fw(1,k_full+1),q_fw(2,k_full+1),q_fw(3,k_full+1),q_fw(4,k_full+1),q_fw(5,k_full+1));
            if P_J(10) < 1e-6 && k > 10
                
                q1_min = q_fw(1,k_full+1); q2_min = q_fw(2,k_full+1); q3_min = q_fw(3,k_full+1); q4_min = q_fw(4,k_full+1); q5_min = q_fw(5,k_full+1); 
%                 q1_plus = q1(:,1); q2_plus = q2(:,1); q3_plus = q3(:,1); q4_plus = q4(:,1); q5_plus = q5(:,1);
                dq1_min = dq_fw(1,k_full+1); dq2_min = dq_fw(2,k_full+1); dq3_min = dq_fw(3,k_full+1); dq4_min = dq_fw(4,k_full+1); dq5_min = dq_fw(5,k_full+1); 
%                 dq1_plus = dq1(:,1); dq2_plus = dq2(:,1); dq3_plus = dq3(:,1); dq4_plus = dq4(:,1); dq5_plus = dq5(:,1);
                rf = rootfinder('rf','newton',struct('x',state_plus_sol,'g',f_heelStrike_error(dq1_min,dq2_min,dq3_min,dq4_min,dq5_min,state_plus_sol(6),state_plus_sol(7),state_plus_sol(8),state_plus_sol(9),state_plus_sol(10),q1_min,q2_min,q3_min,q4_min,q5_min,state_plus_sol(1),state_plus_sol(2),state_plus_sol(3),state_plus_sol(4),state_plus_sol(5))),struct('abstol',1e-16));
                solution_state_plus = rf(zeros(10,1),[]);
                solution_state_plus = full(solution_state_plus);
                q_fw(:,k_full+2) = solution_state_plus(1:5);
                dq_fw(:,k_full+2) = solution_state_plus(6:10);
                contactInfo(gaitcylce,:) = [k+1   P_J(9)   P_J(10)];
                k_full = k_full + 2;
                break
            end
            k_full = k_full + 1;
        end
    end
    robustSim(i_sim).q_fw = q_fw;
    robustSim(i_sim).dq_fw = dq_fw;
    robustSim(i_sim).T_FB_fw = T_FB_fw;

    for k = 1:N_robust
        Kk = 0*K_sol(:,(k-1)*10+1:10*k);
        T_noise = 0*randn(5,1); % No motor noise (for now)
        rf = rootfinder('rf','newton',struct('x',ddq_sol,'g',f_eq_SysDyn_Error(T1_sol(k) + T_noise(1),T2_sol(k) + T_noise(2),T3_sol(k) + T_noise(3),T4_sol(k) + T_noise(4),T5_sol(k) + T_noise(5),ddq_sol(1),ddq_sol(2),ddq_sol(3),ddq_sol(4),ddq_sol(5),dq_fw(1,k),dq_fw(2,k),dq_fw(3,k),dq_fw(4,k),dq_fw(5,k),q_fw(1,k),q_fw(2,k),q_fw(3,k),q_fw(4,k),q_fw(5,k),refTraj_sol(:,k),Kk)),struct('abstol',1e-16));
        solution_u = rf(zeros(5,1),[]);
        solution_u = full(solution_u);
        q_fw(:,k+1) = q_fw(:,k) + dt_robust*dq_fw(:,k);
        dq_fw(:,k+1) = dq_fw(:,k) + dt_robust*solution_u;
    end
    nominalSim(i_sim).q_fw = q_fw;
    nominalSim(i_sim).dq_fw = dq_fw;
end
animationRobustvsNominal_vMultipleCycles(refSim,nominalSim,robustSim,nr_sim,dt_robust,l1,l2,l3,l4,l5,0,contactInfo)