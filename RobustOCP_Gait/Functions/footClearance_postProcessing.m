
close all; clc;
import casadi.*
load('resultAdapted_dt0.01_slow_Sensory_0.00_Motor_1.mat');
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
indicesP = NaN(10,10);
vecindicesP = NaN(55,1);
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
% Foot clearance

% Segment angles
q1_MX = MX.sym('q1_MX',1); q2_MX = MX.sym('q2_MX',1); q3_MX = MX.sym('q3_MX',1); q4_MX = MX.sym('q4_MX',1); q5_MX = MX.sym('q5_MX',1);
P_J_MX = JointPos(l1,l2,l3,l4,l5,q1_MX,q2_MX,q3_MX,q4_MX,q5_MX);
footClearance_MX = P_J_MX(10);
dfootClearance_dX_MX = jacobian(footClearance_MX,[q1_MX q2_MX q3_MX q4_MX q5_MX]);
f_dfootClearance_dX_MX = Function('f_dfootClearance_dX_MX',{q1_MX, q2_MX, q3_MX, q4_MX, q5_MX},{dfootClearance_dX_MX});



footClearance = NaN(81,1);
footClearance_error = NaN(81,1);
relativeJointPosMat = NaN(5,81);
for k = 1:81
    P_J = JointPos(l1,l2,l3,l4,l5,result.q1_sol(k),result.q2_sol(k),result.q3_sol(k),result.q4_sol(k),result.q5_sol(k));
    footClearance(k)= P_J(10);
    P_Jfull(:,k) = P_J;
        Lk = result.L_sol(:,k);
    Lmatk = tril(ones(10)).*Lk(indicesP);
    Pmatk = Lmatk*Lmatk';
        dfootClearance_dX = f_dfootClearance_dX_MX(result.q1_sol(k),result.q2_sol(k),result.q3_sol(k),result.q4_sol(k),result.q5_sol(k));
        footClearance_error(k) = sqrt(full(9*dfootClearance_dX*Pmatk(1:5,1:5)*dfootClearance_dX'));
        relativeJointPosMat(:,k) = relativeJointPos(result.q1_sol(k),result.q2_sol(k),result.q3_sol(k),result.q4_sol(k),result.q5_sol(k));
end

figure(10)
plot(footClearance_error); hold on;
plot(footClearance); hold on;
plot(pi/180*relativeJointPosMat(5,:));hold on;
plot(pi/180*relativeJointPosMat(2,:));

figure(2)
plot([180/pi*relativeJointPosMat(2,:) 180/pi*relativeJointPosMat(5,:)]);
