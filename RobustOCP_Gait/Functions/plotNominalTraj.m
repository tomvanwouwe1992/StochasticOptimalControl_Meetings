close all

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


cs = linspecer(10,'sequential');

for k = 1:2
    if k == 1
        load(['resultAdapted_dt0.01_slow_Sensory_0.00_Motor_1.mat']);
        %                 load(['resultNominal.mat']);
        
        stepTime = 0.8;
        strideTime = 2*stepTime;
        time = [0:0.01:stepTime];
        N_robust = 80;
        T_std = NaN(5,N_robust+1);
        PvecRel_sol = NaN(55,N_robust+1);
        q_std = NaN(5,N_robust+1);
        dq_std = NaN(5,N_robust+1);
        qRel_std = NaN(5,N_robust+1);
        dqRel_std = NaN(5,N_robust+1);
    elseif k == 2
        load(['resultAdapted_dt0.01_slow1.mat']);
        stepTime = 0.8;
        strideTime = 2*stepTime;
        time = [0:0.01:stepTime];
        N_robust = 80;
        T_std = NaN(5,N_robust+1);
        PvecRel_sol = NaN(55,N_robust+1);
        q_std = NaN(5,N_robust+1);
        dq_std = NaN(5,N_robust+1);
        qRel_std = NaN(5,N_robust+1);
        dqRel_std = NaN(5,N_robust+1);
    end
    q1_sol = result.q1_sol; q2_sol = result.q2_sol; q3_sol = result.q3_sol; q4_sol = result.q4_sol; q5_sol = result.q5_sol;
    dq1_sol = result.dq1_sol; dq2_sol = result.dq2_sol; dq3_sol = result.dq3_sol; dq4_sol = result.dq4_sol; dq5_sol = result.dq5_sol;
    ddq1_sol = result.ddq1_sol; ddq2_sol = result.ddq2_sol; ddq3_sol = result.ddq3_sol; ddq4_sol = result.ddq4_sol; ddq5_sol = result.ddq5_sol;
    T1_sol = result.T1_sol; T2_sol = result.T2_sol; T3_sol = result.T3_sol; T4_sol = result.T4_sol; T5_sol = result.T5_sol;
    
    if k == 1
        
        P_J = JointPos(l1,l2,l3,l4,l5,q1_sol,q2_sol,q3_sol,q4_sol,q5_sol);
        P_J = spline(time,P_J,0:0.005:0.8);
    end
    L_sol = result.L_sol;
    K_sol = result.K_sol;
    for j = 1:N_robust+1
        Lvec_sol = L_sol(:,j);
        Lmat_sol = tril(ones(10)).*Lvec_sol(indicesP);
        Pmat_sol = Lmat_sol*Lmat_sol';
        Pvec_sol(:,j) = Pmat_sol(vecindicesP);
        Kk = K_sol(:,(j-1)*10+1:10*j);
        %         cost_expectedEffort = cost_expectedEffort + sum(diag((Kk*feedbackTransform)*Pmat_sol*(Kk*feedbackTransform)'))*dt;
        T_std(:,j) = sqrt(diag((Kk*feedbackTransform)*Pmat_sol*(Kk*feedbackTransform)'));
        PmatRel_sol = feedbackTransform*Pmat_sol*feedbackTransform';
        PvecRel_sol(:,j) = PmatRel_sol(vecindicesP);
        q_std(:,j) = sqrt(diag(Pmat_sol(1:5,1:5)));
        dq_std(:,j) = sqrt(diag(Pmat_sol(6:10,6:10)));
        qRel_std(:,j) = sqrt(diag(PmatRel_sol(1:5,1:5)));
        dqRel_std(:,j) = sqrt(diag(PmatRel_sol(6:10,6:10)));
    end
    
    
    %     if k == 1
    %         P_J_1 = JointPos(l1,l2,l3,l4,l5,q1_sol,q2_sol,q3_sol,q4_sol,q5_sol);
    %         P_J_1 = spline(time,P_J_1,0:0.005:0.8);
    %     end
    %
    %     if k == 10
    %         P_J_10 = JointPos(l1,l2,l3,l4,l5,q1_sol,q2_sol,q3_sol,q4_sol,q5_sol);
    %         P_J_10 = spline(time,P_J_10,0:0.005:0.8);
    %     end
    
    qrel = relativeJointPos(q1_sol,q2_sol,q3_sol,q4_sol,q5_sol);
    
    figure(1)
    
    subplot(3,2,1)
    plot([time time(1:end) + time(end) + 0.00000001],180/pi*[qrel(1,:)';zeros(N_robust+1,1)],'LineWidth',1,'Color',cs(k,:)); hold on
    xlabel('time [s]')
    ylabel('angle [°]')
    title('Ankle')
    xlim([0 strideTime]);
    ylim([-35 20])
    
    subplot(3,2,2)
    plot([time time(1:end) + time(end) + 0.00000001],180/pi*[qRel_std(1,:)';zeros(N_robust+1,1)],'LineWidth',1,'Color',cs(k,:)); hold on
    xlabel('time [s]')
    ylabel('angle [°]')
    title('Ankle std')
    xlim([0 strideTime]);
    
    subplot(3,2,3)
    plot([time time(2:end) + time(end)],180/pi*[qrel(2,:)';qrel(5,2:end)'],'LineWidth',1,'Color',cs(k,:)); hold on
    xlabel('time [s]')
    ylabel('angle [°]')
    title('Knee')
    xlim([0 strideTime]);
    ylim([-105 10])
    
    
    subplot(3,2,4)
    plot([time time(2:end) + time(end)],180/pi*[qRel_std(2,:)';qRel_std(5,2:end)'],'LineWidth',1,'Color',cs(k,:)); hold on
    xlabel('time [s]')
    ylabel('angle [°]')
    title('Knee std')
    xlim([0 strideTime]);
    
    
    subplot(3,2,5)
    plot([time time(2:end) + time(end)],-180/pi*[qrel(3,:)';qrel(4,2:end)'],'LineWidth',1,'Color',cs(k,:)); hold on
    xlabel('time [s]')
    ylabel('angle [°]')
    title('Hip')
    xlim([0 strideTime]);
    ylim([-60 50])
    
    subplot(3,2,6)
    plot([time time(2:end) + time(end)],180/pi*[qRel_std(3,:)';qRel_std(4,2:end)'],'LineWidth',1,'Color',cs(k,:)); hold on
    xlabel('time [s]')
    ylabel('angle [°]')
    title('Hip std')
    xlim([0 strideTime]);
    
    
    figure(2)
    scatter(result.motorNoise_continuousVAR,result.J_torque,'d','MarkerFaceColor',cs(k,:),'MarkerEdgeColor',cs(k,:)); hold on;
    scatter(result.motorNoise_continuousVAR,result.J_expectedEffort,'MarkerFaceColor',cs(k,:),'MarkerEdgeColor',cs(k,:)); hold on;
    xlabel('Motor noise')
    ylabel('Cost')
    
    figure(3)
    subplot(3,2,1)
    plot([time time(2:end) + time(end) ],[T1_sol NaN*T5_sol(2:end) ],'LineWidth',1,'Color',cs(k,:));hold on;
    ylabel('torque [Nm]')
    xlabel('time [s]')
    title('Ankle')
    xlim([0 strideTime]);
    
    subplot(3,2,2)
    plot([time time(2:end) + time(end) ],[T_std(1,:) NaN*T_std(5,2:end) ],'LineWidth',1,'Color',cs(k,:));hold on;
    ylabel('torque [Nm]')
    xlabel('time [s]')
    title('Ankle')
    xlim([0 strideTime]);
    
    subplot(3,2,3)
    plot([time time(2:end) + time(end)],[T2_sol T5_sol(2:end)],'LineWidth',1,'Color',cs(k,:));hold on;
    title('Knee')
    ylabel('torque [Nm]')
    xlabel('time [s]')
    xlim([0 strideTime]);
    
    subplot(3,2,4)
    plot([time time(2:end) + time(end)],[T_std(2,:) T_std(5,2:end)],'LineWidth',1,'Color',cs(k,:));hold on;
    title('Knee')
    ylabel('torque [Nm]')
    xlabel('time [s]')
    xlim([0 strideTime]);
    
    subplot(3,2,5)
    plot([time time(2:end) + time(end)],[T3_sol T4_sol(2:end)],'LineWidth',1,'Color',cs(k,:));hold on;
    title('Knee')
    ylabel('torque [Nm]')
    xlabel('time [s]')
    xlim([0 strideTime]);
    
    subplot(3,2,6)
    plot([time time(2:end) + time(end)],[T_std(3,:) T_std(4,2:end)],'LineWidth',1,'Color',cs(k,:));hold on;
    title('Hip')
    ylabel('torque [Nm]')
    xlabel('time [s]')
    xlim([0 strideTime]);
end
figure(1)
legend('0.2', '0.5', '1', '2', '5', '10', '20', '30', '40', '50')

figure(2)
legend('Expected feedforward cost','Expected feedback cost')

figure(3)
legend('0.2', '0.5', '1', '2', '5', '10', '20', '30', '40', '50')

animation(P_J',0.005)

animation_compareGaits(P_J_1',P_J_10',0.005)