index = 1;

q1_sol = q_scale(1).*output(index:index + N_robust)'; index = index + N_robust + 1; result.q1_sol = q1_sol;
q2_sol = q_scale(2).*output(index:index + N_robust)'; index = index + N_robust + 1; result.q2_sol = q2_sol;
q3_sol = q_scale(3).*output(index:index + N_robust)'; index = index + N_robust + 1; result.q3_sol = q3_sol;
q4_sol = q_scale(4).*output(index:index + N_robust)'; index = index + N_robust + 1; result.q4_sol = q4_sol;
q5_sol = q_scale(5).*output(index:index + N_robust)'; index = index + N_robust + 1; result.q5_sol = q5_sol;

dq1_sol = dq_scale(1).*output(index:index + N_robust)'; index = index + N_robust + 1; result.dq1_sol = dq1_sol;
dq2_sol = dq_scale(2).*output(index:index + N_robust)'; index = index + N_robust + 1; result.dq2_sol = dq2_sol;
dq3_sol = dq_scale(3).*output(index:index + N_robust)'; index = index + N_robust + 1; result.dq3_sol = dq3_sol;
dq4_sol = dq_scale(4).*output(index:index + N_robust)'; index = index + N_robust + 1; result.dq4_sol = dq4_sol;
dq5_sol = dq_scale(5).*output(index:index + N_robust)'; index = index + N_robust + 1; result.dq5_sol = dq5_sol;

ddq1_sol = ddq_scale(1).*output(index:index + N_robust)'; index = index + N_robust + 1; result.ddq1_sol = ddq1_sol;
ddq2_sol = ddq_scale(2).*output(index:index + N_robust)'; index = index + N_robust + 1; result.ddq2_sol = ddq2_sol;
ddq3_sol = ddq_scale(3).*output(index:index + N_robust)'; index = index + N_robust + 1; result.ddq3_sol = ddq3_sol;
ddq4_sol = ddq_scale(4).*output(index:index + N_robust)'; index = index + N_robust + 1; result.ddq4_sol = ddq4_sol;
ddq5_sol = ddq_scale(5).*output(index:index + N_robust)'; index = index + N_robust + 1; result.ddq5_sol = ddq5_sol;

T1_sol = T_scale(1).*output(index:index + N_robust)'; index = index + N_robust + 1; result.T1_sol = T1_sol;
T2_sol = T_scale(2).*output(index:index + N_robust)'; index = index + N_robust + 1; result.T2_sol = T2_sol;
T3_sol = T_scale(3).*output(index:index + N_robust)'; index = index + N_robust + 1; result.T3_sol = T3_sol;
T4_sol = T_scale(4).*output(index:index + N_robust)'; index = index + N_robust + 1; result.T4_sol = T4_sol;
T5_sol = T_scale(5).*output(index:index + N_robust)'; index = index + N_robust + 1; result.T5_sol = T5_sol;


refTraj_sol = refTraj_scale.*reshape(output(index:index + 10*(N_robust+1) - 1),10,N_robust+1); index = index + 10*(N_robust + 1); result.refTraj_sol = refTraj_sol;

K_sol = repmat(K_scale,1,N_robust+1).*reshape(output(index:index + 5*10*(N_robust+1) - 1),5,10*(N_robust+1)); index = index + 5*10*(N_robust + 1); result.K_sol = K_sol;

udx_sol = repmat(udx_scale,1,N_robust+1).*reshape(output(index: index + 10*10*(N_robust+1) - 1),10,10*(N_robust+1)); index = index + 10*10*(N_robust + 1); result.udx_sol = udx_sol;

M_sol = reshape(output(index: index + 10*10*N_robust - 1),10,10*N_robust); index = index + 10*10*N_robust; result.M_sol = M_sol;

L_sol = Lscale_vec.*reshape(output(index: index + 55*(N_robust+1) - 1),55,N_robust+1); index = index + 55*(N_robust + 1); result.L_sol = L_sol;

udw_sol = repmat(udw_scale,1,N_robust+1).*reshape(output(index: index + 10*5*(N_robust+1) - 1),10,5*(N_robust+1)); index = index + 10*5*(N_robust + 1); result.udw_sol = udw_sol;

L_precontact_sol = Lscale(vecindicesP).*output(index: index + 55 - 1); index = index + 55; result.L_precontact_sol = L_precontact_sol;

L_postcontact_sol = Lscale(vecindicesP).*output(index: index + 55 - 1); index = index + 55; result.L_postcontact_sol = L_postcontact_sol;

sigma_postcontact_sol = reshape(output(index: index + 10*21 - 1),10,21); index = index + 10*21; result.sigma_postcontact_sol = sigma_postcontact_sol;

initVarPos_sol = Pscale_diag(1:5).*output(index:index + 5 - 1); index = index + 5; result.initVarPos_sol = initVarPos_sol;

initVarVel_sol = Pscale_diag(6:10).*output(index:index + 5 - 1); index = index + 5; result.initVarVel_sol = initVarVel_sol;

initCovar_sol = Pscale_vec(covarianceIndices).*output(index:index + 45 - 1); index = index + 45; result.initCovar_sol = initCovar_sol;