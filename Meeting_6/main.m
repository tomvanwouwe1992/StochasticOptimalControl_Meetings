wM_std = [0.05]; %0.01 0.025 0.05 0.075 0.1
wPq_std = [3e-4]; % 6e-4];% 1.2e-3];% 2.4e-3];
wPqdot_std = [2.4e-3];% 4.8e-3];% 9.6e-3];
listing = dir();



forceField = 0;

target = 'CIRCLE';
saveName = ['result_' target '_forceField_' num2str(forceField) '_' num2str(wM_std) '_' num2str(wPq_std) '_' num2str(wPqdot_std) '.mat'];
[result,succes] = OCP_ReachingTorqueDriven(target,forceField,wM_std,wPq_std,wPqdot_std,[]);
if succes == true
    save(saveName,'result');
end


target = 'BAR';
saveName = ['result_' target '_forceField_' num2str(forceField) '_' num2str(wM_std) '_' num2str(wPq_std) '_' num2str(wPqdot_std) '.mat'];
% guessName = ['result_CIRCLE_forceField_' num2str(forceField) '_' num2str(0.05) '_' num2str(wPq_std) '_' num2str(wPqdot_std) '.mat'];

[result,succes] = OCP_ReachingTorqueDriven(target,forceField,wM_std,wPq_std,wPqdot_std,[]);
if succes == true
    save(saveName,'result');
end


target = 'OBSTACLE';
saveName = ['result_' target '_forceField_' num2str(forceField) '_' num2str(wM_std) '_' num2str(wPq_std) '_' num2str(wPqdot_std) '.mat'];
% guessName = ['result_' 'CIRCLE' '_forceField_' num2str(forceField) '_' num2str(wM_std) '_' num2str(wPq_std) '_' num2str(wPqdot_std) '.mat'];

[result,succes] = OCP_ReachingTorqueDriven(target,forceField,wM_std,wPq_std,wPqdot_std,[]);
if succes == true
   save(saveName,'result');
end



forceField = 200;

target = 'CIRCLE';
saveName = ['result_' target '_forceField_' num2str(forceField) '_' num2str(wM_std) '_' num2str(wPq_std) '_' num2str(wPqdot_std) '.mat'];
% guessName = ['result_' target '_forceField_' num2str(forceField) '_' num2str(0.05) '_' num2str(wPq_std) '_' num2str(wPqdot_std) '.mat'];

[result,succes] = OCP_ReachingTorqueDriven(target,forceField,wM_std,wPq_std,wPqdot_std,[]);
if succes == true
    save(saveName,'result');
end

