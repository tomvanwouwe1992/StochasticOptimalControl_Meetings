clear all; close all; clc;
dt = 0.01;
D = 0.5; %(step-length)
stepTime = 0.8;
% motorNoise_continuousVAR_vec = [0.01 0.05 0.1 0.2 0.5 1 2 5 10 20 50];
motorNoise_continuousVAR_vec = [0.2 0.5 1 2 5 10 20 30 40 50];


% Note: example here, we input a converged solution to the same problem
% setup (no sensory noise, and same motor noise)
for k = 1:1

    motorNoise_continuousVAR = 1*motorNoise_continuousVAR_vec(k);
    sensoryNoise_continuousVAR = 0*0.01*[pi/180*1*ones(5,1); pi/180*3*ones(5,1)].^2;

    IGName =    [pwd '\Results\resultAdapted_dt0.01_slow_Sensory_0.00_Motor_1.mat'];
    saveName =  [pwd '\Results\resultAdapted_dt0.01_slow_Sensory_0.00_Motor_1.mat'];
    diaryName =  [pwd '\Results\diaryAdapted_dt0.01_slow_Sensory_0.00_Motor_1.txt'];
    
	mainOCP_sensorimotorNoise_withSensory_fcn(saveName,IGName,diaryName,dt,sensoryNoise_continuousVAR,motorNoise_continuousVAR,D,stepTime);
end



