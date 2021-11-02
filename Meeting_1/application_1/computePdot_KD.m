function Pdot = computePdot_KD(KD,m,P,sigma_W)
A = [0 1; KD(1)/m KD(2)/m];
C = [0;1/m];
Pdot = A*P + P*A' + C*sigma_W*C';
Pdot = [Pdot(1,1);Pdot(2,1);Pdot(2,2)];