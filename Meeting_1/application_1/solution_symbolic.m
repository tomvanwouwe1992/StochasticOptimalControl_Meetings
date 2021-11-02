syms P11 P21 P12 P22 K D m sigma 'real'
A = [0 1; K/m D/m];
B = [0; 0];
C = [0; 1/m];
P = [P11 P12; P21 P22];
Pdot = A*P + P*A' + C*sigma*C';