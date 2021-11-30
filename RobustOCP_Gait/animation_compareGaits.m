function [] = animation_compareGaits(P1,P2,dt)
figure(10)
P1x1 = [zeros(size(P1,1),1) P1(:,[1 3 5])];
P1y1 = [zeros(size(P1,1),1) P1(:,[2 4 6])]; 

P2x1 = [zeros(size(P2,1),1) P2(:,[1 3 5])];
P2y1 = [zeros(size(P2,1),1) P2(:,[2 4 6])]; 


P1x2 = [P1(:,[3 7 9])];
P1y2 = [P1(:,[4 8 10])]; 

P2x2 = [P2(:,[3 7 9])];
P2y2 = [P2(:,[4 8 10])]; 

Hl=line(P1x1(1,:), P1y1(1,:)); hold on
H2=line(P1x2(1,:), P1y2(1,:)); hold on

H3=line(P2x1(1,:), P2y1(1,:)); hold on
H4=line(P2x2(1,:), P2y2(1,:)); hold on

axis equal
H3=handle(H3);
H3.Color='r';
H4=handle(H4);
H4.Color='r';
 xlim([-0.6 5])
 ylim([-0.1 1.7])

for i = 1:10
for j=1:size(P1,1)
    Hl.XData=P1x1(j,:);
    Hl.YData=P1y1(j,:);
    H2.XData=P1x2(j,:);
    H2.YData=P1y2(j,:);    
    
    H3.XData=P2x1(j,:);
    H3.YData=P2y1(j,:);
    H4.XData=P2x2(j,:);
    H4.YData=P2y2(j,:);   
    pause(dt)
%     if mod(j,4) == 0 || j ==1 
%         scatter(Px1(j,1),Py1(j,1),1,'b');
%         scatter(Px1(j,2),Py1(j,2),1,'b');
%         scatter(Px1(j,4),Py1(j,4),1,'g');
%         scatter(Px2(j,1),Py2(j,1),1,'g');
%         scatter(Px2(j,2),Py2(j,2),1,'r');
%         scatter(Px2(j,3),Py2(j,3),1,'r');
%     end
end
P1x1 = P1x1 + 0.5;
P1x2 = P1x2 + 0.5;
P2x1 = P2x1 + 0.5;
P2x2 = P2x2 + 0.5;
end