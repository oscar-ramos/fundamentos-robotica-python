function rr_plot(q,l1,l2)
% RR_PLOT - Plot a 2D RR robot.
%    rr_plot(Q,L1,L1) plots a 2D revolute-revolute robot with link lengths
%    L1, L2 and where the joint angles are Q=[q1;q1].

q1=q(1); q2=q(2);

x1=l1*cos(q1);
y1=l1*sin(q1);
x2=x1+l2*cos(q1+q2);
y2=y1+l2*sin(q1+q2);

line([0 x1 x2],[0 y1 y2],'LineWidth',2), axis equal, grid on, hold on
plot(x1,y1,'o','LineWidth',2,'Color','r')
plot(0,0,'o','LineWidth',2,'Color','r')