function X = rr_fkine(q, l1, l2)
% RR_FKINE - Kinematics of a 2D RR robot.
%    [x;y] = rr_fkine(Q,L1,L1) computes the forward kinematics for a 2D RR
%    robot with link lengths L1, L2 and joint angles Q=[q1;q1].

x = l1*cos(q(1)) + l2*cos(q(1)+q(2));
y = l1*sin(q(1)) + l2*sin(q(1)+q(2));

X=[x;y];