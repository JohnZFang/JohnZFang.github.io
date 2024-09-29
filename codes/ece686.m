% A = diag([1.3, 1.1, 1.1]);
% B = [0.7 0;
%     0.7 0.7;
%     0 0.7];
% DeltaA  =  [-0.0102 -0.001 -0.0183; 0.0043 -0.0002 0.0167; -0.0032 0.0051 -0.0022];
% DeltaB  =  [-0.0072 0.0013; 0.0021 0.0065; 0.0019 0.0016];
% Q = 10*eye(3);
% R  =  eye(2);
% n = 3;
% m = 2;
% sum = diag([1, 1]);
% esp = 
% e = 
% rho = 
% epsilon1 = 0.035;
% epsilon2 = 0.01;
% Phi_1 = 1.4;
% Phi_2 = 0.8;
% x_0 = [-1.4 1.5 -1.8]';
% b = 1;
% a = 0;
% tau_1 = 0;
% tau_2 = 0;
clc;
clear all;
A = [0.81 0.15;
     0.31 0.7];
C2 = [2 2.5];
gamma1 = 2.5;

Q = diag([0.09 0.09]);
R = 0.05;
M = 40;
mu = 0.9;
k = M;
P1 = cell(41,1);
P2 = cell(41,1);
L_star = cell(41,1);
L = cell(41,1);
W = cell(41,1);
P1{41} = 0;
P2{41} = 0;
L_star{41} = 0;
L{41} = 0;
W{41} = 0;

A_L = cell(41,1);
A_W = cell(41,1);
Delta = cell(41,1);
Lambda = cell(41,1);
normP1 = [];
normP2 = [];
for i=40:-1:1
    A_L{i} = A+mu*L_star{i+1}*C2;
    A_W{i} = A+W{i+1};
    Delta{i} = eye(2)-1/(gamma1)^2*P1{i+1};
    Lambda{i} = R+C2*P2{i+1}*C2';
    L_star{i} = -A_W{i}*P2{i+1}*C2'/(Lambda{i});
    W{i} = ((gamma1)^2*P1{i+1}-eye(2))\A_L{i};
    %L{i} = -A*P2{i+1}*C2'/(Lambda{i});
    P1{i} = A_L{i}'*P1{i+1}*A_L{i}+eye(2)+(mu-mu^2)*C2'*L_star{i+1}'*P1{i+1}*L_star{i+1}*C2+1/(gamma1^2)*A_L{i}'*P1{i+1}*inv(Delta{i})*P1{i+1}*A_L{i};%whether it's L_star or L needs to verify
    P2{i} = A_W{i}*P2{i+1}*A_W{i}'+Q-mu*A_W{i}*P2{i+1}*C2'*inv(Lambda{i})*C2*P2{i+1}*A_W{i}';
    normP1(i) = norm(P1{i}-P1{i+1},2);
    normP2(i) = norm(P2{i}-P2{i+1},2);
end
%fig1 = plot(1:40,normP1);
% fig2 = plot(1:40,normP2);
% sigma = [0.9;1.2;0.5;1.3;0.6;1.5;0.8;2.7;1;0.7];
% delta = 0.5;
% T = ones(10,1)*16;
L = 300;
x = cell(1,L);
x{1} = [1;1];
X_est = cell(1,L);
X_est{1} = [1;1];
K = cell2mat(L_star);
K = [-0.2564; -0.2070];
for t=1:L-1
    x{t+1} = A*x{t}+abs(0.016*sin(0.1*t))*ones(2,1)+mvnrnd([0;0], Q)';
    y{t} = C2*x{t}+mvnrnd(0, R);
    X_est{t+1} = (A+mu*K(1:2,1)*C2)*X_est{t}-mu*K(1:2,1)*y{t};
end
x = cell2mat(x);
X_est = cell2mat(X_est);
plot(1:L, x(1,:),'-b',1:L,X_est(1,:),'-g');%1:T, newy,'-ko', 
legend('true x1','estimation');
title('state x1');
plot(1:L, x(2,:),'-b',1:L,X_est(2,:),'-g');%1:T, newy,'-ko', 
legend('true x2','estimation');
title('state x2');
