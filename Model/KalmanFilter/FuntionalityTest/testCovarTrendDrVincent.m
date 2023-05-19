clc
clear all 
load SS_SOC50_Ts1.mat

A = sys.A;
B = sys.B;
C=sys.C(1,:);
N=40000;
R = 0.1;
q = .001;

Q = q*eye(size(A,1));
%Q = q*B*B';

P = idare(A',C',Q,R);
[P_infty  ,~,~] =  dare(A', C', Q, R);
ssCov=C*P*C';
ssCov1=C*P_infty*C';

P = 1e-9*eye(size(A,1));
for i=1:N,
 %   P = inv(inv(P)+C*inv(R)*C');
    P = P - P*C'*inv(C*P*C'+R)*C*P'; 
    P = A*P*A' + Q;
%    P = (P+P')/2;
    Cov(i) = C*P*C';
end

figure
semilogy(1:N,Cov,'LineWidth',2)
hold on
semilogy([1 N],[ssCov ssCov],'--','LineWidth',2)
title(['Cell Voltage Q=',num2str(q),' R=',num2str(R)])
