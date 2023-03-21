%% Testing Rob's Spline Function
clear all; close all; clc;
%% Test Vector
X_vec = 0.05:0.001:.28;

%%
tic
for i = 1:length(X_vec)
    D_aniso_1 = [0.109310762, 1.41E-13;...
                 0.121886336, 6.11E-13;...
                 0.158645707, 1.22E-13;...
                 0.195405079, 1.29E-13;...
                 0.231197098, 1.89E-13;...
                 0.272793229 ,4.15E-15];

    af1 = fit(D_aniso_1(:,1),D_aniso_1(:,2),'smoothingspline');

    %%
    D_t(i) = af1(X_vec(i));
end
toc

%% From af1
% for i = 1:length(X_vec)
%     if X_vec(i) < af1.p.breaks(1)
%         D_t_test(i) = NaN;
%     elseif X_vec(i) >= af1.p.breaks(1) && X_vec(i) < af1.p.breaks(2)
%         D_t_test(i) = MyPolyfit(af1.p.coefs(1,:),X_vec(i));
%     elseif X_vec(i) >= af1.p.breaks(2) && X_vec(i) < af1.p.breaks(3)
%         D_t_test(i) = MyPolyfit(af1.p.coefs(2,:),X_vec(i));
%     elseif X_vec(i) >= af1.p.breaks(3) && X_vec(i) < af1.p.breaks(4)
%         D_t_test(i) = MyPolyfit(af1.p.coefs(3,:),X_vec(i));
%     elseif X_vec(i) >= af1.p.breaks(4) && X_vec(i) < af1.p.breaks(5)
%         D_t_test(i) = MyPolyfit(af1.p.coefs(4,:),X_vec(i));
%     elseif X_vec(i) >= af1.p.breaks(5) && X_vec(i) < af1.p.breaks(6)
%         D_t_test(i) = MyPolyfit(af1.p.coefs(5,:),X_vec(i));
%     else
%         D_t_test(i) = NaN;
%     end
% end

%% From af1
tic
for i = 1:length(X_vec)
    if X_vec(i) < af1.p.breaks(1)
        D_t_test(i) = NaN;
    elseif X_vec(i) > af1.p.breaks(6)
        D_t_test(i) = NaN;
    else
        idx = find(X_vec(i) > af1.p.breaks,1);
        D_t_test(i) = MyPolyfit(af1.p.coefs(idx,:),X_vec(i));
    end
end
toc
%% Plot
figure
hold on
semilogy(X_vec,D_t)
semilogy(X_vec,D_t,'o')


%%
function f = MyPolyfit(coef,x)
p1 = coef(1);
p2 = coef(2);
p3 = coef(3);
p4 = coef(4);

f = p1*x^3 + p2*x^2 + p3*x + p4;

end