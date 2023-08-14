%% Asymptotic Pre Calcs
% Get the Asymptotic Kalman Gain and Error Covariance
% Provide SIM to know if state or input process noise
% Provide sys, the discrete time system
%%
function [K_infty, P_infty] = AsymptoticPreCalcs(FLAG,SIM,sys)
A_DT = sys.A;
B_DT = sys.B;
% C_DT = sys.C;
C_DT = sys.C(1,:);


[N_measur, ~] = size(C_DT);
R = SIM.R_0 * eye(N_measur);

switch FLAG.QMode
    case 1 % Input Q
        Q = SIM.Qi;
        Q_matrix = B_DT*Q*B_DT';
        % Q_matrix = Q_matrix + SIM.Q_Add * eye(size(Q_matrix));
    case 2 % State Q
        Q_matrix = SIM.Qs;
end
[P_infty  ,~,~] =  dare(A_DT', C_DT', Q_matrix, R);
K_infty = P_infty*C_DT'*inv(C_DT*P_infty*C_DT' + R);

end