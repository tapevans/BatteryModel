%% Asymptotic Pre Calcs
% Get the Asymptotic Kalman Gain and Error Covariance
% Provide SIM to know if state or input process noise
% Provide sys, the discrete time system
%%
function [K_infty, P_infty] = AsymptoticPreCalcs(FLAG,SIM,sys)
A_DT = sys.A;
B_DT = sys.B;
C_DT = sys.C;


[N_measur, ~] = size(sys.C);
R = SIM.R_0 * eye(N_measur);

switch FLAG.QMode
    case 1 % Input Q
        Q = SIM.Q3i;
        [P_infty  ,~,~] =  dare(A_DT', C_DT', B_DT*Q*B_DT' ,R);
    case 2 % State Q
        Q = SIM.Q3s;
        [P_infty  ,~,~] =  dare(A_DT', C_DT', Q, R);
end

K_infty = P_infty*C_DT'*inv(C_DT*P_infty*C_DT' + R);

end