%% Ho-Kalman Model Reduction
% This code takes the impulse response of a system and converts it into a 
% reduced order model.
% Inputs:
    % Impulse Data
    % Tolerance or Number of States
% Output:
    % State-Space matricies for the reduced model

function [A_r,B_r,C_r,D_r] = HoKalmanReduction(g_k,r)
% Create Hankle Matrix (H)
    H = myHankle(g_k);

% Extract ROM
    % SVD
    [Nrows, Ncolms] = size(H);
    N_outputs = Nrows/Ncolms;
    N_in = 1; %%% !!! Hardcoded but always true for batteries
    
    [U,S,V] = svd(H);
%     r = rank(S);
%     r = 3;
    
    U_colm = U(:,1:r);
    S_colm = S(1:r,1:r);
    V_colm = V(:,1:r);
    
%     if localFLAG.BalRed % If using balanced reduction
        S_sqrt = S_colm.^0.5;
        obsv_r = U_colm*S_sqrt;
        cont_r = S_sqrt*V_colm';
%     else % If using normal reduction
%         obsv_r = U_colm;
%         cont_r = S_colm*V_colm';
%     end
    
% Determine A B C D
    % C_r
    C_r = obsv_r(1:N_outputs,:   );
    
    % B_r
    B_r = cont_r(:        ,N_in);
    
    % A_r
    P   = obsv_r(1:end-N_outputs,:);
    P_p = obsv_r(N_outputs+1:end,:);
    threshold = 1e-10;
    P_inv = pinv(P,threshold);
    A_r = P_inv*P_p;
    
    % D_r ~ !!!!! I'm assuming this is correct
    D_r = zeros(N_outputs,N_in);


end