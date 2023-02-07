%% Get LC Parameters
%
%
%
function [LC_param] = getLCParam(SIM,N,P)
%% Convert Inputs
%     AN    = sim.AN;
%     CA    = sim.CA;
%     SEP   = sim.SEP;
%     EL    = sim.EL;
%     SIM   = sim.SIM;
%     CONS  = sim.CONS;
%     P     = sim.P;
%     N     = sim.N;
%     FLAG  = sim.FLAG;
%     PROPS = sim.PROPS;

%%
% % i_user_sim = SIM.i_user_amp;
% i_user_sim = 1; %%%%%%%%%% Guessing here but I think this is right
% SV_IC = SIM.SV_IC;
% 
% [A,B,C,D] = getSSfromSV(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,SV_IC,i_user_sim);

%% Get SS
[A,B] = getCTAandB(SIM,N,P);


%% Calculate Parameters for Battery System
    % Pseudo-inverse of Mass
        [U, S, V] = svd(SIM.Mass);
        threshold = 1e-7;
        S_cross = pinv(S,threshold);
        M_cross = V*S_cross*U';
    
    % Calculate differential variable dynamics
        A_cross = M_cross*A;
        [~ , ol_poles] = eig(A_cross);
        fastest_diff_pole = min(diag(real(ol_poles)));
    
    % Calculate Null Space matricies
        r = rank(S_cross);
%         r = 2;
        
        U_null = U(:,r+1:end);
        V_null = V(:,r+1:end);
        U_NV_NT = U_null*V_null';
%         U_NV_NT = zeros(5);
%         U_NV_NT(2,2) = -1;
%         U_NV_NT(5,5) =  1;
    
    % Calculate algebraic poles
        LC_k = 1e1;
        diag_vec = zeros(N.states,1);
        algb_poles = LC_k*fastest_diff_pole*ones(length(SIM.algb_idx),1);
        diag_vec(SIM.algb_idx) = algb_poles;
        K = diag(diag_vec);

        
%% Out
LC_param.M_cross = M_cross;
LC_param.U_NV_NT = U_NV_NT;
LC_param.K       = K;

end