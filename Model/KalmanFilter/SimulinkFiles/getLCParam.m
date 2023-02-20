%% Get LC Parameters
%
%
%
function [LC_param] = getLCParam(SIM,N,P,FLAG)

filename = getImpulseFilename(FLAG);
simsys = load(filename);
i_user_sim = -1;
[A,B,C,D] = getSSfromSV(simsys.AN , simsys.CA , simsys.SEP , simsys.EL , simsys.SIM , simsys.CONS , simsys.P , simsys.N , simsys.FLAG , simsys.PROPS, simsys.SIM.SV_IC , i_user_sim);
Mass = simsys.SIM.M;
algb_idx = simsys.SIM.algb_idx;

%% Calculate Parameters for Battery System
    % Pseudo-inverse of Mass
        [U, S, V] = svd(Mass);
        threshold = 1e-7;
        S_cross = pinv(S,threshold);
        M_cross = V*S_cross*U';
    
    % Calculate differential variable dynamics
        A_cross = M_cross*A;
        [~ , ol_poles] = eig(A_cross);
        fastest_diff_pole = min(diag(real(ol_poles)));
    
    % Calculate Null Space matricies
        r = rank(S_cross);
        
        U_null  = U(:,r+1:end);
        V_null  = V(:,r+1:end);
        U_NV_NT = U_null*V_null';
    
    % Calculate algebraic poles
        LC_k = 1e1;
        diag_vec = zeros(N.states,1);
        algb_poles = LC_k*fastest_diff_pole*ones(length(algb_idx),1);
        diag_vec(algb_idx) = algb_poles;
        K = -diag(diag_vec);

        
%% Out
    LC_param.M_cross = M_cross;
    LC_param.U_NV_NT = U_NV_NT;
    LC_param.K       = K;

end

%%%%%%%%%%%%%%%%%%% OLD %%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Get SS
% % switch FLAG.SlinkModel
% %     case 1 % 3 state plant
% %         N.states = 3;
% %         FLAG.C_mode = 4;
% %         [A,B,C,D] = getAll_SS3(SIM,N,P,FLAG);
% %         Mass = SIM.Mass3;
% %         algb_idx = SIM.algb_idx3;
% %     case 2 % 3 state plant with input noise
% %         N.states = 3;
% %         FLAG.C_mode = 4;
% %         [A,B,C,D] = getAll_SS3(SIM,N,P,FLAG);
% %         Mass = SIM.Mass3;
% %         algb_idx = SIM.algb_idx3;
% %     case 3 % 3 state plant with state noise
% %         N.states = 3;
% %         FLAG.C_mode = 4;
% %         [A,B,C,D] = getAll_SS3(SIM,N,P,FLAG);
% %         Mass = SIM.Mass3;
% %         algb_idx = SIM.algb_idx3;
% % end
%     N.states = 4;
%     FLAG.C_mode = 5;
%     [A,B,C,D] = getAll_SS(SIM,N,P,FLAG);
%     Mass = SIM.Mass;
%     algb_idx = [];
