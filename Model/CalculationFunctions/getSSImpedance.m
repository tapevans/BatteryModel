%% getSSImpedance
% This function is used to calculate the state space matricies and its
% associated impedance.
function [A,B,C,D,Z_results] = getSSImpedance(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS)
%,sys_imp,sys_exp,sys_exp_mod
%% Initializing Inputs
SV      = SIM.SV_IC;
i_user  = SIM.i_user;
M       = SIM.M;
freq    = SIM.freq;


%% Solve for State Space (SS) Matricies
[A,B,C,D] = getSSMatricies(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,SV,i_user);


%% Calculate Impedance
Z = zeros(length(freq),1);
for i = 1:length(freq)
    s      = freq(i)*(1i);
    Z_temp = C * ((s*M - A)\B) + D;
    Z(i)   = Z_temp(1); %%%%%%%% Hard-coded for cell voltage
end
multiple = -SIM.A_c^-1;
Z_new = multiple*Z;


%% Combine results to a single variable
    % Calculations
    Z_mag = sqrt(real(Z_new).^2 + imag(Z_new).^2);
    Z_dB = 20*log10(Z_mag);
    Z_angle_deg = asind(imag(Z_new)./Z_mag);
    
    % Save
    Z_results(:,P.SS.omega)    = freq';
    Z_results(:,P.SS.Z_mag)    = Z_mag;
    Z_results(:,P.SS.Z_Re)     = real(Z_new);
    Z_results(:,P.SS.Z_Im)     = imag(Z_new);
    Z_results(:,P.SS.Z_dB)     = Z_dB;
    Z_results(:,P.SS.Z_ps_deg) = Z_angle_deg;


% %% Save the state-space system object
% % Implicit
%     sys_imp = multiple*dss(A,B,C,D,M);
% 
% % Explicit
%     sys_exp = ss(sys_imp,'explicit');
% 
% % Explicit Modified
%     % Pseudo-inverse of Mass
%         [U, S, V] = svd(SIM.M);
%         threshold = 1e-7;
%         S_cross = pinv(S,threshold);
%         M_cross = V*S_cross*U';
%     
%     % Calculate differential variable dynamics
%         A_cross = M_cross*A;
%         [~ , ol_poles] = eig(A_cross);
%         fastest_diff_pole = min(diag(real(ol_poles)));
%     
%     % Calculate Null Space matricies
%         r = rank(S_cross);
%         
%         U_null = U(:,r+1:end);
%         V_null = V(:,r+1:end);
%         U_NV_NT = U_null*V_null';
%     
%     % Calculate algebraic poles
%         LC_k = 10;
%         diag_vec = zeros(N.N_SV_tot,1);
%         algb_poles = LC_k*fastest_diff_pole*ones(length(SIM.algb_idx),1);
%         diag_vec(SIM.algb_idx) = algb_poles;
%         K = diag(diag_vec);
% 
%     % Calculate modified A and B
%         A_mod = M_cross*A - K*U_NV_NT*A;
%         B_mod = M_cross*B - K*U_NV_NT*B;
% 
% 
%     % Make SS object
%         sys_exp_mod = multiple*ss(A_mod,B_mod,C,D);


end