%% Variable Estimator
%
%
%%
function [x_hat,K_k,P_k_pre] = VarEstimator(sys,FLAG,SIM,P,x_hat_0,u,z)
%% Perform Pre-Calcs
[K_infty, P_infty] = AsymptoticPreCalcs(FLAG,SIM,sys);

u_k = reshape(u,1,[]);
z_k = reshape(z,1,[]);


%% Initialize Variables
N_steps = length(u);
[N_measur, N_states] = size(sys.C);

x_var     = zeros(N_states,N_steps);   % Estimator States
x_var_pre = zeros(N_states,N_steps);   % Estimator States predict phase
y_tilde_k = zeros(N_measur,N_steps);   % Estimation Error
S_k       = zeros(N_measur,N_measur,N_steps); % Error Covariance
K_k       = zeros(N_states,N_measur,N_steps); % Kalman Gains
P_k       = zeros(N_states,N_states,N_steps); % Error Covariance
P_k_pre   = zeros(N_states,N_states,N_steps); % Error Covariance predict phase

confidence = 0.5; % Range from [0,1]
P_k_pre(:,:,1) = confidence*eye(N_states);
P_k(:,:,1) = confidence*eye(N_states);
% P_k(:,:,1)     = P_infty;
% P_k_pre(:,:,1) = P_infty;

K_k(:,:,1) = K_infty;

A_DT = sys.A;
B_DT = sys.B;
C_DT = sys.C;

x_var(:,1) = x_hat_0;

switch FLAG.QMode
    case 1 % Input Q
        Q = SIM.Qi;
    case 2 % State Q
        Q = SIM.Qs;
end
R = SIM.R_0 * eye(N_measur);
% SIM.R = (diag(R))';

% if FLAG.UseWrongIC
%     z_init = SIM.y_0_FOM_Offset(P.cell_voltage,1);
% else
    z_init = SIM.y_0_FOM(P.cell_voltage,1);
% end

% z_init = SIM.y_0_FOM(P.cell_voltage,1);


reset_steps = SIM.ResetStep; % Reset P back to the original P


%% Run Estimation
switch FLAG.QMode
    case 1 % Input Q
        for i = 2:N_steps
        % Predict Phase
            x_var_pre(:,i) = A_DT * x_var(:,i-1)  +  B_DT * u_k(:,i-1);
            P_k_pre(:,:,i) = A_DT * P_k(:,:,i-1) * A_DT'  +  B_DT*Q*B_DT';
            %i
            %mod(i,100)
            if mod(i,reset_steps) == 0
                P_k_pre(:,:,i) = P_k_pre(:,:,1);
                P_k(:,:,i)     = P_k(:,:,1);
                K_k(:,:,i)     = K_k(:,:,1);
            end
        
    %     % Actuate System
    %         x_EST_sys(:,i) = A_DT * x_EST_sys(:,i-1)  +  B_DT * u_k(i-1)  +  w_k(:,i-1);
    %         z_k(:,i) = C * x_EST_sys(:,i)  +  v_k(:,i);
    
        % Update (Correction) Phase
            y_tilde_k(:,i) = z_k(:,i) - (C_DT * x_var_pre(:,i) + z_init);
            S_k(:,:,i) = C_DT * P_k_pre(:,:,i) * C_DT' + R;
            K_k(:,:,i) = P_k_pre(:,:,i) * C_DT' * inv(S_k(:,:,i));
            x_var(:,i) = x_var_pre(:,i) + K_k(:,:,i) * y_tilde_k(:,i);
            P_k(:,:,i) = (eye(N_states) - K_k(:,:,i) * C_DT) * P_k_pre(:,:,i);
        end
    case 2 % State Q
        for i = 2:N_steps
        % Predict Phase
            x_var_pre(:,i) = A_DT * x_var(:,i-1)  +  B_DT * u_k(:,i-1);
            P_k_pre(:,:,i) = A_DT * P_k(:,:,i-1) * A_DT'  +  Q;
        
    %     % Actuate System
    %         x_EST_sys(:,i) = A_DT * x_EST_sys(:,i-1)  +  B_DT * u_k(i-1)  +  w_k(:,i-1);
    %         z_k(:,i) = C * x_EST_sys(:,i)  +  v_k(:,i);
    
        % Update (Correction) Phase
            y_tilde_k(:,i) = z_k(:,i) - C_DT * x_var_pre(:,i);
            S_k(:,:,i) = C_DT * P_k_pre(:,:,i) * C' + R;
            K_k(:,:,i) = P_k_pre(:,:,i) * C_DT' * inv(S_k(:,:,i));
            x_var(:,i) = x_var_pre(:,i) + K_k(:,:,i) * y_tilde_k(:,i);
            P_k(:,:,i) = (eye(N_states) - K_k(:,:,i) * C_DT) * P_k_pre(:,:,i);
        end
end

x_hat = x_var;

% %%%% Converging Test
% error_P_k_pre = P_k_pre - P_infty;
% for i = 1:N_steps
%     norm_error_P_k_pre(:,i) = norm(error_P_k_pre(:,:,i));
% end
% 
% error_P_k = P_k - P_infty;
% for i = 1:N_steps
%     norm_error_P_k(:,i) = norm(error_P_k(:,:,i));
% end
% 
% error_K_k = K_k - K_infty;
% for i = 1:N_steps
%     norm_error_K_k(:,i) = norm(error_K_k(:,:,i));
% end
% 
% figure
% title('Convergence to Asymptotic Pre')
% yyaxis left
% plot(1:1:N_steps,norm_error_P_k_pre,'Linewidth',2)
% ylabel('|Error P_k|')
% 
% yyaxis right
% plot(1:1:N_steps,norm_error_K_k,'Linewidth',2)
% ylabel('|Error K_k|')
% 
% xlabel('Number of Discrete Steps')
% % xlim([1,NS])
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% title('Convergence to Asymptotic')
% yyaxis left
% plot(1:1:N_steps,norm_error_P_k,'Linewidth',2)
% ylabel('|Error P_k|')
% 
% yyaxis right
% plot(1:1:N_steps,norm_error_K_k,'Linewidth',2)
% ylabel('|Error K_k|')
% 
% xlabel('Number of Discrete Steps')
% % xlim([1,NS])
end