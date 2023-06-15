%% Variable Estimator
%
%
%%
function [x_hat,K_k,P_k_pre] = VarEstimator(sys,FLAG,SIM,P,x_hat_02,u,z)
    A_DT = sys.A;
    B_DT = sys.B;
    C_DT    = sys.C;
    C_DT_CV = sys.C(P.cell_voltage,:);

    [N_measur, N_states] = size(C_DT_CV);
    [N_measur_all, ~]    = size(C_DT);


%% Perform Pre-Calcs
    if FLAG.DoPreCalc
        [K_infty, P_infty] = AsymptoticPreCalcs(FLAG,SIM,sys);
    else
        K_infty = 0.5*ones(N_states,N_measur);
        P_infty = 0.5*eye(N_states);
    end
    K_infty = K_infty(:,P.cell_voltage);
    
    u_k = reshape(u,1,[]);
    z_k = reshape(z,1,[]);


%% Initialize Variables
    N_steps = length(u);
    
    x_var     = zeros(N_states,N_steps);   % Estimator States
    x_var_pre = zeros(N_states,N_steps);   % Estimator States predict phase
    y_tilde_k = zeros(N_measur,N_steps);   % Estimation Error
    S_k       = zeros(N_measur,N_measur,N_steps); % Error Covariance
    K_k       = zeros(N_states,N_measur,N_steps); % Kalman Gains
    P_k       = zeros(N_states,N_states,N_steps); % Error Covariance
    P_k_pre   = zeros(N_states,N_states,N_steps); % Error Covariance predict phase


%% Calculate an Initial Error Covariance
    % sigma = 3e-1;
    % sigma = 1e-4;
    % sigma = 1e-5;
    sigma = 1e-7;
    y_des = sigma * ones(N_measur_all,1);
    
    % Try using inverse
    x_hat_0_test = C_DT \ y_des;

    % Calculate P_k|k
    error = zeros(size(x_hat_0_test)) - x_hat_0_test;

    P_kk = nan(length(x_hat_0_test));
    for i = 1:length(x_hat_0_test)
        for j = 1:length(x_hat_0_test)
            [P_kk(i,j)] = calcPopCovar(error(i,:), error(j,:)); %%%%%%May be able to change this so we don't need the calcPopCovar here
        end
    end


%% Initialize Kalman Gain and Error Covariance
% Kalman Gain    
    % sensor_confidence = 1; % Range from [0,1]
    sensor_confidence = 0; % Range from [0,1]
    K_k_0     (:,:,1) = sensor_confidence * ones(N_states , N_measur);

% Error Covariance
    % % state_variance    = 1e-1; 
    % state_variance    = 1; 
    % P_k_0     (:,:,1) = state_variance*eye(N.states);
    % P_k_pre_0 (:,:,1) = state_variance*eye(N.states);

    P_k_0     (:,:,1) = P_kk;
    P_k_pre_0 (:,:,1) = P_kk;


% confidence = 1.0; % Range from [0,1]
% P_k_pre(:,:,1) = confidence*eye(N_states);
% P_k(:,:,1) = confidence*eye(N_states);
% % P_k(:,:,1)     = P_infty;
% % P_k_pre(:,:,1) = P_infty;
% 
% % K_k(:,:,1) = K_infty;
% % K_k(:,:,1) = 2*K_infty;
% K_k(:,:,1) = eye(size(K_infty));


%%
    % x_var(:,1) = x_hat_0;
    x_var(:,1) = x_hat_0_test;


%%
    switch FLAG.QMode
        case 1 % Input Q
            Q = SIM.Qi;
            Q_matrix = B_DT*Q*B_DT';
            Q_matrix = Q_matrix + SIM.Q_Add * eye(size(Q_matrix));
        case 2 % State Q
            Q_matrix = SIM.Qs;
    end
    R = SIM.R_0 * eye(N_measur);
    % SIM.R = (diag(R))';
    
    if FLAG.UseWrongIC_y
        z_init = SIM.y_0_FOM_Offset(P.cell_voltage,1);
    else
        z_init = SIM.y_0_FOM(P.cell_voltage,1);
    end
    
    % z_init = SIM.y_0_FOM(P.cell_voltage,1);
    
    % reset_steps = SIM.ResetStep; % Reset P back to the original P


%% Run Estimation
    switch FLAG.QMode
        case 1 % Input Q
            for i = 2:N_steps
            % Predict Phase
                x_var_pre(:,i) = A_DT * x_var(:,i-1)  +  B_DT * u_k(:,i-1);
                % P_k_pre(:,:,i) = A_DT * P_k(:,:,i-1) * A_DT'  +  B_DT*Q*B_DT';
                P_k_pre(:,:,i) = A_DT * P_k(:,:,i-1) * A_DT'  +  Q_matrix;
        
            % Update (Correction) Phase
                y_tilde_k(:,i) = z_k(:,i) - (C_DT_CV * x_var_pre(:,i) + z_init);
                S_k(:,:,i) = C_DT_CV * P_k_pre(:,:,i) * C_DT_CV' + R;
                K_k(:,:,i) = P_k_pre(:,:,i) * C_DT_CV' * inv(S_k(:,:,i));
                x_var(:,i) = x_var_pre(:,i) + K_k(:,:,i) * y_tilde_k(:,i);
                P_k(:,:,i) = (eye(N_states) - K_k(:,:,i) * C_DT_CV) * P_k_pre(:,:,i);
            end
        case 2 % State Q
            for i = 2:N_steps
            % Predict Phase
                x_var_pre(:,i) = A_DT * x_var(:,i-1)  +  B_DT * u_k(:,i-1);
                P_k_pre(:,:,i) = A_DT * P_k(:,:,i-1) * A_DT'  +  Q;
        
            % Update (Correction) Phase
                y_tilde_k(:,i) = z_k(:,i) - C_DT_CV * x_var_pre(:,i);
                S_k(:,:,i) = C_DT_CV * P_k_pre(:,:,i) * C_DT_CV' + R;
                K_k(:,:,i) = P_k_pre(:,:,i) * C_DT_CV' * inv(S_k(:,:,i));
                x_var(:,i) = x_var_pre(:,i) + K_k(:,:,i) * y_tilde_k(:,i);
                P_k(:,:,i) = (eye(N_states) - K_k(:,:,i) * C_DT_CV) * P_k_pre(:,:,i);
            end
    end

    x_hat = x_var;


%% Test Plots
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
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% plot(1:1:N_steps,y_tilde_k,'Linewidth',2)

end

%% 
function [Covar] = calcPopCovar(error_x, error_y)
    Covar = (error_x*error_y')/length(error_x);
end


%%
            %i
            %mod(i,100)
            % if mod(i,reset_steps) == 0
            %     P_k_pre(:,:,i) = P_k_pre(:,:,1);
            %     P_k(:,:,i)     = P_k(:,:,1);
            %     K_k(:,:,i)     = K_k(:,:,1);
            % end


