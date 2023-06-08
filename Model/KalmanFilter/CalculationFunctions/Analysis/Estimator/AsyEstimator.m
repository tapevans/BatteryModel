%% Asymptotic Estimator
%
%
%%
function [x_hat] = AsyEstimator(sys,FLAG,SIM,P,x_hat_0,u,z)
%% Perform Pre-Calcs
if FLAG.DoPreCalc
    [K_infty, ~] = AsymptoticPreCalcs(FLAG,SIM,sys);
else
    % Put something here, len(sys)
end

u_k = reshape(u,1,[]);
z_k = reshape(z,1,[]);


%% Initialize Variables
    A_DT = sys.A;
    B_DT = sys.B;
    C_DT    = sys.C;
    C_DT_CV = sys.C(P.cell_voltage,:);

    [N_measur, N_states] = size(C_DT_CV);
    [N_measur_all, ~]    = size(C_DT);

    K_infty = K_infty(:,P.cell_voltage);

    N_steps = length(u);
    
    x_asy     = zeros(N_states,N_steps);   % Estimator States
    x_asy_pre = zeros(N_states,N_steps);   % Estimator States predict phase
    y_tilde_k = zeros(N_measur,N_steps);   % Estimation Error


%% Calculate an Initial Error Covariance
    % sigma = 3e-1;
    % sigma = 1e-4;
    % sigma = 1e-5;
    sigma = 1e-7;
    y_des = sigma * ones(N_measur_all,1);
    
    % Try using inverse
    x_hat_0_test = C_DT \ y_des;

    % % Calculate P_k|k
    % error = zeros(size(x_hat_0_test)) - x_hat_0_test;
    % 
    % P_kk = nan(length(x_hat_0_test));
    % for i = 1:length(x_hat_0_test)
    %     for j = 1:length(x_hat_0_test)
    %         [P_kk(i,j)] = calcPopCovar(error(i,:), error(j,:)); %%%%%%May be able to change this so we don't need the calcPopCovar here
    %     end
    % end


%%
    % x_asy(:,1) = x_hat_0;
    x_asy(:,1) = x_hat_0_test;
        
    % if FLAG.UseWrongIC
    %     z_init = SIM.y_0_FOM_Offset(P.cell_voltage,1);
    % else
        z_init = SIM.y_0_FOM(P.cell_voltage,1);
    % end
    
    %z_init = SIM.y_0_FOM(P.cell_voltage,1);
    initial_z = C_DT_CV * x_hat_0 + z_init;


%% Run Estimation
    for i = 2:N_steps
        % Predict Phase
            x_asy_pre(:,i) = A_DT * x_asy(:,i-1) + B_DT * u_k(:,i-1);
        
    %     % Actuate System
    %         x_EST_sys(:,i) = A_DT * x_EST_sys(:,i-1)  +  B_DT * u_k(i-1)  +  w_k(:,i-1);
    %         z_k(:,i) = C * x_EST_sys(:,i)  +  v_k(:,i);
    
        % Update (Correction) Phase
            y_tilde_k(:,i) = z_k(:,i) - ( C_DT_CV * x_asy_pre(:,i) + z_init);
            x_asy(:,i) = x_asy_pre(:,i) + K_infty * y_tilde_k(:,i);
    end
    
    x_hat = x_asy;

end