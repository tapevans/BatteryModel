%% Asymptotic Estimator
%
%
%%
function [x_hat] = AsyEstimator(sys,FLAG,SIM,x_hat_0,u,z)
%% Perform Pre-Calcs
[K_infty, ~] = AsymptoticPreCalcs(FLAG,SIM,sys);

u_k = reshape(u.value,1,[]);
z_k = reshape(z.value,1,[]);


%% Initialize Variables
N_steps = length(u.time);
[N_measur, N_states] = size(sys.C);

x_asy     = zeros(N_states,N_steps);   % Estimator States
x_asy_pre = zeros(N_states,N_steps);   % Estimator States predict phase
y_tilde_k = zeros(N_measur,N_steps);   % Estimation Error

A_DT = sys.A;
B_DT = sys.B;
C_DT = sys.C;

x_asy(:,1) = x_hat_0;


%% Run Estimation
for i = 2:N_steps
    % Predict Phase
        x_asy_pre(:,i) = A_DT * x_asy(:,i-1) + B_DT * u_k(:,i-1);
    
%     % Actuate System
%         x_EST_sys(:,i) = A_DT * x_EST_sys(:,i-1)  +  B_DT * u_k(i-1)  +  w_k(:,i-1);
%         z_k(:,i) = C * x_EST_sys(:,i)  +  v_k(:,i);

    % Update (Correction) Phase
        y_tilde_k(:,i) = z_k(:,i) - C_DT * x_asy_pre(:,i);
        x_asy(:,i) = x_asy_pre(:,i) + K_infty * y_tilde_k(:,i);
end

x_hat = x_asy;

end