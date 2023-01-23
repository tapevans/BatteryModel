%% Test Discrete Time in Simulink
%
%
clear all; close all; clc;

%% 
% Input current will be a square wave
% (Maybe change to interpolate instead of hold after) <-- Did this
u_pulses = [0,1,0,0,-1,0,0,1,0,-1,0,0,1,0,0,-1,0,0,1,0,-1,0,0,1,0,0,-1,0,0,1,0,-1,0,0,1,0,0,-1,0,0,1,0,-1,0];
Ts_pulse = 1;
Ts_sample = 1;

N_steps = length(u_pulses) - 1;
t_vec = 0:1:N_steps;
t_vec_a = t_vec*Ts_pulse;

u_t = [t_vec_a' , u_pulses'];
t_final = t_vec_a(end);

multiple = 100;
t_vec_fine = 0:1/multiple:N_steps;
t_vec_fine = t_vec_fine*Ts_pulse;

Q_0 = 1e-6;
R_0 = 1e-4;

%%
A = 1e-3*eye(2);
% A = [1e-3 0
%      0    1e-4];
B = 1*ones(2,1);
% C = 1*ones(1,2);
% C = [1 0];
C = [1 0
     0 1];
D = 0;

[N_meas , N_states]  = size(C);
N_inputs = 1;

sys_CT = ss(A,B,C,D);

x_0 = zeros(N_states,1);

[t_soln,y_soln] = ode45(@(t,x)odefun(t,x,A,B,u_t),t_vec_fine,x_0);

% figure
% hold on
% plot(t_soln,y_soln)
% plot(t_soln,interp1(u_t(:,1),u_t(:,2),t_soln))

%% DT
SS_DT = c2d(sys_CT,Ts_sample);
A_DT = SS_DT.A;
B_DT = SS_DT.B;
% A_DT = expm(A*Ts);
% B_DT = inv(A)*(A_DT - eye(1))*B; % This only works if A is invertable


%% Testing Just P calc in Simulink
Q = Q_0*eye(N_inputs);
w_k = (chol(Q,'lower')*randn(N_inputs,1e6))'; % N_steps*100000

R = R_0*eye(N_meas);
v_k = (chol(R,'lower')*randn(N_meas,1e6))'; % N_steps*100000


    w_k = w_k';
    v_k = v_k';

% [P_infty  ,~,~] =  dare(A_DT', C', B_DT*Q*B_DT' ,R);
[P_infty,~,~] = idare(A_DT',C',B_DT*Q*B_DT',R,0,eye(2));

% Initialize Solution Variables
    x_EST_sys = zeros(N_states,N_steps);   % Real System (Plant) States
    x_asy     = zeros(N_states,N_steps);   % Estimator States
    x_asy_pre = zeros(N_states,N_steps);   % Estimator States predict phase
    x_var     = zeros(N_states,N_steps);   % Estimator States
    x_var_pre = zeros(N_states,N_steps);   % Estimator States predict phase
    z_k       = zeros(N_meas,N_steps);   % Measurement from Real System
    y_tilde_k = zeros(N_meas,N_steps);   % Estimation Error
    S_k       = zeros(N_meas,N_meas,N_steps); % Error Covariance
    K_k       = zeros(N_states,N_meas,N_steps); % Kalman Gains
    P_k       = zeros(N_states,N_states,N_steps); % Error Covariance
    P_k_pre   = zeros(N_states,N_states,N_steps); % Error Covariance predict phase
            
% Set Initial Error Covariance (Identity if confident in model, 0 if not)
    confidence = 0.5; % Range from [0,1]
    P_k(:,:,1) = confidence*eye(N_states);
%     P_k_pre(:,:,1) = confidence*eye(N_states);

% Set Initial Estimator State Estimates
%         x_hat_0 = [0,0,0,0]';
%     x_hat_0 = [0.1,0.01,0.5,0.01]';
    
%     x_asy(:,1) = x_hat_0;
%     x_var(:,1) = x_hat_0;

for i = 2:N_steps
    % Predict Phase
        x_var_pre(:,i) = A_DT * x_var(:,i-1)  +  B_DT * u_pulses(:,i-1);
        P_k_pre(:,:,i) = A_DT * P_k(:,:,i-1) * A_DT'  +  B_DT*Q*B_DT';
    
    % Actuate System
        x_EST_sys(:,i) = A_DT * x_EST_sys(:,i-1)  +  B_DT * (u_pulses(:,i-1)  +  w_k(:,i-1));
        z_k(:,i) = C * x_EST_sys(:,i)  +  v_k(:,i);

    % Update (Correction) Phase
        y_tilde_k(:,i) = z_k(:,i) - C * x_var_pre(:,i);
        S_k(:,:,i) = C * P_k_pre(:,:,i) * C' + R;
        K_k(:,:,i) = P_k_pre(:,:,i) * C' * inv(S_k(:,:,i));
        x_var(:,i) = x_var_pre(:,i) + K_k(:,:,i) * y_tilde_k(:,i);
        P_k(:,:,i) = (eye(N_states) - K_k(:,:,i) * C) * P_k_pre(:,:,i);
end

%% Creating Trigger Signal
% Want a pulse at every 1 second but I need two zeros before each pulse
% So I need a sample at every 1/3 sec
delta = 0.5;
trigger_t = 0:delta:10;
idx = 3:2:length(trigger_t);
trigger_pulse = zeros(size(trigger_t));
trigger_pulse(idx) = 1;
TriggerSignal = [trigger_t' trigger_pulse'];


% delta = 1/3;
% trigger_t = 0:delta:10;
% idx = 4:3:length(trigger_t);
% trigger_pulse = zeros(size(trigger_t));
% trigger_pulse(idx) = 1;
% TriggerSignal = [trigger_t' trigger_pulse'];


%%
function [x_dot] = odefun(t,x,A,B,u)
    u_k = interp1(u(:,1),u(:,2),t);
    x_dot = A*x + B*u_k;
end