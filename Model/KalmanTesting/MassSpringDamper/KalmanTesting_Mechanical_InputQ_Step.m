%% KalmanTesting - Mechanical System - Disturbance on the Input
%
%
%
clear all; %close all; clc;
%% FLAGS
% Simulations
    FLAG.ODE         = 1;
    FLAG.SS_CT       = 1;
    FLAG.SS_DT       = 1;
    FLAG.SS_DT_Noise = 1;
    FLAG.Estimator   = 1;
        FLAG.RunAsymtotic = 1;
        FLAG.RunVariable  = 1;
    FLAG.SIMULINK         = 0;

%     FLAG.UseSavedDataForCodedPlant    = 0;
%     FLAG.UseSavedDataForSimulinkPlant = 0;

% Input Force
    % 1) Step
    % 2) Harmonic
    FLAG.InputMode = 1;
        SIM.F_in = 1; % Force Input (Amplitude)
        SIM.fq = 1; % [Hz]

% Plots
    FLAG.PLOT.x1Compare = 0;
    FLAG.PLOT.v1Compare = 0;
    FLAG.PLOT.x2Compare = 0;
    FLAG.PLOT.v2Compare = 0;
    FLAG.PLOT.ODE         = 0;
    FLAG.PLOT.SS_CT       = 0;
    FLAG.PLOT.SS_DT       = 0;
    FLAG.PLOT.SS_DT_Noise = 0;
    FLAG.PLOT.Estimator   = 1;
        FLAG.PLOT.RunAsymtotic = 1;
        FLAG.PLOT.RunVariable  = 0;
    FLAG.PLOT.SIMULINK.ODE         = 0;
    FLAG.PLOT.SIMULINK.SS_CT       = 0;
    FLAG.PLOT.SIMULINK.SS_DT       = 0;
    FLAG.PLOT.SIMULINK.SS_CT_Noise = 0;
    FLAG.PLOT.SIMULINK.Estimator   = 0;
        FLAG.PLOT.SIMULINK.RunAsymtotic = 0;
        FLAG.PLOT.SIMULINK.RunVariable  = 0;


%% System Parameters
SIM.k1 = .75;  % Spring
SIM.k2 = .75;  % Spring
SIM.b  = 10;  % Damper
SIM.m1 = 5; % Mass 1
SIM.m2 = 5; % Mass 2

% Noise
    Q_0 = 1e-5; % Process Noise
    R_0 = 1e1; % Measurement Noise

Ts = 5; %[s]

% t_final = 50;
% t_final = 100;
t_final = 300;


%% Pointers
i = 1;
P.x1 = i; i = i+1;
P.v1 = i; i = i+1;
P.x2 = i; i = i+1;
P.v2 = i; i = i+1;


%% Create System A and B
[A_CT,B_CT] = getSSMat(SIM);


%% Run Step for ODE
if FLAG.ODE
    %     [x1 v1 x2 v2]
    x_0 = [0 ,0 ,0 ,0 ]'; 
    [t_ODE,x_ODE] = ode45(@(t,x)myFun(t,x,SIM,P,FLAG,A_CT,B_CT),[0,t_final],x_0);
end


%% Test Observability
% Measure x1
C = [1 0 0 0];
Ob = obsv(A_CT,C);
r = rank(Ob);
disp(['Measuring x1, the rank is ' num2str(r)])

% % Measure v1
% C = [0 1 0 0];
% Ob = obsv(A_CT,C);
% r = rank(Ob);
% disp(['Measuring v1, the rank is ' num2str(r)])

% % Measure x2
% C = [0 0 1 0];
% Ob = obsv(A_CT,C);
% r = rank(Ob);
% disp(['Measuring x2, the rank is ' num2str(r)])

% % Measure v2
% C = [0 0 0 1];
% Ob = obsv(A_CT,C);
% r = rank(Ob);
% disp(['Measuring v2, the rank is ' num2str(r)])

% C = eye(4);

[N_meas , N_states]  = size(C);
N_inputs = 1;


%% Make SS CT
if FLAG.SS_CT
    [r,c] = size(C);
    D = zeros(r,1);
    SS_CT = ss(A_CT,B_CT,C,D);
end


%% Step SS CT
if FLAG.SS_CT
    [y_SS_CT,t_SS_CT,x_SS_CT] = step(SS_CT,t_final);
end


%% Make SS DT
% % Manually
% A_DT = expm(A_CT*Ts);
% B_DT = inv(A_CT)*(A_DT - eye(4))*B_CT; % This only works if A is invertable

% Using c2d
if FLAG.SS_DT
    SS_DT = c2d(SS_CT,Ts);
    A_DT = SS_DT.A;
    B_DT = SS_DT.B;
    C_DT = SS_DT.C; % Shouldn't change from CT -> DT, will change from FOM to ROM
end


%% Step SS DT
if FLAG.SS_DT
    [y_SS_DT,t_SS_DT,x_SS_DT] = step(SS_DT,t_final);
end


%% Try SS DT with Noise
if FLAG.SS_DT_Noise
    N_steps = ceil((t_final/Ts)+1);
    % Make input vector (not not since constant)
        u_k = SIM.F_in * ones(N_inputs,N_steps);

    % Make noise vectors
        mu = 0;

        % Process Noise
        Q = Q_0*eye(N_inputs);
%         Q = 1e-1*diag([1e-4 1e-6 1e-4 1e-6]);
%         Q = zeros(4);
        w_k = (chol(Q,'lower')*randn(N_inputs,1e6))'; % N_steps*100000
%         cov(w_k)
%         [mean(w_k(:,1)) mean(w_k(:,2)) mean(w_k(:,3)) mean(w_k(:,4))]
%         [(w_k(:,1))'*(w_k(:,1))    (w_k(:,2))'*(w_k(:,2))    (w_k(:,3))'*(w_k(:,3))    (w_k(:,4))'*(w_k(:,4))]

%         w_k = (sqrt(Q) * randn(4,1e6))'; % N_steps*100000
%         cov(w_k)
%         [mean(w_k(:,1)) mean(w_k(:,2)) mean(w_k(:,3)) mean(w_k(:,4))]
%         [(w_k(:,1))'*(w_k(:,1))    (w_k(:,2))'*(w_k(:,2))    (w_k(:,3))'*(w_k(:,3))    (w_k(:,4))'*(w_k(:,4))]

        % Measurement Noise
        R = R_0*eye(N_meas);
        v_k = (chol(R,'lower')*randn(N_meas,1e6))'; % N_steps*100000
%         cov(v_k)
%         [mean(v_k(:,1))] % mean(v_k(:,2)) mean(v_k(:,3)) mean(v_k(:,4))]
%         [(v_k(:,1))'*(v_k(:,1))]%    (v_k(:,2))'*(v_k(:,2))    (v_k(:,3))'*(w_k(:,3))    (w_k(:,4))'*(w_k(:,4))]

    % Create for loop
    w_k = w_k';
    v_k = v_k';
    t_SS_DT_n = (0:Ts:t_final)';
    x_SS_DT_n = zeros(N_states,N_steps);

    for i = 1:N_steps-1
        x_SS_DT_n(:,i+1) = A_DT*x_SS_DT_n(:,i) + B_DT*( u_k(:,i) + w_k(:,i) );
    end
    x_SS_DT_n = x_SS_DT_n';

end


%%
if FLAG.Estimator
    % Set Noise
        %Q = 1e-1*diag([1e-4 1e-6 1e-4 1e-6]);
%         Q = 1-3;
%         w_k = (chol(Q,'lower')*randn(1,1e6))';
%         w_k = w_k';

%         R = 1e-5;
%         v_k = (chol(R,'lower')*randn(1,1e6))';
%         v_k = v_k';

    % Set input vector
        u_k = SIM.F_in * ones(N_inputs,N_steps);

    % Set time vector
        t_asy = (0:Ts:t_final)';
        t_var = (0:Ts:t_final)';

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
        P_k_pre(:,:,1) = confidence*eye(N_states);
    
    % Set Initial Estimator State Estimates
        x_hat_0 = [0,0,0,0]';
%         x_hat_0 = [0.1,0.01,0.5,0.01]';
        x_asy(:,1) = x_hat_0;
        x_var(:,1) = x_hat_0;

    % Run Variable Estimator
    if FLAG.RunVariable
        for i = 2:N_steps
            % Predict Phase
                x_var_pre(:,i) = A_DT * x_var(:,i-1)  +  B_DT * u_k(:,i-1);
                P_k_pre(:,:,i) = A_DT * P_k(:,:,i-1) * A_DT'  +  B_DT*Q*B_DT';
            
            % Actuate System
                x_EST_sys(:,i) = A_DT * x_EST_sys(:,i-1)  +  B_DT * (u_k(:,i-1)  +  w_k(:,i-1));
                z_k(:,i) = C * x_EST_sys(:,i)  +  v_k(:,i);

            % Update (Correction) Phase
                y_tilde_k(:,i) = z_k(:,i) - C * x_var_pre(:,i);
                S_k(:,:,i) = C * P_k_pre(:,:,i) * C' + R;
                K_k(:,:,i) = P_k_pre(:,:,i) * C' * inv(S_k(:,:,i));
                x_var(:,i) = x_var_pre(:,i) + K_k(:,:,i) * y_tilde_k(:,i);
                P_k(:,:,i) = (eye(N_states) - K_k(:,:,i) * C) * P_k_pre(:,:,i);
        end
    end
    x_var = x_var';
    
    % Run Asymptotic Estimator
    if FLAG.RunAsymtotic
            % Calculate K_infty and P_infty
        %         [P_infty_idare ,~,~] = idare(A_DT , C , Q,R,0,eye(4));
        %         [P_infty_idaret,~,~] = idare(A_DT', C', Q,R,0,eye(4));
%                 [P_infty  ,~,~] =  dare(A_DT', C', B_DT*Q*B_DT' ,R);
                [P_infty  ,~,~] =  idare(A_DT', C', B_DT*Q*B_DT' ,R);
                K_infty = P_infty*C'*inv(C*P_infty*C' + R);

        for i = 2:N_steps
            % Predict Phase
                x_asy_pre(:,i) = A_DT * x_asy(:,i-1) + B_DT * u_k(:,i-1);
            
            % Actuate System
                x_EST_sys(:,i) = A_DT * x_EST_sys(:,i-1)  +  B_DT * (u_k(i-1)  +  w_k(:,i-1));
                z_k(:,i) = C * x_EST_sys(:,i)  +  v_k(:,i);

            % Update (Correction) Phase
                y_tilde_k(:,i) = z_k(:,i) - C * x_asy_pre(:,i);
                x_asy(:,i) = x_asy_pre(:,i) + K_infty * y_tilde_k(:,i);
        end
    end
    x_asy = x_asy';
    x_EST_sys = x_EST_sys';
    
    % kalman function output
        [~,K_fnc,P_fnc] = kalman(SS_DT,Q,R);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some Post Calcs
if FLAG.Estimator && FLAG.RunAsymtotic && FLAG.RunVariable
    P_k_last = P_k(:,:,end);
    P_k_pre_last = P_k_pre(:,:,end);


    disp('Real Time calculated P_k converged to ')
    disp(num2str(P_k_pre(:,:,end)))
    
    disp(newline)
    disp('P_\infty is calculated to be ')
    disp(num2str(P_infty(:,:)) )
    
    error_P_k = P_k_pre - P_infty;
    for i = 1:N_steps
        norm_error_P_k(:,i) = norm(error_P_k(:,:,i));
    end

    disp(newline)
    disp('Real Time Kalman gain K_k converged to ')
    disp( num2str(K_k(:,:,end)) )
    disp(newline)
    disp('Asymptotic Kalman gain is calculated to be ')
    disp(num2str(K_infty(:,:)))

    error_K_k = K_k - K_infty;
    for i = 1:N_steps
        norm_error_K_k(:,i) = norm(error_K_k(:,:,i));
    end

%     % Error for first NS steps
%     NS = 200;
%         figure
%         title('Convergence to Asymptotic')
%         yyaxis left
%         plot(1:1:N_steps,norm_error_P_k,'Linewidth',2)
%         ylabel('|Error P_k|')
%     
%         yyaxis right
%         plot(t_var,norm_error_K_k,'Linewidth',2)
%         ylabel('|Error K_k|')
%     
%         xlabel('Number of Discrete Steps')
%         xlim([1,NS])
% %         exportgraphics(gca,'PandKConverge.png','Resolution',1000)
%     
%     % Error for after 40 steps
%         figure
%         yyaxis left
%         plot(1:1:N_steps,norm_error_P_k,'Linewidth',2)
%         ylabel('|Error P_k|')
%     
%         yyaxis right
%         plot(t_var,norm_error_K_k,'Linewidth',2)
%         ylabel('|Error K_k|')
%         
%         xlabel('Number of Discrete Steps')
%         xlim([NS,Inf])

    % Covariance
        steady_state_step = 250;

        % Error Calc (Deviation)
        for i = 1:N_states
            error_var(:,i) = x_var(steady_state_step:end,i) - x_EST_sys(steady_state_step:end,i);
            error_asy(:,i) = x_asy(steady_state_step:end,i) - x_EST_sys(steady_state_step:end,i);
        end
        error_var = error_var';
        error_asy = error_asy';
        
        % Covar calc
        covar_var = nan(N_states);
        covar_asy = nan(N_states);
        for i = 1:N_states
            for j = 1:N_states
                [covar_var(i,j)] = calcPopCovar(error_var(i,:), error_var(j,:));
                [covar_asy(i,j)] = calcPopCovar(error_asy(i,:), error_asy(j,:));
            end
        end

        % Compare with CPC^T
        C_tilde = eye(N_states);

        for i = 1:N_states
            CPCT_calc(:,i) = C_tilde(i,:) * P_infty * C_tilde(i,:)' ;%+ R(i,i);
        end
        %CPCT_calc_var(:,i)

        disp(newline)
        disp('CPC^T calc is')
        disp(num2str(CPCT_calc))
        
        disp(newline)
        disp('Cov(x_variable) calc is')
        disp(num2str(covar_var))
        disp('With difference of')
%         covar_var-CPCT_calc
        disp(num2str(covar_var-P_infty))
        disp(newline)
        disp('Ratio of error')
        disp(num2str(abs((covar_var-P_infty)./P_infty)))


        disp(newline)
        disp('Cov(x_asymptotic) calc is')
        disp(num2str(covar_asy))
        disp('With difference of')
%         covar_asy-CPCT_calc))
        disp(num2str(covar_asy-P_infty))
        disp(newline)
        disp('Ratio of error')
        disp(num2str(abs((covar_asy-P_infty)./P_infty)))


%% SVD Analysis
% Observability Measurable
%     obsv_r_meas = getObservability(A_DT,C_DT(1,:));
    obsv_r_meas = obsv(A_DT,C_DT(1,:));
    [U_Orm,S_Orm,V_Orm] = svd(obsv_r_meas);

% C_tilde (This normally would be of the ROM)
    C_r = eye(N_states);

% Loop through C_tilde
    sing_val = nan(1,N_states);
    for mm = 1:N_states
        C_tilde = C_r(mm,:);

        % SVD of \tilde{C}
        [~ , ~ , V_tilde] = svd(C_tilde);

        % SVD of Augmented Observability
        obsv_aug = obsv_r_meas*V_tilde(:,1); %!!!!!! Do I need to check that V_tilde is in the column space if it is only 1 row?
        [~ , S_hat , ~] = svd(obsv_aug);

        sing_val(mm) = S_hat(1,1);
    end
    sing_val_norm = sing_val / S_Orm(1,1);



















end
















%% Simulink
if FLAG.SIMULINK
    k1 = SIM.k1;  % Spring
    k2 = SIM.k2;  % Spring
    b  = SIM.b;  % Damper
    m1 = SIM.m1; % Mass 1
    m2 = SIM.m2; % Mass 2
    F_in = SIM.F_in;
    Q = Q_0*eye(N_inputs); % Process Noise
    R = R_0*eye(N_meas);   % Measurement Noise
    Iden = eye(N_states);
    confidence = 0.5;
    P_k_0 = confidence*eye(N_states);
    x_0_est = [0;0;0;0];

%     N_samples_per_Ts_sample = 3;
%     delta = Ts/N_samples_per_Ts_sample;
%     trigger_t = 0:delta:t_final;
%     idx = (1 + N_samples_per_Ts_sample):N_samples_per_Ts_sample:length(trigger_t);
%     trigger_pulse = zeros(size(trigger_t));
%     trigger_pulse(idx) = 1;
%     TriggerSignal = [trigger_t' trigger_pulse'];

    [P_infty  ,~,~] =  dare(A_DT', C', B_DT*Q*B_DT' ,R);
    K_infty = P_infty*C'*inv(C*P_infty*C' + R);

    mdl = 'MSD_Overall';
    load_system(mdl)
    in = Simulink.SimulationInput(mdl);
    in = in.setModelParameter('StartTime','0','StopTime',num2str(t_final));
    
    mdlWks = get_param(in,'ModelWorkspace');
    assignin(mdlWks,'x0'   ,x_0)
    assignin(mdlWks,'k1'   ,SIM.k1)
    assignin(mdlWks,'k2'   ,SIM.k2)
    assignin(mdlWks,'b'    ,SIM.b)
    assignin(mdlWks,'m1'   ,SIM.m1)
    assignin(mdlWks,'m2'   ,SIM.m2)
    assignin(mdlWks,'F_in' ,SIM.F_in)
    assignin(mdlWks,'A_CT' ,A_CT)
    assignin(mdlWks,'B_CT' ,B_CT)
    assignin(mdlWks,'C'    ,C)
    assignin(mdlWks,'Ts'   ,Ts)
    assignin(mdlWks,'Q'    ,Q)
    assignin(mdlWks,'R'    ,R)
    assignin(mdlWks,'Iden' ,Iden)
    assignin(mdlWks,'P_k_0',P_k_0)
    assignin(mdlWks,'x_0_est',x_0_est)
%     assignin(mdlWks,'P_infty',P_infty)
    assignin(mdlWks,'K_infty',K_infty)
    
    out = sim(in);

    save_system
    close_system

%% Post Simulink
x_ODE_Slink = out.x_ODE.signals.values;
t_ODE_Slink = out.x_ODE.time;

for i = 1:length(out.x_CT_SS.time)
    x_SS_CT_Slink(i,:) = (out.x_CT_SS.signals.values(:,:,i))';
end
t_SS_CT_Slink = out.x_CT_SS.time;

for i = 1:length(out.x_DT_SS.time)
    x_SS_DT_Slink(i,:) = (out.x_DT_SS.signals.values(:,:,i))';
end
t_SS_DT_Slink = out.x_DT_SS.time;

for i = 1:length(out.x_CT_SS_N.time)
    x_SS_CT_N_Slink(i,:) = (out.x_CT_SS_N.signals.values(:,:,i))';
end
t_SS_CT_N_Slink = out.x_CT_SS_N.time;

% Calculate mean and covar of noise
w_k_Slink = reshape(out.w_k.signals.values,1,[]);
[mu_x] = calcMean(w_k_Slink);
[error_x] = calcError(w_k_Slink,mu_x);
[Covar_xx] = calcPopCovar(error_x, error_x);
% Q_0

v_k_Slink = reshape(out.v_k.signals.values,1,[]);
[mu_x] = calcMean(v_k_Slink);
[error_x] = calcError(v_k_Slink,mu_x);
[Covar_xx] = calcPopCovar(error_x, error_x);
% R_0

% P_infty
    out.P_k_k1.signals.values(:,:,end);
%     P_infty
    out.P_k_k.signals.values(:,:,end);

    out.K_k.signals.values(:,:,end);
%     K_infty


for i = 1:length(out.x_hat_insideEST.time)
    x_hat_Slink(i,:) = (out.x_hat_insideEST.signals.values(:,:,i))';
end
t_hat_Slink = out.x_hat_insideEST.time;

% Asymptotic
for i = 1:length(out.x_hat_insideEST_asy.time)
    x_hat_asy_Slink(i,:) = (out.x_hat_insideEST_asy.signals.values(:,:,i))';
end
t_hat_asy_Slink = out.x_hat_insideEST_asy.time;

%% CPCT+R Calc
time = 149;
idx_hat   = find(t_hat_Slink == time);
idx_plant = find(t_SS_CT_N_Slink == time);
error = x_hat_Slink(idx_hat+1:end,:) - x_SS_CT_N_Slink(idx_plant:end-1,:); 
error = error';

covar_var_Slink = nan(N_states);
for i = 1:N_states
    for j = 1:N_states
        [covar_var_Slink(i,j)] = calcPopCovar(error(i,:), error(j,:));
%         [covar_asy(i,j)] = calcPopCovar(error_asy(i,:), error_asy(j,:));
    end
end


disp('Compare simulink results')
disp('P_calc')
disp( num2str(covar_var_Slink) )

disp(newline)
disp('P_infty')
disp( num2str(P_infty) )

disp(newline)
disp('Error between them')
disp( num2str(covar_var_Slink-P_infty) )

disp(newline)
disp('Ratio of error')
disp( num2str(abs((covar_var_Slink-P_infty)/P_infty) ))

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting
% Position 1
    if FLAG.PLOT.x1Compare
        figure
        hold on
        title('Step Response Position 1')
        xlabel('Time (s)')
        ylabel('Position')
    
        if FLAG.PLOT.SIMULINK.ODE
            plot(t_ODE_Slink,x_ODE_Slink(:,P.x1),'o','LineWidth',2,'DisplayName','ODE Slink')
        end
        if FLAG.PLOT.SIMULINK.SS_CT
            plot(t_SS_CT_Slink,x_SS_CT_Slink(:,P.x1),'o','LineWidth',2,'DisplayName','SS CT Slink')
        end
        if FLAG.PLOT.SS_DT
            plot(t_SS_DT,x_SS_DT(:,P.x1),'o','LineWidth',5,'DisplayName','SS DT')
        end
        if FLAG.PLOT.SS_CT
            plot(t_SS_CT,x_SS_CT(:,P.x1),'LineWidth',2,'DisplayName','SS CT')
        end
        if FLAG.PLOT.ODE
            plot(t_ODE,x_ODE(:,P.x1),'LineWidth',2,'DisplayName','ODE')
        end
        if FLAG.PLOT.SIMULINK.SS_DT
            plot(t_SS_DT_Slink,x_SS_DT_Slink(:,P.x1),'o','LineWidth',2,'DisplayName','SS DT Slink')
        end
        if FLAG.PLOT.SIMULINK.SS_CT_Noise
            plot(t_SS_CT_N_Slink,x_SS_CT_N_Slink(:,P.x1),'o','LineWidth',2,'DisplayName','SS CT Noise Slink')
        end
        if FLAG.PLOT.SS_DT_Noise
            plot(t_SS_DT_n,x_SS_DT_n(:,P.x1),'o','LineWidth',2,'DisplayName','DT with Noise')
        end
        if FLAG.PLOT.Estimator
            if FLAG.PLOT.RunAsymtotic
                plot(t_asy,x_asy(:,P.x1),'o','LineWidth',5,'DisplayName','Asymp')
            end
            if FLAG.PLOT.RunVariable
                plot(t_var,x_var(:,P.x1),'o','LineWidth',2,'DisplayName','Variable K_k')
            end
            plot(t_asy,x_EST_sys(:,P.x1),'LineWidth',2,'DisplayName','Plant')
        end
        if FLAG.PLOT.SIMULINK.Estimator
            if FLAG.PLOT.SIMULINK.RunAsymtotic
                plot(t_hat_asy_Slink,x_hat_asy_Slink(:,P.x1),'o','LineWidth',5,'DisplayName','Asymp Slink')
            end
            if FLAG.PLOT.SIMULINK.RunVariable
                plot(t_hat_Slink,x_hat_Slink(:,P.x1),'o','LineWidth',2,'DisplayName','Variable k Slink')
            end
            plot(t_asy,x_EST_sys(:,P.x1),'LineWidth',2,'DisplayName','Plant Slink')
        end
        lgn = legend;
        lgn.Location = 'southeast';
        xlim([0,t_final])
    end

% Velocity 1
    if FLAG.PLOT.v1Compare
        figure
        hold on
        title('Step Response Velocity 1')
        xlabel('Time (s)')
        ylabel('Velocity')
    
        if FLAG.PLOT.SIMULINK.ODE
            plot(t_ODE_Slink,x_ODE_Slink(:,P.v1),'o','LineWidth',2,'DisplayName','ODE Slink')
        end
        if FLAG.PLOT.ODE
            plot(t_ODE,x_ODE(:,P.v1),'LineWidth',2,'DisplayName','ODE')
        end
        if FLAG.PLOT.SIMULINK.SS_CT
            plot(t_SS_CT_Slink,x_SS_CT_Slink(:,P.v1),'o','LineWidth',2,'DisplayName','SS CT Slink')
        end
        if FLAG.PLOT.SS_CT
            plot(t_SS_CT,x_SS_CT(:,P.v1),'LineWidth',2,'DisplayName','SS CT')
        end
        if FLAG.PLOT.SS_DT
            plot(t_SS_DT,x_SS_DT(:,P.v1),'o','LineWidth',5,'DisplayName','SS DT')
        end
        if FLAG.PLOT.SIMULINK.SS_DT
            plot(t_SS_DT_Slink,x_SS_DT_Slink(:,P.v1),'o','LineWidth',2,'DisplayName','SS DT Slink')
        end
        if FLAG.PLOT.SIMULINK.SS_CT_Noise
            plot(t_SS_CT_N_Slink,x_SS_CT_N_Slink(:,P.v1),'o','LineWidth',2,'DisplayName','SS CT Noise Slink')
        end
        if FLAG.PLOT.SS_DT_Noise
            plot(t_SS_DT_n,x_SS_DT_n(:,P.v1),'o','LineWidth',2,'DisplayName','DT with Noise')
        end
        if FLAG.PLOT.Estimator
            if FLAG.PLOT.RunAsymtotic
                plot(t_asy,x_asy(:,P.v1),'o','LineWidth',5,'DisplayName','Asymp')
            end
            if FLAG.PLOT.RunVariable
                plot(t_var,x_var(:,P.v1),'o','LineWidth',2,'DisplayName','Variable k')
            end
            plot(t_asy,x_EST_sys(:,P.v1),'LineWidth',2,'DisplayName','Plant')
        end
        if FLAG.PLOT.SIMULINK.Estimator
            if FLAG.PLOT.SIMULINK.RunAsymtotic
                plot(t_hat_asy_Slink,x_hat_asy_Slink(:,P.v1),'o','LineWidth',5,'DisplayName','Asymp Slink')
            end
            if FLAG.PLOT.SIMULINK.RunVariable
                plot(t_hat_Slink,x_hat_Slink(:,P.v1),'o','LineWidth',2,'DisplayName','Variable k Slink')
            end
            plot(t_asy,x_EST_sys(:,P.v1),'LineWidth',2,'DisplayName','Plant Slink')
        end
        lgn = legend;
        xlim([0,t_final])
    end

% Position 2
    if FLAG.PLOT.x2Compare
        figure
        hold on
        title('Step Response Position 2')
        xlabel('Time (s)')
        ylabel('Position')
    
        if FLAG.PLOT.SIMULINK.ODE
            plot(t_ODE_Slink,x_ODE_Slink(:,P.x2),'o','LineWidth',2,'DisplayName','ODE Slink')
        end
        if FLAG.PLOT.ODE
            plot(t_ODE,x_ODE(:,P.x2),'LineWidth',2,'DisplayName','ODE')
        end
        if FLAG.PLOT.SIMULINK.SS_CT
            plot(t_SS_CT_Slink,x_SS_CT_Slink(:,P.x2),'o','LineWidth',2,'DisplayName','SS CT Slink')
        end
        if FLAG.PLOT.SS_CT
            plot(t_SS_CT,x_SS_CT(:,P.x2),'LineWidth',2,'DisplayName','SS CT')
        end
        if FLAG.PLOT.SS_DT
            plot(t_SS_DT,x_SS_DT(:,P.x2),'o','LineWidth',5,'DisplayName','SS DT')
        end
        if FLAG.PLOT.SIMULINK.SS_DT
            plot(t_SS_DT_Slink,x_SS_DT_Slink(:,P.x2),'o','LineWidth',2,'DisplayName','SS DT Slink')
        end
        if FLAG.PLOT.SIMULINK.SS_CT_Noise
            plot(t_SS_CT_N_Slink,x_SS_CT_N_Slink(:,P.x2),'o','LineWidth',2,'DisplayName','SS CT Noise Slink')
        end
        if FLAG.PLOT.SS_DT_Noise
            plot(t_SS_DT_n,x_SS_DT_n(:,P.x2),'o','LineWidth',2,'DisplayName','DT with Noise')
        end
        if FLAG.PLOT.Estimator
            if FLAG.PLOT.RunAsymtotic
                plot(t_asy,x_asy(:,P.x2),'o','LineWidth',5,'DisplayName','Asymp')
            end
            if FLAG.PLOT.RunVariable
                plot(t_var,x_var(:,P.x2),'o','LineWidth',2,'DisplayName','Variable k')
            end
            plot(t_asy,x_EST_sys(:,P.x2),'LineWidth',2,'DisplayName','Plant')
        end
        if FLAG.PLOT.SIMULINK.Estimator
            if FLAG.PLOT.SIMULINK.RunAsymtotic
                plot(t_hat_asy_Slink,x_hat_asy_Slink(:,P.x2),'o','LineWidth',5,'DisplayName','Asymp Slink')
            end
            if FLAG.PLOT.SIMULINK.RunVariable
                plot(t_hat_Slink,x_hat_Slink(:,P.x2),'o','LineWidth',2,'DisplayName','Variable k Slink')
            end
            plot(t_asy,x_EST_sys(:,P.x2),'LineWidth',2,'DisplayName','Plant Slink')
        end
        lgn = legend;
        lgn.Location = 'southeast';
        xlim([0,t_final])
    end

% Velocity 2
    if FLAG.PLOT.v2Compare
        figure
        hold on
        title('Step Response Velocity 2')
        xlabel('Time (s)')
        ylabel('Velocity')
    
        if FLAG.PLOT.SIMULINK.ODE
            plot(t_ODE_Slink,x_ODE_Slink(:,P.v2),'o','LineWidth',2,'DisplayName','ODE Slink')
        end
        if FLAG.PLOT.ODE
            plot(t_ODE,x_ODE(:,P.v2),'LineWidth',2,'DisplayName','ODE')
        end
        if FLAG.PLOT.SIMULINK.SS_CT
            plot(t_SS_CT_Slink,x_SS_CT_Slink(:,P.v2),'o','LineWidth',2,'DisplayName','SS CT Slink')
        end
        if FLAG.PLOT.SS_CT
            plot(t_SS_CT,x_SS_CT(:,P.v2),'LineWidth',2,'DisplayName','SS CT')
        end
        if FLAG.PLOT.SS_DT
            plot(t_SS_DT,x_SS_DT(:,P.v2),'o','LineWidth',5,'DisplayName','SS DT')
        end
        if FLAG.PLOT.SIMULINK.SS_DT
            plot(t_SS_DT_Slink,x_SS_DT_Slink(:,P.v2),'o','LineWidth',2,'DisplayName','SS DT Slink')
        end
        if FLAG.PLOT.SIMULINK.SS_CT_Noise
            plot(t_SS_CT_N_Slink,x_SS_CT_N_Slink(:,P.v2),'o','LineWidth',2,'DisplayName','SS CT Noise Slink')
        end
        if FLAG.PLOT.SS_DT_Noise
            plot(t_SS_DT_n,x_SS_DT_n(:,P.v2),'o','LineWidth',2,'DisplayName','DT with Noise')
        end
        if FLAG.PLOT.Estimator
            if FLAG.PLOT.RunAsymtotic
                plot(t_asy,x_asy(:,P.v2),'o','LineWidth',5,'DisplayName','Asymp')
            end
            if FLAG.PLOT.RunVariable
                plot(t_var,x_var(:,P.v2),'o','LineWidth',2,'DisplayName','Variable k')
            end
            plot(t_asy,x_EST_sys(:,P.v2),'LineWidth',2,'DisplayName','Plant')
        end
        if FLAG.PLOT.SIMULINK.Estimator
            if FLAG.PLOT.SIMULINK.RunAsymtotic
                plot(t_hat_asy_Slink,x_hat_asy_Slink(:,P.v2),'o','LineWidth',5,'DisplayName','Asymp Slink')
            end
            if FLAG.PLOT.SIMULINK.RunVariable
                plot(t_hat_Slink,x_hat_Slink(:,P.v2),'o','LineWidth',2,'DisplayName','Variable k Slink')
            end
            plot(t_asy,x_EST_sys(:,P.v2),'LineWidth',2,'DisplayName','Plant Slink')
        end
        lgn = legend;
        xlim([0,t_final])
    end

%% Arrange Figures
FigArrange = 0;
if FigArrange == 1
    fig = gcf;
    NumFig = fig.Number;

    Ncol = 3;
    
    for i = 1:NumFig
        f = figure(i);
        k = mod(i-1,Ncol);
        row = mod(fix((i-1)/Ncol),2);
        if row == 0
            r = 575;
%             r = 540;
        elseif row == 1
            r = 62;
        end
        f.Position = [k*575+15 r 560 420];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
%% ODE Function
function [x_dot] = myFun(t,x,SIM,P,FLAG,A,B)
% if FLAG.DOSinInput
%     F = SIM.F_in*sin(SIM.omega*t);
% else
%     F = SIM.F_in;
% end

% Input Force
u = SIM.F_in;

% Derivative 
x_dot = A*x + B*u;


end

%% A and B SS CT
% function [A,B,C,D] = getSSMat(SIM)
function [A,B] = getSSMat(SIM)
k1 = SIM.k1; % Spring
k2 = SIM.k2; % Spring
b  = SIM.b;  % Damper
m1 = SIM.m1; % Mass 1
m2 = SIM.m2; % Mass 2

A = [ 0                 1     0     0
     -(k1/m1 + k2/m1)  -b/m1  k2/m1 0
      0                 0     0     1
      k2/m2             0    -k2/m2 0];

B = [0
     1/m1
     0
     0];
% 
% C = [1 0 0 0
%      0 0 1 0];
% 
% D = [0
%      0];

end
function [mu] = calcMean(x)
    mu = mean(x);
end
function [error] = calcError(x,mu)
% difference between x and the expected value (mu)
    error = x-mu;
end
function [Covar] = calcPopCovar(error_x, error_y)
    Covar = (error_x*error_y')/length(error_x);
end
function [Covar] = calcSampleCovar(error_x, error_y)
    Covar = (error_x*error_y')/(length(error_x)-1);
end
