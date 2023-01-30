%% Battery Kalman Testing
%
%
%
clear all; close all; clc;
% States of Interest
    % Cell Voltage
    % Delta Phi   @AN/SEP
    % Temperature @AN/SEP
    % C_Liion     @AN/SEP
    % X_surf      @AN/SEP
    % i_Far       @AN/SEP
    % eta         @AN/SEP

%% Parameters
step_filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanTesting\P2DSystem\KalmanTest_KPCont_KalmanTestStepSOC50.mat';
impl_filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanTesting\P2DSystem\KalmanTest_KPCont_KalmanTestDTImpulseTs1.0SOC50.mat';
sine_filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanTesting\P2DSystem\KalmanTest_EIS_SIN_w0.41888_SOC50.mat';

r  = 10; % Rank(states) of the ROM

% % Simulation inputs
%     % Modes
%     % 1) Step input
%     % 2) Sinusoidal
    FLAG.InputType = 2;
% %         SIM.V_s = 1;                 % Input amplitude
% %         SIM.fq  = 1/((SIM.R_1*SIM.C_1));
% % %         SIM.fq  = 1/20;
% %         SIM.Npoints_step = 60;
% %         SIM.Npoints_sine = 1000;

% What's being measured
        % 1) Measure Cell Voltage
        FLAG.Measure = 1;

% Estimator
% Noise
    Q_0 = 1e-6;
    R_0 = 1e-6;
    N.inputs = 1;
    N.measur = 1;


%% Code used to test different frequencies for sine wave
% Ns [=] samples/cycle
% Ts [=] samples/sec
% fq [=] cycles/sec
% Ns = 15; Ts = 1; fq = Ts/Ns; figure; plot((0:Ts:1/fq),sin((2*pi*fq)*(0:Ts:1/fq)),'-o')

%% Step Data
sim_data.step = load(step_filename);


%% Impulse Data
sim_data.impl = load(impl_filename);


%% Harmonic Data
sim_data.sine = load(sine_filename);


%% Create ROM
% Local Variable Pointers
    i = 1;
    localP.cell_voltage  = i; i = i + 1;
    localP.delta_phi     = i; i = i + 1;
    localP.C_Liion       = i; i = i + 1;
    localP.X_Li_surf     = i; i = i + 1;
    localP.i_Far         = i; i = i + 1;
    localP.eta           = i; i = i + 1;
    
    localP.numP = i-1;

% Extract g_k
    idx_vec = 1:sim_data.impl.SIM.TsMultiple:length(sim_data.impl.t_soln);

    i_user_DT = sim_data.impl.i_user(idx_vec);
    pulse_idx = find(i_user_DT ~= 0 );

    add = 0; % Additional offset of data

    g_k = zeros(localP.numP , length(idx_vec)-(pulse_idx-2+add)); %%%%%%%%% -2 is hardcoded
    IC  = zeros(localP.numP,1);
    
    % Cell Voltage
        cell_voltage_DT = sim_data.impl.cell_voltage(idx_vec);
        g_k_temp = (cell_voltage_DT(pulse_idx-1+add:end))';
        IC(localP.cell_voltage) = g_k_temp(1);
        g_k(localP.cell_voltage,:) = g_k_temp - g_k_temp(1);

    % Delta Phi
        del_phi_DT = sim_data.impl.del_phi(idx_vec,sim_data.impl.N.N_CV_AN);
        g_k_temp = (del_phi_DT(pulse_idx-1+add:end))';
        IC(localP.delta_phi) = g_k_temp(1);
        g_k(localP.delta_phi,:) = g_k_temp - g_k_temp(1);

    % C_Li^+
        C_Liion_DT = sim_data.impl.C_Liion(idx_vec,sim_data.impl.N.N_CV_AN);
        g_k_temp = (C_Liion_DT(pulse_idx-1+add:end))';
        IC(localP.C_Liion) = g_k_temp(1);
        g_k(localP.C_Liion,:) = g_k_temp - g_k_temp(1);

    % X_surf
        X_Li_surf_DT = sim_data.impl.X_Li_surf(idx_vec,sim_data.impl.N.N_CV_AN);
        g_k_temp = (X_Li_surf_DT(pulse_idx-1+add:end))';
        IC(localP.X_Li_surf) = g_k_temp(1);
        g_k(localP.X_Li_surf,:) = g_k_temp - g_k_temp(1);

    % i_Far
        i_Far_DT = sim_data.impl.i_Far(idx_vec);
        g_k_temp = (i_Far_DT(pulse_idx-1+add:end))';
        IC(localP.i_Far) = g_k_temp(1);
        g_k(localP.i_Far,:) = g_k_temp - g_k_temp(1);

    % eta
        eta_DT = sim_data.impl.eta(idx_vec);
        g_k_temp = (eta_DT(pulse_idx-1+add:end))';
        IC(localP.eta) = g_k_temp(1);
        g_k(localP.eta,:) = g_k_temp - g_k_temp(1);

    [A_r,B_r,C_r,D_r] = HoKalmanReduction(g_k,r);
    new_obsv_rank = rank(obsv(A_r,C_r))
    r

    
%% ROM Step Response from SS
    Ts = sim_data.impl.SIM.Ts; % [s]    
    SS_ROM = ss(A_r,B_r,C_r,D_r,Ts);
    t_final_step = sim_data.step.t_soln(end);
    [y_SS_ROM,t_SS_ROM,x_SS_ROM] = step(SS_ROM,t_final_step);
    
    C_r_meas = C_r(1,:);

%% get LC Parameters for battery
    [LC_param]     = getLCParam(sim_data.impl);
    [OutputMatrix,P.OM] = getOutputMatrix(sim_data.impl.N  ,sim_data.impl.P , sim_data.impl.AN);
    

%% Determine Estimator parameters
    Q = Q_0*eye(N.inputs); % Process Noise
    R = R_0*eye(N.measur); % Measurement Noise
    [P_infty  ,~,~] =  dare(A_r', C_r_meas', B_r*Q*B_r' ,R);
    K_infty = P_infty*C_r_meas'*inv(C_r_meas*P_infty*C_r_meas' + R);


%% Run Simulink Step with Noise
% Input Signal
InputSignal.t_final   = 100;
if FLAG.InputType == 1 % Step
    InputSignal.step_time = 2*Ts;
    InputSignal.step_amp  = 1;
else
    InputSignal.omega     = sim_data.sine.SIM.freq;
    InputSignal.sine_amp  = 1;
end


% Noise
    Noise.Q = Q_0*eye(N.inputs); % Process Noise
    Noise.R = R_0*eye(N.measur); % Measurement Noise
%     Noise.Ts_noise = 0.5;
    Noise.Ts_noise = 10;
    Noise.pow_Q = Noise.Ts_noise * Noise.Q;
    Noise.pow_R = Noise.Ts_noise * Noise.R;

% Simulink Parameters
    confidence = 0.5;
    P_k_0 = confidence*eye(r);
    Iden = eye(r);
%     x0_est = C_r\IC;
    x0_est = zeros(length(A_r),1);
    IC_est = C_r*x0_est;    

%     x0_ROM = C_r\IC;
    x0_ROM = zeros(length(A_r),1);

    Slink.C_plant = OutputMatrix(1,:);
    Slink.A_r = A_r;
    Slink.B_r = B_r;
    Slink.C_r_meas = C_r_meas;
    Slink.x0_ROM = x0_ROM;
    Slink.x0_est = x0_est;
    Slink.P_k_0 = P_k_0;
    Slink.Iden = Iden;
    Slink.K_infty = K_infty;

% Call Simulink
    disp('Running Sine')
    [out] = runSimulink(InputSignal,Noise,Slink,LC_param,sim_data.impl,FLAG);


%% Post Calcs
% Plant
    t_slink_plant = out.x_LC_Battery.time;
    if length(size(out.x_LC_Battery.signals.values))>2
        for i = 1:length(t_slink_plant)
            x_slink_plant(i,:) = reshape(out.x_LC_Battery.signals.values(:,:,i),1,[]);
        end
    else
        x_slink_plant = out.x_LC_Battery.signals.values;
    end
    y_slink_plant = (OutputMatrix*x_slink_plant')';

% Plant w/ Noise
    t_slink_plant_noise = out.x_LC_Battery_noise.time;
    if length(size(out.x_LC_Battery_noise.signals.values))>2
        for i = 1:length(t_slink_plant_noise)
            x_slink_plant_noise(i,:) = reshape(out.x_LC_Battery_noise.signals.values(:,:,i),1,[]);
        end
    else
        x_slink_plant_noise = out.x_LC_Battery_noise.signals.values;
    end
    y_slink_plant_noise = (OutputMatrix*x_slink_plant_noise')';

% Variable Estimator
    t_slink_est = out.x_hat.time;
    if length(size(out.x_hat.signals.values))>2
        for i = 1:length(t_slink_est)
            x_slink_hat(i,:) = reshape(out.x_hat.signals.values(:,:,i),1,[]);
        end
    else
        x_slink_hat = out.x_hat.signals.values;
    end
    y_slink_est = (C_r*x_slink_hat')';

% Asymtoptic Estimator
    t_slink_est_asy = out.x_hat_asy.time;
    if length(size(out.x_hat_asy.signals.values))>2
        for i = 1:length(t_slink_est_asy)
            x_slink_hat_asy(i,:) = reshape(out.x_hat_asy.signals.values(:,:,i),1,[]);
        end
    else
        x_slink_hat_asy = out.x_hat_asy.signals.values;
    end
    y_slink_est_asy = (C_r*x_slink_hat_asy'+IC)';

% ROM Model
    t_slink_ROM = out.x_ROM.time;
    if length(size(out.x_ROM.signals.values))>2
        for i = 1:length(t_slink_ROM)
            x_slink_ROM(i,:) = reshape(out.x_ROM.signals.values(:,:,i),1,[]);
        end
    else
        x_slink_ROM = out.x_ROM.signals.values;
    end
    y_slink_ROM = (C_r*x_slink_ROM' +IC)';
    
% P_infty for each desired output
P_infty_des = C_r*P_infty*C_r';


%% Save Results
save('ResultsSine.mat')


%% Plot Results
P2D_plots