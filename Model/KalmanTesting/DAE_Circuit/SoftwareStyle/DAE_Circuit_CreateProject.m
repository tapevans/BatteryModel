%% Create DAE Circuit Project

%% Filename
% Project Folder

%% State to Measure
%  1) V_1
%  2) V_2
%  3) V_3
%  4) All
    FLAG.C_mode = 4;

%     FLAG.NStates = 3; % Use the 5 or 3 state system
%% Noise
% Q Modes
%  1) Input Q
%  2) State Q
    FLAG.QMode = 1;

%% Simulation
    SIM.Ts = 1e-6; % Sampling Rate (When to save outputs)

%     SIM.t_final = 6e-6;
    SIM.t_final = 6e-4;
%     SIM.t_final = 5e-3;
    % SIM.t_final = 5*max(SIM.C_1*SIM.R_1,SIM.C_4*SIM.R_4);


    SIM.V_init        = 0;
    SIM.V_step_height = 1;


%% Call Inputs
[SIM,FLAG] = inputs(SIM,FLAG);


%% Make Save Filename
%     Q_str = num2str(Q_0);
%     R_str = num2str(R_0);
%     Ts_str = num2str(Ts);
%     C_str = num2str(C_mode);
%     switch FLAG.QMode
%         case 1 % Input Q
%             Noise_str = 'InputQ';
%         case 2 % State Q
%             Noise_str = 'StateQ';
%     end
%     
%     filename = ['Slink_Step_' Noise_str '_Q' Q_str '_R' R_str '_Ts' Ts_str '_C' C_str '.mat'];
%     filepath = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanTesting\DAE_Circuit\PlantOutputData\Results';
%     save_filename = [filepath filesep filename];




