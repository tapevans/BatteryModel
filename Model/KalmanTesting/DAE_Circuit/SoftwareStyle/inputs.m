%% Inputs
function [SIM,FLAG] = inputs(SIM,FLAG)
%% System Properties
    SIM.C_1 = 1e-6;
    SIM.C_4 = 1e-4;
    
    SIM.R_1 = 2;
    SIM.R_2 = 2;
    SIM.R_3 = 2;
    SIM.R_4 = 2;
    

%% Noise Properties
    SIM.tc  = 1e-7; % [s],Correlation Time (Noise Sampling Rate)
    SIM.Q_0 = 1e-6; % Process Noise
    SIM.R_0 = 1e-6; % Measurement Noise

end