%% Create Modes of Operation file

%%
% Modes of Operation (MO)
% 0) Simulation Ended
% 1) Constant Current
% 2) Constant Voltage
% 3) Relaxation

function [MO_File] = getControlProfile()
j = 0;
% ---- Template ----
% j = j + 1;
% MO_File(j).MO         = ;   % Mode of Operation (MO)
% MO_File(j).CorD       = ''; % Charge or Discharge
% MO_File(j).C_rate     = ;   % C-rate (for CC) ABSOLUTE VALUE
% MO_File(j).Volt_ref   = ;   % [V] Reference voltage (for CV)
% MO_File(j).Volt_lim   = ;   % [V] Voltage Limit (for CC)
% MO_File(j).Time_lim   = ;   % [s] Time Limit (for all MO)
% MO_File(j).delta_tol  = ;   % Change in SV must be below this (For Relax)
% MO_File(j).C_rate_sat = ;   % Saturation C-rate (for CV) ABSOLUTE VALUE
% _______________________________________________________________________


%% # Chg Dcgh Cycles
for jj = 1:1
    % Constant Current Charge C/5
    j = j + 1;
    MO_File(j).MO        = 1;
    MO_File(j).CorD      = 'C';
    MO_File(j).C_rate    = 1/5;
    MO_File(j).Volt_ref  = [];   % [V]
    MO_File(j).Volt_lim  = 4.2;  % [V]
    MO_File(j).Time_lim  = 3600*5;  % [s]
    MO_File(j).delta_tol = [];

    % Constant Current Discharge C/5
    j = j + 1;
    MO_File(j).MO        = 1;
    MO_File(j).CorD      = 'D';
    MO_File(j).C_rate    = 1/5;
    MO_File(j).Volt_ref  = [];   % [V]
    MO_File(j).Volt_lim  = 3.2;  % [V]
    MO_File(j).Time_lim  = 3600*5;  % [s]
    MO_File(j).delta_tol = [];
end

%% GITT
    % % Relaxation
    % j = j + 1;
    % MO_File(j).MO         = 3;
    % MO_File(j).CorD       = 'C';
    % MO_File(j).C_rate     = 0;
    % MO_File(j).Volt_ref   = []; % [V]
    % MO_File(j).Volt_lim   = 4.4;  % [V]
    % MO_File(j).Time_lim   = 5; % [s]
    % MO_File(j).delta_tol  = 1e-10;
    % MO_File(j).C_rate_sat = 2;
    % 
    % for jj = 1:34
    % % Constant Current Charge C/17.25
    % j = j + 1;
    % MO_File(j).MO        = 1;
    % MO_File(j).CorD      = 'C';
    % MO_File(j).C_rate    = 1/17.25;
    % MO_File(j).Volt_ref  = [];
    % MO_File(j).Volt_lim  = 4.4;
    % MO_File(j).Time_lim  = 3600*0.5;
    % MO_File(j).delta_tol = [];
    % 
    % % Relaxation
    % j = j + 1;
    % MO_File(j).MO         = 3;
    % MO_File(j).CorD       = 'C';
    % MO_File(j).C_rate     = 0;
    % MO_File(j).Volt_ref   = []; % [V]
    % MO_File(j).Volt_lim   = 4.4;  % [V]
    % MO_File(j).Time_lim   = 3600*2; % [s]
    % MO_File(j).delta_tol  = 1e-10;
    % MO_File(j).C_rate_sat = 2;
    % end
    % 
    % % for jj = 1:34
    % % % Constant Current Charge 1C
    % % j = j + 1;
    % % MO_File(j).MO        = 1;
    % % MO_File(j).CorD      = 'D';
    % % MO_File(j).C_rate    = 1/17.25;
    % % MO_File(j).Volt_ref  = [];
    % % MO_File(j).Volt_lim  = 3.2;
    % % MO_File(j).Time_lim  = 3600*0.5;
    % % MO_File(j).delta_tol = [];
    % % 
    % % % Relaxation
    % % j = j + 1;
    % % MO_File(j).MO         = 3;
    % % MO_File(j).CorD       = 'C';
    % % MO_File(j).C_rate     = 0;
    % % MO_File(j).Volt_ref   = []; % [V]
    % % MO_File(j).Volt_lim   = 3.2;  % [V]
    % % MO_File(j).Time_lim   = 3600*2; % [s]
    % % MO_File(j).delta_tol  = 1e-10;
    % % MO_File(j).C_rate_sat = 2;
    % % end

%% Stair Step with intial and ending Relax
    % % Relaxation
    % j = j + 1;
    % MO_File(j).MO         = 3;
    % MO_File(j).CorD       = 'C';
    % MO_File(j).C_rate     = 0;
    % MO_File(j).Volt_ref   = []; % [V]
    % MO_File(j).Volt_lim   = 4.4;  % [V]
    % MO_File(j).Time_lim   = 5*60; %200; % [s]
    % MO_File(j).delta_tol  = 1e-10;
    % MO_File(j).C_rate_sat = 2;
    % 
    % % Constant Current Charge 1C
    % j = j + 1;
    % MO_File(j).MO        = 1;
    % MO_File(j).CorD      = 'C';
    % MO_File(j).C_rate    = 1;
    % MO_File(j).Volt_ref  = [];
    % MO_File(j).Volt_lim  = 3.9;
    % MO_File(j).Time_lim  = 5*60;%3600*10;
    % MO_File(j).delta_tol = [];
    % 
    % % Constant Current Charge C/2
    % j = j + 1;
    % MO_File(j).MO        = 1;
    % MO_File(j).CorD      = 'C';
    % MO_File(j).C_rate    = 1/2;
    % MO_File(j).Volt_ref  = [];
    % MO_File(j).Volt_lim  = 4.0;
    % MO_File(j).Time_lim  = 4*5*60;%3600*100;
    % MO_File(j).delta_tol = [];
    % 
    % % Constant Current Charge C/3
    % j = j + 1;
    % MO_File(j).MO        = 1;
    % MO_File(j).CorD      = 'C';
    % MO_File(j).C_rate    = 1/3;
    % MO_File(j).Volt_ref  = [];
    % MO_File(j).Volt_lim  = 4.1;
    % MO_File(j).Time_lim  = 4*5*60;%3600*100;
    % MO_File(j).delta_tol = [];
    % 
    % % % Constant Current Charge C/20 (Step Small Charge)
    % % j = j + 1;
    % % MO_File(j).MO        = 1;
    % % MO_File(j).CorD      = 'C';
    % % MO_File(j).C_rate    = 1/20;
    % % MO_File(j).Volt_ref  = [];
    % % MO_File(j).Volt_lim  = 4.4;
    % % MO_File(j).Time_lim  = 12*5*60;%1000;
    % % MO_File(j).delta_tol = [];
    % 
    % % Constant Current Discharge C/20 (Step Small Discharge)
    % j = j + 1;
    % MO_File(j).MO        = 1;
    % MO_File(j).CorD      = 'D';
    % MO_File(j).C_rate    = 1/20;
    % MO_File(j).Volt_ref  = [];
    % MO_File(j).Volt_lim  = 3.4;
    % MO_File(j).Time_lim  = 12*5*60;%1000;
    % MO_File(j).delta_tol = [];
    % 
    % % Relaxation
    % j = j + 1;
    % MO_File(j).MO         = 3;
    % MO_File(j).CorD       = 'C';
    % MO_File(j).C_rate     = 0;
    % MO_File(j).Volt_ref   = []; % [V]
    % MO_File(j).Volt_lim   = 4.2;  % [V]
    % MO_File(j).Time_lim   = 5*60; % [s]
    % MO_File(j).delta_tol  = 1e-10;
    % MO_File(j).C_rate_sat = 2;


%% Stair Step, No Relax
    % % Constant Current Charge 1C
    % j = j + 1;
    % MO_File(j).MO        = 1;
    % MO_File(j).CorD      = 'C';
    % MO_File(j).C_rate    = 1;
    % MO_File(j).Volt_ref  = [];
    % MO_File(j).Volt_lim  = 3.7;
    % MO_File(j).Time_lim  = 3600*1;
    % MO_File(j).delta_tol = [];
    % 
    % % Constant Current Charge C/2
    % j = j + 1;
    % MO_File(j).MO        = 1;
    % MO_File(j).CorD      = 'C';
    % MO_File(j).C_rate    = 1/2;
    % MO_File(j).Volt_ref  = [];
    % MO_File(j).Volt_lim  = 3.9;
    % MO_File(j).Time_lim  = 3600*1;
    % MO_File(j).delta_tol = [];
    % 
    % % Constant Current Charge C/3
    % j = j + 1;
    % MO_File(j).MO        = 1;
    % MO_File(j).CorD      = 'C';
    % MO_File(j).C_rate    = 1/3;
    % MO_File(j).Volt_ref  = [];
    % MO_File(j).Volt_lim  = 4.1;
    % MO_File(j).Time_lim  = 3600*1;
    % MO_File(j).delta_tol = [];

%% Charge to 50% SOC
    % % Relaxation
    % j = j + 1;
    % MO_File(j).MO         = 3;
    % MO_File(j).CorD       = 'C';
    % MO_File(j).C_rate     = 0;
    % MO_File(j).Volt_ref   = []; % [V]
    % MO_File(j).Volt_lim   = 4.4;  % [V]
    % MO_File(j).Time_lim   = 10; % [s]
    % MO_File(j).delta_tol  = 1e-10;
    % MO_File(j).C_rate_sat = 2;
    % 
    % % Constant Current Charge C/20
    % j = j + 1;
    % MO_File(j).MO        = 1;
    % MO_File(j).CorD      = 'C';
    %     % MO_File(j).C_rate    = 1/20; % 0.05
    %     % MO_File(j).C_rate    = 1/10; % 0.10
    %     % MO_File(j).C_rate    = 1/5; % 0.20
    %     % MO_File(j).C_rate    = 1/2; % 0.50
    %     % MO_File(j).C_rate    = 1; % 1.00
    %     MO_File(j).C_rate    = 2; % 2.00
    % MO_File(j).Volt_ref  = [];
    %     % MO_File(j).Volt_lim  = 3.553062401138020; % 10% SOC equilibrium voltage 
    %     % MO_File(j).Volt_lim  = 3.635333656536661; % 25% SOC equilibrium voltage 
    %     % MO_File(j).Volt_lim  = 3.749054804993171; % 50% SOC equilibrium voltage 
    %     % MO_File(j).Volt_lim  = 3.954758227144845; % 75% SOC equilibrium voltage 
    %     MO_File(j).Volt_lim  = 4.082762550630444; % 90% SOC equilibrium voltage 
    % 
    % MO_File(j).Time_lim  = 100*3600;
    % MO_File(j).delta_tol = [];
    % 
    % % Relaxation
    % j = j + 1;
    % MO_File(j).MO         = 3;
    % MO_File(j).CorD       = 'C';
    % MO_File(j).C_rate     = 0;
    % MO_File(j).Volt_ref   = []; % [V]
    % MO_File(j).Volt_lim   = 4.4;  % [V]
    % MO_File(j).Time_lim   = 4*3600; % [s]
    % MO_File(j).delta_tol  = 1e-10;
    % MO_File(j).C_rate_sat = 2;
    % 
    % % Constant Current Charge C/20 after Relaxing
    % j = j + 1;
    % MO_File(j).MO        = 1;
    % MO_File(j).CorD      = 'C';
    % MO_File(j).C_rate    = 1/20;
    % MO_File(j).Volt_ref  = [];
    % MO_File(j).Volt_lim  = 4.5;
    % MO_File(j).Time_lim  = 4*3600;
    % MO_File(j).delta_tol = [];


%% Impulse Simulation
    % initial_relax_samples = 2;
    % Ts = 1.0;   % [s], discrete time sampling rate
    % N_meas = 300-1-initial_relax_samples; % Number of discrete values to measure
    % N_meas = 100;
    % 
    % % Relaxation
    % j = j + 1;
    % MO_File(j).MO         = 3;
    % MO_File(j).CorD       = 'D';
    % MO_File(j).C_rate     = 0;
    % MO_File(j).Volt_ref   = []; % [V]
    % MO_File(j).Volt_lim   = 4.4;  % [V]
    % MO_File(j).Time_lim   = initial_relax_samples*Ts; % [s]
    % MO_File(j).delta_tol  = 1e-10;
    % MO_File(j).C_rate_sat = 2;
    % 
    % % Constant Current Charge C/20
    % j = j + 1;
    % MO_File(j).MO        = 1;
    % MO_File(j).CorD      = 'D';
    % MO_File(j).C_rate    = 0.040843474405010; % Results in i = 1 A/m^2
    % MO_File(j).Volt_ref  = [];
    % MO_File(j).Volt_lim  = 4.2;
    % MO_File(j).Time_lim  = Ts;
    % MO_File(j).delta_tol = [];
    % 
    % % Relaxation
    % j = j + 1;
    % MO_File(j).MO         = 3;
    % MO_File(j).CorD       = 'D';
    % MO_File(j).C_rate     = 0;
    % MO_File(j).Volt_ref   = []; % [V]
    % MO_File(j).Volt_lim   = 4.4;  % [V]
    % MO_File(j).Time_lim   = N_meas*Ts; % [s]
    % MO_File(j).delta_tol  = 1e-10;
    % MO_File(j).C_rate_sat = 2;


%% Impulse Simulation (Longer)
    % % %(steps at t^-) Didn't change anything
    % % Sampling Time
    %     N_t_s = 25;
    %     T_s_min = -1; % T_s = 10^(T_s_min)
    %     T_s_max =  1; % T_s = 10^(T_s_max)
    %     T_s_vec = logspace(T_s_min,T_s_max,N_t_s);
    %     Ts = T_s_vec(25);   % [s], discrete time sampling rate
    %     % Ts = 1.0;   % [s], discrete time sampling rate
    % 
    % initial_relax_samples = 2;
    % N_meas = 300-1-initial_relax_samples; % Number of discrete values to measure
    % N_meas = 800;
    % 
    % % Relaxation
    % j = j + 1;
    % MO_File(j).MO         = 3;
    % MO_File(j).CorD       = 'D';
    % MO_File(j).C_rate     = 0;
    % MO_File(j).Volt_ref   = []; % [V]
    % MO_File(j).Volt_lim   = 4.4;  % [V]%4.4
    % MO_File(j).Time_lim   = initial_relax_samples*Ts; % [s]
    % MO_File(j).delta_tol  = 1e-10;
    % MO_File(j).C_rate_sat = 2;
    % 
    % % Constant Current 
    % j = j + 1;
    % MO_File(j).MO        = 1;
    % MO_File(j).CorD      = 'D';
    % MO_File(j).C_rate    = 0.040843474405010; % Results in i = 1 A/m^2
    % MO_File(j).Volt_ref  = [];
    % MO_File(j).Volt_lim  = 5.2;%4.4
    % MO_File(j).Time_lim  = Ts;
    % MO_File(j).delta_tol = [];
    % 
    % % Relaxation
    % j = j + 1;
    % MO_File(j).MO         = 3;
    % MO_File(j).CorD       = 'D';
    % MO_File(j).C_rate     = 0;
    % MO_File(j).Volt_ref   = []; % [V]
    % MO_File(j).Volt_lim   = 5.4;  % [V]%4.4
    % MO_File(j).Time_lim   = N_meas*Ts; % [s]
    % MO_File(j).delta_tol  = 1e-10;
    % MO_File(j).C_rate_sat = 2;


%% Step with Initial relax
    % initial_relax_samples = 2;
    % Ts = 10.0;   % [s], discrete time sampling rate
    % N_meas = 600; % Number of discrete values to measure
    % 
    % % Relaxation
    % j = j + 1;
    % MO_File(j).MO         = 3;
    % MO_File(j).CorD       = 'D';
    % MO_File(j).C_rate     = 0;
    % MO_File(j).Volt_ref   = []; % [V]
    % MO_File(j).Volt_lim   = 4.4;  % [V]
    % MO_File(j).Time_lim   = initial_relax_samples*Ts; % [s]
    % MO_File(j).delta_tol  = 1e-10;
    % MO_File(j).C_rate_sat = 2;
    % 
    % % Constant Current Charge C/20
    % j = j + 1;
    % MO_File(j).MO        = 1;
    % MO_File(j).CorD      = 'D';
    % MO_File(j).C_rate    = 0.040843474405010; % Results in i = 1 A/m^2
    % MO_File(j).Volt_ref  = [];
    % MO_File(j).Volt_lim  = 3.8;%4.2 (Max), 3.8 (Min)
    % MO_File(j).Time_lim  = N_meas*Ts;
    % MO_File(j).delta_tol = [];

%%
    % Constant Current Charge C/20
    % j = j + 1;
    % MO_File(j).MO        = 1;
    % MO_File(j).CorD      = 'C';
    % MO_File(j).C_rate    = 1;
    % MO_File(j).Volt_ref  = [];
    % MO_File(j).Volt_lim  = 3.75;
    % MO_File(j).Time_lim  = 3600*20;
    % MO_File(j).delta_tol = [];
    
    % % Constant Current Discharge 1C/3
    % j = j + 1;
    % MO_File(j).MO        = 1;
    % MO_File(j).CorD      = 'D';
    % MO_File(j).C_rate    = 1/20;
    % MO_File(j).Volt_ref  = [];
    % MO_File(j).Volt_lim  = 4.0;
    % MO_File(j).Time_lim  = 3600*20;
    % MO_File(j).delta_tol = [];
    % 
    % % Constant Current Discharge 1C/3
    % j = j + 1;
    % MO_File(j).MO        = 1;
    % MO_File(j).CorD      = 'D';
    % MO_File(j).C_rate    = 1/50;
    % MO_File(j).Volt_ref  = [];
    % MO_File(j).Volt_lim  = 4.0;
    % MO_File(j).Time_lim  = 3600*50;
    % MO_File(j).delta_tol = [];

%% Voltage Step
    % % Constant Voltage @ 4.0V
    % j = j + 1;
    % MO_File(j).MO         = 2;
    % MO_File(j).CorD       = [];
    % MO_File(j).C_rate     = [];
    % MO_File(j).Volt_ref   = 4.0; % [V]
    % MO_File(j).Volt_lim   = [];  % [V]
    % MO_File(j).Time_lim   = 0.5; % [s]
    % MO_File(j).delta_tol  = 1e-10;
    % MO_File(j).C_rate_sat = 2;
    % 
    % % Constant Voltage @ 3.98V
    % j = j + 1;
    % MO_File(j).MO         = 2;
    % MO_File(j).CorD       = [];
    % MO_File(j).C_rate     = [];
    % MO_File(j).Volt_ref   = 3.98;  % [V]
    % MO_File(j).Volt_lim   = [];   % [V]
    % MO_File(j).Time_lim   = 0.05; % [s]
    % MO_File(j).delta_tol  = 1e-10;
    % MO_File(j).C_rate_sat = 2;
    
    % % Constant Voltage @ 4.0V
    % j = j + 1;
    % MO_File(j).MO         = 2;
    % MO_File(j).CorD       = [];
    % MO_File(j).C_rate     = [];
    % MO_File(j).Volt_ref   = 4.0;  % [V]
    % MO_File(j).Volt_lim   = [];   % [V]
    % MO_File(j).Time_lim   = 60; % [s]
    % MO_File(j).delta_tol  = 1e-10;
    % MO_File(j).C_rate_sat = 2;

%% CC to CV
    % % Constant Current Charge 
    % j = j + 1;
    % MO_File(j).MO        = 1;
    % MO_File(j).CorD      = 'C';
    % MO_File(j).C_rate    = 2;
    % MO_File(j).Volt_ref  = [];
    % MO_File(j).Volt_lim  = 4.2;%3.749054804993171;
    % MO_File(j).Time_lim  = 3600*1;
    % MO_File(j).delta_tol = [];
    % 
    % % Constant Voltage @ 50% SOC
    % j = j + 1;
    % MO_File(j).MO         = 2;
    % MO_File(j).CorD       = [];
    % MO_File(j).C_rate     = [];
    % MO_File(j).Volt_ref   = 4.2;%3.749054804993171;  % [V]
    % MO_File(j).Volt_lim   = [];   % [V]
    % MO_File(j).Time_lim   = 2000; % [s]
    % MO_File(j).delta_tol  = 1e-10;
    % MO_File(j).C_rate_sat = 2;
    
    % % Constant Current Charge 1C
    % j = j + 1;
    % MO_File(j).MO        = 1;
    % MO_File(j).CorD      = 'D';
    % MO_File(j).C_rate    = 1;
    % MO_File(j).Volt_ref  = [];
    % MO_File(j).Volt_lim  = 3.44;
    % MO_File(j).Time_lim  = 3600*1;
    % MO_File(j).delta_tol = [];
end