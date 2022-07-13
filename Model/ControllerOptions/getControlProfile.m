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
% MO_File(j).MO         = ;   % Mode of Operation
% MO_File(j).CorD       = ''; % Charge or Discharge
% MO_File(j).C_rate     = ;   % C-rate (for CC) ABSOLUTE VALUE
% MO_File(j).Volt_ref   = ;   % [V] Reference voltage (for CV)
% MO_File(j).Volt_lim   = ;   % [V] Voltage Limit (for CC)
% MO_File(j).Time_lim   = ;   % [s] Time Limit (for all MO)
% MO_File(j).delta_tol  = ;   % Change in SV must be below this (For Relax)
% MO_File(j).C_rate_sat = ;   % Saturation C-rate (for CV) ABSOLUTE VALUE
% _______________________________________________________________________

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
% for jj = 1:34
% % Constant Current Charge 1C
% j = j + 1;
% MO_File(j).MO        = 1;
% MO_File(j).CorD      = 'D';
% MO_File(j).C_rate    = 1/17.25;
% MO_File(j).Volt_ref  = [];
% MO_File(j).Volt_lim  = 3.2;
% MO_File(j).Time_lim  = 3600*0.5;
% MO_File(j).delta_tol = [];
% 
% % Relaxation
% j = j + 1;
% MO_File(j).MO         = 3;
% MO_File(j).CorD       = 'C';
% MO_File(j).C_rate     = 0;
% MO_File(j).Volt_ref   = []; % [V]
% MO_File(j).Volt_lim   = 3.2;  % [V]
% MO_File(j).Time_lim   = 3600*2; % [s]
% MO_File(j).delta_tol  = 1e-10;
% MO_File(j).C_rate_sat = 2;
% end

%% CC
% % Constant Current Charge C/2
% j = j + 1;
% MO_File(j).MO        = 1;
% MO_File(j).CorD      = 'C';
% MO_File(j).C_rate    = 1/2;
% MO_File(j).Volt_ref  = [];
% MO_File(j).Volt_lim  = 3.9;
% MO_File(j).Time_lim  = 3600*1;
% MO_File(j).delta_tol = [];

% % Constant Current Charge C/3
% j = j + 1;
% MO_File(j).MO        = 1;
% MO_File(j).CorD      = 'C';
% MO_File(j).C_rate    = 1/3;
% MO_File(j).Volt_ref  = [];
% MO_File(j).Volt_lim  = 4.1;
% MO_File(j).Time_lim  = 3600*1;
% MO_File(j).delta_tol = [];

%% Voltage Step
% Constant Voltage @ 4.0V
j = j + 1;
MO_File(j).MO         = 2;
MO_File(j).CorD       = [];
MO_File(j).C_rate     = [];
MO_File(j).Volt_ref   = 4.0; % [V]
MO_File(j).Volt_lim   = [];  % [V]
MO_File(j).Time_lim   = 0.5; % [s]
MO_File(j).delta_tol  = 1e-10;
MO_File(j).C_rate_sat = 2;

% Constant Voltage @ 3.98V
j = j + 1;
MO_File(j).MO         = 2;
MO_File(j).CorD       = [];
MO_File(j).C_rate     = [];
MO_File(j).Volt_ref   = 3.98;  % [V]
MO_File(j).Volt_lim   = [];   % [V]
MO_File(j).Time_lim   = 0.05; % [s]
MO_File(j).delta_tol  = 1e-10;
MO_File(j).C_rate_sat = 2;

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
% % Constant Current Charge 1C
% j = j + 1;
% MO_File(j).MO        = 1;
% MO_File(j).CorD      = 'C';
% MO_File(j).C_rate    = 1;
% MO_File(j).Volt_ref  = [];
% MO_File(j).Volt_lim  = 4.2;
% MO_File(j).Time_lim  = 3600*1;
% MO_File(j).delta_tol = [];
% 
% % Constant Voltage @ 4.0V
% j = j + 1;
% MO_File(j).MO         = 2;
% MO_File(j).CorD       = [];
% MO_File(j).C_rate     = [];
% MO_File(j).Volt_ref   = 4.2;  % [V]
% MO_File(j).Volt_lim   = [];   % [V]
% MO_File(j).Time_lim   = 10; % [s]
% MO_File(j).delta_tol  = 1e-10;
% MO_File(j).C_rate_sat = 2;
end