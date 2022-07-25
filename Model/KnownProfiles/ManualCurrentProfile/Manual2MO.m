%% Make MO file from manual plating profile
% 
% 
% 
clear all;
close all;
clc;
%%
filename = 'F:\TylerFiles\GitHubRepos\BatteryModelingExtras\OLD MODELS\WORKINGP2D_BeforeMovingRepos\BatteryModel\BatchMode_DAE\Results\Manual_Profile_for_Plating_100steps_1000Iter_0.01tol\Manual_Profile_for_Plating_100steps_1000Iter_0.01tol_CurrentProfile.mat';
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModelingExtras\OLD MODELS\WORKINGP2D_BeforeMovingRepos\BatteryModel\BatchMode_DAE\Results\Manual_Profile_for_Plating_100steps_1000Iter_0.01tol\Manual_Profile_for_Plating_100steps_1000Iter_0.01tol.mat';
filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\Semi_Explicit_Test\Semi_Explicit_Test_SOC25_ManCurrProf_100steps_1000Iter_0.01tol_CurrentProfile_Output.mat';
load(filename)

save_filename = 'ManualMOSOC25.mat';

Cell_Cap = 24.48371531970692;

%%
j = 0;
% Relaxation
j = j + 1;
MO_File(j).MO         = 3;
MO_File(j).CorD       = 'C';
MO_File(j).C_rate     = 0;
MO_File(j).Volt_ref   = []; % [V]
MO_File(j).Volt_lim   = 4.2;  % [V]
MO_File(j).Time_lim   = 10; % [s]
MO_File(j).delta_tol  = 1e-10;
MO_File(j).C_rate_sat = 2;

% for i = 1:length(current_region)
for i = 1:length(region_current_vec)
    % Constant Current Charge
    j = j + 1;
    MO_File(j).MO        = 1;
    MO_File(j).CorD      = 'C';
%     MO_File(j).C_rate    = -current_region(i)/Cell_Cap;
    MO_File(j).C_rate    = abs(region_current_vec(i))/Cell_Cap;
    MO_File(j).Volt_ref  = [];
    MO_File(j).Volt_lim  = 4.2;
    MO_File(j).Time_lim  = 36;
    MO_File(j).delta_tol = [];
end

%% Save File
save(save_filename,"MO_File")
