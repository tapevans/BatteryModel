clear all
close all
clc
%%
% Test Initialization
%%% Mode 1 ---- Polarization ----
SIM.SOC_start = 50;   % [%], Initial state of charge
SIM.SimMode = 1;
SIM.C_rate = 1/20;
SIM.ChargeOrDischarge = 1;

%%% Mode 2 ---- Harmonic Perturbation ----

%%% Mode 3 ---- State Space EIS ----
% SIM.SimMode = 3;
% SIM.freq      = 1e-2; % [rad/s], frequency of the sin wave 
% SIM.SOC_start = 95;   % [%], Initial state of charge

%%% Mode 4 ---- Known BC Profile Controller ----

%%% Mode 5 ---- MOO Controller ----

%%% Mode 6 ---- Simulink ----

%%% Mode 7 ---- Manual Profile ----


%%%
[AN,CA,SEP,EL,SIM,N,FLAG] = batt_inputs(SIM);
%
[AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS] = batt_init(AN,CA,SEP,EL,SIM,N,FLAG);

cell_voltage = SIM.SV_IC(341) - SIM.SV_IC(3)
%% Test Governing Eqns Output
t = 2;
i_user = 14;
SV = SIM.SV_IC;
dSVdt = batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user);

%% Run Single Battery Batch
clear all
close all
clc

CreateProject
%
RunSimulations


%% Post-Processing
filename = 'F:\TylerFiles\GitHubRepos\p2d-model\BatteryModel\BatchMode_DAE\Results\Half_Cell_Test\Half_Cell_Test_halfCA_Polar_1.00C_D.mat';
postProcessing(filename)

%% Plot Single Results
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\KBCP_Mode_Test\KBCP_Mode_Test_KPCont_Profile_CCChg4.2_CCDchg3.4SOC0.mat';
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\KBCP_Mode_Test\KBCP_Mode_Test_KPCont_Profile_CCChg4.2_CCDchg3.4_C3SOC0.mat';
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\MoleInitFix\MoleInitFix_KPCont_Profile_CC_Test_3Step_wRelaxSOC0.mat';
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\TestingForTyrone\TestingForTyrone_KPCont_Profile_CC_Test_3Step_wRelaxSOC50.mat';
filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\Semi_Explicit_Test\Semi_Explicit_Test_Polar_1.00C_D.mat';
plotfcn(filename)

%%
clear all
close all
clc

% Load Results
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\Final_Lui_Wiley_Model\Final_Lui_Wiley_Model_KPCont_Profile_CV_Test_1SmallStepSOC81.93.mat';

% filename = 'F:\TylerFiles\GitHubRepos\BatteryModelingExtras\DataToTyrone\2022_07_13_Lui_SS\StairStepSim.mat';
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModelingExtras\DataToTyrone\2022_07_11_Lui_SS\50SOC.mat';

filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\Semi_Explicit_Test\Semi_Explicit_Test_SS_EIS_SOC95.mat';
load(filename)


%% Change Parameters of an existing sim
clear all; close all; clc;
%Thinking the simulations are struggling to make it past the knee of the
%voltage curve so I'm thinking I'll adjust the final time
filename = 'F:\TylerFiles\GitHubRepos\p2d-model\BatteryModel\BatchMode_DAE\Results\Carderock_C_dl_ANx3_C_dl_CAx3\C_dl_ANx3_C_dl_CAx3_Polar_0.33C_D.mat';
% filename = 'F:\TylerFiles\GitHubRepos\p2d-model\BatteryModel\BatchMode_DAE\Results\Carderock_C_dl_ANx3\C_dl_ANx3_Polar_0.33C_D.mat';
% filename = 'F:\TylerFiles\GitHubRepos\p2d-model\BatteryModel\BatchMode_DAE\Results\Carderock_C_dl_CAx3\C_dl_CAx3_Polar_0.33C_D.mat';
load(filename)
t_final = 8100;
SIM.tspan = [0, t_final];

save(filename,'AN','CA','SEP','EL','SIM','CONS','P','N','FLAG')


%% Comparing Step Response to Continuous Model
clear all;
close all;
clc;
SS_sys = load('F:\TylerFiles\GitHubRepos\BatteryModelingExtras\DataToTyrone\2022_07_15_Lui_SS\5SOC.mat');
tFinal = 1e-6;
figure
[y,tOut] = step(SS_sys.sys,tFinal);
SS_sys.cell_voltage = y(:,1);
SS_sys.del_phi      = y(:,2);
SS_sys.TempK        = y(:,3);
SS_sys.C_Liion      = y(:,4);
SS_sys.X_surf       = y(:,5);

filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\TestingForTyrone\TestingForTyrone_KPCont_Profile_CC_Test_StepResponse_0.4CSOC50.mat';
model = load(filename);
idx = find(model.t_soln == 1);
idx = idx(1);
close all;

% Cell Voltage
figure
hold on
plot(model.t_soln(idx:end) - model.t_soln(idx) , model.cell_voltage(idx:end) , 'ob', 'LineWidth',2,'DisplayName','Model')
plot(tOut , SS_sys.cell_voltage + model.cell_voltage(idx,1), 'k' , 'LineWidth',2,'DisplayName','SS')
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Cell Voltage')
xlim([0,1e-6])
lgn = legend;

% Del_phi
figure
hold on
plot(model.t_soln(idx:end) - model.t_soln(idx), model.del_phi(idx:end,model.N.N_CV_AN) , 'ob', 'LineWidth',2,'DisplayName','Model')
plot(tOut , SS_sys.del_phi + model.del_phi(idx,model.N.N_CV_AN) , 'k' , 'LineWidth',2,'DisplayName','SS')
xlabel('Time (s)')
ylabel('Voltage (V)')
title('\Delta \phi')
xlim([0,1e-6])
lgn = legend;

% Temperature
figure
hold on
plot(model.t_soln(idx:end) - model.t_soln(idx) , model.TemperatureK(idx:end,model.N.N_CV_AN) , 'ob', 'LineWidth',2,'DisplayName','Model')
plot(tOut , SS_sys.TempK + model.TemperatureK(idx,model.N.N_CV_AN) , 'k' , 'LineWidth',2,'DisplayName','SS')
xlabel('Time (s)')
ylabel('Temperature (K)')
title('Temperature')
xlim([0,1e-6])
lgn = legend;

% C_Liion
figure
hold on
plot(model.t_soln(idx:end) - model.t_soln(idx) , model.C_Liion(idx:end,model.N.N_CV_AN) , 'ob', 'LineWidth',2,'DisplayName','Model')
plot(tOut , SS_sys.C_Liion + model.C_Liion(idx,model.N.N_CV_AN) , 'k' , 'LineWidth',2,'DisplayName','SS')
xlabel('Time (s)')
ylabel('Concentration (kmol/m^3)')
title('C_{Li^+}')
xlim([0,1e-6])
lgn = legend;

% X_surf
figure
hold on
plot(model.t_soln(idx:end) - model.t_soln(idx) , model.X_Li_surf(idx:end,model.N.N_CV_AN) , 'ob', 'LineWidth',2,'DisplayName','Model')
plot(tOut , SS_sys.X_surf + model.X_Li_surf(idx,model.N.N_CV_AN) , 'k' , 'LineWidth',2,'DisplayName','SS')
xlabel('Time (s)')
ylabel('Mole Frac (-)')
title('X_{surf}')
xlim([0,1e-6])
lgn = legend;
