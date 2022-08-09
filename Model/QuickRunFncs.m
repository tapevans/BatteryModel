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

%% Load Results
clear all
close all
clc

filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\TestingForTyrone\TestingForTyrone_KPCont_SingleStep200sRelaxSOC80.mat';
load(filename)

%% Make SS from SV during simulation
% Load results from a simulation using the previous section, then set the
% index in the next line from which time you'd like to pull SV from
SIM.SV_IC = SV_soln(320,:)';
SIM.freq = logspace(-1,11,101);
SIM.i_user = -SIM.Cell_Cap/20;
N.N_In = 1;
N.N_Out = 5;

SIM.OutputMatrix = zeros(N.N_Out , N.N_SV_tot);
    j = 0;
    % Cell Voltage
    j = j+1;
    idx_phi_ed_AN = P.phi_ed;
    i = N.N_CV_CA(end);
    index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
    idx_phi_ed_CA = index_offset + P.phi_ed;
    SIM.OutputMatrix(j,idx_phi_ed_AN) = -1;
    SIM.OutputMatrix(j,idx_phi_ed_CA) =  1;
    % @AN/SEP
    i = N.N_CV_AN(end);
    index_offset = (i-1)*N.N_SV_AN;
    % Delta Phi   @AN/SEP
    j = j+1;
    idx = index_offset + P.del_phi;
    SIM.OutputMatrix(j,idx) =  1;
    % Temperature @AN/SEP
    j = j+1;
    idx = index_offset + P.T;
    SIM.OutputMatrix(j,idx) = 1;
    % C_Liion     @AN/SEP
    j = j+1;
    idx = index_offset + P.C_Liion;
    SIM.OutputMatrix(j,idx) = 1;
    % X_surf      @AN/SEP
    j = j+1;
    idx = index_offset + P.C_Li_surf_AN;
    SIM.OutputMatrix(j,idx) = 1/AN.C_Li_max;

i = 1;
P.SS.omega    = i; i = i + 1;
P.SS.Z_mag    = i; i = i + 1;
P.SS.Z_Re     = i; i = i + 1;
P.SS.Z_Im     = i; i = i + 1;
P.SS.Z_dB     = i; i = i + 1;
P.SS.Z_ps_deg = i; i = i + 1;

[A,B,C,D,Z_results] = getSSImpedance(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS);

save_filename = [num2str(SIM.SOC_start) 'SOC.mat'];
multiple = -SIM.A_c^-1;
sys = multiple*ss(A,B,C,D,'E',SIM.M);
save(save_filename, 'sys')



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

% Model
filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\TestingForTyrone\TestingForTyrone_KPCont_SingleStep200sRelaxSOC80.mat';
model = load(filename);
idx = find(model.t_soln == 200);
idx = idx(1);
OneCRate = model.SIM.Cell_Cap;

% Linear SS Model
% SS_sys = load('F:\TylerFiles\GitHubRepos\BatteryModelingExtras\DataToTyrone\2202_08_09_Lui_SS\SS80SOC.mat');
% SS_sys = load('F:\TylerFiles\GitHubRepos\BatteryModelingExtras\DataToTyrone\2202_08_09_Lui_SS\79.9SOC.mat');
% SS_sys = load('F:\TylerFiles\GitHubRepos\BatteryModelingExtras\DataToTyrone\2202_08_09_Lui_SS\79.95SOC.mat');
SS_sys = load('F:\TylerFiles\GitHubRepos\BatteryModelingExtras\DataToTyrone\2202_08_09_Lui_SS\80SOCfromC3Sim.mat');
tFinal = 13*60;
Amplitude = OneCRate/3;
opt = stepDataOptions;
opt.StepAmplitude = Amplitude;

% figure
[y,tOut] = step(SS_sys.sys,tFinal,opt);
SS_sys.cell_voltage = y(:,1);
SS_sys.del_phi      = y(:,2);
SS_sys.TempK        = y(:,3);
SS_sys.C_Liion      = y(:,4);
SS_sys.X_surf       = y(:,5);

% close all;
% Cell Voltage (min)
figure
hold on
plot((model.t_soln(idx:end) - model.t_soln(idx))/60 , model.cell_voltage(idx:end) , 'or', 'LineWidth',2,'DisplayName','Model')
plot((tOut)/60 , SS_sys.cell_voltage + model.cell_voltage(idx,1), 'k' , 'LineWidth',2,'DisplayName','SS')
xlabel('Time (min)')
ylabel('Voltage (V)')
title('Cell Voltage')
xlim([0,tFinal/60])
lgn = legend;

% Del_phi (min)
figure
hold on
plot((model.t_soln(idx:end) - model.t_soln(idx))/60, model.del_phi(idx:end,model.N.N_CV_AN) , 'or', 'LineWidth',2,'DisplayName','Model')
plot((tOut)/60 , SS_sys.del_phi + model.del_phi(idx,model.N.N_CV_AN) , 'k' , 'LineWidth',2,'DisplayName','SS')
xlabel('Time (min)')
ylabel('Voltage (V)')
title('\Delta \phi')
xlim([0,(tFinal)/60])
lgn = legend;

% % Cell Voltage (s)
% figure
% hold on
% plot(model.t_soln(idx:end) - model.t_soln(idx) , model.cell_voltage(idx:end) , 'or', 'LineWidth',2,'DisplayName','Model')
% plot(tOut , SS_sys.cell_voltage + model.cell_voltage(idx,1), 'k' , 'LineWidth',2,'DisplayName','SS')
% xlabel('Time (min)')
% ylabel('Voltage (V)')
% title('Cell Voltage')
% xlim([0,tFinal])
% lgn = legend;
% 
% % Del_phi (s)
% figure
% hold on
% plot(model.t_soln(idx:end) - model.t_soln(idx), model.del_phi(idx:end,model.N.N_CV_AN) , 'or', 'LineWidth',2,'DisplayName','Model')
% plot(tOut , SS_sys.del_phi + model.del_phi(idx,model.N.N_CV_AN) , 'k' , 'LineWidth',2,'DisplayName','SS')
% xlabel('Time (min)')
% ylabel('Voltage (V)')
% title('\Delta \phi')
% xlim([0,tFinal])
% lgn = legend;

% % Temperature
% figure
% hold on
% plot(model.t_soln(idx:end) - model.t_soln(idx) , model.TemperatureK(idx:end,model.N.N_CV_AN) , 'or', 'LineWidth',2,'DisplayName','Model')
% plot(tOut , SS_sys.TempK + model.TemperatureK(idx,model.N.N_CV_AN) , 'k' , 'LineWidth',2,'DisplayName','SS')
% xlabel('Time (s)')
% ylabel('Temperature (K)')
% title('Temperature')
% xlim([0,tFinal])
% lgn = legend;

% C_Liion
figure
hold on
plot(model.t_soln(idx:end) - model.t_soln(idx) , model.C_Liion(idx:end,model.N.N_CV_AN) , 'or', 'LineWidth',2,'DisplayName','Model')
plot(tOut , SS_sys.C_Liion + model.C_Liion(idx,model.N.N_CV_AN) , 'k' , 'LineWidth',2,'DisplayName','SS')
xlabel('Time (s)')
ylabel('Concentration (kmol/m^3)')
title('C_{Li^+}')
xlim([0,tFinal])
lgn = legend;
% 
% % X_surf
% figure
% hold on
% plot(model.t_soln(idx:end) - model.t_soln(idx) , model.X_Li_surf(idx:end,model.N.N_CV_AN) , 'or', 'LineWidth',2,'DisplayName','Model')
% plot(tOut , SS_sys.X_surf + model.X_Li_surf(idx,model.N.N_CV_AN) , 'k' , 'LineWidth',2,'DisplayName','SS')
% xlabel('Time (s)')
% ylabel('Mole Frac (-)')
% title('X_{surf}')
% xlim([0,tFinal])
% lgn = legend;
