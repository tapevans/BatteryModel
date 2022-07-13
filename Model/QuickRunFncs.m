clear all
close all
clc
%%
% Test Initialization
SIM.ChargeOrDischarge = 1;
SIM.SOC_start = 95;
% SIM.ChargeOrDischarge = -1; %Charge
% SIM.SOC_start = 5;
SIM.SimMode = 1;
SIM.C_rate = 1/20;
[AN,CA,SEP,EL,SIM,N,FLAG] = batt_inputs(SIM);
%
[AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS] = batt_init(AN,CA,SEP,EL,SIM,N,FLAG);

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
filename = 'F:\TylerFiles\GitHubRepos\p2d-model\BatteryModel\BatchMode_DAE\Results\KBCP_Mode_Test\KBCP_Mode_Test_KPCont_Profile_CC_Test_3StepSOC0.mat';
plotfcn(filename)

%%
clear all
close all
clc

%% Load Results
filename = 'F:\TylerFiles\GitHubRepos\p2d-model\BatteryModel\BatchMode_DAE\Results\KBCP_Mode_Test\KBCP_Mode_Test_KPCont_Profile_CC_Test_3StepSOC0.mat';
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModelingExtras\DataToTyrone\2022_07_11_Lui_SS\95SOC.mat';
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