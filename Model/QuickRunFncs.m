clear all
close all
clc
%%
% Test Initialization
%%% Mode 1 ---- Polarization ----
    % SIM.SOC_start = 50;   % [%], Initial state of charge
    % SIM.SimMode = 1;
    % SIM.C_rate = 1/20;
    % SIM.ChargeOrDischarge = 1;

%%% Mode 2 ---- Harmonic Perturbation ----

%%% Mode 3 ---- State Space EIS ----
    SIM.SimMode = 3;
    SIM.freq      = 1e-2; % [rad/s], frequency of the sin wave 
    SIM.SOC_start = 95;   % [%], Initial state of charge

%%% Mode 4 ---- Known BC Profile Controller ----

%%% Mode 5 ---- MOO Controller ----

%%% Mode 6 ---- Simulink ----

%%% Mode 7 ---- Manual Profile ----


%%%
[AN,CA,SEP,EL,SIM,N,FLAG] = batt_inputs(SIM);
%
[AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS] = batt_init(AN,CA,SEP,EL,SIM,N,FLAG);

% cell_voltage = SIM.SV_IC(341) - SIM.SV_IC(3)

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
clc; close all;
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\KBCP_Mode_Test\KBCP_Mode_Test_KPCont_Profile_CCChg4.2_CCDchg3.4SOC0.mat';
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\KBCP_Mode_Test\KBCP_Mode_Test_KPCont_Profile_CCChg4.2_CCDchg3.4_C3SOC0.mat';
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\MoleInitFix\MoleInitFix_KPCont_Profile_CC_Test_3Step_wRelaxSOC0.mat';
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\TestingForTyrone\TestingForTyrone_KPCont_Profile_CC_Test_3Step_wRelaxSOC50.mat';
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\Semi_Explicit_Test\Semi_Explicit_Test_Polar_1.00C_D.mat';
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


