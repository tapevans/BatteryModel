clear all
close all
clc


%
    inputHandle             = @batt_inputs_NREL;
    SIM.ANEqPotentialHandle = @E_eqGraphite_NREL;
    SIM.ANi_oHandle         = @i_oC6_NREL;
    SIM.ANsigmaHandle       = @sigmaC6_NREL;
    SIM.AND_oHandle         = @D_o_Graphite_NREL;
    SIM.CAEqPotentialHandle = @E_eqNMC_NREL;
    SIM.CAi_oHandle         = @i_oNMC_NREL;  
    SIM.CAsigmaHandle       = @sigmaNMC_NREL;
    SIM.CAD_oHandle         = @D_o_NMC532_NREL;
    SIM.ELtf_numHandle      = @transferenceNumber_NREL;
    SIM.ELActivityHandle    = @activity_NREL;
    SIM.ELD_o_Li_ionHandle  = @D_oLiion_NREL;
    SIM.ELkappaHandle       = @kappa_NREL;

% Test Initialization
%%% Mode 1 ---- Polarization ----
    SIM.SOC_start = 10;   % [%], Initial state of charge
    SIM.SimMode = 1;
    SIM.C_rate = 1/20;
    % SIM.C_rate = 0;
    SIM.ChargeOrDischarge = 1;

%%% Mode 2 ---- Harmonic Perturbation ----

%%% Mode 3 ---- State Space EIS ----
    % SIM.SimMode = 3;
    % SIM.freq      = 1e-2; % [rad/s], frequency of the sin wave 
    % SOC_start = 90
    % SIM.SOC_start = SOC_start;   % [%], Initial state of charge
    
%%% Mode 4 ---- Known BC Profile Controller ----

%%% Mode 5 ---- MOO Controller ----

%%% Mode 7 ---- Manual Profile ----

%%% Mode 8 ---- PRBS ----



%
[AN,CA,SEP,EL,SIM,N,FLAG] = inputHandle(SIM);
%
[AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS] = batt_init(AN,CA,SEP,EL,SIM,N,FLAG);

% cell_voltage = SIM.SV_IC(341) - SIM.SV_IC(3)
% % % x_surf = SIM.SV_IC(17)/AN.C_Li_max

%%% Mode 3
    % cell_voltage = SIM.OutputAtEquil(1)

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


%% Uncheck Post-Processing Complete
clear all; close all; clc;
filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\TestNewInputFileNoise\Test_KPCont_StairStepNoRelaxSOC0.mat';
load(filename,'postProcessComplete')
postProcessComplete = 0;
save(filename)


%% Post-Processing
clear all; close all; clc;
filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\OptiFinal\Standard_Polar_0.25C_C.mat';
postProcessing(filename)


%% Plot Single Results
clc; 
clear all;
close all;

% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\COETest\444All_Constant_Polar_1.00C_C.mat';
filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\TestMLI\test_MultiLevelInput_SOC_init50_pmSOC1_maxCrate1.mat';
plotfcn(filename)


%% Arrange Figures
FigArrange = 1;
if FigArrange == 1
    fig = gcf;
    NumFig = fig.Number;
    
    Ncol = 3;
    
    for i = 1:NumFig
        f = figure(i);
        k = mod(i-1,Ncol);
        row = mod(fix((i-1)/Ncol),2);
        if row == 0
            r = 575;
%             r = 540;
        elseif row == 1
            r = 62;
        end
        f.Position = [k*575+15 r 560 420];
    end
end


%% Load Results
clear all
close all
clc

filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\TestImpedanceContributions\Standard_SS_EIS_SOC95.mat';


load(filename)

% Resave with PostProcessingNotComplete
    % postProcessComplete = 0;
    % save(filename)

%% Create a SS representation from a SV during a simulation

SV = SV_soln(end,:)';
% SIM.I_user_amp = SIM.C_rate*SIM.Cell_Cap; %%% This doesn't work cause this is based on the MO file
SIM.I_user_amp = (1/20)*SIM.Cell_Cap; 
SIM.i_user_amp = SIM.I_user_amp/SIM.A_c;
SIM.i_user = SIM.i_user_amp;
i_user = SIM.i_user;

[SS_A,SS_B,SS_C,SS_D] = getSSfromSV(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,SV,i_user);

%%
%%%%%%%%!!!!!!!!!!!!!!!! May just get rid of this
multiple = -SIM.A_c^-1;
C = C(1,:); % Cell_Voltage is the first row
D = 0;
sys_E = multiple*ss(A,B,C,D,'E',SIM.M);
sys_min = minreal(sys_E);
LQGsys = sys_min;
    N_output = 1;

    rho = 1e0;
    Q = rho * (LQGsys.C' * LQGsys.C) + 1e-10*eye(length(LQGsys.A));
    R = 1e0; % Inputs is length 1, makes R (1x1) %%%%!!!! Needs better comment
    QXU = blkdiag(Q,R);
    QXU = eye(size(QXU));

    QWV = eye(size(QXU)); % Needs to be positive definite

    rho_I = 1e0;
    QI = rho_I*eye(N_output);


    reg_SS = lqg(LQGsys,QXU,QWV,QI);


%% Make SS from SV during simulation
% % Load results from a simulation using the previous section, then set the
% % index in the next line from which time you'd like to pull SV from
% SIM.SV_IC = SV_soln(end,:)';
% SIM.freq = logspace(-1,11,101);
% SIM.i_user = -SIM.Cell_Cap/20;
% 
% i = 1;
%     P.OM.cell_volt = i; i = i + 1;
%     P.OM.del_phi   = i; i = i + 1;
%     P.OM.temp      = i; i = i + 1;
%     P.OM.C_Liion   = i; i = i + 1;
%     P.OM.X_surf    = i; i = i + 1;
%     P.OM.i_Far     = i; i = i + 1;
%     P.OM.eta       = i; i = i + 1;
% 
% N.N_In = 1;
% N.N_Out = length(fieldnames(P.OM));
% 
%     SIM.OutputMatrix = zeros(N.N_Out , N.N_SV_tot);
%         % Cell Voltage
%             idx_phi_ed_AN = P.phi_ed;
% 
%             i = N.N_CV_CA(end);
%             index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
%             idx_phi_ed_CA = index_offset + P.phi_ed;
% 
%             SIM.OutputMatrix(P.OM.cell_volt,idx_phi_ed_AN) = -1;
%             SIM.OutputMatrix(P.OM.cell_volt,idx_phi_ed_CA) =  1;
%         % @AN/SEP
%             i = N.N_CV_AN(end);
%             index_offset = (i-1)*N.N_SV_AN;
%         % Delta Phi   @AN/SEP
%             idx = index_offset + P.del_phi;
%             SIM.OutputMatrix(P.OM.del_phi,idx) =  1;
%         % Temperature @AN/SEP
%             idx = index_offset + P.T;
%             SIM.OutputMatrix(P.OM.temp,idx) = 1;
%         % C_Liion     @AN/SEP
%             idx = index_offset + P.C_Liion;
%             SIM.OutputMatrix(P.OM.C_Liion,idx) = 1;
%         % X_surf      @AN/SEP
%             idx = index_offset + P.C_Li_surf_AN;
%             SIM.OutputMatrix(P.OM.X_surf,idx) = 1/AN.C_Li_max;
%         % i_Far      @AN/SEP
%             idx = index_offset + P.i_PS;
%             SIM.OutputMatrix(P.OM.i_Far,idx) = 1;
%         % Eta      @AN/SEP
%             idx = index_offset + P.V_2;
%             SIM.OutputMatrix(P.OM.eta,idx) = 1;
%             idx = index_offset + P.V_1;
%             SIM.OutputMatrix(P.OM.eta,idx) = -1;
% 
% i = 1;
% P.SS.omega    = i; i = i + 1;
% P.SS.Z_mag    = i; i = i + 1;
% P.SS.Z_Re     = i; i = i + 1;
% P.SS.Z_Im     = i; i = i + 1;
% P.SS.Z_dB     = i; i = i + 1;
% P.SS.Z_ps_deg = i; i = i + 1;
% 
% [A,B,C,D,~] = getSSImpedance(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS);
% 
% % save_filename = [num2str(SIM.SOC_start) 'SOC.mat'];
% % save_filename = ['80SOC_SimEndState.mat'];
% save_filename = ['LongTau'];
% multiple = -SIM.A_c^-1;
% sys = multiple*ss(A,B,C,D,'E',SIM.M);
% postProcessComplete = 1;
% clear SIM
% SIM.SimMode = 0;
% save(save_filename, 'sys','postProcessComplete','SIM')


%% Change Parameters of an existing sim
clear all; close all; clc;
%Thinking the simulations are struggling to make it past the knee of the
%voltage curve so I'm thinking I'll adjust the final time
filename = 'F:\TylerFiles\GitHubRepos\p2d-model\BatteryModel\BatchMode_DAE\Results\Carderock_C_dl_ANx3_C_dl_CAx3\C_dl_ANx3_C_dl_CAx3_Polar_0.33C_D.mat';
% filename = 'F:\TylerFiles\GitHubRepos\p2d-model\BatteryModel\BatchMode_DAE\Results\Carderock_C_dl_ANx3\C_dl_ANx3_Polar_0.33C_D.mat';
% filename = 'F:\TylerFiles\GitHubRepos\p2d-model\BatteryModel\BatchMode_DAE\Results\Carderock_C_dl_CAx3\C_dl_CAx3_Polar_0.33C_D.mat';
load(filename)
% t_final = 8100;
% SIM.tspan = [0, t_final];



save(filename,'AN','CA','SEP','EL','SIM','CONS','P','N','FLAG')



%% Testing Cap Calc
%%%%%%%%!!!!!!!!!!!!!!!! May just get rid of this
% OneC = SIM.Cell_Cap;
% C2   = OneC/2;
% C3   = OneC/3;
% C20  = OneC/20;
% 
% idx_1C   = find(I_user == -OneC,1,'first');
% idx_C2   = find(I_user == -C2,1,'first');
% idx_C3   = find(abs(I_user - (-1*C3))<1e-3,1,'first');
% idx_C20  = find(abs(I_user - (-1*C20))<1e-3,1,'first');
% idx_Rest = find(abs(I_user - (-1*C20))<1e-3,1,'last') + 1;
% 
% t_1C  = t_soln(idx_C2)   - t_soln(idx_1C);
% t_C2  = t_soln(idx_C3)   - t_soln(idx_C2);
% t_C3  = t_soln(idx_C20)  - t_soln(idx_C3);
% t_C20 = t_soln(idx_Rest) - t_soln(idx_C20);
% 
% Cap_1C  = OneC*t_1C /3600;
% Cap_C2  = C2  *t_C2 /3600;
% Cap_C3  = C3  *t_C3 /3600;
% Cap_C20 = C20 *t_C20/3600;
% Cap_tot = Cap_1C + Cap_C2 + Cap_C3 + Cap_C20;
% 
% SOC_1C  = 100*Cap_1C  / SIM.Cell_Cap;
% SOC_C2  = 100*Cap_C2  / SIM.Cell_Cap;
% SOC_C3  = 100*Cap_C3  / SIM.Cell_Cap;
% SOC_C20 = 100*Cap_C20 / SIM.Cell_Cap;
% SOC_tot = SOC_1C + SOC_C2 + SOC_C3 + SOC_C20;
% 
% SOC_error = SIM.SOC_start - (SOC(end) - SOC_tot)
