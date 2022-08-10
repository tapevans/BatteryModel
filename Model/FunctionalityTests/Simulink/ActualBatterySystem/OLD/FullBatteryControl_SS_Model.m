%% Test Battery Controller
clear all; close all; clc;

%% Flags
loc_FLAG.Proportional = 1;
loc_FLAG.LowPass      = 1;

%% Load File Parameters
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\FunctionalityTests\Simulink\ActualBatterySystem\Semi_Explicit_Test_KPCont_Profile_CV_Test_1SmallStepSOC81.93.mat';
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\Semi_Explicit_Test\Semi_Explicit_Test_KPCont_Profile_CV_Test_1SmallStep_IC4.0VSOC100.mat';
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\FunctionalityTests\Simulink\ActualBatterySystem\Semi_Explicit_Test_SS_EIS_SOC81.93.mat';
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\FunctionalityTests\Simulink\ActualBatterySystem\MassIdentity_Test_SS_EIS_SOC80.46.mat';
filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\FunctionalityTests\Simulink\ActualBatterySystem\MassIdentity_Test_SS_EIS_SOC50.mat';
load(filename);

%% Battery System
multiple = -SIM.A_c^-1;
% sys = multiple*ss(A,B,C,D,'E',SIM.M);
sys = multiple*ss(A,B,C,D);
% syse = ss(sys,'explicit');
% tf_sys = tf(sys);
% num_sys = cell2mat(tf_sys.Numerator(1));
% den_sys = cell2mat(tf_sys.Denominator(1));

% bodeplot(sys,syse,'g--')
% figure
% step(sys)

%% Controller Parameters
controller_tf = tf(1,1);
if loc_FLAG.Proportional
%     K_p = 1;
    K_p = 275;
    controller_tf = controller_tf * tf(K_p,1);
else
    K_p = 0; % This is here for filenaming
end
if loc_FLAG.LowPass
%     cutoff_fq = 7e9;
    cutoff_fq = 7e1;
    controller_tf = controller_tf * tf([cutoff_fq],[1 , cutoff_fq]);
else
    cutoff_fq = 0;
end
[num,den] = tfdata(controller_tf);
num = cell2mat(num);
den = cell2mat(den);

%% Controller Saturation
CrateSat = 2;
CurrentLimit = CrateSat*SIM.Cell_Cap;
UpperLimit =  CurrentLimit;
LowerLimit = -CurrentLimit;

%% Convert Function Handles into Strings
AN.EqPotentialHandle = func2str(AN.EqPotentialHandle);
AN.i_oHandle         = func2str(AN.i_oHandle);
AN.sigmaHandle       = func2str(AN.sigmaHandle);
AN.D_oHandle         = func2str(AN.D_oHandle);

CA.EqPotentialHandle = func2str(CA.EqPotentialHandle);
CA.i_oHandle         = func2str(CA.i_oHandle);
CA.sigmaHandle       = func2str(CA.sigmaHandle);
CA.D_oHandle         = func2str(CA.D_oHandle);

EL.tf_numHandle      = func2str(EL.tf_numHandle);
EL.ActivityHandle    = func2str(EL.ActivityHandle);
EL.D_o_Li_ionHandle  = func2str(EL.D_o_Li_ionHandle);
EL.kappaHandle       = func2str(EL.kappaHandle);

SIM = rmfield(SIM,'fsolve_options');
% SIM = rmfield(SIM,'ControllerHandle');
%% Simulation Parameters
t_final = 0.015; % [s]

%% Initial Conditions
SS_IC = SIM.SV_IC_Static;

%% Run Simulink
mdl = 'Simple_Controller_SS_Model';
set_param('Simple_Controller_SS_Model','AlgebraicLoopSolver','Auto')
in = Simulink.SimulationInput(mdl);
in = in.setModelParameter('StartTime','0','StopTime',num2str(t_final));
% activeConfigObj = getActiveConfigSet(mdl);
% set_param(activeConfigObj,'RelTol','1e-6')
% activeConfigObj
out = sim(in);

%% Save Results
save_filename = ['F:\TylerFiles\GitHubRepos\BatteryModel\Model\FunctionalityTests\Simulink\ActualBatterySystem\Kp' num2str(K_p) 'omega_o' num2str(cutoff_fq) '.mat'];


%% Plot Results
sub_off = 0;
% sub_off = 0.04;
% Time vs Error
figure
plot(out.tout,out.Error,'LineWidth',2)
title('Error')
xlabel('Time (s)')
ylabel('Error (V)')
xlim([0,out.tout(end)-sub_off])

% Time vs Voltage
figure
plot(out.tout,out.cell_voltage,'LineWidth',2)
title('Voltage')
xlabel('Time (s)')
ylabel('Voltage (V)')
xlim([0,out.tout(end)-sub_off])

% Time vs Current
figure
plot(out.tout,out.i_user,'LineWidth',2)
title('Current')
xlabel('Time (s)')
ylabel('Current (A)')
xlim([0,out.tout(end)-sub_off])

% Time vs Voltage and Current
figure
title('Cell Voltage and Current vs Time')
xlabel('Time (s)')

yyaxis left
plot(out.tout , out.cell_voltage , 'Linewidth' , 2 )
ylabel('Cell Voltage (V)')
xlim([0,out.tout(end)-sub_off])

yyaxis right
plot(out.tout , out.i_user, 'Linewidth' , 2 )
ylabel('Current (A)')

%% Arrange Figures
FigArrange = 1;
fig = gcf;
NumFig = fig.Number;
if FigArrange == 1
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