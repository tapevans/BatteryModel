%% DAE Circuit Quick Run Function
clear all; close all; clc;

%% Subdirecties to Include
% Script's filepath
[current_file_path,~,~] = fileparts(mfilename('fullpath'));

% Include all folders
addpath(genpath(current_file_path)); 
% genpath creates a string with all folders and subfolders in a given directory
% addpath then adds all of them to the current workspace

%% Test Inputs
DAE_Circuit_CreateProject
%
[SIM,N,P,FLAG] = init(SIM,FLAG);


%% Test Run Step ode 3
    FLAG.StateMode = 3;
    FLAG.InputType = 3; % 1) Ramp, %2) Sine, 3) Step
    [t_soln3, x_soln3, z_soln3] = runODE(SIM,N,P,FLAG);
    
    % figure
    % hold on
    % plot(t_soln3,x_soln3(:,P.V_1),'LineWidth',2,'DisplayName','V_1')
    % plot(t_soln3,x_soln3(:,P.V_2),'LineWidth',2,'DisplayName','V_2')
    % plot(t_soln3,x_soln3(:,P.V_3),'LineWidth',2,'DisplayName','V_3')

%% Test Run Step ode 5
    FLAG.StateMode = 5;
    FLAG.InputType = 1;
    [t_soln5, x_soln5, z_soln5] = runODE(SIM,N,P,FLAG);
    
    % figure
    % hold on
    % plot(t_soln5,x_soln5(:,P.V_1),'LineWidth',2,'DisplayName','V_1')
    % plot(t_soln5,x_soln5(:,P.V_2),'LineWidth',2,'DisplayName','V_2')
    % plot(t_soln5,x_soln5(:,P.V_3),'LineWidth',2,'DisplayName','V_3')

    % Compare ode5 to ode3
    % figure
    % hold on
    % plot(t_soln5,x_soln5(:,P.V_1),'LineWidth',2,'DisplayName','V_1')
    % plot(t_soln5,x_soln5(:,P.V_2),'LineWidth',2,'DisplayName','V_2')
    % plot(t_soln5,x_soln5(:,P.V_3),'LineWidth',2,'DisplayName','V_3')
    % plot(t_soln3+2*SIM.Ts,x_soln3(:,P.V_1),'ko','LineWidth',2,'DisplayName','V_1')
    % plot(t_soln3+2*SIM.Ts,x_soln3(:,P.V_2),'ko','LineWidth',2,'DisplayName','V_2')
    % plot(t_soln3+2*SIM.Ts,x_soln3(:,P.V_3),'ko','LineWidth',2,'DisplayName','V_3')


%% Test getting SS systems
FLAG.C_mode = 4;
[sys_CT3 , sys_DT3 , sys_CT5 , sys_DT5] = getSS_System(SIM,N,P,FLAG);

%% Test 3CT SS by comparing step to ode
FLAG.InputType = 3; %No Ramp
[InputSignal] = getInputSignal(SIM,N,P,FLAG);
% creating the ss system in matlab, the states are reduced, need to fix the IC
% This method hasn't always worked
y_0 = SIM.x_0_3;
x_red = sys_CT3.C\y_0;

[y_CT3,t_CT3,x_CT3] = lsim(sys_CT3,InputSignal(:,2),InputSignal(:,1),x_red);

% [y_CT3,t_CT3,x_CT3] = step(sys_CT3,SIM.t_final_sim); Doesn't do IC which
% is needed for 3 state system

%     figure
%     hold on
%     plot(t_CT3,y_CT3(:,P.V_1),'o','LineWidth',2,'DisplayName','V_1')
%     plot(t_CT3,y_CT3(:,P.V_2),'o','LineWidth',2,'DisplayName','V_2')
%     plot(t_CT3,y_CT3(:,P.V_3),'o','LineWidth',2,'DisplayName','V_3')
%     plot(t_soln3,x_soln3(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
%     plot(t_soln3,x_soln3(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
%     plot(t_soln3,x_soln3(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')
%     title('Compare CT3')
%     

%% Test 3DT SS by comparing step to ode
    % [y_DT3,t_DT3,x_DT3] = step(sys_DT3,SIM.t_final_sim);
    FLAG.InputType = 3; %No Ramp
    [InputSignal] = getInputSignal(SIM,N,P,FLAG);
    % Fix for actual sampling time
    InputSignal = InputSignal((1:10:length(InputSignal)) , :);
    y_0 = SIM.x_0_3;
    x_red = sys_CT3.C\y_0;
    
    [y_DT3,t_DT3,x_DT3] = lsim(sys_DT3,InputSignal(:,2),InputSignal(:,1),x_red);

%     figure
%     hold on
%     plot(t_DT3,y_DT3(:,P.V_1),'o','LineWidth',2,'DisplayName','V_1')
%     plot(t_DT3,y_DT3(:,P.V_2),'o','LineWidth',2,'DisplayName','V_2')
%     plot(t_DT3,y_DT3(:,P.V_3),'o','LineWidth',2,'DisplayName','V_3')
%     plot(t_soln3,x_soln3(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
%     plot(t_soln3,x_soln3(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
%     plot(t_soln3,x_soln3(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')
%     title('Compare DT3')

%% Test 5CT SS by comparing step to ode
    [y_CT5,t_CT5,x_CT5] = step(sys_CT5,SIM.t_final_sim);
     % Only work with zero IC

%     figure
%     hold on
%     plot(t_CT5,y_CT5(:,P.V_1),'o','LineWidth',2,'DisplayName','V_1')
%     plot(t_CT5,y_CT5(:,P.V_2),'o','LineWidth',2,'DisplayName','V_2')
%     plot(t_CT5,y_CT5(:,P.V_3),'o','LineWidth',2,'DisplayName','V_3')
%     plot(t_soln3,x_soln3(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
%     plot(t_soln3,x_soln3(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
%     plot(t_soln3,x_soln3(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')
%     title('Compare CT5')
%     

%% Test 5DT SS by comparing step to ode
    [y_DT5,t_DT5,x_DT5] = step(sys_DT5,SIM.t_final_sim);
     % Only work with zero IC     

%     figure
%     hold on
%     plot(t_DT5,y_DT5(:,P.V_1),'o','LineWidth',2,'DisplayName','V_1')
%     plot(t_DT5,y_DT5(:,P.V_2),'o','LineWidth',2,'DisplayName','V_2')
%     plot(t_DT5,y_DT5(:,P.V_3),'o','LineWidth',2,'DisplayName','V_3')
%     plot(t_soln3,x_soln3(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
%     plot(t_soln3,x_soln3(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
%     plot(t_soln3,x_soln3(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')
%     title('Compare DT5')


%% Test Simulink
% Simulink Models to run
% 1) 3 state plant
% 2) 3 state plant with input noise
% 3) 3 state plant with state noise
FLAG.C_mode = 4;
FLAG.SlinkModel = 2;
switch FLAG.SlinkModel
    case 1 % 3 state plant
        [t,x,z,u,w,v] = RunSimulink(SIM,N,P,FLAG);
    case 2 % 3 state plant with input noise
        [t,x,z,u,w,v] = RunSimulink(SIM,N,P,FLAG);%_IN
    case 3 % 3 state plant with state noise
        [t,x,z,u,w,v] = RunSimulink(SIM,N,P,FLAG);%_SN
end

figure
hold on
plot(x.time,x.value(P.V_1,:),'o','LineWidth',2,'DisplayName','V_1')
plot(x.time,x.value(P.V_2,:),'o','LineWidth',2,'DisplayName','V_2')
plot(x.time,x.value(P.V_3,:),'o','LineWidth',2,'DisplayName','V_3')
plot(t_soln3,x_soln3(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
plot(t_soln3,x_soln3(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
plot(t_soln3,x_soln3(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')
% plot(t_soln5,x_soln5(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
% plot(t_soln5,x_soln5(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
% plot(t_soln5,x_soln5(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')

%% Save Simulink Outputs for Estimation
switch FLAG.SlinkModel
    case 1 % 3 state plant
        type_str = 'Plant3_Step';
        Q_str = num2str(0);
        R_str = num2str(0);
    case 2 % 3 state plant with input noise
        type_str = 'Plant3_IQ_Step';
        Q_str = num2str(SIM.Q_0);
        R_str = num2str(SIM.R_0);
    case 3 % 3 state plant with state noise
        type_str = 'Plant3_SQ_Step';
        Q_str = num2str(SIM.Q_0);
        R_str = num2str(SIM.R_0);
end

Ts_str = num2str(SIM.Ts);
switch FLAG.C_mode
    case 1
        C_str = num2str(FLAG.C_mode);
    case 2
        C_str = num2str(FLAG.C_mode);
    case 3
        C_str = num2str(FLAG.C_mode);
    case 4
        C_str = num2str('All');
end


Folder = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanTesting\DAE_Circuit\SoftwareStyle\Results';
filename = [type_str '_Q' Q_str '_R' R_str '_Ts' Ts_str '_C' C_str '.mat'];
overall_filename = [Folder filesep filename];
save(overall_filename, "z", "x","u" ,"w" ,"v" , "SIM", "FLAG", "sys_DT5", "sys_DT3","sys_CT3", "sys_CT5")


%% Test Estimation
[x_hat_asy, x_hat_var, y_hat_asy, y_hat_var, P_infty, K_infty, K_k, P_k_pre] = PerformEstimation(overall_filename);

    figure
    hold on
    plot(x.time,y_hat_asy(P.V_1,:),'o','LineWidth',2,'DisplayName','V_1')
    plot(x.time,y_hat_asy(P.V_2,:),'o','LineWidth',2,'DisplayName','V_2')
    plot(x.time,y_hat_asy(P.V_3,:),'o','LineWidth',2,'DisplayName','V_3')
    plot(t_CT3,y_CT3(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
    plot(t_CT3,y_CT3(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
    plot(t_CT3,y_CT3(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')
    title('Compare Asy Estimation')

    figure
    hold on
    plot(x.time,y_hat_var(P.V_1,:),'o','LineWidth',2,'DisplayName','V_1')
    plot(x.time,y_hat_var(P.V_2,:),'o','LineWidth',2,'DisplayName','V_2')
    plot(x.time,y_hat_var(P.V_3,:),'o','LineWidth',2,'DisplayName','V_3')
    plot(t_CT3,y_CT3(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
    plot(t_CT3,y_CT3(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
    plot(t_CT3,y_CT3(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')
    title('Compare Var Estimation')

    CPCT_calc = sys_CT3.C * P_infty * sys_CT3.C';
    disp('CPCT_calc=')
    disp(num2str(CPCT_calc))
    disp(newline)
    disp('P_k_pre final value')
    disp(num2str(sys_CT3.C *P_k_pre(:,:,end)* sys_CT3.C'))


%% Test Covariance Calc
idx = 1000;
[Pcalc_var] = getP_calc(z.value,y_hat_var,idx)
[Pcalc_asy] = getP_calc(z.value,y_hat_asy,idx)
% CPcalcCT_var = sys_CT3.C * Pcalc_var * sys_CT3.C';
% CPcalcCT_asy = sys_CT3.C * Pcalc_asy * sys_CT3.C';
% disp(newline)
% disp(num2str(CPcalcCT_var))
% disp(newline)
% disp(num2str(CPcalcCT_asy))

%% Test  Run Sine 3

%% Test  Run Sine 5





%% Test Save File

%% Test Plot
