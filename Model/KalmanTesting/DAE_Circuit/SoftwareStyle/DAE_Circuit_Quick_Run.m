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

%% Test 3DT SS by comparing step to ode
    % [y_DT3,t_DT3,x_DT3] = step(sys_DT3,SIM.t_final_sim);
    FLAG.InputType = 3; %No Ramp
    [InputSignal] = getInputSignal(SIM,N,P,FLAG);
    % Fix for actual sampling time
    InputSignal = InputSignal((1:10:length(InputSignal)) , :);
    y_0 = SIM.x_0_3;
    x_red = sys_CT3.C\y_0;
    
    [y_DT3,t_DT3,x_DT3] = lsim(sys_DT3,InputSignal(:,2),InputSignal(:,1),x_red);


%% Test 5CT SS by comparing step to ode
    [y_CT5,t_CT5,x_CT5] = step(sys_CT5,SIM.t_final_sim);
     % Only work with zero IC

%% Test 5DT SS by comparing step to ode
    [y_DT5,t_DT5,x_DT5] = step(sys_DT5,SIM.t_final_sim);
     % Only work with zero IC     

%% Test Run Step ode 3
    FLAG.StateMode = 3;
    FLAG.InputType = 3;
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


%% Test Simulink
% Simulink Models to run
% 1) 3 state plant
% 2) 3 state plant with input noise
% 3) 3 state plant with state noise
FLAG.SlinkModel = 2;
switch FLAG.SlinkModel
    case 1 % 3 state plant
        [t_Slink_3Plant,x_Slink_3Plant,z_Slink_3Plant] = RunSimulink(SIM,N,P,FLAG);
    case 2 % 3 state plant with input noise
        [t_Slink_3Plant,x_Slink_3Plant,z_Slink_3Plant] = RunSimulink(SIM,N,P,FLAG);%_IN
    case 3 % 3 state plant with state noise
        [t_Slink_3Plant,x_Slink_3Plant,z_Slink_3Plant] = RunSimulink(SIM,N,P,FLAG);%_SN
end

%% Compare Plots
    % figure
    % hold on
    % plot(t_soln5,x_soln5(:,P.V_1),'LineWidth',2,'DisplayName','V_1')
    % plot(t_soln5,x_soln5(:,P.V_2),'LineWidth',2,'DisplayName','V_2')
    % plot(t_soln5,x_soln5(:,P.V_3),'LineWidth',2,'DisplayName','V_3')
    % plot(t_soln3+2*SIM.Ts,x_soln3(:,P.V_1),'ko','LineWidth',2,'DisplayName','V_1')
    % plot(t_soln3+2*SIM.Ts,x_soln3(:,P.V_2),'ko','LineWidth',2,'DisplayName','V_2')
    % plot(t_soln3+2*SIM.Ts,x_soln3(:,P.V_3),'ko','LineWidth',2,'DisplayName','V_3')

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
%     figure
%     hold on
%     plot(t_DT3,y_DT3(:,P.V_1),'o','LineWidth',2,'DisplayName','V_1')
%     plot(t_DT3,y_DT3(:,P.V_2),'o','LineWidth',2,'DisplayName','V_2')
%     plot(t_DT3,y_DT3(:,P.V_3),'o','LineWidth',2,'DisplayName','V_3')
%     plot(t_soln3,x_soln3(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
%     plot(t_soln3,x_soln3(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
%     plot(t_soln3,x_soln3(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')
%     title('Compare DT3')
%     
%     figure
%     hold on
%     plot(t_DT5,y_DT5(:,P.V_1),'o','LineWidth',2,'DisplayName','V_1')
%     plot(t_DT5,y_DT5(:,P.V_2),'o','LineWidth',2,'DisplayName','V_2')
%     plot(t_DT5,y_DT5(:,P.V_3),'o','LineWidth',2,'DisplayName','V_3')
%     plot(t_soln3,x_soln3(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
%     plot(t_soln3,x_soln3(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
%     plot(t_soln3,x_soln3(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')
%     title('Compare DT5')

figure
hold on
plot(x_Slink_3Plant.time,x_Slink_3Plant.value(P.V_1,:),'o','LineWidth',2,'DisplayName','V_1')
plot(x_Slink_3Plant.time,x_Slink_3Plant.value(P.V_2,:),'o','LineWidth',2,'DisplayName','V_2')
plot(x_Slink_3Plant.time,x_Slink_3Plant.value(P.V_3,:),'o','LineWidth',2,'DisplayName','V_3')
plot(t_soln5,x_soln5(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
plot(t_soln5,x_soln5(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
plot(t_soln5,x_soln5(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')


%% Test  Run Sine 3

%% Test  Run Sine 5





%% Test Save File

%% Test Plot
