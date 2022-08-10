%% Testing Simulink 
% Tests:
% * ^Transfer function against ODE solver
% * ^String to function handle and then used inside another function
% * ^Passing a string
% * ^Passing structs
% * IC
% * Time Data for reference input
% * Controller Algorithms

%% Notes
% Mode of Operation:
% 1: CC
% 2: CV

clear all; close all; clc;

%% Inputs
SIM.R = 1;
SIM.C = 1;

% SIM.ProfileFnH = @Profile_CC;
% SIM.ProfileFnH = @Profile_CC_Steps;
% SIM.ProfileFnH = @Profile_CV;
SIM.ProfileFnH = @Profile_CV_Steps;

filename = 'Test_OverallSys.mat';

%% Create Project
% Initialize Parameters % tfinal, SV_0, 
% SIM.tfinal = (5*SIM.R*SIM.C)*3;
SIM.tfinal = 65;

% Save
save(filename,'SIM');

clearvars -except filename
%% Run Simulation
% Load File
load(filename)

% Convert Fcn Handles to string
SIM.ProfileFnH = func2str(SIM.ProfileFnH);

% Set necessary parameters for Simulink Workspace
num = [SIM.R];
den = [SIM.R*SIM.C, 1];

% Call Simulink
mdl = 'PuttingItAllTogether';
in = Simulink.SimulationInput(mdl);
in = in.setModelParameter('StartTime','0','StopTime',num2str(SIM.tfinal));
out = sim(in);

% Save Solution
save(filename,'SIM','out')

%% Post-Processing
tsoln = out.tout;
MO    = out.MO;
I_ref = out.I_ref;
V_ref = out.V_ref;
V_out = out.V_out;
I_out = out.I_out;

% Save Solution
save(filename,'SIM','out','tsoln','MO','I_ref','V_ref','V_out','I_out')

%% Plots
% t vs V_out
    figure
    plot(tsoln,V_out,'LineWidth',2)
    xlabel('Time (s)')
    ylabel('Voltage (V)')
    title('Time vs V_{out}')

% t vs V and I
    figure
    title('Voltage and Current vs Time')
    xlabel('Time (s)')

    yyaxis left
    plot(tsoln , V_out , 'Linewidth' , 2 )
    ylabel('V_{out}')

    yyaxis right
    plot(tsoln , I_out, 'Linewidth' , 2 )
    ylabel('I_{out}')

% t vs V_out and V_ref
    figure
    hold on
    plot(tsoln,V_ref,'k','LineWidth',2)
    plot(tsoln,V_out,'r','LineWidth',2)
    xlabel('Time (s)')
    ylabel('Voltage')
    title('Vout and Vref')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% OLD STUFF %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Properties 
% I_in = 2; % [A]
% R    = 1; % [Ohm]
% C    = 1; % [F]

% tfinal = 5*R*C;
% SIM.fh_str = '@getMO';

%% Initialize
% V_0 = 0; % [V], Initial Voltage

%% Make TF sys
% num = [R];
% den = [R*C, 1];
% sys = tf(num,den);

%% Run ODE Soln
% SOLN = ode45(@(t,V) myFun(t,V,R,C,I_in), [0,tfinal],V_0);

%% Run tf soln
% opt = stepDataOptions('InputOffset',0,'StepAmplitude',I_in);
% [V,tOut] = step(sys,tfinal,opt);

%% Run Simulink
% mdl = 'JustSys';
% in = Simulink.SimulationInput(mdl);
% in = in.setModelParameter('StartTime','0','StopTime',num2str(tfinal));
% out = sim(in);

% mdl = 'ModeSelectorTest';
% in = Simulink.SimulationInput(mdl);
% in = in.setModelParameter('StartTime','0','StopTime',num2str(1));
% out = sim(in);

%% Plot Soln
% figure
% hold on
% plot(SOLN.x,SOLN.y,'k','LineWidth',2)
% plot(tOut,V,'ro')
% plot(out.simout.Time,out.simout.Data,'b*')


%% Function Handle
% function dVdt = myFun(t,V,R,C,I_in)
%     dVdt = (I_in - V/R)/C;
% end
