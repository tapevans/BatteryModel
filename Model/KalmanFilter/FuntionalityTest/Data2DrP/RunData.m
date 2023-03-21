%% Run Data
% Data loaded from file
%   time:    1xN vector with discrete time in seconds
%   u:       1xN vector with input signal in A/m^2
%   outputs: 6xN matrix with the desired outputs
%   P:       struct with Pointers that can be used to index the output rows
%   Lables:  struct used for labeling the title and y-axis in plots

clear all; close all; clc;

%% FLAGS
FLAG.PRBSAmp = 1;  % Amplitude of Current
FLAG.SOC     = 50; % State of Charge
FLAG.Tswitch = 100; % PRBS switching time
FLAG.Ts      = 1;  % Sampling Rate
FLAG.PlotData= 1;  % 1 if plot data, 0 if not
FLAG.UseRelax = 1;


%% Import Data
if FLAG.UseRelax 
    filename = ['PRBS_Relax_Amp' num2str(FLAG.PRBSAmp) '_SOC' num2str(FLAG.SOC) '_SwitchingTime' num2str(FLAG.Tswitch) '_Ts' num2str(FLAG.Ts) '.mat'];
else
    filename = ['PRBS_Amp' num2str(FLAG.PRBSAmp) '_SOC' num2str(FLAG.SOC) '_SwitchingTime' num2str(FLAG.Tswitch) '_Ts' num2str(FLAG.Ts) '.mat'];
end
load(filename)

if FLAG.PlotData
    %% Plot Current
    figure
    plot(time,u)
    xlabel('Time [s]')
    ylabel('Current [A/m^2]')


    %% Plot Desired Outputs
    for i = 1:length(fieldnames(P))
        figure
        plot(time,outputs(i,:))
        xlabel('Time [s]')
        ylabel(Labels.unit{i})
        title(Labels.title{i})
    end
end