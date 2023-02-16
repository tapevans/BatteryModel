%% Run DAE Circuit Project
clear all; close all; clc;


%% Subdirecties to Include
% Script's filepath
[current_file_path,~,~] = fileparts(mfilename('fullpath'));

% Include all folders
addpath(genpath(current_file_path));
% genpath creates a string with all folders and subfolders in a given directory
% addpath then adds all of them to the current workspace


%% Analysis to Perform
FLAG.Analysis.NoNoiseCompare    = 1;
FLAG.Analysis.ode3      = 1;
FLAG.Analysis.ode5      = 1;
FLAG.Analysis.SS_CT3    = 1;
FLAG.Analysis.SS_DT3    = 1;
FLAG.Analysis.ROM_Mlab3 = 1;
FLAG.Analysis.ROM_HoKal = 1;
FLAG.Analysis.NoisyPlant        = 1;
FLAG.Analysis.Estimator         = 1;
FLAG.Analysis.Est_Error_calc    = 1;
    FLAG.Analysis.dispResults = 1;
FLAG.Analysis.GenComparData     = 0;
FLAG.Analysis.ComparSVD2Pinf    = 0;

SIM.AnalysisComplete.Initialization             = 0;
SIM.AnalysisComplete.NoNoiseCompare.ode3        = 0;
SIM.AnalysisComplete.NoNoiseCompare.ode5        = 0;
SIM.AnalysisComplete.NoNoiseCompare.SS_CT3      = 0;
SIM.AnalysisComplete.NoNoiseCompare.SS_DT3      = 0;
SIM.AnalysisComplete.NoNoiseCompare.ROM_Mlab3   = 0;
SIM.AnalysisComplete.NoNoiseCompare.ROM_HoKal   = 0;
SIM.AnalysisComplete.NoisyPlant                 = 0;
SIM.AnalysisComplete.Estimator                  = 0;
SIM.AnalysisComplete.Est_Error_calc             = 0;
SIM.AnalysisComplete.GenComparData              = 0;
SIM.AnalysisComplete.ComparSVD2Pinf             = 0;

SIM.Analysis.OverwriteData = 1;
SIM.Analysis.NoisyPlant.OverwriteData = 1;

FLAG.PLOT.PlotResults = 1;

FLAG.TestCase = 1;

%% Subdirecties to Include
% Script's filepath
[current_file_path,~,~] = fileparts(mfilename('fullpath'));

% Include all folders
addpath(genpath(current_file_path));
% genpath creates a string with all folders and subfolders in a given directory
% addpath then adds all of them to the current workspace


%% State to Measure
%  1) V_1
%  2) V_2
%  3) V_3
%  4) All ~ Use all for comparing no noise simulations mostly
%%%FLAG.C_mode = 3;

%!!! C_vec replaces FLAG.C_mode
%     C_vec = [1 2 3 ]; 
    C_vec = [3]; 

%% Sampling Rate
%%%SIM.Ts = 1e-6; % Sampling Rate (When to save outputs)

%!!! Ts_vec replaces SIM.Ts
    %Ts_vec = [1e-7 5e-7 1e-6 5e-6 1e-5 1e-4];
    Ts_vec = [5e-6 ];


%% Noise
% Q Modes
%  1) Input Q
%  2) State Q
FLAG.QMode = 1;

% Noise Properties
% SIM.tc  = SIM.Ts; % [s],Correlation Time (Noise Sampling Rate)
SIM.Q_0 = 1e-1; % Process Noise
SIM.R_0 = 1e-1; % Measurement Noise

%!!! Q_0_vec and R_0_vec replaces SIM.Q_0 and SIM.R_0
    %Q_0_vec = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1];
    %R_0_vec = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1];

    Q_0_vec = [1e-3 ];
    R_0_vec = [1e-3 ];


%% Input Type
% Inputs
%  1) Step
%  2) Harmonic
%  3) Ramp ~ subset of step
%  4) Impulse ~ Not to be used as a simulation
FLAG.InputMode = 1;
SIM.V_init        = 0;
SIM.V_step_height = 4;


%% Simulation t_final
%%%% These are for Step Simulations %%%%%%
%     SIM.t_final = 6e-6;
%     SIM.t_final = 1e-4;
%     SIM.t_final = 6e-4;
    SIM.t_final = 5e-3;
%     SIM.t_final = 5*max(SIM.C_1*SIM.R_1,SIM.C_4*SIM.R_4);

%%%% These are for Harmonic Simulations %%%%%%
%     N_cycles = 10; %% WAAAAAAAY too long
%     SIM.t_final = N_cycles/SIM.fq;

if FLAG.InputMode == 2
    SIM.V_init        = 0;
    SIM.t_final = 0.005; % (5000 * 1e-6) USE FOR SINE, 5000 data points sampled at Ts = 1e-6
    % SIM.fq = SIM.Ts/15; %% WAAAAAAAY too long
    SIM.fq = 400; %(2 cycles over t_final = 5000*Ts, had to calculate manually)
    SIM.V_s = 1; %Amplitude of sine wave
end


%% Mode for Simulink Model
% Simulnk Model to run
%  3) 3 States
%  5) 5 States !!!!!!!!! Not implemented
FLAG.SlinkModel = 3;


%% Estimator Model
% Model for Estimator to use
%  1) Matlab SS_DT3
%  2) Matlab SS_DT5 !!!!! Not implemented, also, is this able to mix with 3 state Slink plant?
%  3) Ho-Kalman
FLAG.EstimatorModel = 3;


%% Call Analysis Main
% Add for loop later for data generation
for QQ = Q_0_vec
    for RR = R_0_vec
        for TT = Ts_vec
            for CC = C_vec
                SIM.Q_0 = QQ; % Process Noise
                SIM.R_0 = RR; % Measurement Noise
                SIM.Ts = TT; % Sampling Rate (When to save outputs)
                FLAG.C_mode = CC;

                SIM.tc  = SIM.Ts; % [s],Correlation Time (Noise Sampling Rate)
                
                %% Make Save Filename
                Q_str = num2str(SIM.Q_0);
                R_str = num2str(SIM.R_0);
                Ts_str = num2str(SIM.Ts);
                C_str = num2str(FLAG.C_mode);
                switch FLAG.QMode
                    case 1 % Input Q
                        Noise_str = '_InputQ';
                    case 2 % State Q
                        Noise_str = '_StateQ';
                end
                switch FLAG.InputMode
                    case 1 % Step
                        Input_str = '_Step';
                    case 2 % Harmonic
                        Input_str = '_Sine';
                        fq_str    = [num2str(SIM.fq) 'Hz'];
                        Input_str = [Input_str fq_str];
                end
                
                if FLAG.TestCase
                    filename = 'TestCase.mat';
                else
                    filename = ['Plant' Input_str Noise_str '_Q' Q_str '_R' R_str '_Ts' Ts_str '_C' C_str  '.mat'];
                end
                filepath = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanTesting\DAE_Circuit\SoftwareStyle\Results';
                SIM.save_filename = [filepath filesep filename];
                
                %% Run Analysis
                AnalysisMain(SIM,FLAG)
            end
        end
    end
end


