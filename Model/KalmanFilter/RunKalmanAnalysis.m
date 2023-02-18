%% Run Kalman Analysis for P2D
clear all; close all; clc


%% Analysis to Perform
    FLAG.Analysis.NoNoiseCompare        = 1;
        FLAG.Analysis.ode           = 1;
        FLAG.Analysis.SS_CT         = 0;
        FLAG.Analysis.SS_DT         = 0;
        FLAG.Analysis.ROM_Mlab      = 0;
        FLAG.Analysis.ROM_HoKal     = 0;
            FLAG.Analysis.PlotImp = 0;
    FLAG.Analysis.NoisyPlant            = 0;
    FLAG.Analysis.Estimator             = 0;
    FLAG.Analysis.Est_Error_calc        = 0;
        FLAG.Analysis.dispResults = 1;
    FLAG.Analysis.GenComparData         = 0;
    FLAG.Analysis.ComparSVD2Pinf        = 0;
        % FLAG.CompareQMode
        %  1) Input Q
        %  2) State Q
        FLAG.CompareQMode = 1;
        FLAG.READ_IN_DATA = 1;

    
    SIM.Analysis.OverwriteData = 1; %% If you try to run the same sim with different Analysis, it loads old flags and won't do new analysis
    SIM.Analysis.NoisyPlant.OverwriteData = 1; %% If changing t_final, need to rerun Slink
    
    FLAG.PLOT.PlotResults = 0;
    
    FLAG.TestCase = 1; % Filename just becomes 'TestCase.mat'. Will still overwrite filenames for Slink and data generation

    SIM.AnalysisComplete.Initialization = 0;


%% Conditions to Run
Q_0_vec = [1e-6];
R_0_vec = [1e-6];
Ts_vec  = [1];
SOC_vec = [50];

N.N_samples = 600;


%% Input Type
    % Inputs
    %  1) Step
    %  2) Harmonic
    %  3) Ramp ~ subset of step
    %  4) Impulse ~ Not to be used as a simulation
    FLAG.InputMode = 1;


%% Call Analysis Main






















