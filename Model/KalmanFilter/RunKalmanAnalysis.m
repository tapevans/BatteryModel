%% Run Kalman Analysis for P2D
clear all; close all; clc

%% Subdirecties to Include
    % Script's filepath
    [current_file_path,~,~] = fileparts(mfilename('fullpath'));

    % Ensure this is the current directory, important for pwd
    cd(current_file_path)
    
    % Include all folders
    addpath(genpath(current_file_path));
    % genpath creates a string with all folders and subfolders in a given directory
    % addpath then adds all of them to the current workspace


%% Analysis to Perform
    FLAG.Analysis.NoNoiseCompare        = 1;
        FLAG.Analysis.ode           = 0;
        FLAG.Analysis.SS_CT         = 0;
        FLAG.Analysis.SS_DT         = 0;
        FLAG.Analysis.OptimalHK     = 1;
            FLAG.r_max = 50;
            FLAG.UseInput_r = 1;
        FLAG.Analysis.ROM_HoKal     = 0;
            FLAG.Analysis.PlotImp = 0;
            FLAG.EST.SepHK        = 0; % 1 if calculate a ROM for each desired variable
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

    
    FLAG.OverwriteData.All      = 0; %% If you try to run the same sim with different Analysis, it loads old flags and won't do new analysis
    FLAG.OverwriteData.Slink    = 0; %% If changing t_final, need to rerun Slink
    FLAG.OverwriteData.GenData  = 0;
    FLAG.OverwriteData.ROMError = 1;
    FLAG.OverwriteData.ROM      = 0; % Don't think I have used this yet
    
    FLAG.PLOT.PlotResults = 0; 
    
    FLAG.TestCase = 1; % Don't think I have used this yet % Filename just becomes 'TestCase.mat'. Will still overwrite filenames for Slink and data generation

%     FLAG.AnalysisComplete.Initialization = 0;


%% Desired Outputs
    FLAG.DesOut.cell_voltage  = 1; % Terminal Voltage
    FLAG.DesOut.delta_phi     = 1; % Electrostatic potential difference between active material and electrolyte @AN/SEP interface
    FLAG.DesOut.i_Far         = 0; % Chemical reaction current density at SEI @AN/SEP interface
    FLAG.DesOut.eta           = 1; % Overpotential at SEI @AN/SEP interface
    FLAG.DesOut.C_Liion       = 1; % Concentration of Li^+ in the electrolyte @AN/SEP interface
    FLAG.DesOut.C_Li          = 1; % Concentration of Li at the surface of active material @AN/SEP interface
    FLAG.DesOut.delta_C_Li    = 1; % Difference of Li concentration between the surface and center of active material particle @AN/SEP interface
    FLAG.DesOut.T             = 0; % Temperature of control volume @AN/SEP interface


%% Filepath to get simulation results
    FLAG.folderpath         = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\ObservabilityTest';
    FLAG.folderpathROMError = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanFilter\Results\ROMErrorData';
    FLAG.folderpathPRBS     = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\PRBS_Sims';

%% Noise
    % Q Modes
    %  1) Input Q
    %  2) State Q ####### Not implemented in Slink
    FLAG.QMode = 1;


%% Estimator Model
    % Model for Estimator to use
    %  1) Matlab SS_DT
    %  2) Ho-Kalman
    FLAG.EstimatorModel = 2;
        
    

%% Conditions to Run
    Q_0_vec = [1e-3];
    R_0_vec = [1e-3];
%     R_0_vec = [1e1];

    % Desired Sampling Times
        N_t_s   = 25;   % **Keep this fix for now
        T_s_min = -1; % T_s = 10^(T_s_min), **Keep this fix for now
        T_s_max =  1; % T_s = 10^(T_s_max), **Keep this fix for now
        Ts_vec = logspace(T_s_min,T_s_max,N_t_s);
%         Ts_vec  = [10];

    SOC_vec = 0:1:100;    
%     SOC_vec = [50];

    FLAG.N_samples = 600;


%% Input Type
    % Inputs
    %  1) Step
    %  2) Harmonic ########Not implemented
    %  3) Ramp ~ subset of step ##### Not implemented since I get step data from my P2D results
    %  4) Impulse ~ Not to be used as a simulation ## Not implemented 
    %  5) PRBS
    FLAG.InputMode = 1;

%%%%%%%%%%%%% 
% What's available:
% PRBSAmp: 1
% Tswitch: 0.1 , 1 , 10 , 100 
% SOC:     25 , 50 , 75
%(Can't do Tswitch=0.1 when doing Ho-Kalman since I don't have Ts=0.01 Impulse Response, unless SamplesPerSwitch == 1)
%(Can't do Tswitch=100 when doing Ho-Kalman and using SamplesPerSwitch == 1 since I don't have Ts=100 Impulse Response)
%!!!(Can fix these issues if I implement a DT Impulse simulation instead of calling existing simulations)
if FLAG.InputMode == 5
    FLAG.PRBSAmp = 1; 
    FLAG.Tswitch = 100; 
    %FLAG.SamplesPerSwitch = 10;
    %FLAG.N_samples
end


%% Offset IC for ROM
    FLAG.offsetROM_IC_Rel = 1e-1; % [%], Relative offset for the initial conditions for ROM Estimator
    FLAG.offsetROM_IC_Abs = 1e-6; % [%], Relative offset for the initial conditions for ROM Estimator

%% Call Analysis Main
    for QQ = Q_0_vec
        for RR = R_0_vec
            for TT = Ts_vec
                for SS = SOC_vec
                    FLAG.Q_0 = QQ; % Process Noise
                    FLAG.R_0 = RR; % Measurement Noise
                    FLAG.Ts = TT;  % Sampling Rate (When to save outputs)
                    FLAG.SOC = SS; % State of Charge
                    disp(['Q: ' num2str(FLAG.Q_0) ' R: ' num2str(FLAG.R_0) ' Ts: ' num2str(FLAG.Ts) ' SOC: ' num2str(FLAG.SOC) ] )

                    AnalysisMain(FLAG)
                end
            end
        end
    end


%% Generate Comparison Data Call Main
%     Q_0 = 1e-6; R_0 = 1e-6; % TL % Completed
%     Q_0 = 1e-6; R_0 = 1e1;  % BL
%     Q_0 = 1e1;  R_0 = 1e1;  % BR
%     Q_0 = 1e1;  R_0 = 1e-6; % TR
%     Q_0 = 1e-3; R_0 = 1e-3; % Middle

%     FLAG.Analysis.NoisyPlant            = 0;
%     FLAG.Analysis.GenComparData         = 1;
%     
%     Q_vec = [1e-6 1e1 1e1 1e-3];
%     R_vec = [1e1 1e1 1e-6 1e-3];
%     Noise_Vec = [Q_vec;R_vec];
%     Ts_vec = logspace(T_s_min,T_s_max,N_t_s);
%     SOC_vec = 0:1:100; 
% 
%         for i = 1:length(Noise_Vec)
%             for TT = Ts_vec
%                 for SS = SOC_vec
%                     FLAG.Q_0 = Noise_Vec(1,i); % Process Noise
%                     FLAG.R_0 = Noise_Vec(2,i); % Measurement Noise
%                     FLAG.Ts = TT;  % Sampling Rate (When to save outputs)
%                     FLAG.SOC = SS; % State of Charge
%                 disp(['Q: ' num2str(FLAG.Q_0) ' R: ' num2str(FLAG.R_0) ' Ts: ' num2str(FLAG.Ts) ' SOC: ' num2str(FLAG.SOC) ] )
% 
%                     AnalysisMain(FLAG)
%                 end
%             end
%         end



%% Generate Slink Data Call Main
%     FLAG.Analysis.NoisyPlant            = 1;
%     FLAG.Analysis.GenComparData         = 0;
%     
%     Q_vec = [1e1  1e1 1e1];
%     R_vec = [1e1 1e-3 1e-6];
%     Noise_Vec = [Q_vec;R_vec];
%     [r,c] = size(Noise_Vec);
%     Ts_vec = logspace(T_s_min,T_s_max,N_t_s);
%     Ts_vec = Ts_vec([25]); % [1,7,13, 20,25]
%     SOC_vec = [25]; 
% 
%     for i = 1:c
%         for TT = Ts_vec
%             for SS = SOC_vec
%                 FLAG.Q_0 = Noise_Vec(1,i); % Process Noise
%                 FLAG.R_0 = Noise_Vec(2,i); % Measurement Noise
%                 FLAG.Ts = TT;  % Sampling Rate (When to save outputs)
%                 FLAG.SOC = SS; % State of Charge
%                 disp(['Q: ' num2str(FLAG.Q_0) ' R: ' num2str(FLAG.R_0) ' Ts: ' num2str(FLAG.Ts) ' SOC: ' num2str(FLAG.SOC) ] )
% 
%                 AnalysisMain(FLAG)
%             end
%         end
%     end















