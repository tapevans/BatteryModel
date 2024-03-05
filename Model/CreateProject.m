%% Create Project
% This file is used to create a project folder. In this folder is a file of
% battery parameters used, which simulations are to be performed, and a
% data file (.mat) for each simulation to be performed.
%
% Simulations generated with this file will use the values in batt_inputs
% for other simulation parameters such as temperature, voltage limits, how
% deep to discharge, and etc...
%
% File names
    % ---- Polarization ----
        % 'battery_name'_Polar_'Crate'C_(CorD)
%    
    % ---- Harmonic Perturbation ----
        % 'battery_name'_EIS_SIN_w'omega'_SOC'SOC'
%    
    % ---- State Space EIS ----
        % 'battery_name'_SS_EIS_SOC'SOC'
%    
    % ---- Known BC Profile Controller ----
        % 'battery_name'_KPCont_'KBCPProfileFilename'_SOC'KBSOC'
%        
    % ---- MOO Controller ----
        % 'battery_name'_MOOCont_'ControllerName'
%    
    % ---- Manual Current Profile ----
        % For platingRefinement
            % 'battery_name'_ManCurrProf_ 'MCP.N_regions' steps_ 'MCP.max_iterations' Iter_ 'MCP.tol_Delta_phi' tol 
        % For other
            % 'battery_name'_ManCurrProf_'ManualProfName'
%    
    % ---- PRBS ----
        % 'battery_name'_PRBS_Amp'PRBS_Amp'_SOC'PRBS_SOC'_SwitchingTime'PRBS_Tswitch'
%
    % ---- EIS from Stiching PRBS ---- 
        % 'battery_name'_PRBS_Amp'PRBS_Amp'_SOC'EIS_PRBS.SOC'_SwitchingTime'PRBS_Tswitch'
        % 'battery_name'_PRBS_EIS_SOC'EIS_PRBS.SOC'
%    
    % ---- EIS Ho-Kalman ----
        % 'battery_name'_EIS_HoKalman_SOC'HK_SOC'_SamplingTime'HK_Ts'
        
        
%%
% clear all; close all; clc;


%% Subdirecties to Include
    % Script's filepath
        [current_file_path,~,~] = fileparts(mfilename('fullpath'));
    
    % Include all folders
    addpath(genpath(current_file_path)); 
        % genpath creates a string with all folders and subfolders in a
        % given directory addpath then adds all of them to the current
        % workspace


%% Inputs
% InputFile
    %   -1) batt_inputs_Manual
    %    0) batt_inputs
    %    1) batt_inputs_Wiley_Lui
    %    2) batt_inputs_HZ_A
    %    3) batt_inputs_HZ_B
    %    4) batt_inputs_NREL
    %    5) batt_inputs_ToddKingston
    %    6) batt_inputs_MPC
    InputFile = 5;

    FLAG_local.folder_overwrite = 0; % 1 if delete folder if it already exists
    FLAG_local.folder_add       = 1; % 1 if just want to add simulations to folder
        % if folder_overwrite and folder_add are 0, then a new name should
        % be used for the folder
    
    FLAG_local.sim_overwrite    = 1; % 1 if older simulation is deleted and new one is created

% Test Degradation
    % folder_name  = 'TK_Test';
    % battery_name = 'Test_CTRG_Both';

    % folder_name  = 'TK_DelV_DelT_Tests';
    % % battery_name = 'Test';
    % % battery_name = 'TestAMLoss';

    folder_name  = 'TK_CyclingDegradation';
    % % battery_name = 'Iso20_ConstantCTRG_AN';
    % % battery_name = 'Iso20_ConstantCTRG_CA';
    % % battery_name = 'Iso20_ConstantCTRG_Both';
    % % battery_name = 'Iso20_ConstantAML_AN';
    % % battery_name = 'Iso20_ConstantAML_CA';
    % % battery_name = 'Iso20_ConstantAML_Both';
    % % battery_name = 'Iso20_CombinedCTRG_AML_AN';
    % % battery_name = 'Iso20_CombinedCTRG_AML_CA';
    % battery_name = 'Iso20_CombinedCTRG_AML_Both';


%% Simulations
% ---- Polarization ----
% Positive is discharge, Negative is charge
    C_rates      = [];
    % C_rates      = [0]; 
    % C_rates      = [0.040843474405010]; % Results in 1A/m^2
    % C_rates      = [-1/5 -1/2 -1 -1.5 -2 -5];
    % C_rates      = [-1/20 -1/5 -1/2 -1 -2 -5];
    % C_rates      = [-1]; 
    % C_rates      = [-5]; 
    % C_rates      = -(1:5); 
    % C_rates      = [1/20]; 
    % C_rates      = [1/20 -1/20]; 
    % C_rates      = [-1/4];
    % C_rates      = [1/20 1/10 1/3 1 2]; 
    % C_rates      = [-1/20  -1 -2]; 

% ---- Harmonic Perturbation ----
% [rad/s], frequency of the sin wave
    EIS_SIN_freq = [];
    % EIS_SIN_freq = [2*pi*0.0666666666666667];
    % EIS_SIN_freq = [1e2 1e-2];
    % EIS_SIN_freq = logspace(-2,1,31);
    % EIS_SIN_freq = 50 * (2*pi);

    exp_min = -3;
    exp_max =  5;
    exp_diff = exp_max - exp_min;
    % EIS_SIN_freq = logspace(exp_min ,exp_max ,exp_diff*10+1);

        % [%], Initial state of charge
        EIS_SOC      = [50];  

% ---- State Space EIS ----
    SS_SOC = [];
    % SS_SOC = 0:1:100;
    % SS_SOC = [5, 10, 25, 50, 75, 90, 95];
    % SS_SOC = [80.46];
    % SS_SOC = [50];
    
    % Desired frequency [rad s^-1] for impedance results (Starting with known Hz)
        % SS_freq = [];
        % SS_freq = logspace(-1,11,101);
        % SS_freq = (logspace(-2,6,75) *(2*pi));
        exp_min = -2;
        exp_max =  6;
        exp_diff = exp_max - exp_min;
        SS_freq = logspace(exp_min ,exp_max ,exp_diff*10+1);
        SS_omega = SS_freq*2*pi; % [rad/s]
    
    % % Desired frequency [rad s^-1] for impedance results
    %     % SS_omega = [];
    %     % SS_omega = logspace(-1,11,101);
    %     % SS_freq = (logspace(-2,6,75) *(2*pi));
    %     exp_min = -3;
    %     exp_max =  5;
    %     exp_diff = exp_max - exp_min;
    %     SS_omega = logspace(exp_min ,exp_max ,exp_diff*10+1);
        
        
% ---- Known BC Profile Controller ----
    KBCP   = 1;
        KBCPProfileOverwrite = 0;
        % KBCPProfileFilename = 'StairStepNoRelax';
        % KBCPProfileFilename = 'GITT';
        % KBCPProfileFilename = 'Chg_Dchg_1Cycles';
        % KBCPProfileFilename = 'Chg_Dchg_1.5Cycles';
        % KBCPProfileFilename = 'Chg_Dchg_2Cycles';
        % KBCPProfileFilename = 'Chg_Dchg_4Cycles';
        % KBCPProfileFilename = 'Chg_Dchg_10Cycles';
        % KBCPProfileFilename = 'Chg_Dchg_17Cycles';
        % KBCPProfileFilename = 'Chg_Dchg_20Cycles';

    % Initial SOC
        % KBSOC = 81.93;
        % KBSOC = 100; 
        KBSOC = 0;


% ---- MOO Controller ----
    MOO = 0;
        ControllerName = '';


% ---- Manual Current Profile ----
    ManCurrProfile = 0; % 1 if want to make a manual current profile simulation from what is currently in makeCurrentProfile.m and input.m
        MCP.plating_refine = 1; % 1 if to use the plating refinement naming structure; if 0, then ManProfileName will be used
            MCP.N_regions      = 100;  % Number of steps used in the current profile
            MCP.max_iterations = 1000; % Number of refinements used in the while loop
            MCP.tol_Delta_phi  = 0.01; % Goal for the largest delta_phi
            
        MCP.UseExistingProfile = 0; % 1 if using file found at MCP.Existing_Profile_filepath
            MCP.Existing_Profile_filepath = 'F:\TylerFiles\GitHubRepos\BatteryModelingExtras\OLD MODELS\WORKINGP2D_BeforeMovingRepos\BatteryModel\BatchMode_DAE\Results\Final_Lui_Wiley_Model\Final_Lui_Wiley_Model_ManCurrProf_100steps_1000Iter_0.03tol_CurrentProfile_Output.mat';
        MCP.ManProfileName = ''; % Use this if MCP.plating_refine = 0;
        

% ---- PRBS ---- 
    doPRBS = 0;
        PRBS_Amp     = 1;  % [A/m^2], Amplitude of PRBS Signal
        PRBS_SOC     = 50; % [%],     State of Charge
        PRBS_Tswitch = 1;  % [s],     Switching Time
        % PRBS_t_Report = 1/200; % Inverse of number of samples per Tswitch
        % %%%%%%%%%% Do I need this^^^^

        AddIntermediateRelaxTime = 0;
            NumTsRelax = 2;
            NumZeroCrossingUntilNextRelax = 5;

        MakeLongPRBSSignal = 0;
            DesiredLength  = 4e3;

        UseODEextend = 0; %%%%%%%% When do I need to do this?
            % NumCuts          = 3;
            NumCuts          = 40;
            REDUCESOLN       = 1;
            SaveIntermediate = 1;


% ---- EIS from Stiching PRBS ---- 
    getEIS_PRBS = 0;
        EIS_PRBS.Amp     = 1;  % [A/m^2], Amplitude of PRBS Signal
        EIS_PRBS.SOC     = 50; % [%],     State of Charge
        EIS_PRBS.Tswitch = [0.1 , 1 , 10];  % [s],     Switching Time
        EIS_PRBS.t_Report = 1/200; % Inverse of number of samples per Tswitch
        EIS_PRBS.AddIntermediateRelaxTime = 0;
            EIS_PRBS.NumTsRelax = 2;
            EIS_PRBS.NumZeroCrossingUntilNextRelax = 5;

        EIS_PRBS.MakeLongPRBSSignal = 0;
            EIS_PRBS.DesiredLength  = 100;
            % EIS_PRBS.DesiredLength  = 255;

        EIS_PRBS.UseODEextend = 0;
            % NumCuts          = 3;
            EIS_PRBS.NumCuts          = 40;
            EIS_PRBS.REDUCESOLN       = 1;
            EIS_PRBS.SaveIntermediate = 1;

            exp_min = -3;
            exp_max =  5;
            exp_diff = exp_max - exp_min;
        EIS_PRBS.freq = logspace(exp_min ,exp_max ,exp_diff*10+1);


% ---- EIS Ho-Kalman ----
    getEIS_HoKalman = 0;
        HK_Ts  = 1;  % [s], Sampling Time
        HK_SOC = 50; % [%], State of Charge
            exp_min = -3;
            exp_max =  5;
            exp_diff = exp_max - exp_min;
        HK_freq = logspace(exp_min ,exp_max ,exp_diff*10+1);
        HK_nSamples = 800; % [], Number of relaxation samples (2 inital relax, 1 pulse, nSamples relax)
        % HK_nSamples = 5; % [], Number of relaxation samples (2 inital relax, 1 pulse, nSamples relax)


%% Create Folder    
% Make Folder Path
    save_file_path = [current_file_path filesep 'Results' filesep folder_name];
    folder_exist = isfolder(save_file_path);
    if folder_exist
        if FLAG_local.folder_overwrite
            disp('Deleting folder with the same name and all of its contents then creating a new folder with the same name')
            rmdir(save_file_path, 's')
            mkdir(save_file_path)
        end
        if ~FLAG_local.folder_overwrite && ~FLAG_local.folder_add
            disp('Creating new folder name')
            count = 1;
            good_new_name = 0;
            while ~good_new_name
                new_name = [save_file_path '_(' num2str(count) ')'];
                % folder_exist = isfolder(new_name);
                if isfolder(new_name)
                    count = count + 1;
                else
                    mkdir(new_name)
                    save_file_path = new_name;
                    good_new_name = 1;
                end
            end
        else
            disp('Folder exist. Adding new simulations to existing folder')
        end
    else
        mkdir(save_file_path)
    end


%% Create Simulation Files
%% ---- Polarization ----
for i = 1:length(C_rates)% -1 if Charge, 1 if Discharge 
    if C_rates(i) < 0
        CorD = 'C';
        SIM.ChargeOrDischarge = -1;
        SIM.SOC_start = 5;
    elseif C_rates(i) > 0
        CorD = 'D';
        SIM.ChargeOrDischarge = 1;
        SIM.SOC_start = 95;
        % SIM.SOC_start = 100;
    else
        CorD = 'D';
        SIM.ChargeOrDischarge = 1;
        SIM.SOC_start = 10;
    end
    cRateTxt = sprintf( '%.2f' , abs(C_rates(i)) ) ;
    filename = [ battery_name , '_Polar_' , cRateTxt , 'C_' , CorD , '.mat'];
    
    % Check if sim already exist
    if isfile([save_file_path filesep filename])
        if FLAG_local.sim_overwrite
            disp('Deleting previous sim')% delete old file
            delete([save_file_path filesep filename])
        else
            disp('Simulation already exists')
        end        
    end
    
    % Create sim if it doesn't exist
    if ~isfile([save_file_path filesep filename])
        disp('Making sim')
        SIM.results_filename = filename;
        SIM.SimMode          = 1;
        SIM.C_rate           = abs(C_rates(i));
        [SIM,inputHandle]    = getFunctionHandles(InputFile,SIM);
        
        % Call Inputs
            [AN,CA,SEP,EL,SIM,N,FLAG] = inputHandle(SIM);
        
        % Call Init
            [AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS] = batt_init(AN,CA,SEP,EL,SIM,N,FLAG);
        
        % Save Simulation File
            postProcessComplete = 0;
            save([save_file_path filesep filename],'AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete')
        
        % Clear Variables
            clear AN CA SEP EL SIM CONS P N FLAG PROPS
    end
end


%% ---- Harmonic Perturbation ----
for i = 1:length(EIS_SIN_freq)
    for j = 1:length(EIS_SOC)
        filename = [ battery_name , '_EIS_SIN_w' , num2str(EIS_SIN_freq(i)) , '_SOC' , num2str(EIS_SOC(j)), '.mat' ];
        if isfile([save_file_path filesep filename])
            if FLAG_local.sim_overwrite
                disp('Deleting previous sim')% delete old file
                delete([save_file_path filesep filename])
            else
                disp('Simulation already exists')
            end
        end
        
        % Create sim if it doesn't exist
        if ~isfile([save_file_path filesep filename])
            SIM.results_filename = filename;
            SIM.SimMode          = 2;
            SIM.freq             = EIS_SIN_freq(i);
            SIM.SOC_start        = EIS_SOC(j);
            [SIM,inputHandle]    = getFunctionHandles(InputFile,SIM);

            % Call Inputs
                [AN,CA,SEP,EL,SIM,N,FLAG] = inputHandle(SIM);

            % Call Init
                [AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS] = batt_init(AN,CA,SEP,EL,SIM,N,FLAG);

            % Save File
                postProcessComplete = 0;
                save([save_file_path filesep filename],'AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete')
            
            % Clear Variables
                clear AN CA SEP EL SIM CONS P N FLAG PROPS
        end
    end
end


%% ---- State Space EIS ----
for i = 1:length(SS_SOC)
    filename = [ battery_name , '_SS_EIS_SOC' , num2str(SS_SOC(i)), '.mat' ];
    if isfile([save_file_path filesep filename])
        if FLAG_local.sim_overwrite
            disp('Deleting previous sim')% delete old file
            delete([save_file_path filesep filename])
        else
            disp('Simulation already exists')
        end
    end
    
    % Create sim if it doesn't exist
    if ~isfile([save_file_path filesep filename])
        SIM.results_filename = filename;
        SIM.SimMode          = 3;
        SIM.omega            = SS_omega;
        SIM.SOC_start        = SS_SOC(i);
        [SIM,inputHandle]    = getFunctionHandles(InputFile,SIM);

        % Call Inputs
            [AN,CA,SEP,EL,SIM,N,FLAG] = inputHandle(SIM);

        % Call Init
            [AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS] = batt_init(AN,CA,SEP,EL,SIM,N,FLAG);

        % Save File
            postProcessComplete = 0;
            save([save_file_path filesep filename],'AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete')
        
        % Clear Variables
            clear AN CA SEP EL SIM CONS P N FLAG PROPS
    end
end


%% ---- Known BC Profile ----
if KBCP    
    % Check if controller profile already exist
    profile_filename = [current_file_path filesep 'KnownProfiles' filesep KBCPProfileFilename '.mat'];
    
    if isfile(profile_filename)
        if KBCPProfileOverwrite
            disp('Deleting previous profile')% delete old file
            delete(profile_filename)
        else
            disp('Profile already exists')
        end        
    end
    
    % Create KBCP profile if it doesn't exist and load into SIM
    if isfile(profile_filename)
        disp('Using Existing Profile')
        MO_File = load(profile_filename);
        SIM.Controller_MO_File = MO_File.MO_File;
    else
        [MO_File] = getControlProfile();
        save(profile_filename,'MO_File');
        SIM.Controller_MO_File = MO_File;
    end
    
    % Check if simulation already exist    
    filename = [ battery_name , '_KPCont_' , KBCPProfileFilename , 'SOC' , num2str(KBSOC) ,  '.mat'];
    
    if isfile([save_file_path filesep filename])
        if FLAG_local.sim_overwrite
            disp('Deleting previous sim')% delete old file
            delete([save_file_path filesep filename])
        else
            disp('Simulation already exists')
        end        
    end
    
    % Create simulation if it doesn't exist
    if ~isfile([save_file_path filesep filename])
        disp('Making sim')
        SIM.results_filename = filename;
        SIM.SimMode          = 4;
        SIM.SOC_start        = KBSOC;
        [SIM,inputHandle]    = getFunctionHandles(InputFile,SIM);
        
        % Call Inputs
            [AN,CA,SEP,EL,SIM,N,FLAG] = inputHandle(SIM);
        
        % Call Init
            [AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS] = batt_init(AN,CA,SEP,EL,SIM,N,FLAG);
        
        % Save Simulation File
            postProcessComplete = 0;
            save([save_file_path filesep filename],'AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete')
        
        % Clear Variables
            clear AN CA SEP EL SIM CONS P N FLAG PROPS
    end
end


%% ---- MOO Controller ----
if MOO
    
end


%% ---- Manual Current Profile ----
if ManCurrProfile
    % Select naming convention
    if MCP.plating_refine
        filename = [ battery_name , '_ManCurrProf_' ,num2str(MCP.N_regions), 'steps_' ,num2str(MCP.max_iterations), 'Iter_' ,num2str(MCP.tol_Delta_phi), 'tol' , '.mat' ];
    else
        filename = [ battery_name , '_ManCurrProf_' , MCP.ManProfileName , '.mat' ];
    end
    
    % Check if file already exists
    if isfile([save_file_path filesep filename])
        if FLAG_local.sim_overwrite
            disp('Deleting previous sim')% delete old file
            delete([save_file_path filesep filename])
        else
            disp('Simulation already exists')
        end
    end
    
    % Create sim if it doesn't exist
    if ~isfile([save_file_path filesep filename])
        SIM.N_regions      = MCP.N_regions;
        SIM.tol_Delta_phi  = MCP.tol_Delta_phi; 
        SIM.max_iterations = MCP.max_iterations;
        SIM.plating_refine = MCP.plating_refine;
        SIM.SimMode        = 7;
        [SIM,inputHandle]  = getFunctionHandles(InputFile,SIM);
        
        % Call Inputs
            [AN,CA,SEP,EL,SIM,N,FLAG] = inputHandle(SIM);
        
        % Call Init
            [AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS] = batt_init(AN,CA,SEP,EL,SIM,N,FLAG);
        
        % Make the current profile
            if ~MCP.UseExistingProfile
                profile_save_filepath = [save_file_path filesep filename];
                [profile_save_filepath] = makeCurrentProfile(profile_save_filepath, SIM);
                variables = load(profile_save_filepath);
                SIM.region_time_vec        = variables.region_time_vec;
                SIM.region_current_vec     = variables.region_current_vec;
                SIM.profile_time           = variables.profile_time;
                SIM.profile_current        = variables.profile_current;
                SIM.tspan                  = [0, variables.t_final];
                SIM.Input_Profile_filepath = variables.profile_save_filepath;
            else
                variables = load(MCP.Existing_Profile_filepath);
                SIM.region_time_vec        = variables.region_time_vec;
                SIM.region_current_vec     = variables.region_current_vec;
                SIM.profile_time           = variables.profile_time;
                SIM.profile_current        = variables.profile_current;
                SIM.tspan                  = [0, variables.t_final];
                SIM.Input_Profile_filepath = variables.profile_save_filepath;
            end
        % Save File
            postProcessComplete = 0;
            save([save_file_path filesep filename],'AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete')
        
            % Clear Variables
            clear AN CA SEP EL SIM CONS P N FLAG PROPS
    end
end


%% ---- PRBS ----
if doPRBS
    filename = [ battery_name , '_PRBS' , '_Amp' , num2str(PRBS_Amp) , '_SOC' , num2str(PRBS_SOC) , '_SwitchingTime' , num2str(PRBS_Tswitch) , '.mat'];
    
    % Check if sim already exist
    if isfile([save_file_path filesep filename])
        if FLAG_local.sim_overwrite
            disp('Deleting previous sim')% delete old file
            delete([save_file_path filesep filename])
        else
            disp('Simulation already exists')
        end        
    end
    
    % Create sim if it doesn't exist
    if ~isfile([save_file_path filesep filename])
        disp('Making sim')
        SIM.results_filename = filename;
        SIM.SimMode   = 8;
        SIM.SOC_start = PRBS_SOC;
        SIM.Tswitch   = PRBS_Tswitch;
        % SIM.t_Report  = PRBS_t_Report ;
        SIM.PRBSAmp   = PRBS_Amp;

        SIM.AddIntermediateRelaxTime      = AddIntermediateRelaxTime;
        SIM.NumTsRelax                    = NumTsRelax;
        SIM.NumZeroCrossingUntilNextRelax = NumZeroCrossingUntilNextRelax;

        SIM.MakeLongPRBSSignal = MakeLongPRBSSignal;
        SIM.DesiredLength      = DesiredLength;
        SIM.UseODEextend       = UseODEextend;
        SIM.NumCuts            = NumCuts;
        SIM.REDUCESOLN         = REDUCESOLN;
        SIM.SaveIntermediate   = SaveIntermediate;
        [SIM,inputHandle]      = getFunctionHandles(InputFile,SIM);
        
        % Call Inputs
            [AN,CA,SEP,EL,SIM,N,FLAG] = inputHandle(SIM);
        
        % Call Init
            [AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS] = batt_init(AN,CA,SEP,EL,SIM,N,FLAG);
        
        % Save Simulation File
            postProcessComplete = 0;
            save([save_file_path filesep filename],'AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete')
        
        % Clear Variables
            clear AN CA SEP EL SIM CONS P N FLAG PROPS
    end

end


%% ---- EIS PRBS ----
if getEIS_PRBS
    % Create PRBS Simulations
        PRBS_desTswitch_filenames = cell([]);
        for i = 1:length(EIS_PRBS.Tswitch)
            Tswitch = EIS_PRBS.Tswitch(i);
            filename = [ battery_name , '_PRBS' , '_Amp' , num2str(EIS_PRBS.Amp) , '_SOC' , num2str(EIS_PRBS.SOC) , '_SwitchingTime' , num2str(Tswitch) , '.mat'];
            PRBS_desTswitch_filenames{i,1} = [save_file_path filesep filename];

            % Check if sim already exist
            if isfile([save_file_path filesep filename])
                if FLAG_local.sim_overwrite
                    disp('Deleting previous sim')% delete old file
                    delete([save_file_path filesep filename])
                else
                    disp('Simulation already exists')
                end        
            end

            % Create sim if it doesn't exist
            if ~isfile([save_file_path filesep filename])
                disp('Making sim')
                SIM.results_filename = filename;
                SIM.SimMode   = 8;
                SIM.SOC_start = EIS_PRBS.SOC;
                SIM.Tswitch   = Tswitch;
                SIM.t_Report  = EIS_PRBS.t_Report;
                SIM.PRBSAmp   = EIS_PRBS.Amp;

                SIM.AddIntermediateRelaxTime      = EIS_PRBS.AddIntermediateRelaxTime;
                SIM.NumTsRelax                    = EIS_PRBS.NumTsRelax;
                SIM.NumZeroCrossingUntilNextRelax = EIS_PRBS.NumZeroCrossingUntilNextRelax;

                SIM.MakeLongPRBSSignal = EIS_PRBS.MakeLongPRBSSignal;
                SIM.DesiredLength      = EIS_PRBS.DesiredLength;
                SIM.UseODEextend       = EIS_PRBS.UseODEextend;
                SIM.NumCuts            = EIS_PRBS.NumCuts;
                SIM.REDUCESOLN         = EIS_PRBS.REDUCESOLN;
                SIM.SaveIntermediate   = EIS_PRBS.SaveIntermediate;
                [SIM,inputHandle]      = getFunctionHandles(InputFile,SIM);

                % Call Inputs
                    [AN,CA,SEP,EL,SIM,N,FLAG] = inputHandle(SIM);

                % Call Init
                    [AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS] = batt_init(AN,CA,SEP,EL,SIM,N,FLAG);

                % Save Simulation File
                    postProcessComplete = 0;
                    save([save_file_path filesep filename],'AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete')

                % Clear Variables
                    clear AN CA SEP EL SIM CONS P N FLAG PROPS
            end
        end

    % Create Stitching File
        disp('Making Overall Sim')
        filename = [battery_name ,'_PRBS_EIS_SOC' num2str(EIS_PRBS.SOC)];
        postProcessComplete = 0;
        SIM.SimMode = 9;
        i = 1;
        P.SS.omega    = i; i = i + 1;
        P.SS.Z_mag    = i; i = i + 1;
        P.SS.Z_Re     = i; i = i + 1;
        P.SS.Z_Im     = i; i = i + 1;
        P.SS.Z_dB     = i; i = i + 1;
        P.SS.Z_ps_deg = i; i = i + 1;
        FLAG.doPostProcessing = 1;   % 1 if the postprocessing function is performed after a simulation completes
            FLAG.ReduceSolnTime = 0; % 1 if the results that are saved don't use all the points produced by t_soln ######NOT IMPLEMENTED YET
        FLAG.Plot             = 1;
        save([save_file_path filesep filename],'PRBS_desTswitch_filenames','EIS_PRBS','SIM','FLAG','P','postProcessComplete')
end


%% ---- EIS Ho-Kalman ----
if getEIS_HoKalman
    filename = [ battery_name , '_EIS_HoKalman_SOC' , num2str(HK_SOC) , '_SamplingTime' , num2str(HK_Ts) , '.mat'];
    
    % Check if sim already exist
    if isfile([save_file_path filesep filename])
        if FLAG_local.sim_overwrite
            disp('Deleting previous sim')% delete old file
            delete([save_file_path filesep filename])
        else
            disp('Simulation already exists')
        end        
    end
    
    % Create sim if it doesn't exist
    if ~isfile([save_file_path filesep filename])
        disp('Making sim')
        SIM.results_filename = filename;
        SIM.SimMode       = 10;
        SIM.SOC_start     = HK_SOC;
        SIM.Tsample       = HK_Ts;
        SIM.freq          = HK_freq;
        SIM.HK_nSamples   = HK_nSamples;
        [SIM,inputHandle] = getFunctionHandles(InputFile,SIM);
        
        % Call Inputs
            [AN,CA,SEP,EL,SIM,N,FLAG] = inputHandle(SIM);
        
        % Call Init
            [AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS] = batt_init(AN,CA,SEP,EL,SIM,N,FLAG);
        
        % Save Simulation File
            postProcessComplete = 0;
            save([save_file_path filesep filename],'AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete')
        
        % Clear Variables
            clear AN CA SEP EL SIM CONS P N FLAG PROPS
    end

end

%% Function Handles
function [SIM,inputHandle] = getFunctionHandles(InputFile,SIM)
% FLAG.InputFile
    %   -1) batt_inputs_Manual
    %    0) batt_inputs
    %    1) batt_inputs_Wiley_Lui
    %    2) batt_inputs_HZ_A
    %    3) batt_inputs_HZ_B
    %    4) batt_inputs_NREL
    %    5) batt_inputs_ToddKingston
    %    6) batt_inputs_MPC

switch InputFile
    case -1 % batt_input_Manual
        inputHandle             = @batt_inputs_Manual;
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
    case 0 % batt_input
        inputHandle             = @batt_inputs;
        SIM.ANEqPotentialHandle = @E_eqGraphite;
        SIM.ANi_oHandle         = @i_oC6;
        SIM.ANsigmaHandle       = @sigmaC6;
        SIM.AND_oHandle         = @D_o_Graphite;
        SIM.CAEqPotentialHandle = @E_eqNMC;
        SIM.CAi_oHandle         = @i_oC6;  
        SIM.CAsigmaHandle       = @sigmaNMC;
        SIM.CAD_oHandle         = @D_o_NMC532;
        SIM.ELtf_numHandle      = @transferenceNumber;
        SIM.ELActivityHandle    = @activity;
        SIM.ELD_o_Li_ionHandle  = @D_oLiion;
        SIM.ELkappaHandle       = @kappa;
    case 1 % batt_input_Wiley_Lui
        inputHandle             = @batt_inputs_Wiley_Lui;
        SIM.ANEqPotentialHandle = @E_eqGraphite;
        SIM.ANi_oHandle         = @i_oC6;
        SIM.ANsigmaHandle       = @sigmaC6;
        SIM.AND_oHandle         = @D_o_Graphite;
        SIM.CAEqPotentialHandle = @E_eqNMC;
        SIM.CAi_oHandle         = @i_oC6;  
        SIM.CAsigmaHandle       = @sigmaNMC;
        SIM.CAD_oHandle         = @D_o_NMC532;
        SIM.ELtf_numHandle      = @transferenceNumber;
        SIM.ELActivityHandle    = @activity;
        SIM.ELD_o_Li_ionHandle  = @D_oLiion;
        SIM.ELkappaHandle       = @kappa;
    case 2 % batt_input_HZ_A
        inputHandle             = @batt_inputs_HZ_A;
        SIM.ANEqPotentialHandle = @E_eqGraphite_HZ;
        SIM.ANi_oHandle         = @i_oC6_HZ;
        SIM.ANsigmaHandle       = @sigmaC6_HZ;
        SIM.AND_oHandle         = @D_o_Graphite_HZ;
        SIM.CAEqPotentialHandle = @E_eqNMC_HZ;
        SIM.CAi_oHandle         = @i_oC6_HZ;  
        SIM.CAsigmaHandle       = @sigmaNMC_HZ;
        SIM.CAD_oHandle         = @D_o_NMC532_HZ;
        SIM.ELtf_numHandle      = @transferenceNumber_HZ;
        SIM.ELActivityHandle    = @activity_HZ;
        SIM.ELD_o_Li_ionHandle  = @D_oLiion_HZ;
        SIM.ELkappaHandle       = @kappa_HZ;
    case 3 % batt_input_HZ_B
        inputHandle             = @batt_inputs_HZ_B;
        SIM.ANEqPotentialHandle = @E_eqGraphite_HZ;
        SIM.ANi_oHandle         = @i_oC6_HZ;
        SIM.ANsigmaHandle       = @sigmaC6_HZ;
        SIM.AND_oHandle         = @D_o_Graphite_HZ;
        SIM.CAEqPotentialHandle = @E_eqNMC_HZ;
        SIM.CAi_oHandle         = @i_oNMC_HZ;  
        SIM.CAsigmaHandle       = @sigmaNMC_HZ;
        SIM.CAD_oHandle         = @D_o_NMC532_HZ;
        SIM.ELtf_numHandle      = @transferenceNumber_HZ;
        SIM.ELActivityHandle    = @activity_HZ;
        SIM.ELD_o_Li_ionHandle  = @D_oLiion_HZ;
        SIM.ELkappaHandle       = @kappa_HZ;
    case 4 % batt_input_NREL
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
    case 5 % batt_input_ToddKingston
        inputHandle             = @batt_inputs_ToddKingston;
        SIM.ANEqPotentialHandle = @E_eqGraphite_NREL;
        SIM.ANi_oHandle         = @i_oC6_NREL;
        SIM.ANsigmaHandle       = @sigmaC6_NREL;
        SIM.AND_oHandle         = @D_o_Graphite_NREL;
        SIM.CAEqPotentialHandle = @E_eqNMC_NREL;
        SIM.CAi_oHandle         = @i_oNMC_NREL;  
        SIM.CAsigmaHandle       = @sigmaNMC_NREL;
        SIM.CAD_oHandle         = @D_o_NMC532_NREL;
        SIM.ELtf_numHandle      = @transferenceNumber_Landesfeind;
        SIM.ELActivityHandle    = @activity_Landesfeind;
        SIM.ELD_o_Li_ionHandle  = @D_oLiion_Landesfeind;
        SIM.ELkappaHandle       = @kappa_Landesfeind;
    case 6 % batt_input_MPC
        inputHandle             = @batt_inputs_MPC;
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
end
end


%% OOOOOOOLLLLLDDDD Stuff

% folder_name  = 'Final_Lui_Wiley_Model';
% battery_name = 'Final_Lui_Wiley_Model';

% folder_name  = 'TestingForTyrone';
% battery_name = 'TestingForTyrone';

% folder_name  = 'SeminarPres_Nov2022';
% battery_name = 'SeminarPres_Nov2022';

% % folder_name  = 'ObservabilityTest';
% % battery_name = 'ObservabilityTest';

% folder_name  = 'KalmanTest';
% battery_name = 'KalmanTest_JustPlant';

% folder_name  = 'PRBS_Sims';
% battery_name = 'PRBS_Sims';

% folder_name  = 'zeroMeanPRBS_Sims';
% battery_name = 'PRBS_Sims';

% folder_name  = 'zeroMeanPRBS_withRelax_Sims';
% battery_name = 'PRBS_Sims';

% folder_name  = 'LongerImpulse';
% battery_name = 'ObservabilityTest';
% battery_name = 'PRBS_Sims';

% folder_name  = 'LongerZeroMeanPRBS_Sims';
% battery_name = 'PRBS_Sims';

% folder_name  = 'EISCompareForPeter';
% battery_name = 'EISCompareForPeter';

% folder_name  = 'TestODExtend';
% % % battery_name = 'TestODExtend_Baseline';
% % % battery_name = 'TestODExtend_ODEextend_2Cuts';
% % % battery_name = 'TestODExtend_ODEextend_3Cuts';
% % % battery_name = 'TestODExtend_ODEextend_3Cuts_Save';
% % % battery_name = 'TestODExtend_ODEextend_RedSoln_2Cuts';
% % % battery_name = 'TestODExtend_ODEextend_RedSoln_3Cuts';
% % % battery_name = 'TestODExtend_ODEextend_RedSoln_3Cuts_Save';
% battery_name = 'TestODExtend_idata';


        % KBCPProfileFilename = 'DTImpulseTs0.1';               % 1
        % KBCPProfileFilename = 'DTImpulseTs0.121152765862859'; % 2
        % KBCPProfileFilename = 'DTImpulseTs0.146779926762207'; % 3
        % KBCPProfileFilename = 'DTImpulseTs0.177827941003892'; % 4
        % KBCPProfileFilename = 'DTImpulseTs0.215443469003188'; % 5
        % KBCPProfileFilename = 'DTImpulseTs0.261015721568254'; % 6
        % KBCPProfileFilename = 'DTImpulseTs0.316227766016838'; % 7
        % KBCPProfileFilename = 'DTImpulseTs0.383118684955729'; % 8
        % KBCPProfileFilename = 'DTImpulseTs0.464158883361278'; % 9
        % KBCPProfileFilename = 'DTImpulseTs0.562341325190349'; % 10
        % KBCPProfileFilename = 'DTImpulseTs0.681292069057961'; % 11
        % KBCPProfileFilename = 'DTImpulseTs0.825404185268018'; % 12
        % KBCPProfileFilename = 'DTImpulseTs1.0';               % 13
        % KBCPProfileFilename = 'DTImpulseTs1.21152765862859';  % 14
        % KBCPProfileFilename = 'DTImpulseTs1.46779926762207';  % 15
        % KBCPProfileFilename = 'DTImpulseTs1.77827941003892';  % 16
        % KBCPProfileFilename = 'DTImpulseTs2.15443469003188';  % 17
        % KBCPProfileFilename = 'DTImpulseTs2.61015721568254';  % 18
        % KBCPProfileFilename = 'DTImpulseTs3.16227766016838';  % 19
        % KBCPProfileFilename = 'DTImpulseTs3.83118684955729';  % 20
        % KBCPProfileFilename = 'DTImpulseTs4.64158883361278';  % 21
        % KBCPProfileFilename = 'DTImpulseTs5.62341325190349';  % 22
        % KBCPProfileFilename = 'DTImpulseTs6.81292069057961';  % 23
        % KBCPProfileFilename = 'DTImpulseTs8.25404185268018';  % 24
        % KBCPProfileFilename = 'DTImpulseTs10.0';              % 25
        % 
        % KBCPProfileFilename = 'KalmanTestDTImpulseTs1.0';
        % KBCPProfileFilename = 'KalmanTestStep';
        % KBCPProfileFilename = 'KalmanTestStep_JustPlant_NoRelax';
        % KBCPProfileFilename = 'Relax_Step';

        % KBCPProfileFilename = 'CCCV_1.5C';
        % KBCPProfileFilename = 'CCCV_2C';
        % KBCPProfileFilename = '5minImpulseLongResponse';
        % KBCPProfileFilename = 'LongTau_wStep';
    
    % % folder_name  = 'TestNewInputFile';
    % folder_name  = 'TestNewInputFileNoise';
    % % folder_name  = 'TestNewInputFileNoise_SplitSim';
    % battery_name = 'Test';
    
    % folder_name  = 'ThermalGradient';
    % battery_name = 'NoThermalGradient';
    % % battery_name = '4C_A2C_ThermalGradient';
    % % battery_name = '4C_C2A_ThermalGradient';
    
    % folder_name  = 'COETest';
    % battery_name = 'Test123';
    
    % folder_name  = 'QuickTest';
    % battery_name = 'NoisePlots';
    
    % folder_name  = 'TestImpedanceContributions';
    % battery_name = 'Standard';
    
    % folder_name  = 'TestImpedanceSparse';
    % % battery_name = 'Standard';
    % % battery_name = 'WithCOE';
    % battery_name = 'WithCOEWithVarProps';

    % folder_name  = 'SeminarFall2023';
    % % battery_name = 'M_A';
    % % battery_name = 'M_B';
    % battery_name = 'M_C';

    % % folder_name  = 'OptiBase';
    % % folder_name  = 'OptiAll';
    % folder_name  = 'OptiFinal';
    % battery_name = 'Standard';
    % % battery_name = 'VaryProps';