%% Create Project
% This file is used to create a project folder. In this folder is a file of
% battery parameters used, which simulations are to be performed, and a
% data file (.mat) for each simulation to be performed.

% Simulations generated with this file will use the values in batt_inputs
% for other simulation parameters such as temperature, voltage limits, how
% deep to discharge, and etc...

% File names
    % ---- Polarization ----
        % 'battery_name'_Polar_'Crate'C_(CorD)
    
    % ---- Harmonic Perturbation ----
        % 'battery_name'_EIS_SIN_w'freq'_SOC'SOC'
    
    % ---- State Space EIS ----
        % 'battery_name'_SS_EIS_SOC'SOC'
    
    % ---- Known BC Profile Controller ----
        % 'battery_name'_KPCont_'KBCPProfileFilename'_SOC'KBSOC'
        
    % ---- MOO Controller ----
        % 'battery_name'_MOOCont_'ControllerName'
    
    % ---- Manual Current Profile ----
        % For platingRefinement
            % 'battery_name'_ManCurrProf_ 'MCP.N_regions' steps_ 'MCP.max_iterations' Iter_ 'MCP.tol_Delta_phi' tol 
        % For other
            % 'battery_name'_ManCurrProf_'ManualProfName'
            
        
        
%%
clear all; close all; clc;

%% Subdirecties to Include
% Script's filepath
[current_file_path,~,~] = fileparts(mfilename('fullpath'));

% Include all folders
addpath(genpath(current_file_path)); 
% genpath creates a string with all folders and subfolders in a given directory
% addpath then adds all of them to the current workspace

%% Inputs
FLAG_local.folder_overwrite = 0; % 1 if delete folder if it already exists
FLAG_local.folder_add       = 1; % 1 if just want to add simulations to folder
% if folder_overwrite and folder_add are 0, then a new name should be used
% for the folder

FLAG_local.sim_overwrite    = 1; % 1 if older simulation is deleted and new one is created

% folder_name  = 'Final_Lui_Wiley_Model';
% battery_name = 'Final_Lui_Wiley_Model';

folder_name  = 'KBCP_Mode_Test';
battery_name = 'KBCP_Mode_Test';

% ---- Polarization ----
% Positive is discharge, Negative is charge
    C_rates      = []; 
%     C_rates      = [-1/5 -1/2 -1 -1.5 -2 -5]; 
%     C_rates      = [-1]; 
%     C_rates      = [1/20 1/10]; 
%     C_rates      = [1/20];
%     C_rates      = [1/20 1/10 1/3 1 2]; 

% ---- Harmonic Perturbation ----
% [rad/s], frequency of the sin wave
    EIS_SIN_freq = [];
%     EIS_SIN_freq = [1e2 1e-2];
    % EIS_SIN_freq = logspace(-2,1,31);

        % [%], Initial state of charge
        EIS_SOC      = [50];  

% ---- State Space EIS ----
    SS_SOC = [];
%     SS_SOC = [5, 10, 25, 50, 75, 90, 95];
%     SS_SOC = [81.93];
    
%         SS_freq = [];
        SS_freq = logspace(-1,11,101);
%         SS_freq = (logspace(-2,6,75) *(2*pi));
        
% ---- Known BC Profile Controller ----
    KBCP   = 1;
        KBCPProfileOverwrite = 1;
        KBCPProfileFilename = 'Profile_CCChg4.2_CCDchg3.4_C3';
%         KBSOC = 81.93;
        KBSOC = 0;
        
% ---- MOO Controller ----
    MOO = 0;
        ControllerName = '';
        
% ---- Manual Current Profile ----
    ManCurrProfile = 0; % 1 if want to make a manual current profile simulation from what is currently in makeCurrentProfile.m and input.m
        MCP.plating_refine = 1; % 1 if to use the plating refinement naming structure; if 0, then ManProfileName will be used
            MCP.N_regions  = 100;  % Number of steps used in the current profile
            MCP.max_iterations = 1000; % Number of refinements used in the while loop
            MCP.tol_Delta_phi  = 0.01; % Goal for the largest delta_phi
            
        MCP.UseExistingProfile = 1; % 1 if using file found at MCP.Existing_Profile_filepath
            MCP.Existing_Profile_filepath = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\Final_Lui_Wiley_Model\Final_Lui_Wiley_Model_ManCurrProf_100steps_1000Iter_0.03tol_CurrentProfile_Output.mat';
        MCP.ManProfileName = ''; % Use this if MCP.plating_refine = 0;
        
        
        
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
%             folder_exist = isfolder(new_name);
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
%%
% ---- Polarization ----
for i = 1:length(C_rates)% -1 if Charge, 1 if Discharge 
    if C_rates(i) < 0
        CorD = 'C';
        SIM.ChargeOrDischarge = -1;
        SIM.SOC_start = 0;
    else
        CorD = 'D';
        SIM.ChargeOrDischarge = 1;
        SIM.SOC_start = 95;
%         SIM.SOC_start = 100;
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
        SIM.SimMode = 1;
        SIM.C_rate = abs(C_rates(i));
        
        % Call Inputs
        [AN,CA,SEP,EL,SIM,N,FLAG] = batt_inputs(SIM);
        
        % Call Init
        [AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS] = batt_init(AN,CA,SEP,EL,SIM,N,FLAG);
        
        % Save Simulation File
        save([save_file_path filesep filename],'AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS')
        
        % Clear Variables
        clear AN CA SEP EL SIM CONS P N FLAG PROPS
    end
end

%%
% ---- Harmonic Perturbation ----
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
            SIM.SimMode = 2;
            SIM.freq = EIS_SIN_freq(i);
            SIM.SOC_start = EIS_SOC(j);
            % Call Inputs
            [AN,CA,SEP,EL,SIM,N,FLAG] = batt_inputs(SIM);
            % Call Init
            [AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS] = batt_init(AN,CA,SEP,EL,SIM,N,FLAG);
            % Save File
            save([save_file_path filesep filename],'AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS')
            % Clear Variables
            clear AN CA SEP EL SIM CONS P N FLAG PROPS
        end
    end
end

%%
% ---- State Space EIS ----
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
        SIM.SimMode = 3;
        SIM.freq = SS_freq;
        SIM.SOC_start = SS_SOC(i);
        % Call Inputs
        [AN,CA,SEP,EL,SIM,N,FLAG] = batt_inputs(SIM);
        % Call Init
        [AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS] = batt_init(AN,CA,SEP,EL,SIM,N,FLAG);
        % Save File
        save([save_file_path filesep filename],'AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS')
        % Clear Variables
        clear AN CA SEP EL SIM CONS P N FLAG PROPS
    end
end

%%
% ---- Known BC Profile ----
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
        SIM.SimMode = 4;
        SIM.SOC_start = KBSOC;
        
        % Call Inputs
        [AN,CA,SEP,EL,SIM,N,FLAG] = batt_inputs(SIM);
        
        % Call Init
        [AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS] = batt_init(AN,CA,SEP,EL,SIM,N,FLAG);
        
        % Save Simulation File
        save([save_file_path filesep filename],'AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS')
        
        % Clear Variables
        clear AN CA SEP EL SIM CONS P N FLAG PROPS
    end
end

%% ---- MOO Controller ----
if MOO
    
end

%%
% ---- Manual Current Profile ----
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
        SIM.SimMode = 7;
        
        % Call Inputs
            [AN,CA,SEP,EL,SIM,N,FLAG] = batt_inputs(SIM);
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
            save([save_file_path filesep filename],'AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS')
        % Clear Variables
            clear AN CA SEP EL SIM CONS P N FLAG PROPS
    end
end
