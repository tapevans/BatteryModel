%% Adjust Optimal Rank

%%
clear all; close all; clc;

%% Parameters to Run
SOC_vec = 0:5:100;
% SOC_vec = 0;
Ts_vec  = 1;


%% Simulation Type
FLAG.InputMode = 5;
FLAG.PRBSAmp   = 1;
FLAG.Tswitch   = 10;
FLAG.folderpathOptiBatch = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanFilter\Results\Optimal_Rank';
SIM = 'temp';


%% Plot FLAGS
% FLAG.PlotSOCVsRank_Ind     = 1;
% FLAG.PlotSOCVsRank_Com     = 0;
% FLAG.PlotSOCVsRank_Com_All = 0;


%% Get Labels
[RESULTS , N, P] = getPlotLabels();


%% Plot Results
% Loop through SOC
i = 0;
for SS = 1:length(SOC_vec)
    SOC = SOC_vec(SS);
    disp(['SOC: ' num2str(SOC)])
    % Loop through Ts
    for TT = Ts_vec
        FLAG.SOC = SOC;
        FLAG.Ts  = TT;
        %% Check if batch exists
        [batch_filename] = getBatchFilename(SIM,FLAG);
        if isfile(batch_filename)
            load(batch_filename);
            % if exist('adjusted')
                

            % else
                adjusted = 1;
                for OO = 1:N.DesOut
                    % Voltage
                    adjusted_volt = MinAndIDX{OO}.ind.min_v * 10;
                    adj_idx = find(All_ErrorSqSum{OO}(1,:)<adjusted_volt,1);
                    if isempty(adj_idx)
                        adj_idx = MinAndIDX{OO}.ind.idx_volt;
                    end
                    adjustedData{OO}.idx_volt = adj_idx;
                    adjustedData{OO}.min_v    = All_ErrorSqSum{OO}(1,adj_idx);

                    % Output
                    adjusted_outp = MinAndIDX{OO}.ind.min_v * 10;
                    adj_idx = find(All_ErrorSqSum{OO}(2,:)<adjusted_outp,1);
                    if isempty(adj_idx)
                        adj_idx = MinAndIDX{OO}.ind.idx_othe;
                    end
                    adjustedData{OO}.idx_othe = adj_idx;
                    adjustedData{OO}.min_o    = All_ErrorSqSum{OO}(2,adj_idx);
                end
                % Combined
                for OO = 1:N.DesOut
                    % Output
                    adjusted_outp = MinAndIDX{end}.com.min_com(OO) * 10;
                    adj_idx = find(All_ErrorSqSum{end}(OO,:)<adjusted_outp,1);
                    if isempty(adj_idx)
                        adj_idx = MinAndIDX{end}.com.idx_com(OO);
                    end
                    adjustedData{N.DesOut+1}.idx_outp(OO) = adj_idx;
                    adjustedData{N.DesOut+1}.min_o(OO)    = All_ErrorSqSum{end}(OO,adj_idx);
                end

                % Save Results
                save(batch_filename,'All_ErrorSqSum','MinAndIDX','SOC','Ts','RESULTS','P','N','adjustedData','adjusted') 
            % end
        else
            disp("This simulation doesn't exist")
        end
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RESULTS , N, P] = getPlotLabels()
    FLAG.DesOut.cell_voltage  = 1; % Terminal Voltage
    FLAG.DesOut.delta_phi     = 1; % Electrostatic potential difference between active material and electrolyte @AN/SEP interface
    FLAG.DesOut.i_Far         = 0; % Chemical reaction current density at SEI @AN/SEP interface
    FLAG.DesOut.eta           = 1; % Overpotential at SEI @AN/SEP interface
    FLAG.DesOut.C_Liion       = 1; % Concentration of Li^+ in the electrolyte @AN/SEP interface
    FLAG.DesOut.C_Li          = 1; % Concentration of Li at the surface of active material @AN/SEP interface
    FLAG.DesOut.delta_C_Li    = 1; % Difference of Li concentration between the surface and center of active material particle @AN/SEP interface
    FLAG.DesOut.T             = 0; % Temperature of control volume @AN/SEP interface
    
    i = 1;
    if FLAG.DesOut.cell_voltage   % Cell Voltage
        P.cell_voltage = i;
        RESULTS.Labels.title{i} = 'Cell Voltage';
        RESULTS.Labels.unit{i}  = 'Voltage [V]';
        RESULTS.foldername{i}   = 'CellVoltage';
        i = i + 1;
    end
    if FLAG.DesOut.delta_phi      % Delta Phi
        P.delta_phi = i;
        RESULTS.Labels.title{i} = '\Delta \phi';
        RESULTS.Labels.unit{i}  = 'Voltage [V]';
        RESULTS.foldername{i}   = 'DeltaPhi';
        i = i + 1;
    end
    if FLAG.DesOut.i_Far          % i_Far
        P.i_Far = i;
        RESULTS.Labels.title{i} = 'i_{Far}';
        RESULTS.Labels.unit{i}  = 'Current [A/m^2]';
        RESULTS.foldername{i}   = 'iFar';
        i = i + 1;
    end
    if FLAG.DesOut.eta          % eta
        P.eta = i;
        RESULTS.Labels.title{i} = '\eta';
        RESULTS.Labels.unit{i}  = 'Voltage [V]';
        RESULTS.foldername{i}   = 'eta';
        i = i + 1;
    end
    if FLAG.DesOut.C_Liion        % C_Li+
        P.C_Liion = i;
        RESULTS.Labels.title{i} = 'C_{Li^+}';
        RESULTS.Labels.unit{i}  = 'Concentration [kmol/m^3]';
        RESULTS.foldername{i}   = 'C_Liion';
        i = i + 1;
    end
    if FLAG.DesOut.C_Li           % C_Li
        P.C_Li = i;
        RESULTS.Labels.title{i} = 'C_{Li,surf}';
        RESULTS.Labels.unit{i}  = 'Concentration [kmol/m^3]';
        RESULTS.foldername{i}   = 'C_Li_surf';
        i = i + 1;
    end
    if FLAG.DesOut.delta_C_Li     % Delta C_Li
        P.delta_C_Li = i;
        RESULTS.Labels.title{i} = '\Delta C_{Li}';
        RESULTS.Labels.unit{i}  = 'Concentration [kmol/m^3]';
        RESULTS.foldername{i}   = 'DeltaC_Li';
        i = i + 1;
    end
    if FLAG.DesOut.T              % Temperature
        P.T = i;
        RESULTS.Labels.title{i} = 'Temperature';
        RESULTS.Labels.unit{i}  = 'Temperature [K]';
        RESULTS.foldername{i}   = 'Temperature';
        i = i + 1;
    end
    
    N.DesOut = i-1;
end
