%% Plot Opti Error
% Plot the sum of the squared error for each desired output across all
% ranks. Do this for a given simulation, SOC, and sampling rate

% NOTES:
% Tswitch10, Ts=1, 
%    SOC 90,100 didn't run
%    Voltage looks bad on all individual plots
%%
clear all; close all; clc

%% Parameters to Run
SOC_vec = 20;
% SOC_vec = 35:5:100;
Ts_vec  = 1;


%% Simulation Type
FLAG.InputMode = 5;
FLAG.PRBSAmp = 1;
FLAG.Tswitch = 10;
FLAG.folderpathOptiBatch = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanFilter\Results\Optimal_Rank';
SIM = 'temp';


%% Plot FLAGS
FLAG.PlotErrorVsRank_Ind = 0;
FLAG.PlotErrorVsRank_Com = 1;


%% Get Labels
[RESULTS , N, P] = getPlotLabels();


%% Plot Results
% Loop through SOC
for SS = SOC_vec
    % disp(['SOC: ' num2str(SS)])
    % disp(datestr(datetime));
    % Loop through Ts
    for TT = Ts_vec
        FLAG.SOC = SS;
        FLAG.Ts  = TT;
        %% Check if batch exists
        RUNSIM = false;
        [batch_filename] = getBatchFilename(SIM,FLAG);
        if isfile(batch_filename)
            RUNSIM = true;
        else
            RUNSIM = false;
            disp("This simulation doesn't exist")
        end

        %% Plot Results
        if RUNSIM
            % Load File
            data = load(batch_filename);

            [row,col] = size(data.All_ErrorSqSum);
            % Plot Error vs Rank for Each Output using Individual ROM
            if FLAG.PlotErrorVsRank_Ind
                for OO = 1:col-1
                    figure
                    % plot(1:1:length(ind.saveErrorSqSum) , Tswitch10.ind.saveErrorSqSum(2,:),'-k','Linewidth',2)
                    semilogy(1:1:length( data.All_ErrorSqSum{OO}(2,:) ) , data.All_ErrorSqSum{OO}(2,:),'-k','Linewidth',2)
                    hold on
                    plot(data.MinAndIDX{OO}.ind.idx_othe , data.MinAndIDX{OO}.ind.min_o,'bo','Linewidth',2,'MarkerSize',10,'DisplayName',[RESULTS.Labels.title{OO} ' Minimum'])
                    plot(data.MinAndIDX{OO}.ind.idx_volt , data.All_ErrorSqSum{OO}(2,data.MinAndIDX{OO}.ind.idx_volt),'ro','Linewidth',2,'MarkerSize',10,'DisplayName','Voltage Minimum')
                    xlabel('Number of States')
                    ylabel(RESULTS.Labels.unit{OO})
                    title([RESULTS.Labels.title{OO} ' \Sigma error^2, SOC: ' num2str(SS) ' ROM: Individual'])
                    xlim([1,length(data.All_ErrorSqSum{OO}(2,:))])
                    % ylim([0 1e-1])
                    lgn = legend;
                end
            end

            % Plot Error vs Rank for Each Output using Combined ROM
            if FLAG.PlotErrorVsRank_Com
                for OO = 1:col-1
                    figure
                    % plot(1:1:length(ind.saveErrorSqSum) , Tswitch10.ind.saveErrorSqSum(2,:),'-k','Linewidth',2)
                    semilogy(1:1:length( data.All_ErrorSqSum{end}(OO,:) ) , data.All_ErrorSqSum{end}(OO,:),'-k','Linewidth',2)
                    hold on
                    plot(data.MinAndIDX{end}.com.idx_com(OO) , data.MinAndIDX{end}.com.min_com(OO),'bo','Linewidth',2,'MarkerSize',10,'DisplayName',[RESULTS.Labels.title{OO} ' Minimum'])
                    plot(data.MinAndIDX{end}.com.idx_com(1)  , data.All_ErrorSqSum{end}(OO,data.MinAndIDX{end}.com.idx_com(1)) ,'ro','Linewidth',2,'MarkerSize',10,'DisplayName','Voltage Minimum')
                    xlabel('Number of States')
                    ylabel(RESULTS.Labels.unit{OO})
                    title([RESULTS.Labels.title{OO} ' \Sigma error^2, SOC: ' num2str(SS) ' ROM: Combined'])
                    xlim([1,length(data.All_ErrorSqSum{OO}(2,:))])
                    % ylim([0 1e-1])
                    lgn = legend;
                end
            end
        end
    end
end


%% Arrange Figures
FigArrange = 1;
if FigArrange == 1
    fig = gcf;
    NumFig = fig.Number;

    Ncol = 3;
    
    for i = 1:NumFig
        f = figure(i);
        k = mod(i-1,Ncol);
        row = mod(fix((i-1)/Ncol),2);
        if row == 0
            r = 575;
%             r = 540;
        elseif row == 1
            r = 62;
        end
        f.Position = [k*575+15 r 560 420];
    end
end

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
