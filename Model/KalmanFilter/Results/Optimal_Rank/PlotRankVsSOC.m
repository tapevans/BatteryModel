%% Plot rank across SOC
% Plot the rank at which a minimum occured for a desired variable for a
% range of SOC.

clear all; close all; clc
%% Parameters to Run
% SOC_vec = [50,55];
SOC_vec = 0:5:100;
Ts_vec  = 1;


%% Simulation Type
FLAG.InputMode = 5;
FLAG.PRBSAmp = 1;
FLAG.Tswitch = 10;
FLAG.folderpathOptiBatch = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanFilter\Results\Optimal_Rank';
SIM = 'temp';


%% Plot FLAGS
FLAG.PlotSOCVsRank_Ind = 0;
FLAG.PlotSOCVsRank_Com = 1;
FLAG.PlotSOCVsRank_Com_All = 1;


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
            % Load File
            i = i+1;
            data = load(batch_filename);
            MyData(i).name = num2str(SOC);
            MyData(i).SOC  = SOC;
            [~,col] = size(data.MinAndIDX);
            for OO = 1:col-1
                ind_idx_voltage(OO) = data.MinAndIDX{OO}.ind.idx_volt;
                ind_idx_output(OO)  = data.MinAndIDX{OO}.ind.idx_othe;
            end
            for OO = col
                com_idx_output = data.MinAndIDX{OO}.com.idx_com;
            end
            MyData(i).ind_idx_voltage = ind_idx_voltage;
            MyData(i).ind_idx_output  = ind_idx_output;
            MyData(i).com_idx_output  = com_idx_output;
        else
            disp("This simulation doesn't exist")
        end
    end
end

%% Plot Results
    % Plot SOC vs Rank for Each Output using Individual ROM
    if FLAG.PlotSOCVsRank_Ind
        for i = 1:length(MyData)
            % Create SOC Vector
            SOC_vec_new(1,i) = MyData(i).SOC;
    
            % Create idx Vector
            idx_vec_volt(:,i) = MyData(i).ind_idx_voltage';
            idx_vec_outp(:,i) = MyData(i).ind_idx_output';
        end
        for OO = 1:N.DesOut
            figure
            hold on
            plot(SOC_vec_new , idx_vec_volt(OO,:),'-ko','Linewidth',2,'MarkerSize',5,'DisplayName',[RESULTS.Labels.title{1}  ' Rank'])
            plot(SOC_vec_new , idx_vec_outp(OO,:),'-ro','Linewidth',2,'MarkerSize',2,'DisplayName',[RESULTS.Labels.title{OO} ' Rank'])
            yline(10,'--')
            xlabel('State of Charge')
            ylabel('Minimum Rank')
            title([RESULTS.Labels.title{OO} ' Minimum Rank vs SOC, ROM:Individual'])
            xlim([ min(SOC_vec_new) , max(SOC_vec_new) ])
            % ylim([0 1e-1])
            lgn = legend;
        end
    end
    
    % Plot SOC vs Rank for Each Output using Combined ROM
    if FLAG.PlotSOCVsRank_Com
        for i = 1:length(MyData)
            % Create SOC Vector
            SOC_vec_new(1,i) = MyData(i).SOC;
    
            % Create idx Vector
            idx_vec_outp(:,i) = MyData(i).com_idx_output';
        end
        for OO = 1:N.DesOut
            figure
            hold on
            plot(SOC_vec_new , idx_vec_outp(OO,:),'-ro','Linewidth',2,'MarkerSize',2,'DisplayName',[RESULTS.Labels.title{OO} ' Rank'])
            yline(10,'--')
            xlabel('State of Charge')
            ylabel('Minimum Rank')
            title([RESULTS.Labels.title{OO} ' Minimum Rank vs SOC, ROM:Combined'])
            xlim([ min(SOC_vec_new) , max(SOC_vec_new) ])
            % ylim([0 1e-1])
            lgn = legend;
        end
    end

    % All Combined on the same plot
    if FLAG.PlotSOCVsRank_Com_All
        figure
        hold on
        for OO = 1:N.DesOut
            plot(SOC_vec_new , idx_vec_outp(OO,:),'-o','Linewidth',2,'MarkerSize',2,'DisplayName',[RESULTS.Labels.title{OO} ' Rank'])
        end
        xlim([ min(SOC_vec_new) , max(SOC_vec_new) ])
        % ylim([0 1e-1])
        lgn = legend;
        xlabel('State of Charge')
        ylabel('Minimum Rank')
        title(['All Minimum Rank vs SOC, ROM:Combined'])
        yline(10,'--')
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
