%% Magnitude of CPCT vs SOC and Ts
clear all; close all; clc;


%% Define Noise to compare
    % Q_0 = 1e-6; R_0 = 1e-6; % TL
    % Q_0 = 1e-6; R_0 = 1e1;  % BL
    % Q_0 = 1e1;  R_0 = 1e1;  % BR
    % Q_0 = 1e1;  R_0 = 1e-6; % TR
    Q_0 = 1e-3; R_0 = 1e-3; % Middle
    % Q_0 = 1e-3; R_0 = 1e-6; % 

    FLAG.QMode = 1; % 1) Input Q   2) State Q


%% Local FLAGS
    localFLAG.Overwrite.MyData   = 0;
    
    localFLAG.IDVorCOM                    = 0; % 1 if individual (IDV), 0 if combined (COM)
    localFLAG.Plot.NotNormalized          = 0;
    localFLAG.Plot.Normalized2CellVoltage = 0;
    localFLAG.Plot.Normalized2MaxCPCT     = 1;
    
    localFLAG.Plot.cell_voltage  = 1; % Terminal Voltage
    localFLAG.Plot.delta_phi     = 1; % Electrostatic potential difference between active material and electrolyte @AN/SEP interface
    localFLAG.Plot.i_Far         = 0; % Chemical reaction current density at SEI @AN/SEP interface
    localFLAG.Plot.eta           = 1; % Overpotential at SEI @AN/SEP interface
    localFLAG.Plot.C_Liion       = 1; % Concentration of Li^+ in the electrolyte @AN/SEP interface
    localFLAG.Plot.C_Li          = 1; % Concentration of Li at the surface of active material @AN/SEP interface
    localFLAG.Plot.delta_C_Li    = 1; % Difference of Li concentration between the surface and center of active material particle @AN/SEP interface
    localFLAG.Plot.T             = 0; % Temperature
    
    if localFLAG.Plot.NotNormalized
        cmin = -14;
        cmax = -7;
    elseif localFLAG.Plot.Normalized2CellVoltage
        cmin = 0;
        cmax = 3;
    elseif localFLAG.Plot.Normalized2MaxCPCT 
        cmin = -6;
        cmax =  0;
    end
    
    colormapping = 'jet';


%% Make string for filenames to read in
    NoiseStr = ['Q' num2str(Q_0) '_R' num2str(R_0)];
    % if localFLAG.IDVorCOM
    %     IDVorCOM = '_IDV';
    % else
    %     IDVorCOM = '_COM';
    % end


%% Filepath 
% Assumes starting in 'KalmanFilter' folder
    if FLAG.QMode == 1 % InputQ
        oldFolder = cd([pwd filesep 'Results' filesep 'ComparisonData' filesep 'InputQ']);
    else
        oldFolder = cd([pwd filesep 'Results' filesep 'ComparisonData' filesep 'StateQ']);
    end


%% Check if Data has already been consolidated
    data_filename = ['MyData_' NoiseStr '.mat'];
    % data_filename = ['MyData_' NoiseStr IDVorCOM '.mat'];
    ReadData = false;
    if isfile(data_filename)
        if localFLAG.Overwrite.MyData
            ReadData = true;
            delete(data_filename)
        end
    else
        ReadData = true;
    end
    

%% Read in Data
    if ReadData
        % Get Comparison filenames
            list = dir(['*' NoiseStr '*']);
            % list = dir(['*' NoiseStr IDVorCOM '*']);
            num_files = length(list);
            for j = 1:num_files % Creates a cell array with all simulations' full path name
                if ~exist('sim_filenames')
                    sim_filenames{1} = [pwd filesep list(j).name];
                else
                    sim_filenames{end+1,1} = [pwd filesep list(j).name];
                end
            end
    
        % Loop through all files and save data to a struct
            for i = 1:num_files
                sim = load(sim_filenames{i});
                MyData(i).name      = sim_filenames{i};
                MyData(i).Q0        = sim.Q_0;
                MyData(i).R0        = sim.R_0;
                MyData(i).Ts        = sim.Ts;
                MyData(i).SOC       = sim.SOC;
                MyData(i).SVD       = sim.sing_val_norm;
                MyData(i).CPCT      = sim.CPCT;
                MyData(i).dot       = sim.matrixOFdot;
                MyData(i).deg       = sim.matrixOFdot_deg;
                MyData(i).S_Orm     = sim.S_Orm;
                MyData(i).sing_val  = sim.sing_val;
                MyData(i).P         = sim.P     ;
            end
    
        % Save struct
            save(data_filename, 'MyData')
    
        %Go back to oldFolder
            cd(oldFolder);
    else
        load(data_filename);
        cd(oldFolder);
    end

    N_sims = length(MyData);


%% Add Normalized Values
    for i = 1:N_sims
        for OO = 1:length(MyData(i).CPCT)
            MyData(i).CPCT_diag{OO}       = diag(MyData(i).CPCT{OO});
            MyData(i).CPCT_normalized{OO} = MyData(i).CPCT_diag{OO} / MyData(i).CPCT_diag{OO}(MyData(i).P.cell_voltage); % Normalized WRT voltage CPCT
            MyData(i).CPCT_normMAX{OO}    = MyData(i).CPCT_diag{OO} / max(MyData(i).CPCT_diag{OO}); % Normalized to the largest CPCT
        end
    end


%% Axis Variables
% State of Charge
    SOC_vec = 0:1:100;
    % SOC_vec = 0:5:100;

% Desired Sampling Rate
    N_t_s   = 25;   % **Keep this fix for now
    T_s_min = -1; % T_s = 10^(T_s_min), **Keep this fix for now
    T_s_max =  1; % T_s = 10^(T_s_max), **Keep this fix for now
    Ts_vec = logspace(T_s_min,T_s_max,N_t_s);

    % SOC_vec = [50,55];
    % Ts_vec  = Ts_vec([12,13]);


%% Initialize variables for Contour plots
    T_s_linear = log10(Ts_vec);
    [X,Y] = meshgrid(T_s_linear , SOC_vec);


%% Match MyData idx with plot location
    IDX_Mat = nan(size(X));
    for SS = 1:length(SOC_vec)
        for TT = 1:length(Ts_vec)
            SOC = SOC_vec(SS);
            Ts  = Ts_vec(TT);
            for i = 1:N_sims
                if MyData(i).SOC == SOC && MyData(i).Ts == Ts
                    IDX_Mat(SS,TT) = i;
                end
            end
        end
    end


%% Save CPCT Data for plotting
n_Outs = length(fieldnames(MyData(1).P));
    for SS = 1:length(SOC_vec)
        for TT = 1:length(Ts_vec)
            if localFLAG.IDVorCOM % Individual
                for OO = 1:n_Outs
                    Z_CPCT(SS,TT,OO)         = log10( MyData(IDX_Mat(SS,TT)).CPCT_diag{OO}(2) );
                    Z_CPCT_norm(SS,TT,OO)    = log10( MyData(IDX_Mat(SS,TT)).CPCT_normalized{OO}(2) );
                    Z_CPCT_normMax(SS,TT,OO) = log10( MyData(IDX_Mat(SS,TT)).CPCT_normMAX{OO}(2) );
                end
                if SS == 1 && TT == 1
                    Z_CPCT_norm(SS,TT,1) = 1.00000000000000000000000000001;
                end
            else % Combined
                for OO = 1:n_Outs
                    Z_CPCT(SS,TT,OO)         = log10( MyData(IDX_Mat(SS,TT)).CPCT_diag{end}(OO) );
                    Z_CPCT_norm(SS,TT,OO)    = log10( MyData(IDX_Mat(SS,TT)).CPCT_normalized{end}(OO) );
                    Z_CPCT_normMax(SS,TT,OO) = log10( MyData(IDX_Mat(SS,TT)).CPCT_normMAX{end}(OO) );
                end
                if SS == 1 && TT == 1
                    Z_CPCT_norm(SS,TT,1) = 1.00000000000000000000000000001;
                end
            end
        end
    end
    

%% Plot Titles
% Title Vec
    i = 1;
    if isfield(MyData(1).P,'cell_voltage')   % Cell Voltage
        P.cell_voltage = i; 
        RESULTS.Labels.title{i} = 'Cell Voltage';
        RESULTS.Labels.unit{i}  = 'Voltage [V]';
        i = i + 1;       
    end
    if isfield(MyData(1).P,'delta_phi')      % Delta Phi
        P.delta_phi = i; 
        RESULTS.Labels.title{i} = '\Delta \phi';
        RESULTS.Labels.unit{i}  = 'Voltage [V]';
        i = i + 1;    
    end
    if isfield(MyData(1).P,'i_Far')          % i_Far
        P.i_Far = i;
        RESULTS.Labels.title{i} = 'i_{Far}';
        RESULTS.Labels.unit{i}  = 'Current [A/m^2]';
        i = i + 1;
    end
    if isfield(MyData(1).P,'eta')          % eta
        P.eta = i;
        RESULTS.Labels.title{i} = '\eta';
        RESULTS.Labels.unit{i}  = 'Voltage [V]';
        i = i + 1;
    end
    if isfield(MyData(1).P,'C_Liion')        % C_Li+
        P.C_Liion = i; 
        RESULTS.Labels.title{i} = 'C_{Li^+}';
        RESULTS.Labels.unit{i}  = 'Concentration [kmol/m^3]';
        i = i + 1;
    end
    if isfield(MyData(1).P,'C_Li')          % C_Li
        P.C_Li = i;
        RESULTS.Labels.title{i} = 'C_{Li,surf}';
        RESULTS.Labels.unit{i}  = 'Concentration [kmol/m^3]';
        i = i + 1;
    end
    if isfield(MyData(1).P,'delta_C_Li')    % Delta C_Li
        P.delta_C_Li = i;
        RESULTS.Labels.title{i} = '\Delta C_{Li}';
        RESULTS.Labels.unit{i}  = 'Concentration [kmol/m^3]';
        i = i + 1;
    end
    if isfield(MyData(1).P,'T')              % Temperature
        P.T = i;
        RESULTS.Labels.title{i} = 'Temperature';
        RESULTS.Labels.unit{i}  = 'Temperature [K]';
        i = i + 1;
    end

    idx = find(NoiseStr == '_');
    NoiseStr(idx) = ' ';

%% Plots
%% Contour not normalized
if localFLAG.Plot.NotNormalized
    if localFLAG.Plot.cell_voltage  && isfield(P,'cell_voltage')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT(:,:,P.cell_voltage));
        title([RESULTS.Labels.title{P.cell_voltage} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.delta_phi     && isfield(P,'delta_phi')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT(:,:,P.delta_phi));
        title([RESULTS.Labels.title{P.delta_phi} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.i_Far         && isfield(P,'i_Far')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT(:,:,P.i_Far));
        title([RESULTS.Labels.title{P.i_Far} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.eta           && isfield(P,'eta')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT(:,:,P.eta));
        title([RESULTS.Labels.title{P.eta} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.C_Liion       && isfield(P,'C_Liion')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT(:,:,P.C_Liion));
        title([RESULTS.Labels.title{P.C_Liion} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.C_Li          && isfield(P,'C_Li')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT(:,:,P.C_Li));
        title([RESULTS.Labels.title{P.C_Li} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.delta_C_Li    && isfield(P,'delta_C_Li')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT(:,:,P.delta_C_Li));
        title([RESULTS.Labels.title{P.delta_C_Li} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.T             && isfield(P,'T')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT(:,:,P.T));
        title([RESULTS.Labels.title{P.T} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
end


%% Contour normalized to Cell Voltage
if localFLAG.Plot.Normalized2CellVoltage
    if localFLAG.Plot.cell_voltage  && isfield(P,'cell_voltage')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT_norm(:,:,P.cell_voltage));
        title([RESULTS.Labels.title{P.cell_voltage} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
%         colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.delta_phi     && isfield(P,'delta_phi')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT_norm(:,:,P.delta_phi));
        title([RESULTS.Labels.title{P.delta_phi} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
%         colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        colormap(f, flipud(colormap(colormapping)))
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.i_Far         && isfield(P,'i_Far')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT_norm(:,:,P.i_Far));
        title([RESULTS.Labels.title{P.i_Far} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.eta           && isfield(P,'eta')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT_norm(:,:,P.eta));
        title([RESULTS.Labels.title{P.eta} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.C_Liion       && isfield(P,'C_Liion')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT_norm(:,:,P.C_Liion));
        title([RESULTS.Labels.title{P.C_Liion} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.C_Li          && isfield(P,'C_Li')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT_norm(:,:,P.C_Li));
        title([RESULTS.Labels.title{P.C_Li} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.delta_C_Li    && isfield(P,'delta_C_Li')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT_norm(:,:,P.delta_C_Li));
        title([RESULTS.Labels.title{P.delta_C_Li} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.T             && isfield(P,'T')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT_norm(:,:,P.T));
        title([RESULTS.Labels.title{P.T} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
end


%% Countour normalized to max CPCT value
if localFLAG.Plot.Normalized2MaxCPCT
    if localFLAG.Plot.cell_voltage  && isfield(P,'cell_voltage')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT_normMax(:,:,P.cell_voltage));
        title([RESULTS.Labels.title{P.cell_voltage} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.delta_phi     && isfield(P,'delta_phi')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT_normMax(:,:,P.delta_phi));
        title([RESULTS.Labels.title{P.delta_phi} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.i_Far         && isfield(P,'i_Far')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT_normMax(:,:,P.i_Far));
        title([RESULTS.Labels.title{P.i_Far} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.eta           && isfield(P,'eta')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT_normMax(:,:,P.eta));
        title([RESULTS.Labels.title{P.eta} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.C_Liion       && isfield(P,'C_Liion')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT_normMax(:,:,P.C_Liion));
        title([RESULTS.Labels.title{P.C_Liion} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.C_Li          && isfield(P,'C_Li')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT_normMax(:,:,P.C_Li));
        title([RESULTS.Labels.title{P.C_Li} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.delta_C_Li    && isfield(P,'delta_C_Li')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT_normMax(:,:,P.delta_C_Li));
        title([RESULTS.Labels.title{P.delta_C_Li} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
    if localFLAG.Plot.T             && isfield(P,'T')
        f = figure;
        [~,h] = contourf(X,Y,Z_CPCT_normMax(:,:,P.T));
        title([RESULTS.Labels.title{P.T} ' ' NoiseStr])
        xlabel('log_{10}(T_s) [s]')
        ylabel('SOC [%]')
        colorbar
        colormap(f, flipud(colormap(colormapping)))
        % clim([cmin cmax])
        h.LineStyle = 'none';
    end
end


%% Arrange Figures
FigArrange = 1;
if FigArrange == 1
    Ncol = 3;

    fig = gcf;
    NumFig = fig.Number;

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