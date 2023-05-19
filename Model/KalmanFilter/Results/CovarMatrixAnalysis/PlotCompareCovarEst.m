%% Compare Estimation Covariance Matrix Analysis
FLAG.CLEARDATA = 0;
if FLAG.CLEARDATA
    clearvars -EXCEPT FLAG
end
    close all; clc;


%% Switch to current directory
    [current_file_path,~,~] = fileparts(mfilename('fullpath'));

    % Ensure this is the current directory, important for pwd
    cd(current_file_path)


%% Local FLAGS
    localFLAG.Overwrite.MyData   = 0;

    localFLAG.PLOT.VaryQ   = 1;
    localFLAG.PLOT.VaryR   = 1;
    localFLAG.PLOT.VarySOC = 0;
    localFLAG.PLOT.VaryTs  = 0;

    Static.Q   = 1e-3; % [1e-3 1e-2 1e-1 1e0 1e1]
    Static.R   = 1e-1; % [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1]
    Static.SOC = 50;
    Static.Ts  = 1;
    
    colormapping = 'cool';
    % colormapping = 'copper';

    localFLAG.Plot.cell_voltage  = 1; % Terminal Voltage
    localFLAG.Plot.delta_phi     = 1; % Electrostatic potential difference between active material and electrolyte @AN/SEP interface
    localFLAG.Plot.eta           = 1; % Overpotential at SEI @AN/SEP interface
    localFLAG.Plot.C_Liion       = 1; % Concentration of Li^+ in the electrolyte @AN/SEP interface
    localFLAG.Plot.C_Li          = 1; % Concentration of Li at the surface of active material @AN/SEP interface
    localFLAG.Plot.delta_C_Li    = 1; % Difference of Li concentration between the surface and center of active material particle @AN/SEP interface

    localFLAG.xLim  = 3;
    % 1) Short
    % 2) Mid
    % 3) Long
    % 4) Custom
    xLim_custom = 1500;
    

%% Check if Data has already been consolidated
    data_filename = 'MyData.mat';
    ReadData = false;
    if isfile(data_filename)
        if localFLAG.Overwrite.MyData
            ReadData = true;
            delete(data_filename)
        end
    else
        ReadData = true;
    end
    % ReadData = true;


%% Read in Data
    if ReadData
        % Get Comparison filenames
            list = dir('*.mat');
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
                MyData(i).name          = sim_filenames{i};
                MyData(i).QQ            = sim.QQ;
                MyData(i).RR            = sim.RR;
                % MyData(i).SIM           = sim.SIM;
                MyData(i).FLAG          = sim.FLAG;
                MyData(i).N             = sim.N;
                MyData(i).P             = sim.P;
                % MyData(i).sys           = sim.sys;
                % MyData(i).P_k_pre_red   = sim.P_k_pre_red;
                % MyData(i).S_k_red       = sim.S_k_red;
                % MyData(i).K_k_red       = sim.K_k_red;
                % MyData(i).P_k_red       = sim.P_k_red;
                % MyData(i).CPCT_red      = sim.CPCT_red;
                % MyData(i).CPCT_pre_red  = sim.CPCT_pre_red;
                % MyData(i).CK_k_CV_red   = sim.CK_k_CV_red;
                % MyData(i).CK_k_red      = sim.CK_k_red;
                MyData(i).CPCT_pre_diag = sim.CPCT_pre_diag;
                MyData(i).CPCT_diag     = sim.CPCT_diag;
                % MyData(i).CK            = sim.CK;
                % MyData(i).CK_infty_CV   = sim.CK_infty_CV;
                % MyData(i).CK_infty      = sim.CK_infty;
                MyData(i).CP_inftyCT    = sim.CP_inftyCT;
                MyData(i).RESULTS       = sim.RESULTS;
            end
    
        % Save struct
            save(data_filename, 'MyData')
    else
        if FLAG.CLEARDATA
            load(data_filename);
        end
    end

    N_sims = length(MyData);


%% Get Desired Variable IDX
    IDX_DesVar = [];
    if localFLAG.Plot.cell_voltage
        IDX_DesVar = [IDX_DesVar , MyData(1).P.cell_voltage];
    end
    if localFLAG.Plot.delta_phi
        IDX_DesVar = [IDX_DesVar , MyData(1).P.delta_phi];
    end
    if localFLAG.Plot.eta
        IDX_DesVar = [IDX_DesVar , MyData(1).P.eta];
    end
    if localFLAG.Plot.C_Liion
        IDX_DesVar = [IDX_DesVar , MyData(1).P.C_Liion];
    end
    if localFLAG.Plot.C_Li
        IDX_DesVar = [IDX_DesVar , MyData(1).P.C_Li];
    end
    if localFLAG.Plot.delta_C_Li
        IDX_DesVar = [IDX_DesVar , MyData(1).P.delta_C_Li];
    end


%% Get IDX
    IDX_Q   = [];
    IDX_R   = [];
    IDX_SOC = [];
    IDX_Ts  = [];

    for j = 1:N_sims
        if MyData(j).QQ == Static.Q
            IDX_Q   = [IDX_Q , j];
        end
        if MyData(j).RR == Static.R
            IDX_R   = [IDX_R , j];
        end
        if MyData(j).FLAG.SOC == Static.SOC
            IDX_SOC = [IDX_SOC , j];
        end
        if MyData(j).FLAG.Ts == Static.Ts
            IDX_Ts  = [IDX_Ts , j];
        end
    end


%% Get colormap    
    figure(1);
    N_DataPoints = max(length(IDX_Q) , length(IDX_R));
    cmap         = colormap(colormapping);
    [r,c]        = size(cmap);
    IDX_Color    = linspace(1,r,N_DataPoints);
    IDX_Color    = floor(IDX_Color);
    close 1 % closes this figure


%% Get xLim
    switch localFLAG.xLim 
        case 1 % 1) Short
            xLim = MyData(1).N.steps_single;
        case 2 % 2) Mid
            xLim = MyData(1).N.steps_mid;
        case 3 % 3) Long
            xLim = MyData(1).N.steps_final;
        case 4
            xLim = xLim_custom;
    end


%% Plot Vary Q
if localFLAG.PLOT.VaryQ
    % Get combined IDX
    IDX = intersect(IDX_R,IDX_SOC);
    IDX = intersect(IDX  ,IDX_Ts );

    % Sort By Q
    for i = 1:length(IDX)
        Q_sort(i) = MyData(IDX(i)).QQ;
    end
    [~ , newI] = sort(Q_sort,'descend');
    IDX = IDX(newI);

    

    % Plot Results
    for j = IDX_DesVar
        figure
        hold on
        for i = 1:length(IDX)
            plot( MyData(IDX(i)).N.sample_vec , MyData(IDX(i)).CPCT_diag(j,:) , '-' , 'Color' , cmap(IDX_Color(i),:), 'LineWidth',2 , 'DisplayName' , ['Q = ' num2str(MyData(IDX(i)).QQ) ] )
        end
        ax = gca;
        ax.YScale = 'log';
        lgn = legend;
        lgn.Location = 'best';
        xlim([1 , xLim])
        title( [MyData(IDX(1)).RESULTS.Labels.title{j} ', SOC = ' num2str(MyData(IDX(1)).FLAG.SOC) ', T_s = ' num2str(MyData(IDX(1)).FLAG.Ts) , ', R = ' , num2str(MyData(IDX(1)).RR)])
        xlabel('Samples')
        ylabel('Covariance (CP_kC^T)')
    end
end


%% Plot Vary R
if localFLAG.PLOT.VaryR
    % Get combined IDX
    IDX = intersect(IDX_Q,IDX_SOC);
    IDX = intersect(IDX  ,IDX_Ts );

    % Sort By R
    for i = 1:length(IDX)
        R_sort(i) = MyData(IDX(i)).RR;
    end
    [~ , newI] = sort(R_sort,'descend');
    IDX = IDX(newI);

    % Plot Results
    for j = IDX_DesVar
        figure
        hold on
        for i = 1:length(IDX)
            plot( MyData(IDX(i)).N.sample_vec , MyData(IDX(i)).CPCT_diag(j,:) , '-' , 'Color' , cmap(IDX_Color(i),:), 'LineWidth',2 , 'DisplayName' , ['R = ' num2str(MyData(IDX(i)).RR) ] )
        end
        ax = gca;
        ax.YScale = 'log';
        lgn = legend;
        lgn.Location = 'best';
        xlim([1 , xLim])
        title( [MyData(IDX(1)).RESULTS.Labels.title{j}, ', SOC = ' , num2str(MyData(IDX(1)).FLAG.SOC) , ', T_s = ' , num2str(MyData(IDX(1)).FLAG.Ts) , ', Q = ' , num2str(MyData(IDX(1)).QQ) ])
        xlabel('Samples')
        ylabel('Covariance (CP_kC^T)')
    end
end


%% Plot Vary SOC
% if localFLAG.PLOT.VarySOC
%     % Get combined IDX
%     IDX = union(IDX_R,IDX_Q  );
%     IDX = union(IDX  ,IDX_Ts );
% 
%     % Get colormap
% 
%     % Plot Results
%     for j = IDX_DesVar
%         figure
%         hold on
%         for i = 1:length(IDX)
%             plot( MyData(i).N.sample_vec , MyData(i).CPCT_diag(j,:) , '-' , 'LineWidth',2 , 'DisplayName' , MyData(i).FLAG.SOC ) % , 'Color' , ColorVec(j,:)
%         end
%         lgn = legend;
%         title( [MyData(IDX(1)).RESULTS.Labels.title{j} ', Ts' num2str(MyData(IDX(1)).FLAG.Ts) , ', Q' , num2str(MyData(IDX(1)).QQ) , ', R' , num2str(MyData(IDX(1)).RR)])
%     end
% end
% 
% 
%% Plot Vary Ts
% if localFLAG.PLOT.VaryTs
%     % Get combined IDX
%     IDX = union(IDX_R,IDX_SOC);
%     IDX = union(IDX  ,IDX_Q  );
% 
%     % Get colormap
% 
%     % Plot Results
%     for j = IDX_DesVar
%         figure
%         hold on
%         for i = 1:length(IDX)
%             plot( MyData(i).N.sample_vec , MyData(i).CPCT_diag(j,:) , '-' , 'LineWidth',2 , 'DisplayName' , MyData(i).FLAG.Ts ) % , 'Color' , ColorVec(j,:)
%         end
%         lgn = legend;
%         title( [MyData(IDX(1)).RESULTS.Labels.title{j} ', SOC' num2str(MyData(IDX(1)).FLAG.SOC) ', Q' , num2str(MyData(IDX(1)).QQ) , ', R' , num2str(MyData(IDX(1)).RR)])
%     end
% end


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




%% OOOOOOOOOOOOOOOOOOOOOOOLLLLLLLLLLLLLLLLLLLLLLLDDDDDD
    
    % if localFLAG.Plot.NotNormalized
    %     cmin = -12;
    %     cmax = -7;
    % elseif localFLAG.Plot.Normalized2CellVoltage
    %     cmin = 0;
    %     cmax = 1;
    % elseif localFLAG.Plot.Normalized2MaxCPCT 
    %     cmin = -6;
    %     cmax =  0;
    % elseif localFLAG.Plot.Normalized2T10Delta 
    %     cmin = -12;
    %     cmax = -7;
    % end
