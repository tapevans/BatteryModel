%% Compare Data
% Is their a relationship between the singular value and P_infty

clear all; close all; clc;

%% FLAGS
    FLAG.READ_IN_DATA = 0;
    
    FLAG.ALLDATA      = 0;
    FLAG.ALLDATA_spec = 0;
    FLAG.NOISE        = 0;
    FLAG.NOISE_Cm     = 0; % plots for each Cm
    FLAG.ALLNOISE     = 0;
    FLAG.Ts           = 0;
    FLAG.C_Meas       = 0;
    FLAG.AngleComp    = 1;
    
    FLAG.UseNormCPC      = 0;
    FLAG.UseMAXNormCPC   = 1;
    
    FLAG.exportFigs   = 0;
        FLAG.FigAllData = 0;
        FLAG.FigNoise   = 0;
        FLAG.FigTs      = 0;
        FLAG.FigCm      = 0;
    

%% Read in data
if FLAG.READ_IN_DATA
    % Get filenames
        oldFolder = cd('F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanTesting\DAE_Circuit\SoftwareStyle\Results\ComparisonData');
        list = dir('*.mat*');
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
            MyData(i).C         = sim.C_idx;
            MyData(i).SVD       = sim.sing_val_norm;
            MyData(i).CPC       = sim.CPCT_calc;
            MyData(i).deg       = sim.matrixOFdot_deg;
            MyData(i).S_Orm     = sim.S_Orm;
            MyData(i).sing_val  = sim.sing_val;
        end

    % Save struct
        save('MyData.mat', 'MyData')

    %Go back to oldFolder
        cd(oldFolder);
else
    oldFolder = cd('F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanTesting\DAE_Circuit\SoftwareStyle\Results\ComparisonData');
    load('MyData.mat')
    cd(oldFolder);
end

%% Metrics of comparison
    Q_0_vec = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1]; % Process Noise
    R_0_vec = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1]; % Measurement Noise
%     Q_0_vec = [1e-6 1e-5 ]; % Process Noise
%     R_0_vec = [1e-3 1e-2 1e-1 1e0 1e1]; % Measurement Noise
    Ts_vec = [1e-7 5e-7 1e-6 5e-6 1e-5 1e-4];
    C_vec   = [1 2 3];

%% Create a high level idx
N_sims = length(MyData);

Q_idx = zeros(N_sims,1);
R_idx = zeros(N_sims,1);
T_idx = zeros(N_sims,1);
C_idx = zeros(N_sims,1);

for i = 1:N_sims
    for QQ = Q_0_vec
        if MyData(i).Q0 == QQ
            Q_idx(i) = 1;
        end
    end
    for RR = R_0_vec
        if MyData(i).R0 == RR
            R_idx(i) = 1;
        end
    end
    for TT = Ts_vec
        if MyData(i).Ts == TT
            T_idx(i) = 1;
        end
    end
    for CC = C_vec
        if MyData(i).C == CC
            C_idx(i) = 1;
        end
    end
end

combined_highLevel = Q_idx & R_idx & C_idx & T_idx;

%% Add Normalization of P_infty
    
    for i = 1:N_sims
        MyData(i).CPC_normalized = MyData(i).CPC / MyData(i).CPC(MyData(i).C);
        MyData(i).CPC_normMAX    = MyData(i).CPC / max(MyData(i).CPC);
    end

%% Add Marker Style data
    color_vec = {'g','r','b','c','m','#D95319'};
    markr_vec = {'o','square','^','hexagram'};
    
    for i = 1:N_sims
        MyData(i).shape = markr_vec{MyData(i).C};
        switch MyData(i).Ts
            case Ts_vec(1)
                MyData(i).color = color_vec{1};
            case Ts_vec(2)
                MyData(i).color = color_vec{2};
            case Ts_vec(3)
                MyData(i).color = color_vec{3};
            case Ts_vec(4)
                MyData(i).color = color_vec{4};
            case Ts_vec(5)
                MyData(i).color = color_vec{5};
            case Ts_vec(6)
                MyData(i).color = color_vec{6};
            otherwise
                MyData(i).color = 'w';
        end
    end

%% Make an ALL DATA Vector
if FLAG.ALLDATA
    svd_n_data_total = [];
    cpct_data_total  = [];
    % Store Data in vector for analysis
        for j = 1:N_sims
            svd_n_data_total = [svd_n_data_total, MyData(j).SVD];
            if FLAG.UseNormCPC
                cpct_data_total  = [cpct_data_total,  MyData(j).CPC_normalized ];
            elseif FLAG.UseMAXNormCPC
                cpct_data_total  = [cpct_data_total,  MyData(j).CPC_normMAX ];
            else
                cpct_data_total  = [cpct_data_total,  MyData(j).CPC ];
            end
            
        end

    % Calculate the covariance between the data
        [mu_x]    = calcMean( svd_n_data_total);
        [error_x] = calcError(svd_n_data_total , mu_x);
%         if FLAG.UseNormCPC
%             [mu_y]    = calcMean( cpct_data_total);
%             [error_y] = calcError(cpct_data_total , mu_y);
%         else
            [mu_y]    = calcMean( log10(cpct_data_total));
            [error_y] = calcError(log10(cpct_data_total) , mu_y);
%         end
        [Covar] = calcPopCovar(error_x, error_y);
    
    % Scatter plot of the data
        figure
        hold on
        if FLAG.UseNormCPC
            scatter(svd_n_data_total , log10(cpct_data_total),'filled','r','MarkerEdgeColor','k')
            yline(0)
        else
            scatter(svd_n_data_total , log10(cpct_data_total),'filled','r','MarkerEdgeColor','r')
            yline(0)
        end
        title(['All Data, Covariance: ' num2str(Covar)])
        if FLAG.UseNormCPC
            ylabel('log10(Error Covariance Normalized)')
        else
            ylabel('log10(Error Covariance)')
        end
        xlabel('Normalized Singular Values (Of SV)')
        if FLAG.exportFigs && FLAG.FigAllData
            exportgraphics(gca,'AllDataCorr_InputQ.png','Resolution',800)
        end
end


%% Make an ALL DATA Vector For specific Parameters
if FLAG.ALLDATA_spec
    Q_0_vec_spec = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1]; % Process Noise
    R_0_vec_spec = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1]; % Measurement Noise
    Ts_vec_spec  = [1e-7 5e-7 1e-6 5e-6 1e-5 1e-4];
    C_vec_spec   = [1 2 3];

    Q_idx = zeros(N_sims,1);
    R_idx = zeros(N_sims,1);
    C_idx = zeros(N_sims,1);
    T_idx = zeros(N_sims,1);
    for QQ = Q_0_vec_spec
        for RR = R_0_vec_spec
            for  CC = C_vec_spec
                for TT = Ts_vec_spec
                    for j = 1:N_sims
                        if MyData(j).Q0 == QQ
                            Q_idx(j) = 1;
                        end
                        if MyData(j).R0 == RR
                            R_idx(j) = 1;
                        end
                        if MyData(j).C == CC
                            C_idx(j) = 1;
                        end
                        if MyData(j).Ts == TT
                            T_idx(j) = 1;
                        end
                    end
                end
            end
        end
    end

    combined = Q_idx & R_idx & C_idx & T_idx;
    idx = find(combined == 1);

    idx_test = find((combined & combined_highLevel) == 1);

    svd_n_data = [];
    cpct_data  = [];
    % Store Data in vector for analysis
        for j = idx'
            svd_n_data = [svd_n_data, MyData(j).SVD];
            if FLAG.UseNormCPC
                cpct_data  = [cpct_data,  MyData(j).CPC_normalized ];
            elseif FLAG.UseMAXNormCPC
                cpct_data  = [cpct_data,  MyData(j).CPC_normMAX ];
            else
                cpct_data  = [cpct_data,  MyData(j).CPC ];
            end
            
        end

    % Calculate the covariance between the data
        [mu_x]    = calcMean( svd_n_data);
        [error_x] = calcError(svd_n_data , mu_x);
%         if FLAG.UseNormCPC
%             [mu_y]    = calcMean( cpct_data);
%             [error_y] = calcError(cpct_data , mu_y);
%         else
            [mu_y]    = calcMean( log10(cpct_data));
            [error_y] = calcError(log10(cpct_data) , mu_y);
%         end
        [Covar] = calcPopCovar(error_x, error_y);
    
    % Scatter plot of the data
        figure
        hold on
        if FLAG.UseNormCPC
            scatter(svd_n_data , log10(cpct_data),'filled','r','MarkerEdgeColor','k')
            yline(0)
        else
            scatter(svd_n_data , log10(cpct_data),'filled','r','MarkerEdgeColor','r')
            yline(0)
        end
        title(['All Data, Covariance: ' num2str(Covar)])
        if FLAG.UseNormCPC
            ylabel('log10(Error Covariance Normalized)')
        else
            ylabel('log10(Error Covariance)')
        end
        xlabel('Normalized Singular Values (Of SV)')
        if FLAG.exportFigs && FLAG.FigAllData
            exportgraphics(gca,'AllDataCorr_InputQ.png','Resolution',800)
        end
end


%% Sort by Noise  
if FLAG.NOISE
    Q_0 = 1e-6;
    R_0 = 1e1;
    Q_idx = zeros(N_sims,1);
    R_idx = zeros(N_sims,1);
    for j = 1:N_sims
        if MyData(j).Q0 == Q_0
            Q_idx(j) = 1;
        end
        if MyData(j).R0 == R_0
            R_idx(j) = 1;
        end
    end
    combined = Q_idx & R_idx;
%     idx = find(combined == 1);
    idx = find((combined & combined_highLevel) == 1);

    svd_n_data = [];
    cpct_data  = [];

    % Store Data in vector for analysis
    for j = idx'
        svd_n_data = [svd_n_data, MyData(j).SVD];
        if FLAG.UseNormCPC
            cpct_data  = [cpct_data,  MyData(j).CPC_normalized ];
        elseif FLAG.UseMAXNormCPC
            cpct_data  = [cpct_data,  MyData(j).CPC_normMAX ];
        else
            cpct_data  = [cpct_data,  MyData(j).CPC ];
        end
    end

% Calculate the covariance between the data
    [mu_x] = calcMean(svd_n_data);
    [error_x] = calcError(svd_n_data,mu_x);
%     if FLAG.UseNormCPC
%         [mu_y]    = calcMean( cpct_data);
%         [error_y] = calcError(cpct_data , mu_y);
%     else
        [mu_y]    = calcMean( log10(cpct_data));
        [error_y] = calcError(log10(cpct_data) , mu_y);
%     end
    [Covar] = calcPopCovar(error_x, error_y);

    % Scatter plot of the data (Colors are Ts, Cm are shapes)
        figure
        hold on
        for j = idx'
            if FLAG.UseNormCPC
                plot(MyData(j).SVD , log10(MyData(j).CPC_normalized),'LineStyle','none','MarkerSize',10,'Marker',MyData(j).shape,'MarkerFaceColor',MyData(j).color,'MarkerEdgeColor','k')
            elseif FLAG.UseMAXNormCPC
                plot(MyData(j).SVD , log10(MyData(j).CPC_normMAX),'LineStyle','none','MarkerSize',10,'Marker',MyData(j).shape,'MarkerFaceColor',MyData(j).color,'MarkerEdgeColor','k')
            else
                plot(MyData(j).SVD , log10(MyData(j).CPC),'LineStyle','none','MarkerSize',20    ,'Marker',MyData(j).shape,'MarkerFaceColor',MyData(j).color,'MarkerEdgeColor','k')
            end
        end
        yline(0)
        title(['Q = ' num2str(Q_0) '  R = ' num2str(R_0) '  Covariance: ' num2str(Covar) ])
        if FLAG.UseNormCPC
            ylabel('log10(Error Covariance Normalized)')
        else
            ylabel('log10(Error Covariance)')
        end
        xlabel('Normalized Singular Values (Of SV)')

        % Custom Legend
        for k = 1:length(color_vec)
            qw{k} = plot(nan,'LineStyle','none' ,'Marker','o','MarkerFaceColor',color_vec{k} ,'MarkerEdgeColor','k', 'DisplayName',[num2str(Ts_vec(k)) ' s']);
        end
        for k = length(color_vec)+1:length(color_vec)+length(markr_vec)
            offset = length(color_vec);
            qw{k} = plot(nan,'LineStyle','none', 'Marker',markr_vec{k-offset},'MarkerFaceColor','k','MarkerEdgeColor','k', 'DisplayName',['C_m = ' num2str(C_vec(k-offset))]);
        end
        lgn = legend([qw{:}]);
        lgn.NumColumns = 2;
        lgn.Orientation = 'horizontal';
        lgn.Location    = 'northeast';

        if FLAG.exportFigs && FLAG.FigNoise
            exportgraphics(gca,['SingleNoiseCorr' '_Q' num2str(Q_0) '_R' num2str(R_0) '_InputQ.png'],'Resolution',800)
        end
end


%% Sort by Noise and Plot by Cm
if FLAG.NOISE_Cm
    Q_0 = 1e-6;
    R_0 = 1e1;
    for  i = C_vec
        Q_idx = zeros(N_sims,1);
        R_idx = zeros(N_sims,1);
        C_idx = zeros(N_sims,1);
        for j = 1:N_sims
            if MyData(j).Q0 == Q_0
                Q_idx(j) = 1;
            end
            if MyData(j).R0 == R_0
                R_idx(j) = 1;
            end
            if MyData(j).C == i
                C_idx(j) = 1;
            end
        end
        combined = Q_idx & R_idx & C_idx;
%         idx = find(combined == 1);
        idx = find((combined & combined_highLevel) == 1);

        svd_n_data = [];
        cpct_data  = [];

        % Store Data in vector for analysis
        for j = idx'
            svd_n_data = [svd_n_data, MyData(j).SVD];
            if FLAG.UseNormCPC
                cpct_data  = [cpct_data,  MyData(j).CPC_normalized ];
            elseif FLAG.UseMAXNormCPC
                cpct_data  = [cpct_data,  MyData(j).CPC_normMAX ];
            else
                cpct_data  = [cpct_data,  MyData(j).CPC ];
            end
        end

        % Calculate the covariance between the data
        [mu_x] = calcMean(svd_n_data);
        [error_x] = calcError(svd_n_data,mu_x);
%         if FLAG.UseNormCPC
%             [mu_y]    = calcMean( cpct_data);
%             [error_y] = calcError(cpct_data , mu_y);
%         else
            [mu_y]    = calcMean( log10(cpct_data));
            [error_y] = calcError(log10(cpct_data) , mu_y);
%         end
        [Covar] = calcPopCovar(error_x, error_y);

        % Scatter plot of the data (Colors are Ts, Cm are shapes)
        figure
        hold on
        for j = idx'
            if FLAG.UseNormCPC
                plot(MyData(j).SVD , log10(MyData(j).CPC_normalized),'LineStyle','none','MarkerSize',10,'Marker',MyData(j).shape,'MarkerFaceColor',MyData(j).color,'MarkerEdgeColor','k')
            elseif FLAG.UseMAXNormCPC
                plot(MyData(j).SVD , log10(MyData(j).CPC_normMAX),'LineStyle','none','MarkerSize',10,'Marker',MyData(j).shape,'MarkerFaceColor',MyData(j).color,'MarkerEdgeColor','k')
            else
                plot(MyData(j).SVD , log10(MyData(j).CPC),'LineStyle','none','MarkerSize',20    ,'Marker',MyData(j).shape,'MarkerFaceColor',MyData(j).color,'MarkerEdgeColor','k')
            end
        end
        yline(0)
        title(['Q = ' num2str(Q_0) '  R = ' num2str(R_0) '  C_m = ' num2str(i) '   Covariance: ' num2str(Covar) ])
        if FLAG.UseNormCPC
            ylabel('log10(Error Covariance Normalized)')
        else
            ylabel('log10(Error Covariance)')
        end
        xlabel('Normalized Singular Values (Of SV)')
        
        % Custom Legend
        qw = {};
        for k = 1:length(color_vec)
            qw{k} = plot(nan,'LineStyle','none' ,'Marker','o','MarkerFaceColor',color_vec{k} ,'MarkerEdgeColor','k', 'DisplayName',[num2str(Ts_vec(k)) ' s']);
        end
%         for k = length(color_vec)+1:length(color_vec)+length(markr_vec)
%             offset = length(color_vec);
%             qw{k} = plot(nan,'LineStyle','none', 'Marker',markr_vec{k-offset},'MarkerFaceColor','k','MarkerEdgeColor','k', 'DisplayName',['C_m = ' num2str(C_vec(k-offset))]);
%         end
        lgn = legend([qw{:}]);
        lgn.NumColumns = 1;
        lgn.Orientation = 'horizontal';
        lgn.Location    = 'northeast';

        if FLAG.exportFigs && FLAG.FigNoise
            exportgraphics(gca,['SingleNoiseCmCorr' '_Q' num2str(Q_0) '_R' num2str(R_0) '_C' num2str(i) '_InputQ.png'],'Resolution',800)
        end
    end
end

%% Sort by All Noise  
if FLAG.ALLNOISE
    covar_out = [];
    for QQ = Q_0_vec
        for RR = R_0_vec 
            Q_0 = QQ;
            R_0 = RR;
            Q_idx = zeros(N_sims,1);
            R_idx = zeros(N_sims,1);
            for j = 1:N_sims
                if MyData(j).Q0 == Q_0
                    Q_idx(j) = 1;
                end
                if MyData(j).R0 == R_0
                    R_idx(j) = 1;
                end
            end
            combined = Q_idx & R_idx;
%             idx = find(combined == 1);
            idx = find((combined & combined_highLevel) == 1);

            % Store Data in vector for analysis
                svd_n_data = [];
                cpct_data  = [];
                for j = idx'
                    svd_n_data = [svd_n_data, MyData(j).SVD];
                    if FLAG.UseNormCPC
                        cpct_data  = [cpct_data,  MyData(j).CPC_normalized ];
                    elseif FLAG.UseMAXNormCPC
                        cpct_data  = [cpct_data,  MyData(j).CPC_normMAX ];
                    else
                        cpct_data  = [cpct_data,  MyData(j).CPC ];
                    end
                end

            % Calculate the covariance between the data
        [mu_x] = calcMean(svd_n_data);
        [error_x] = calcError(svd_n_data,mu_x);
%         if FLAG.UseNormCPC
%             [mu_y]    = calcMean( cpct_data);
%             [error_y] = calcError(cpct_data , mu_y);
%         else
            [mu_y]    = calcMean( log10(cpct_data));
            [error_y] = calcError(log10(cpct_data) , mu_y);
%         end
        [Covar] = calcPopCovar(error_x, error_y);

        % Scatter plot of the data (Colors are Ts, Cm are shapes)
                figure
                hold on
                for j = idx'
                    if FLAG.UseNormCPC
                        plot(MyData(j).SVD , log10(MyData(j).CPC_normalized),'LineStyle','none','MarkerSize',10,'Marker',MyData(j).shape,'MarkerFaceColor',MyData(j).color,'MarkerEdgeColor','k')
                    elseif FLAG.UseMAXNormCPC
                        plot(MyData(j).SVD , log10(MyData(j).CPC_normMAX),'LineStyle','none','MarkerSize',10,'Marker',MyData(j).shape,'MarkerFaceColor',MyData(j).color,'MarkerEdgeColor','k')
                    else
                        plot(MyData(j).SVD , log10(MyData(j).CPC),'LineStyle','none','MarkerSize',10,'Marker',MyData(j).shape,'MarkerFaceColor',MyData(j).color,'MarkerEdgeColor','k')
                    end
                end
                yline(0)
                title(['Q = ' num2str(Q_0) '  R = ' num2str(R_0) '  Covariance: ' num2str(Covar) ])
                if FLAG.UseNormCPC
                    ylabel('log10(Error Covariance Normalized)')
                else
                    ylabel('log10(Error Covariance)')
                end
                xlabel('Normalized Singular Values (Of SV)')
%                 if FLAG.exportFigs
%                     %exportgraphics(gca,'All_NoiseCorr.png','Resolution',800)
%                 end

            covar_out(end+1,1) = Covar;
        end
    end
    covar_out = reshape(covar_out,length(R_0_vec),length(Q_0_vec));
end

%% Sort by Ts
if FLAG.Ts
    covar_out_Ts = [];
    for i = Ts_vec
        Ts_idx = zeros(N_sims,1);
        for j = 1:N_sims
            if MyData(j).Ts == i
                Ts_idx(j) = 1;
            end
        end
%         idx = find(Ts_idx == 1);
        idx = find((Ts_idx & combined_highLevel) == 1);

        % Store Data in vector for analysis
        svd_n_data = [];
        cpct_data  = [];
        for j = idx'
            svd_n_data = [svd_n_data, MyData(j).SVD];
            if FLAG.UseNormCPC
                cpct_data  = [cpct_data,  MyData(j).CPC_normalized ];
            elseif FLAG.UseMAXNormCPC
                cpct_data  = [cpct_data,  MyData(j).CPC_normMAX ];
            else
                cpct_data  = [cpct_data,  MyData(j).CPC ];
            end
        end

        % Calculate the covariance between the data
        [mu_x] = calcMean(svd_n_data);
        [error_x] = calcError(svd_n_data,mu_x);
%         if FLAG.UseNormCPC
%             [mu_y]    = calcMean( cpct_data);
%             [error_y] = calcError(cpct_data , mu_y);
%         else
            [mu_y]    = calcMean( log10(cpct_data));
            [error_y] = calcError(log10(cpct_data) , mu_y);
%         end
        [Covar] = calcPopCovar(error_x, error_y);

        % Scatter plot of the data (Colors are Ts, Cm are shapes)
        figure
        hold on
        for j = idx'
            if FLAG.UseNormCPC
                plot(MyData(j).SVD , log10(MyData(j).CPC_normalized),'LineStyle','none','MarkerSize',10,'Marker',MyData(j).shape,'MarkerFaceColor',MyData(j).color,'MarkerEdgeColor','k')
            elseif FLAG.UseMAXNormCPC
                plot(MyData(j).SVD , log10(MyData(j).CPC_normMAX),'LineStyle','none','MarkerSize',10,'Marker',MyData(j).shape,'MarkerFaceColor',MyData(j).color,'MarkerEdgeColor','k')
            else
                plot(MyData(j).SVD , log10(MyData(j).CPC),'LineStyle','none','MarkerSize',10,'Marker',MyData(j).shape,'MarkerFaceColor',MyData(j).color,'MarkerEdgeColor','k')
            end
        end
        yline(0)
        title(['T_s = ' num2str(i) '  Covariance: ' num2str(Covar) ])
        if FLAG.UseNormCPC
            ylabel('log10(Error Covariance Normalized)')
        else
            ylabel('log10(Error Covariance)')
        end
        xlabel('Normalized Singular Values (Of SV)')
        if FLAG.exportFigs && FLAG.FigTs
            exportgraphics(gca,['TsCorr' '_Ts' num2str(i) '_InputQ.png'],'Resolution',800)
        end

        covar_out_Ts(end+1,1) = Covar;
    end
end


%% Sort by Measured State
if FLAG.C_Meas
    covar_out_C_meas = [];
    for i = C_vec
        C_idx = zeros(N_sims,1);
        for j = 1:N_sims
            if MyData(j).C == i
                C_idx(j) = 1;
            end
        end
%         idx = find(C_idx == 1);
        idx = find((C_idx & combined_highLevel) == 1);

        % Store Data in vector for analysis
        svd_n_data = [];
        cpct_data  = [];
        for j = idx'
            svd_n_data = [svd_n_data, MyData(j).SVD];
            if FLAG.UseNormCPC
                cpct_data  = [cpct_data,  MyData(j).CPC_normalized ];
            elseif FLAG.UseMAXNormCPC
                cpct_data  = [cpct_data,  MyData(j).CPC_normMAX ];
            else
                cpct_data  = [cpct_data,  MyData(j).CPC ];
            end
        end

        % Calculate the covariance between the data
        [mu_x] = calcMean(svd_n_data);
        [error_x] = calcError(svd_n_data,mu_x);
%         if FLAG.UseNormCPC
%             [mu_y]    = calcMean( cpct_data);
%             [error_y] = calcError(cpct_data , mu_y);
%         else
            [mu_y]    = calcMean( log10(cpct_data));
            [error_y] = calcError(log10(cpct_data) , mu_y);
%         end
        [Covar] = calcPopCovar(error_x, error_y);

        % Scatter plot of the data (Colors are Ts, Cm are shapes)
        figure
        hold on
        for j = idx'
            if FLAG.UseNormCPC
                plot(MyData(j).SVD , log10(MyData(j).CPC_normalized),'LineStyle','none','MarkerSize',10,'Marker',MyData(j).shape,'MarkerFaceColor',MyData(j).color,'MarkerEdgeColor','k')
            elseif FLAG.UseMAXNormCPC
                plot(MyData(j).SVD , log10(MyData(j).CPC_normMAX),'LineStyle','none','MarkerSize',10,'Marker',MyData(j).shape,'MarkerFaceColor',MyData(j).color,'MarkerEdgeColor','k')
            else
                plot(MyData(j).SVD , log10(MyData(j).CPC)    ,'LineStyle','none','MarkerSize',10,'Marker',MyData(j).shape,'MarkerFaceColor',MyData(j).color,'MarkerEdgeColor','k')
            end
        end
        yline(0)
        title(['C_m = ' num2str(i) '  Covariance: ' num2str(Covar) ])
        if FLAG.UseNormCPC
            ylabel('log10(Error Covariance Normalized)')
        else
            ylabel('log10(Error Covariance)')
        end
        xlabel('Normalized Singular Values (Of SV)')
        if FLAG.exportFigs  && FLAG.FigCm
            exportgraphics(gca,['Cm' num2str(i) 'Corr.png'],'Resolution',800)
        end

        covar_out_C_meas(end+1,1) = Covar;
    end
end


%% Angle Comparison
if FLAG.AngleComp
    makeAngleComparison(MyData,FLAG)
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




%% Functions
function [mu] = calcMean(x)
mu = mean(x);
end
function [error] = calcError(x,mu)
% difference between x and the expected value (mu)
error = x-mu;
end
function [Covar] = calcPopCovar(error_x, error_y)
Covar = (error_x*error_y')/length(error_x);
end
