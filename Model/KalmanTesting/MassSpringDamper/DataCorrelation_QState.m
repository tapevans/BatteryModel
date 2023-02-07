%% Data Correlation
% Using QStates data
clear all; close all; clc;

%% FLAGS
FLAG.READ_IN_DATA = 0;

FLAG.ALLDATA      = 0;
FLAG.NOISE        = 1;
FLAG.Ts           = 0;
FLAG.C_Meas       = 0;

%% Read in data
if FLAG.READ_IN_DATA
    % Get filenames
    oldFolder = cd([pwd filesep 'Step_CodedPlantResultsQState']);
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
        %     if ~exist('MyData')
        %         MyData(i).name = sim_filenames{i};
        %         MyData(i).Q0  = Q_0;
        %         MyData(i).R0  = R_0;
        %         MyData(i).Ts  = Ts;
        %         MyData(i).C   = C_idx;
        %         MyData(i).SVD = sing_val_norm;
        %         MyData(i).CPC = CPCT_calc;
        %     else
        MyData(i).name = sim_filenames{i};
        MyData(i).Q0  = sim.Q_0;
        MyData(i).R0  = sim.R_0;
        MyData(i).Ts  = sim.Ts;
        MyData(i).C   = sim.C_idx;
        MyData(i).SVD = sim.sing_val_norm;
        MyData(i).CPC = sim.CPCT_calc;
        %     end
    end

    % Save struct
    save('MyData.mat', 'MyData')

    %Go back to oldFolder
    cd(oldFolder);
else
    oldFolder = cd([pwd filesep 'Step_CodedPlantResultsQState']);
    load('MyData.mat')
    cd(oldFolder);
end
%% Sort how to compare
Q_0_vec = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1]; % Process Noise
R_0_vec = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1]; % Measurement Noise
Ts_vec  = [.01 .05 .1 .5 1 5]; %[s]
C_vec   = [1 2 3 4];

N_sims = length(MyData);

%% Make an ALL DATA Vector
if FLAG.ALLDATA
    svd_n_data_total = [];
    cpct_data_total  = [];
    for j = 1:N_sims
        svd_n_data_total = [svd_n_data_total, MyData(j).SVD];
        cpct_data_total  = [cpct_data_total,  MyData(j).CPC ];
    end

    figure
    scatter(svd_n_data_total , log10(cpct_data_total))
    title('All Data')
    xlabel('Normalized Singular Values (Of SV)')
    ylabel('log10(Error Covariance)')

    [mu_x] = calcMean(svd_n_data_total);
    [error_x] = calcError(svd_n_data_total,mu_x);
    [mu_y] = calcMean(log10(cpct_data_total));
    [error_y] = calcError(log10(cpct_data_total),mu_y);
    [Covar_total] = calcPopCovar(error_x, error_y)

end

%% Sort by Noise
% Q_0 = 1e-6;
% R_0 = 1e-6;
% Q_idx = zeros(N_sims,1);
% R_idx = zeros(N_sims,1);
% for j = 1:N_sims
%     if MyData(j).Q0 == Q_0
%         Q_idx(j) = 1;
%     end
%     if MyData(j).R0 == R_0
%         R_idx(j) = 1;
%     end
% end
% combined = Q_idx & R_idx;
% idx = find(combined == 1);
%
% svd_n_data = [];
% cpct_data  = [];
% for j = idx
%     svd_n_data = [svd_n_data, MyData(j).SVD];
%     cpct_data  = [cpct_data,  MyData(j).CPC ];
% end
%
% [mu_x] = calcMean(svd_n_data);
% [error_x] = calcError(svd_n_data,mu_x);
% [mu_y] = calcMean(log10(cpct_data));
% [error_y] = calcError(log10(cpct_data),mu_y);
% [Covar] = calcPopCovar(error_x, error_y);
%
% figure
% scatter(svd_n_data , log10(cpct_data),'fill','k')
% title(['Q = ' num2str(Q_0) '  R = ' num2str(R_0) '  Covariance: ' num2str(Covar) ])
% xlabel('Normalized Singular Values (Of SV)')
% ylabel('log10(Error Covariance)')

if FLAG.NOISE
    covar_out = [];
    for QQ = Q_0_vec
        for RR = R_0_vec % Q_0_vec
            %         Q_0 = 1e-6;
            %         R_0 = 1e-6;
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
            idx = find(combined == 1);

            svd_n_data = [];
            cpct_data  = [];
            for j = idx
                svd_n_data = [svd_n_data, MyData(j).SVD];
                cpct_data  = [cpct_data,  MyData(j).CPC ];
            end

            [mu_x] = calcMean(svd_n_data);
            [error_x] = calcError(svd_n_data,mu_x);
            [mu_y] = calcMean(log10(cpct_data));
            [error_y] = calcError(log10(cpct_data),mu_y);
            [Covar] = calcPopCovar(error_x, error_y);

            figure
            scatter(svd_n_data , log10(cpct_data),'filled','r','MarkerEdgeColor','k')
            title(['Q = ' num2str(Q_0) '  R = ' num2str(R_0) '  Covariance: ' num2str(Covar) ])
            xlabel('Normalized Singular Values (Of SV)')
            ylabel('log10(Error Covariance)')

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
        idx = find(Ts_idx == 1);

        svd_n_data = [];
        cpct_data  = [];
        for j = idx
            svd_n_data = [svd_n_data, MyData(j).SVD];
            cpct_data  = [cpct_data,  MyData(j).CPC ];
        end

        [mu_x] = calcMean(svd_n_data);
        [error_x] = calcError(svd_n_data,mu_x);
        [mu_y] = calcMean(log10(cpct_data));
        [error_y] = calcError(log10(cpct_data),mu_y);
        [Covar] = calcPopCovar(error_x, error_y);

        figure
        scatter(svd_n_data , log10(cpct_data),'fill','k')
        title(['T_s = ' num2str(i) '  Covariance: ' num2str(Covar) ])
        xlabel('Normalized Singular Values (Of SV)')
        ylabel('log10(Error Covariance)')

        covar_out_Ts(end+1,1) = Covar;
    end
end


%% Sort by Measured
if FLAG.C_Meas
    covar_out_C_meas = [];
    for i = C_vec
        C_idx = zeros(N_sims,1);
        for j = 1:N_sims
            if MyData(j).C == i
                C_idx(j) = 1;
            end
        end
        idx = find(C_idx == 1);

        svd_n_data = [];
        cpct_data  = [];
        for j = idx
            svd_n_data = [svd_n_data, MyData(j).SVD];
            cpct_data  = [cpct_data,  MyData(j).CPC ];
        end

        [mu_x] = calcMean(svd_n_data);
        [error_x] = calcError(svd_n_data,mu_x);
        [mu_y] = calcMean(log10(cpct_data));
        [error_y] = calcError(log10(cpct_data),mu_y);
        [Covar] = calcPopCovar(error_x, error_y);

        figure
        scatter(svd_n_data , log10(cpct_data),'fill','k')
        title(['C_m = ' num2str(i) '  Covariance: ' num2str(Covar) ])
        xlabel('Normalized Singular Values (Of SV)')
        ylabel('log10(Error Covariance)')

        covar_out_C_meas(end+1,1) = Covar;
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

