%% Plot Single Estimation Covariance Matrix Analysis
    clear all; close all; clc;


%% Switch to current directory
    [current_file_path,~,~] = fileparts(mfilename('fullpath'));

    % Ensure this is the current directory, important for pwd
    cd(current_file_path)


%% Load File
    QQ = 1e-3;    RR = 1e-1;    SOC = 50;    Ts = 1;

    filename = ['CovarEst_Combined_SOC' num2str(SOC) '_Ts' num2str(Ts) '_Q' num2str(QQ) '_R' num2str(RR) '.mat'];
    
    load(filename)


%% Parameters
    FLAG.PLOT.KALMAN_GAIN                   = 0;
        FLAG.PLOT.KALMAN.CV         = 0;
        FLAG.PLOT.KALMAN.AllDes     = 1;
        FLAG.KalmanLim  = 4;
            xLim_custom = 1500;
            % 1) Short
            % 2) Mid
            % 3) Long
            % 4) Custom
            switch FLAG.KalmanLim 
                case 1
                    KalmanLim = N.steps_single;
                case 2
                    KalmanLim = N.steps_mid;
                case 3
                    KalmanLim = N.steps_final;
                case 4
                    KalmanLim = xLim_custom;
            end
        FLAG.PLOT.KALMAN.AllDes_olap= 1; % Overlapping
        FLAG.PLOT.KALMAN.Logy   = 0;

    FLAG.PLOT.COVAR                         = 1;
    % Overlapping
    FLAG.PLOT.COVAR_short           = 1;
    FLAG.PLOT.COVAR_mid             = 1;
    FLAG.PLOT.COVAR_long            = 1;
        FLAG.PLOT.SUB           = 1;

    % Individual Variables
    FLAG.PLOT.COVAR_short_desVar    = 0;
    FLAG.PLOT.COVAR_mid_desVar      = 0;
    FLAG.PLOT.COVAR_long_desVar     = 0;
        FLAG.UseLogY            = 1;
        FLAG.UseCPC_pre         = 1;

    FLAG.WithAsy    = 1;

    ColorVec = [0      0.4470 0.7410
                0.8500 0.3250 0.0980
                0.4940 0.1840 0.5560
                0.4660 0.6740 0.1880
                0.3010 0.7450 0.9330
                0.6350 0.0780 0.1840];


%% Kalman Gain
if FLAG.PLOT.KALMAN_GAIN
    % Kalman Gain (CK) of just voltage
    if FLAG.PLOT.KALMAN.CV        
        figure
        plot( N.sample_vec , CK_CV , '-b' , 'LineWidth',2 )
        if FLAG.WithAsy
            yL = yline(CK_infty_CV   , 'k--', 'LineWidth',2);
            yL.Alpha = 1;
        end
        title([RESULTS.Labels.title{P.cell_voltage} ' Kalman Gain'])
        xlabel('Samples')
        ylabel('C K_k')
        xlim([1 , KalmanLim])
    end

    % Kalman Gain (CK) of all desired variables
    if FLAG.PLOT.KALMAN.AllDes 
        for j = 1:N.measur_all
            figure
            % plot( N.sample_vec , reshape(CK_k_red(j,:,:) , 1 , []) , '-b' , 'LineWidth',2 )
            plot( N.sample_vec , CK(j,:) , '-b' , 'LineWidth',2 )
            if FLAG.WithAsy
                yL = yline(CK_infty(j,1)   , 'k--', 'LineWidth',2);
                yL.Alpha = 1;
            end
            title([RESULTS.Labels.title{j} ' Kalman Gain'])
            xlabel('Samples')
            ylabel('C K_k')
            xlim([1 , KalmanLim])
        end
    end

    % Kalman Gain (CK) of all desired variables overlapping
    if FLAG.PLOT.KALMAN.AllDes_olap 
        figure
        hold on
        for j = 1:N.measur_all
            plot( N.sample_vec , reshape(CK_k_red(j,:,:) , 1 , []) , '-b' , 'LineWidth',2, 'Color' , ColorVec(j,:) , 'DisplayName' , RESULTS.Labels.title{j}  )
        end
        title('Kalman Gain')
        xlabel('Samples')
        ylabel('C K_k')
        xlim([1 , KalmanLim])
        lgn = legend;
        lgn.Location = 'northwest';

        if FLAG.PLOT.KALMAN.Logy
            ax = gca;
            ax.YScale = 'log';
        end
    end
end


%% Output Covariance Matrix
if FLAG.PLOT.COVAR
    % Covariance (CPCT) of each variable for the short time, overlapped
    if FLAG.PLOT.COVAR_short
        figure
        hold on
        for j = 1:N.measur_all
            plot( N.sample_vec , CPCT_diag(j,:) , '-' , 'LineWidth',2, 'Color' , ColorVec(j,:) , 'DisplayName' , RESULTS.Labels.title{j} )
        end
        lgn = legend;
        title(['CP_kC^T with Asymptotic' , '  Q = ' num2str(QQ) , '  R = ' num2str(RR) ])
        xlabel('Samples')
        ylabel('Covariance')
        % lgn.Location = 'eastoutside';
        lgn.Location = 'northwest';
        xlim([1 , N.steps_single])
        if FLAG.PLOT.SUB
            ax2 = axes('Position',[0.5 0.2 .35 .45]);
            box on;
            hold on
            for j = 1:N.measur_all
                plot( N.sample_vec , CPCT_diag(j,:) , '-' , 'LineWidth',2, 'Color' , ColorVec(j,:) , 'DisplayName' , RESULTS.Labels.title{j} )
            end
            xlim([1 , N.steps_single])
            ylim([0 , 2e-8])
            ax2.XTick =  ax2.XTick(1:2:5);
        end
    end

    % Covariance (CPCT) of each variable for the mid time, overlapped
    if FLAG.PLOT.COVAR_mid
        figure
        hold on
        for j = 1:N.measur_all
            plot( N.sample_vec , CPCT_diag(j,:) , '-' , 'LineWidth',2, 'Color' , ColorVec(j,:) , 'DisplayName' , RESULTS.Labels.title{j} )
        end
        lgn = legend;
        title(['CP_kC^T with Asymptotic' , '  Q = ' num2str(QQ) , '  R = ' num2str(RR) ])
        xlabel('Samples')
        ylabel('Covariance')
        lgn.Location = 'southwest';
        xlim([1 , N.steps_mid])
        if FLAG.PLOT.SUB
            ax2 = axes('Position',[0.395 0.4 .35 .45]);
            box on;
            hold on
            for j = 1:N.measur_all
                plot( N.sample_vec , CPCT_diag(j,:) , '-' , 'LineWidth',2, 'Color' , ColorVec(j,:) , 'DisplayName' , RESULTS.Labels.title{j} )
            end
            xlim([2000 , N.steps_mid])
            ylim([0 , 1e-7])
        end
    end

    % Covariance (CPCT) of each variable for the long time, overlapped
    if FLAG.PLOT.COVAR_long
        figure
        hold on
        for j = 1:N.measur_all
            plot( N.sample_vec , CPCT_diag(j,:) , '-' , 'LineWidth',2, 'Color' , ColorVec(j,:) , 'DisplayName' , RESULTS.Labels.title{j} )
        end
        lgn = legend;
        title(['CP_kC^T with Asymptotic' , '  Q = ' num2str(QQ) , '  R = ' num2str(RR) ])
        xlabel('Samples')
        ylabel('Covariance')
        lgn.Location = 'east';
        xlim([1 , N.steps_final])
        if FLAG.PLOT.SUB
            ax2 = axes('Position',[0.2 0.22 .35 .45]);
            box on;
            hold on
            for j = 1:N.measur_all
                plot( N.sample_vec , CPCT_diag(j,:) , '-' , 'LineWidth',2, 'Color' , ColorVec(j,:) , 'DisplayName' , RESULTS.Labels.title{j} )
            end
            ax2.YScale = 'log';
            ax2.XLim = [1  , N.steps_final];
            % ax2.XLim = [1e4  , N.steps_final];
            % ax2.XLim = [1e4  , 4e4];
            ax2.YLim = [1e-11 , 1e-4];
        end
    end

    % Covariance (CPCT) of each variable for the short time, individual
    if FLAG.PLOT.COVAR_short_desVar
        for j = 1:N.measur_all
            figure
            hold on        
            if FLAG.UseCPC_pre
                plot( N.sample_vec , CPCT_pre_diag(j,:) , '-' , 'LineWidth',2)
            else
                plot( N.sample_vec , CPCT_diag(j,:)     , '-' , 'LineWidth',2)
            end
            % if FLAG.UseLogY
            %     ax = gca;
            %     ax.YScale = 'log';
            % end
            if FLAG.WithAsy
                yL = yline(CP_inftyCT(j,j)              , '--', 'LineWidth',2);
                yL.Alpha = 1;
            end
            title([RESULTS.Labels.title{j} , '  Q = ' num2str(QQ) , '  R = ' num2str(RR) ])
            xlabel('Samples')
            ylabel('Covariance')
            xlim([1 , N.steps_single])
        end
    end

    % Covariance (CPCT) of each variable for the mid time, individual
    if FLAG.PLOT.COVAR_mid_desVar
        for j = 1:N.measur_all
            figure
            hold on        
            if FLAG.UseCPC_pre
                plot( N.sample_vec , CPCT_pre_diag(j,:) , '-' , 'LineWidth',2)
            else
                plot( N.sample_vec , CPCT_diag(j,:)     , '-' , 'LineWidth',2)
            end
            % if FLAG.UseLogY
            %     ax = gca;
            %     ax.YScale = 'log';
            % end
            if FLAG.WithAsy
                yL = yline(CP_inftyCT(j,j)              , '--', 'LineWidth',2);
                yL.Alpha = 1;
            end
            title([RESULTS.Labels.title{j} , '  Q = ' num2str(QQ) , '  R = ' num2str(RR) ])
            xlabel('Samples')
            ylabel('Covariance')
            xlim([1 , N.steps_mid])
        end
    end

    % Covariance (CPCT) of each variable for the long time, individual
    if FLAG.PLOT.COVAR_long_desVar
        for j = 1:N.measur_all
            figure
            hold on        
            if FLAG.UseCPC_pre
                plot( N.sample_vec , CPCT_pre_diag(j,:) , '-' , 'LineWidth',2);
            else
                plot( N.sample_vec , CPCT_diag(j,:)     , '-' , 'LineWidth',2);
            end
            if FLAG.UseLogY
                ax = gca;
                ax.YScale = 'log';
            end
            if FLAG.WithAsy
                yL = yline(CP_inftyCT(j,j) , '--','Asymptotic', 'LineWidth',2);
                yL.Alpha = 1;
            end
            title([RESULTS.Labels.title{j} , '  Q = ' num2str(QQ) , '  R = ' num2str(RR) ])
            xlabel('Samples')
            ylabel('Covariance')
            xlim([1 , N.steps_final])
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