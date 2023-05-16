%% Perform Estimation Covariance Matrix Analysis

function [] = CovarEstimationAnalysis(SIM,FLAG,P,sys_HK,RESULTS)
%% Other Parameters
    FLAG.folderpathCovarEst = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanFilter\Results\CovarMatrixAnalysis';

    FLAG.OverwriteResults = 0;
    FLAG.SAVEResults      = 1;
    FLAG.PLOT.Any         = 0;
    FLAG.PLOT.KALMAN_GAIN = 0;
    FLAG.PLOT.COVAR       = 1;
        FLAG.PLOT.SUB     = 1;

    ColorVec = [0      0.4470 0.7410
                0.8500 0.3250 0.0980
                0.4940 0.1840 0.5560
                0.4660 0.6740 0.1880
                0.3010 0.7450 0.9330
                0.6350 0.0780 0.1840];


%% Number of Steps
    N.steps_single = 3000; % Samples plotted every 1 sample, after this only every N.interval samples are saved
    N.steps_mid    = 10000; % For plotting
    % N.steps_final  = 100000; % Last sample
    N.steps_final  = 4e4; % Last sample
    N.interval     = 100;
    N.sample_vec     = [1:1:N.steps_single , N.steps_single+N.interval : N.interval : N.steps_final];
    N.steps = length(N.sample_vec);


%% System Variables
    sys = sys_HK{end}; % This is hardcoded for combined system for right now

    A_DT    = sys.A;
    B_DT    = sys.B;
    C_DT    = sys.C;
    C_DT_CV = sys.C(P.cell_voltage,:);
    
    [N.measur, N.states] = size(C_DT_CV);
    [N.measur_all, ~]    = size(C_DT);


%% Asymptotic Calculations
    [K_infty, P_infty] = AsymptoticPreCalcs(FLAG,SIM,sys);
     K_infty           = K_infty(:,P.cell_voltage);


%% Noise Matrix
    Q_vec = [1e-3 1e-2 1e-1 1e0 1e1];
    R_vec = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1];


%% Initialize Variables
    % S_k_0      = zeros(N.measur     , N.measur     , N.steps); % Pre-fit Error Covariance (difference between measurement and estimation)
    % % K_k_0      = zeros(N.states     , N.measur     , N.steps); % Kalman Gains
    % K_k_0      =  ones(N.states     , N.measur     , N.steps); % Kalman Gains
    % P_k_0      = zeros(N.states     , N.states     , N.steps); % Error Covariance Update Phase
    % P_k_pre_0  = zeros(N.states     , N.states     , N.steps); % Error Covariance Predict Phase
    % CPCT_0     = zeros(N.measur_all , N.measur_all , N.steps); % Error Covariance Predict Phase
    % CPCT_pre_0 = zeros(N.measur_all , N.measur_all , N.steps); % Error Covariance Predict Phase
    % CK_k_0     = zeros(N.measur     , N.measur     , N.steps); % Error Covariance Predict Phase

    S_k_0      = zeros(N.measur     , N.measur     , N.steps_final); % Pre-fit Error Covariance (difference between measurement and estimation)
    % K_k_0      = zeros(N.states     , N.measur     , N.steps_final); % Kalman Gains
    K_k_0      =  ones(N.states     , N.measur     , N.steps_final); % Kalman Gains
    P_k_0      = zeros(N.states     , N.states     , N.steps_final); % Error Covariance Update Phase
    P_k_pre_0  = zeros(N.states     , N.states     , N.steps_final); % Error Covariance Predict Phase
    CPCT_0     = zeros(N.measur_all , N.measur_all , N.steps_final); % Error Covariance Predict Phase
    CPCT_pre_0 = zeros(N.measur_all , N.measur_all , N.steps_final); % Error Covariance Predict Phase
    CK_k_0     = zeros(N.measur     , N.measur     , N.steps_final); % Error Covariance Predict Phase

    confidence        = 1.0; % Range from [0,1]
    P_k_0(:,:,1)      = confidence*eye(N.states);
    P_k_pre_0(:,:,1)  = confidence*eye(N.states);
    CPCT_0(:,:,1)     = C_DT * P_k_0(:,:,1) * C_DT'; 
    CPCT_pre_0(:,:,1) = C_DT * P_k_pre_0(:,:,1) * C_DT'; 
    CK_k_0(:,:,1)     = C_DT_CV * K_k_0(:,:,1);


%% Loop Through Noise Vectors
    for QQ = Q_vec
        for RR = R_vec
            disp(['Q= ' num2str(QQ) '  R= ' num2str(RR)])
            %% Check if Results Exist
            RUNSIM = false;
            
            save_filename = getCovarEstimationAnalysis(FLAG , QQ , RR);
            if isfile(save_filename)
                if FLAG.OverwriteResults
                    RUNSIM = true;
                    delete(save_filename)
                end
            else
                RUNSIM = true;
            end

            %% Run Simulation
            if RUNSIM
                Q_matrix = B_DT*QQ*B_DT';
                S_k      = S_k_0;
                K_k      = K_k_0;
                P_k      = P_k_0;
                P_k_pre  = P_k_pre_0;
                CPCT     = CPCT_0;
                CPCT_pre = CPCT_pre_0;
                CK_k     = CK_k_0;

                %% Run Estimation
                for i = 2:N.steps_final
                    % Predict Phase
                    P_k_pre(:,:,i)  = A_DT * P_k(:,:,i-1) * A_DT'  +  Q_matrix;
                    CPCT_pre(:,:,i) = C_DT * P_k_pre(:,:,i) * C_DT';

                    % Update (Correction) Phase
                    S_k(:,:,i)  = C_DT_CV * P_k_pre(:,:,i) * C_DT_CV' + RR;
                    K_k(:,:,i)  = P_k_pre(:,:,i) * C_DT_CV' * inv(S_k(:,:,i));
                    P_k(:,:,i)  = (eye(N.states) - K_k(:,:,i) * C_DT_CV) * P_k_pre(:,:,i);
                    CPCT(:,:,i) = C_DT * P_k(:,:,i) * C_DT';
                    CK_k(:,:,i) = C_DT_CV * K_k(:,:,i);
                end


                %% Get idx of actual samples
                S_k_red      = S_k(:,:,N.sample_vec);
                K_k_red      = K_k(:,:,N.sample_vec);
                P_k_red      = P_k(:,:,N.sample_vec);
                P_k_pre_red  = P_k_pre(:,:,N.sample_vec);
                CPCT_red     = CPCT(:,:,N.sample_vec);
                CPCT_pre_red = CPCT_pre(:,:,N.sample_vec);
                CK_k_red     = CK_k(:,:,N.sample_vec);


                %% Get diag of CPCT
                for i = 2:N.steps
                    CPCT_pre_diag(:,i) = diag( CPCT_pre_red(:,:,i) );
                    CPCT_diag(:,i)     = diag( CPCT_red(:,:,i) );
                end


                %% Reshape CK_k
                CK       = reshape(CK_k_red,1,[]);
                CK_infty = C_DT_CV * K_infty;

                CP_inftyCT = C_DT * P_infty * C_DT';


                %% Save Results
                if FLAG.SAVEResults
                    save(save_filename, 'QQ' , 'RR' , 'SIM' , 'FLAG' , 'N' , 'P' , 'sys' , 'S_k_red' , 'K_k_red' , 'P_k_red' , 'P_k_pre_red' , 'CPCT_red' , 'CPCT_pre_red' , 'CK_k_red' , 'CPCT_pre_diag' , 'CPCT_diag' , 'CK' , 'CK_infty' , 'CP_inftyCT' )
                end


                %% Plot Results
                if FLAG.PLOT.Any
                    if FLAG.PLOT.KALMAN_GAIN
                        % Kalman Gain (CK) of just voltage
                        figure
                        plot( N.sample_vec , CK , '-k' , 'LineWidth',2 )
                        title('Cell Voltage Kalman Gain')
                        xlabel('Samples')
                        ylabel('C K_k')
    
    
                        % Kalman Gain (CK) of just voltage with Asymptotic
                        figure
                        hold on
                        plot( N.sample_vec , CK , '-k' , 'LineWidth',2 )
                        yline(CK_infty, 'LineWidth',2,'Color','b')
                        title('Cell Voltage Kalman Gain')
                        xlabel('Samples')
                        ylabel('C K_k')
                    end
    
                    if FLAG.PLOT.COVAR
                        % % Plot Each CPCT with Asymptotic
                        %     for j = 1:N.measur_all
                        %         figure
                        %         hold on
                        %         plot( N.sample_vec , CPCT_diag(j,:) , '-' , 'LineWidth',2 )
                        %         yL = yline(CP_inftyCT(j,j)       , '--', 'LineWidth',2, 'Color','k');
                        %         yL.Alpha = 1;
                        %         title(RESULTS.Labels.title{j})
                        %     end
    
    
                        % Covariance (CPCT) of each variable with respect to steps
                        figure
                        hold on
                        for j = 1:N.measur_all
                            plot( N.sample_vec , CPCT_diag(j,:) , '-' , 'LineWidth',2, 'Color' , ColorVec(j,:) , 'DisplayName' , RESULTS.Labels.title{j} )
                            % yL = yline(CP_inftyCT(j,j)       , '--', 'LineWidth',2, 'Color' , ColorVec(j,:) , 'DisplayName' , [RESULTS.Labels.title{j} '_\infty']);
                            % yL.Alpha = 1;
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
    
    
                        % Covariance (CPCT) of each variable with respect to steps
                        figure
                        hold on
                        for j = 1:N.measur_all
                            plot( N.sample_vec , CPCT_diag(j,:) , '-' , 'LineWidth',2, 'Color' , ColorVec(j,:) , 'DisplayName' , RESULTS.Labels.title{j} )
                            % yL = yline(CP_inftyCT(j,j)       , '--', 'LineWidth',2, 'Color' , ColorVec(j,:) , 'DisplayName' , [RESULTS.Labels.title{j} '_\infty']);
                            % yL.Alpha = 1;
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
    
    
                        % Covariance (CPCT) of each variable with respect to steps
                        figure
                        hold on
                        for j = 1:N.measur_all
                            plot( N.sample_vec , CPCT_diag(j,:) , '-' , 'LineWidth',2, 'Color' , ColorVec(j,:) , 'DisplayName' , RESULTS.Labels.title{j} )
                            % yL = yline(CP_inftyCT(j,j)       , '--', 'LineWidth',2, 'Color' , ColorVec(j,:) , 'DisplayName' , [RESULTS.Labels.title{j} '_\infty']);
                            % yL.Alpha = 1;
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
    
                        % % Covariance (CPCT) of each variable with respect to steps
                        %     figure
                        %     hold on
                        %     for j = 1:N.measur_all
                        %         plot( N.sample_vec , CPCT_diag(j,:) , '-' , 'LineWidth',2 , 'DisplayName' , RESULTS.Labels.title{j} )
                        %     end
                        %     lgn = legend;
                        %     title('CP_kC^T')
                        %     xlabel('Samples')
                        %     ylabel('Covariance')
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




end