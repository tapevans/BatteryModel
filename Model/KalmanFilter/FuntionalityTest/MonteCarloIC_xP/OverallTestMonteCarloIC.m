%% Overall Test for Monte Carlo Initialization
    clear all; close all; clc;


%% FLAGS
    FLAG.UseSingleMeasurement = 1;
        FLAG.MeasurementToUse = 3;
    
% Plots
    FLAG.PLOT.ANY         = 1;
    FLAG.PLOT.InputSignal = 0;
    FLAG.PLOT.NoisyPlant  = 0;
    FLAG.PLOT.AsyEst      = 1;
    FLAG.PLOT.VarEst      = 1;


%% Parameters
    % omega       = 2*pi/(3600*24); % GeoSynchronous
    R_0         = 1e-1; % Measurement Noise
    Q_0         = 1e0; % Process Noise
    Ts          = 1e-4;    % Sampling Rate
    Tswitch     = 10;   % Switching Time on PRBS
    PRBS_Amp    = 1;    % Amplitude of the PRBS signal
    IC_multiple = 1.2;    % What to multiply P_infty by for the initial P_k|k-1


%% ROM
%%%%%%%%%%%%%% Spacecraft, unstable in DT
    % A = [0          1         0  0
    %      3*omega^2  0         0  2*omega
    %      0          0         0  1
    %      0          -2*omega  0  0];
    % B = [0 0 
    %      1 0
    %      0 0 
    %      0 1];
    % 
    % C = [1 0 0 0
    %      0 0 1 0];
    % 
    % [N_outputs , N_states] = size(C);
    % [~         , N_inputs] = size(B);
    % 
    % D = zeros(N_outputs,N_inputs);

%%%%%%%%%%%%%% Problem 7.34, Unstable
    % A = [-2  1
    %       1  0];
    % B = [1
    %      0];
    % 
    % C = [1 2
    %      1 0
    %      0 1];
    % 
    % [N_outputs , N_states] = size(C);
    % [~         , N_inputs] = size(B);
    % 
    % D = zeros(N_outputs,N_inputs);

%%%%%%%%%%%%%% Problem 7.49, Unstable
    A = [-0.4    0  -0.01
            1    0      0
         -1.4  9.8  -0.02];
    B = [6.3
           0
         9.8];

    % C = [0 0 1];
    C = eye(3);

    [N_outputs , N_states] = size(C);
    [~         , N_inputs] = size(B);

    D = zeros(N_outputs,N_inputs);


%% Get Input Vector
    time = 0:Ts:2;
    N_time_steps = length(time);
    u_k = zeros(N_inputs , N_time_steps);

% Plot Time 
    if FLAG.PLOT.ANY && FLAG.PLOT.InputSignal
        figure
        plot(time,u_k,'LineWidth',2)
        xlabel('Time')
        ylabel('Input Thrust')
        title('Input Signal')
    end


%% Get DT ROM
    sys   = ss(A,B,C,D);
    sysDT = c2d(sys,Ts);

    % eigVal = eig(sysDT.A)


%% Initialize Q and R
% Process Noise
    Q = Q_0 * eye(N_states);
    
% Measurement Noise
    if FLAG.UseSingleMeasurement
        N_meas_outputs = 1;
        meas_vec = [FLAG.MeasurementToUse];
    else
        N_meas_outputs = N_outputs;
        meas_vec = (1:N_outputs)';
    end
    R = R_0 * eye(N_meas_outputs);
    
% Split the measured matrix
    C_meas = sysDT.C(meas_vec,:);


%% Get P_infty
    % [P_infty  ,~,~] =  dare(sysDT.A', sysDT.C', Q, R);
    % K_infty = P_infty*sysDT.C'*inv(sysDT.C*P_infty*sysDT.C' + R);
    [P_infty ,~,~] =  dare(sysDT.A', C_meas', Q, R);
    K_infty         = P_infty * C_meas' * inv(C_meas * P_infty * C_meas' + R);


%% Simulate Plant
% Noise vectors
    w_k = mvnrnd(zeros(1,N_states)       , Q , N_time_steps);   w_k = w_k';
    v_k = mvnrnd(zeros(1,N_meas_outputs) , R , N_time_steps);   v_k = v_k';

    % covar_p = data2Covar(w_k)
    % covar_m = data2Covar(v_k)
    % mean(w_k,2)
    % mean(v_k,2)

% Initialize Simulation Variables
    x_plant = nan(N_states , N_time_steps);
    x_plant(:,1) = zeros(N_states , 1);

% Simulate Plant
    for i = 1:N_time_steps-1
        x_plant(:,i+1) = sysDT.A * x_plant(:,i) + sysDT.B * u_k(:,i) + w_k(:,i);
    end

    y_plant      = sysDT.C * x_plant;
    y_plant_meas = y_plant(meas_vec,:) + v_k;

% Plot Outputs
    if FLAG.PLOT.ANY && FLAG.PLOT.NoisyPlant
        figure
        hold on
        for i = 1:N_meas_outputs
            plot(time,y_plant_meas(i,:),'LineWidth',2)
        end
        xlabel('Time (s)')
        ylabel('Output')
        title('Noisy Plant Results')
    end


%% Initialize Estimation
% Split the measured matrix
    C_meas = sysDT.C(meas_vec,:);

% Get inflated state error covariance
    P_0 = P_infty*IC_multiple;

% Get states initial distribution
    % rng('default')  % For reproducibility
    n = 1000;    
    x_distribution = mvnrnd( zeros(1,N_states) , P_0 , n ); 
    x_distribution = x_distribution'; 

% Choose a random state vector as x_IC
    idx_IC = randi(n,1,1);
    x_IC   = x_distribution(:,idx_IC);


%% Run Asymptotic Estimation
   x_asy     = zeros(N_states      ,N_time_steps);   % Estimator States
   x_asy_pre = zeros(N_states      ,N_time_steps);   % Estimator States predict phase
   y_tilde_k = zeros(N_meas_outputs,N_time_steps);   % Estimation Error 

   x_asy(:,1)     = x_IC;
   x_asy_pre(:,1) = x_IC;

   for i = 2:N_time_steps
        % Predict Phase
            x_asy_pre(:,i) = sysDT.A * x_asy(:,i-1) + sysDT.B * u_k(:,i-1);
    
        % Update (Correction) Phase
            y_tilde_k(:,i) = y_plant_meas(:,i) - ( C_meas  * x_asy_pre(:,i) );
            x_asy(:,i)     = x_asy_pre(:,i) + K_infty * y_tilde_k(:,i);
   end
   y_est_asy_pre = sysDT.C * x_asy_pre;
   y_est_asy     = sysDT.C * x_asy;


%% Run Variable Estimation
% Initialize
    x_var     = zeros(N_states                       ,N_time_steps); % Estimator States
    x_var_pre = zeros(N_states                       ,N_time_steps); % Estimator States predict phase
    y_tilde_k = zeros(N_meas_outputs                 ,N_time_steps); % Estimation Error
    S_k       = zeros(N_meas_outputs ,N_meas_outputs ,N_time_steps); % Error Covariance
    K_k       = zeros(N_states       ,N_meas_outputs ,N_time_steps); % Kalman Gains
    P_k       = zeros(N_states       ,N_states       ,N_time_steps); % Error Covariance
    P_k_pre   = zeros(N_states       ,N_states       ,N_time_steps); % Error Covariance predict phase
    
    P_k(:,:,1)     = P_0;
    P_k_pre(:,:,1) = P_0; 
    
    x_var(:,1)     = x_IC;
    x_var_pre(:,1) = x_IC;
    
% Run Estimation
    for i = 2:N_time_steps
    % Predict Phase
        x_var_pre(:,i) = sysDT.A * x_var(:,i-1)  +  sysDT.B * u_k(:,i-1);
        P_k_pre(:,:,i) = sysDT.A * P_k(:,:,i-1) * sysDT.A'  +  Q;
    
    % Update (Correction) Phase
        y_tilde_k(:,i) = y_plant_meas(:,i) - C_meas * x_var_pre(:,i);
        S_k(:,:,i)     = C_meas * P_k_pre(:,:,i) * C_meas' + R;
        K_k(:,:,i)     = P_k_pre(:,:,i) * C_meas' * inv(S_k(:,:,i));
        x_var(:,i)     = x_var_pre(:,i) + K_k(:,:,i) * y_tilde_k(:,i);
        P_k(:,:,i)     = (eye(N_states) - K_k(:,:,i) * C_meas) * P_k_pre(:,:,i);
    end
    y_est_var_pre = sysDT.C * x_var_pre;
    y_est_var     = sysDT.C * x_var;

    disp('Final P_k')
    % disp(num2str(P_k(:,:,end)))
    disp(num2str(P_k_pre(:,:,end)))
    disp(newline)
    disp('P_infty')
    disp(num2str(P_infty))

    disp(newline)
    disp('Final Post Meas P_k')
    disp(num2str(P_k(:,:,end)))
    disp(newline)
    disp('Post P_infty')
    P_infty_post = (eye(N_states) - K_infty * C_meas) * P_infty;
    disp(num2str(P_infty_post))


%% Perform Error Calculations
    error_abs_states_asy = abs(x_asy - x_plant);
    error_abs_states_var = abs(x_var - x_plant);
    
    error_abs_outputs_asy = abs(y_est_asy - y_plant);
    error_abs_outputs_var = abs(y_est_var - y_plant);

    error_states_asy = x_asy - x_plant;
    error_states_var = x_var - x_plant;
    
    error_outputs_asy = y_est_asy - y_plant;
    error_outputs_var = y_est_var - y_plant;

    start_idx = floor(N_time_steps*7/8);

    % covar_es_asy = data2Covar(error_states_asy(:,start_idx:end))
    % covar_es_var = data2Covar(error_states_var(:,start_idx:end))  
    % 
    % covar_eo_asy = data2Covar(error_outputs_asy(:,start_idx:end))
    % covar_eo_var = data2Covar(error_outputs_var(:,start_idx:end))

    for i = 1:N_time_steps
        for j = 1:N_states
            std_var(j,i) = P_k(j,j,i)^0.5;
        end
    end

    xconf      = [time               flip(time)];
    % yconf_1sig = [x_plant+  std_var  flip(x_plant-  std_var)];
    % yconf_2sig = [x_plant+2*std_var  flip(x_plant-2*std_var)];
    % yconf_3sig = [x_plant+3*std_var  flip(x_plant-3*std_var)];
    yconf_1sig = [x_plant+  std_var  zeros(3,length(x_plant))];
    yconf_2sig = [x_plant+2*std_var  zeros(3,length(x_plant))];
    yconf_3sig = [x_plant+3*std_var  zeros(3,length(x_plant))];


%% Plot Performance
if FLAG.PLOT.ANY
    % % Asymptotic Estimator performance
    % for i = 1:N_outputs
    %     figure
    %     hold on
    %     plot(time,y_est_asy(i,:),'ro','LineWidth',2,'DisplayName','Asymptotic')
    %     plot(time,y_plant(i,:)  ,'k-','LineWidth',2,'DisplayName','Plant')
    %     xlabel('Time')
    %     ylabel('Output')
    %     title('Compare Asymptotic Estimator to Plant')
    %     lgn = legend;
    % end
    % 
    % % Variable Estimator performance
    % for i = 1:N_outputs
    %     figure
    %     hold on
    %     plot(time,y_est_var(i,:),'go','LineWidth',2,'DisplayName','Variable')
    %     plot(time,y_plant(i,:)  ,'k-','LineWidth',2,'DisplayName','Plant')
    %     xlabel('Time')
    %     ylabel('Output')
    %     title('Compare Variable Estimator to Plant')
    %     lgn = legend;
    % end

    % % Both Estimator performance
    % for i = 1:N_outputs
    %     figure
    %     hold on
    %     plot(time,y_est_asy(i,:),'ro','LineWidth',2,'DisplayName','Asymptotic')
    %     plot(time,y_est_var(i,:),'go','LineWidth',2,'DisplayName','Variable')
    %     plot(time,y_plant(i,:)  ,'k-','LineWidth',2,'DisplayName','Plant')
    %     xlabel('Time')
    %     ylabel('Output')
    %     title('Compare Variable Estimator to Plant')
    %     lgn = legend;
    % end
    % 
    % % Absolute Value of the state error for Asymptotic
    % figure
    % hold on
    % plot(time(2:end),error_abs_states_asy(1,2:end),'r','LineWidth',2,'DisplayName','State 1')
    % plot(time(2:end),error_abs_states_asy(2,2:end),'b','LineWidth',2,'DisplayName','State 2')
    % plot(time(2:end),error_abs_states_asy(3,2:end),'g','LineWidth',2,'DisplayName','State 3')
    % xlabel('Time')
    % title('Asymptotic states error')
    % lgn = legend;
    % 
    % % Absolute Value of the state error for Variable
    % figure
    % hold on
    % plot(time(2:end),error_abs_states_var(1,2:end),'r','LineWidth',2,'DisplayName','State 1')
    % plot(time(2:end),error_abs_states_var(2,2:end),'b','LineWidth',2,'DisplayName','State 2')
    % plot(time(2:end),error_abs_states_var(3,2:end),'g','LineWidth',2,'DisplayName','State 3')
    % xlabel('Time')
    % title('Variable states error')
    % lgn = legend;

    % State Confidence Intervals for Variable Estimator
    for i = 1:N_states
        figure
        hold on
        p3 = fill(xconf , yconf_3sig(i,:),'r');
        p3.FaceColor = [1 0.8 0.8];      
        p3.EdgeColor = 'none'; 
        % p2 = fill(xconf , yconf_2sig(i,:),'b');
        % p2.FaceColor = [0.8 0.8 1];      
        % p2.EdgeColor = 'none'; 
        % p1 = fill(xconf , yconf_1sig(i,:),'k');
        % p1.FaceColor = [0.5 0.5 0.5];      
        % p1.EdgeColor = 'none'; 
        plot(time,error_abs_states_var(i,:),'r')
        xlabel('Time')
        ylabel('State Error')
        title(['State ' num2str(i) ' Absolute Error with 3\sigma CI'])
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
