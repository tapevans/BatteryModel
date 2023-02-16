%% Get Ho-Kalman ROM
function [sys] = getHoKalmanROM(SIM,N,P,FLAG)
%% Get Impulse Response
    FLAG.StateMode = 5;
    FLAG.InputMode = 4;
    
    N_samples = 50;
    SIM.time_rest = SIM.Ts;
    SIM.time_ramp = SIM.Ts/10;
    SIM.t_final   = (N_samples + 1)*SIM.Ts;
    SIM.t_vec_imp = 0:SIM.Ts/10:SIM.t_final;
    SIM.t_vec = SIM.t_vec_imp;

    SIM.V_init = 0;
    SIM.V_imp = 1;
%     SIM.V_init = SIM.V_init;    
%     SIM.V_imp = SIM.V_step_height;
    
    %SIM.V_imp = 1/SIM.Ts;

    [InputSignal] = getInputSignal(SIM,N,P,FLAG);

%     figure
%     plot(InputSignal(:,1),InputSignal(:,2),'-o')

    %% Run Impulse Response
    %% Get IC
    SIM.V_init = 0;
    SIM.V_step = 1;
    % Initial Steady State
    R_tot = SIM.R_1 + SIM.R_2 + SIM.R_3 + SIM.R_4;
    i_PS  = SIM.V_init / R_tot;
    V_1   = SIM.V_init - SIM.R_1 * i_PS;
    V_2   = V_1        - SIM.R_2 * i_PS;
    V_3   = V_2        - SIM.R_3 * i_PS;
    delV_1= SIM.V_init - V_1;
    delV_2= V_1        - V_2;
    delV_3= V_2        - V_3;
    delV_4= V_3        - 0;

    old_inputMode = FLAG.InputMode;
    FLAG.InputMode = 1;
    switch FLAG.InputMode
        case 1 % Step
            % Initial Step
            V_eff = SIM.V_step - delV_1 -delV_4;
            R_eff = SIM.R_2 + SIM.R_3;
            i_eff = V_eff/R_eff;
            delV_2_eff = i_eff*SIM.R_2;

            % 3 State System
            SIM.x_0_3        = zeros(P.V_3,1);
            SIM.x_0_3(P.V_1) = SIM.V_step - delV_1;
            SIM.x_0_3(P.V_3) = delV_4;
            SIM.x_0_3(P.V_2) = (SIM.V_step - delV_1) - delV_2_eff ;

            % 5 State System
            SIM.x_0_5 = zeros(5,1);

            SIM.x_0_5(P.V_s)  = SIM.V_init;
            SIM.x_0_5(P.i_PS) = i_PS;
            SIM.x_0_5(P.V_1)  = V_1;
            SIM.x_0_5(P.V_2)  = V_2;
            SIM.x_0_5(P.V_3)  = V_3;
        case 2 %%%%%%%%%%%%%%%%%%%%%%%%%%Add more for inital DC and DC offset
            % 3 State System
            SIM.x_0_3 = zeros(P.V_3,1);
            % 5 State System
            SIM.x_0_5 = zeros(5,1);
    end
    FLAG.InputMode = old_inputMode;


    [t_imp, x_imp, z_imp] = runODE(SIM,N,P,FLAG);

%     figure
%     hold on
%     plot(t_imp,x_imp(:,P.V_1),'-ok')
%     plot(t_imp,x_imp(:,P.V_2),'-ob')
%     plot(t_imp,x_imp(:,P.V_3),'-og')
%     title('ODE impulse Response')

%% Sort Data for g_k
    idx_vec = 11:10:length(SIM.t_vec);

    g_k = nan(P.V_3,N_samples+1);

    g_k_temp = z_imp(idx_vec,1:3)';
    g_k = g_k_temp - g_k_temp(:,1);


%     figure
%     hold on
%     plot(t_imp,x_imp(:,P.V_1),'-ok')
% %     plot(t_imp(idx_vec,:),g_k(P.V_1,:)+g_k_temp(P.V_1,1),'or','MarkerFaceColor','r')
%     plot(t_imp(idx_vec,:),g_k(P.V_1,:)+SIM.x_0_5(P.V_1,1),'or','MarkerFaceColor','r')
% %     ylim([0,SIM.V_imp+SIM.V_init])
%     title('Data Used in Ho-Kalman Reduction V_1')
% 
%     figure
%     hold on
%     plot(t_imp,x_imp(:,P.V_2),'-ob')
% %     plot(t_imp(idx_vec,:),g_k(P.V_2,:)+g_k_temp(P.V_2,1),'or','MarkerFaceColor','r')
%     plot(t_imp(idx_vec,:),g_k(P.V_2,:)+SIM.x_0_5(P.V_2,1),'or','MarkerFaceColor','r')
% %     ylim([0,SIM.V_imp+SIM.V_init])
%     title('Data Used in Ho-Kalman Reduction V_2')
% 
%     figure
%     hold on
%     plot(t_imp,x_imp(:,P.V_3),'-og')
% %     plot(t_imp(idx_vec,:),g_k(P.V_3,:)+g_k_temp(P.V_3,1),'or','MarkerFaceColor','r')
%     plot(t_imp(idx_vec,:),g_k(P.V_3,:)+SIM.x_0_5(P.V_3,1),'or','MarkerFaceColor','r')
% %     ylim([0,SIM.V_imp+SIM.V_init])
%     title('Data Used in Ho-Kalman Reduction V_3')

%     figure
%     hold on
%     plot(t_imp,x_imp(:,P.V_s),'-ob')
%     plot(t_imp(idx_vec,:),x_imp(idx_vec,P.V_s),'or','MarkerFaceColor','r')
% 
%     figure
%     hold on
%     plot(t_imp,x_imp(:,P.i_PS),'-ob')
%     plot(t_imp(idx_vec,:),x_imp(idx_vec,P.i_PS),'or','MarkerFaceColor','r')



%% Get Hankle Matrix
    H = myHankle(g_k);

%% Convert to SS DT sys
% SVD of Hankle
    [Nrows, Ncolms] = size(H);
    N_outputs = Nrows/Ncolms;
    N_in = 1; %%% !!! Hardcoded but always true for batteries
    
    [U,S,V] = svd(H);
    r = rank(S,1e-8);
%     r = 5;
    
    U_colm = U(:,1:r);
    S_colm = S(1:r,1:r);
    V_colm = V(:,1:r);
    % U_null = U(:,r+1:end);
    % S_null = S(r+1:end,r+1:end);
    % V_null = V(:,r+1:end);
    
%     if localFLAG.BalRed % If using balanced reduction
        S_sqrt = S_colm.^0.5;
        obsv_r = U_colm*S_sqrt;
        cont_r = S_sqrt*V_colm';
%     else % If using normal reduction
%         obsv_r = U_colm;
%         cont_r = S_colm*V_colm';
%     end

% C_r
    C_r = obsv_r(1:N_outputs,:   );
    
% B_r
    B_r = cont_r(:        ,N_in);
    
% A_r
    P   = obsv_r(1:end-N_outputs,:);
    P_p = obsv_r(N_outputs+1:end,:);
    threshold = 1e-7;
    P_inv = pinv(P,threshold);
    A_r = P_inv*P_p;

% D_r ~ !!!!! I'm assuming this is correct
    D_r = zeros(N_outputs,N_in);


%% Return the system
sys = ss(A_r, B_r, C_r, D_r, SIM.Ts);

end