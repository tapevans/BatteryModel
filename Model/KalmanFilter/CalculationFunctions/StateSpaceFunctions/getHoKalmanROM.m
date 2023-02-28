%% Get Ho-Kalman ROM
function [sys] = getHoKalmanROM(SIM,N,P,FLAG,RESULTS)
%% Load Impulse Data
    filename = getImpulseFilename(FLAG);
    simsys = load(filename);


%% Sort Data for g_k
    % Extract g_k
    idx_vec = 1:simsys.SIM.TsMultiple:length(simsys.t_soln);

    i_user_DT = simsys.i_user(idx_vec);
    pulse_idx = find(i_user_DT ~= 0 );
    
    z_imp = SIM.OutputMatrix * simsys.SV_soln';
    g_k_temp = z_imp(:,idx_vec(pulse_idx-1:end));
    IC  = z_imp(:,1);
    g_k = g_k_temp - g_k_temp(:,1);

    t_imp = simsys.t_soln;


%% Plot Impulse Response
    if FLAG.Analysis.PlotImp
        for i = 1:N.DesOut
        figure
        hold on
        plot(t_imp           ,z_imp(i,:)               ,'-ok')
        plot(t_imp(idx_vec,:),g_k(i,:)+IC(i,1),'or','MarkerFaceColor','r')
        title(['Data Used in Ho-Kalman Reduction ' RESULTS.Labels.title{i}])
        end
    end


%% Get Hankle Matrix
    H = myHankle(g_k);


%% Convert to SS DT sys
% SVD of Hankle
    [Nrows, Ncolms] = size(H);
    N_outputs = Nrows/Ncolms;
    N_in = 1; %%% !!! Hardcoded but always true for batteries
    
    [U,S,V] = svd(H);
    r = rank(S,1.15e-7);
%     r = 18;
    
    U_colm = U(:,1:r);
    S_colm = S(1:r,1:r);
    V_colm = V(:,1:r);

    % Balanced Reduction
    S_sqrt = S_colm.^0.5;
    obsv_r = U_colm*S_sqrt;
    cont_r = S_sqrt*V_colm';

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


%%%%%%%%%%%%%%%%%%%%%%%%%% OOOOOOOOOOOLD %%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Get Impulse Response
%     FLAG.InputMode = 4;
%     
%     N_samples = 50;
%     SIM.time_rest = SIM.Ts;
%     SIM.time_ramp = SIM.Ts/10;
%     SIM.t_final   = (N_samples + 1)*SIM.Ts;
%     SIM.t_vec_imp = 0:SIM.Ts/10:SIM.t_final;
%     SIM.t_vec = SIM.t_vec_imp;
% 
%     SIM.F_init = 0;
%     SIM.F_imp = 1;
% %     SIM.F_init = SIM.F_init;    
% %     SIM.F_imp = SIM.F_step_height;
%     
%     %SIM.F_imp = 1/SIM.Ts;
% 
%     [InputSignal] = getInputSignal(SIM,N,P,FLAG);
% 
%     if FLAG.Analysis.PlotImp
%         figure
%         plot(InputSignal(:,1),InputSignal(:,2),'-o')
%         title('Impulse Signal')
%     end

    %% Run Impulse Response
%     SIM.x_0 = [0 ,0 ,0 ,0 ]';
%     
%     [t_imp, x_imp, z_imp] = runODE(SIM,N,P,FLAG);
% 
%     if FLAG.Analysis.PlotImp
%         figure
%         hold on
%         plot(t_imp,x_imp(:,P.x1),'-ok','DisplayName','x_1')
%         plot(t_imp,x_imp(:,P.v1),'-ob','DisplayName','v_1')
%         plot(t_imp,x_imp(:,P.x2),'-og','DisplayName','x_2')
%         plot(t_imp,x_imp(:,P.v2),'-oc','DisplayName','v_2')
%         title('ODE impulse Response')
%         lgn = legend;
%         lgn.Location = 'northeast';
%     end

% %     plot(t_imp(idx_vec,:),g_k(P.x1,:)+g_k_temp(P.x1,1),'or','MarkerFaceColor','r')
% %     plot(t_imp(idx_vec,:),g_k(P.v1,:)+g_k_temp(P.v1,1),'or','MarkerFaceColor','r')
% %     plot(t_imp(idx_vec,:),g_k(P.x2,:)+g_k_temp(P.x2,1),'or','MarkerFaceColor','r')
% %     plot(t_imp(idx_vec,:),g_k(P.v2,:)+g_k_temp(P.v2,1),'or','MarkerFaceColor','r')

%%
%         figure
%         hold on
%         plot(t_imp,x_imp(:,P.v1),'-ob')
%         plot(t_imp(idx_vec,:),g_k(P.v1,:)+SIM.x_0(P.v1,1),'or','MarkerFaceColor','r')
%         title('Data Used in Ho-Kalman Reduction v_1')
%     
%         figure
%         hold on
%         plot(t_imp,x_imp(:,P.x2),'-og')
%         plot(t_imp(idx_vec,:),g_k(P.x2,:)+SIM.x_0(P.x2,1),'or','MarkerFaceColor','r')
%         title('Data Used in Ho-Kalman Reduction x_2')
%     
%         figure
%         hold on
%         plot(t_imp,x_imp(:,P.v2),'-oc')
%         plot(t_imp(idx_vec,:),g_k(P.v2,:)+SIM.x_0(P.v2,1),'or','MarkerFaceColor','r')
%         title('Data Used in Ho-Kalman Reduction v_2')


%%
%     g_k = zeros(numP , length(idx_vec)-(pulse_idx-2)); %%%%%%%%% -2 is hardcoded



%     idx_vec = 11:10:length(SIM.t_vec);
% 
%     g_k = nan(P.v2,N_samples+1);
% 
%     g_k_temp = z_imp(idx_vec,:)';
%     g_k = g_k_temp - g_k_temp(:,1);
