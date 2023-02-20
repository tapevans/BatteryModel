%% Create HoKalman SS from Impulse

clear all; close all; clc;

%%
% FLAGS
% States to estimate
localFLAG.cell_voltage  = 1; % Terminal Voltage
localFLAG.delta_phi     = 1; % Electrostatic potential difference between active material and electrolyte @AN/SEP interface
localFLAG.C_Liion       = 1; % Concentration of Li^+ in the electrolyte @AN/SEP interface
localFLAG.C_Li          = 1; % Concentration of Li at the surface of active material @AN/SEP interface
localFLAG.delta_C_Li    = 1; % Difference of Li concentration between the surface and center of active material particle @AN/SEP interface
localFLAG.i_Far         = 1; % Chemical reaction current density at SEI @AN/SEP interface
localFLAG.T             = 1; % Temperature of control volume @AN/SEP interface

%% Set pointers for variables
i = 1;
if localFLAG.cell_voltage   % Cell Voltage
    %     PLOT_STR.title{i} = 'Cell Voltage';
    %
    %     idx = find(PLOT_STR.title{i} ==' ');
    %     KalmanStr{i} = PLOT_STR.title{i};
    %     KalmanStr{i}(idx) = '';

    localP.cell_voltage  = i; i = i + 1;
end
if localFLAG.delta_phi      % Delta Phi
    %     PLOT_STR.title{i} = '\Delta \phi';
    %
    %     idx = find(PLOT_STR.title{i} ==' ');
    %     KalmanStr{i} = PLOT_STR.title{i};
    %     KalmanStr{i}(idx) = '';
    %
    %
    localP.delta_phi     = i; i = i + 1;
end
if localFLAG.C_Liion        % C_Li+
    %     PLOT_STR.title{i} = 'C_{Li^+}';
    %
    %     idx = find(PLOT_STR.title{i} ==' ');
    %     KalmanStr{i} = PLOT_STR.title{i};
    %     KalmanStr{i}(idx) = '';
    %
    %
    localP.C_Liion       = i; i = i + 1;
end
if localFLAG.C_Li           % C_Li
    %     PLOT_STR.title{i} = 'C_{Li,surf}';
    %
    %     idx = find(PLOT_STR.title{i} ==' ');
    %     KalmanStr{i} = PLOT_STR.title{i};
    %     KalmanStr{i}(idx) = '';
    %
    %
    localP.C_Li          = i; i = i + 1;
end
if localFLAG.delta_C_Li     % Delta C_Li
    %     PLOT_STR.title{i} = '\Delta C_{Li}';
    %
    %     idx = find(PLOT_STR.title{i} ==' ');
    %     KalmanStr{i} = PLOT_STR.title{i};
    %     KalmanStr{i}(idx) = '';
    %
    %
    localP.delta_C_Li    = i; i = i + 1;
end
if localFLAG.i_Far          % i_Far
    %     PLOT_STR.title{i} = 'i_{Far}';
    %
    %     idx = find(PLOT_STR.title{i} ==' ');
    %     KalmanStr{i} = PLOT_STR.title{i};
    %     KalmanStr{i}(idx) = '';
    %
    %
    localP.i_Far         = i; i = i + 1;
end
if localFLAG.T              % Temperature
    %     PLOT_STR.title{i} = 'Temperature';
    %
    %     idx = find(PLOT_STR.title{i} ==' ');
    %     KalmanStr{i} = PLOT_STR.title{i};
    %     KalmanStr{i}(idx) = '';
    %
    %
    localP.T             = i; i = i + 1;
end

numP = i-1;
%%
% Desired Sampling Times
N_t_s   = 25;   % **Keep this fix for now
T_s_min = -1; % T_s = 10^(T_s_min), **Keep this fix for now
T_s_max =  1; % T_s = 10^(T_s_max), **Keep this fix for now
Ts_vec = logspace(T_s_min,T_s_max,N_t_s);
% Ts_vec  = 1;

% Split up like this for Github pushing
% SOC_vec = 0:1:25;
% SOC_vec = 26:1:49;
% SOC_vec = 51:1:75;
% SOC_vec = 76:1:100;
% SOC_vec = 50;

SOC_vec = 0:1:100;


Ts_vec_str = {'0.1',
    '0.121152765862859',
    '0.146779926762207',
    '0.177827941003892',
    '0.215443469003188',
    '0.261015721568254',
    '0.316227766016838',
    '0.383118684955729',
    '0.464158883361278',
    '0.562341325190349',
    '0.681292069057961',
    '0.825404185268018',
    '1.0',
    '1.21152765862859',
    '1.46779926762207',
    '1.77827941003892',
    '2.15443469003188',
    '2.61015721568254',
    '3.16227766016838',
    '3.83118684955729',
    '4.64158883361278',
    '5.62341325190349',
    '6.81292069057961',
    '8.25404185268018',
    '10.0'};

% Ts_vec_str = {'1.0'};

for SS = 1:length(SOC_vec)
    for TT = 1:length(Ts_vec)
        %% Get Filename
        SOC = SOC_vec(SS);
        Ts = Ts_vec(TT);
        filepath = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\ObservabilityTest';
        filename = ['ObservabilityTest_KPCont_DTImpulseTs' num2str(Ts_vec_str{TT}) 'SOC' num2str(SOC) '.mat'];
%         load([filepath filesep filename])
        %% Load Impulse Data
        % Load Simulation
        simsys = load([filepath filesep filename]);
%         simsys = load([filepath filesep filename],'t_soln','i_user','SIM','N','cell_voltage','del_phi','C_Liion','C_Li','C_Li_surf','i_Far','TemperatureC');

        % Extract g_k
        idx_vec = 1:simsys.SIM.TsMultiple:length(simsys.t_soln);

        i_user_DT = simsys.i_user(idx_vec);
        pulse_idx = find(i_user_DT ~= 0 );

        g_k = zeros(numP , length(idx_vec)-(pulse_idx-2)); %%%%%%%%% -2 is hardcoded
        IC  = zeros(numP,1);

        if localFLAG.cell_voltage   % Cell Voltage
            cell_voltage_DT = simsys.cell_voltage(idx_vec);
            g_k_temp = (cell_voltage_DT(pulse_idx-1:end))';
            IC(localP.cell_voltage) = g_k_temp(1);
            g_k(localP.cell_voltage,:) = g_k_temp - g_k_temp(1);
        end
        if localFLAG.delta_phi      % Delta Phi
            del_phi_DT = simsys.del_phi(idx_vec,simsys.N.N_CV_AN);
            g_k_temp = (del_phi_DT(pulse_idx-1:end))';
            IC(localP.delta_phi) = g_k_temp(1);
            g_k(localP.delta_phi,:) = g_k_temp - g_k_temp(1);
        end
        if localFLAG.C_Liion        % C_Li+
            C_Liion_DT = simsys.C_Liion(idx_vec,simsys.N.N_CV_AN);
            g_k_temp = (C_Liion_DT(pulse_idx-1:end))';
            IC(localP.C_Liion) = g_k_temp(1);
            g_k(localP.C_Liion,:) = g_k_temp - g_k_temp(1);
        end
        if localFLAG.C_Li           % C_Li_surf
            C_Li_surf_DT = simsys.C_Li_surf(idx_vec,simsys.N.N_CV_AN);
            g_k_temp = (C_Li_surf_DT(pulse_idx-1:end))';
            IC(localP.C_Li) = g_k_temp(1);
            g_k(localP.C_Li,:) = g_k_temp - g_k_temp(1);
        end
        if localFLAG.delta_C_Li     % Delta C_Li
            C_Li_DT = simsys.C_Li(:,simsys.N.N_CV_AN,idx_vec);
            Del_C_Li = C_Li_DT(simsys.N.N_R_AN , :,:) - C_Li_DT(1 , :,:);
            Del_C_Li = reshape(Del_C_Li,1,[]);
            g_k_temp = (Del_C_Li(pulse_idx-1:end))';
            IC(localP.delta_C_Li) = g_k_temp(1);
            g_k(localP.delta_C_Li,:) = g_k_temp - g_k_temp(1);
        end
        if localFLAG.i_Far          % i_Far
            i_Far_DT = simsys.i_Far(idx_vec);
            g_k_temp = (i_Far_DT(pulse_idx-1:end))';
            IC(localP.i_Far) = g_k_temp(1);
            g_k(localP.i_Far,:) = g_k_temp - g_k_temp(1);
        end
        if localFLAG.T              % Temperature
            TemperatureC_DT = simsys.TemperatureC(idx_vec);
            g_k_temp = (TemperatureC_DT(pulse_idx-1:end))';
            IC(localP.T) = g_k_temp(1);
            g_k(localP.T,:) = g_k_temp - g_k_temp(1);
        end


        %% Plot Impulse Data
        if false
            figure
            hold on
            plot(t_soln         ,simsys.cell_voltage         ,'ok')
            %                 plot(t_soln(idx_vec),simsys.cell_voltage(idx_vec),'or','MarkerFaceColor','r')
            %                 plot(t_soln(idx_vec(3:end)),g_k(localP.cell_voltage,:),'or','MarkerFaceColor','r')
            plot(t_soln(idx_vec(3:end)),g_k(localP.cell_voltage,:)+simsys.cell_voltage(1),'or','MarkerFaceColor','r')
            title('Impulse Data: Cell Voltage')

            figure
            hold on
            plot(t_soln         ,simsys.del_phi(:,simsys.N.N_CV_AN)      ,'ok')
            plot(t_soln(idx_vec),simsys.del_phi(idx_vec,simsys.N.N_CV_AN),'or','MarkerFaceColor','r')
            title('Impulse Data: \Delta \phi')

            figure
            hold on
            plot(t_soln         ,simsys.C_Liion(:,simsys.N.N_CV_AN)      ,'ok')
            plot(t_soln(idx_vec),simsys.C_Liion(idx_vec,simsys.N.N_CV_AN),'or','MarkerFaceColor','r')
            title('Impulse Data: C_{Li^+}')

            figure
            hold on
            plot(t_soln         ,simsys.C_Li_surf(:,simsys.N.N_CV_AN)      ,'ok')
            plot(t_soln(idx_vec),simsys.C_Li_surf(idx_vec,simsys.N.N_CV_AN),'or','MarkerFaceColor','r')
            title('Impulse Data: C_{Li,surf}')

            %                 figure
            %                 hold on
            %                 plot(t_soln         ,C_Li(:,simsys.N.N_CV_AN,:)         ,'ok')
            %                 plot(t_soln(idx_vec),C_Li(:,simsys.N.N_CV_AN,idx_vec),'or','MarkerFaceColor','r')
            %                 title('Impulse Data: \Delta C_{Li}')

            figure
            hold on
            plot(t_soln         ,simsys.i_Far(:,simsys.N.N_CV_AN)      ,'-ok')
            plot(t_soln(idx_vec),simsys.i_Far(idx_vec,simsys.N.N_CV_AN),'or','MarkerFaceColor','r')
            title('Impulse Data: i_{Far}')

            %                 figure
            %                 hold on
            %                 plot(t_soln         ,simsys.TemperatureC         ,'ok')
            %                 plot(t_soln(idx_vec),simsys.TemperatureC(idx_vec),'or','MarkerFaceColor','r')
            %                 title('Impulse Data: Temperature')

        end
        %% Get Hankle Matrix
        H = myHankle(g_k);

        %% Convert to SS DT sys
        % SVD of Hankle
        [Nrows, Ncolms] = size(H);
        N_outputs = Nrows/Ncolms;
        N_in = 1; %%% !!! Hardcoded but always true for batteries

        [U,S,V] = svd(H);
        r = rank(S,1e-7);
        %     r = 5;

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
        sys = ss(A_r, B_r, C_r, D_r, Ts);

        %% Save System
        save_filepath = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanFilter\DATA\SS_HoKalman';
        save_filename = ['SS_HK_SOC' num2str(SOC) '_Ts' num2str(Ts) '.mat'];
        save([save_filepath filesep save_filename],'sys','IC')

    end
end