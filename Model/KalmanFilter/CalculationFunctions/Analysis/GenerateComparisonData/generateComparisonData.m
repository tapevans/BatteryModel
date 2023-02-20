%% Generate Comparison Data
function generateComparisonData(SIM,FLAG,N,P,RESULTS) % ,ESTIMATOR,ERROR_CALC
%% Check if the results exist
    RUNSIM = false;
    save_filename = getGenDataFilename(FLAG);
    if isfile(save_filename)
        if FLAG.OverwriteData.GenData
            RUNSIM = true;
            delete(save_filename)
    %     else
    %         Plant = load(save_filename);
    %         RESULTS.Slink_plant = Plant.Plant_Data;
        end
    else
        RUNSIM = true;
    end

if RUNSIM
    %% CPCT Calc
    %     if length(fieldnames(ERROR_CALC)) == 0
    %         CPCT_calc = NA;
    %     else
    %         CPCT_calc = ERROR_CALC.CPC_infty_diag;
    %     end

    %% Get ROM
    switch FLAG.EstimatorModel
        case 1 % Matlab SS_DT
            [~ , est_sys] = getSS_System(SIM,N,P,FLAG);
        case 2 % Ho-Kalman
            FLAG.Analysis.PlotImp = 0;
            [est_sys] = getHoKalmanROM(SIM,N,P,FLAG,RESULTS);
    end
    C_des = est_sys.C;
    C_m = est_sys.C(P.cell_voltage,:);


    %% CPCT from DARE
    [~, P_infty] = AsymptoticPreCalcs(FLAG,SIM,est_sys);
    CPCT = C_des * P_infty * C_des';


    %% SVD Analysis
    obsv_r_meas = obsv(est_sys.A,C_m);
    [~,S_Orm,V_Orm] = svd(obsv_r_meas);

    [N_outputs_ROM,N_states_ROM] = size(C_des);

    % Loop through C_tilde
    sing_val = nan(1,N_outputs_ROM);
    for mm = 1:N_outputs_ROM
        C_tilde = C_des(mm,:);

        % SVD of \tilde{C}
        [~ , ~ , V_tilde] = svd(C_tilde);

        % SVD of Augmented Observability
        obsv_aug = obsv_r_meas*V_tilde(:,1); %!!!!!! Do I need to check that V_tilde is in the column space if it is only 1 row?
        [~ , S_hat , V_hat] = svd(obsv_aug);

        sing_val(mm) = S_hat(1,1);
    end
    sing_val_norm = sing_val / S_Orm(1,1);


    %% Angle Analysis
    matrixOFdot = zeros(N_outputs_ROM,N_states_ROM);
    matrixOFdot_deg = zeros(N_outputs_ROM,N_states_ROM);
    for i = 1:N_outputs_ROM
        for j = 1:N_states_ROM
            matrixOFdot(i,j) = dot(  C_des(i,:)/norm(C_des(i,:))  ,  V_Orm(:,j)  );
            matrixOFdot_deg(i,j) = acosd(dot(  C_des(i,:)/norm(C_des(i,:))  ,  V_Orm(:,j)  ));
            if matrixOFdot_deg(i,j) > 90
                matrixOFdot_deg(i,j) = 180 - matrixOFdot_deg(i,j) ;
            end
        end
    end


    %% Save Results
    Q_0 = SIM.Q_0;
    R_0 = SIM.R_0;
    Ts  = SIM.Ts;
    SOC = SIM.SOC;
    %
    %     Q_str    = num2str(Q_0);
    %     R_str    = num2str(R_0);
    %     Ts_str   = num2str(Ts);
    %     SOC_str  = num2str(SOC);
    %
    %     filename = ['CompData_Q' Q_str '_R' R_str '_Ts' Ts_str '_C' SOC_str '.mat'];
    %     if FLAG.QMode == 1 % Input Q
    %         path_str = [pwd filesep 'Results' filesep 'ComparisonData' filesep 'InputQ'];
    %     else % State Q
    %         path_str = [pwd filesep 'Results' filesep 'ComparisonData' filesep 'StateQ'];
    %     end
    %     save_filename = [path_str filesep filename];
    save(save_filename, 'Q_0', 'R_0' , 'Ts', 'SOC', 'sing_val_norm', 'CPCT','matrixOFdot_deg','matrixOFdot',"S_Orm","sing_val",'C_des','P')
end


end