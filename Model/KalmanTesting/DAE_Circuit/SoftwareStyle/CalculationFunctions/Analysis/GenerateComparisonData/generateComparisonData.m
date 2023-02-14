%% Generate Comparison Data
function generateComparisonData(SIM,FLAG,N,P,RESULTS,ESTIMATOR,ERROR_CALC)
%% %%%%%%%%%%%%%%%Need something to check if data file already exists!!!!!!!!!!!!!!

%% CPCT Calc
    CPCT_calc = ERROR_CALC.CPC_infty_diag;

 %% SVD Analysis
    % Observability Measurable
    %     obsv_r_meas = getObservability(A_DT,C_DT(1,:));

switch FLAG.C_mode
    case 1
        C_DT = ESTIMATOR.sys.C(1,:);
    case 2
        C_DT = ESTIMATOR.sys.C(2,:);
    case 3
        C_DT = ESTIMATOR.sys.C(3,:);
    case 4
        C_DT = ESTIMATOR.sys.C;
end


    obsv_r_meas = obsv(ESTIMATOR.sys.A,C_DT);
    [~,S_Orm,V_Orm] = svd(obsv_r_meas);

    % C_tilde (This normally would be of the ROM)
    C_r = ESTIMATOR.C_des_ROM;

    [N_outputs,N_states] = size(C_r);

    % Loop through C_tilde
    sing_val = nan(1,N_outputs);
    for mm = 1:N_outputs
        C_tilde = C_r(mm,:);

        % SVD of \tilde{C}
        [~ , ~ , V_tilde] = svd(C_tilde);

        % SVD of Augmented Observability
        obsv_aug = obsv_r_meas*V_tilde(:,1); %!!!!!! Do I need to check that V_tilde is in the column space if it is only 1 row?
        [~ , S_hat , V_hat] = svd(obsv_aug);

        sing_val(mm) = S_hat(1,1);
    end
    sing_val_norm = sing_val / S_Orm(1,1);


%     matrixOFdot = zeros(N_outputs,N_states);
    matrixOFdot_deg = zeros(N_outputs,N_states);
    for i = 1:N_outputs
        for j = 1:N_states
%             matrixOFdot(i,j) = dot(  C_r(i,:)/norm(C_r(i,:))  ,  V_Orm(:,j)  );
            matrixOFdot_deg(i,j) = acosd(dot(  C_r(i,:)/norm(C_r(i,:))  ,  V_Orm(:,j)  ));
            if matrixOFdot_deg(i,j) > 90
                matrixOFdot_deg(i,j) = 180 - matrixOFdot_deg(i,j) ;
            end
        end
    end


%% Save Results
    Q_0 = SIM.Q_0;
    R_0 = SIM.R_0;
    Ts = SIM.Ts;
    C_idx = FLAG.C_mode;
    
    Q_str = num2str(Q_0);
    R_str = num2str(R_0);
    Ts_str = num2str(Ts);
    C_str = num2str(C_idx);
    
    filename = ['CompData_Q' Q_str '_R' R_str '_Ts' Ts_str '_C' C_str '.mat'];
    path_str = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanTesting\DAE_Circuit\SoftwareStyle\Results\ComparisonData';
    save_filename = [path_str filesep filename];
    save(save_filename, 'Q_0', 'R_0' , 'Ts', 'C_idx', 'sing_val_norm', 'CPCT_calc','matrixOFdot_deg',"S_Orm","sing_val")

end