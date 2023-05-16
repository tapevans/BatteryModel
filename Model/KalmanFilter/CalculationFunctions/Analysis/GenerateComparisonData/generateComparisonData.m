%% Generate Comparison Data
function generateComparisonData(SIM,FLAG,N,P,RESULTS) 
%% Check if the results exist
    RUNSIM = false;
    save_filename = getGenDataFilename(FLAG);
    if isfile(save_filename)
        if FLAG.OverwriteData.GenData
            RUNSIM = true;
            delete(save_filename)
            disp('Deleting existing GenData')
        else
            disp('GenData already exists')
        end
    else
        RUNSIM = true;
    end


%% Run data generation if needed
    if RUNSIM
        %% Get ROM
        switch FLAG.EstimatorModel
            case 1
                [~ , sys_DT] = getSS_System(SIM,N,P,FLAG);
                est_sys_tot{1} = sys_DT;
            case 2
                [HK_sys] = getHoKalmanROM(SIM,N,P,FLAG,RESULTS);
                % if FLAG.EST.SepHK
                %     est_sys_tot = HK_sys(1:N.DesOut);
                % else
                %     est_sys_tot{1} = HK_sys{end};
                % end
                est_sys_tot = HK_sys;
        end
    

        %% Loop through all ROMs
        for OO = 1:length(est_sys_tot)
            %% Get Single ROM
            est_sys = est_sys_tot{OO};
            % C_des   = est_sys.C;
            C_des_temp = est_sys.C;
            C_m        = est_sys.C(P.cell_voltage,:);


            %% CPCT from DARE
            [~, P_infty] = AsymptoticPreCalcs(FLAG,SIM,est_sys);
            % CPCT = C_des * P_infty * C_des';
            CPCT{OO} = C_des_temp * P_infty * C_des_temp';
    
    
            %% SVD Analysis
            obsv_r_meas = obsv(est_sys.A,C_m);
            [~,S_Orm_temp,V_Orm] = svd(obsv_r_meas);
    
            [N_outputs_ROM,N_states_ROM] = size(C_des_temp);
    
            % Loop through C_tilde
            sing_val_temp = nan(1,N_outputs_ROM);
            for mm = 1:N_outputs_ROM
                C_tilde = C_des_temp(mm,:);
    
                % SVD of \tilde{C}
                [~ , ~ , V_tilde] = svd(C_tilde);
    
                % SVD of Augmented Observability
                obsv_aug = obsv_r_meas*V_tilde(:,1); %!!!!!! Do I need to check that V_tilde is in the column space if it is only 1 row?
                [~ , S_hat , V_hat] = svd(obsv_aug);
    
                sing_val_temp(mm) = S_hat(1,1);
            end
            % sing_val_norm = sing_val / S_Orm(1,1);
            sing_val_norm{OO} = sing_val_temp / S_Orm_temp(1,1);
    
    
            %% Angle Analysis
            % matrixOFdot     = zeros(N_outputs_ROM,N_states_ROM);
            % matrixOFdot_deg = zeros(N_outputs_ROM,N_states_ROM);
            matrixOFdot_temp     = zeros(N_outputs_ROM,N_states_ROM);
            matrixOFdot_deg_temp = zeros(N_outputs_ROM,N_states_ROM);
            for i = 1:N_outputs_ROM
                for j = 1:N_states_ROM
                    matrixOFdot_temp(i,j)     = dot(  C_des_temp(i,:)/norm(C_des_temp(i,:))  ,  V_Orm(:,j)  );
                    matrixOFdot_deg_temp(i,j) = acosd(dot(  C_des_temp(i,:)/norm(C_des_temp(i,:))  ,  V_Orm(:,j)  ));
                    if matrixOFdot_deg_temp(i,j) > 90
                        matrixOFdot_deg_temp(i,j) = 180 - matrixOFdot_deg_temp(i,j) ;
                    end
                end
            end

            matrixOFdot_deg{OO} = matrixOFdot_deg_temp;
            matrixOFdot{OO}     = matrixOFdot_temp;
            S_Orm{OO}           = S_Orm_temp;
            sing_val{OO}        = sing_val_temp;
            C_des{OO}           = C_des_temp;
        end

        %% PRBS Results
            OutputMatrix = SIM.OutputMatrix;
            [r,~] = size(RESULTS.ode.z_soln);
            for OO = 1:r
                ODEProfile{OO} = RESULTS.ode.z_soln(OO,:); % Can do IC and delta from max-min later
            end
        

        %% Save Results
            Q_0 = SIM.Q_0;
            R_0 = SIM.R_0;
            Ts  = SIM.Ts;
            SOC = SIM.SOC;
        
            save(save_filename, 'Q_0', 'R_0' , 'Ts', 'SOC', 'CPCT', 'sing_val_norm' ,'matrixOFdot_deg','matrixOFdot',"S_Orm","sing_val",'C_des','P','OutputMatrix','ODEProfile')
    end
end

%%%% OOOOOOOOOOOOOOOOOOOOOOOOOOLLLLLLD

% ,ESTIMATOR,ERROR_CALC


        %% CPCT Calc
        %     if length(fieldnames(ERROR_CALC)) == 0
        %         CPCT_calc = NA;
        %     else
        %         CPCT_calc = ERROR_CALC.CPC_infty_diag;
        %     end


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