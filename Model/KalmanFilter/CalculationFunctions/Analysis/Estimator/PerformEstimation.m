%% Estimator

function [ESTIMATOR,RESULTS] = PerformEstimation(plant_filename,SIM,FLAG,N,P,RESULTS)
%% Load Plant Data
    if ~FLAG.UseROMAsPlant % Load from Slink
        if isfile(plant_filename)
            Plant = load(plant_filename);
        else
            disp('Data for this plant has not run in Simulink')
        end
        
        % Plant Data
        u = Plant.Plant_Data.i_user.i_user_soln;
        u = [u(2:end) ; u(end)]; %% Fix the delay
        z = Plant.Plant_Data.zN.z_soln;
        t = Plant.Plant_Data.i_user.t_soln;
    else
        if FLAG.Tswitch == 100 && FLAG.UseROMAsPlant % Troublesome case
            [t, u , ~ , z_all, z_cell_voltage] = ROMplantData(SIM,FLAG,N,P,RESULTS);
        else
            [t, u , x , z_all, z_cell_voltage] = ROMplantData(SIM,FLAG,N,P,RESULTS);
        end
    end
    
    z_init = SIM.y_0_FOM;


%% Get Estimator Model
%  1) Matlab SS_DT
%  2) Ho-Kalman
    switch FLAG.EstimatorModel
        case 1
            [~ , sys_DT] = getSS_System(SIM,N,P,FLAG);
            est_sys_tot{1} = sys_DT;
        case 2
            [HK_sys] = getHoKalmanROM(SIM,N,P,FLAG,RESULTS);
            if FLAG.EST.SepHK
                est_sys_tot = HK_sys(1:N.DesOut);
            else
                est_sys_tot{1} = HK_sys{end};
            end
    end

    oldQ = SIM.Qi;
    SIM.Qi = SIM.Qi * SIM.LargeQMultiply;
for OO = 1:length(est_sys_tot)
    est_sys = est_sys_tot{OO};
    C_des_ROM = est_sys.C;
    ESTIMATOR.C_des_ROM{OO} = C_des_ROM;
    ESTIMATOR.sys{OO} = est_sys;

    % C_m = est_sys.C(P.cell_voltage,:);
    % D_m = est_sys.D(P.cell_voltage,:);
    % est_sys = ss(est_sys.A , est_sys.B , C_m , D_m , SIM.Ts);

    if FLAG.UseWrongIC
        random = randi([-9,9],length(est_sys.A),1);
        x_hat_0 = random * FLAG.offsetROM;
    else
        x_hat_0 = zeros(length(est_sys.A),1);
    end
    
    if FLAG.UseROMAsPlant
        z = z_cell_voltage{OO}(P.cell_voltage,:);
    end
    

%% Do Asymptotic Estimation
    if FLAG.DoAsy
        [x_hat_asy]   = AsyEstimator(est_sys,FLAG,SIM,P,x_hat_0,u,z);
        if FLAG.EstimatorModel == 1 || FLAG.EST.SepHK == 0
            y_hat_asy_ALL = C_des_ROM * x_hat_asy + z_init;
        else
            y_hat_asy_ALL = C_des_ROM * x_hat_asy + z_init([P.cell_voltage,OO],1);
        end
        y_hat_asy = y_hat_asy_ALL(P.cell_voltage,:);
        
        RESULTS.EST.ASY.x_soln{OO}      = x_hat_asy;
        RESULTS.EST.ASY.z_soln{OO}      = y_hat_asy;
        RESULTS.EST.ASY.z_soln_ALL{OO}  = y_hat_asy_ALL;
    end

%% Do Variable Estimation
    if FLAG.DoVar
        if FLAG.Tswitch == 100 && FLAG.UseROMAsPlant
            [x_hat_var, ~, ~] = VarEstimator(est_sys,FLAG,SIM,P,x_hat_0,u,z);
        else
            [x_hat_var, ESTIMATOR.K_k{OO}, ESTIMATOR.P_k_pre{OO}] = VarEstimator(est_sys,FLAG,SIM,P,x_hat_0,u,z);
        end

        if FLAG.EstimatorModel == 1 || FLAG.EST.SepHK == 0
            y_hat_var_ALL = C_des_ROM * x_hat_var + z_init;
        else
            y_hat_var_ALL = C_des_ROM * x_hat_var + z_init([P.cell_voltage,OO],1);
        end
        y_hat_var = y_hat_var_ALL(P.cell_voltage,:);
        
        if FLAG.Tswitch == 100 && FLAG.UseROMAsPlant % Troublesome case
        else
            RESULTS.EST.VAR.x_soln{OO}      = x_hat_var;
        end
        RESULTS.EST.VAR.z_soln{OO}      = y_hat_var;
        RESULTS.EST.VAR.z_soln_ALL{OO}  = y_hat_var_ALL;
    end  


%% Do Asymptotic Calcs
    if FLAG.DoPreCalc
        [ESTIMATOR.K_infty{OO}, ESTIMATOR.P_infty{OO}] = AsymptoticPreCalcs(FLAG,SIM,est_sys);
    else
        if FLAG.Tswitch == 100 && FLAG.UseROMAsPlant % Troublesome case
        else
            ESTIMATOR.K_infty{OO} = ESTIMATOR.K_k{OO}(:,:,end);
            ESTIMATOR.P_infty{OO} = ESTIMATOR.P_k_pre{OO}(:,:,end);
        end
    end  
    
end    


%% Something
    SIM.Qi = oldQ;
    RESULTS.EST.ASY.t_soln   = t;
    RESULTS.EST.VAR.t_soln   = t;
    RESULTS.EST.PLANT.t_soln = t;


%% Save Plant Results
    if ~FLAG.UseROMAsPlant % Load from Slink
        if FLAG.EstimatorModel == 1 || FLAG.EST.SepHK == 0
            if FLAG.Tswitch == 100 && FLAG.UseROMAsPlant % Troublesome case
            else
                RESULTS.EST.PLANT.x_soln{1}     = Plant.Plant_Data.xN.x_soln;
            end
            RESULTS.EST.PLANT.z_soln_ALL{1} = SIM.OutputMatrix * Plant.Plant_Data.xN.x_soln;
            RESULTS.EST.PLANT.z_soln{1}     = Plant.Plant_Data.zN.z_soln;
        else
            plant_z_all = SIM.OutputMatrix * Plant.Plant_Data.xN.x_soln;
            for OO = 1:length(est_sys_tot)
                if FLAG.Tswitch == 100 && FLAG.UseROMAsPlant
                else
                    RESULTS.EST.PLANT.x_soln{OO}     = Plant.Plant_Data.xN.x_soln;
                end
                RESULTS.EST.PLANT.z_soln_ALL{OO} = plant_z_all([P.cell_voltage,OO],:);
                RESULTS.EST.PLANT.z_soln{OO}     = Plant.Plant_Data.zN.z_soln;
            end
        end
    else
        if FLAG.Tswitch == 100 && FLAG.UseROMAsPlant % Troublesome case
        else
            RESULTS.EST.PLANT.x_soln     = x;
        end
        
        RESULTS.EST.PLANT.z_soln_ALL     = z_all;
        
        for OO = 1:length(est_sys_tot)
            RESULTS.EST.PLANT.z_soln{OO} = z_cell_voltage{OO}(P.cell_voltage,:);
        end
    end


end