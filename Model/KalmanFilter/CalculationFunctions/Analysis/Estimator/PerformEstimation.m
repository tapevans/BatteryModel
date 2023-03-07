%% Estimator

function [ESTIMATOR,RESULTS] = PerformEstimation(plant_filename,SIM,FLAG,N,P,RESULTS)
%% Load Plant Data
    if ~FLAG.UseROMAsPlant % Load from Slink
        Plant = load(plant_filename);
        
        % Plant Data
        u = Plant.Plant_Data.i_user.i_user_soln;
        z = Plant.Plant_Data.zN.z_soln;
    else
        [u , z] = ROMplantData(SIM,FLAG,N,P,RESULTS);
    end
    
%     if FLAG.UseWrongIC
%         z_init = SIM.y_0_FOM_Offset;
%     else
        z_init = SIM.y_0_FOM;
%     end


%% Get Estimator Model
%  1) Matlab SS_DT
%  2) Ho-Kalman
    switch FLAG.EstimatorModel
        case 1
            [~ , sys_DT] = getSS_System(SIM,N,P,FLAG);
            est_sys_tot{1} = sys_DT;
        case 2
            [HK_sys] = getHoKalmanROM(SIM,N,P,FLAG);
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

    C_m = est_sys.C(P.cell_voltage,:);
    D_m = est_sys.D(P.cell_voltage,:);
    est_sys = ss(est_sys.A , est_sys.B , C_m , D_m , SIM.Ts);

    if FLAG.UseWrongIC
        random = randi([-9,9],length(est_sys.A),1);
        x_hat_0 = random * FLAG.offsetROM;
    else
        x_hat_0 = zeros(length(est_sys.A),1);
    end


%% Do Asymptotic Calcs
    [ESTIMATOR.K_infty{OO}, ESTIMATOR.P_infty{OO}] = AsymptoticPreCalcs(FLAG,SIM,est_sys);
    

%% Do Asymptotic Estimation
    [x_hat_asy]   = AsyEstimator(est_sys,FLAG,SIM,P,x_hat_0,u,z);
    if OO == 1
        y_hat_asy     = est_sys.C * x_hat_asy + z_init(P.cell_voltage,1);
    end
    if FLAG.EstimatorModel == 1 || FLAG.EST.SepHK == 0
        y_hat_asy_ALL = C_des_ROM * x_hat_asy + z_init;
    else
        y_hat_asy_ALL = C_des_ROM * x_hat_asy + z_init([1,OO],1);
    end
    
    if  OO == 1
        % RESULTS.EST.ASY.t_soln      = Plant.Plant_Data.i_user.t_soln;
        RESULTS.EST.ASY.t_soln      = u;
        RESULTS.EST.ASY.x_soln      = x_hat_asy;
        RESULTS.EST.ASY.z_soln      = y_hat_asy;
        RESULTS.EST.ASY.z_soln_ALL  = y_hat_asy_ALL;
    else
        %RESULTS.EST.ASY.t_soln      = Plant.Plant_Data.i_user.t_soln;
        %RESULTS.EST.ASY.x_soln(OO,:)      = x_hat_asy; % If I need this, I have to think more on how to store the data
        %RESULTS.EST.ASY.z_soln      = y_hat_asy;
        RESULTS.EST.ASY.z_soln_ALL(OO,:)  = y_hat_asy_ALL(2,:);
    end


%% Do Variable Estimation
    [x_hat_var, ESTIMATOR.K_k{OO}, ESTIMATOR.P_k_pre{OO}] = VarEstimator(est_sys,FLAG,SIM,P,x_hat_0,u,z);
    if OO == 1
        y_hat_var     = est_sys.C * x_hat_var + z_init(P.cell_voltage,1);
    end
    if FLAG.EstimatorModel == 1 || FLAG.EST.SepHK == 0
        y_hat_var_ALL = C_des_ROM * x_hat_var + z_init;
    else
        y_hat_var_ALL = C_des_ROM * x_hat_var + z_init([1,OO],1);
    end
    
    if  OO == 1
        %RESULTS.EST.VAR.t_soln      = Plant.Plant_Data.i_user.t_soln;
        RESULTS.EST.VAR.t_soln      = u;
        RESULTS.EST.VAR.x_soln      = x_hat_var;
        RESULTS.EST.VAR.z_soln      = y_hat_var;
        RESULTS.EST.VAR.z_soln_ALL  = y_hat_var_ALL;
    else
        %RESULTS.EST.VAR.t_soln      = Plant.Plant_Data.i_user.t_soln;
        %RESULTS.EST.VAR.x_soln(OO,:)      = x_hat_asy; % If I need this, I have to think more on how to store the data
        %RESULTS.EST.VAR.z_soln      = y_hat_asy;
        RESULTS.EST.VAR.z_soln_ALL(OO,:)  = y_hat_var_ALL(2,:);
    end
    
    
end    
SIM.Qi = oldQ;


%% Convert Plant true states into ROM states
%     [U, S, V]   = svd(C_des_ROM);
%     threshold   = 1e-7;
%     S_cross     = pinv(S,threshold);
%     C_des_cross = V*S_cross*U';
% 
%     x_plant_ROM = C_des_cross * Plant.Plant_Data.zNA.z_soln;
%     
%     RESULTS.EST.PLANT.t_soln = Plant.Plant_Data.i_user.t_soln;
%     RESULTS.EST.PLANT.x_soln = x_plant_ROM;
%     RESULTS.EST.PLANT.z_soln = Plant.Plant_Data.zN.z_soln;
%     RESULTS.EST.PLANT.z_soln_ALL = Plant.Plant_Data.zNA.z_soln;


    %RESULTS.EST.PLANT.t_soln     = Plant.Plant_Data.i_user.t_soln;
    RESULTS.EST.PLANT.t_soln      = u;
    RESULTS.EST.PLANT.x_soln     = Plant.Plant_Data.xN.x_soln;
    RESULTS.EST.PLANT.z_soln_ALL = SIM.OutputMatrix * Plant.Plant_Data.xN.x_soln;
    RESULTS.EST.PLANT.z_soln     = Plant.Plant_Data.zN.z_soln;


end

%% OLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reduce the outputs to C_mode
%     switch FLAG.C_mode
%         case N.states
%             C_m = est_sys.C;
%             D_m = est_sys.D;
%             est_sys = ss(est_sys.A,est_sys.B,C_m,D_m,SIM.Ts);
%         otherwise
%     end

% Initial Conditions
%     y_0 = SIM.x_0;
%     x_hat_0 = C_des_ROM\y_0;
%     x_hat_0(3) = 5.25;
%     x_hat_0(4) = -0.25;

% y_hat_var_ALL = C_des_ROM * x_hat_var + z_init;
% RESULTS.EST.VAR.t_soln      = Plant.Plant_Data.i_user.t_soln;
%     RESULTS.EST.VAR.x_soln      = x_hat_var;
%     RESULTS.EST.VAR.z_soln      = y_hat_var;
%     RESULTS.EST.VAR.z_soln_ALL  = y_hat_var_ALL;