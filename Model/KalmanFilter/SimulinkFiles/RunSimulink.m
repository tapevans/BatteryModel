%% Run Simulink
function [RESULTS] = RunSimulink(SIM,N,P,FLAG,RESULTS)
%% Check if the simulation results exist
    RUNSIM = false;
    [slink_filename] = getSlinkfilename(FLAG,SIM);
    if isfile(slink_filename)
        if FLAG.OverwriteData.Slink
            RUNSIM = true;
            delete(slink_filename)
        else
            Plant = load(slink_filename);
            RESULTS.Slink_plant = Plant.Plant_Data;
        end
    else
        RUNSIM = true;
    end


%% Run Simulation
if RUNSIM
%% Add Battery Model to path
addpath(genpath('F:\TylerFiles\GitHubRepos\BatteryModel\Model'))


%% Change Working Directory
    [current_file_path,~,~] = fileparts(mfilename('fullpath'));
    oldFolder = cd(current_file_path); % Change to the folder that contains the Simulink files. 
    % This is mostly so the cache files don't save to the main workspace




%% Calculate LC parameters
    [LC_param] = getLCParam(SIM,N,P,FLAG);
    M_cross = LC_param.M_cross;
    U_NV_NT = LC_param.U_NV_NT;
    K_plant = LC_param.K;


%% Get other input parameters
    SV_IC = SIM.x_0_FOM;
    switch FLAG.QMode
        case 1 % Input Q
            mdl = 'Plant_Overall';
            pow_Q = SIM.pow_Qi;
        case 2 % State Q
            mdl = 'Plant_Overall_StateQ'; %### Not implemented
            pow_Q = SIM.pow_Qs;
    end
    C_des = SIM.OutputMatrix;

    % Estimator System
    switch FLAG.EstimatorModel
        case 1
            [~ , sys_DT] = getSS_System(SIM,N,P,FLAG);
            est_sys = sys_DT;
        case 2
            [est_sys] = getHoKalmanROM(SIM,N,P,FLAG);
    end
    A_ROM = est_sys.A;
    B_ROM = est_sys.B;
    C_des_ROM = est_sys.C;

    C_m     = C_des(P.cell_voltage,:);
    C_m_ROM = C_des_ROM(P.cell_voltage,:);

    x_0_ROM = zeros(length(A_ROM),1);


%% Get sim Structs
filename = getImpulseFilename(FLAG);
simsys = load(filename);


%% Function Handle Conversion
    simsys.AN.EqPotentialHandle = func2str(simsys.AN.EqPotentialHandle);
    simsys.AN.i_oHandle         = func2str(simsys.AN.i_oHandle);
    simsys.AN.sigmaHandle       = func2str(simsys.AN.sigmaHandle);
    simsys.AN.D_oHandle         = func2str(simsys.AN.D_oHandle);
    
    simsys.CA.EqPotentialHandle = func2str(simsys.CA.EqPotentialHandle);
    simsys.CA.i_oHandle         = func2str(simsys.CA.i_oHandle);
    simsys.CA.sigmaHandle       = func2str(simsys.CA.sigmaHandle);
    simsys.CA.D_oHandle         = func2str(simsys.CA.D_oHandle);
    
    simsys.EL.tf_numHandle      = func2str(simsys.EL.tf_numHandle);
    simsys.EL.ActivityHandle    = func2str(simsys.EL.ActivityHandle);
    simsys.EL.D_o_Li_ionHandle  = func2str(simsys.EL.D_o_Li_ionHandle);
    simsys.EL.kappaHandle       = func2str(simsys.EL.kappaHandle);
    
    simsys.SIM = rmfield(simsys.SIM,'fsolve_options');
    simsys.SIM = rmfield(simsys.SIM,'ControllerHandle');
    simsys.SIM = rmfield(simsys.SIM,'Controller_MO_File');


%% Run Simulink
    load_system(mdl)
    in = Simulink.SimulationInput(mdl);
    in = in.setModelParameter('StartTime','0','StopTime',num2str(SIM.t_final_Relax_Step));
    
    % Assign variable values by modifying the workspace
    mdlWks = get_param(in,'ModelWorkspace');
    
    assignin(mdlWks,'InputSignal' ,SIM.InputSignal)
    assignin(mdlWks,'AN'          ,simsys.AN)
    assignin(mdlWks,'CA'          ,simsys.CA)
    assignin(mdlWks,'SEP'         ,simsys.SEP)
    assignin(mdlWks,'EL'          ,simsys.EL)
    assignin(mdlWks,'SIM'         ,simsys.SIM)
    assignin(mdlWks,'CONS'        ,simsys.CONS)
    assignin(mdlWks,'P'           ,simsys.P)
    assignin(mdlWks,'N'           ,simsys.N)
    assignin(mdlWks,'PROPS'       ,simsys.PROPS)
    assignin(mdlWks,'FLAG'        ,simsys.FLAG)
    assignin(mdlWks,'C_m'         ,C_m)
    assignin(mdlWks,'C_des'       ,C_des)
    assignin(mdlWks,'A_ROM'       ,A_ROM)
    assignin(mdlWks,'B_ROM'       ,B_ROM)
    assignin(mdlWks,'C_m_ROM'     ,C_m_ROM)
    assignin(mdlWks,'C_des_ROM'   ,C_des_ROM)
    assignin(mdlWks,'x_0_ROM'     ,x_0_ROM)
    assignin(mdlWks,'M_cross'     ,M_cross)
    assignin(mdlWks,'U_NV_NT'     ,U_NV_NT)
    assignin(mdlWks,'K'           ,K_plant)
    assignin(mdlWks,'SV_IC'       ,SV_IC)
    assignin(mdlWks,'pow_Qi'      ,pow_Q)
    assignin(mdlWks,'pow_R'       ,SIM.pow_R)
    assignin(mdlWks,'Ts'          ,SIM.Ts)
    assignin(mdlWks,'tc'          ,SIM.tc)

    tic
    out = sim(in);
    toc
%     save_system
%     close_system

    
%% Save Simulink Variables
    %%%%%%%%%%%% States of Ideal Plant %%%%%%%%%%%% 
    if length(size(out.x_plant.signals.values)) >2
        for i = 1:length(out.x_plant.time)
            x_plant.value(:,i) = reshape(out.x_plant.signals.values(:,:,i),[],1);
        end
    else
        x_plant.value = out.x_plant.signals.values;
    end
    x_plant.time = out.x_plant.time;
    RESULTS.Slink_plant.x.t_soln = x_plant.time;
    RESULTS.Slink_plant.x.x_soln = x_plant.value;

    %%%%%%%%%%%% Measurement of Ideal Plant %%%%%%%%%%%% 
    if length(size(out.z_plant.signals.values)) >2
        for i = 1:length(out.z_plant.time)
            z_plant.value(:,i) = reshape(out.z_plant.signals.values(:,:,i),[],1);
        end
    else
        z_plant.value = out.z_plant.signals.values;
    end
    z_plant.time = out.z_plant.time;
    RESULTS.Slink_plant.z.t_soln = z_plant.time;
    RESULTS.Slink_plant.z.z_soln = z_plant.value;

    %%%%%%%%%%%% Measurement of All Outputs of Ideal Plant %%%%%%%%%%%% 
    if length(size(out.z_plant_All.signals.values)) >2
        for i = 1:length(out.z_plant_All.time)
            z_plant_All.value(:,i) = reshape(out.z_plant_All.signals.values(:,:,i),[],1);
        end
    else
        z_plant_All.value = out.z_plant_All.signals.values;
    end
    z_plant_All.time = out.z_plant_All.time;
    RESULTS.Slink_plant.zA.t_soln = z_plant_All.time;
    RESULTS.Slink_plant.zA.z_soln = z_plant_All.value;

    %%%%%%%%%%%% States of Input Noise Plant %%%%%%%%%%%% 
    if length(size(out.x_plant_IN.signals.values)) >2
        for i = 1:length(out.x_plant_IN.time)
            x_plant_IN.value(:,i) = reshape(out.x_plant_IN.signals.values(:,:,i),[],1);
        end
    else
        x_plant_IN.value = out.x_plant_IN.signals.values;
    end
    x_plant_IN.time = out.x_plant_IN.time;
    RESULTS.Slink_plant.xN.t_soln = x_plant_IN.time;
    RESULTS.Slink_plant.xN.x_soln = x_plant_IN.value;
    
    %%%%%%%%%%%% Measurement of Input Noise Plant %%%%%%%%%%%% 
    if length(size(out.z_plant_IN.signals.values)) >2
        for i = 1:length(out.z_plant_IN.time)
            z_plant_IN.value(:,i) = reshape(out.z_plant_IN.signals.values(:,:,i),[],1);
        end
    else
        z_plant_IN.value = out.z_plant_IN.signals.values;
    end
    z_plant_IN.time = out.z_plant_IN.time;
    RESULTS.Slink_plant.zN.t_soln = z_plant_IN.time;
    RESULTS.Slink_plant.zN.z_soln = z_plant_IN.value;
    
    %%%%%%%%%%%% Measurement of All Outputs of Input Noise Plant %%%%%%%%%%%% 
    if length(size(out.z_plant_All_IN.signals.values)) >2
        for i = 1:length(out.z_plant_All_IN.time)
            z_plant_All_IN.value(:,i) = reshape(out.z_plant_All_IN.signals.values(:,:,i),[],1);
        end
    else
        z_plant_All_IN.value = out.z_plant_All_IN.signals.values;
    end
    z_plant_All_IN.time = out.z_plant_All_IN.time;
    RESULTS.Slink_plant.zNA.t_soln = z_plant_All_IN.time;
    RESULTS.Slink_plant.zNA.z_soln = z_plant_All_IN.value;
            
    %%%%%%%%%%%% Process Noise Applied to Input Noise Plant %%%%%%%%%%%% 
    if length(size(out.w_k_IN.signals.values)) >2
        for i = 1:length(out.w_k_IN.time)
            w_k_IN.value(:,i) = reshape(out.w_k_IN.signals.values(:,:,i),[],1);
        end
    else
        w_k_IN.value = out.w_k_IN.signals.values;
    end
    w_k_IN.time = out.w_k_IN.time;
    RESULTS.Slink_plant.w.t_soln = w_k_IN.time;
    RESULTS.Slink_plant.w.w_soln = w_k_IN.value;
    
    %%%%%%%%%%%% Measurement Noise Applied to Input Noise Plant %%%%%%%%%%%% 
    if length(size(out.v_k_IN.signals.values)) >2
        for i = 1:length(out.v_k_IN.time)
            v_k_IN.value(:,i) = reshape(out.v_k_IN.signals.values(:,:,i),[],1);
        end
    else
        v_k_IN.value = out.v_k_IN.signals.values;
    end
    v_k_IN.time = out.v_k_IN.time;
    RESULTS.Slink_plant.v.t_soln = v_k_IN.time;
    RESULTS.Slink_plant.v.v_soln = v_k_IN.value;
    
    %%%%%%%%%%%% Pure Input Signal to All Plants %%%%%%%%%%%% 
    if length(size(out.i_user_in.signals.values)) >2
        for i = 1:length(out.i_user_in.time)
            i_user_in.value(:,i) = reshape(out.i_user_in.signals.values(:,:,i),[],1);
        end
    else
        i_user_in.value = out.i_user_in.signals.values;
    end
    i_user_in.time = out.i_user_in.time;
    RESULTS.Slink_plant.i_user.t_soln  = i_user_in.time;
    RESULTS.Slink_plant.i_user.i_user_soln = i_user_in.value;

    %%%%%%%%%%%% States of Input Noise ROM %%%%%%%%%%%% 
    if length(size(out.x_ROM_IN.signals.values)) >2
        for i = 1:length(out.x_ROM_IN.time)
            x_ROM_IN.value(:,i) = reshape(out.x_ROM_IN.signals.values(:,:,i),[],1);
        end
    else
        x_ROM_IN.value = out.x_ROM_IN.signals.values;
    end
    x_ROM_IN.time = out.x_ROM_IN.time;
    RESULTS.Slink_plant.xNR.t_soln = x_ROM_IN.time;
    RESULTS.Slink_plant.xNR.x_soln = x_ROM_IN.value;
    
    %%%%%%%%%%%% Measurement of Input Noise ROM %%%%%%%%%%%% 
    if length(size(out.z_ROM_IN.signals.values)) >2
        for i = 1:length(out.z_ROM_IN.time)
            z_ROM_IN.value(:,i) = reshape(out.z_ROM_IN.signals.values(:,:,i),[],1);
        end
    else
        z_ROM_IN.value = out.z_ROM_IN.signals.values;
    end
    z_ROM_IN.time = out.z_ROM_IN.time;
    RESULTS.Slink_plant.zNR.t_soln = z_ROM_IN.time;
    RESULTS.Slink_plant.zNR.z_soln = z_ROM_IN.value;
    
    %%%%%%%%%%%% Measurement of All Outputs of Input Noise ROM %%%%%%%%%%%% 
    if length(size(out.z_ROM_All_IN.signals.values)) >2
        for i = 1:length(out.z_ROM_All_IN.time)
            z_ROM_All_IN.value(:,i) = reshape(out.z_ROM_All_IN.signals.values(:,:,i),[],1);
        end
    else
        z_ROM_All_IN.value = out.z_ROM_All_IN.signals.values;
    end
    z_ROM_All_IN.time = out.z_ROM_All_IN.time;
    RESULTS.Slink_plant.zNRA.t_soln = z_ROM_All_IN.time;
    RESULTS.Slink_plant.zNRA.z_soln = z_ROM_All_IN.value;

    
%% Return to Old Working Directory
    cd(oldFolder);


end
end

%%%%%%%%%%%%%%%%%%%%%%%%%% OOOOOOOOOOOLD %%%%%%%%%%%%%%%%%%%%%%%%%% 
%% InputSignalMode
%  1) From code
%  2) From function input
%     FLAG.InputSignalMode = 2;


%     old_C_mode = FLAG.C_mode;
%     FLAG.C_mode = 5;
%     N.states = 4;
%     [A,B,C,~] = getAll_SS(SIM,N,P,FLAG);
%     [InputSignal] = getInputSignal(SIM,N,P,FLAG);


%     FLAG.C_mode = old_C_mode;

    % Get Measured for both FOM and ROM   
%     switch FLAG.C_mode
%         case N.states+1
%             C_m     = C_des;
%             C_m_ROM = C_des_ROM;
%         otherwise
%     end
            
%     y_0 = C_des*SIM.x_0;

%     assignin(mdlWks,'A'           ,A)
%     assignin(mdlWks,'B'           ,B)