%% Run Simulink
function [RESULTS] = RunSimulink(SIM,N,P,FLAG,RESULTS)
%% Check if the simulation results exist
RUNSIM = false;
[slink_filename] = getSlinkfilename(FLAG,SIM);
if isfile(slink_filename)
    if SIM.Analysis.NoisyPlant.OverwriteData
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
%% Change Working Directory
    [current_file_path,~,~] = fileparts(mfilename('fullpath'));
    oldFolder = cd(current_file_path); % Change to the folder that contains the Simulink files. 
    % This is mostly so the cache files don't save to the main workspace


%% InputSignalMode
%  1) From code
%  2) From function input
    FLAG.InputSignalMode = 2;


%% Calculate LC parameters
    [LC_param] = getLCParam(SIM,N,P,FLAG);
    M_cross = LC_param.M_cross;
    U_NV_NT = LC_param.U_NV_NT;
    K_plant = LC_param.K;


%% Get other input parameters
    switch FLAG.SlinkModel
        case 3 % 3 state plant
            old_C_mode = FLAG.C_mode;
            FLAG.C_mode = 4;
            N.states = 3;
            SV_IC = SIM.x_0_3;
            [A,B,C,~] = getAll_SS3(SIM,N,P,FLAG);
            [InputSignal] = getInputSignal(SIM,N,P,FLAG);
            switch FLAG.QMode
                case 1 % Input Q
                    mdl = 'Circuit_3State_Overall';
                    pow_Q = SIM.pow_Q3i;
                case 2 % State Q
                    mdl = 'Circuit_3State_Overall_StateQ';
                    pow_Q = SIM.pow_Q3s;
            end
            FLAG.C_mode = old_C_mode;
        case 5 % 5 state plant
            old_C_mode = FLAG.C_mode;
            FLAG.C_mode = 4;
            N.states = 5;
            SV_IC = SIM.x_0_5;
            [A,B,C,~] = getAll_SS5(SIM,N,P,FLAG);
            [InputSignal] = getInputSignal(SIM,N,P,FLAG);
            switch FLAG.QMode
                case 1 % Input Q
                    mdl = 'Circuit_5State_Overall';
                    pow_Q = SIM.pow_Q5i;
                case 2 % State Q
                    mdl = 'Circuit_5State_Overall_StateQ';
                    pow_Q = SIM.pow_Q5s;
            end
            FLAG.C_mode = old_C_mode;
    end

    C_des = C;
    switch FLAG.C_mode
        case 1
            C_m = C_des(1,:);
        case 2
            C_m = C_des(2,:);
        case 3
            C_m = C_des(3,:);
        case 4
            C_m = C_des;
    end

    % Estimator System
    old_C_mode = FLAG.C_mode;
    FLAG.C_mode = 4;

    switch FLAG.EstimatorModel
        case 1
            [~ , sys_DT3 , ~ , ~] = getSS_System(SIM,N,P,FLAG);
            est_sys = sys_DT3;
        case 2
            [~ , ~ , ~ , sys_DT5] = getSS_System(SIM,N,P,FLAG);
            est_sys = sys_DT5;
        case 3
            % Get Ho-Kalman ROM
    end

    FLAG.C_mode = old_C_mode;

    A_ROM = est_sys.A;
    B_ROM = est_sys.B;
    C_des_ROM = est_sys.C;
    switch FLAG.C_mode
        case 1
            C_m_ROM = C_des_ROM(1,:);
        case 2
            C_m_ROM = C_des_ROM(2,:);
        case 3
            C_m_ROM = C_des_ROM(3,:);
        case 4
            C_m_ROM = C_des_ROM;
    end
    
    switch FLAG.SlinkModel
        case 3 % 3 state plant
            y_0 = SIM.x_0_3;
        case 5
            y_0 = SIM.x_0_5;
    end
    
    x_0_ROM = C_des_ROM\y_0;


%% Run Simulink
    load_system(mdl)
    in = Simulink.SimulationInput(mdl);
    in = in.setModelParameter('StartTime','0','StopTime',num2str(SIM.t_final_sim));
    
    % Assign variable values by modifying the workspace
    mdlWks = get_param(in,'ModelWorkspace');
    
    assignin(mdlWks,'InputSignal' ,InputSignal)
    assignin(mdlWks,'SIM'         ,SIM)
    assignin(mdlWks,'N'           ,N)
    assignin(mdlWks,'P'           ,P)
    assignin(mdlWks,'FLAG'        ,FLAG)
    assignin(mdlWks,'A'           ,A)
    assignin(mdlWks,'B'           ,B)
    assignin(mdlWks,'C_m'         ,C_m)
    assignin(mdlWks,'C_des'       ,C_des)
    assignin(mdlWks,'A_ROM'       ,A_ROM)
    assignin(mdlWks,'B_ROM'       ,B_ROM)
    assignin(mdlWks,'C_m_ROM'     ,C_m_ROM)
    assignin(mdlWks,'C_des_ROM'   ,C_des_ROM)
    assignin(mdlWks,'x_0_ROM'     ,x_0_ROM)
    assignin(mdlWks,'M_cross'     ,M_cross)
    assignin(mdlWks,'U_NV_NT'     ,U_NV_NT)
    assignin(mdlWks,'K_plant'     ,K_plant)
    assignin(mdlWks,'SV_IC'       ,SV_IC)
    assignin(mdlWks,'pow_Q3i'     ,pow_Q)
    assignin(mdlWks,'pow_R'       ,SIM.pow_R)
    assignin(mdlWks,'Ts'          ,SIM.Ts)
    assignin(mdlWks,'tc'          ,SIM.tc)

    tic
    out = sim(in);
    toc
    save_system
    close_system

%% Save Simulink Variables
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

    %%%%%%%%%%%%%%%%%%%%%%
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

    %%%%%%%%%%%%%%%%%%%%%%
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

    %%%%%%%%%%%%%%%%%%%%%%
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
    
    %%%%%%%%%%%%%%%%%%%%%%
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
    
    %%%%%%%%%%%%%%%%%%%%%%
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
            
    %%%%%%%%%%%%%%%%%%%%%%
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
    
    %%%%%%%%%%%%%%%%%%%%%%
    if length(size(out.v_k_IN.signals.values)) >2
        for i = 1:length(out.v_k_IN.time)
            v_k_IN.value(:,i) = reshape(out.v_k_IN.signals.values(:,:,i),[],1);
        end
    else
        v_k_IN.value = out.v_k_IN.signals.values;
    end
    v_k_IN.time = out.v_k_IN.time;
    RESULTS.Slink_plant.v.t_soln = v_k_IN.time;
    RESULTS.Slink_plant.v.w_soln = v_k_IN.value;
    
    %%%%%%%%%%%%%%%%%%%%%%
    if length(size(out.V_in.signals.values)) >2
        for i = 1:length(out.V_in.time)
            V_in.value(:,i) = reshape(out.V_in.signals.values(:,:,i),[],1);
        end
    else
        V_in.value = out.V_in.signals.values;
    end
    V_in.time = out.V_in.time;
    RESULTS.Slink_plant.Vs.t_soln  = V_in.time;
    RESULTS.Slink_plant.Vs.Vs_soln = V_in.value;

    %%%%%%%%%%%%%%%%%%%%%%
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
    
    %%%%%%%%%%%%%%%%%%%%%%
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
    
    %%%%%%%%%%%%%%%%%%%%%%
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