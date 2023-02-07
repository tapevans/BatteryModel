%% Run Simulink
function [t,x,z,u,w,v] = RunSimulink(SIM,N,P,FLAG)
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


%% Run Simulink - %%%%%%%%%%%%%%%%%% Have stuff for different sim modes
    switch FLAG.SlinkModel
        case 1 % 3 state plant
            N.states = 3;
            SV_IC = SIM.x_0_3;
            mdl = 'Circuit_3State_Overall';
            FLAG.C_mode = 4;
            [A,B,C,~] = getAll_SS3(SIM,N,P,FLAG);
            FLAG.InputType = 3; %No Ramp
            [InputSignal] = getInputSignal(SIM,N,P,FLAG);
            

        case 2 % 3 state plant with input noise
            N.states = 3;
            mdl = 'Circuit_3State_Overall';
            SV_IC = SIM.x_0_3;
            FLAG.C_mode = 4;
            [A,B,C,~] = getAll_SS3(SIM,N,P,FLAG);
            FLAG.InputType = 3; %No Ramp
            [InputSignal] = getInputSignal(SIM,N,P,FLAG);

        case 3 % 3 state plant with state noise
            N.states = 3;
            mdl = 'Circuit_3State_Overall';
            SV_IC = SIM.x_0_3;
            FLAG.C_mode = 4;
            [A,B,C,~] = getAll_SS3(SIM,N,P,FLAG);
            FLAG.InputType = 3; %No Ramp
            [InputSignal] = getInputSignal(SIM,N,P,FLAG);
    end

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
    assignin(mdlWks,'C_m'         ,C)
    assignin(mdlWks,'M_cross'     ,M_cross)
    assignin(mdlWks,'U_NV_NT'     ,U_NV_NT)
    assignin(mdlWks,'K_plant'     ,K_plant)
    assignin(mdlWks,'SV_IC'       ,SV_IC)
    assignin(mdlWks,'pow_Q3i'     ,SIM.pow_Q3i)
    assignin(mdlWks,'pow_R'       ,SIM.pow_R)
    assignin(mdlWks,'Ts'          ,SIM.Ts)
    assignin(mdlWks,'tc'          ,SIM.tc)

    tic
    out = sim(in);
    toc
%     save_system
%     close_system

%% Save Simulink Variables

    if length(size(out.x_plant.signals.values)) >2
        for i = 1:length(out.x_plant.time)
            x_plant.value(:,i) = reshape(out.x_plant.signals.values(:,:,i),[],1);
        end
    else
        x_plant.value = out.x_plant.signals.values;
    end
    x_plant.time = out.x_plant.time;
    
    if length(size(out.z_plant.signals.values)) >2
        for i = 1:length(out.z_plant.time)
            z_plant.value(:,i) = reshape(out.z_plant.signals.values(:,:,i),[],1);
        end
    else
        z_plant.value = out.z_plant.signals.values;
    end
    z_plant.time = out.z_plant.time;

    if length(size(out.x_plant_IN.signals.values)) >2
        for i = 1:length(out.x_plant_IN.time)
            x_plant_IN.value(:,i) = reshape(out.x_plant_IN.signals.values(:,:,i),[],1);
        end
    else
        x_plant_IN.value = out.x_plant_IN.signals.values;
    end
    x_plant_IN.time = out.x_plant_IN.time;
    
    if length(size(out.z_plant_IN.signals.values)) >2
        for i = 1:length(out.z_plant_IN.time)
            z_plant_IN.value(:,i) = reshape(out.z_plant_IN.signals.values(:,:,i),[],1);
        end
    else
        z_plant_IN.value = out.z_plant_IN.signals.values;
    end
    z_plant_IN.time = out.z_plant_IN.time;
                
    if length(size(out.w_k_IN.signals.values)) >2
        for i = 1:length(out.w_k_IN.time)
            w_k_IN.value(:,i) = reshape(out.w_k_IN.signals.values(:,:,i),[],1);
        end
    else
        w_k_IN.value = out.w_k_IN.signals.values;
    end
    w_k_IN.time = out.w_k_IN.time;
    
    if length(size(out.v_k_IN.signals.values)) >2
        for i = 1:length(out.v_k_IN.time)
            v_k_IN.value(:,i) = reshape(out.v_k_IN.signals.values(:,:,i),[],1);
        end
    else
        v_k_IN.value = out.v_k_IN.signals.values;
    end
    v_k_IN.time = out.v_k_IN.time;
    
    if length(size(out.V_in.signals.values)) >2
        for i = 1:length(out.V_in.time)
            V_in.value(:,i) = reshape(out.V_in.signals.values(:,:,i),[],1);
        end
    else
        V_in.value = out.V_in.signals.values;
    end
    V_in.time = out.V_in.time;

    switch FLAG.SlinkModel
        case 1 % 3 state plant
            t = out.tout;
            x = x_plant;
            z = z_plant;
            u = V_in;
            w = 0;
            v = 0;
        case 2 % 3 state plant with input noise
            t = out.tout;
            x = x_plant_IN;
            z = z_plant_IN;
            u = V_in;
            w = w_k_IN;
            v = v_k_IN;
        case 3 % 3 state plant with input noise
            t = out.tout;
            x = x_plant_SN;
            z = z_plant_SN;
            u = V_in;
            w = w_k_SN;
            v = v_k_SN;
    end
%% Return to Old Working Directory
cd(oldFolder);
end