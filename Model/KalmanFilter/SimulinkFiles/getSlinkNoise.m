%% Test Getting Noise for Slink
function [InputSig_w , InputSig_v] = getSlinkNoise(SIM,N,P,FLAG)
%% Add Battery Model to path
    addpath(genpath('F:\TylerFiles\GitHubRepos\BatteryModel\Model'))


%% Change Working Directory
    [current_file_path,~,~] = fileparts(mfilename('fullpath'));
    oldFolder = cd(current_file_path); % Change to the folder that contains the Simulink files. 
    % This is mostly so the cache files don't save to the main workspace


%% Define Variables
    switch FLAG.QMode
        case 1 % Input Q
            pow_Q = SIM.pow_Qi;
        case 2 % State Q
            pow_Q = SIM.pow_Qs;
    end

 %% Run Simulink
    mdl = 'getNoise';
    load_system(mdl)
    in = Simulink.SimulationInput(mdl);
    if FLAG.Tswitch == 100 && FLAG.UseROMAsPlant
        in = in.setModelParameter('StartTime','0','StopTime',num2str( SIM.InputSignal(1001,1) ) );
    else
        in = in.setModelParameter('StartTime','0','StopTime',num2str( SIM.InputSignal(end,1) ) );
    end
    
    
    % Assign variable values by modifying the workspace
    mdlWks = get_param(in,'ModelWorkspace');   


    assignin(mdlWks,'pow_Qi'      ,pow_Q)
    assignin(mdlWks,'pow_R'       ,SIM.pow_R)
    assignin(mdlWks,'Ts'          ,SIM.Ts)
    assignin(mdlWks,'tc'          ,SIM.tc)


%     tic
    out = sim(in);
%     toc
    save_system
    close_system


%%
    %%%%%%%%%%%% Process Noise Applied to Input Noise Plant %%%%%%%%%%%% 
    if length(size(out.w_k_IN.signals.values)) >2
        for i = 1:length(out.w_k_IN.time)
            w_k_IN.value(:,i) = reshape(out.w_k_IN.signals.values(:,:,i),[],1);
        end
    else
        w_k_IN.value = out.w_k_IN.signals.values;
    end
    w_k_IN.time = out.w_k_IN.time;
%     RESULTS.Slink_plant.w.t_soln = w_k_IN.time;
%     RESULTS.Slink_plant.w.w_soln = w_k_IN.value;
    
    %%%%%%%%%%%% Measurement Noise Applied to Input Noise Plant %%%%%%%%%%%% 
    if length(size(out.v_k_IN.signals.values)) >2
        for i = 1:length(out.v_k_IN.time)
            v_k_IN.value(:,i) = reshape(out.v_k_IN.signals.values(:,:,i),[],1);
        end
    else
        v_k_IN.value = out.v_k_IN.signals.values;
    end
    v_k_IN.time = out.v_k_IN.time;
%     RESULTS.Slink_plant.v.t_soln = v_k_IN.time;
%     RESULTS.Slink_plant.v.v_soln = v_k_IN.value;


%% Add Zero Initial 
    zero_off_time_vec = [0; v_k_IN.time+SIM.Ts];
    zero_off_w_vec    = [0; w_k_IN.value'];
    zero_off_v_vec    = [0; v_k_IN.value'];


%% Make New Vector
if FLAG.Tswitch == 100 && FLAG.UseROMAsPlant
    for i = 1:10
        if i == 1
            t_vec = [ zero_off_time_vec ];
            w_vec = [ zero_off_w_vec ];
            v_vec = [ zero_off_v_vec ];
        else
            t_vec = [t_vec; v_k_IN.time+t_vec(end)+SIM.Ts];
            w_vec = [w_vec; w_k_IN.value'];
            v_vec = [v_vec; v_k_IN.value'];
        end
    end
    InputSig_w = [t_vec , w_vec];
    InputSig_v = [t_vec , v_vec];
else
    t_vec_interp = 0;
    w_vec_interp = 0;
    v_vec_interp = 0;
    
    for i = 2:length(zero_off_time_vec)
        % Minus
        t_vec_interp(end+1,1) = zero_off_time_vec(i) - SIM.Ts/20;
        w_vec_interp(end+1,1) = zero_off_w_vec(i-1);
        v_vec_interp(end+1,1) = zero_off_v_vec(i-1);
    
        % k
        t_vec_interp(end+1,1) = zero_off_time_vec(i);
        w_vec_interp(end+1,1) = ( zero_off_w_vec(i-1) + zero_off_w_vec(i) )/2;
        v_vec_interp(end+1,1) = ( zero_off_v_vec(i-1) + zero_off_v_vec(i) )/2;
    
        % Plus
        t_vec_interp(end+1,1) = zero_off_time_vec(i) + SIM.Ts/20;
        w_vec_interp(end+1,1) = zero_off_w_vec(i);
        v_vec_interp(end+1,1) = zero_off_v_vec(i);
    end
    
    InputSig_w = [t_vec_interp , w_vec_interp];
    InputSig_v = [t_vec_interp , v_vec_interp];
end

%% Plot
% % Process Noise
%     figure
%     hold on
%     plot(w_k_IN.time      ,w_k_IN.value  ,'o','DisplayName' ,'Slink')
%     plot(zero_off_time_vec,zero_off_w_vec,'ro','DisplayName','Zero Offset')
%     plot(t_vec_interp     ,w_vec_interp  ,'k','DisplayName' ,'Interpolated')
%     title('Process Noise')
%     lgn = legend;
%     lgn.Location = 'best';
% 
% % Measurement Noise
%     figure
%     hold on
%     plot(v_k_IN.time      ,v_k_IN.value  ,'o','DisplayName' ,'Slink')
%     plot(zero_off_time_vec,zero_off_v_vec,'ro','DisplayName','Zero Offset')
%     plot(t_vec_interp     ,v_vec_interp  ,'k','DisplayName' ,'Interpolated')
%     title('Measurement Noise')
%     lgn = legend;
%     lgn.Location = 'best';


    
%% Return to Old Working Directory
    cd(oldFolder);


end