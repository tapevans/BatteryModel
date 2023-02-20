%% Test Getting Noise for Slink
%     Q_0_vec = [1e-6];
%     R_0_vec = [1e-6];
%     Ts_vec  = [1];
%     SOC_vec = [50];
% 
%     FLAG.Q_0 = QQ; % Process Noise
%     FLAG.R_0 = RR; % Measurement Noise
%     FLAG.Ts = TT;  % Sampling Rate (When to save outputs)
%     FLAG.SOC = SS; % State of Charge
% 
%     SIM.Q_0 = FLAG.Q_0;
%     SIM.R_0 = FLAG.R_0;
%     SIM.Ts  = FLAG.Ts;
%     SIM.tc  = FLAG.Ts;

%% Change Working Directory !!!!!!!!!!!(Remove after testing)
    [current_file_path,~,~] = fileparts(mfilename('fullpath'));
    oldFolder = cd(current_file_path); % Change to the folder that contains the Simulink files. 
    % This is mostly so the cache files don't save to the main workspace

%% Define Variables !!!!!!!!!!!(Remove after testing)
    SIM.Q_0 = [1e-6];
    SIM.R_0 = [1e-6];
    SIM.Ts  = [1];
    SIM.SOC = [50];

    FLAG.N_samples = 600;
    SIM.t_final_Relax_Step = (FLAG.N_samples+2)*SIM.Ts;
    SIM.t_vec_Relax_Step = 0:SIM.Ts:SIM.t_final_Relax_Step;

    N.inputs = 1;
    N.measur = 1;
    N.states = 355;

    SIM.tc = SIM.Ts;

    SIM.Qi = SIM.Q_0 * eye(N.inputs);
    
    Q       = SIM.Q_0 * eye(N.states);
    SIM.Qs = (diag(Q));
    
    R     = SIM.R_0 * eye(N.measur);
    SIM.R = (diag(R));
    
    SIM.pow_Qi = SIM.tc * SIM.Qi;
    SIM.pow_Qs = SIM.tc * SIM.Qs;
    SIM.pow_R  = SIM.tc * SIM.R;

    pow_Q = SIM.pow_Qi;

 %% Run Simulink
    mdl = 'getNoise';
    load_system(mdl)
    in = Simulink.SimulationInput(mdl);
    in = in.setModelParameter('StartTime','0','StopTime',num2str(SIM.t_final_Relax_Step));
    
    % Assign variable values by modifying the workspace
    mdlWks = get_param(in,'ModelWorkspace');   


    assignin(mdlWks,'pow_Qi'      ,pow_Q)
    assignin(mdlWks,'pow_R'       ,SIM.pow_R)
    assignin(mdlWks,'Ts'          ,SIM.Ts)
    assignin(mdlWks,'tc'          ,SIM.tc)


    tic
    out = sim(in);
    toc

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


%% Add Zero Initial 
    zero_off_time_vec = [0; v_k_IN.time+SIM.Ts];
    zero_off_w_vec    = [0; w_k_IN.value'];
    zero_off_v_vec    = [0; v_k_IN.value'];

%% Make New Vector
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

%% Plot
% Process Noise
    figure
    hold on
    plot(w_k_IN.time      ,w_k_IN.value  ,'o','DisplayName' ,'Slink')
    plot(zero_off_time_vec,zero_off_w_vec,'ro','DisplayName','Zero Offset')
    plot(t_vec_interp     ,w_vec_interp  ,'k','DisplayName' ,'Interpolated')
    title('Process Noise')
    lgn = legend;
    lgn.Location = 'best';

% Measurement Noise
    figure
    hold on
    plot(v_k_IN.time      ,v_k_IN.value  ,'o','DisplayName' ,'Slink')
    plot(zero_off_time_vec,zero_off_v_vec,'ro','DisplayName','Zero Offset')
    plot(t_vec_interp     ,v_vec_interp  ,'k','DisplayName' ,'Interpolated')
    title('Measurement Noise')
    lgn = legend;
    lgn.Location = 'best';
