%% Run DAE Circuit Plant
% Run the DAE Circuit model with Linear Control in the Simulink framework 
% to obtain noisy measurements for estimation.
% This will first use step data

% Run model at a short sampling rate so I can extract longer sampling
% rates from the same data.

clear all; close all; clc;
%% Parameters to run
% FLAG.InputSignalMode
%  1) From code
%  2) From function input
    FLAG.InputSignalMode = 2;

% Q Modes
%  1) Input Q
%  2) State Q
    FLAG.QMode = 1;

% State to Measure
%  1) V_1
%  2) V_2
%  3) V_3
    C_mode = 1;

tc = 1e-7;
% Q_0 = 1e-5;
% R_0 = 1e-5;
Q_0 = 0e-5;
R_0 = 0e-5;
Ts = 1e-6;

% t_final = 5e-3;
t_final = 4e-6;
t_final = 8e-6;


%% Make Save Filename
    Q_str = num2str(Q_0);
    R_str = num2str(R_0);
    Ts_str = num2str(Ts);
    C_str = num2str(C_mode);
    switch FLAG.QMode
        case 1 % Input Q
            Noise_str = 'InputQ';
        case 2 % State Q
            Noise_str = 'StateQ';
    end
    
    filename = ['Slink_Step_' Noise_str '_Q' Q_str '_R' R_str '_Ts' Ts_str '_C' C_str '.mat'];
    filepath = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanTesting\DAE_Circuit\PlantOutputData\Results';
    save_filename = [filepath filesep filename];


%% System Parameters
    SIM.C_1 = 1e-6;
    SIM.C_4 = 1e-4;
    
    SIM.R_1 = 2;
    SIM.R_2 = 2;
    SIM.R_3 = 2;
    SIM.R_4 = 2;


% Pointers
    i = 1;
    P.V_s   = i; i = i+1;
    P.deltaV= i; i = i+1;%P.V_1   = i; i = i+1;
    P.V_2   = i; i = i+1;
    P.V_3   = i; i = i+1;
    P.i_PS  = i; i = i+1;

% C Matrix
    C_des = eye(5);
    switch C_mode
        case 1
            C_m = C_des(1,:);
        case 2
            C_m = C_des(1,:);
        case 3
            C_m = C_des(1,:);
    end
    
    [N.measur , N.states]  = size(C_m);
     N.inputs = 1;

% Mass Matrix
    Mass = zeros(N.states);

    Mass(P.V_s   ,P.deltaV) =  SIM.C_1;
    Mass(P.deltaV,P.deltaV) = -SIM.C_1;
    Mass(P.V_3   ,P.V_3   ) =  SIM.C_4;

%     Mass(P.V_s,P.V_1) = -SIM.C_1;

    SIM.Mass = Mass;

    SIM.algb_idx = [P.V_2    , P.i_PS];
    SIM.diff_idx = [P.deltaV , P.V_3, P.V_s]; %P.V_1

% Initial Conditions
    SIM.V_init        = 0;
        V_step_height = 1;
    SIM.V_step        = SIM.V_init + V_step_height;

    SIM.time_rest = 8*tc+Ts;
    SIM.time_ramp = 2*tc;
    
    SIM.x_0 = zeros(N.states,1);

    SIM.x_0(P.V_s)  = SIM.V_init;
    SIM.x_0(P.i_PS) = SIM.V_init / ( SIM.R_2 + SIM.R_3 ); % SIM.R_1 +  + SIM.R_4
    SIM.x_0(P.deltaV) = (SIM.x_0(P.i_PS) * SIM.R_1 );
    V_1 = SIM.x_0(P.V_s) - SIM.x_0(P.deltaV);
    SIM.x_0(P.V_2)    = - (SIM.x_0(P.i_PS) * SIM.R_2  - V_1);
    SIM.x_0(P.V_3)    = - (SIM.x_0(P.i_PS) * SIM.R_3  - SIM.x_0(P.V_2));

    %SIM.x_0(P.V_1)  = - (SIM.x_0(P.i_PS) * SIM.R_1  - SIM.x_0(P.V_s));

% Time constants
    tau_min = SIM.R_1 * SIM.C_1;
    tau_max = SIM.R_4 * SIM.C_4;

%% Generate Input Signal
% Rest for initial 8*tc + Ts
% Ramp from (Ts + 8*tc) to (Ts + 10*tc)
    t_final_sim = t_final + 2*Ts;
    t_vec = 0:tc:t_final_sim;
    for i = 1:length(t_vec)
        t = t_vec(i);
        [V_in(i)] = getInput(t,SIM);
    end
    InputSignal = [t_vec' V_in'];


%% Test DAE
%     tspan = [0,t_final_sim];
%     Tol.Rel = 1e-4;
%     Tol.Abs = 1e-7;
%     options_DAE = odeset('RelTol' ,Tol.Rel,      ...
%                          'AbsTol' ,Tol.Abs,      ...
%                          'Mass'   ,Mass);
%     V_in = SIM.V_init;
%     FLAG.InputSignalMode = 1;
% 
%     [t_DAE,x_DAE] = ode15s(@(t,SV) DAECircGovnEqn(t,SV,SIM,N,P,FLAG,V_in),tspan,SIM.x_0,options_DAE);
%     z_DAE = (C_m*x_DAE')';
% 
%     % Calc V1
%     V_1 = x_DAE(:,P.V_s) - x_DAE(:,P.deltaV);
%     
% 
%     figure
%     hold on
% %     plot(t_DAE,x_DAE(:,P.V_1))
%     plot(t_DAE,V_1)
%     plot(t_DAE,x_DAE(:,P.V_2))
%     plot(t_DAE,x_DAE(:,P.V_3))


% %% Test A,B
%     [A,B] = getCTAandB(SIM,N,P);
%     C = eye(5,5);
%     D = zeros(5,1);
%     
%     sys = dss(A,B,C,D,Mass);
%     sys = ss(sys,'explicit');
%     [y_sys,t_sys,x_sys] = step(sys,6e-6);
% 
%     V_1 = y_sys(:,P.V_s) - y_sys(:,P.deltaV);
% 
%     figure
%     hold on
%     plot(t_sys,V_1)
%     plot(t_sys,y_sys(:,P.V_2))
%     plot(t_sys,y_sys(:,P.V_3))




%% Calculate Noise Power
    switch FLAG.QMode
        case 1 % Input Q
            Q = Q_0 * eye(N.inputs);
        case 2 % State Q
            Q = Q_0 * eye(N.states);
            Q = (diag(Q))';
    end
    R = R_0 * eye(N.measur);
    R = (diag(R))';
    
    pow_Q = tc * Q;
    pow_R = tc * R;


%% SS System
[A,B] = getCTAandB(SIM,N,P);


%% Calculate LC parameters
[LC_param] = getLCParam(SIM,N,P);
M_cross = LC_param.M_cross;
U_NV_NT = LC_param.U_NV_NT;
K_plant = LC_param.K;


%% Run Simulink
    SV_IC = SIM.x_0;

    switch FLAG.QMode
        case 1 % Input Q
            mdl = 'Circuit_5State_Overall';
        case 2 % State Q
            mdl = 'Circuit_5State_Overall_StateQ';
    end

    load_system(mdl)
    in = Simulink.SimulationInput(mdl);
    in = in.setModelParameter('StartTime','0','StopTime',num2str(t_final_sim));
    
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

if length(size(out.w_k.signals.values)) >2
    for i = 1:length(out.w_k.time)
        w_k.value(:,i) = reshape(out.w_k.signals.values(:,:,i),[],1);
    end
else
    w_k.value = out.w_k.signals.values;
end
w_k.time = out.w_k.time;

if length(size(out.v_k.signals.values)) >2
    for i = 1:length(out.v_k.time)
        v_k.value(:,i) = reshape(out.v_k.signals.values(:,:,i),[],1);
    end
else
    v_k.value = out.v_k.signals.values;
end
v_k.time = out.v_k.time;

if length(size(out.V_in.signals.values)) >2
    for i = 1:length(out.V_in.time)
        V_in.value(:,i) = reshape(out.V_in.signals.values(:,:,i),[],1);
    end
else
    V_in_sim.value = out.V_in.signals.values;
end
V_in_sim.time = out.V_in.time;



%% Save File
save(save_filename,'Q_0','R_0','C_mode','Ts','FLAG','N','P','SIM','x_plant','z_plant','w_k','v_k','V_in_sim','InputSignal','A','B')


%% Check solution plot
figure
hold on
% plot(x_plant.time,x_plant.value(P.V_1,:))
plot(x_plant.time,x_plant.value(P.V_2,:))
plot(x_plant.time,x_plant.value(P.V_3,:))

% figure
% hold on
% plot(x_plant.time,x_plant.value(P.V_s,:))
% plot(x_plant.time,x_plant.value(P.i_PS,:))

%% Run DAE check
% Test input voltage logic
%     t_vec = 0:tc:11*tc;
%     for i = 1:length(t_vec)
%         t = t_vec(i);
%         if t <= SIM.time_rest
%             V_in(i) = SIM.V_init;
%         elseif t > SIM.time_rest && t < (SIM.time_rest + SIM.time_ramp)
%             V_in(i) = (SIM.V_step - SIM.V_init)/SIM.time_ramp*(t-SIM.time_ramp) + SIM.V_init;
%         else
%             V_in(i) = SIM.V_step;
%         end
%     end
%     plot(t_vec,V_in)

% for i = 1:length(t_vec)
%     t = t_vec(i);
%     [V_in(i)] = getInput(t,SIM);
% end
% plot(t_vec,V_in)
    

%     tspan = [0,t_final_sim];
%     Tol.Rel = 1e-4;
%     Tol.Abs = 1e-7;
%     options_DAE = odeset('RelTol' ,Tol.Rel,      ...
%                          'AbsTol' ,Tol.Abs,      ...
%                          'Mass'   ,Mass);
%     V_in = SIM.V_init;
%     FLAG.InputSignalMode = 1;
% 
%     [t_DAE,x_DAE] = ode15s(@(t,SV) DAECircGovnEqn(t,SV,SIM,N,P,FLAG,V_in),tspan,SIM.x_0,options_DAE);
%     z_DAE = (C_m*x_DAE')';
% 
%     figure
%     hold on
% %     plot(t_DAE,x_DAE(:,P.V_1))
%     plot(t_DAE,x_DAE(:,P.V_2))
%     plot(t_DAE,x_DAE(:,P.V_3))

%     figure
%     hold on
%     plot(t_DAE,x_DAE(:,P.V_s))
%     plot(t_DAE,x_DAE(:,P.i_PS))

%% Test A,B
%     [A,B] = getCTAandB(SIM,N,P);
%     C = eye(3,5);
%     D = zeros(3,1);
%     
%     sys = dss(A,B,C,D,Mass);
%     sys = ss(sys,'explicit');
%     [y_sys,t_sys,x_sys] = step(sys,6e-6);
% 
%     figure
%     hold on
%     plot(t_sys,y_sys(:,1))
%     plot(t_sys,y_sys(:,2))
%     plot(t_sys,y_sys(:,3))

%% Arrange Figures
FigArrange = 1;
if FigArrange == 1
    fig = gcf;
    NumFig = fig.Number;

    Ncol = 3;
    
    for i = 1:NumFig
        f = figure(i);
        k = mod(i-1,Ncol);
        row = mod(fix((i-1)/Ncol),2);
        if row == 0
            r = 575;
%             r = 540;
        elseif row == 1
            r = 62;
        end
        f.Position = [k*575+15 r 560 420];
    end
end
