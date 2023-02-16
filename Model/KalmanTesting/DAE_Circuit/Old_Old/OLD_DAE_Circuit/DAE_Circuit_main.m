%% DAE Circuit
%
%
%
clear all; close all; clc;
%% Parameters
% Filename
    filename = 'Test.mat';

% System
    SIM.C_1 = 1e-6;
    SIM.C_4 = 1e-4;
%     SIM.C_4 = 1e-6;
    
    SIM.R_1 = 2;
    SIM.R_2 = 2;
    SIM.R_3 = 2;
    SIM.R_4 = 2;
    
% Noise
    SIM.Q_0 = 1e-2;
    SIM.R_0 = 1e-4;

% Pointers
    i = 1;
    P.V_1  = i; i = i+1;
    P.V_2  = i; i = i+1;
    P.V_3  = i; i = i+1;

    
%% Simulations
% Which simulations to run
    FLAG.DAE              = 1;
    FLAG.Simulink_NoNoise = 1;   
    FLAG.Simulink_Noise   = 1;   
    FLAG.Estimator        = 0;
        FLAG.RunAsymtotic = 0;
        FLAG.RunVariable  = 0;

    FLAG.Simulink_Noise_LC= 0;   
    FLAG.Estimator_LC     = 0;
        FLAG.RunAsymtotic_LC= 0;
        FLAG.RunVariable_LC = 0;

% Simulation inputs
    % Modes
    % 1) Step input
    % 2) Sinusoidal
    FLAG.InputType = 1;
        SIM.V_s = 1;                 % Input amplitude
        SIM.fq  = 1;
    t_final = 6e-4;               % Simulation end time
%     t_final = 6e-6;               % Simulation end time
    SIM.Ts  = SIM.C_1*SIM.R_1/10; % Sampling Rate
    
    % What's being measured
        % 1) Measure V1
        % 2) Measure V2
        % 3) Measure V3
        % 4) All
        FLAG.Measure = 1;


%% Plots
    FLAG.PLOT.V1Compare = 1;
    FLAG.PLOT.V2Compare = 1;
    FLAG.PLOT.V3Compare = 1;
    
    FLAG.PLOT.DAE              = 1;
    FLAG.PLOT.Simulink_NoNoise = 1;
    FLAG.PLOT.Simulink_Noise   = 1;
    FLAG.PLOT.Estimator        = 0;
        FLAG.PLOT.RunAsymtotic = 0;
        FLAG.PLOT.RunVariable  = 0;

    FLAG.PLOT.Simulink_Noise_LC= 0;
    FLAG.PLOT.Estimator_LC     = 0;
        FLAG.PLOT.RunAsymtotic_LC = 0;
        FLAG.PLOT.RunVariable_LC  = 0;


%% Initialization
    % Mass
        Mass = zeros(P.V_3);
        
        Mass(P.V_1 , P.V_1 ) =  SIM.C_1;
        Mass(P.V_3 , P.V_3 ) =  SIM.C_4;

    % Time Vector    
        N_steps = ceil((t_final/SIM.Ts)+1);

    % Input signal
    
    % Measured Matrix
        switch FLAG.Measure
            case 1
                C = [1 0 0];
            case 2
                C = [0 1 0];
            case 3
                C = [0 0 1];
            case 4
                C = eye(3);
        end

    % Parameters
        [N.measur , N.states]  = size(C);
        N.inputs = 1;

    % Covariance
        Q = SIM.Q_0*eye(N.inputs);
        R = SIM.R_0*eye(N.measur);
    
    % Initial Conditions
        SIM.x_0 = zeros(N.states,1);
    
        SIM.x_0(P.V_1) = SIM.V_s;
        SIM.x_0(P.V_2) = (SIM.R_3/(SIM.R_2 + SIM.R_3))*SIM.V_s;

        SIM.x_0_2D = SIM.x_0([P.V_1,P.V_3],:);


%% DAE Simulation
if FLAG.DAE
    tspan = [0,t_final ];
    Tol.Rel = 1e-4;
    Tol.Abs = 1e-7;
    options = odeset('RelTol' ,Tol.Rel,      ...
                     'AbsTol' ,Tol.Abs,      ...
                     'Mass'   ,Mass);
    
    [t_DAE,x_DAE] = ode15s(@(t,SV) odefun_RedCirc(t,SV,P,SIM),tspan,SIM.x_0,options);
    z_DAE = (C*x_DAE')';
end


%% Simulink No Noise
if FLAG.Simulink_NoNoise
    % Covariance
        Q = 0*eye(N.inputs);
        R = 0*eye(N.measur);

    mdl = 'DAE_Circuit_Sim';
    load_system(mdl)
    in = Simulink.SimulationInput(mdl);
    in = in.setModelParameter('StartTime','0','StopTime',num2str(t_final));
    
    mdlWks = get_param(in,'ModelWorkspace');
    assignin(mdlWks,'C_1'    ,SIM.C_1)
    assignin(mdlWks,'C_4'    ,SIM.C_4)
    assignin(mdlWks,'R_1'    ,SIM.R_1)
    assignin(mdlWks,'R_2'    ,SIM.R_2)
    assignin(mdlWks,'R_3'    ,SIM.R_3)
    assignin(mdlWks,'R_4'    ,SIM.R_4)
    assignin(mdlWks,'Vs'     ,SIM.V_s)
    assignin(mdlWks,'Q'      ,Q)
    assignin(mdlWks,'R'      ,R)
    assignin(mdlWks,'x_0_2D' ,SIM.x_0_2D)
    assignin(mdlWks,'C'      ,C)
    assignin(mdlWks,'Ts'     ,SIM.Ts)

    out_NN = sim(in);

    t_Slink_NN   = out_NN.tout;
    x_Slink_NN   = out_NN.x;
    z_Slink_NN   = out_NN.z;
    un_Slink_NN  = out_NN.Vs_n;
    y_Slink_NN   = out_NN.y;
    w_k_Slink_NN = out_NN.w_k;
    v_k_Slink_NN = out_NN.v_k;

end


%% Simulink with Noise
if FLAG.Simulink_Noise
    % Covariance
        Q = SIM.Q_0*eye(N.inputs);
        R = SIM.R_0*eye(N.measur);

    mdl = 'DAE_Circuit_Sim';
    load_system(mdl)
    in = Simulink.SimulationInput(mdl);
    in = in.setModelParameter('StartTime','0','StopTime',num2str(t_final));
    
    mdlWks = get_param(in,'ModelWorkspace');
    assignin(mdlWks,'C_1'    ,SIM.C_1)
    assignin(mdlWks,'C_4'    ,SIM.C_4)
    assignin(mdlWks,'R_1'    ,SIM.R_1)
    assignin(mdlWks,'R_2'    ,SIM.R_2)
    assignin(mdlWks,'R_3'    ,SIM.R_3)
    assignin(mdlWks,'R_4'    ,SIM.R_4)
    assignin(mdlWks,'Vs'     ,SIM.V_s)
    assignin(mdlWks,'Q'      ,Q)
    assignin(mdlWks,'R'      ,R)
    assignin(mdlWks,'x_0_2D' ,SIM.x_0_2D)
    assignin(mdlWks,'C'      ,C)
    assignin(mdlWks,'Ts'     ,SIM.Ts)

    out_wN = sim(in);

    t_Slink_wN   = out_wN.tout;
    un_Slink_wN  = out_wN.Vs_n;
    x_Slink_wN   = out_wN.x;
    y_Slink_wN   = out_wN.y;
    z_Slink_wN   = out_wN.z;
    w_k_Slink_wN = out_wN.w_k;
    v_k_Slink_wN = out_wN.v_k;

end


%% Estimator
if FLAG.Estimator
    % Initialize Solution Variables
%         x_EST_sys = zeros(N_states,N_steps);   % Real System (Plant) States
%         x_asy     = zeros(N_states,N_steps);   % Estimator States
%         x_asy_pre = zeros(N_states,N_steps);   % Estimator States predict phase
%         x_var     = zeros(N_states,N_steps);   % Estimator States
%         x_var_pre = zeros(N_states,N_steps);   % Estimator States predict phase
%         z_k       = zeros(N_meas,N_steps);   % Measurement from Real System
%         y_tilde_k = zeros(N_meas,N_steps);   % Estimation Error
%         S_k       = zeros(N_meas,N_meas,N_steps); % Error Covariance
%         K_k       = zeros(N_states,N_meas,N_steps); % Kalman Gains
%         P_k       = zeros(N_states,N_states,N_steps); % Error Covariance
%         P_k_pre   = zeros(N_states,N_states,N_steps); % Error Covariance predict phase
end

%% Save Results
    save(filename)


%% Plot Results
    DAE_Circuit_PlotFnc(filename)


%% Helper Functions
%% DAE Governing Equations
function [dSVdt] = odefun_RedCirc(t,SV,P,SIM)
dSVdt = zeros(size(SV));

dSVdt(P.V_1,1) = -(SV(P.V_1)-SIM.V_s  )/SIM.R_1...
                 -(SV(P.V_1)-SV(P.V_2))/SIM.R_2;

dSVdt(P.V_2,1) =  (SV(P.V_2)-SV(P.V_1))/SIM.R_2...
                 +(SV(P.V_2)-SV(P.V_3))/SIM.R_3;

dSVdt(P.V_3,1) = -(SV(P.V_3)-SV(P.V_2))/SIM.R_3...
                 -(SV(P.V_3)-0        )/SIM.R_4;
end
