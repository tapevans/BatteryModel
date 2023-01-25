%% DAE Circuit Main
%
%
%
clear all; close all; clc;
%% Inputs
% Plant
    SIM.C_1 = 1e-6;
    SIM.C_4 = 1e-4;
%     SIM.C_4 = 1e-6;
    
    SIM.R_1 = 2;
    SIM.R_2 = 2;
    SIM.R_3 = 2;
    SIM.R_4 = 2;

% Simulation inputs
    % Modes
    % 1) Step input
    % 2) Sinusoidal
    FLAG.InputType = 1;
        SIM.V_s = 1;                 % Input amplitude
        SIM.fq  = 1/(100*(SIM.R_1*SIM.C_1));

    t_final = 6e-4;               % Simulation end time
%     t_final = 6e-6;               % Simulation end time
%     t_final = 5*max(SIM.C_1*SIM.R_1,SIM.C_4*SIM.R_4);
%     t_final = 1e-1;

% What's being measured
        % 1) Measure V1
        % 2) Measure V2
        % 3) Measure V3
        % 4) All
        FLAG.Measure = 1;

% Estimator
% Noise
    SIM.Q_0 = 1e-2;
    SIM.R_0 = 1e-4;

% Sampling Rate
    % Modes
    % 1) Known
    % 2) Samples/Cycle
    % 3) Samples/Tau
    FLAG.SampleType = 1;
        SIM.Ts_known = 1e-7;
        SIM.SpC = 10;
        SIM.SpT = 10;


%% Initialize
% Pointers
    i = 1;
    P.V_1  = i; i = i+1;
    P.V_2  = i; i = i+1;
    P.V_3  = i; i = i+1;

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
    [N.measur , N.states]  = size(C);
     N.inputs = 1;

% Initial Conditions
    switch FLAG.InputType
        case 1 % Step
            SIM.x_0 = zeros(N.states,1);
%             SIM.x_0(P.V_1) = SIM.V_s;
%             SIM.x_0(P.V_2) = (SIM.R_3/(SIM.R_2 + SIM.R_3))*SIM.V_s;
        case 2 % Sine
            SIM.x_0 = zeros(N.states,1);
            %SIM.x_0(P.V_1) = SIM.V_s;
            %SIM.x_0(P.V_2) = (SIM.R_3/(SIM.R_2 + SIM.R_3))*SIM.V_s;
    end 
    
% DAE Parameters
    Mass = eye(N.states);
    Mass(P.V_2,P.V_2) = 0;

    tspan = [0,t_final];
    Tol.Rel = 1e-4;
    Tol.Abs = 1e-7;
    options_DAE = odeset('RelTol' ,Tol.Rel,      ...
                         'AbsTol' ,Tol.Abs,      ...
                         'Mass'   ,Mass);

% ODE Parameters
    Mass_ODE = Mass([P.V_1,P.V_3],[P.V_1,P.V_3]);
    tspan = [0,t_final];
    Tol.Rel = 1e-4;
    Tol.Abs = 1e-7;
    options_ODE = odeset('RelTol' ,Tol.Rel,      ...
                         'AbsTol' ,Tol.Abs,      ...
                         'Mass'   ,Mass_ODE);

    SIM.x_0_2D = SIM.x_0([P.V_1,P.V_3],:);

% Input signal
    switch FLAG.InputType
        case 1 % Step
            SIM.omega = 0;
        case 2 % Sine
            SIM.omega = 2*pi*SIM.fq;
    end

% Sampling
    switch FLAG.SampleType
        case 1 % Known
            SIM.Ts = SIM.Ts_known;
        case 2 % Samples/Cycle
            SIM.Ts = 1/SIM.fq/SIM.SpC;
        case 3 % Samples/Tau
            tau_min = min(SIM.R_1*SIM.C_1 , SIM.R_4*SIM.C_4);
            SIM.Ts  = tau_min/SIM.SpT;
    end


%% Run DAE 
[t_DAE,x_DAE] = ode15s(@(t,SV) daefun_RedCirc(t,SV,P,SIM,FLAG),tspan,SIM.x_0,options_DAE);
z_DAE = (C*x_DAE')';

figure
hold on
plot(t_DAE,x_DAE(:,1))
plot(t_DAE,x_DAE(:,2))
plot(t_DAE,x_DAE(:,3))


%% Run ODE
[t_ODE,x_ODE] = ode45(@(t,SV) odefun_RedCirc(t,SV,P,SIM,FLAG),tspan,SIM.x_0_2D,options_ODE);
% V_2 calc
    [V_2] = getV2(x_ODE(:,1),x_ODE(:,2),SIM);
    x_ODE = [x_ODE(:,1),V_2,x_ODE(:,2)];
    z_ODE = (C*x_ODE')';

% figure
% hold on
% plot(t_ODE,x_ODE(:,1))
% plot(t_ODE,x_ODE(:,2))
% plot(t_ODE,x_ODE(:,3))

% figure
% hold on
% plot(t_ODE,z_ODE,'o')
% plot(t_DAE,z_DAE,'Linewidth',2)


%% Quick Observability Check
[A,B] = getSSMat(SIM);
Mass = zeros(P.V_3);
Mass(P.V_1 , P.V_1 ) =  SIM.C_1;
Mass(P.V_3 , P.V_3 ) =  SIM.C_4;
D = zeros(N.measur,N.inputs);

% Measure V_1
C = [1 0 0];
sys_SS = dss(A,B,C,D,Mass);
% sys_SS = ss(sys_SS,'explicit');
Ob = obsv(sys_SS);
r = rank(Ob);
disp(['Measuring V_1, the rank is ' num2str(r)])

% Measure V_2
C = [0 1 0];
sys_SS = dss(A,B,C,D,Mass);
% sys_SS = ss(sys_SS,'explicit');
Ob = obsv(sys_SS);
r = rank(Ob);
disp(['Measuring V_2, the rank is ' num2str(r)])

% Measure V_3
C = [0 0 1];
sys_SS = dss(A,B,C,D,Mass);
sys_SS = ss(sys_SS,'explicit');
Ob = obsv(sys_SS);
r = rank(Ob);
disp(['Measuring V_3, the rank is ' num2str(r)])

% If finding OBSV using just A and C, measuring V_2 is the only rank
% deficient obsv
% (F:\TylerFiles\GitHubRepos\BatteryModel\Model\FunctionalityTests\ObservabilitySensitivity\CircuitAnalysis_ROM.m)

% t_final_step = 6e-4;
% [y_SS_ROM,t_SS_ROM,x_SS_ROM] = step(sys_SS,t_final_step);
% figure
% hold on
% plot(t_SS_ROM,x_SS_ROM(:,1))
% plot(t_SS_ROM,x_SS_ROM(:,2))
% % plot(t_SS_ROM,x_SS_ROM(:,3))
% 
% x_SS_ROM_converted = ((sys_SS.C*x_SS_ROM') + [0.69;0.35;0])'; %
% 
% figure
% hold on
% plot(t_SS_ROM,x_SS_ROM_converted(:,1))
% plot(t_SS_ROM,x_SS_ROM_converted(:,2))
% plot(t_SS_ROM,x_SS_ROM_converted(:,3))


%% Find a ROM using Ho-Kalman
% Impulse Response
    FLAG.InputType = 3;
    
    % DAE Parameters
        tspan = 0:SIM.Ts/10:t_final;
        SIM.x_0_imp = zeros(N.states,1);
%         SIM.x_0_imp(P.V_1) = 1/SIM.Ts;
%         SIM.x_0_imp(P.V_2) = (SIM.R_3/(SIM.R_2 + SIM.R_3))*(1/SIM.Ts);
%         SIM.x_0_imp(P.V_1) = 1;
%         SIM.x_0_imp(P.V_2) = (SIM.R_3/(SIM.R_2 + SIM.R_3))*(1);
    
        Mass = eye(N.states);
        Mass(P.V_2,P.V_2) = 0;
        Tol.Rel = 1e-4;
        Tol.Abs = 1e-7;
        options_DAE = odeset('RelTol' ,Tol.Rel,      ...
                             'AbsTol' ,Tol.Abs,      ...
                             'Mass'   ,Mass);
            
    [t_impulse,x_impulse] = ode15s(@(t,SV) daefun_RedCirc(t,SV,P,SIM,FLAG),tspan,SIM.x_0_imp,options_DAE);
    
    figure
    hold on
    plot(t_impulse,x_impulse(:,1),'o')
    plot(t_impulse,x_impulse(:,2),'o')
    plot(t_impulse,x_impulse(:,3),'o')
    xlim([0 6e-4])

% Create Hankle Matrix
    idx = 1:10:length(t_impulse); %+10
    g_k = x_impulse(idx,:)';
    g_k = g_k(:,1:100);
    g_k = g_k - g_k(:,1);
    H = myHankle(g_k);

% Extract ROM
    % SVD
    [Nrows, Ncolms] = size(H);
    N_outputs = Nrows/Ncolms;
    N_in = N.inputs; %%% !!! Hardcoded but always true for batteries
    
    [U,S,V] = svd(H);
    r = rank(S);
    r = 5;
    
    U_colm = U(:,1:r);
    S_colm = S(1:r,1:r);
    V_colm = V(:,1:r);
    % U_null = U(:,r+1:end);
    % S_null = S(r+1:end,r+1:end);
    % V_null = V(:,r+1:end);
    
%     if localFLAG.BalRed % If using balanced reduction
        S_sqrt = S_colm.^0.5;
        obsv_r = U_colm*S_sqrt;
        cont_r = S_sqrt*V_colm';
%     else % If using normal reduction
%         obsv_r = U_colm;
%         cont_r = S_colm*V_colm';
%     end
    
    % ABC
    % C_r
    C_r = obsv_r(1:N_outputs,:   );
    
    % B_r
    B_r = cont_r(:        ,N_in);
    
    % A_r
    P   = obsv_r(1:end-N_outputs,:);
    P_p = obsv_r(N_outputs+1:end,:);
    threshold = 1e-10;
    P_inv = pinv(P,threshold);
    A_r = P_inv*P_p;
    
    % D_r ~ !!!!! I'm assuming this is correct
    D_r = zeros(N_outputs,N_in);


%% ROM Step Response
    SS_ROM = ss(A_r,B_r,C_r,D_r,SIM.Ts);
    t_final_step = 6e-4;
    [y_SS_ROM,t_SS_ROM,x_SS_ROM] = step(SS_ROM,t_final_step);
    
    x_SS_ROM_converted = ((C_r*x_SS_ROM') )'; %+ SIM.x_0

    figure
    hold on
    plot(t_SS_ROM,x_SS_ROM_converted(:,1))
    plot(t_SS_ROM,x_SS_ROM_converted(:,2))
    plot(t_SS_ROM,x_SS_ROM_converted(:,3))
%     ylim([0,1])


%% Simulink
% if FLAG.Simulink_NoNoise
%     % Covariance
%         Q = 0*eye(N.inputs);
%         R = 0*eye(N.measur);
% 
%     mdl = 'DAE_Circuit_Sim';
%     load_system(mdl)
%     in = Simulink.SimulationInput(mdl);
%     in = in.setModelParameter('StartTime','0','StopTime',num2str(t_final));
%     
%     mdlWks = get_param(in,'ModelWorkspace');
%     assignin(mdlWks,'C_1'    ,SIM.C_1)
%     assignin(mdlWks,'C_4'    ,SIM.C_4)
%     assignin(mdlWks,'R_1'    ,SIM.R_1)
%     assignin(mdlWks,'R_2'    ,SIM.R_2)
%     assignin(mdlWks,'R_3'    ,SIM.R_3)
%     assignin(mdlWks,'R_4'    ,SIM.R_4)
%     assignin(mdlWks,'Vs'     ,SIM.V_s)
%     assignin(mdlWks,'Q'      ,Q)
%     assignin(mdlWks,'R'      ,R)
%     assignin(mdlWks,'x_0_2D' ,SIM.x_0_2D)
%     assignin(mdlWks,'C'      ,C)
%     assignin(mdlWks,'Ts'     ,SIM.Ts)
% 
%     out_NN = sim(in);
% 
%     t_Slink_NN   = out_NN.tout;
%     x_Slink_NN   = out_NN.x;
%     z_Slink_NN   = out_NN.z;
%     un_Slink_NN  = out_NN.Vs_n;
%     y_Slink_NN   = out_NN.y;
%     w_k_Slink_NN = out_NN.w_k;
%     v_k_Slink_NN = out_NN.v_k;
% 
% end



%% Helper Functions
%% DAE Governing Equations
function [dSVdt] = daefun_RedCirc(t,SV,P,SIM,FLAG)
% Input Voltage
    switch FLAG.InputType
        case 1 % Step
            u = SIM.V_s;
        case 2 % Sine
            u = SIM.V_s*sin(SIM.omega*t);
        case 3 % Impulse
            if t <= SIM.Ts
%                 u = 1/SIM.Ts;
                u = 1;
            else
                u = 0;
            end
    end

% Derivative of Input
    dVsdt = SIM.V_s*SIM.omega*cos(SIM.omega*t);

% State Derivatives
dSVdt = zeros(size(SV));

dSVdt(P.V_1,1) = -( 1/(SIM.R_1*SIM.C_1) +  1/(SIM.R_2*SIM.C_1) )*SV(P.V_1)...
                 +(1/(SIM.R_2*SIM.C_1))                         *SV(P.V_2)...
                 +(1/(SIM.R_1*SIM.C_1))                         *u...
                 +dVsdt;

dSVdt(P.V_2,1) =  (-1/SIM.R_2)              *SV(P.V_1)...
                 +( 1/SIM.R_2 +  1/SIM.R_3 )*SV(P.V_2)...
                 +(-1/SIM.R_3)              *SV(P.V_3);

dSVdt(P.V_3,1) =  ( 1/(SIM.R_3*SIM.C_4))                        *SV(P.V_2)...
                 -( 1/(SIM.R_3*SIM.C_4) +  1/(SIM.R_4*SIM.C_4) )*SV(P.V_3);
end


%% ODE Governing Equations
function [dSVdt] = odefun_RedCirc(t,SV,P,SIM,FLAG)
% Input Voltage
    switch FLAG.InputType
        case 1 % Step
            u = SIM.V_s;
        case 2 % Sine
            u = SIM.V_s*sin(SIM.omega*t);
        case 3 % Impulse
            if t < SIM.Ts
                u = 1/SIM.Ts;
            else
                u = 0;
            end
    end

% Derivative of Input
    dVsdt = SIM.V_s*SIM.omega*cos(SIM.omega*t);

% State Derivatives
    dSVdt = zeros(size(SV));

% V_2 Calculation
    [V_2] = getV2(SV(P.V_1),SV(P.V_2),SIM); %V_2 is V_3 in this case

dSVdt(P.V_1,1) = -( 1/(SIM.R_1*SIM.C_1) +  1/(SIM.R_2*SIM.C_1) )*SV(P.V_1)...
                 +(1/(SIM.R_2*SIM.C_1))                         *V_2...
                 +(1/(SIM.R_1*SIM.C_1))                         *u...
                 +dVsdt;


dSVdt(P.V_2,1) =  ( 1/(SIM.R_3*SIM.C_4))                        *V_2...
                 -( 1/(SIM.R_3*SIM.C_4) +  1/(SIM.R_4*SIM.C_4) )*SV(2);
end


%% V_2 calc
function [V_2] = getV2(V_1,V_3,SIM)
Cons = (1/SIM.R_2 + 1/SIM.R_3)^(-1);

V_2 = (1/SIM.R_2)*Cons*V_1...
     +(1/SIM.R_3)*Cons*V_3;
end

%% SS Matricies
function [A,B] = getSSMat(SIM)
% This is only valid for step function
A = [-(1/SIM.R_1+1/SIM.R_2)  1/SIM.R_2             0
     -1/SIM.R_2             (1/SIM.R_2+1/SIM.R_3) -1/SIM.R_3
      0                      1/SIM.R_3            -(1/SIM.R_3 + 1/SIM.R_4)];

B = [1/SIM.R_1 ; 0 ; 0];
end







