%% DAE Circuit Main
%
%
%
clear all; close all; clc;
%% Inputs
% Plant
%     SIM.C_1 = 1e-6;
%     SIM.C_4 = 1e-4;
% %     SIM.C_4 = 1e-6;
%     
%     SIM.R_1 = 2;
%     SIM.R_2 = 2;
%     SIM.R_3 = 2;
%     SIM.R_4 = 2;

    SIM.C_1 = 1e1;
    SIM.C_4 = 1e3;
    
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
        SIM.fq  = 1/((SIM.R_1*SIM.C_1));
        SIM.ramp_time = 1e-2;
%         SIM.fq  = 1/20;
        SIM.Npoints_step = 60;
        SIM.Npoints_sine = 1000;

%     t_final = 6e-4;               % Simulation end time
%     t_final = 10e-6;               % Simulation end time
%     t_final = 5*max(SIM.C_1*SIM.R_1,SIM.C_4*SIM.R_4);
    t_final = 5e1;
%     t_final = 100;
%     t_final = 10000;

% What's being measured
        % 1) Measure V1
        % 2) Measure V2
        % 3) Measure V3
        % 4) All
        FLAG.Measure = 3;

% Estimator
% Noise
    SIM.Q_0 = 1e-1;
    SIM.R_0 = 1e-3;

% Sampling Rate
    % Modes
    % 1) Known
    % 2) Samples/Cycle
    % 3) Samples/Tau
    FLAG.SampleType = 1;
%         SIM.Ts_known = 1e-7;
        SIM.Ts_known = 1e0;
% %         SIM.Ts_known = 1e1;
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
SIM.x_0
[t_DAE,x_DAE] = ode15s(@(t,SV) daefun_RedCirc(t,SV,P,SIM,FLAG),tspan,SIM.x_0,options_DAE);
z_DAE = (C*x_DAE')';

figure
hold on
plot(t_DAE,x_DAE(:,1))
plot(t_DAE,x_DAE(:,2))
plot(t_DAE,x_DAE(:,3))
title('DAE Results')


%% Run ODE
SIM.x_0_2D
[t_ODE,x_ODE] = ode45(@(t,SV) odefun_RedCirc(t,SV,P,SIM,FLAG),tspan,SIM.x_0_2D,options_ODE);
% V_2 calc
    [V_2] = getV2(x_ODE(:,1),x_ODE(:,2),SIM);
    x_ODE = [x_ODE(:,1),V_2,x_ODE(:,2)];
    z_ODE = (C*x_ODE')';

figure
hold on
plot(t_ODE,x_ODE(:,1))
plot(t_ODE,x_ODE(:,2))
plot(t_ODE,x_ODE(:,3))
title('ODE Results')

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
    old_input_type = FLAG.InputType;
    FLAG.InputType = 3;
    
    % DAE Parameters
        tspan = 0:SIM.Ts/10:t_final;
        SIM.x_0_imp = zeros(N.states,1);
%         SIM.x_0_imp(P.V_1) = 1/SIM.Ts;
%         SIM.x_0_imp(P.V_2) = (SIM.R_3/(SIM.R_2 + SIM.R_3))*(1/SIM.Ts);
%         SIM.x_0_imp(P.V_1) = SIM.V_s;
%         SIM.x_0_imp(P.V_2) = (SIM.R_3/(SIM.R_2 + SIM.R_3))*(SIM.V_s);
    
        Mass = eye(N.states);
        Mass(P.V_2,P.V_2) = 0;
        Tol.Rel = 1e-4;
        Tol.Abs = 1e-7;
        options_DAE = odeset('RelTol' ,Tol.Rel,      ...
                             'AbsTol' ,Tol.Abs,      ...
                             'Mass'   ,Mass);

    SIM.x_0_imp
    [t_impulse,x_impulse] = ode15s(@(t,SV) daefun_RedCirc(t,SV,P,SIM,FLAG),tspan,SIM.x_0_imp,options_DAE);
    
    figure
    hold on
    plot(t_impulse,x_impulse(:,1),'o')
    plot(t_impulse,x_impulse(:,2),'o')
    plot(t_impulse,x_impulse(:,3),'o')
%     xlim([0 6e-4])
    title('Impulse Response')

% Create Hankle Matrix
    idx = 1:10:length(t_impulse); %+10
    g_k = x_impulse(idx,:)';
%     g_k = g_k(:,1:50);
    g_k = g_k - g_k(:,1);
    H = myHankle(g_k);

% Extract ROM
    % SVD
    [Nrows, Ncolms] = size(H);
    N_outputs = Nrows/Ncolms;
    N_in = N.inputs; %%% !!! Hardcoded but always true for batteries
    
    [U,S,V] = svd(H);
    r = rank(S);
    r = 4;
    
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

    switch FLAG.Measure
        case 1
            C_r_meas = C_r(1,:);
        case 2
            C_r_meas = C_r(2,:);
        case 3
            C_r_meas = C_r(3,:);
        case 4
            C_r_meas = C_r;
    end
    N.states_ROM = r;
    x0_ROM = zeros(length(A_r),1);

%% ROM Step Response
    SS_ROM = ss(A_r,B_r,C_r,D_r,SIM.Ts);
%     t_final_step = 6e-4;
    t_final_step = t_final;
    [y_SS_ROM,t_SS_ROM,x_SS_ROM] = step(SS_ROM,t_final_step);
    
    x_SS_ROM_converted = ((C_r*x_SS_ROM')+ SIM.x_0 )'; %

    figure
    hold on
    plot(t_SS_ROM,x_SS_ROM_converted(:,1))
    plot(t_SS_ROM,x_SS_ROM_converted(:,2))
    plot(t_SS_ROM,x_SS_ROM_converted(:,3))
    title('ROM Step Response')
    %     ylim([0,1])

FLAG.InputType = old_input_type;

%% Simulink
% if FLAG.Simulink_NoNoise
%     t_final = 100;

    % Input Signal
    switch FLAG.InputType
    case 1 % Step
        t_input = (linspace(0,t_final,SIM.Npoints_step))';
        signal  = ones(length(t_input),1);
        InputSignal = [t_input , signal];
    case 2 % Sine
        t_input = (linspace(0,t_final,SIM.Npoints_sine))';
        signal  = SIM.V_s*sin(SIM.omega*t_input);
        InputSignal = [t_input , signal];
    end 

    % Covariance
    Q = SIM.Q_0*eye(N.inputs); % Process Noise
    R = SIM.R_0*eye(N.measur); % Measurement Noise
    confidence = 0.5;
    P_k_0 = confidence*eye(N.states_ROM);
    Iden = eye(N.states_ROM);

%     x0_est = SIM.x_0;
    x0_est = zeros(length(A_r),1);
%     x0_est = C_r\SIM.x_0;

    [P_infty  ,~,~] =  dare(A_r', C_r_meas', B_r*Q*B_r' ,R);
    K_infty = P_infty*C_r_meas'*inv(C_r_meas*P_infty*C_r_meas' + R);

    % Model Parameters
    mdl = 'Circuit_Overall';
    load_system(mdl)
    in = Simulink.SimulationInput(mdl);
    in = in.setModelParameter('StartTime','0','StopTime',num2str(t_final));
    
    mdlWks = get_param(in,'ModelWorkspace');
    assignin(mdlWks,'C1'    ,SIM.C_1)
    assignin(mdlWks,'C4'    ,SIM.C_4)
    assignin(mdlWks,'R1'    ,SIM.R_1)
    assignin(mdlWks,'R2'    ,SIM.R_2)
    assignin(mdlWks,'R3'    ,SIM.R_3)
    assignin(mdlWks,'R4'    ,SIM.R_4)
    assignin(mdlWks,'Vs'    ,SIM.V_s)
    assignin(mdlWks,'F_in'  ,SIM.V_s)
    assignin(mdlWks,'Q'     ,Q)
    assignin(mdlWks,'R'     ,R)
    assignin(mdlWks,'x0_2D' ,SIM.x_0_2D)
    assignin(mdlWks,'x0_ROM',x0_ROM)
    assignin(mdlWks,'C'     ,C)
    assignin(mdlWks,'Ts'    ,SIM.Ts)
    assignin(mdlWks,'omega' ,SIM.omega)
    assignin(mdlWks,'A_ROM' ,A_r)
    assignin(mdlWks,'B_ROM' ,B_r)
    assignin(mdlWks,'C_ROM' ,C_r_meas)
    assignin(mdlWks,'InputSignal' ,InputSignal)
    assignin(mdlWks,'A_est' ,A_r)
    assignin(mdlWks,'B_est' ,B_r)
    assignin(mdlWks,'C_est' ,C_r_meas)
    assignin(mdlWks,'x0_est',x0_est)
    assignin(mdlWks,'P_k_0',P_k_0)
    assignin(mdlWks,'Iden' ,Iden)
    assignin(mdlWks,'K_infty',K_infty)


    out = sim(in);
    save_system
    close_system
%%

    t_ODE_Slink  = out.x_ODE.time;
    x_ODE_Slink  = out.x_ODE.signals.values;

    t_ODE_N_Slink  = out.x_ODE_N.time;
    x_ODE_N_Slink  = out.x_ODE_N.signals.values;

    t_ROM_Slink  = out.x_ROM.time;
    x_ROM_Slink  = out.x_ROM.signals.values;
    x_ROM_Slink_converted = (C_r*x_ROM_Slink')'; % Only works because C = eye. Would need a transformation matrix otherwise

    t_EST_Slink  = out.x_hat.time;
    if length(size(out.x_hat.signals.values))>2
        for i = 1:length(out.x_hat.time)
            x_EST_Slink(i,:) = (out.x_hat.signals.values(:,:,i))';
        end
    else
        x_EST_Slink  = out.x_hat.signals.values;
    end
    x_EST_Slink_converted = (C_r*x_EST_Slink')'; % Only works because C = eye. Would need a transformation matrix otherwise

    t_EST_asy_Slink  = out.x_hat_asy.time;
    if length(size(out.x_hat_asy.signals.values))>2
        for i = 1:length(out.x_hat_asy.time)
            x_EST_asy_Slink(i,:) = (out.x_hat_asy.signals.values(:,:,i))';
        end
    else
        x_EST_asy_Slink  = out.x_hat_asy.signals.values;
    end
    x_EST_asy_Slink_converted = (C_r*x_EST_asy_Slink')'; % Only works because C = eye. Would need a transformation matrix otherwise

figure
hold on
% plot(t_DAE,x_DAE(:,1))
% plot(t_ODE_Slink,x_ODE_Slink(:,1),'o','Linewidth',3)
plot(t_ODE_N_Slink,x_ODE_N_Slink(:,1),'ko')
plot(t_ROM_Slink,x_ROM_Slink_converted(:,1),'go')
% plot(t_EST_Slink,x_EST_Slink_converted(:,1),'ro','Linewidth',2)
plot(t_EST_asy_Slink,x_EST_asy_Slink_converted(:,1),'bo')
title('Comparison Voltage 1')

figure
hold on
% plot(t_DAE,x_DAE(:,2))
% plot(t_ODE_Slink,x_ODE_Slink(:,2),'o','Linewidth',3)
plot(t_ODE_N_Slink,x_ODE_N_Slink(:,2),'ko')
plot(t_ROM_Slink,x_ROM_Slink_converted(:,2),'go')
% plot(t_EST_Slink,x_EST_Slink_converted(:,2),'ro','Linewidth',2)
plot(t_EST_asy_Slink,x_EST_asy_Slink_converted(:,2),'bo')
title('Comparison Voltage 2')

figure
hold on
% plot(t_DAE,x_DAE(:,3))
% plot(t_ODE_Slink,x_ODE_Slink(:,3),'o','Linewidth',3)
plot(t_ODE_N_Slink,x_ODE_N_Slink(:,3),'ko')
plot(t_ROM_Slink,x_ROM_Slink_converted(:,3),'go')
% plot(t_EST_Slink,x_EST_Slink_converted(:,3),'ro','Linewidth',2)
plot(t_EST_asy_Slink,x_EST_asy_Slink_converted(:,3),'bo')
title('Comparison Voltage 3')

% end

%% Arrange Figures
FigArrange = 1;
fig = gcf;
NumFig = fig.Number;
if FigArrange == 1
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

%% Helper Functions
%% DAE Governing Equations
function [dSVdt] = daefun_RedCirc(t,SV,P,SIM,FLAG)
% Input Voltage
    switch FLAG.InputType
        case 1 % Step
            if t < SIM.ramp_time
                u = SIM.V_s/SIM.ramp_time*t;
            else
                u = SIM.V_s;
            end
        case 2 % Sine
            u = SIM.V_s*sin(SIM.omega*t);
        case 3 % Impulse
            if t < SIM.ramp_time
%                 u = SIM.V_s/SIM.ramp_time*t;
                u = (1/SIM.Ts)/SIM.ramp_time*t;
            elseif t >=SIM.ramp_time && t < SIM.Ts
                u = 1/SIM.Ts;
            else
                u = 0;
            end
    end

% Derivative of Input
    switch FLAG.InputType
        case 1 % Step
            if t < SIM.ramp_time
                dVsdt = SIM.V_s/SIM.ramp_time;
            else
                dVsdt = 0;
            end
        case 2 % Sine
            dVsdt = SIM.V_s*SIM.omega*cos(SIM.omega*t);
        case 3 % Impulse
            if t < SIM.ramp_time
%                 dVsdt = SIM.V_s/SIM.ramp_time;
                dVsdt = (1/SIM.Ts)/SIM.ramp_time;
            else
                dVsdt = 0;
            end
    end


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
            if t < SIM.ramp_time
                u = SIM.V_s/SIM.ramp_time*t;
            else
                u = SIM.V_s;
            end
        case 2 % Sine
            u = SIM.V_s*sin(SIM.omega*t);
        case 3 % Impulse
            if t < SIM.ramp_time
%                 u = SIM.V_s/SIM.ramp_time*t;
                u = (1/SIM.Ts)/SIM.ramp_time*t;
            elseif t >=SIM.ramp_time && t < SIM.Ts
                u = 1/SIM.Ts;
            else
                u = 0;
            end
    end

% Derivative of Input
    switch FLAG.InputType
        case 1 % Step
            if t < SIM.ramp_time
                dVsdt = SIM.V_s/SIM.ramp_time;
            else
                dVsdt = 0;
            end
        case 2 % Sine
            dVsdt = SIM.V_s*SIM.omega*cos(SIM.omega*t);
        case 3 % Impulse
            if t < SIM.ramp_time
%                 dVsdt = SIM.V_s/SIM.ramp_time;
                dVsdt = (1/SIM.Ts)/SIM.ramp_time;
            else
                dVsdt = 0;
            end
    end
    

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
