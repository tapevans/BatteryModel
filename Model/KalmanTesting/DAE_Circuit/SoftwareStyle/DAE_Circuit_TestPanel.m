%% DAE Circuit Quick Run Function
clear all; close all; clc;

%% Test Inputs
%% State to Measure
%  1) V_1
%  2) V_2
%  3) V_3
%  4) All ~ Use all for comparing no noise simulations mostly
FLAG.C_mode = 4;

%% Sampling Rate
SIM.Ts = 1e-6; % Sampling Rate (When to save outputs)

%% Noise
% Q Modes
%  1) Input Q
%  2) State Q
FLAG.QMode = 1;

% Noise Properties
SIM.tc  = SIM.Ts; % [s],Correlation Time (Noise Sampling Rate)
SIM.Q_0 = 1e-6; % Process Noise
SIM.R_0 = 1e-6; % Measurement Noise


%% Input Type
% Inputs
%  1) Step
%  2) Harmonic
%  3) Ramp ~ subset of step
FLAG.InputMode = 1;
SIM.V_init        = 0;
SIM.V_step_height = 1;
% SIM.V_init        = 4;
% SIM.V_step_height = 4;
%         SIM.fq = SIM.Ts/15; %% WAAAAAAAY too long
SIM.fq = 400;
SIM.V_s = 1; %Amplitude of sine wave


%% Simulation t_final
%     SIM.t_final = 6e-6;
    SIM.t_final = 1e-4;
%     SIM.t_final = 6e-4;
% SIM.t_final = 5e-3;
% SIM.t_final = 5*max(SIM.C_1*SIM.R_1,SIM.C_4*SIM.R_4);

%     N_cycles = 10; %% WAAAAAAAY too long
%     SIM.t_final = N_cycles/SIM.fq;

if FLAG.InputMode == 2
    SIM.t_final = 0.005; % (5000 * 1e-6)
end




%%
[SIM,FLAG] = inputs(SIM,FLAG);
[SIM,N,P,FLAG,RESULTS] = init(SIM,FLAG);

%% Test Run Step ode 3
if true
    FLAG.StateMode = 3;
%         FLAG.InputMode = 3; % 1) Step 2) Sine 3) Ramp
    [t_soln3, x_soln3, z_soln3] = runODE(SIM,N,P,FLAG);

    %     figure
    %     hold on
    %     plot(t_soln3,x_soln3(:,P.V_1),'LineWidth',2,'DisplayName','V_1')
    %     plot(t_soln3,x_soln3(:,P.V_2),'LineWidth',2,'DisplayName','V_2')
    %     plot(t_soln3,x_soln3(:,P.V_3),'LineWidth',2,'DisplayName','V_3')
end


%% Test Run Step ode 5
if true
    FLAG.StateMode = 5;
%         FLAG.InputMode = 3; % 1) Step 2) Sine, 3) Ramp
    [t_soln5, x_soln5, z_soln5] = runODE(SIM,N,P,FLAG);

    %     figure
    %     hold on
    %     plot(t_soln5,x_soln5(:,P.V_1),'LineWidth',2,'DisplayName','V_1')
    %     plot(t_soln5,x_soln5(:,P.V_2),'LineWidth',2,'DisplayName','V_2')
    %     plot(t_soln5,x_soln5(:,P.V_3),'LineWidth',2,'DisplayName','V_3')
    %
    % % Compare ode5 to ode3
    %     figure
    %     hold on
    %     plot(t_soln5,x_soln5(:,P.V_1),'o','LineWidth',2,'DisplayName','V_1')
    %     plot(t_soln5,x_soln5(:,P.V_2),'o','LineWidth',2,'DisplayName','V_2')
    %     plot(t_soln5,x_soln5(:,P.V_3),'o','LineWidth',2,'DisplayName','V_3')
    %     plot(t_soln3,x_soln3(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
    %     plot(t_soln3,x_soln3(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
    %     plot(t_soln3,x_soln3(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')

end


%% Test getting SS systems
FLAG.C_mode = 4;
[sys_CT3 , sys_DT3 , sys_CT5 , sys_DT5] = getSS_System(SIM,N,P,FLAG);

%% Test 3CT SS by comparing step to ode
if false
    % FLAG.InputMode = 2; % 1) Step 2) Sine 3) Ramp
    [InputSignal] = getInputSignal(SIM,N,P,FLAG);
    % creating the ss system in matlab, the states are reduced, need to fix the IC
    % This method hasn't always worked
    y_0 = SIM.x_0_3;
    x_red = sys_CT3.C\y_0;

    [y_CT3,t_CT3,x_CT3] = lsim(sys_CT3,InputSignal(:,2),InputSignal(:,1),x_red);

    % [y_CT3,t_CT3,x_CT3] = step(sys_CT3,SIM.t_final_sim); Doesn't do IC which
    % is needed for 3 state system

    figure
    hold on
    plot(t_CT3,y_CT3(:,P.V_1),'o','LineWidth',2,'DisplayName','V_1')
    plot(t_CT3,y_CT3(:,P.V_2),'o','LineWidth',2,'DisplayName','V_2')
    plot(t_CT3,y_CT3(:,P.V_3),'o','LineWidth',2,'DisplayName','V_3')
    plot(t_soln3,x_soln3(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
    plot(t_soln3,x_soln3(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
    plot(t_soln3,x_soln3(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')
    title('Compare CT3')
    %
end


%% Test 3DT SS by comparing step to ode
if false
    % [y_DT3,t_DT3,x_DT3] = step(sys_DT3,SIM.t_final_sim);
    FLAG.InputMode = 2; % 1) Step 2) Sine 3) Ramp
    [InputSignal] = getInputSignal(SIM,N,P,FLAG);

    y_0 = SIM.x_0_3;
    x_red = sys_CT3.C\y_0;

    [y_DT3,t_DT3,x_DT3] = lsim(sys_DT3,InputSignal(:,2),InputSignal(:,1),x_red);

    figure
    hold on
    plot(t_DT3,y_DT3(:,P.V_1),'o','LineWidth',2,'DisplayName','V_1')
    plot(t_DT3,y_DT3(:,P.V_2),'o','LineWidth',2,'DisplayName','V_2')
    plot(t_DT3,y_DT3(:,P.V_3),'o','LineWidth',2,'DisplayName','V_3')
    plot(t_soln3,x_soln3(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
    plot(t_soln3,x_soln3(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
    plot(t_soln3,x_soln3(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')
    title('Compare DT3')
end


%% Test 5CT SS by comparing step to ode
%     [y_CT5,t_CT5,x_CT5] = step(sys_CT5,SIM.t_final_sim);
% Only work with zero IC

%     figure
%     hold on
%     plot(t_CT5,y_CT5(:,P.V_1),'o','LineWidth',2,'DisplayName','V_1')
%     plot(t_CT5,y_CT5(:,P.V_2),'o','LineWidth',2,'DisplayName','V_2')
%     plot(t_CT5,y_CT5(:,P.V_3),'o','LineWidth',2,'DisplayName','V_3')
%     plot(t_soln3,x_soln3(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
%     plot(t_soln3,x_soln3(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
%     plot(t_soln3,x_soln3(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')
%     title('Compare CT5')
%

%% Test 5DT SS by comparing step to ode
%     [y_DT5,t_DT5,x_DT5] = step(sys_DT5,SIM.t_final_sim);
% Only work with zero IC

%     figure
%     hold on
%     plot(t_DT5,y_DT5(:,P.V_1),'o','LineWidth',2,'DisplayName','V_1')
%     plot(t_DT5,y_DT5(:,P.V_2),'o','LineWidth',2,'DisplayName','V_2')
%     plot(t_DT5,y_DT5(:,P.V_3),'o','LineWidth',2,'DisplayName','V_3')
%     plot(t_soln3,x_soln3(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
%     plot(t_soln3,x_soln3(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
%     plot(t_soln3,x_soln3(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')
%     title('Compare DT5')


%% Test Converting Plant states into ROM IC
if false
    C_des       = sys_DT3.C;
    [U, S, V]   = svd(C_des);
    threshold   = 1e-7;
    S_cross     = pinv(S,threshold);
    C_des_cross = V*S_cross*U';

    x_plant_ROM = C_des_cross * z_soln3';

    % figure
    % hold on
    % plot(t_soln3,x_plant_ROM(1,:),'o','DisplayName','x_{plant ROM,1}' )
    % plot(t_soln3,x_plant_ROM(2,:),'o','DisplayName','x_{plant ROM,2}' )
    % plot(t_CT3,x_CT3(:,P.V_1),'k','LineWidth',2,'DisplayName','x_{ROM,1}')
    % plot(t_CT3,x_CT3(:,P.V_2),'k','LineWidth',2,'DisplayName','x_{ROM,2}')
    % lgn = legend;
end


%% Test SS DT using C_des matches all SS DT using C_m
if false
    FLAG.InputMode = 1; % 1) Step 2) Sine 3) Ramp
    [InputSignal] = getInputSignal(SIM,N,P,FLAG);

    % Using C_des to create SS DT
    FLAG.C_mode = 4;
    [A,B,C,D] = getAll_SS3(SIM,N,P,FLAG);
    sys = dss(A,B,C,D,SIM.Mass3);
    sys_CT3 = ss(sys,'explicit');
    sys_DT3 = c2d(sys_CT3,SIM.Ts);

    y_0 = SIM.x_0_3;
    x_red = sys_CT3.C\y_0

    % [y_DT3,t_DT3,x_DT3] = lsim(sys_DT3,InputSignal(:,2),InputSignal(:,1),x_red);

    % % Using C_m from C_des to create sys
    % for i = 1:3
    % idx = i;
    % FLAG.C_mode = idx;
    % [A,B,C,D] = getAll_SS3(SIM,N,P,FLAG);
    % sys = dss(A,B,C,D,SIM.Mass3);
    % sys_CT3 = ss(sys,'explicit');
    % sys_DT3 = c2d(sys_CT3,SIM.Ts);
    %
    % y_0 = SIM.x_0_3(idx);
    % x_red = sys_CT3.C\y_0
    %
    % [y_DT3_CM,t_DT3_CM,x_DT3_CM] = lsim(sys_DT3,InputSignal(:,2),InputSignal(:,1),x_red);
    %
    %     figure
    %     hold on
    %     plot(t_DT3,y_DT3(:,P.V_1),'o','LineWidth',2,'DisplayName','V_1')
    %     plot(t_DT3,y_DT3(:,P.V_2),'o','LineWidth',2,'DisplayName','V_2')
    %     plot(t_DT3,y_DT3(:,P.V_3),'o','LineWidth',2,'DisplayName','V_3')
    %     plot(t_DT3_CM,y_DT3_CM(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
    %     title('Compare C_{des} to C_{m} Outputs')
    %
    %
    %     figure
    %     hold on
    %     plot(t_DT3,x_DT3(:,P.V_1),'o','LineWidth',2,'DisplayName','X_1')
    %     plot(t_DT3,x_DT3(:,P.V_2),'o','LineWidth',2,'DisplayName','X_2')
    % %     plot(t_DT3,x_DT3(:,P.V_3),'o','LineWidth',2,'DisplayName','X_3')
    %     plot(t_DT3_CM,x_DT3_CM(:,P.V_1),'k','LineWidth',2,'DisplayName','X_1')
    %     plot(t_DT3_CM,x_DT3_CM(:,P.V_2),'k','LineWidth',2,'DisplayName','X_2')
    %     title('Compare C_{des} to C_{m} States')
    % end

    % Using C_m from C_des to create sys
    for i = 1:3
        idx = i;
        FLAG.C_mode = idx;
        C_m = sys_DT3.C(idx,:);
        y_DT3_CM = C_m*x_DT3';
        x_DT3_CM = x_DT3;
        t_DT3_CM = t_DT3;
        % [y_DT3_CM,t_DT3_CM,x_DT3_CM] = lsim(sys_DT3,InputSignal(:,2),InputSignal(:,1),x_red);

        %     figure
        %     hold on
        %     plot(t_DT3,y_DT3(:,P.V_1),'o','LineWidth',2,'DisplayName','V_1')
        %     plot(t_DT3,y_DT3(:,P.V_2),'o','LineWidth',2,'DisplayName','V_2')
        %     plot(t_DT3,y_DT3(:,P.V_3),'o','LineWidth',2,'DisplayName','V_3')
        %     plot(t_DT3_CM,y_DT3_CM(P.V_1,:),'k','LineWidth',2,'DisplayName','V_1')
        %     title('Compare C_{des} to C_{m} Outputs')


        %     figure
        %     hold on
        %     plot(t_DT3,x_DT3(:,P.V_1),'o','LineWidth',2,'DisplayName','X_1')
        %     plot(t_DT3,x_DT3(:,P.V_2),'o','LineWidth',2,'DisplayName','X_2')
        %     plot(t_DT3_CM,x_DT3_CM(:,P.V_1),'k','LineWidth',2,'DisplayName','X_1')
        %     plot(t_DT3_CM,x_DT3_CM(:,P.V_2),'k','LineWidth',2,'DisplayName','X_2')
        %     title('Compare C_{des} to C_{m} States')
    end

    % If I want to find the des outputs, I have to make the SS DT using all
    % outputs. Then in the plants, I have to set C_m to the state I am
    % currently measuring.
end


%% Test Simulink
if false
    % % Simulink Models to run
    % % 3) 3 state plant
    % % 5) 5 state plant %%!!! Not implemented yet
    % FLAG.C_mode = 4;
    % FLAG.SlinkModel = 3;
    % [RESULTS] = RunSimulink(SIM,N,P,FLAG,RESULTS);

    %%
    % figure
    % hold on
    % % plot(RESULTS.Slink_plant.x.t_soln , RESULTS.Slink_plant.x.x_soln(P.V_1,:),'o','LineWidth',2,'DisplayName','V_1')
    % % plot(RESULTS.Slink_plant.x.t_soln , RESULTS.Slink_plant.x.x_soln(P.V_2,:),'o','LineWidth',2,'DisplayName','V_2')
    % % plot(RESULTS.Slink_plant.x.t_soln , RESULTS.Slink_plant.x.x_soln(P.V_3,:),'o','LineWidth',2,'DisplayName','V_3')
    % plot(RESULTS.Slink_plant.xN.t_soln , RESULTS.Slink_plant.xN.x_soln(P.V_1,:),'o','LineWidth',2,'DisplayName','V_1')
    % plot(RESULTS.Slink_plant.xN.t_soln , RESULTS.Slink_plant.xN.x_soln(P.V_2,:),'o','LineWidth',2,'DisplayName','V_2')
    % plot(RESULTS.Slink_plant.xN.t_soln , RESULTS.Slink_plant.xN.x_soln(P.V_3,:),'o','LineWidth',2,'DisplayName','V_3')
    %
    % % plot(t_soln3,x_soln3(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
    % % plot(t_soln3,x_soln3(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
    % % plot(t_soln3,x_soln3(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')
    % plot(t_soln5,x_soln5(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
    % plot(t_soln5,x_soln5(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
    % plot(t_soln5,x_soln5(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')

    %% Save Simulink Outputs for Estimation
    % switch FLAG.SlinkModel
    %     case 1 % 3 state plant
    %         type_str = 'Plant3_Step';
    %         Q_str = num2str(0);
    %         R_str = num2str(0);
    %     case 2 % 3 state plant with input noise
    %         type_str = 'Plant3_IQ_Step';
    %         Q_str = num2str(SIM.Q_0);
    %         R_str = num2str(SIM.R_0);
    %     case 3 % 3 state plant with state noise
    %         type_str = 'Plant3_SQ_Step';
    %         Q_str = num2str(SIM.Q_0);
    %         R_str = num2str(SIM.R_0);
    % end
    %
    % Ts_str = num2str(SIM.Ts);
    % switch FLAG.C_mode
    %     case 1
    %         C_str = num2str(FLAG.C_mode);
    %     case 2
    %         C_str = num2str(FLAG.C_mode);
    %     case 3
    %         C_str = num2str(FLAG.C_mode);
    %     case 4
    %         C_str = num2str('All');
    % end
    %
    %
    % Folder = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanTesting\DAE_Circuit\SoftwareStyle\Results';
    % filename = [type_str '_Q' Q_str '_R' R_str '_Ts' Ts_str '_C' C_str '.mat'];
    % overall_filename = [Folder filesep filename];
    % save(overall_filename, "z", "x","u" ,"w" ,"v" , "SIM", "FLAG", "sys_DT5", "sys_DT3","sys_CT3", "sys_CT5")
end

%% Test Estimation
if false
    % [x_hat_asy, x_hat_var, y_hat_asy, y_hat_var, P_infty, K_infty, K_k, P_k_pre] = PerformEstimation(overall_filename);
    %
    %     figure
    %     hold on
    %     plot(x.time,y_hat_asy(P.V_1,:),'o','LineWidth',2,'DisplayName','V_1')
    %     plot(x.time,y_hat_asy(P.V_2,:),'o','LineWidth',2,'DisplayName','V_2')
    %     plot(x.time,y_hat_asy(P.V_3,:),'o','LineWidth',2,'DisplayName','V_3')
    %     plot(t_CT3,y_CT3(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
    %     plot(t_CT3,y_CT3(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
    %     plot(t_CT3,y_CT3(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')
    %     title('Compare Asy Estimation')
    %
    %     figure
    %     hold on
    %     plot(x.time,y_hat_var(P.V_1,:),'o','LineWidth',2,'DisplayName','V_1')
    %     plot(x.time,y_hat_var(P.V_2,:),'o','LineWidth',2,'DisplayName','V_2')
    %     plot(x.time,y_hat_var(P.V_3,:),'o','LineWidth',2,'DisplayName','V_3')
    %     plot(t_CT3,y_CT3(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
    %     plot(t_CT3,y_CT3(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
    %     plot(t_CT3,y_CT3(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')
    %     title('Compare Var Estimation')
    %
    %     CPCT_calc = sys_CT3.C * P_infty * sys_CT3.C';
    %     disp('CPCT_calc=')
    %     disp(num2str(CPCT_calc))
    %     disp(newline)
    %     disp('P_k_pre final value')
    %     disp(num2str(sys_CT3.C *P_k_pre(:,:,end)* sys_CT3.C'))


    %% Test Covariance Calc
    % idx = 1000;
    % [Pcalc_var] = getP_calc(z.value,y_hat_var,idx)
    % [Pcalc_asy] = getP_calc(z.value,y_hat_asy,idx)


    % CPcalcCT_var = sys_CT3.C * Pcalc_var * sys_CT3.C';
    % CPcalcCT_asy = sys_CT3.C * Pcalc_asy * sys_CT3.C';
    % disp(newline)
    % disp(num2str(CPcalcCT_var))
    % disp(newline)
    % disp(num2str(CPcalcCT_asy))
end

%% Test Ho-Kalman
if true
    [sys_HK] = getHoKalmanROM(SIM,N,P,FLAG);

    % IC
    y_0 = SIM.x_0_5(1:3,1);
    x_red = sys_HK.C\y_0;

    % u_k
    [InputSignal] = getInputSignal(SIM,N,P,FLAG);

    % Simulate
    [y_HK,t_HK,x_HK] = lsim(sys_HK,InputSignal(:,2),InputSignal(:,1),x_red);

    % Results
    figure
    hold on
    plot(t_HK,y_HK(:,P.V_1),'b','LineWidth',2,'DisplayName','V_1')
    plot(t_HK,y_HK(:,P.V_2),'r','LineWidth',2,'DisplayName','V_2')
    plot(t_HK,y_HK(:,P.V_3),'y','LineWidth',2,'DisplayName','V_3')
    title('Just Ho-Kalman Results')
    lgn = legend;


    % Compare 3 State ODE
    figure
    hold on
    plot(t_HK,y_HK(:,P.V_1),'o','LineWidth',2,'DisplayName','V_1')
    plot(t_HK,y_HK(:,P.V_2),'o','LineWidth',2,'DisplayName','V_2')
    plot(t_HK,y_HK(:,P.V_3),'o','LineWidth',2,'DisplayName','V_3')
    plot(t_soln3,z_soln3(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
    plot(t_soln3,z_soln3(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
    plot(t_soln3,z_soln3(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')
    title('Compare Ho-Kalman to ode3')


    % Compare 5 State ODE
    figure
    hold on
    plot(t_HK,y_HK(:,P.V_1),'o','LineWidth',2,'DisplayName','V_1')
    plot(t_HK,y_HK(:,P.V_2),'o','LineWidth',2,'DisplayName','V_2')
    plot(t_HK,y_HK(:,P.V_3),'o','LineWidth',2,'DisplayName','V_3')
    plot(t_soln5,z_soln5(:,P.V_1),'k','LineWidth',2,'DisplayName','V_1')
    plot(t_soln5,z_soln5(:,P.V_2),'k','LineWidth',2,'DisplayName','V_2')
    plot(t_soln5,z_soln5(:,P.V_3),'k','LineWidth',2,'DisplayName','V_3')
    title('Compare Ho-Kalman to ode5')
end


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
