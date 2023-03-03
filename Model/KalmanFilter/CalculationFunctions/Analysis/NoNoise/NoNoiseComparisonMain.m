%% No Noise Comparison Main

function [SIM,FLAG,N,P,RESULTS] = NoNoiseComparisonMain(SIM,FLAG,N,P,RESULTS)
%% Change Input Signal to be evenly spaced 
t_vec  = (0 : SIM.Ts : SIM.InputSignal(end,1))';
i_user = interp1(SIM.InputSignal(:,1) , SIM.InputSignal(:,2) , t_vec);

NewInputSignal = [t_vec , i_user];

%% ODE
if FLAG.Analysis.ode
    [RESULTS.ode.t_soln, RESULTS.ode.x_soln, RESULTS.ode.z_soln, ~, ~] = getNoNoiseODE(SIM,FLAG,N,P);
end


%% State Space Continuous Time
if FLAG.Analysis.SS_CT   
    [sys_CT , ~ ] = getSS_System(SIM,N,P,FLAG);

    % Initial Conditions for the SS
        x_red = zeros(length(sys_CT.A),1);

    % Run SS Simulation
        %[RESULTS.SS_CT.z_soln,RESULTS.SS_CT.t_soln,RESULTS.SS_CT.x_soln] = lsim(sys_CT,SIM.InputSignal(:,2),SIM.InputSignal(:,1),x_red);
        [RESULTS.SS_CT.z_soln,RESULTS.SS_CT.t_soln,RESULTS.SS_CT.x_soln] = lsim(sys_CT,NewInputSignal(:,2),NewInputSignal(:,1),x_red);

    % Add Initial Offset back to SS
        RESULTS.SS_CT.z_soln = RESULTS.SS_CT.z_soln + SIM.y_0_FOM'; 
end


%% State Space Discrete Time
if FLAG.Analysis.SS_DT   
    [~ , sys_DT] = getSS_System(SIM,N,P,FLAG);

    % Initial Conditions for the SS
        x_red = zeros(length(sys_DT.A),1);

    % Run SS Simulation
%         [RESULTS.SS_DT.z_soln , RESULTS.SS_DT.t_soln , RESULTS.SS_DT.x_soln] = lsim(sys_DT,SIM.InputSignal(:,2),SIM.InputSignal(:,1),x_red);
        [RESULTS.SS_DT.z_soln , RESULTS.SS_DT.t_soln , RESULTS.SS_DT.x_soln] = lsim(sys_DT,NewInputSignal(:,2),NewInputSignal(:,1),x_red);
    
    % Add Initial Offset back to SS
        RESULTS.SS_DT.z_soln = RESULTS.SS_DT.z_soln + SIM.y_0_FOM';
end


%% Optimal Ho-Kalman
if FLAG.Analysis.OptimalHK
    oldFLAGUseInput_r = FLAG.UseInput_r;
    FLAG.UseInput_r = 1;

    tic
    HK_ErrorCalc(SIM,N,P,FLAG,RESULTS);
    toc

    FLAG.UseInput_r = oldFLAGUseInput_r;
end


%% Ho-Kalman
if FLAG.Analysis.ROM_HoKal
    oldFLAGUseInput_r = FLAG.UseInput_r;
    FLAG.UseInput_r = 0;
    
    [sys_HK,~] = getHoKalmanROM(SIM,N,P,FLAG,RESULTS);

    FLAG.UseInput_r = oldFLAGUseInput_r;

    if FLAG.EST.SepHK
        for i = 1:N.DesOut  %length(sys_HK)-1
            % Initial Conditions for the SS
            x_red = zeros(length(sys_HK{i}.A),1);

            % Run Simulation
            [z_soln,t_soln,~] = lsim(sys_HK{i},NewInputSignal(:,2),NewInputSignal(:,1),x_red);
            if i == 1
                RESULTS.ROM_HoKal.t_soln      = t_soln;
            end
            RESULTS.ROM_HoKal.z_soln(:,i) = z_soln(:,2);
        end
    else
        % Initial Conditions for the SS
        x_red = zeros(length(sys_HK{end}.A),1);

        % Run Simulation
        [RESULTS.ROM_HoKal.z_soln,RESULTS.ROM_HoKal.t_soln,RESULTS.ROM_HoKal.x_soln] = lsim(sys_HK{end},NewInputSignal(:,2),NewInputSignal(:,1),x_red);

    end
    % Add Initial Offset back to SS
    RESULTS.ROM_HoKal.z_soln = RESULTS.ROM_HoKal.z_soln + SIM.y_0_FOM';

    %old Stuff
%             [z_soln,t_soln,~] = lsim(sys_HK{i},SIM.InputSignal(:,2),SIM.InputSignal(:,1),x_red);
%                 RESULTS.ROM_HoKal.z_soln(:,1) = z_soln(:,1);
            %RESULTS.ROM_HoKal.x_soln
    %         [RESULTS.ROM_HoKal.z_soln,RESULTS.ROM_HoKal.t_soln,RESULTS.ROM_HoKal.x_soln] = lsim(sys_HK{end},SIM.InputSignal(:,2),SIM.InputSignal(:,1),x_red);
end




end


%%%%%%%%%%%%%%%%%%%%%%%%%% OOOOOOOOOOOLD %%%%%%%%%%%%%%%%%%%%%%%%%% 
%   
%% 
% SOC = SIM.SOC;
% Ts  = SIM.Ts;



%     filepath = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanFilter\DATA\ODE_Step';
%     filename = ['Step_SOC' num2str(SOC) '_Ts' num2str(Ts) '.mat'];
%     sim = load([filepath filesep filename]);
%     RESULTS.ode.t_soln = sim.t_sim_vec;
%     RESULTS.ode.x_soln = sim.x;
%     RESULTS.ode.z_soln = sim.z;    

%%
 %&& ~SIM.AnalysisComplete.NoNoiseCompare.SS_CT
% && ~SIM.AnalysisComplete.NoNoiseCompare.SS_DT
% && ~SIM.AnalysisComplete.NoNoiseCompare.ROM_HoKal



 %% 
% if FLAG.Analysis.ROM_Mlab && ~SIM.AnalysisComplete.NoNoiseCompare.ROM_Mlab
%     [~ , sys_DT ] = getSS_System(SIM,N,P,FLAG);
%     
%     FLAG.InputType = 3; % 1) Step 2) Sine 3) Ramp
%     [InputSignal] = getInputSignal(SIM,N,P,FLAG);
%         
%     y_0 = SIM.x_0;
%     x_red = sys_DT.C\y_0;
%     
%     [RESULTS.ROM_Mlab.z_soln , RESULTS.ROM_Mlab.t_soln , RESULTS.ROM_Mlab.x_soln] = lsim(sys_DT,InputSignal(:,2),InputSignal(:,1),x_red);
% 
% 
%     SIM.AnalysisComplete.NoNoiseCompare.ROM_Mlab = 1;
% end

%%


%     FLAG.InputType = 3; % 1) Step 2) Sine 3) Ramp
%     [InputSignal] = getInputSignal(SIM,N,P,FLAG);

%     y_0 = SIM.x_0;
%     x_red = sys_CT.C\SIM.y_0_FOM;


%     SIM.AnalysisComplete.NoNoiseCompare.SS_CT = 1;

%%
    
%     FLAG.InputType = 3; % 1) Step 2) Sine 3) Ramp
%     [InputSignal] = getInputSignal(SIM,N,P,FLAG);
        
%     y_0 = SIM.x_0;
%     x_red = sys_DT.C\SIM.y_0_FOM;


%     SIM.AnalysisComplete.NoNoiseCompare.SS_DT = 1;

%%


%     % IC
%     y_0 = SIM.x_0;
%     x_red = sys_HK.C\y_0;

%     % u_k
%     FLAG.InputType = 3; % 1) Step 2) Sine 3) Ramp
%     [InputSignal] = getInputSignal(SIM,N,P,FLAG);
%     SIM.AnalysisComplete.NoNoiseCompare.ROM_HoKal = 1;