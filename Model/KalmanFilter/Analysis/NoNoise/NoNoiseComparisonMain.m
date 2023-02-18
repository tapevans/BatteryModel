%% No Noise Comparison Main

function [SIM,FLAG,RESULTS] = NoNoiseComparisonMain(SIM,FLAG)%,N,P,RESULTS   N,P,
%% 
SOC = SIM.SOC;
Ts  = SIM.Ts;
if FLAG.Analysis.ode %&& ~SIM.AnalysisComplete.NoNoiseCompare.ode
    filepath = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanFilter\DATA\ODE_Step';
    filename = ['Step_SOC' num2str(SOC) '_Ts' num2str(Ts) '.mat'];
    sim = load([filepath filesep filename]);
    RESULTS.ode.t_soln = sim.t_sim_vec;
    RESULTS.ode.x_soln = sim.x;
    RESULTS.ode.z_soln = sim.z;    
end
% %% 
% if FLAG.Analysis.SS_CT    && ~SIM.AnalysisComplete.NoNoiseCompare.SS_CT
%     [sys_CT , ~ ] = getSS_System(SIM,N,P,FLAG);
% 
%     FLAG.InputType = 3; % 1) Step 2) Sine 3) Ramp
%     [InputSignal] = getInputSignal(SIM,N,P,FLAG);
% 
%     y_0 = SIM.x_0;
%     x_red = sys_CT.C\y_0;
% 
%     [RESULTS.SS_CT.z_soln,RESULTS.SS_CT.t_soln,RESULTS.SS_CT.x_soln] = lsim(sys_CT,InputSignal(:,2),InputSignal(:,1),x_red);
% 
%     SIM.AnalysisComplete.NoNoiseCompare.SS_CT = 1;
% end
% %% 
% if FLAG.Analysis.SS_DT    && ~SIM.AnalysisComplete.NoNoiseCompare.SS_DT
%     [~ , sys_DT] = getSS_System(SIM,N,P,FLAG);
%     
%     FLAG.InputType = 3; % 1) Step 2) Sine 3) Ramp
%     [InputSignal] = getInputSignal(SIM,N,P,FLAG);
%         
%     y_0 = SIM.x_0;
%     x_red = sys_DT.C\y_0;
%     
%     [RESULTS.SS_DT.z_soln , RESULTS.SS_DT.t_soln , RESULTS.SS_DT.x_soln] = lsim(sys_DT,InputSignal(:,2),InputSignal(:,1),x_red);
% 
% 
%     SIM.AnalysisComplete.NoNoiseCompare.SS_DT = 1;
% end
% %% 
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
% if FLAG.Analysis.ROM_HoKal% && ~SIM.AnalysisComplete.NoNoiseCompare.ROM_HoKal
%     [sys_HK] = getHoKalmanROM(SIM,N,P,FLAG);
% 
%     % IC
%     y_0 = SIM.x_0;
%     x_red = sys_HK.C\y_0;
% 
%     % u_k
%     FLAG.InputType = 3; % 1) Step 2) Sine 3) Ramp
%     [InputSignal] = getInputSignal(SIM,N,P,FLAG);
% 
%     % Simulate
%     [RESULTS.ROM_HoKal.z_soln,RESULTS.ROM_HoKal.t_soln,RESULTS.ROM_HoKal.x_soln] = lsim(sys_HK,InputSignal(:,2),InputSignal(:,1),x_red);
%     SIM.AnalysisComplete.NoNoiseCompare.ROM_HoKal = 1;
% end
%% 


end