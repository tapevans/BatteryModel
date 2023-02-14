%% No Noise Comparison Main

function [SIM,FLAG,N,P,RESULTS] = NoNoiseComparisonMain(SIM,FLAG,N,P,RESULTS)
%% 
if FLAG.Analysis.ode3 && ~SIM.AnalysisComplete.NoNoiseCompare.ode3
    FLAG.StateMode = 3; % 3 or 5 depending on which system to use
    FLAG.InputType = 3; % 1) Step 2) Sine 3) Ramp
    [RESULTS.ode3.t_soln, RESULTS.ode3.x_soln, RESULTS.ode3.z_soln] = runODE(SIM,N,P,FLAG);

    SIM.AnalysisComplete.NoNoiseCompare.ode3 = 1;
end
%% 
if FLAG.Analysis.ode5    &&   ~SIM.AnalysisComplete.NoNoiseCompare.ode5
    FLAG.StateMode = 5; % 3 or 5 depending on which system to use
    FLAG.InputType = 1; % 1) Step 2) Sine 3) Ramp
    [RESULTS.ode5.t_soln, RESULTS.ode5.x_soln, RESULTS.ode5.z_soln] = runODE(SIM,N,P,FLAG);

    SIM.AnalysisComplete.NoNoiseCompare.ode5 = 1;

end
%% 
if FLAG.Analysis.SS_CT3    && ~SIM.AnalysisComplete.NoNoiseCompare.SS_CT3
    [sys_CT3 , ~ , ~ , ~] = getSS_System(SIM,N,P,FLAG);

    FLAG.InputType = 1; % 1) Step 2) Sine 3) Ramp
    [InputSignal] = getInputSignal(SIM,N,P,FLAG);

    y_0 = SIM.x_0_3;
    x_red = sys_CT3.C\y_0;

    [RESULTS.SS_CT3.z_soln,RESULTS.SS_CT3.t_soln,RESULTS.SS_CT3.x_soln] = lsim(sys_CT3,InputSignal(:,2),InputSignal(:,1),x_red);

    SIM.AnalysisComplete.NoNoiseCompare.SS_CT3 = 1;
end
%% 
if FLAG.Analysis.SS_DT3    && ~SIM.AnalysisComplete.NoNoiseCompare.SS_DT3
    [~ , sys_DT3 , ~ , ~] = getSS_System(SIM,N,P,FLAG);
    
    FLAG.InputType = 1; % 1) Step 2) Sine 3) Ramp
    [InputSignal] = getInputSignal(SIM,N,P,FLAG);
        
    y_0 = SIM.x_0_3;
    x_red = sys_DT3.C\y_0;
    
    [RESULTS.SS_DT3.z_soln , RESULTS.SS_DT3.t_soln , RESULTS.SS_DT3.x_soln] = lsim(sys_DT3,InputSignal(:,2),InputSignal(:,1),x_red);


    SIM.AnalysisComplete.NoNoiseCompare.SS_DT3 = 1;
end
%% 
if FLAG.Analysis.ROM_Mlab3 && ~SIM.AnalysisComplete.NoNoiseCompare.ROM_Mlab3
    [~ , sys_DT3 , ~ , ~] = getSS_System(SIM,N,P,FLAG);
    
    FLAG.InputType = 1; % 1) Step 2) Sine 3) Ramp
    [InputSignal] = getInputSignal(SIM,N,P,FLAG);
        
    y_0 = SIM.x_0_3;
    x_red = sys_DT3.C\y_0;
    
    [RESULTS.ROM_Mlab3.z_soln , RESULTS.ROM_Mlab3.t_soln , RESULTS.ROM_Mlab3.x_soln] = lsim(sys_DT3,InputSignal(:,2),InputSignal(:,1),x_red);


    SIM.AnalysisComplete.NoNoiseCompare.ROM_Mlab3 = 1;
end
%% 
if FLAG.Analysis.ROM_HoKal && ~SIM.AnalysisComplete.NoNoiseCompare.ROM_HoKal
    %%%%%%%%%%%%%%%%%%%%
    SIM.AnalysisComplete.NoNoiseCompare.ROM_HoKal = 1;
end
%% 


end