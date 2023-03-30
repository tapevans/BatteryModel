%% No Noise Comparison Main

function [SIM,FLAG,N,P,RESULTS] = NoNoiseComparisonMain(SIM,FLAG,N,P,RESULTS)
%% Change Input Signal to be evenly spaced 
    t_vec  = (0 : SIM.Ts : SIM.InputSignal(end,1))';
    i_user = interp1(SIM.InputSignal(:,1) , SIM.InputSignal(:,2) , t_vec);
    
    NewInputSignal = [t_vec , i_user]; % I think this ensures proper DT signal

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
            [RESULTS.SS_DT.z_soln , RESULTS.SS_DT.t_soln , RESULTS.SS_DT.x_soln] = lsim(sys_DT,NewInputSignal(:,2),NewInputSignal(:,1),x_red);
        
        % Add Initial Offset back to SS
            RESULTS.SS_DT.z_soln = RESULTS.SS_DT.z_soln + SIM.y_0_FOM';
    end


%% Optimal Ho-Kalman
    if FLAG.Analysis.OptimalHK
        oldFLAGUseInput_r = FLAG.UseInput_r;
        FLAG.UseInput_r = 1;
    
        %tic
        %HK_OptiStates(SIM,N,P,FLAG,RESULTS);
        % HK_OptiStatesTest(SIM,N,P,FLAG,RESULTS);
        %toc
    
        % batchOptiRank(SIM,FLAG,N,P,RESULTS)

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
            [RESULTS.ROM_HoKal.z_soln , RESULTS.ROM_HoKal.t_soln , RESULTS.ROM_HoKal.x_soln] = lsim(sys_HK{end},NewInputSignal(:,2),NewInputSignal(:,1),x_red);
    
        end
        % Add Initial Offset back to SS
        RESULTS.ROM_HoKal.z_soln = RESULTS.ROM_HoKal.z_soln + SIM.y_0_FOM';
    end
end