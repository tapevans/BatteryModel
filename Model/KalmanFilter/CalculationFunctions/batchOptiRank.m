%% Optimal Rank
%
function batchOptiRank(SIM,FLAG,N,P,RESULTS)
%% Initialize
% SOC_vec = [50];
SOC_vec = 0:5:100;
Ts_vec  = [1];
rankMax = 20;

%% Set FLAGS
FLAG.OverwriteData.OptiBatch = 1;
FLAG.folderpathOptiBatch = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanFilter\Results\Optimal_Rank';
FLAG.PRBSAmp = 1;
FLAG.Tswitch = 10;
FLAG.UseInput_r = 1;


% Loop through SOC
for SS = SOC_vec
    disp(['SOC: ' num2str(SS)])
    disp(datestr(datetime));
    % Loop through Ts
    for TT = Ts_vec
        FLAG.SOC = SS;
        FLAG.Ts  = TT;
        %% Check if batch exists
        RUNSIM = false;
        [batch_filename] = getBatchFilename(SIM,FLAG);
        if isfile(batch_filename)
            if FLAG.OverwriteData.OptiBatch
                disp('Overwriting Batch File')
                RUNSIM = true;
                delete(batch_filename)
            else
                disp('Optimal Batch File Exists')
            end
        else
            RUNSIM = true;
        end

        %% Create Batch Data
        if RUNSIM
            % Get ODE simulation
            [t_soln_ODE, x_soln_ODE, z_soln_ODE, i_user_time_ODE, i_user_value_ODE] = getNoNoiseODE(SIM,FLAG,N,P);

            % DT Input Signal
            if TT == SIM.Ts
                NewInputSignal = SIM.InputSignal;
            else
                disp('Wrong Ts')
                break
            end

            % Loop through Desired Outputs
            for OO = 1:N.DesOut+1
                % Loop through Ranks
                for RR = 1:rankMax
                    SIM.Input_r = RR;
                    % get HK
                        [sys_HK,~] = getHoKalmanROM(SIM,N,P,FLAG,RESULTS);
                    
                    % Simulate ROM
                        if OO < N.DesOut+1
                            x_red = zeros(length(sys_HK{OO}.A),1);
                            [z_soln, t_soln, ~] = lsim(sys_HK{OO} , NewInputSignal(:,2) , NewInputSignal(:,1) , x_red);
                            % if OO == 1
                            individualOutputs.t_soln = t_soln;
                            % end
                            individualOutputs.z_soln = z_soln;
                        else
                            x_red = zeros(length(sys_HK{end}.A),1);
                            [combinedOutputs.z_soln, combinedOutputs.t_soln, combinedOutputs.x_soln] = lsim(sys_HK{end} , NewInputSignal(:,2) , NewInputSignal(:,1),x_red);
                        end
                    
                    % Add Initial Off set
                        if OO ~= N.DesOut+1
                            individualOutputs.z_soln = individualOutputs.z_soln + SIM.y_0_FOM([1,OO],1)';
                        else
                            combinedOutputs.z_soln   = combinedOutputs.z_soln   + SIM.y_0_FOM';
                        end

                    % Solve for ODE at DT
                        if OO ~= N.DesOut+1
                            z_interpODE(1,:) = interp1(t_soln_ODE, z_soln_ODE(1,:)  , NewInputSignal(:,1));
                            z_interpODE(2,:) = interp1(t_soln_ODE, z_soln_ODE(OO,:) , NewInputSignal(:,1));
                        else
                            for OOO = 1:N.DesOut
                                z_interpODE(OOO,:) = interp1(t_soln_ODE, z_soln_ODE(OOO,:) , NewInputSignal(:,1));
                            end
                        end

                    % Calculate Error
                        if OO < N.DesOut+1
                        % Individual
                            individualOutputs.Error      =     individualOutputs.z_soln'  - z_interpODE;
                            individualOutputs.ErrorSQ    =    (individualOutputs.Error).^2;
                            individualOutputs.ErrorSQSum = sum(individualOutputs.ErrorSQ,2);
                            ind.saveErrorSqSum(:,RR)     =     individualOutputs.ErrorSQSum;
                        else
                        % Combined
                            combinedOutputs.Error      =     combinedOutputs.z_soln'    - z_interpODE;
                            combinedOutputs.ErrorSQ    =    (combinedOutputs.Error).^2;
                            combinedOutputs.ErrorSQSum = sum(combinedOutputs.ErrorSQ,2);
                            com.saveErrorSqSum(:,RR)   =     combinedOutputs.ErrorSQSum;
                        end
                end

                % Save Error for each Output
                if  OO < N.DesOut+1
                    All_ErrorSqSum{OO} = ind.saveErrorSqSum;
                else
                    All_ErrorSqSum{OO} = com.saveErrorSqSum;
                end
                
                % Determine the idx (rank) of each minimum error
                if OO < N.DesOut+1
                    [min_v , idx_volt] = min(ind.saveErrorSqSum(1,:));
                    [min_o , idx_othe] = min(ind.saveErrorSqSum(2,:));
                    MinAndIDX{OO}.ind.min_v    = min_v;
                    MinAndIDX{OO}.ind.idx_volt = idx_volt;
                    MinAndIDX{OO}.ind.min_o    = min_o;
                    MinAndIDX{OO}.ind.idx_othe = idx_othe;
                else
                    [min_com , idx_com] = min(com.saveErrorSqSum,[],2);
                    MinAndIDX{OO}.com.min_com  = min_com;
                    MinAndIDX{OO}.com.idx_com  = idx_com;
                end
            end

            % Save Results
            save(batch_filename,'All_ErrorSqSum','MinAndIDX','SS','TT','RESULTS','P','N')
        end
        clear z_interpODE All_ErrorSqSum MinAndIDX individualOutputs combinedOutputs ind com
    end
end

end