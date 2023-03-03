%% ROM Error Calc
% Trying to determine if there is a trend between how many states to use for
% the ROM and Ts, SOC, or Desired Variable
function HK_ErrorCalc(SIM,N,P,FLAG,RESULTS)
%% Check if file exists
    RUNSIM = false;
    [ROMError_filename] = getROMErrorFilename(FLAG);
    if isfile(ROMError_filename)
        if FLAG.OverwriteData.ROMError
            disp('Overwriting ROM Error File')
            RUNSIM = true;
            delete(ROMError_filename)
        else
            disp('ROM Error File Exists') %%%%!! May add something here later to check if certain des outputs are a part of the saved data
        end
    else
        RUNSIM = true;
    end
    if RUNSIM
        FLAG.Analysis.PlotImp = 0;
        %% Remove fields from RESULTS
        if isfield(RESULTS,'ROM_HoKal')
            RESULTS = rmfield(RESULTS,'ROM_HoKal');
        end
        if isfield(RESULTS,'SS_DT')
            RESULTS = rmfield(RESULTS,'SS_DT');
        end
        if isfield(RESULTS,'SS_CT')
            RESULTS = rmfield(RESULTS,'SS_CT');
        end
        if isfield(RESULTS,'ode')
            RESULTS = rmfield(RESULTS,'ode');
        end
        if isfield(RESULTS,'foldername')
            RESULTS = rmfield(RESULTS,'foldername');
        end

        %% Initialize Save Variables
        singularValue      = nan(1 ,SIM.r_max ,N.DesOut+1);

        error.ind.step     = nan(1 ,SIM.r_max ,N.DesOut+1);
        error.ind.T10      = nan(1 ,SIM.r_max ,N.DesOut+1);
        error.ind.T100     = nan(1 ,SIM.r_max ,N.DesOut+1);
        error.ind.AllSims  = nan(1 ,SIM.r_max ,N.DesOut+1);
        error.com.step     = nan(1 ,SIM.r_max ,N.DesOut+1);
        error.com.T10      = nan(1 ,SIM.r_max ,N.DesOut+1);
        error.com.T100     = nan(1 ,SIM.r_max ,N.DesOut+1);
        error.com.AllSims  = nan(1 ,SIM.r_max ,N.DesOut+1);

        optimal_r.ind.step    = nan(N.DesOut,1);
        optimal_r.ind.T10     = nan(N.DesOut,1);
        optimal_r.ind.T100    = nan(N.DesOut,1);
        optimal_r.ind.AllSims = nan(N.DesOut,1);
        optimal_r.com.step    = nan(N.DesOut,1);
        optimal_r.com.T10     = nan(N.DesOut,1);
        optimal_r.com.T100    = nan(N.DesOut,1);
        optimal_r.com.AllSims = nan(N.DesOut,1);

        optimal_r.AllOutputs.ind.step    = nan(1,1);
        optimal_r.AllOutputs.ind.T10     = nan(1,1);
        optimal_r.AllOutputs.ind.T100    = nan(1,1);
        optimal_r.AllOutputs.ind.AllSims = nan(1,1);
        optimal_r.AllOutputs.com.step    = nan(1,1);
        optimal_r.AllOutputs.com.T10     = nan(1,1);
        optimal_r.AllOutputs.com.T100    = nan(1,1);
        optimal_r.AllOutputs.com.AllSims = nan(1,1);


        %% Load Comparison Data
        % Step
            FLAG.InputMode = 1;
            [step.t_soln,       step.x_soln,       step.z_soln,       i_user_time, i_user_value] = getNoNoiseODE(SIM,FLAG,N,P);
            % Input Signal
                t_final_Relax_Step = (FLAG.N_samples+2)*SIM.Ts;
                t_vec_Relax_Step = 0:SIM.Ts:t_final_Relax_Step;
                % Add a time right after the step
                t_vec_Relax_Step = [t_vec_Relax_Step(1:3) , t_vec_Relax_Step(3)+SIM.Ts/10 , t_vec_Relax_Step(4:end)];
                input_load = ones(size(SIM.t_vec_Relax_Step));
                input_load(1) = 0;
                input_load(2) = 0;
                input_load(3) = 0;            
                if FLAG.SOC < 91
                    input_load = -1*input_load;
                end
                %t_vec  = (0 : SIM.Ts : i_user_time(end))';
                %i_user = interp1(i_user_time , i_user_value , t_vec);
                %step.NewInputSignal = [t_vec , i_user];
                InputSignal = [t_vec_Relax_Step' input_load'];

                t_vec  = (0 : SIM.Ts : SIM.InputSignal(end,1))';
                i_user = interp1(InputSignal(:,1) , InputSignal(:,2) , t_vec);
                step.NewInputSignal = [t_vec i_user];

        % Tswitch10
            FLAG.InputMode = 5;
            FLAG.PRBSAmp   = 1;
            FLAG.Tswitch   = 10;
            [Tswitch10.t_soln,  Tswitch10.x_soln,  Tswitch10.z_soln,  i_user_time, i_user_value] = getNoNoiseODE(SIM,FLAG,N,P);
            % Input Signal 
                t_vec  = (0 : SIM.Ts : i_user_time(end))';
                i_user = interp1(i_user_time , i_user_value , t_vec);
                Tswitch10.NewInputSignal = [t_vec , i_user];

        % Tswitch100
            FLAG.InputMode = 5;
            FLAG.PRBSAmp   = 1;
            FLAG.Tswitch   = 100;
            [Tswitch100.t_soln, Tswitch100.x_soln, Tswitch100.z_soln, i_user_time, i_user_value] = getNoNoiseODE(SIM,FLAG,N,P);
            % Input Signal 
                t_vec  = (0 : SIM.Ts : i_user_time(end))';
                i_user = interp1(i_user_time , i_user_value , t_vec);
                Tswitch100.NewInputSignal = [t_vec , i_user];


        %% Loop through all r
        for RR = 1:SIM.r_max
            %% Create ROM
            SIM.Input_r = RR;
            [sys_HK,singVals] = getHoKalmanROM(SIM,N,P,FLAG,RESULTS);
            

            %% Get Singular Values
            for OO = 1:N.DesOut+1
                singularValue(1,RR,OO) = singVals(1,OO);
            end


            %% Get ROM simulation data for each test for all desired outputs
            for OO = 1:N.DesOut+1
                %% Step Data
                if OO < N.DesOut+1
                    x_red = zeros(length(sys_HK{OO}.A),1);
                    [z_soln, t_soln, ~] = lsim(sys_HK{OO},step.NewInputSignal(:,2),step.NewInputSignal(:,1),x_red);
                    if OO == 1
                        step.individualOutputs.t_soln   = t_soln;
                    end
                    step.individualOutputs.z_soln(:,OO) = z_soln(:,2);
                else
                    x_red = zeros(length(sys_HK{end}.A),1);
                    [step.combinedOutputs.z_soln, step.combinedOutputs.t_soln, step.combinedOutputs.x_soln] = lsim(sys_HK{end},step.NewInputSignal(:,2),step.NewInputSignal(:,1),x_red);
                end
                    
                %% Tswitch10
                if OO < N.DesOut+1
                    x_red = zeros(length(sys_HK{OO}.A),1);
                    [z_soln, t_soln, ~] = lsim(sys_HK{OO},Tswitch10.NewInputSignal(:,2),Tswitch10.NewInputSignal(:,1),x_red);
                    if OO == 1
                        Tswitch10.individualOutputs.t_soln   = t_soln;
                    end
                    Tswitch10.individualOutputs.z_soln(:,OO) = z_soln(:,2);
                else
                    x_red = zeros(length(sys_HK{end}.A),1);
                    [Tswitch10.combinedOutputs.z_soln, Tswitch10.combinedOutputs.t_soln, Tswitch10.combinedOutputs.x_soln] = lsim(sys_HK{end},Tswitch10.NewInputSignal(:,2),Tswitch10.NewInputSignal(:,1),x_red);
                end
    
                %% Tswitch100
                if OO < N.DesOut+1
                    x_red = zeros(length(sys_HK{OO}.A),1);
                    [z_soln, t_soln, ~] = lsim(sys_HK{OO},Tswitch100.NewInputSignal(:,2),Tswitch100.NewInputSignal(:,1),x_red);
                    if OO == 1
                        Tswitch100.individualOutputs.t_soln   = t_soln;
                    end
                    Tswitch100.individualOutputs.z_soln(:,OO) = z_soln(:,2);
                else
                    x_red = zeros(length(sys_HK{end}.A),1);
                    [Tswitch100.combinedOutputs.z_soln, Tswitch100.combinedOutputs.t_soln, Tswitch100.combinedOutputs.x_soln] = lsim(sys_HK{end},Tswitch100.NewInputSignal(:,2),Tswitch100.NewInputSignal(:,1),x_red);
                end
            end
            % Add Initial Offset back to SS
            step.individualOutputs.z_soln = step.individualOutputs.z_soln + SIM.y_0_FOM';
            step.combinedOutputs.z_soln   = step.combinedOutputs.z_soln   + SIM.y_0_FOM';

            Tswitch10.individualOutputs.z_soln = Tswitch10.individualOutputs.z_soln + SIM.y_0_FOM';
            Tswitch10.combinedOutputs.z_soln   = Tswitch10.combinedOutputs.z_soln   + SIM.y_0_FOM';

            Tswitch100.individualOutputs.z_soln = Tswitch100.individualOutputs.z_soln + SIM.y_0_FOM';
            Tswitch100.combinedOutputs.z_soln   = Tswitch100.combinedOutputs.z_soln   + SIM.y_0_FOM';


            %% Plot data comparison
%             % Step Individual
%                 plotCompare(step.t_soln,       step.z_soln,       step.individualOutputs.t_soln,       step.individualOutputs.z_soln,       N.DesOut, RESULTS)
%             % Step Combined
%                 plotCompare(step.t_soln,       step.z_soln,       step.combinedOutputs.t_soln,         step.combinedOutputs.z_soln,         N.DesOut, RESULTS)
%             % Tswitch10 Individual
%                 plotCompare(Tswitch10.t_soln,  Tswitch10.z_soln,  Tswitch10.individualOutputs.t_soln,  Tswitch10.individualOutputs.z_soln,  N.DesOut, RESULTS)
%             % Tswitch10 Combined
%                 plotCompare(Tswitch10.t_soln,  Tswitch10.z_soln,  Tswitch10.combinedOutputs.t_soln,    Tswitch10.combinedOutputs.z_soln,    N.DesOut, RESULTS)
%             % Tswitch100 Individual
%                 plotCompare(Tswitch100.t_soln, Tswitch100.z_soln, Tswitch100.individualOutputs.t_soln, Tswitch100.individualOutputs.z_soln, N.DesOut, RESULTS)
%             % Tswitch100 Combined
%                 plotCompare(Tswitch100.t_soln, Tswitch100.z_soln, Tswitch100.combinedOutputs.t_soln,   Tswitch100.combinedOutputs.z_soln,   N.DesOut, RESULTS)

            %% Shift Data
            % Step Individual
                step.individualOutputs.t_soln = step.individualOutputs.t_soln(1:end-2);
                step.individualOutputs.z_soln = step.individualOutputs.z_soln(3:end,:);
            % Step Combined
                step.combinedOutputs.t_soln = step.combinedOutputs.t_soln(1:end-2);
                step.combinedOutputs.z_soln = step.combinedOutputs.z_soln(3:end,:);
            % Tswitch10 Individual
                Tswitch10.individualOutputs.t_soln = Tswitch10.individualOutputs.t_soln(1:end-2);
                Tswitch10.individualOutputs.z_soln = Tswitch10.individualOutputs.z_soln(3:end,:);
            % Tswitch10 Combined
                Tswitch10.combinedOutputs.t_soln = Tswitch10.combinedOutputs.t_soln(1:end-2);
                Tswitch10.combinedOutputs.z_soln = Tswitch10.combinedOutputs.z_soln(3:end,:);
            % Tswitch100 Individual
                Tswitch100.individualOutputs.t_soln = Tswitch100.individualOutputs.t_soln(1:end-2);
                Tswitch100.individualOutputs.z_soln = Tswitch100.individualOutputs.z_soln(3:end,:);
            % Tswitch100 Combined
                Tswitch100.combinedOutputs.t_soln = Tswitch100.combinedOutputs.t_soln(1:end-2);
                Tswitch100.combinedOutputs.z_soln = Tswitch100.combinedOutputs.z_soln(3:end,:);


            %% Plot Shifted data comparison
%             % Step Individual
%                 plotCompare(step.t_soln,       step.z_soln,       step.individualOutputs.t_soln,       step.individualOutputs.z_soln,       N.DesOut, RESULTS)
%             % Step Combined
%                 plotCompare(step.t_soln,       step.z_soln,       step.combinedOutputs.t_soln,         step.combinedOutputs.z_soln,         N.DesOut, RESULTS)
%             % Tswitch10 Individual
%                 plotCompare(Tswitch10.t_soln,  Tswitch10.z_soln,  Tswitch10.individualOutputs.t_soln,  Tswitch10.individualOutputs.z_soln,  N.DesOut, RESULTS)
%             % Tswitch10 Combined
%                 plotCompare(Tswitch10.t_soln,  Tswitch10.z_soln,  Tswitch10.combinedOutputs.t_soln,    Tswitch10.combinedOutputs.z_soln,    N.DesOut, RESULTS)
%             % Tswitch100 Individual
%                 plotCompare(Tswitch100.t_soln, Tswitch100.z_soln, Tswitch100.individualOutputs.t_soln, Tswitch100.individualOutputs.z_soln, N.DesOut, RESULTS)
%             % Tswitch100 Combined
%                 plotCompare(Tswitch100.t_soln, Tswitch100.z_soln, Tswitch100.combinedOutputs.t_soln,   Tswitch100.combinedOutputs.z_soln,   N.DesOut, RESULTS)


            %% Calculate Interpolated ODE Data
            for OO = 1:N.DesOut
                % Step Data
                    step.z_interpODE(OO,:) = interp1(step.t_soln, step.z_soln(OO,:) ,step.NewInputSignal(:,1));

                % Tswitch10
                    Tswitch10.z_interpODE(OO,:) = interp1(Tswitch10.t_soln, Tswitch10.z_soln(OO,:) ,Tswitch10.NewInputSignal(:,1));

                % Tswitch100
                    Tswitch100.z_interpODE(OO,:) = interp1(Tswitch100.t_soln, Tswitch100.z_soln(OO,:) ,Tswitch100.NewInputSignal(:,1));
            end


            %% Calculate Error at each time step
            % Step Individual
                step.individualOutputs.Error       = abs(step.individualOutputs.z_soln'       - step.z_interpODE(:,1:end-2));
            % Step Combined
                step.combinedOutputs.Error         = abs(step.combinedOutputs.z_soln'         - step.z_interpODE(:,1:end-2));
            % Tswitch10 Individual
                Tswitch10.individualOutputs.Error  = abs(Tswitch10.individualOutputs.z_soln'  - Tswitch10.z_interpODE(:,1:end-2));
            % Tswitch10 Combined
                Tswitch10.combinedOutputs.Error    = abs(Tswitch10.combinedOutputs.z_soln'    - Tswitch10.z_interpODE(:,1:end-2));
            % Tswitch100 Individual
                Tswitch100.individualOutputs.Error = abs(Tswitch100.individualOutputs.z_soln' - Tswitch100.z_interpODE(:,1:end-2));
            % Tswitch100 Combined
                Tswitch100.combinedOutputs.Error   = abs(Tswitch100.combinedOutputs.z_soln'   - Tswitch100.z_interpODE(:,1:end-2));


            %% Calculate Relative Error at each time step
            % Step Individual
                step.individualOutputs.RelError       = step.individualOutputs.Error       ./ step.z_interpODE(:,1:end-2);
            % Step Combined
                step.combinedOutputs.RelError         = step.combinedOutputs.Error         ./ step.z_interpODE(:,1:end-2);
            % Tswitch10 Individual
                Tswitch10.individualOutputs.RelError  = Tswitch10.individualOutputs.Error  ./ Tswitch10.z_interpODE(:,1:end-2);
            % Tswitch10 Combined
                Tswitch10.combinedOutputs.RelError    = Tswitch10.combinedOutputs.Error    ./ Tswitch10.z_interpODE(:,1:end-2);
            % Tswitch100 Individual
                Tswitch100.individualOutputs.RelError = Tswitch100.individualOutputs.Error ./ Tswitch100.z_interpODE(:,1:end-2);
            % Tswitch100 Combined
                Tswitch100.combinedOutputs.RelError   = Tswitch100.combinedOutputs.Error   ./ Tswitch100.z_interpODE(:,1:end-2);


            %% Calculate Average Relative Error for each desired output
            % Step Individual
                step.individualOutputs.AvgRelError       = (sum(step.individualOutputs.RelError,2)      / length(step.individualOutputs.RelError));
            % Step Combined
                step.combinedOutputs.AvgRelError         = (sum(step.combinedOutputs.RelError,2)        / length(step.combinedOutputs.RelError));
            % Tswitch10 Individual
                Tswitch10.individualOutputs.AvgRelError  = (sum(Tswitch10.individualOutputs.RelError,2) / length(Tswitch10.individualOutputs.RelError));
            % Tswitch10 Combined
                Tswitch10.combinedOutputs.AvgRelError    = (sum(Tswitch10.combinedOutputs.RelError,2)   / length(Tswitch10.combinedOutputs.RelError));
            % Tswitch100 Individual
                Tswitch100.individualOutputs.AvgRelError = (sum(Tswitch100.individualOutputs.RelError,2)/ length(Tswitch100.individualOutputs.RelError));
            % Tswitch100 Combined
                Tswitch100.combinedOutputs.AvgRelError   = (sum(Tswitch100.combinedOutputs.RelError,2)  / length(Tswitch100.combinedOutputs.RelError));


            %% Calculate Combined Average Relative Error for each desired output
            % Average across desired outputs (the extra entry)
            % Step Individual
                step.individualOutputs.CombinedAvgRelError       = (sum(step.individualOutputs.AvgRelError)      / length(step.individualOutputs.AvgRelError));
            % Step Combined
                step.combinedOutputs.CombinedAvgRelError         = (sum(step.combinedOutputs.AvgRelError)        / length(step.combinedOutputs.AvgRelError));
            % Tswitch10 Individual
                Tswitch10.individualOutputs.CombinedAvgRelError  = (sum(Tswitch10.individualOutputs.AvgRelError) / length(Tswitch10.individualOutputs.AvgRelError));
            % Tswitch10 Combined
                Tswitch10.combinedOutputs.CombinedAvgRelError    = (sum(Tswitch10.combinedOutputs.AvgRelError)   / length(Tswitch10.combinedOutputs.AvgRelError));
            % Tswitch100 Individual
                Tswitch100.individualOutputs.CombinedAvgRelError = (sum(Tswitch100.individualOutputs.AvgRelError)/ length(Tswitch100.individualOutputs.AvgRelError));
            % Tswitch100 Combined
                Tswitch100.combinedOutputs.CombinedAvgRelError   = (sum(Tswitch100.combinedOutputs.AvgRelError)  / length(Tswitch100.combinedOutputs.AvgRelError));


            %% Calculate Combined Average Relative Error
            % Each desired outputs, average across all sims
            combinedSims.individualOutputs.AvgRelError = (step.individualOutputs.AvgRelError + Tswitch10.individualOutputs.AvgRelError + Tswitch100.individualOutputs.AvgRelError)/3;
            combinedSims.combinedOutputs.AvgRelError   = (step.combinedOutputs.AvgRelError   + Tswitch10.combinedOutputs.AvgRelError   + Tswitch100.combinedOutputs.AvgRelError)  /3;


            %% Calculate Combined Combined Average Relative Error
            % Average across desired outputs (the extra entry)
            com_com_ind = sum(combinedSims.individualOutputs.AvgRelError) / length(combinedSims.individualOutputs.AvgRelError);
            com_com_com = sum(combinedSims.combinedOutputs.AvgRelError)   / length(combinedSims.combinedOutputs.AvgRelError);


            %% Fill Save Matricies
            for OO = 1:N.DesOut
                error.ind.step   (1 ,RR ,OO) = step.individualOutputs.AvgRelError(OO);
                error.ind.T10    (1 ,RR ,OO) = Tswitch10.individualOutputs.AvgRelError(OO);
                error.ind.T100   (1 ,RR ,OO) = Tswitch100.individualOutputs.AvgRelError(OO);
                error.ind.AllSims(1 ,RR ,OO) = combinedSims.individualOutputs.AvgRelError(OO);
                error.com.step   (1 ,RR ,OO) = step.combinedOutputs.AvgRelError(OO);
                error.com.T10    (1 ,RR ,OO) = Tswitch10.combinedOutputs.AvgRelError(OO);
                error.com.T100   (1 ,RR ,OO) = Tswitch100.combinedOutputs.AvgRelError(OO);
                error.com.AllSims(1 ,RR ,OO) = combinedSims.combinedOutputs.AvgRelError(OO);
            end
            %Average across outputs
            error.ind.step   (1 ,RR ,N.DesOut+1) = step.individualOutputs.CombinedAvgRelError;
            error.ind.T10    (1 ,RR ,N.DesOut+1) = Tswitch10.individualOutputs.CombinedAvgRelError;
            error.ind.T100   (1 ,RR ,N.DesOut+1) = Tswitch100.individualOutputs.CombinedAvgRelError;
            error.ind.AllSims(1 ,RR ,N.DesOut+1) = com_com_ind;
            error.com.step   (1 ,RR ,N.DesOut+1) = step.combinedOutputs.CombinedAvgRelError;
            error.com.T10    (1 ,RR ,N.DesOut+1) = Tswitch10.combinedOutputs.CombinedAvgRelError;
            error.com.T100   (1 ,RR ,N.DesOut+1) = Tswitch100.combinedOutputs.CombinedAvgRelError;
            error.com.AllSims(1 ,RR ,N.DesOut+1) = com_com_com;

        step = rmfield(step,'individualOutputs');
        step = rmfield(step,'combinedOutputs');
        Tswitch10 = rmfield(Tswitch10,'individualOutputs');
        Tswitch10 = rmfield(Tswitch10,'combinedOutputs');
        Tswitch100 = rmfield(Tswitch100,'individualOutputs');
        Tswitch100 = rmfield(Tswitch100,'combinedOutputs');
        end
        %% Find Optimal r
        % Each Desired Output
        for OO = 1:N.DesOut
            [~,optimal_r.ind.step(OO)]    = min(error.ind.step(1 ,: ,OO)   );
            [~,optimal_r.ind.T10(OO)]     = min(error.ind.T10(1 ,: ,OO)    );
            [~,optimal_r.ind.T100(OO)]    = min(error.ind.T100(1 ,: ,OO)   );
            [~,optimal_r.ind.AllSims(OO)] = min(error.ind.AllSims(1 ,: ,OO));
            [~,optimal_r.com.step(OO)]    = min(error.com.step(1 ,: ,OO)   );
            [~,optimal_r.com.T10(OO)]     = min(error.com.T10(1 ,: ,OO)    );
            [~,optimal_r.com.T100(OO)]    = min(error.com.T100(1 ,: ,OO)   );
            [~,optimal_r.com.AllSims(OO)] = min(error.com.AllSims(1 ,: ,OO));
        end

        % Across All Outputs
        [~,optimal_r.AllOutputs.ind.step]    = min(error.ind.step   (1 ,: ,N.DesOut+1));
        [~,optimal_r.AllOutputs.ind.T10]     = min(error.ind.T10    (1 ,: ,N.DesOut+1));
        [~,optimal_r.AllOutputs.ind.T100]    = min(error.ind.T100   (1 ,: ,N.DesOut+1));
        [~,optimal_r.AllOutputs.ind.AllSims] = min(error.ind.AllSims(1 ,: ,N.DesOut+1));
        [~,optimal_r.AllOutputs.com.step]    = min(error.com.step   (1 ,: ,N.DesOut+1));
        [~,optimal_r.AllOutputs.com.T10]     = min(error.com.T10    (1 ,: ,N.DesOut+1));
        [~,optimal_r.AllOutputs.com.T100]    = min(error.com.T100   (1 ,: ,N.DesOut+1));
        [~,optimal_r.AllOutputs.com.AllSims] = min(error.com.AllSims(1 ,: ,N.DesOut+1));


        %% Save Results for this SOC and Ts
        save(ROMError_filename,'singularValue','error','optimal_r','P','RESULTS')

    end
end


%% Plot Comparison Function
function plotCompare(ODE_T, ODE_Z, ROM_T, ROM_Z, DesOut,RESULTS)
for i = 1:DesOut
    figure
    hold on
    plot(ROM_T , ROM_Z(:,i),'ro','LineWidth',2,'DisplayName','ROM Ho-Kal')
    plot(ODE_T , ODE_Z(i,:),'k' ,'LineWidth',2,'DisplayName','ODE')
    lgn = legend;
    lgn.Location = 'best';
    title([RESULTS.Labels.title{i} ' No Noise Comparison'])
    xlabel('Time [s]')
    ylabel(RESULTS.Labels.unit{i})
end

end

