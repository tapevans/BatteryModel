%% RunSimulations
% This file will run all simulations in a given project folder

%% Change to this script's working directory
[filepath,~,~] = fileparts(mfilename('fullpath'));
cd(filepath)

%% List of Project Folders
i = 1;
Project_Folder{i} = 'zeroMeanPRBS_withRelax_Sims';   i = i+1;
% Project_Folder{i} = 'zeroMeanPRBS_Sims';   i = i+1;
% Project_Folder{i} = 'PRBS_Sims';   i = i+1;
% Project_Folder{i} = 'LongerImpulse';   i = i+1;
% Project_Folder{i} = 'KalmanTest';   i = i+1;
% Project_Folder{i} = 'ObservabilityTest';   i = i+1;
% Project_Folder{i} = 'SeminarPres_Nov2022';   i = i+1;
% Project_Folder{i} = 'TestingForTyrone';   i = i+1;
% Project_Folder{i} = 'Final_Lui_Wiley_Model';   i = i+1;

%%
num_Proj = length(Project_Folder);

%% Make a list of all sim file names with full path name
for i = 1:num_Proj
    % Make a list of simulations in the project folder
    oldFolder = cd([pwd filesep 'Results' filesep Project_Folder{i}]);
    list = dir('*.mat*');
    num_files = length(list);
    for j = 1:num_files % Creates a cell array with all simulations' full path name
        if ~exist('sim_filenames')
            sim_filenames{1} = [pwd filesep list(j).name];
        else
            sim_filenames{end+1,1} = [pwd filesep list(j).name];
        end
    end
    %Go back to oldFolder
    cd(oldFolder);
end
num_sim_files = length(sim_filenames);


%% Run all the simulations
for i = 1:num_sim_files
    disp(' ')
    disp(['Performing Simulation ' num2str(i) '/' num2str(num_sim_files)])
    disp(datestr(datetime));
    load(sim_filenames{i},'SIM')
    if ~(SIM.SimMode == 0)
        load(sim_filenames{i},'postProcessComplete')
    else
        postProcessComplete = 1;
    end
    
    %% Check if it has already ran
    if ~postProcessComplete
        load(sim_filenames{i})
        %% Run simulation
        tSimStart = tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ---- Polarization ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if SIM.SimMode == 1 
            % Simulation Parameters
            Tol.Abs = 1E-7;
            Tol.Rel = 1E-7;

            events = @(t,SV) batt_events(t,SV,SIM,P,N,FLAG);

            options = odeset('RelTol' ,Tol.Rel,      ...
                             'AbsTol' ,Tol.Abs,      ...
                             'Mass'   ,SIM.M,        ...
                             'Events' ,events);%,       ...
                            %'MaxStep',1e2);
            
            i_user = 0;
            for k = 1:length(SIM.tspan)-1
                if k == 1
                    tspan = [SIM.tspan(k), SIM.tspan(k+1)-1e-8];
                    SV_IC = SIM.SV_IC;
                    SOLN = ode15s(@(t,SV)batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user),tspan,SV_IC,options);
                    if FLAG.SaveSolnDiscreteTime
                        new_tfinal = SOLN.x(end);
                        save_time = (0:SIM.SaveTimeStep:new_tfinal)';
                        t_soln = save_time;
                        SV_soln = (deval(SOLN,save_time))';
                    else
                        t_soln  = SOLN.x';
                        SV_soln = SOLN.y';
                    end
                else
                    %New IC
                    SV_IC = SV_soln(end,:);
                    tspan = [SIM.tspan(k), SIM.tspan(k+1)-1e-8];
                    SOLN = ode15s(@(t,SV)batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user),tspan,SV_IC,options);
                    if FLAG.SaveSolnDiscreteTime
                        new_tfinal = SOLN.x(end);
                        save_time = (0:SIM.SaveTimeStep:new_tfinal)';
                        t_soln_int = save_time;
                        SV_soln_int = (deval(SOLN,save_time))';
                    else
                        t_soln_int  = SOLN.x';
                        SV_soln_int = SOLN.y';
                    end
                    t_soln  = [t_soln ; t_soln_int ];
                    SV_soln = [SV_soln; SV_soln_int];
                end
            end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ---- Harmonic Perturbation ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif SIM.SimMode == 2 
            % Simulation Parameters
            Tol.Abs = 1E-7;
            Tol.Rel = 1E-7;

            events = @(t,SV) batt_events(t,SV,SIM,P,N,FLAG);

            options = odeset('RelTol' ,Tol.Rel,      ...
                             'AbsTol' ,Tol.Abs,      ...
                             'Mass'   ,SIM.M,        ...
                             'Events' ,events);%,       ...
                            %'MaxStep',1e2);
            
            i_user = 0;
            tspan = SIM.tspan;
            SV_IC = SIM.SV_IC;
            SOLN = ode15s(@(t,SV)batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user),tspan,SV_IC,options);
            if FLAG.SaveSolnDiscreteTime
                new_tfinal = SOLN.x(end);
                save_time = (0:SIM.SaveTimeStep:new_tfinal)';
                t_soln = save_time;
                SV_soln = (deval(SOLN,save_time))';
            else
                t_soln  = SOLN.x';
                SV_soln = SOLN.y';
            end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ---- State Space EIS ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif SIM.SimMode == 3 
            [A,B,C,D,Z_results] = getSSImpedance(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS);
            %,sys_imp,sys_exp,sys_exp_mod
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ---- Known BC Profile Controller ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif SIM.SimMode == 4 
            % Initialize Save Variables
                t_soln      = [];
                SV_soln     = [];
                i_user_soln = [];
                mode_soln   = [];
                step_soln   = [];
                
            % Loop through all of the steps
            for k = 1:length(SIM.Controller_MO_File)
                SIM.current_MO_step = k;
                MO = SIM.Controller_MO_File(SIM.current_MO_step).MO;
                
                % Set SV_IC
                    if k == 1
                        SV_IC = SIM.SV_IC;
                    else
                        SV_IC = SV_soln(end,:)';
                    end
                
                % Determine i_user and tspan
                    if MO == 1 %CC
                        if SIM.Controller_MO_File(SIM.current_MO_step).CorD == 'C'
                            i_user = -SIM.Controller_MO_File(SIM.current_MO_step).C_rate*SIM.Cell_Cap/SIM.A_c;
                        else
                            i_user =  SIM.Controller_MO_File(SIM.current_MO_step).C_rate*SIM.Cell_Cap/SIM.A_c;
                        end
                        tfinal = SIM.Controller_MO_File(SIM.current_MO_step).Time_lim;
                        tspan = [0,tfinal];
                    elseif MO == 2 % CV
                        tspan_vec = [0,SIM.Controller_MO_File(SIM.current_MO_step).Time_lim];
                        if isempty(i_user_soln)
                            i_user = 0;
                        else
                            i_user = i_user_soln(end);
                        end
                    elseif MO == 3 % Relaxation
                        i_user = 0;
                        tfinal = SIM.Controller_MO_File(SIM.current_MO_step).Time_lim;
                        tspan = [0,tfinal];
                    end
                
                % Simulation Parameters
                % Needs to be inside loop since the voltage limits change with each type
                Tol.Abs = 1E-7;
                Tol.Rel = 1E-7;

                events = @(t,SV) batt_events(t,SV,SIM,P,N,FLAG);

                options_CC = odeset('RelTol' ,Tol.Rel,      ...
                                    'AbsTol' ,Tol.Abs,      ...
                                    'Mass'   ,SIM.M,        ...
                                    'Events' ,events);%,       ...
%                                     'MaxStep',1e0);%

                options_CV = odeset('RelTol' ,Tol.Rel,      ...
                                    'AbsTol' ,Tol.Abs,      ...
                                    'Mass'   ,SIM.M);

                % call ODE
                    if k~=1
                        clear SOLN
                    end
                    if MO == 1 %CC
                        SOLN = ode15s(@(t,SV)batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user),tspan,SV_IC,options_CC);
                    elseif MO == 2 % CV
                        [time_CV, SV_CV, current_CV, solver_CV] = runCVHold(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,SV_IC,i_user);
                        
                        SOLN.x = time_CV;
                        SOLN.y = SV_CV;
                        SOLN.c = current_CV;
                        SOLN.solver = solver_CV;

                    elseif MO == 3 % Relaxation
                        SOLN = ode15s(@(t,SV)batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user),tspan,SV_IC,options_CC);
                    end
                
                % Save solution
                    if MO ~= 2 % CC, Relax
                        if FLAG.SaveSolnDiscreteTime
                            new_tfinal  = SOLN.x(end);
                            save_time   = (0:SIM.SaveTimeStep:new_tfinal)';
                            t_soln_int  = save_time;
                            SV_soln_int = (deval(SOLN,save_time))';
                        else
                            t_soln_int  = SOLN.x';
                            SV_soln_int = SOLN.y';
                        end
                        i_user_soln_int = i_user*ones(length(t_soln_int),1);
                    else % CV
                        t_soln_int      = SOLN.x'; 
                        SV_soln_int     = SOLN.y';
                        i_user_soln_int = SOLN.c';
                    end
                    
                    if ~isempty(t_soln) % Adjust time values for continuing sims
                        t_soln_int = t_soln_int + t_soln(end);
                    end
                    
                    mode_soln_int = MO * ones(length(t_soln_int),1);
                    step_soln_int = k  * ones(length(t_soln_int),1);
                    
                    if k == 1
                        idx = 1;
                    else
                        idx = 2; % t^+ simulations
                        %idx = 1; % t^- simulations
                        %disp(idx)
                    end
                    
                    % t^+ simulations
                    t_soln      = [t_soln      ; t_soln_int(idx:end)     ];
                    SV_soln     = [SV_soln     ; SV_soln_int(idx:end,:)  ];
                    i_user_soln = [i_user_soln ; i_user_soln_int(idx:end)];
                    mode_soln   = [mode_soln   ; mode_soln_int(idx:end)  ];
                    step_soln   = [step_soln   ; step_soln_int(idx:end)  ];
                    
%                     % t^- simulations
%                     t_soln      = [t_soln(1:end-1,:)      ; t_soln_int(idx:end)     ];
%                     SV_soln     = [SV_soln(1:end-1,:)     ; SV_soln_int(idx:end,:)  ];
%                     i_user_soln = [i_user_soln(1:end-1,:) ; i_user_soln_int(idx:end)];
%                     mode_soln   = [mode_soln(1:end-1,:)   ; mode_soln_int(idx:end)  ];
%                     step_soln   = [step_soln(1:end-1,:)   ; step_soln_int(idx:end)  ];
            end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ---- MOO Controller ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif SIM.SimMode == 5
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ---- Manual Current Profile ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif SIM.SimMode == 7 
            % Simulation Parameters
            Tol.Abs = 1E-7;
            Tol.Rel = 1E-7;

            events = @(t,SV) batt_events(t,SV,SIM,P,N,FLAG);

            options = odeset('RelTol' ,Tol.Rel,      ...
                             'AbsTol' ,Tol.Abs,      ...
                             'Mass'   ,SIM.M,        ...
                             'Events' ,events);%,       ...
                            %'MaxStep',1e2);
            
            counter = 0;
            continue_Refinement = 1;
            while continue_Refinement
                counter = counter + 1;
                if counter == 1
                    disp(['Iteration ' num2str(counter)])
                else
                    disp(['Iteration ' num2str(counter) ,' sum_of_refinement_vec = ' num2str(SIM.sum_of_refinement_vec(counter - 1))])
                end
                %% Run a simulation
                i_user = 0;
                tspan = SIM.tspan;
                SV_IC = SIM.SV_IC;
                [t_soln,SV_soln,~,~,~] = ode15s(@(t,SV)batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user),tspan,SV_IC,options);

                %% Analyze the results
                r = max( N.N_SV_AN , N.N_SV_CA );
                N_t_steps = length(t_soln);
                SV      = zeros(r , N.N_CV_tot, N_t_steps );
                del_phi = zeros( N_t_steps , N.N_CV_tot );

                for j = 1:length(t_soln)
                    SV( : , : , j )  = SV1Dto2D( SV_soln( j , : ) , N , P , FLAG );
% !!!!!!!!!!!!!!!!                    SV( : , : , i ) = addPhiEl2SV(SV_temp,P,N);
                    del_phi( j , : ) = SV( P.del_phi , : , j);
                end

                %% Modify the charge profile
                if FLAG.Optimize_Profile
                    [continue_Refinement, SIM.profile_current , SIM.region_current_vec , sum_of_refinement] = plating_Refinement(del_phi, t_soln, SIM, N);
                    SIM.sum_of_refinement_vec(counter) = sum_of_refinement;
                else
                    continue_Refinement = 0;
                end

                %% Max number of iterations
                if counter == SIM.max_iterations
                    continue_Refinement = 0;
                end
                
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ---- PRBS ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif SIM.SimMode == 8
            % Simulation Parameters
            Tol.Abs = 1E-7;
            Tol.Rel = 1E-7;

            events = @(t,SV) batt_events(t,SV,SIM,P,N,FLAG);

            options = odeset('RelTol' ,Tol.Rel,      ...
                             'AbsTol' ,Tol.Abs,      ...
                             'Mass'   ,SIM.M);%,        ...
                             %'Events' ,events);%,       ...
                             %'MaxStep',1e2);
            
            i_user = 0;
            SOLN = ode15s(@(t,SV)batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user),SIM.tspan,SIM.SV_IC,options);

            if FLAG.SaveSolnDiscreteTime
                new_tfinal = SOLN.x(end);
                save_time = (0:SIM.SaveTimeStep:new_tfinal)';
                t_soln = save_time;
                SV_soln = (deval(SOLN,save_time))';
            else
                t_soln  = SOLN.x';
                SV_soln = SOLN.y';
            end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ---- Data Files ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif SIM.SimMode == 0 
            disp('This file is not a simulation file')
        else 
            % Maybe don't need else
            disp('Not a recognized simulation mode')
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SIM.tSimEnd = toc(tSimStart);
        
        
        %% Save results
        SIZE_SV_soln = whos('SV_soln');
        if SIM.SimMode == 0
            % Just a data file, don't save
        elseif SIM.SimMode == 3  % ---- State Space EIS ----
            postProcessComplete = 1;
            save(sim_filenames{i},'AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','A','B','C','D','Z_results','postProcessComplete') %,'sys_imp','sys_exp','sys_exp_mod'
            if FLAG.SaveSystemForEst % Save System to be used in Estimator
                multiple = -SIM.A_c^-1;
                sys = multiple*ss(A,B,C,D,'E',SIM.M);
                OutputAtEquil = SIM.OutputAtEquil;
                filename = [num2str(SIM.SOC_start) 'SOC'];
                save(filename,'sys','OutputAtEquil')
            end
        elseif SIM.SimMode == 4 % ---- Known BC Profile Controller ----
            save(sim_filenames{i},'t_soln','SV_soln','i_user_soln','mode_soln','step_soln','AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete','SOLN')
        elseif SIM.SimMode == 5 % ---- MOO Controller ----
            save(sim_filenames{i},'t_soln','SV_soln','i_user_soln','AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete')
        elseif SIM.SimMode == 8 % ---- PRBS ----
            if SIZE_SV_soln.bytes > 1e9 % 1 GB = 1e9 bytes
                save(sim_filenames{i},'t_soln','SV_soln','AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete','SOLN','-v7.3')
            else
                save(sim_filenames{i},'t_soln','SV_soln','AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete','SOLN')
            end
        else
            if SIZE_SV_soln.bytes > 1e9 % 1 GB = 1e9 bytes
                save(sim_filenames{i},'t_soln','SV_soln','AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete','-v7.3')
            else
                save(sim_filenames{i},'t_soln','SV_soln','AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete')
            end
        end
        
        %% Perform post-processing (Also re-saves data)
        if SIM.SimMode == 0 % ---- Data Files ----
            %Don't do anything to the data file
        else
            if SIM.SimMode ~= 3 % ---- State Space EIS ----
                if FLAG.doPostProcessing
                    disp('Performing Post-Processing')
                    postProcessing(sim_filenames{i});
                end
            end
        end
        
        %% Plot
        if exist('FLAG','var') && isfield(FLAG,'Plot')
            if FLAG.Plot
                plotfcn(sim_filenames{i});
            end
        end
    else
        disp('Simulation has already been analyzed')
    end
    clearvars -except sim_filenames i k num_sim_files
end
disp('Finished all simulations')
