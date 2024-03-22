%% RunSimulations
% This file will run all simulations in a given project folder

%% Change to this script's working directory
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    cd(filepath)

    
%% List of Project Folders
    i = 1;
    % Project_Folder{i} = 'TK_TestDeg';   i = i+1;
    % Project_Folder{i} = 'TK_DelV_DelT_Tests';   i = i+1;
    % Project_Folder{i} = 'TK_CyclingDegradation';   i = i+1;
    % Project_Folder{i} = 'TK_Test';   i = i+1;
    Project_Folder{i} = 'TK_GeneralComparison';   i = i+1;


%% Parameters
    num_Proj = length(Project_Folder);


%% Make a list of all sim file names with full path name
    for i = 1:num_Proj
        % Make a list of simulations in the project folder
            oldFolder = cd([pwd filesep 'Results' filesep Project_Folder{i}]);
            list = dir('*.mat*');
            num_files = length(list);
            sim_filenames = {};
            for j = 1:num_files % Creates a cell array with all simulations' full path name
                sim_filenames{end+1,1} = [pwd filesep list(j).name];
            end
    
        %Go back to oldFolder
            cd(oldFolder);
    end
    num_sim_files = length(sim_filenames);
    
    previousSimMode = 0; % Used to clear SimMode from persistent memory
    clear batt_GovEqn

%% Run all the simulations
for i = 1:num_sim_files
    disp(' ')
    disp(['Performing Simulation ' num2str(i) '/' num2str(num_sim_files)])
    disp(['Performing Simulation ' sim_filenames{i}])
    disp(datestr(datetime));
    simStart = tic;
    load(sim_filenames{i},'SIM')
    if ~(SIM.SimMode == 0)
        load(sim_filenames{i},'postProcessComplete')
    else
        postProcessComplete = 1;
    end
    
    %% Check if it has already ran
    if ~postProcessComplete
        clear batt_GovEqn
        load(sim_filenames{i})
        %% Run simulation
        tSimStart = tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ---- Polarization ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if SIM.SimMode == 1 
            % if SIM.SimMode ~= previousSimMode
            %     clear batt_GovEqn
            %     previousSimMode = SIM.SimMode;
            % end

            % Simulation Parameters
            Tol.Abs = 1E-7;
            Tol.Rel = 1E-7;

            events = @(t,SV) batt_events(t,SV,SIM,P,N,FLAG);

            options = odeset('RelTol' ,Tol.Rel,      ...
                             'AbsTol' ,Tol.Abs,      ...
                             'Mass'   ,SIM.M,        ...
                             'Events' ,events,       ...
                             'MaxStep',1e0);%);%
                if isfield(SIM,'JPattern')
                    options.JPattern = SIM.JPattern;
                end
            
            i_user = 0;
            tspan = [SIM.profile_time(1), SIM.profile_time(end)];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ---- Harmonic Perturbation ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif SIM.SimMode == 2 
            % if SIM.SimMode ~= previousSimMode
            %     clear batt_GovEqn
            %     previousSimMode = SIM.SimMode;
            % end

            % Simulation Parameters
            Tol.Abs = 1E-7;
            Tol.Rel = 1E-7;

            events = @(t,SV) batt_events(t,SV,SIM,P,N,FLAG);

            options = odeset('RelTol' ,Tol.Rel,      ...
                             'AbsTol' ,Tol.Abs,      ...
                             'Mass'   ,SIM.M,        ...
                             'Events' ,events);%,       ...
                             %'MaxStep',1e0);%
                if isfield(SIM,'JPattern')
                    options.JPattern = SIM.JPattern;
                end
            
            i_user = 0;
            tspan = [SIM.profile_time(1), SIM.profile_time(end)];
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
            % if SIM.SimMode ~= previousSimMode
            %     clear batt_GovEqn
            %     previousSimMode = SIM.SimMode;
            % end
            [A,B,C,D,Z_results] = getSSImpedance(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS);
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ---- Known BC Profile Controller ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif SIM.SimMode == 4 
            % if SIM.SimMode ~= previousSimMode
            %     clear batt_GovEqn
            %     previousSimMode = SIM.SimMode;
            % end
            % Initialize Save Variables
                t_soln      = [];
                SV_soln     = [];
                i_user_soln = [];
                mode_soln   = [];
                step_soln   = [];
                
            % Loop through all of the steps
            for k = 1:length(SIM.Controller_MO_File)
                clear batt_GovEqn
                k
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
                        SIM.profile_time    = SIM.profile_time_cell{SIM.current_MO_step};
                        SIM.profile_current = SIM.profile_current_cell{SIM.current_MO_step};
                        i_user = SIM.profile_current(1); % Just initializes the variable
                        tspan = [SIM.profile_time(1), SIM.profile_time(end)];
                    elseif MO == 2 % CV
                        % tspan_vec = [0,SIM.Controller_MO_File(SIM.current_MO_step).Time_lim];
                        if isempty(i_user_soln)
                            i_user = 0;
                        else
                            i_user = i_user_soln(end);
                        end
                        SIM.profile_time    = SIM.profile_time_cell{SIM.current_MO_step};
                        SIM.profile_current = SIM.profile_current_cell{SIM.current_MO_step};
                        tspan = [SIM.profile_time(1), SIM.profile_time(end)];
                    elseif MO == 3 % Relaxation
                        i_user = 0;
                        SIM.profile_time    = SIM.profile_time_cell{SIM.current_MO_step};
                        SIM.profile_current = SIM.profile_current_cell{SIM.current_MO_step};
                        tspan = [SIM.profile_time(1), SIM.profile_time(end)];
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
                                    % 'MaxStep',1e0);%
                    if isfield(SIM,'JPattern')
                        options.JPattern = SIM.JPattern;
                    end

                options_CV = odeset('RelTol' ,Tol.Rel,      ...
                                    'AbsTol' ,Tol.Abs,      ...
                                    'Mass'   ,SIM.M);
                    if isfield(SIM,'JPattern')
                        options.JPattern = SIM.JPattern;
                    end

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
            end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ---- MOO Controller ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif SIM.SimMode == 5
            % if SIM.SimMode ~= previousSimMode
            %     clear batt_GovEqn
            %     previousSimMode = SIM.SimMode;
            % end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ---- Manual Current Profile ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif SIM.SimMode == 7 
            % if SIM.SimMode ~= previousSimMode
            %     clear batt_GovEqn
            %     previousSimMode = SIM.SimMode;
            % end

            % Simulation Parameters
            Tol.Abs = 1E-7;
            Tol.Rel = 1E-7;

            events = @(t,SV) batt_events(t,SV,SIM,P,N,FLAG);

            options = odeset('RelTol' ,Tol.Rel,      ...
                             'AbsTol' ,Tol.Abs,      ...
                             'Mass'   ,SIM.M,        ...
                             'Events' ,events);%,       ...
                            %'MaxStep',1e2);
                if isfield(SIM,'JPattern')
                    options.JPattern = SIM.JPattern;
                end
            
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
                tspan = SIM.tspan; % This is created in CreateProject.m
                SV_IC = SIM.SV_IC;
                [t_soln,SV_soln,~,~,~] = ode15s(@(t,SV)batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user),tspan,SV_IC,options);

                %% Analyze the results
                r = max( N.N_SV_AN , N.N_SV_CA );
                N_t_steps = length(t_soln);
                SV      = zeros(r , N.N_CV_tot, N_t_steps );
                del_phi = zeros( N_t_steps , N.N_CV_tot );

                for j = 1:length(t_soln)
                    % SV( : , : , j ) = SV1Dto2D( SV_soln( j , : ) , N , P , FLAG );
                    SV( : , : , j ) = SV1Dto2D_short(SV_soln( j , : ), SIM.SV_nan, N.IDX_1Dto2D);
                    % SV( : , : , i ) = addPhiEl2SV(SV_temp,P,N); %!!!!!!!!!!!!!!!!
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
            % if SIM.SimMode ~= previousSimMode
            %     clear batt_GovEqn
            %     previousSimMode = SIM.SimMode;
            % end

            % Simulation Parameters
            Tol.Abs = 1E-7;
            Tol.Rel = 1E-7;

            events = @(t,SV) batt_events(t,SV,SIM,P,N,FLAG);

            options = odeset('RelTol' ,Tol.Rel,      ...
                             'AbsTol' ,Tol.Abs,      ...
                             'Mass'   ,SIM.M,        ...
                             'MaxStep',0.99*SIM.Tswitch);
                             %'Events' ,events);%,       ...
                if isfield(SIM,'JPattern')
                    options.JPattern = SIM.JPattern;
                end

            i_user = 0;
            
            % TestingSim = 1;
            % if TestingSim
                N_steps = length(SIM.profile_time_Noisy2D);
                % First Step
                    tspan  = SIM.profile_time_Noisy2D(1,:);
                    i_user = SIM.profile_current_Noisy2D(1,1);
                    SV_IC = SIM.SV_IC;
                    SOLN = ode15s(@(t,SV)batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user),tspan,SV_IC,options);
                    if FLAG.SaveSolnDiscreteTime
                        new_tfinal = SOLN.x(end);
                        save_time  = (0:SIM.SaveTimeStep:new_tfinal)';
                        t_soln     = save_time;
                        SV_soln    = (deval(SOLN,save_time))';
                    else
                        t_soln  = SOLN.x';
                        SV_soln = SOLN.y';
                    end
                        i_user_soln = i_user * ones(length(t_soln),1);

                % Other Steps
                for k = 2:N_steps
                    tspan  = SIM.profile_time_Noisy2D(k,:) - SIM.profile_time_Noisy2D(k,1);
                    i_user = SIM.profile_current_Noisy2D(k,1);
                    
                    % New IC
                        SV_IC = SV_soln(end,:);
                        SOLN = ode15s(@(t,SV)batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user),tspan,SV_IC,options);
                        if FLAG.SaveSolnDiscreteTime
                            new_tfinal      = SOLN.x(end);
                            save_time       = (0:SIM.SaveTimeStep:new_tfinal)';
                            t_soln_int      = save_time;
                            SV_soln_int     = (deval(SOLN,save_time))';
                        else
                            t_soln_int      = SOLN.x';
                            SV_soln_int     = SOLN.y';
                        end
                        i_user_soln_int = i_user * ones(length(t_soln_int),1);
                        t_soln_int = t_soln_int + SIM.profile_time_Noisy2D(k,1);
                        t_soln_int(1) = t_soln_int(1) + 1e-6 ;

                        t_soln      = [t_soln      ; t_soln_int ];
                        SV_soln     = [SV_soln     ; SV_soln_int];
                        i_user_soln = [i_user_soln ; i_user_soln_int];
                end
            % else
            %     %%%%%%%%%% Testing splitting sim %%%%%%%%%%
            %     if SIM.UseODEextend
            %         split_time = linspace(0 , SIM.profile_time(end) , SIM.NumCuts+1);
            % 
            %         % Save the initial parameters
            %         if SIM.SaveIntermediate
            %             SIM.SimMode = 0;
            %             postProcessComplete = 1;
            %             save([sim_filenames{i}(1:end-4) '_InputParams.mat'],'AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete')
            %             postProcessComplete = 0;
            %             SIM.SimMode = 8;
            %         end
            % 
            %         % Run First Cut
            %         SOLN = ode15s(@(t,SV)batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user),split_time(1,1:2) , SIM.SV_IC , options);
            %         SOLN_Next = SOLN;
            % 
            %         % Save First Cut
            %         if SIM.SaveIntermediate
            %             SIM.SimMode = 0;
            %             postProcessComplete = 1;
            %             save([sim_filenames{i}(1:end-4) '_Cut' num2str(1) '.mat'],'SOLN_Next','postProcessComplete','SIM')
            %             postProcessComplete = 0;
            %             SIM.SimMode = 8;
            %         end
            % 
            %         % Loop through the rest of the Cuts
            %         for j = 2:SIM.NumCuts
            %             SOLN_temp = SOLN_Next;
            %             if SIM.REDUCESOLN
            %                 n_points = length(SOLN_temp.x);
            %                 y_IDX    = [n_points-1:n_points];
            %                 SOLN_temp.x = SOLN_Next.x(1,y_IDX);
            %                 SOLN_temp.y = SOLN_Next.y(:,y_IDX);
            %                 SOLN_temp.idata.kvec = SOLN_Next.idata.kvec(1,y_IDX);
            %                 SOLN_temp.idata.dif3d = SOLN_Next.idata.dif3d(:,:,y_IDX);
            %             end
            %             SOLN_Next = odextend(SOLN_temp , @(t,SV)batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user), split_time(1,j+1));
            %             % SOLN_Next = rmfield( SOLN_Next , 'idata' );
            % 
            %             % Save Cuts
            %             if SIM.SaveIntermediate
            %                 SIM.SimMode = 0;
            %                 postProcessComplete = 1;
            %                 save([sim_filenames{i}(1:end-4) '_Cut' num2str(j) '.mat'],'SOLN_Next','postProcessComplete','SIM')
            %                 postProcessComplete = 0;
            %                 SIM.SimMode = 8;
            %             end
            % 
            %             % Combine SOLN
            %             if SIM.REDUCESOLN
            %                 IDX = find(SOLN_Next.x> SOLN.x(end),1);
            %                 SOLN.x = [SOLN.x , SOLN_Next.x(1,IDX:end) ];
            %                 SOLN.y = [SOLN.y , SOLN_Next.y(:,IDX:end) ];
            %             else
            %                 SOLN = SOLN_Next;
            %             end
            %         end
            %     else
            %         SIM.tspan = [0, SIM.profile_time(end)];
            %         SOLN = ode15s(@(t,SV)batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user),SIM.tspan,SIM.SV_IC,options);
            %     end
            % 
            %     % %%%%% When having to run smaller time steps %%%%%%%%%%%Testing!!!!!!!!!!!!!!!!!!!!!!!!!!
            %     %     old_tFinal = SIM.profile_time(end);
            %     %     t_intermed = 0.99*SIM.Tswitch;
            %     % 
            %     %     %Run first step
            %     %         options = odeset('RelTol' ,Tol.Rel,      ...
            %     %                          'AbsTol' ,Tol.Abs,      ...
            %     %                          'Mass'   ,SIM.M,        ...
            %     %                          'MaxStep',0.01*SIM.Tswitch);
            %     %                          %'Events' ,events);%,       ...
            %     % 
            %     %         SIM.tspan = [0,t_intermed]; 
            %     %         i_user = 0;
            %     %         SOLN = ode15s(@(t,SV)batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user),SIM.tspan,SIM.SV_IC,options);
            %     % 
            %     %     %Run the rest
            %     %         IC_new = SOLN.y(:,end);
            %     %         options = odeset('RelTol' ,Tol.Rel,      ...
            %     %                          'AbsTol' ,Tol.Abs,      ...
            %     %                          'Mass'   ,SIM.M,        ...
            %     %                          'MaxStep',0.99*SIM.Tswitch);
            %     %                          %'Events' ,events);%,       ...
            %     % 
            %     %         % SIM.tspan = [0,t_intermed]; 
            %     %         % i_user = 0;
            %     % 
            %     %         SOLN = odextend(SOLN, @(t,SV)batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user) , old_tFinal , IC_new , options);
            % 
            %     if FLAG.SaveSolnDiscreteTime
            %         new_tfinal = SOLN.x(end);
            %         save_time = (0:SIM.SaveTimeStep:new_tfinal)';
            %         t_soln = save_time;
            %         SV_soln = (deval(SOLN,save_time))';
            %     else
            %         t_soln  = SOLN.x';
            %         SV_soln = SOLN.y';
            %     end
            % end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ---- EIS from Stitching PRBS ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif SIM.SimMode == 9
            % if SIM.SimMode ~= previousSimMode
            %     clear batt_GovEqn
            %     previousSimMode = SIM.SimMode;
            % end

            % Load in PRBS results
                NumPRBS = length(PRBS_desTswitch_filenames);
                % Results = struct(NumPRBS , 1);
                for j = 1:NumPRBS
                    data = load(PRBS_desTswitch_filenames{j});
                    Results(j).t            = data.t_soln;
                    Results(j).cell_voltage = data.cell_voltage;
                    Results(j).inputCurrent = data.SIM.A_c * interp1(data.SIM.profile_time , data.SIM.profile_current, data.t_soln); % (data.SIM.A_c*) converts from A/m^2 to A
                    Results(j).Tswitch      = data.SIM.Tswitch;
                end
                PRBSLength = data.SIM.DesiredLength;

            % Identify initial state-space models
                % Inititalize
                    dofilter = 0; %%%%%%%%%Hardcoded
                    sys                   = cell(NumPRBS , 1);
                    info_all              = cell(NumPRBS , 1);
                    smoothdat             = cell(NumPRBS , 1);
                    Update_Freq_all       = cell(NumPRBS , 1);
                    Ts_all                = cell(NumPRBS , 1);

                for j = 1:NumPRBS
                   Ts = Results(j).t(2)-Results(j).t(1);
                   freq = 1/Ts;
                   [sysc_out , sys_out , dat_out , structuredat_out]  = IdentificationProcess(  Results(j).cell_voltage,-Results(j).inputCurrent, Results(j).t, freq, dofilter);
                   sys{j}                   = sys_out;
                   info_all{j}.delayCount   = 1;
                   smoothdat{j}             = dat_out;
                   Update_Freq_all{j}       = freq;
                   Ts_all{j}                = Ts;

                   % Results(j).sys                 = sys_out;
                   % Results(j).info_all.delayCount = 1;
                   % Results(j).smoothdat           = dat_out;
                   % Results(j).Update_Freq_all     = freq;
                   % Results(j).Ts_all              = Ts;
                end

            % Stitch togehter state-space models
                Tswitch_vec = nan(NumPRBS,1);
                for j = 1:NumPRBS
                    Tswitch_vec(j) = Results(j).Tswitch;
                end
                Tswitch_min = min(Tswitch_vec);
                Tswitch_max = max(Tswitch_vec);

                wtot = logspace(log10((pi/(Tswitch_min/100))) , log10((2*pi/(PRBSLength*100*Tswitch_max))) , 50)';% rad/s
                % if isempty(EIS_PRBS.freq)
                %     wtot = logspace(log10((pi/(Tswitch_min/100))) , log10((2*pi/(PRBSLength*100*Tswitch_max))) , 50)';% rad/s
                % else
                %     wtot = EIS_PRBS.freq;
                % end

                [re_est,re_est_var,im_est,im_est_var,re,im,var_re,var_im] = FrequencyEstimateProcess(sys,wtot,Update_Freq_all,PRBSLength);
                % EISresultsdisplay

            % fit of all data to continuous time model (G)
                Iterations   = 4E+3;   % 4E+3, 500
                [xout,np,nz] = CTSysID_EIS_Estimate(sys,info_all,smoothdat,Iterations);
                [G,~]        = pzderiv_tot(xout,np,nz,smoothdat,sys);
 
            % Extract EIS from continuous time model
                % [re_est,re_est_var,im_est,im_est_var,nomsys] = FrequencyEstimate_tot(xout,np,nz,P,wtot);
               
            % Plot of continuous time results
                % refineresultsdisplay

            % Get Impedance
                if ~isempty(EIS_PRBS.freq)
                    [mag,phase,~] = bode(G,EIS_PRBS.freq);
                    [re ,im   ,~] = nyquist(G,EIS_PRBS.freq);
                    mag   = reshape(mag  ,[],1);
                    phase = reshape(phase,[],1);
                    re    = reshape(re   ,[],1);
                    im    = reshape(im   ,[],1);
                    Z_results(:,P.SS.omega)    = EIS_PRBS.freq';
                    Z_results(:,P.SS.Z_mag)    = mag;
                    Z_results(:,P.SS.Z_Re)     = re;
                    Z_results(:,P.SS.Z_Im)     = im;
                    Z_results(:,P.SS.Z_dB)     = 20*log10(mag);
                    Z_results(:,P.SS.Z_ps_deg) = phase;
                else
                    Z_results = [];
                end

                postProcessComplete = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ---- EIS Ho-Kalman ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif SIM.SimMode == 10
            % if SIM.SimMode ~= previousSimMode
            %     clear batt_GovEqn
            %     previousSimMode = SIM.SimMode;
            % end
            
        % Get impulse response
            % Simulation Parameters
            Tol.Abs = 1E-7;
            Tol.Rel = 1E-7;

            events = @(t,SV) batt_events(t,SV,SIM,P,N,FLAG);

            options = odeset('RelTol' ,Tol.Rel,      ...
                             'AbsTol' ,Tol.Abs,      ...
                             'Mass'   ,SIM.M,        ...
                             'Events' ,events,       ...
                             'MaxStep',0.5*SIM.Tsample);%
                if isfield(SIM,'JPattern')
                    options.JPattern = SIM.JPattern;
                end
            
            i_user = 0;
            tspan = [0, SIM.profile_time(end)];
            SV_IC = SIM.SV_IC;
            SOLN_CT = ode15s(@(t,SV)batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user),tspan,SV_IC,options);
            if FLAG.SaveSolnDiscreteTime
                new_tfinal = SOLN_CT.x(end);
                save_time = (0:SIM.SaveTimeStep:new_tfinal)';
                t_soln_CT = save_time;
                SV_soln_CT = (deval(SOLN_CT,save_time))';
            else
                t_soln_CT  = SOLN_CT.x';
                SV_soln_CT = SOLN_CT.y';
            end

            % DT Solution
                SOLN       = SOLN_CT;
                new_tfinal = SOLN.x(end);
                save_time  = (0:SIM.SaveTimeStep:new_tfinal)';
                t_soln     = save_time;
                SV_soln    = (deval(SOLN,save_time))';


        % Some how reduce the data to just cell voltage response
            SIM.TsMultiple = 5;
            t_test = 0 : (SIM.Tsample/SIM.TsMultiple) : SIM.profile_time(end);
            i_user = i_user_calc( t_test , SIM);
            z_imp = SIM.OutputMatrix * SV_soln';
            z_imp = z_imp(P.OM.cell_volt,:);

            freq    = SIM.freq;

        % Get full-order state-space from Ho-Kalman
            r = 355;
            [sys_FOM , ~] = getHoKalman(t_soln , i_user , z_imp , r , SIM , FLAG);

        % Get full-order impedence at desired frequencies
            M = eye(size(sys_FOM.A));
            [Z_results_FOM] = getImpedanceFromSSSystem(sys_FOM , M , freq , SIM , P);

        % Get reduced-order state-space from Ho-Kalman
            r = 10;
            [sys_ROM , ~] = getHoKalman(t_soln , i_user , z_imp , r , SIM , FLAG);

        % Get reduced-order impedence at desired frequencies
            M = eye(size(sys_ROM.A));
            [Z_results_ROM] = getImpedanceFromSSSystem(sys_ROM , M , freq , SIM , P);


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

        % ---- Data Files ----
        if SIM.SimMode == 0
            % Just a data file, don't save

        % ---- State Space EIS ----
        elseif SIM.SimMode == 3  
            postProcessComplete = 1;
            save(sim_filenames{i},'AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','A','B','C','D','Z_results','postProcessComplete') %,'sys_imp','sys_exp','sys_exp_mod'
            if FLAG.SaveSystemForEst % Save System to be used in Estimator
                multiple = -SIM.A_c^-1;
                sys = multiple*ss(A,B,C,D,'E',SIM.M);
                OutputAtEquil = SIM.OutputAtEquil;
                filename = [num2str(SIM.SOC_start) 'SOC'];
                save(filename,'sys','OutputAtEquil')
            end

        % ---- Known BC Profile Controller ----
        elseif SIM.SimMode == 4 
            if SIZE_SV_soln.bytes > 1e9
                save(sim_filenames{i},'t_soln','SV_soln','i_user_soln','mode_soln','step_soln','AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete','SOLN','-v7.3')
            else
                save(sim_filenames{i},'t_soln','SV_soln','i_user_soln','mode_soln','step_soln','AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete','SOLN')
            end

        % ---- MOO Controller ----
        elseif SIM.SimMode == 5 
            save(sim_filenames{i},'t_soln','SV_soln','i_user_soln','AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete')

        % ---- PRBS ----
        elseif SIM.SimMode == 8 
            if SIZE_SV_soln.bytes > 1e9 % 1 GB = 1e9 bytes
                save(sim_filenames{i},'t_soln','SV_soln','i_user_soln','AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete','SOLN','-v7.3')
            else
                save(sim_filenames{i},'t_soln','SV_soln','i_user_soln','AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete','SOLN')
            end

        % ---- EIS PRBS ----
        elseif SIM.SimMode == 9
            save(sim_filenames{i},'G','Z_results','EIS_PRBS','SIM','P','PRBS_desTswitch_filenames','Results','postProcessComplete')
            % if SIZE_SV_soln.bytes > 1e9 % 1 GB = 1e9 bytes
            %     save(sim_filenames{i},'G','Z_results','EIS_PRBS','SIM','PRBS_desTswitch_filenames','Results','postProcessComplete','SOLN','-v7.3')
            % else
            %     save(sim_filenames{i},'G','Z_results','EIS_PRBS','SIM','PRBS_desTswitch_filenames','Results','postProcessComplete','SOLN')
            % end

        % ---- EIS Ho-Kalman ----
        elseif SIM.SimMode == 10
            % postProcessComplete = 1;
            save(sim_filenames{i},'t_soln','SV_soln','t_soln_CT','SV_soln_CT','AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','SOLN','SOLN_CT','sys_FOM','sys_ROM','Z_results_FOM','Z_results_ROM','postProcessComplete')
            
        else
            if SIZE_SV_soln.bytes > 1e9 % 1 GB = 1e9 bytes
                save(sim_filenames{i},'t_soln','SV_soln','SOLN','AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete','-v7.3')
            else
                save(sim_filenames{i},'t_soln','SV_soln','SOLN','AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete')
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
    simEnd = toc(simStart);
    disp(['Simulation took:' num2str(simEnd) 's'])
    clearvars -except sim_filenames i k num_sim_files previousSimMode
end
disp('Finished all simulations')

