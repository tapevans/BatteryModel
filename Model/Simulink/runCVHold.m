function [time_CV, SV_CV, current_CV, solver_CV] = runCVHold(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,SV_IC,i_user_IC)
%% runCVHold
% runCVHold is a part of the KBCP simulation mode that handles the constant
% voltage case.
%
% This function performs the constant voltage profile by creating a PI
% controller and feeding the parameters into a Simulink model specifically
% made for CV mode.
%
% Inputs:
%  * Battery parameter structs
%  * SV_IC: State vector typically taken from the end of the constant
%           current cycle
%  * i_user_IC: [A/m^2], Load current flux used during the previous cycle
%
% Outputs:
%  * time_CV:    [s],     Simulation time response
%  * SV_soln:       ,     Simulation state vector
%  * current_CV: [A/m^2], Simulation current flux response


%% Constants
    t_r = 0.1; w_c = 2/t_r; % Rise time [s]; Cutoff frequency [1/s]
    
    t_final_cutoff = 1000; % [s], Time that is reported back to the solver if the initial current flux is less than C/20


%% Change Working Directory
    [current_file_path,~,~] = fileparts(mfilename('fullpath'));
    oldFolder = cd(current_file_path); % Change to the folder that contains the Simulink files. This is mostly so the cache files don't save to the main workspace


%% Parse Incoming Data
    % Run Time
        t_final = SIM.Controller_MO_File(SIM.current_MO_step).Time_lim;
    
    % Reference Voltage (Set Point)
        CV_ref  = SIM.Controller_MO_File(SIM.current_MO_step).Volt_ref;


%% Create Plant for Tuning PID
    [A,B,C,D] = getSSfromSV(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,SV_IC,i_user_IC);
    C = C(1,:); 
    D = 0;


%% Tune PID
    plant = ss(A,B,C,D,'E',SIM.M);
    
    Controller = pidtune(plant, 'PI', w_c);
    % Extract controller gains
        K_P = Controller.Kp;
        K_I = Controller.Ki;
        K_D = 0;


%% Calculate Parameters for Battery System
    % Pseudo-inverse of Mass
        [U, S, V] = svd(SIM.M);
        threshold = 1e-7;
        S_cross = pinv(S,threshold);
        M_cross = V*S_cross*U';
    
    % Calculate differential variable dynamics
        A_cross = M_cross*A;
        [~ , ol_poles] = eig(M_cross*A);
        fastest_diff_pole = min(diag(real(ol_poles)));
    
    % Calculate Null Space matricies
        r = rank(S_cross);
        
        U_null = U(:,r+1:end);
        V_null = V(:,r+1:end);
        U_NV_NT = U_null*V_null';
    
    % Calculate algebraic poles
        LC_k = 10;
        diag_vec = zeros(N.N_SV_tot,1);
        algb_poles = LC_k*fastest_diff_pole*ones(length(SIM.algb_idx),1);
        diag_vec(SIM.algb_idx) = algb_poles;
        K_algb = diag(diag_vec);


%% Function Handle Conversion
    AN.EqPotentialHandle = func2str(AN.EqPotentialHandle);
    AN.i_oHandle         = func2str(AN.i_oHandle);
    AN.sigmaHandle       = func2str(AN.sigmaHandle);
    AN.D_oHandle         = func2str(AN.D_oHandle);
    
    CA.EqPotentialHandle = func2str(CA.EqPotentialHandle);
    CA.i_oHandle         = func2str(CA.i_oHandle);
    CA.sigmaHandle       = func2str(CA.sigmaHandle);
    CA.D_oHandle         = func2str(CA.D_oHandle);
    
    EL.tf_numHandle      = func2str(EL.tf_numHandle);
    EL.ActivityHandle    = func2str(EL.ActivityHandle);
    EL.D_o_Li_ionHandle  = func2str(EL.D_o_Li_ionHandle);
    EL.kappaHandle       = func2str(EL.kappaHandle);
    
    SIM = rmfield(SIM,'fsolve_options');
    SIM = rmfield(SIM,'ControllerHandle');
    SIM = rmfield(SIM,'Controller_MO_File');


%% Run Simulation
    % Create Model Object
        mdl = 'battery_system_CV_PID';
        load_system(mdl)

        in = Simulink.SimulationInput(mdl);
        in = in.setModelParameter('StartTime','0','StopTime',num2str(t_final));
    
    % Assign variable values by modifying the workspace
        mdlWks = get_param(in,'ModelWorkspace');
        assignin(mdlWks,'CV_ref'    ,CV_ref)
        assignin(mdlWks,'C'         ,C)
        assignin(mdlWks,'i_user_IC' ,i_user_IC)
        assignin(mdlWks,'K_P'       ,K_P)
        assignin(mdlWks,'K_I'       ,K_I)
        assignin(mdlWks,'K_D'       ,K_D)
        assignin(mdlWks,'M_cross'   ,M_cross)
        assignin(mdlWks,'U_NV_NT'   ,U_NV_NT)
        assignin(mdlWks,'K_algb'    ,K_algb)
        assignin(mdlWks,'SV_IC'     ,SV_IC)
        assignin(mdlWks,'AN'        ,AN)
        assignin(mdlWks,'CA'        ,CA)
        assignin(mdlWks,'SEP'       ,SEP)
        assignin(mdlWks,'EL'        ,EL)
        assignin(mdlWks,'SIM'       ,SIM)
        assignin(mdlWks,'CONS'      ,CONS)
        assignin(mdlWks,'P'         ,P)
        assignin(mdlWks,'N'         ,N)
        assignin(mdlWks,'FLAG'      ,FLAG)
        assignin(mdlWks,'PROPS'     ,PROPS)
        if FLAG.SaveSolnDiscreteTime
            SampTime = SIM.SaveTimeStep;
        else
            SampTime = -1;
        end
        assignin(mdlWks,'SampTime'  ,SampTime)
        

    % Run Simulation
        out = sim(in);
    

%% Post-Processing
    % Time
        time_CV = out.tout_mine';

            
    % SV_soln
        SV_CV = out.SV';
    
    % Current Flux
        out.i_user = reshape(out.i_user,1,[]);
        current_CV = out.i_user;

    % Solver
        solver_CV = out.SimulationMetadata.ModelInfo.SolverInfo.Solver;

%% Remove Data below C/20
    C_rate = 1/20;
    current_limit = C_rate*SIM.Cell_Cap/SIM.A_c;
    if abs(i_user_IC) <= current_limit % If the CC mode before was less than C/20
        idx = find(time_CV > t_final_cutoff,1);
        if ~isempty(idx)
            time_CV    = time_CV(1:idx);
            SV_CV      = SV_CV(:,1:idx);
            current_CV = current_CV(1:idx);
        end
    else % End condition of the CV mode
        idx = find(abs(current_CV) < current_limit,1);
        if ~isempty(idx)
            time_CV    = time_CV(1:idx);
            SV_CV      = SV_CV(:,1:idx);
            current_CV = current_CV(1:idx);
        end
    end


%% Return to Old Working Directory
    cd(oldFolder);


end