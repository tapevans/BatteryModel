%% Run CC Profile for Just simulink Model
clear all; close all; clc;


%% Load Initial Results
% Load ode solution to obtain tfinal
% ode_filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\SeminarPres_Nov2022\SeminarPres_Nov2022_Polar_0.05C_C.mat';
% ode_filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\SeminarPres_Nov2022\SeminarPres_Nov2022_Polar_0.20C_C.mat';
% ode_filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\SeminarPres_Nov2022\SeminarPres_Nov2022_Polar_0.50C_C.mat';
% ode_filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\SeminarPres_Nov2022\SeminarPres_Nov2022_Polar_1.00C_C.mat';
% ode_filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\SeminarPres_Nov2022\SeminarPres_Nov2022_Polar_2.00C_C.mat';
% ode_filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\SeminarPres_Nov2022\SeminarPres_Nov2022_Polar_5.00C_C.mat';

% ode_SIM = load(ode_filename);
% t_final = ode_SIM.t_soln(end);
% clear ode_SIM

% Load Simulink Init files
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\SeminarPres_Nov2022\SeminarPres_Nov2022_SimulinkChargeCC_Polar_0.05C_C.mat';
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\SeminarPres_Nov2022\SeminarPres_Nov2022_SimulinkChargeCC_Polar_0.20C_C.mat';
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\SeminarPres_Nov2022\SeminarPres_Nov2022_SimulinkChargeCC_Polar_0.50C_C.mat';
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\SeminarPres_Nov2022\SeminarPres_Nov2022_SimulinkChargeCC_Polar_1.00C_C.mat';
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\SeminarPres_Nov2022\SeminarPres_Nov2022_SimulinkChargeCC_Polar_2.00C_C.mat';
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\SeminarPres_Nov2022\SeminarPres_Nov2022_SimulinkChargeCC_Polar_5.00C_C.mat';

filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanTesting\P2DSystem\TestingPlant\KalmanTest_JustPlant_Polar_0.04C_D.mat'; % Jan 28,2023

load(filename)
t_final = 50;

%% Change Working Directory
    [current_file_path,~,~] = fileparts(mfilename('fullpath'));
    oldFolder = cd(current_file_path); % Change to the folder that contains the Simulink files. 
    % This is mostly so the cache files don't save to the main workspace


%% Calculate SS from SV_IC
i_user_sim = SIM.i_user_amp;
SV_IC = SIM.SV_IC;

[A,B,C,D] = getSSfromSV(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,SV_IC,i_user_sim);

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
        K = diag(diag_vec);


%% Output Matrix
% Pointers for Output Matrix
    i = 1;
    P.OM.cell_volt = i; i = i + 1;
    P.OM.del_phi   = i; i = i + 1;
%     P.OM.temp      = i; i = i + 1;
    P.OM.C_Liion   = i; i = i + 1;
    P.OM.X_surf    = i; i = i + 1;
    P.OM.i_Far     = i; i = i + 1;
    P.OM.eta       = i; i = i + 1;
    
    N.N_Out = length(fieldnames(P.OM));

    N.N_In  = 1;
        % I_user
    % Outputs
        % Cell Voltage
        % Delta Phi   @AN/SEP
        % Temperature @AN/SEP
        % C_Liion     @AN/SEP
        % X_surf      @AN/SEP
        % i_Far       @AN/SEP
        % eta         @AN/SEP
    OutputMatrix = zeros(N.N_Out , N.N_SV_tot);
        % Cell Voltage
            idx_phi_ed_AN = P.phi_ed;

            i = N.N_CV_CA(end);
            index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
            idx_phi_ed_CA = index_offset + P.phi_ed;

            OutputMatrix(P.OM.cell_volt,idx_phi_ed_AN) = -1;
            OutputMatrix(P.OM.cell_volt,idx_phi_ed_CA) =  1;
        % @AN/SEP
            i = N.N_CV_AN(end);
            index_offset = (i-1)*N.N_SV_AN;
        % Delta Phi   @AN/SEP
            idx = index_offset + P.del_phi;
            OutputMatrix(P.OM.del_phi,idx) =  1;
%         % Temperature @AN/SEP
%             idx = index_offset + P.T;
%             OutputMatrix(P.OM.temp,idx) = 1;
        % C_Liion     @AN/SEP
            idx = index_offset + P.C_Liion;
            OutputMatrix(P.OM.C_Liion,idx) = 1;
        % X_surf      @AN/SEP
            idx = index_offset + P.C_Li_surf_AN;
            OutputMatrix(P.OM.X_surf,idx) = 1/AN.C_Li_max;
        % i_Far      @AN/SEP
            idx = index_offset + P.i_PS;
            OutputMatrix(P.OM.i_Far,idx) = 1;
        % Eta      @AN/SEP
            idx = index_offset + P.V_2;
            OutputMatrix(P.OM.eta,idx) = 1;
            idx = index_offset + P.V_1;
            OutputMatrix(P.OM.eta,idx) = -1;


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
%     SIM = rmfield(SIM,'ControllerHandle');
%     SIM = rmfield(SIM,'Controller_MO_File');


%% Run Simulation
    % Create Model Object
        mdl = 'battery_system';
        load_system(mdl)

        in = Simulink.SimulationInput(mdl);
        in = in.setModelParameter('StartTime','0','StopTime',num2str(t_final));

    % Run Simulation
        tSimStart = tic;
        out = sim(in);
        SIM.tSimEnd = toc(tSimStart);


%% Post-Processing
    % Time
        t_soln = out.tout;
            
    % SV_soln
        n_row = length(t_soln);
        SV_soln = zeros(n_row,N.N_SV_tot);
        for i = 1:n_row
            SV_soln(i,:) = reshape(out.SV(:,:,i),[1,N.N_SV_tot]);
        end

%% Function Handle Conversion
    AN.EqPotentialHandle = str2func(AN.EqPotentialHandle);
    AN.i_oHandle         = str2func(AN.i_oHandle);
    AN.sigmaHandle       = str2func(AN.sigmaHandle);
    AN.D_oHandle         = str2func(AN.D_oHandle);
    
    CA.EqPotentialHandle = str2func(CA.EqPotentialHandle);
    CA.i_oHandle         = str2func(CA.i_oHandle);
    CA.sigmaHandle       = str2func(CA.sigmaHandle);
    CA.D_oHandle         = str2func(CA.D_oHandle);
    
    EL.tf_numHandle      = str2func(EL.tf_numHandle);
    EL.ActivityHandle    = str2func(EL.ActivityHandle);
    EL.D_o_Li_ionHandle  = str2func(EL.D_o_Li_ionHandle);
    EL.kappaHandle       = str2func(EL.kappaHandle);

%% Save Results
save(filename,'t_soln','SV_soln','AN','CA','SEP','EL','SIM','CONS','P','N','FLAG','PROPS','postProcessComplete')


%% Run other postProcessing
postProcessing(filename) % This saves the file again


%% Plot
% plotfcn(filename);
% 
% 
SIM.tSimEnd

