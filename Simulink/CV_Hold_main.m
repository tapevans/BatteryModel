%% Test Battery In Simulink
clear all; close all; clc;

%% Parameters
    p_OS = 1e-6; % [%], Percent Overshoot
    t_s  = 2e-4; % [s], Settling Time
%     tfinal = 30*60; % [s], Final Simulation Time
    tfinal = 1; % [s], Final Simulation Time
    
    model_filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\TestingSimulink\TestingSimulink_KPCont_SingleStepTo50SOCSOC0.mat';
    equil_filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\TestingSimulink\TestingSimulink_SS_EIS_SOC50.mat';
    save_filename  = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\TestingSimulink\SimulinkResults.mat';

    FLAG.PostProcess = 0;
    FLAG.PLOT = 0;

%% Update Working Directory
    [current_file_path,~,~] = fileparts(mfilename('fullpath'));
    cd(current_file_path)

% Add Model Folder
    idx = find(current_file_path == '\',1,'last');
    model_folder = [ current_file_path(1:idx) 'Model'];
    addpath(genpath(model_folder));

% Remove Functionality Test Folder
    rmpath('F:\TylerFiles\GitHubRepos\BatteryModel\Model\FunctionalityTests')

%% Load Model File
model = load(model_filename);
SV_IC = model.SV_soln(end,:)';
clear model;

%% Load Equilibrium State File
load(equil_filename);

% %% Get SS Matricies
% [A,B,C,D,~] = getSSImpedance(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS);
C = C(1,:); % Cell_Voltage is the first row
D = 0;

%% Calculate Algebraic Control Value
% Pseudo-inverse of Mass
    [U, S, V] = svd(SIM.M);
    threshold = 1e-7;
    S_cross = pinv(S,threshold);
    M_cross = V*S_cross*U';

% Calculate differential variable dynamics
    A_cross = M_cross*A;
    [ol_vec , ol_poles] = eig(M_cross*A);
    fastest_diff_pole = min(diag(real(ol_poles)));
%%%     ol_poles_diag = diag(real(ol_poles));

% Calculate Null Space matricies
    r = rank(S_cross); 
    U_colm = U(:,1:r);
    V_colm = V(:,1:r);
    U_null = U(:,r+1:end);
    V_null = V(:,r+1:end);
    
    U_NV_NT = U_null*V_null';

% Calculate algebraic poles
    LC_k = 10;
    diag_vec = zeros(N.N_SV_tot,1);
    algb_poles = LC_k*fastest_diff_pole*ones(length(SIM.algb_idx),1);
    diag_vec(SIM.algb_idx) = algb_poles;
    K_algb = diag(diag_vec);

% Calculate New System A,B
    A_comb = M_cross*A - K_algb*U_NV_NT*A;
%     B_comb = M_cross*B + U_NV_NT*B;
    B_comb = M_cross*B;
    
%     [sys_eig_vec ,sys_poles] = eig(A_comb);
%     sys_poles_diag = diag(real(sys_poles));
%     find(sys_poles_diag>0)
% %     find(sys_poles_diag>0)
%     reduced_sys_eig_vec = real(sys_eig_vec);
%     for i = 1:length(reduced_sys_eig_vec)
%         for j = 1:length(reduced_sys_eig_vec)
%             if abs(reduced_sys_eig_vec(i,j))<1e-5
%                 reduced_sys_eig_vec(i,j) = 0;
%             else
%                 reduced_sys_eig_vec(i,j) = reduced_sys_eig_vec(i,j);
%             end
%         end
%     end

%% Calculate Closed Loop Pole Dynamics
    zeta = sqrt( ((log(p_OS/100))^2)/ ( pi^2 + (log(p_OS/100))^2  )   );
    omega_n = 4.6/( zeta*t_s );

    pole_p = -omega_n*(zeta + 1i*sqrt(1-zeta^2));
    pole_m = -omega_n*(zeta - 1i*sqrt(1-zeta^2));
    fast_pole = -5*zeta*omega_n;
%     fast_pole_vec = linspace(0,1,(N.N_SV_tot-2))';
    fast_pole_vec = -(0:1:N.N_SV_tot-3)';
    fast_pole_vec = fast_pole_vec + fast_pole;

    cl_poles = [pole_p ; pole_m ; fast_pole_vec];

%% Calculate Controller Matricies
% % Determine the Controller Values
%     K = place(A_comb,B_comb,cl_poles);
    K = place(A,B,cl_poles);
    Co = ctrb(A,B);
    rank(Co);
    unco = length(A_comb) - rank(Co);

% Determine N_x and N_u (See notes or control textbook)
    zero_vec = zeros(N.N_SV_tot,1);
%     a = [A_comb , B_comb; C , D];
    a = [A , B; C , D];
    b = [zeros(N.N_SV_tot,1) ;1];
    N_vec = a\b;

    N_x = N_vec(1:N.N_SV_tot,1);
    N_u = N_vec(N.N_SV_tot+1,1);


%% Other Constants
CV = C*SV_IC;

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
% SIM = rmfield(SIM,'ControllerHandle');

%% Run Simulation
mdl = 'battery_system_CV';
in = Simulink.SimulationInput(mdl);
in = in.setModelParameter('StartTime','0','StopTime',num2str(tfinal));

tStart = tic;
% out = sim(in);
t_solve = toc(tStart);


%% Post-Processing
if FLAG.PostProcess

end

%% Plotting
if FLAG.PLOT

end


