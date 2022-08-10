%% Test Battery Controller
clear all; close all; clc;

%% Flags
loc_FLAG.Proportional = 1;
loc_FLAG.LowPass      = 0;

%% Load File Parameters
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\FunctionalityTests\Simulink\ActualBatterySystem\Semi_Explicit_Test_KPCont_Profile_CV_Test_1SmallStepSOC81.93.mat';
filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\Semi_Explicit_Test\Semi_Explicit_Test_KPCont_Profile_CV_Test_1SmallStep_IC4.0VSOC100.mat';
load(filename);

%% Controller Parameters
controller_tf = tf(1,1);
if loc_FLAG.Proportional
    K_p = 275;
    controller_tf = controller_tf * tf(K_p,1);
else
    K_p = 0;
end
if loc_FLAG.LowPass
    cutoff_fq = 7e9;
    controller_tf = controller_tf * tf([cutoff_fq],[1 , cutoff_fq]);
else
    cutoff_fq = 0;
end
[num,den] = tfdata(controller_tf);
num = cell2mat(num);
den = cell2mat(den);

%% Convert Function Handles into Strings
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
%% Simulation Parameters
t_final = 1; % [s]

% Output Matrix (Just looking at cell voltage)
idx_phi_ed_AN = P.phi_ed;

i = N.N_CV_CA(end);
index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
idx_phi_ed_CA = index_offset + P.phi_ed;

OutputMatrix = zeros(1,length(SIM.SV_IC));
OutputMatrix(idx_phi_ed_AN) = -1;
OutputMatrix(idx_phi_ed_CA) =  1;

%% Calculate the inverse of the mass matrix
% % M*y' = f(t,y)
% % Since GovnEqn = f, then I need to solve for y' before sending it through
% % the integrator block
% % y' = inv(M)*f
% mass_inv = pinv(SIM.M); % Psudo-inverse since mass is singular

% In Semi-Explicit form, the differential equations are solved separately
% from the algebraic. 
% u' = f1(t,u,v); 0 = f2(t,u,v)
% Need a mass matrix just for the diffEqs and then an index matrix to
% seperate diff and algebraic equations
mass_inv = SIM.M_DiffEq_Inv;

%% Initial Conditions
SV_IC = SV_soln(end,:)';
u_IC = SV_IC(SIM.Diff_idx);
v_IC = SV_IC(SIM.Algb_idx);
% i_user_IC = (SIM.Cell_Cap/50)/SIM.A_c;
i_user_IC = 0;
% u_IC = SIM.SV_IC(SIM.Diff_idx);
% v_IC = SIM.SV_IC(SIM.Algb_idx);

%% Run Simulink
mdl = 'Simple_Controller';
set_param('Simple_Controller','AlgebraicLoopSolver','Auto')
in = Simulink.SimulationInput(mdl);
in = in.setModelParameter('StartTime','0','StopTime',num2str(t_final));
% activeConfigObj = getActiveConfigSet(mdl);
% set_param(activeConfigObj,'RelTol','1e-6')
% activeConfigObj
out = sim(in);

%% Save Results



save_filename = ['F:\TylerFiles\GitHubRepos\BatteryModel\Model\FunctionalityTests\Simulink\ActualBatterySystem\Kp' num2str(K_p) 'omega_o' num2str(cutoff_fq) '.mat'];


%% Plot Results


