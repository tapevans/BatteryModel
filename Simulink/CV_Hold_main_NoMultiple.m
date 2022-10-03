%% Test Battery In Simulink
clear all; 
% close all; 
clc;

%% Parameters
    tfinal = 30*60; % [s], Final Simulation Time
    tfinal = 2*60*60; % [s], Final Simulation Time
    tfinal = 3000;
%     tfinal = 10*60; % [s], Final Simulation Time
%     tfinal = 1e-5; % [s], Final Simulation Time

    t_r = 1e0;
    
    model_filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\TestingSimulink\TestingSimulink_KPCont_SingleStepTo50SOCSOC0.mat';
    equil_filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\TestingSimulink\TestingSimulink_SS_EIS_SOC50.mat';
    save_filename  = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\TestingSimulink\SimulinkResults.mat';
%     save_filename  = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\TestingSimulink\SimulinkResults_CC_C20SOC50_CV.mat';

    localFLAG.CalcSS      = 1; % 1 if calculating SS from the end of CC, 0 if using the equilibrium model for SS
    localFLAG.PID         = 1;
    localFLAG.LQR         = 0;
    localFLAG.AddTunedController = 0;
    localFLAG.PostProcess = 1;
    localFLAG.PLOT        = 1;
    localFLAG.SaveResults = 0;

%% Update Working Directory
    [current_file_path,~,~] = fileparts(mfilename('fullpath'));
    cd(current_file_path)

% Add 'Model' Folder
    idx = find(current_file_path == '\',1,'last');
    model_folder = [ current_file_path(1:idx) 'Model'];
    addpath(genpath(model_folder));

% Remove Functionality Test Folder
    rmpath('F:\TylerFiles\GitHubRepos\BatteryModel\Model\FunctionalityTests')

%% Load Model File
model = load(model_filename);
SV_IC = model.SV_soln(end,:)';
SIM   = model.SIM;
N     = model.N;

% multiple = -model.SIM.A_c^-1;
% CC_init  = multiple*model.i_user_soln(end); % [A m^-2], The current used at the last time step of the previous mode. Initializes the CV with this current
CC_init  = model.i_user_soln(end); % [A m^-2], The current used at the last time step of the previous mode. Initializes the CV with this current
% if ~localFLAG.CalcSS
%     clear model;
% end

%% Load Equilibrium State Space File
if ~localFLAG.CalcSS
    equil = load(equil_filename);
    A = equil.A; B = equil.B; C = equil.C; D = equil.D;
    C = C(1,:); % Cell_Voltage is the first row
    D = 0;
    %     syseq_ctrb = ctrb(ss(A,B,C,D,'E',SIM.M))
    CV = C*model.SV_soln(end,:)'; %!!!!!!!Can change this back to SV_IC, reference voltage where the hold occurs
end

%% Get SS Matricies
if localFLAG.CalcSS
AN    = model.AN;
CA    = model.CA;
SEP   = model.SEP;
EL    = model.EL;
SIM   = model.SIM;
CONS  = model.CONS;
P     = model.P;
N     = model.N;
FLAG  = model.FLAG;
PROPS = model.PROPS;


[A,B,C,D] = getSSfromSV(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,SV_IC,CC_init);
C = C(1,:); % Cell_Voltage is the first row
D = 0;
CV = C*SV_IC;

% i = 1;
%     P.OM.cell_volt = i; i = i + 1;
%     P.OM.del_phi   = i; i = i + 1;
%     P.OM.temp      = i; i = i + 1;
%     P.OM.C_Liion   = i; i = i + 1;
%     P.OM.X_surf    = i; i = i + 1;
%     P.OM.i_Far     = i; i = i + 1;
%     P.OM.eta       = i; i = i + 1;

% N.N_In = 1;
% N.N_Out = length(fieldnames(P.OM));

% SIM.OutputMatrix = zeros(N.N_Out , N.N_SV_tot);
%     % Cell Voltage
%         idx_phi_ed_AN = P.phi_ed;
% 
%         i = N.N_CV_CA(end);
%         index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
%         idx_phi_ed_CA = index_offset + P.phi_ed;
% 
%         SIM.OutputMatrix(P.OM.cell_volt,idx_phi_ed_AN) = -1;
%         SIM.OutputMatrix(P.OM.cell_volt,idx_phi_ed_CA) =  1;
%     % @AN/SEP
%         i = N.N_CV_AN(end);
%         index_offset = (i-1)*N.N_SV_AN;
%     % Delta Phi   @AN/SEP
%         idx = index_offset + P.del_phi;
%         SIM.OutputMatrix(P.OM.del_phi,idx) =  1;
%     % Temperature @AN/SEP
%         idx = index_offset + P.T;
%         SIM.OutputMatrix(P.OM.temp,idx) = 1;
%     % C_Liion     @AN/SEP
%         idx = index_offset + P.C_Liion;
%         SIM.OutputMatrix(P.OM.C_Liion,idx) = 1;
%     % X_surf      @AN/SEP
%         idx = index_offset + P.C_Li_surf_AN;
%         SIM.OutputMatrix(P.OM.X_surf,idx) = 1/AN.C_Li_max;
%     % i_Far      @AN/SEP
%         idx = index_offset + P.i_PS;
%         SIM.OutputMatrix(P.OM.i_Far,idx) = 1;
%     % Eta      @AN/SEP
%         idx = index_offset + P.V_2;
%         SIM.OutputMatrix(P.OM.eta,idx) = 1;
%         idx = index_offset + P.V_1;
%         SIM.OutputMatrix(P.OM.eta,idx) = -1;

% i = 1;
% P.SS.omega    = i; i = i + 1;
% P.SS.Z_mag    = i; i = i + 1;
% P.SS.Z_Re     = i; i = i + 1;
% P.SS.Z_Im     = i; i = i + 1;
% P.SS.Z_dB     = i; i = i + 1;
% P.SS.Z_ps_deg = i; i = i + 1;

% SIM.SV_IC  = SV_IC;
% % SIM.freq   = logspace(-1,11,101);
% SIM.i_user = -model.SIM.Cell_Cap/20;


% [A,B,C,D,~] = getSSImpedance(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS);
% C = C(1,:); % Cell_Voltage is the first row
% D = 0;
% CV = C*SV_IC;


%     multiple = -SIM.A_c^-1;
%     sys = multiple*ss(A,B,C,D,'E',SIM.M);
%     sys_min = minreal(sys);
% %     postProcessComplete = 1;
% %     save('Calc_SS','sys','postProcessComplete')
%     
%     A = sys.A; B = sys.B; C = sys.C; D = sys.D;

end


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

%% LQG
if localFLAG.LQR
N_states = N.N_SV_tot;
N_output = 1;

rho = 1e0;
Q_LQG = rho * (sys_min.C' * sys_min.C) ;%+ 1e-10*eye(length(sys.A));
Q_LQG = blkdiag(Q_LQG,1);
R_LQG = 1e0; % Inputs is length 1, makes R (1x1) %%%%!!!! Needs better comment
N_LQG = zeros(length(sys_min.A)+1,1); % (states+1 x inputs)

[K_LQI,~,~] = lqi(sys_min,Q_LQG,R_LQG,N_LQG);

% Kalman Filter
% Modify the system?
LQGsysMOD = sys_min;
LQGsysMOD.OutputName = 'Voltage';
LQGsysMOD.InputName  = 'CurrentWNoise';

Sum = sumblk('CurrentWNoise = Current + w');
sys_LQG = connect(LQGsysMOD,Sum,{'Current','w'},'Voltage');


Qn = rho;
Rn = R_LQG;
Nn = 0;
[kest,~,~] = kalman(sys_LQG,Qn,Rn,Nn);

% Combining Kalman filter and Gains
Controller = lqgtrack(kest,K_LQI);
end

%% PID Controller
if localFLAG.PID

AN    = model.AN;
CA    = model.CA;
SEP   = model.SEP;
EL    = model.EL;
SIM   = model.SIM;
CONS  = model.CONS;
P     = model.P;
N     = model.N;
FLAG  = model.FLAG;
PROPS = model.PROPS;


% Going to use Configuration 1 of the Control System Designer
% multiple = -SIM.A_c^-1;
% plant = multiple*ss(A,B,C,D,'E',SIM.M);
plant = ss(A,B,C,D,'E',SIM.M);

% Control System Designer
    % s = sisoinit(1);
    % s.G.Value = plant;
    % controlSystemDesigner(s)

% PID Tuner
% opt = pidtuneOptions;
% Types = {'PI','PID','PI2','PID2'};
% Types = {'PI','PID'};
% for i = 1:length(Types)
%     allController{i} = pidtune(plant, Types{i});
% end
% Controller_PI   = pidtune(plant,'PI'  );
% Controller_PID  = pidtune(plant,'PID' );
% Controller_PI2  = pidtune(plant,'PI2' );
% Controller_PID2 = pidtune(plant,'PID2');

% pidTuner(plant,allController{2})

wc = 2/t_r;
Controller = pidtune(plant, 'PI', wc);
K_P = Controller.Kp;
K_I = Controller.Ki;
K_D = 0;
% K_D = Controller.Kd;




end
% if localFLAG.AddTunedController
% load('TunedController.mat');
% allController{end+1} = myTunedController;
% end

%% Extract Gains
% for i = 1:length(allController)
%     Kp_vec(i) = allController{i}.Kp;
%     Ki_vec(i) = allController{i}.Ki;
%     Kd_vec(i) = allController{i}.Kd;
% end

%% Convert PID to tf or ss
% for i = 1:length(allController)
% %     allController{i} = tf(allController{i});
% %     allController{i} = ss(allController{i});
% end

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
% mdl = 'battery_system_CV_PID';
mdl = 'battery_system_CV_PID_Block_NoMultiple';
in = Simulink.SimulationInput(mdl);
in = in.setModelParameter('StartTime','0','StopTime',num2str(tfinal));
% for i = 1:length(allController)
%     K_P = Kp_vec(i);
%     K_I = Ki_vec(i);
%     K_D = Kd_vec(i);
%     Controller = allController{i};
    tStart = tic;
    out = sim(in);
    t_solve = toc(tStart);
% end

%% Post-Processing
if localFLAG.PostProcess
%     for i = 1:length(out)
        out.i_user = reshape(out.i_user,[],1);
%         out{i}.i_user = reshape(out{i}.i_user,[],1);
%     end

% Combined CC and CV Results
    adj_time_CV = out.tout + model.t_soln(end);
    total_time = [model.t_soln; adj_time_CV];
    voltage_CC = model.SV_soln(:,341) - model.SV_soln(:,3);
    voltage_CV = out.SV(:,341)        - out.SV(:,3);
    total_voltage = [ voltage_CC ; voltage_CV];
end

%% Plotting
if localFLAG.PLOT
%     myTitle = Types;
%     myTitle{end+1} = 'myTunedController';
%     for i = 1:length(out)
        figure
        yyaxis left
        plot(out.tout,out.SV(:,341)-out.SV(:,3))
        ylabel('Voltage [V]')
        yline(CV,'--')
        
        yyaxis right
        plot(out.tout,out.i_user)
        title(['t_r = ', num2str(t_r), 's'])
        ylabel('Current [A m^{-2}]')

        xlabel('Time [s]')
%     end

    % CC with CV
    figure
    plot(total_time , total_voltage , 'k' , 'LineWidth' , 2)
    title(['Combined CC and CV: t_r = ', num2str(t_r), 's'])
    xlabel('Time [s]')
    ylabel('Voltage [V]')

end

%% OOOOOOLD
%     p_OS = 1e-6; % [%], Percent Overshoot
%     t_s  = 2e-4; % [s], Settling Time



% sorteig = sort(real(eig(A)));
% [vec , poles] = eig(A);
% poles = diag(poles);
% find(poles>0)
% reduced_vec = real(vec);
%     for i = 1:length(reduced_vec)
%         for j = 1:length(reduced_vec)
%             if abs(reduced_vec(i,j))<1e-2
%                 reduced_vec(i,j) = 0;
%             else
%                 reduced_vec(i,j) = reduced_vec(i,j);
%             end
%         end
%     end

%     sorteig_cross = sort(real(eig(A_cross)));
%     ol_poles_diag = diag(real(ol_poles));
%     idx_cross = find(ol_poles_diag>0);
%     reduced_ol_vec = real(ol_vec);
%     for i = 1:length(reduced_ol_vec)
%         for j = 1:length(reduced_ol_vec)
%             if abs(reduced_ol_vec(i,j))<1e-2
%                 reduced_ol_vec(i,j) = 0;
%             else
%                 reduced_ol_vec(i,j) = reduced_ol_vec(i,j);
%             end
%         end
%     end

%     U_colm = U(:,1:r);
%     V_colm = V(:,1:r);



% % Calculate New System A,B
%     A_comb = M_cross*A - K_algb*U_NV_NT*A;
%     B_comb = M_cross*B - K_algb*U_NV_NT*B;
%     %     B_comb = M_cross*B;
    

% % Save the system for troubleshooting
%     multiple = 1;
% %     multiple = -SIM.A_c^-1;
%     sys = multiple*ss(A_comb,B_comb,C,D);
%     postProcessComplete = 1;
%     save('LCNull_SS','sys','postProcessComplete')

% % Calc eig
%     sorteig_comb = sort(real(eig(A_comb)));
%     [sys_eig_vec ,sys_poles] = eig(A_comb);
%     sys_poles_diag = diag(real(sys_poles));
%     idx_comb = find(sys_poles_diag>0);
%     reduced_sys_eig_vec = real(sys_eig_vec);
%     for i = 1:length(reduced_sys_eig_vec)
%         for j = 1:length(reduced_sys_eig_vec)
%             if abs(reduced_sys_eig_vec(i,j))<1e-2
%                 reduced_sys_eig_vec(i,j) = 0;
%             else
%                 reduced_sys_eig_vec(i,j) = reduced_sys_eig_vec(i,j);
%             end
%         end
%     end
%     for i = 1:length(A_comb)
%         if i == 1
%             Controllability_comb = [ B_comb];
%         else
%             Controllability_comb = [ Controllability_comb, A_comb^(i-1)*B_comb];
%         end
%     end
% %     Controllability_comb = [ B_comb A_comb*B_comb A_comb^2*B_comb];
%     rank_comb = rank(Controllability_comb)

% idx_comb'
% sys_poles_diag(idx_comb)'
% D;
%% Calculate Closed Loop Pole Dynamics
%     zeta = sqrt( ((log(p_OS/100))^2)/ ( pi^2 + (log(p_OS/100))^2  )   );
%     omega_n = 4.6/( zeta*t_s );
% 
%     pole_p = -omega_n*(zeta + 1i*sqrt(1-zeta^2));
%     pole_m = -omega_n*(zeta - 1i*sqrt(1-zeta^2));
%     fast_pole = -5*zeta*omega_n;
% %     fast_pole_vec = linspace(0,1,(N.N_SV_tot-2))';
%     fast_pole_vec = -(0:1:N.N_SV_tot-3)';
%     fast_pole_vec = fast_pole_vec + fast_pole;
% 
%     cl_poles = [pole_p ; pole_m ; fast_pole_vec];

%% Calculate Controller Matricies
% sys   = ss(A,B,C,D);  Co    = ctrb(A,B);
% sysr  =  minreal(sys);Co_r  = ctrb(sysr.A ,sysr.B);
% sysrs = sminreal(sys);Co_rs = ctrb(sysrs.A,sysrs.B);

% sys_c   = ss(A_comb,B_comb,C,D);  Co_c    = ctrb(A_comb,B_comb);
% sysr_c  =  minreal(sys_c);        Co_r_c  = ctrb(sysr_c.A ,sysr_c.B);
% sysrs_c = sminreal(sys_c);        Co_rs_c = ctrb(sysrs_c.A,sysrs_c.B);

% step(sys,'r',sysr,'--g',sysrs,'--b')
% step(sys,sys_c)

% % Determine the Controller Values
%     K = place(A_comb,B_comb,cl_poles);
    
% K = place(A,B,cl_poles);
% K_r = place(sysr.A,sysr.B,cl_poles(1:length(sysr.B)));
% K_rs = place(sysrs.A,sysrs.B,cl_poles(1:length(sysrs.B)));
% K_c = place(A_comb,B_comb,cl_poles);
% K_r_c = place(sysr_c.A,sysr_c.B,cl_poles(1:length(sysr_c.B)));
% K_rs_c = place(sysrs_c.A,sysrs_c.B,cl_poles(1:length(sysrs_c.B)));

%     Co = ctrb(A,B); % Controllability Matrix
%     rank(Co);
%     unco = length(A_comb) - rank(Co); % Number of Uncontrollable State

% % Determine N_x and N_u (See notes or control textbook)
%     zero_vec = zeros(N.N_SV_tot,1);
%     a = [A_comb , B_comb; C , D];
% %     a = [A , B; C , D];
%     b = [zeros(N.N_SV_tot,1) ;1];
%     N_vec = a\b;
% 
%     N_x = N_vec(1:N.N_SV_tot,1);
%     N_u = N_vec(N.N_SV_tot+1,1);



% t_s  = 2e-4; % [s], Settling Time
% t_s  = 1; % [s], Settling Time
% sys = ss(A_comb,B_comb,C,D);
% sys_ctrb = ctrb(sys)



% QXU = blkdiag(Q,R);
% QXU = eye(size(QXU));
% 
% QWU = eye(size(QXU)); % Needs to be positive definite
% 
% rho_I = 1e0;
% QI = rho_I*eye(N_output);
% 
% 
% reg_SS = lqg(sys,QXU,QWU,QI);
% A_SS = reg_SS.A;
% B_SS = reg_SS.B;
% C_SS = reg_SS.C;
% D_SS = reg_SS.D;