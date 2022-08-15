%% Test Battery In Simulink
clear all; close all; clc;

%% Load File
filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\TestingSimulink\TestingSimulink_Polar_1.00C_D.mat';
load(filename)
[current_file_path,~,~] = fileparts(mfilename('fullpath'));
cd(current_file_path)
%% ODE Simulation
% Simulation Parameters
    tspan = SIM.tspan;
    Tol.Abs = 1E-7;
    Tol.Rel = 1E-7;
    
    events = @(t,SV) batt_events(t,SV,SIM,P,N,FLAG);
    
    options = odeset('RelTol' ,Tol.Rel,      ...
                     'AbsTol' ,Tol.Abs,      ...
                     'Mass'   ,SIM.M,        ...
                     'Events' ,events);%,       ...
                    %'MaxStep',1e2);
    SV_IC = SIM.SV_IC;
    i_user = 0;
    tStart = tic;
    SOLN = ode15s(@(t,SV)batt_GovEqn_test(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user),tspan,SV_IC,options);
    tODE = toc(tStart);
    t_soln_ode  = SOLN.x';
    SV_soln_ode = SOLN.y';

% From Post Processing
    N_t_steps = length(t_soln_ode);
    idx = 1:1:N_t_steps;
    phi_ed       = zeros( N_t_steps , N.N_CV_tot );
    max_SV = max( N.N_SV_AN , N.N_SV_CA );
    SV     = zeros(max_SV+1 , N.N_CV_tot, N_t_steps );
    for i = 1:N_t_steps
        % Go through the solution vector and reshape every SV to 2D (3D matrix)
        SV_temp    = SV1Dto2D( SV_soln_ode( idx(i) , : ) , N , P , FLAG );
        SV( : , : , i ) = addPhiEl2SV(SV_temp,P,N);
        
        %         TemperatureK( i , : ) = SV( P.T , : , i );
        %         
        %         C_Liion( i , : )                = SV( P.C_Liion              , :              , i );
        %         C_Li( : , : , i )               = SV( P.C_Li:P.C_Li_surf_max , :              , i );
        %         C_Li_surf( i , N.CV_Region_AN ) = SV( P.C_Li_surf_AN         , N.CV_Region_AN , i ); 
        %         C_Li_surf( i , N.CV_Region_CA ) = SV( P.C_Li_surf_CA         , N.CV_Region_CA , i ); 
        
        phi_ed( i , : ) = SV( P.phi_ed , : , i );
        %         phi_el( i , : ) = SV( P.phi_el , : , i );
                
        %         eta_AN       = SV(P.V_2 , N.CV_Region_AN ,i) - SV(P.V_1 , N.CV_Region_AN ,i);
        %         eta_CA       = SV(P.V_2 , N.CV_Region_CA ,i) - SV(P.V_1 , N.CV_Region_CA ,i);
        %         eta( i , : ) = [eta_AN  , NaN(1,N.N_CV_SEP) ,eta_CA  ];
        %         
        %         V_SEI(   i , : ) = SV( P.V_1    , : , i) - SV( P.phi_el , : , i);
        %         del_phi( i , : ) = SV( P.del_phi,:  , i);
    end
    % Cell Voltage
    cell_voltage_ode = phi_ed(:,end) - phi_ed(:,1);

%% Simulink
% Simulink will use Linear Controller Method
% Pseudo-inverse of Mass
    [U, S, V] = svd(SIM.M);
    threshold = 1e-7;
    S_cross = pinv(S,threshold);
    M_cross = V*S_cross*U';

% Calculate the A matrix
    t = 0;
    i_user  = SIM.i_user_amp;
%     i_user  = -SIM.A_c*SIM.i_user_amp; %%%%%%%%%Is this needed?
    TOL.Rel = 1e-3;
    TOL.Abs = 1e-6;
    A = zeros(N.N_SV_tot , N.N_SV_tot);
    dSVdt_init = batt_GovEqn_test(t,SIM.SV_IC,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user);
    for i = 1:N.N_SV_tot
        p     = zeros(N.N_SV_tot,1);
        p(i)  = TOL.Rel * SIM.SV_IC(i) + TOL.Abs;
        SV_p  = SIM.SV_IC + p;
        dSVdt = batt_GovEqn_test(t,SV_p,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user);
            
        A(:,i) = ( dSVdt - dSVdt_init ) / p(i);
    end

% Calculate differential variable dynamics
    A_cross = M_cross*A;
    [ol_vec , ol_poles] = eig(M_cross*A);
%     fastest_diff_pole = min(diag(ol_poles));
    fastest_diff_pole = min(diag(real(ol_poles)));
    ol_poles_diag = diag(real(ol_poles));
    
    reduced_ol_eig_vec = real(ol_vec);
    for i = 1:length(reduced_ol_eig_vec)
        for j = 1:length(reduced_ol_eig_vec)
            if abs(reduced_ol_eig_vec(i,j))<1e-5
                reduced_ol_eig_vec(i,j) = 0;
            else
                reduced_ol_eig_vec(i,j) = reduced_ol_eig_vec(i,j);
            end
        end
    end

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
    K = diag(diag_vec);

% Calculate full system poles (this is a check)
    A_comb = M_cross*A  - K*U_NV_NT*A ;
    [sys_eig_vec ,sys_poles] = eig(A_comb);
    sys_poles_diag = diag(real(sys_poles));
%     find(sys_poles_diag>0)
    reduced_sys_eig_vec = real(sys_eig_vec);
    for i = 1:length(reduced_sys_eig_vec)
        for j = 1:length(reduced_sys_eig_vec)
            if abs(reduced_sys_eig_vec(i,j))<1e-5
                reduced_sys_eig_vec(i,j) = 0;
            else
                reduced_sys_eig_vec(i,j) = reduced_sys_eig_vec(i,j);
            end
        end
    end

%% Make Output matrix
    N.N_In  = 1;
        % I_user
    N.N_Out = 5;
        % Cell Voltage
        % Delta Phi   @AN/SEP
        % Temperature @AN/SEP
        % C_Liion     @AN/SEP
        % X_surf      @AN/SEP
    SIM.OutputMatrix = zeros(N.N_Out , N.N_SV_tot);
    j = 0;
        % Cell Voltage
            j = j+1;
            idx_phi_ed_AN = P.phi_ed;

            i = N.N_CV_CA(end);
            index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
            idx_phi_ed_CA = index_offset + P.phi_ed;

            SIM.OutputMatrix(j,idx_phi_ed_AN) = -1;
            SIM.OutputMatrix(j,idx_phi_ed_CA) =  1;
        % @AN/SEP
            i = N.N_CV_AN(end);
            index_offset = (i-1)*N.N_SV_AN;
        % Delta Phi   @AN/SEP
            j = j+1;
            idx = index_offset + P.del_phi;
            SIM.OutputMatrix(j,idx) =  1;
        % Temperature @AN/SEP
            j = j+1;
            idx = index_offset + P.T;
            SIM.OutputMatrix(j,idx) = 1;
        % C_Liion     @AN/SEP
            j = j+1;
            idx = index_offset + P.C_Liion;
            SIM.OutputMatrix(j,idx) = 1;
        % X_surf      @AN/SEP
            j = j+1;
            idx = index_offset + P.C_Li_surf_AN;
            SIM.OutputMatrix(j,idx) = 1/AN.C_Li_max;

%% Make constant variables used in Simulink
    i_user_sim = SIM.i_user_amp;
    SV_IC = SIM.SV_IC;
    OutputMatrix = SIM.OutputMatrix;

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
% SIM = rmfield(SIM,'ControllerHandle');

%% Run Simulink
t_final = t_soln_ode(end); %[s]
mdl = 'battery_system';
% set_param('Simple_Controller_SS_Model','AlgebraicLoopSolver','Auto')
in = Simulink.SimulationInput(mdl);
in = in.setModelParameter('StartTime','0','StopTime',num2str(t_final));
% activeConfigObj = getActiveConfigSet(mdl);
% set_param(activeConfigObj,'RelTol','1e-6')
% activeConfigObj
tStart = tic;
out = sim(in);
tSim = toc(tStart);

%% Calculations
t_soln = out.tout;
SV_soln = zeros(length(out.tout),N.N_SV_tot);
for i = 1:length(out.tout)
    SV_soln(i,:) = out.SV(:,:,i)';
end

% From Post Processing
    N_t_steps = length(t_soln);
    idx = 1:1:N_t_steps;
    phi_ed       = zeros( N_t_steps , N.N_CV_tot );
    max_SV = max( N.N_SV_AN , N.N_SV_CA );
    SV     = zeros(max_SV+1 , N.N_CV_tot, N_t_steps );
    for i = 1:N_t_steps
        % Go through the solution vector and reshape every SV to 2D (3D matrix)
        SV_temp    = SV1Dto2D( SV_soln( idx(i) , : ) , N , P , FLAG );
        SV( : , : , i ) = addPhiEl2SV(SV_temp,P,N);
        
        %         TemperatureK( i , : ) = SV( P.T , : , i );
        %         
        %         C_Liion( i , : )                = SV( P.C_Liion              , :              , i );
        %         C_Li( : , : , i )               = SV( P.C_Li:P.C_Li_surf_max , :              , i );
        %         C_Li_surf( i , N.CV_Region_AN ) = SV( P.C_Li_surf_AN         , N.CV_Region_AN , i ); 
        %         C_Li_surf( i , N.CV_Region_CA ) = SV( P.C_Li_surf_CA         , N.CV_Region_CA , i ); 
        
        phi_ed( i , : ) = SV( P.phi_ed , : , i );
        %         phi_el( i , : ) = SV( P.phi_el , : , i );
                
        %         eta_AN       = SV(P.V_2 , N.CV_Region_AN ,i) - SV(P.V_1 , N.CV_Region_AN ,i);
        %         eta_CA       = SV(P.V_2 , N.CV_Region_CA ,i) - SV(P.V_1 , N.CV_Region_CA ,i);
        %         eta( i , : ) = [eta_AN  , NaN(1,N.N_CV_SEP) ,eta_CA  ];
        %         
        %         V_SEI(   i , : ) = SV( P.V_1    , : , i) - SV( P.phi_el , : , i);
        %         del_phi( i , : ) = SV( P.del_phi,:  , i);
    end
    % Cell Voltage
    cell_voltage = phi_ed(:,end) - phi_ed(:,1);
    
%% Plot
%Cell Voltage
    figure
    hold on
    plot(t_soln    ,cell_voltage    ,'ro','LineWidth',2,'DisplayName','Simulink')
    plot(t_soln_ode,cell_voltage_ode,'k' ,'LineWidth',2,'DisplayName','ode15s')
    title('Cell Voltage')
    xlabel('Time (s)')
    ylabel('Voltage (V)')
    xlim([0,t_soln(end)])
    lgn = legend;
%     exportgraphics(gcf,'ODEvSIMULINK.png','Resolution',1000)