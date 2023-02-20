%% Initialization File
function [SIM,N,P,FLAG,RESULTS] = init(FLAG)
%% Convert Some FLAG to SIM
    SIM.Q_0 = FLAG.Q_0;
    SIM.R_0 = FLAG.R_0;
    SIM.Ts  = FLAG.Ts;
    SIM.tc  = FLAG.Ts;
    SIM.SOC = FLAG.SOC;


%% Set pointers for variables
    i = 1;
    if FLAG.DesOut.cell_voltage   % Cell Voltage
        P.cell_voltage = i; 
        RESULTS.Labels.title{i} = 'Cell Voltage';
        RESULTS.Labels.unit{i}  = 'Voltage [V]';
        i = i + 1;       
    end
    if FLAG.DesOut.delta_phi      % Delta Phi
        P.delta_phi = i; 
        RESULTS.Labels.title{i} = '\Delta \phi';
        RESULTS.Labels.unit{i}  = 'Voltage [V]';
        i = i + 1;    
    end
    if FLAG.DesOut.i_Far          % i_Far
        P.i_Far = i;
        RESULTS.Labels.title{i} = 'i_{Far}';
        RESULTS.Labels.unit{i}  = 'Current [A/m^2]';
        i = i + 1;
    end
    if FLAG.DesOut.eta          % eta
        P.eta = i;
        RESULTS.Labels.title{i} = '\eta';
        RESULTS.Labels.unit{i}  = 'Voltage [V]';
        i = i + 1;
    end
    if FLAG.DesOut.C_Liion        % C_Li+
        P.C_Liion = i; 
        RESULTS.Labels.title{i} = 'C_{Li^+}';
        RESULTS.Labels.unit{i}  = 'Concentration [kmol/m^3]';
        i = i + 1;
    end
    if FLAG.DesOut.C_Li           % C_Li
        P.C_Li = i;
        RESULTS.Labels.title{i} = 'C_{Li,surf}';
        RESULTS.Labels.unit{i}  = 'Concentration [kmol/m^3]';
        i = i + 1;
    end
    if FLAG.DesOut.delta_C_Li     % Delta C_Li
        P.delta_C_Li = i;
        RESULTS.Labels.title{i} = '\Delta C_{Li}';
        RESULTS.Labels.unit{i}  = 'Concentration [kmol/m^3]';
        i = i + 1;
    end
    if FLAG.DesOut.T              % Temperature
        P.T = i;
        RESULTS.Labels.title{i} = 'Temperature';
        RESULTS.Labels.unit{i}  = 'Temperature [K]';
        i = i + 1;
    end
    
    N.DesOut = i-1;


%% Output Matrix
% Load a file to get pointers
    filename = getSSFilename(FLAG);
    sys = load(filename);

    N.N_In  = sys.N.N_In;
        % I_user
    % Outputs
        % Cell Voltage
        % Delta Phi   @AN/SEP
        % Temperature @AN/SEP
        % C_Liion     @AN/SEP
        % X_surf      @AN/SEP
        % i_Far       @AN/SEP
        % eta         @AN/SEP
    SIM.OutputMatrix = zeros(N.DesOut , sys.N.N_SV_tot);
        if FLAG.DesOut.cell_voltage
        % Cell Voltage
            idx_phi_ed_AN = sys.P.phi_ed;

            i = sys.N.N_CV_CA(end);
            index_offset = (i-1)*sys.N.N_SV_CA + sys.N.N_SV_AN_tot + sys.N.N_SV_SEP_tot;
            idx_phi_ed_CA = index_offset + sys.P.phi_ed;

            SIM.OutputMatrix(P.cell_voltage,idx_phi_ed_AN) = -1; 
            SIM.OutputMatrix(P.cell_voltage,idx_phi_ed_CA) =  1;
        end
        % @AN/SEP
            i = sys.N.N_CV_AN(end);
            index_offset = (i-1)*sys.N.N_SV_AN;
        if FLAG.DesOut.delta_phi
        % Delta Phi   @AN/SEP
            idx = index_offset + sys.P.del_phi;
            SIM.OutputMatrix(P.delta_phi,idx) =  1; 
        end
        if FLAG.DesOut.i_Far
        % i_Far      @AN/SEP
            idx = index_offset + sys.P.i_PS;
            SIM.OutputMatrix(P.i_Far,idx) = 1;  
        end
        if FLAG.DesOut.eta
        % eta       @AN/SEP
            idx = index_offset + sys.P.V_2;
            SIM.OutputMatrix(P.eta,idx) = 1;
            idx = index_offset + sys.P.V_1;
            SIM.OutputMatrix(P.eta,idx) = -1;
        end
        if FLAG.DesOut.C_Liion
        % C_Liion     @AN/SEP
            idx = index_offset + sys.P.C_Liion;
            SIM.OutputMatrix(P.C_Liion,idx) = 1;
        end
        if FLAG.DesOut.C_Li
        % C_Li,surf   @AN/SEP
            idx = index_offset + sys.P.C_Li_surf_AN;
            SIM.OutputMatrix(P.C_Li,idx) = 1;%/AN.C_Li_max;
        end
        if FLAG.DesOut.delta_C_Li
        % Delta C_Li  @AN/SEP (Over entire radius of particle)
            idx = index_offset + sys.P.C_Li_surf_AN; % Surface Node
            SIM.OutputMatrix(P.delta_C_Li,idx) = 1;
            idx = index_offset + sys.P.C_Li;         % Most Interior Node
            SIM.OutputMatrix(P.delta_C_Li,idx) = -1;
        end
        if FLAG.DesOut.T
        % Temperature @AN/SEP
            idx = index_offset + sys.P.T;
            SIM.OutputMatrix(P.T,idx) = 1; 
        end


%% N Variables
    % N.measure is the number of outputs that are being measured
    % N.inputs  is the number of input sources
    % N.states  is the number of states used in the FOM

    N.states = sys.N.N_SV_tot; % FOM
    N.measur = 1;          % Should always be 1 
    N.inputs = sys.N.N_In; % Should always be 1 


%% Calculate Noise Power
    SIM.Qi = SIM.Q_0 * eye(N.inputs);
    
    Q       = SIM.Q_0 * eye(N.states);
    SIM.Qs = (diag(Q));
    
    R     = SIM.R_0 * eye(N.measur);
    SIM.R = (diag(R));
    
    SIM.pow_Qi = SIM.tc * SIM.Qi;
    SIM.pow_Qs = SIM.tc * SIM.Qs;
    SIM.pow_R  = SIM.tc * SIM.R;


%% Get IC from FOM
    filename = getImpulseFilename(FLAG);
    sys = load(filename);
    SIM.y_0_FOM = SIM.OutputMatrix * sys.SIM.SV_IC;
    SIM.x_0_FOM = sys.SIM.SV_IC;
    

%% Input Signal Parameters
% time_rest: The time before any type of input occurs, where the system holds its initial conditions
% time_ramp: ODE simulations need a continuous signal, therefore a step is converted into a quick 
%            ramp. This is the time over which the input source ramps up
% t_final_sim: Overall time for the simulation to fun for. It includes t_final and rest and ramp times
% t_vec: The discretized time vector to return the solution

    SIM.t_final_Relax_Step = (FLAG.N_samples+2)*SIM.Ts;
    SIM.t_vec_Relax_Step = 0:SIM.Ts:SIM.t_final_Relax_Step;
    input_load = ones(size(SIM.t_vec_Relax_Step));
    input_load(1) = 0;
    input_load(2) = 0;
    
    if FLAG.SOC < 91
        input_load = -1*input_load;
    end
    
    SIM.InputSignal = [SIM.t_vec_Relax_Step' input_load'];


end


%%%%%%%%%%%%%%%%%%%%%%%%%% OOOOOOOOOOOLD %%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Set All Functions to Incomplete
%     SIM.AnalysisComplete.Initialization             = 0;
%     SIM.AnalysisComplete.NoNoiseCompare.ode         = 0;
%     SIM.AnalysisComplete.NoNoiseCompare.SS_CT       = 0;
%     SIM.AnalysisComplete.NoNoiseCompare.SS_DT       = 0;
%     SIM.AnalysisComplete.NoNoiseCompare.ROM_Mlab    = 0;
%     SIM.AnalysisComplete.NoNoiseCompare.ROM_HoKal   = 0;
%     SIM.AnalysisComplete.NoisyPlant                 = 0;
%     SIM.AnalysisComplete.Estimator                  = 0;
%     SIM.AnalysisComplete.Est_Error_calc             = 0;
%     SIM.AnalysisComplete.GenComparData              = 0;
%     SIM.AnalysisComplete.ComparSVD2Pinf             = 0;

%%
% 
% F_step: Force after the step has taken place
%     SIM.F_step = SIM.F_init + SIM.F_step_height;

% SIM.time_rest = 8*SIM.tc+SIM.Ts;
%     if FLAG.InputMode == 2
%         SIM.omega = 2*pi*SIM.fq; % [rad/s], Frequency for the sinusoidal input
%     end




%% Mass Matrix
%     SIM.Mass = eye(N.states);


%% Initial Conditions
%     SIM.x_0 = [0 ,0 ,0 ,0 ]';

    
%% Initialize results struct
%     RESULTS = struct();




%     SIM.time_rest = 2*SIM.Ts;
%     SIM.time_ramp = SIM.tc/10;
%     
%     SIM.t_final_sim = SIM.t_final + SIM.time_rest + SIM.time_ramp;
%     SIM.t_vec = 0:SIM.tc:SIM.t_final_sim;