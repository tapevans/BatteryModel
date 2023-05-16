%% Initialization File
function [SIM,N,P,FLAG,RESULTS] = init(FLAG)
%% Convert Some FLAG to SIM
    SIM.Q_0   = FLAG.Q_0;
    SIM.R_0   = FLAG.R_0;
    SIM.Ts    = FLAG.Ts;
    SIM.tc    = FLAG.Ts;
    SIM.SOC   = FLAG.SOC;
    SIM.r_max = FLAG.r_max;
    SIM.LargeQMultiply = FLAG.LargeQMultiply;
    SIM.ResetStep      = FLAG.ResetStep;
    SIM.Q_Add          = FLAG.Q_Add;

    if FLAG.InputMode == 5
        SIM.Tswitch = FLAG.Tswitch;
        SIM.PRBSAmp = FLAG.PRBSAmp;
    end

    SIM.offsetROM_IC_Rel = FLAG.offsetROM_IC_Rel;
    SIM.offsetROM_IC_Abs = FLAG.offsetROM_IC_Abs;


%% Set pointers for variables
    i = 1;
    if FLAG.DesOut.cell_voltage   % Cell Voltage
        P.cell_voltage = i; 
        RESULTS.Labels.title{i} = 'Cell Voltage';
        RESULTS.Labels.unit{i}  = 'Voltage [V]';
        RESULTS.foldername{i}   = 'CellVoltage';
        i = i + 1;       
    end
    if FLAG.DesOut.delta_phi      % Delta Phi
        P.delta_phi = i; 
        RESULTS.Labels.title{i} = '\Delta \phi';
        RESULTS.Labels.unit{i}  = 'Voltage [V]';
        RESULTS.foldername{i}   = 'DeltaPhi';
        i = i + 1;    
    end
    if FLAG.DesOut.i_Far          % i_Far
        P.i_Far = i;
        RESULTS.Labels.title{i} = 'i_{Far}';
        RESULTS.Labels.unit{i}  = 'Current [A/m^2]';
        RESULTS.foldername{i}   = 'iFar';
        i = i + 1;
    end
    if FLAG.DesOut.eta          % eta
        P.eta = i;
        RESULTS.Labels.title{i} = '\eta';
        RESULTS.Labels.unit{i}  = 'Voltage [V]';
        RESULTS.foldername{i}   = 'eta';
        i = i + 1;
    end
    if FLAG.DesOut.C_Liion        % C_Li+
        P.C_Liion = i; 
        RESULTS.Labels.title{i} = 'C_{Li^+}';
        RESULTS.Labels.unit{i}  = 'Concentration [kmol/m^3]';
        RESULTS.foldername{i}   = 'C_Liion';
        i = i + 1;
    end
    if FLAG.DesOut.C_Li           % C_Li
        P.C_Li = i;
        RESULTS.Labels.title{i} = 'C_{Li,surf}';
        RESULTS.Labels.unit{i}  = 'Concentration [kmol/m^3]';
        RESULTS.foldername{i}   = 'C_Li_surf';
        i = i + 1;
    end
    if FLAG.DesOut.delta_C_Li     % Delta C_Li
        P.delta_C_Li = i;
        RESULTS.Labels.title{i} = '\Delta C_{Li}';
        RESULTS.Labels.unit{i}  = 'Concentration [kmol/m^3]';
        RESULTS.foldername{i}   = 'DeltaC_Li';
        i = i + 1;
    end
    if FLAG.DesOut.T              % Temperature
        P.T = i;
        RESULTS.Labels.title{i} = 'Temperature';
        RESULTS.Labels.unit{i}  = 'Temperature [K]';
        RESULTS.foldername{i}   = 'Temperature';
        i = i + 1;
    end
    
    N.DesOut = i-1;


%% Output Matrix
% Load a file to get pointers
    filename = getSSFilename(FLAG);
    sys = load(filename);

    N.N_In  = sys.N.N_In;
        % i_user

    % Outputs
        % Cell Voltage
        % Delta Phi   @AN/SEP
        % i_Far       @AN/SEP
        % eta         @AN/SEP
        % C_Liion     @AN/SEP
        % C_Li        @AN/SEP
        % Delta C_Li  @AN/SEP
        % Temperature @AN/SEP
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
    
    % Wrong IC
    y_0_FOM_Offset = (SIM.y_0_FOM * FLAG.offsetROM_IC_Rel) + FLAG.offsetROM_IC_Abs + SIM.y_0_FOM;
    idx = find(y_0_FOM_Offset == FLAG.offsetROM_IC_Abs);
    for j = idx
        y_0_FOM_Offset(j) = FLAG.offsetROM_IC_Zero;
    end
    SIM.y_0_FOM_Offset = y_0_FOM_Offset;


%% Input Signal Parameters
% time_rest: The time before any type of input occurs, where the system holds its initial conditions
% time_ramp: ODE simulations need a continuous signal, therefore a step is converted into a quick 
%            ramp. This is the time over which the input source ramps up
% t_final_sim: Overall time for the simulation to fun for. It includes t_final and rest and ramp times
% t_vec: The discretized time vector to return the solution
% if ~FLAG.Analysis.GenComparData %%%%%%%Delete Later
    if FLAG.InputMode == 5
        filename = getPRBSFilename(FLAG);
        sysSIM   = load(filename,'SIM');

        t_vec  = (0:SIM.Ts:sysSIM.SIM.profile_time(end))';
        i_user = interp1(sysSIM.SIM.profile_time , sysSIM.SIM.profile_current , t_vec);

        %SIM.InputSignal = [t_vec , i_user];
        SIM.InputSignal = [t_vec , [i_user(2:end);i_user(end)]]; % Proper ZOH
    else
        %%%%%%%%%%%%%%%%%%!!!!!!!!!!! Fix else for proper ZOH
        SIM.t_final_Relax_Step = (FLAG.N_samples+2)*SIM.Ts;
        SIM.t_vec_Relax_Step = 0:SIM.Ts:SIM.t_final_Relax_Step;


        % Add a time right after the step
        SIM.t_vec_Relax_Step = [SIM.t_vec_Relax_Step(1:3) , SIM.t_vec_Relax_Step(3)+SIM.Ts/10 , SIM.t_vec_Relax_Step(4:end)];

        input_load = ones(size(SIM.t_vec_Relax_Step));
        input_load(1) = 0;
        input_load(2) = 0;
        input_load(3) = 0;

        if FLAG.SOC < 91
            input_load = -1*input_load;
        end

        SIM.InputSignal = [SIM.t_vec_Relax_Step' input_load'];
    end
% end
end


%%%%%%%%%%%%%%%%%%%%%%%%%% OOOOOOOOOOOLD %%%%%%%%%%%%%%%%%%%%%%%%%% 
% if 0
%     figure
%     hold on
%     plot(sysSIM2.t_soln        ,sysSIM2.i_user        ,'k','LineWidth' ,2,'DisplayName','ODE')
%     plot(SIM.InputSignal(:,1)  ,SIM.InputSignal(:,2)  ,'ro','LineWidth',5,'DisplayName','DT')
%     plot(SIM.InputSignal2(:,1) ,SIM.InputSignal2(:,2) ,'bo','LineWidth',2,'DisplayName','DT Fix')
%     xlabel('Time [s]')
%     ylabel('Current [A m^{-2}]')
%     title('Compare ODE to DT Input Signal')
%     lgn = legend;
% end