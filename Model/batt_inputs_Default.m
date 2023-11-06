%% batt_inputs
%
%% Battary Model Input File
    % This function is used to define the users inputs. User inputs are
    % * Numerical Parameters
    % * Simulation Modes
    % * Battery Chemistry
    % * Cell Geometry
    % * Overall Simulation Parameters
    %
    % Units used in this model
    % * Mass       : kilogram (kg)
    % * Moles      : kilomole (kmol)
    % * Time       : seconds  (s)
    % * Temperature: Kelvin   (K)
    % * Length     : meter    (m)
    % * Energy     : Joule    (J)
    % * Power      : Watt     (W)
    % * Current    : Amp      (A)
    % * Resistance : Ohms     (Omega)
    %             or Siemens  (S)
%
%% Simulation Modes
    % 1)  Polarization
    % 2)  Harmonic Perturbation
    % 3)  State Space EIS
    % 4)  Known BC Profile Controller
    % 5)  MOO Controller (MOO: Multi-Objective Optimization)
    % 6)  Simulink
    % 7)  Manual Current Profile
    % 8)  PRBS
    % 9)  EIS PRBS Stitching
    % 10) EIS Ho-Kalman
    % 0)  Files that are used strictly for data and not to run a simulation

function [AN,CA,SEP,EL,SIM,N,FLAG] = batt_inputs_Default(SIM)
%% Flags
    FLAG.AN_LI_FOIL = 0; % 1 if the anode   is a Li foil
    FLAG.CA_LI_FOIL = 0; % 1 if the cathode is a Li foil %%%%%Not Implemented
    
    FLAG.R_AN   = 1; % 1 if radial diffusion in anode   active material is considered
    FLAG.R_CA   = 1; % 1 if radial diffusion in cathode active material is considered

    % FLAG.OffDiagOnsager = 1;
    %     FLAG.Soret = 1;
    
    FLAG.COE    = 0; % Cons of Energy (Temperature). 0 if dTdt = 0
        % Thermal Boundary Conditions
        % 1) Known Temperature  T(0,t)      = T_s
        % 2) Known Heat Flux    -k dTdx|x=0 = q''_s
        % 3) Insulated             dTdx|x=0 = 0
        % 4) Convection         -k dTdx|x=0 = h[T_inf - T(0,t)]
            FLAG.T_BC_AN      = 3; % BC applied to the face normal to the AN current collector 
            FLAG.T_BC_CA      = 3; % BC applied to the face normal to the CA current collector 
            FLAG.T_BC_Outside = 3; % BC applied to the surface of each CV, Can only be 3 or 4
    
        FLAG.HEAT_GEN_TOTAL    = 0; % Heat generation total within the control volume
            FLAG.HEAT_GEN_IONIC    = 1; % Heat generation due to ohmic losses from ionic current (el phase)
            FLAG.HEAT_GEN_ELECTRIC = 1; % Heat generation due to ohmic losses from electronic current (ed phase)
            FLAG.HEAT_GEN_RXN      = 1; % Heat generation from current across R_SEI %%% May not be accurate
            % FLAG.HEAT_GEN_CHEM_RXN = 0; % Heat generation from a reaction at the SEI %%%%%%%NOT IMPLEMENTED

    FLAG.InitialThermalGradient = 0;
    FLAG.RampThermalGradient    = 0;
        SIM.RampThermalGradientTime = 3*60; % [s], time to ramp the BC
        % FLAG.TempBC
        %    0) Manual
        %    1) Iso20
        %    2) AN Cold
        %    3) CA Cold
        %    4) Iso18
        %    5) Iso22
        FLAG.TempBC = 5; % Pre-determined temperature BC
            
    FLAG.V_SEI            = 1; % 1 if the overpotential is calculated using V_SEI
    % FLAG.SEI_growth       = 0; % 1 if the SEI grows over time (increased R_SEI) %%%%%%%NOT IMPLEMENTED
    
    FLAG.Bruggeman = 0; % 1 if properties are adjusted for tortuosity
        FLAG.BRUG_ED       = 1; % Apply BRUG to electrode (active material) parameters
        FLAG.BRUG_EL       = 1; % Apply BRUG to electrolyte parameters
        % Specific Parameters to have BRUG applied to
            FLAG.BRUG_sigma    = 1; 
            FLAG.BRUG_kappa    = 1;
            FLAG.BRUG_activity = 0;
            FLAG.BRUG_D_Liion  = 1;
            FLAG.BRUG_tf_num   = 0;
            FLAG.BRUG_D_o_AN   = 0;
            FLAG.BRUG_D_o_CA   = 0;
    
    FLAG.Newman_i_o = 1; % 1 if using  i_o = k_o*F*(C_Li)^alpha_c * (C_max - C_Li)^alpha_a * (C_Liion)^alpha_a; 0 if using a function handle
    FLAG.preDefined_C_max = 1; % 1 if using the C_max found in the input file
    
    FLAG.CONSTANT_PROPS_FROM_HANDLES = 0; % Uses the initial conditions to solve for the properties and uses those throughout the simulation
    FLAG.VARIABLE_PROPS_FROM_HANDLES = 0;
        % FLAG.VARIABLE_sigma    = 0;
        FLAG.VARIABLE_kappa    = 1;
        FLAG.VARIABLE_activity = 0;
        FLAG.VARIABLE_D_Liion  = 1;
        FLAG.VARIABLE_tf_num   = 0;
        FLAG.VARIABLE_D_o_AN   = 1; % 1 if the active material diffusion coefficient is dependent on concentration
        FLAG.VARIABLE_D_o_CA   = 1; % 1 if the active material diffusion coefficient is dependent on concentration
    
    %%%% If CONSTANT_PROPS_FROM_HANDLES or VARIABLE_PROPS_FROM_HANDLES are 1
    %%%% but none of the individual properties are a 1, then all properties are
    %%%% equal to the defaults set below. 
    if ( FLAG.CONSTANT_PROPS_FROM_HANDLES && FLAG.VARIABLE_PROPS_FROM_HANDLES)
        warning('Both FLAG.CONSTANT_PROPS_FROM_HANDLES and FLAG.VARIABLE_PROPS_FROM_HANDLES are 1.')
    end
    
    % % Sampling Time
    %     N_t_s = 25;
    %     T_s_min = -1; % T_s = 10^(T_s_min)
    %     T_s_max =  1; % T_s = 10^(T_s_max)
    %     T_s_vec = logspace(T_s_min,T_s_max,N_t_s);
    %     Ts = T_s_vec(25);   % [s], discrete time sampling rate

    FLAG.AddInputNoise = 0;
        if FLAG.AddInputNoise
            SIM.Ts      = 1;    % [s], Sampling rate of the DT system
            SIM.Q_input = 1e-2; % [-], Input noise covariance matrix
        end
    
    FLAG.SaveSolnDiscreteTime = 0; % 1 if evaluate the ode soln at a given sampling rate
        if FLAG.AddInputNoise % Force DT save
            FLAG.SaveSolnDiscreteTime = 1;
        end
        if FLAG.AddInputNoise % Don't redefine the sampling rate and force to save
            % SIM.Ts           = Ts;                    % [s], Sampling rate of the DT system
        else
            SIM.Ts = 1; % [s], Sampling rate of the DT system
        end
	    SIM.TsMultiple   = 5;                     % Sample faster than desired SaveTimeStep, New SaveTimeStep is Ts/TsMultiple
        SIM.SaveTimeStep = SIM.Ts/SIM.TsMultiple; % [s], Sampling rate of the ode output

    if isfield(SIM,'t_Report') && SIM.SimMode == 9
        if ~isempty(SIM.t_Report)
            FLAG.SaveSolnDiscreteTime = 1;
                SIM.Ts           = SIM.Tswitch * SIM.t_Report;
                SIM.TsMultiple   = 1/SIM.t_Report;
                SIM.SaveTimeStep = SIM.Tswitch * SIM.t_Report;
        end
    end
    
    FLAG.SaveSystemForEst = 0; % 1, if save the system to be used in the estimator ONLY FOR SS EIS (SimMode 3)
    
    FLAG.doPostProcessing = 1;   % 1 if the postprocessing function is performed after a simulation completes
        FLAG.ReduceSolnTime = 0; % 1 if the results that are saved don't use all the points produced by t_soln ######NOT IMPLEMENTED YET
    FLAG.Plot             = 1;   % 1 if the results plot immediately
        FLAG.PlotUserCurrentProfiles = 0;


%% Numerical Parameters
    % N.N_CV_AN   = 3;  % Number of x-direction anode     control volumes (CV) %%%Nodes now. Centered for now
    % N.N_CV_SEP  = 3;  % Number of x-direction seperator control volumes (CV)
    % N.N_CV_CA   = 3;  % Number of x-direction cathode   control volumes (CV)
    % N.N_R_AN    = 3;  % Number of r-direction anode     control volumes (CV) %%%No less than 10 should be used
    % N.N_R_CA    = 3;  % Number of r-direction cathode   control volumes (CV) %%%No less than 10 should be used

    N.N_CV_AN   = 10;  % Number of x-direction anode     control volumes (CV) %%%Nodes now. Centered for now
    N.N_CV_SEP  = 5;   % Number of x-direction seperator control volumes (CV)
    N.N_CV_CA   = 10;  % Number of x-direction cathode   control volumes (CV)
    N.N_R_AN    = 10;  % Number of r-direction anode     control volumes (CV) %%%No less than 10 should be used
    N.N_R_CA    = 10;  % Number of r-direction cathode   control volumes (CV) %%%No less than 10 should be used

    % N.N_CV_AN   = 20;  % Number of x-direction anode     control volumes (CV) %%%Nodes now. Centered for now
    % N.N_CV_SEP  = 10;  % Number of x-direction seperator control volumes (CV)
    % N.N_CV_CA   = 20;  % Number of x-direction cathode   control volumes (CV)
    % N.N_R_AN    = 10;  % Number of r-direction anode     control volumes (CV) %%%No less than 10 should be used
    % N.N_R_CA    = 10;  % Number of r-direction cathode   control volumes (CV) %%%No less than 10 should be used
    

%% fsolve options
    options                   = optimoptions('fsolve');
    options.FunctionTolerance = 1e-15;
    % options.Display           = 'iter';
    options.Display           = 'none';
    SIM.fsolve_options        = options;


%% Setting Simulation Parameters based on mode of operation
%% ---- Polarization ----
    if SIM.SimMode == 1
        SIM.charge_frac         = 1.0;  % How deep do we want to charge/discharge? 
        SIM.t_ramp              = 0;    % [s], Ramp time for load current to go from 0 to i_user
    %%%% Input from batch mode 
        % SIM.C_rate              = 2;	% How many charges per hour, ABSOLUTE VALUE
        % SIM.ChargeOrDischarge   = 1;   % -1 if Charge, 1 if Discharge 
        % SIM.SOC_start           = 50;    % [%], Initial state of charge of the cell 
    end


%% ---- Harmonic Perturbation ----
    if SIM.SimMode == 2
        SIM.C_rate    = 1/5;	  % How many charges per hour, ABSOLUTE VALUE
        SIM.N_cycles  = 10;   % [cycles], number of cycles or full oscillation of the sin wave
        SIM.t_ramp    = 0;
    %%%% Input from batch mode 
        % SIM.freq      = 1e-2; % [rad/s], frequency of the sin wave 
        % SIM.SOC_start = 50;   % [%], Initial state of charge 
    end


%% ---- State Space EIS ----
    if SIM.SimMode == 3
        SIM.C_rate    = 1/20;  % How many charges per hour, ABSOLUTE VALUE
        SIM.t_ramp    = 0;
        
    %%%% Input from batch mode 
        % SIM.omega      = 1e-2; % [rad/s], Frequency of the sin wave 
        % SIM.SOC_start = 50;    % [%],     Initial state of charge
    end


%% ---- Known BC Profile Controller ----
    if SIM.SimMode == 4
        SIM.DiscreteTimeStep    = 0.001; %%%%%%%%%%% IDK if this is ever used
        SIM.ControllerHandle    = @Controller_CV; % Function handle of the controller to call
        SIM.ZeroTime            = 0.01; % [s], How long a CV has to be at min i_user before simulation calls it quits 
    
    %%%% Input from CreateProject
        % SIM.SOC_start           = 81.93;    % [%], Initial state of charge of the cell 81.93 ~ 4.0V
    %%%% Input from controller
        % SIM.C_rate              = 2;	% How many charges per hour, ABSOLUTE VALUE
        % SIM.ChargeOrDischarge   = 1;   % -1 if Charge, 1 if Discharge 
        % SIM.charge_frac         = 1.0;  % How deep do we want to charge/discharge? 
        % SIM.t_ramp              = 0;    % [s], Ramp time for load current to go from 0 to i_user
    
    end


%% ---- MOO Controller ----
    if SIM.SimMode == 5
        SIM.SOC_start           = 81.93;    % [%], Initial state of charge of the cell 81.93 ~ 4.0V
    
    %%%% Input from controller
        % SIM.C_rate              = 2;	% How many charges per hour, ABSOLUTE VALUE
        % SIM.ChargeOrDischarge   = 1;   % -1 if Charge, 1 if Discharge 
        % SIM.charge_frac         = 1.0;  % How deep do we want to charge/discharge? 
        % SIM.t_ramp              = 0;    % [s], Ramp time for load current to go from 0 to i_user
    
    end


%% ---- Manual Profile ----
    if SIM.SimMode == 7
        FLAG.Optimize_Profile       = 1; % 1 if the charge current profile is modified until it reaches the specified tolerance
            % ? Function handle if different types of things we are optimizing for?
            FLAG.Save_Current_Profile   = 1; % 1 if saving the final refined profile to its own .mat file
        
        SIM.increase_percent  = 1.3; % If profile is being optimized for plating, this is measure for how much the profile increases
        SIM.decrease_percent  = 2;   % If profile is being optimized for plating, this is measure for how much the profile decreases
        SIM.SOC_start         = 10;   % [%], Initial state of charge of the cell 
        SIM.C_rate            = 1;  % How many charges per hour, ABSOLUTE VALUE, This is used for initialization of i_user_amp used as a reference current
        SIM.C_rate_min        = 1/20;
        SIM.ramp_time         = 1; % [s],
        SIM.ChargeOrDischarge = -1; % -1 if Charge, 1 if Discharge 
        
        SIM.t_ramp    = 0;
        
    %%%% Inputs loaded from makeCurrentProfile.m (SIM.profile_filepath,SIM)
        % SIM.tol_Delta_phi
        % SIM.max_iterations
        % SIM.N_regions
    
    end


%% ---- PRBS ----
    if SIM.SimMode == 8
        %SIM.PRBSLength     = 100;
        %SIM.PRBSLength     = 77; % New PRBS
        SIM.PRBSLength     = 159; % New PRBS, at this length, the signal ends at a zero mean current
        %SIM.t_ramp_ratio   = 1/20;        % Fraction of the switching time that is used as ramp interpolation
        SIM.t_ramp_ratio   = 1/50;        % Fraction of the switching time that is used as ramp interpolation
        SIM.initial_offset = SIM.Tswitch; % Initial zero offset
        % SIM.initial_offset = 0; % Initial zero offset
    %%%% Input from batch mode 
        % SIM.SOC_start = 50;   % [%], Initial state of charge 
        % SIM.PRBSAmp
        % SIM.Tswitch
    end


%% ---- EIS from Stiching PRBS ----
    if SIM.SimMode == 9
    %     %SIM.PRBSLength     = 100;
    %     %SIM.PRBSLength     = 77; % New PRBS
    %     SIM.PRBSLength     = 159; % New PRBS
    %     %SIM.t_ramp_ratio   = 1/20;        % Fraction of the switching time that is used as ramp interpolation
    %     SIM.t_ramp_ratio   = 1/50;        % Fraction of the switching time that is used as ramp interpolation
    %     SIM.initial_offset = SIM.Tswitch; % Initial zero offset
    % %%%% Input from batch mode 
    % %     SIM.SOC_start = 50;   % [%], Initial state of charge 
    % %     SIM.PRBSAmp
    % %     SIM.Tswitch
    end


%% ---- EIS Ho-Kalman ----
    if SIM.SimMode == 10
        SIM.t_ramp_ratio   = 1/50;        % Fraction of the switching time that is used as ramp interpolation
        
    %%%% Input from batch mode 
        % SIM.SOC_start   = 50;     % [%],     Initial state of charge
        % SIM.Tsample     = HK_Ts;  % [s],     Sampling Time
        % SIM.freq        = 1e-2;   % [rad/s], Frequency 
        % SIM.HK_nSamples = HK_nSamples; % [], Number of Relax Samples
    end


%% Battery Chemistry (Cell Performance)
% ---- Anode ----
    AN.EqPotentialHandle = @E_eqGraphite;
    AN.i_oHandle         = @i_oC6;
    AN.sigmaHandle       = @sigmaC6;
    AN.D_oHandle         = @D_o_Graphite;
    
    AN.k_o      = 3.80589129594505E-09; % [A m^-2],       Exchange current density Rate constant
    AN.alpha_a  = 0.5;                  % [-],            Symmetry factor, annodic 
    AN.alpha_c  = 0.5;                  % [-],            Symmetry factor, cathodic
    AN.C_dl     = 3e-6;                 % [F m^-2],       Double-layer capacitance
    AN.R_SEI    = 0.00733702042141719;  % [Ohm m^2],      Solid electrolyte interface resistance
    AN.sigma    = 100;                  % [S m^-1],       Electrical conductivity (ed phase)
    AN.D_o      = 3E-13;                % [m^2 s^-1],     Solid-state diffusion coefficient
    AN.MW       = (72.0642)*1e0;        % [kg kmol^-1],   Molecular weight of C6
    AN.MW_Lith  = (79.0052)*1e0;        % [kg kmol^-1],   Molecular weight of LiC6
    AN.specCap  = (350)*1e0;            % [Ahr kg^-1],    Specific capacity
    % AN.rho      = (2.266*1000)*1e0;     % [kg m^-3],      Density of the electrode material without Li
    % AN.c_p      = (720)*1e0;            % [J kg^-1 K^-1], Specific heat capacity of active material
    % % AN.lambda = 470;                    % [W m^-1 K^-1],  Thermal conductivity
    % AN.lambda   = 1.04;                 % [W m^-1 K^-1],  Thermal conductivity of active material %sciencedirect.com/science/article/pii/S0735193317300179        
    AN.C_Li_max = 28.544149086751;      % [kmol m^-3],    Max concentration of lithium in the active material
    AN.rho      = (2.266*1000)*1e0; % [kg m^-3],      Testing
    AN.c_p      = (720)*1e0;        % [J kg^-1 K^-1], Testing
    AN.lambda   = 1.04*1e-3;        % [W m^-1 K^-1],  Testing

    % Lithium Foil Properties
        if FLAG.AN_LI_FOIL
            AN.EqPotentialHandle = @E_eqLiFoil;
            AN.i_oHandle         = @i_oLiFoil;
            AN.sigmaHandle       = @sigmaC6;
            AN.D_oHandle         = @D_o_Graphite;
            AN.k_o      = 1E-07;        % [A m^-2],       Exchange current density Rate constant
            AN.alpha_a  = 0.3;          % [-],            Symmetry factor, annodic 
            AN.alpha_c  = 0.7;          % [-],            Symmetry factor, cathodic
            AN.C_dl     = 3e-6;         % [F/m^2],        Double-layer capacitance
            AN.R_SEI    = 1e-2;         % [Ohm m^2],      Solid electrolyte interface resistance
            AN.sigma    = 100;          % [S m^-1],       Electrical conductivity (ed phase)
            AN.D_o      = 3E-13;        % [m^2 s^-1],     Solid-state diffusion coefficient
            AN.rho      = 0.534*1000;   % [kg m^-3],      Density of the electrode material without Li
            AN.c_p      = 3000;         % [J kg^-1 K^-1], Specific heat capacity  of separator material
            AN.lambda   = 71.2;         % [W m^-1 K^-1],  Thermal conductivity
            AN.C_Li_max = AN.rho/AN.MW; % [kmol m^-3], Max concentration of lithium in the active material
        end

% ---- Separator ----
    % SEP.rho    = (2266)*1e0; % [kg m^-3],      Density of separator material
    % SEP.c_p    = (720)*1e0;  % [J kg^-1 K^-1], Specific heat capacity 
    % % SEP.lambda = 470;        % [W m^-1 K^-1],  Thermal conductivity of separator material
    % SEP.lambda = (0.6)*1e0;  % [W m^-1 K^-1],  Thermal conductivity of separator material
    SEP.rho      = (2.266*1000)*1e0; % [kg m^-3],      Testing
    SEP.c_p      = (720)*1e0;        % [J kg^-1 K^-1], Testing
    SEP.lambda   = 1.04*1e-3;        % [W m^-1 K^-1],  Testing

% ---- Cathode ----
    CA.EqPotentialHandle = @E_eqNMC;
    CA.i_oHandle         = @i_oC6;  
    CA.sigmaHandle       = @sigmaNMC;
    CA.D_oHandle         = @D_o_NMC532;

    CA.k_o      = 9.72997022033729E-10; % [A m^-2],       Exchange current density Rate constant
    CA.alpha_a  = 0.5;                  % [-],            Symmetry factor, annodic 
    CA.alpha_c  = 0.5;                  % [-],            Symmetry factor, cathodic
    CA.C_dl     = 1e-6;                 % [F m^-2],       Double-layer capacitance
    CA.R_SEI    = 0.00611911635271628;  % [ohm m^2],      Solid electrolyte interface resistance
    CA.sigma    = 3.8;                  % [S m^-1],       Electrical conductivity (am phase)
    CA.D_o      = 5E-13;                % [m^2 s^-1],     Solid-state diffusion coefficient
    % CA.rho      = (2739.6961)*1e0;      % [kg m^-3],      Density of the electrode material without Li
    % CA.c_p      = (720)*1e0;            % [J kg^-1 K^-1], Specific heat capacity of cathode material
    % CA.lambda   = (0.6)*1e0;            % [W m^-1 K^-1],  Thermal conductivity of cathode material
    CA.C_Li_max = 40.544149086751;      % [kmol m^-3],    Max concentration of lithium in the active material
    CA.rho      = (2.266*1000)*1e0; % [kg m^-3],      Testing
    CA.c_p      = (720)*1e0;        % [J kg^-1 K^-1], Testing
    CA.lambda   = 1.04*1e-3;        % [W m^-1 K^-1],  Testing
    

% ---- Electrolyte ----
    EL.tf_numHandle       = @transferenceNumber;
    EL.ActivityHandle     = @activity;
    EL.D_o_Li_ionHandle   = @D_oLiion;
    EL.kappaHandle        = @kappa;
    
    EL.D_o_Li_ion = 7.5E-11; % [m^2 s^-1],     Li^+ liquid diffusion coefficient
    EL.kappa      = 0.28;    % [S m^-1],       Ionic conductivity
    EL.C          = 1.0;     % [kmol m^-3],    Li concentration
    EL.tf_num     = 0.363;   % [-],            Transference Number 
    EL.Activity   = 1;       % [-],
    % EL.rho        = 2266;    % [kg m^-3],      Density of the electrolyte %%%%%%%%%%Guess value
    % EL.c_p        = 720;     % [J kg^-1 K^-1], Specific heat capacity of electrolyte
    % EL.lambda     = 0.45;    % [W m^-1 K^-1],  Thermal conductivity of electrolyte %sciencedirect.com/science/article/pii/S0735193317300179
    EL.rho      = (2.266*1000)*1e0; % [kg m^-3],      Testing
    EL.c_p      = (720)*1e0;        % [J kg^-1 K^-1], Testing
    EL.lambda   = 1.04*1e-3;        % [W m^-1 K^-1],  Testing
    

%% Cell Geometry
    % c - coin
    % p - pouch
    % s - spiral

    SIM.cell_geo = 'p';

    % ---- Anode ----
        if SIM.cell_geo == 'p'
            AN.del_z        = sqrt(0.1);     % [m], Direction without tab, length of electrode
            AN.del_y        = sqrt(0.1);     % [m], Direction with    tab, height of electrode
        elseif SIM.cell_geo == 'c'
            AN.diam         = 1.5e-2;        % [m], Coin cell diameter
        end
        AN.L            = 100e-6;            % [m], Direction normal to current flow, thickness of electrode
        AN.r_p          = 12.5E-06;          % [m], Outer radius of electrode particle
        AN.eps_ed       = 0.320700992697794; % [-], Volume fraction of active material in the composite electrode
        AN.eps_b        = 0.172;             % [-], Volume fraction of binder
        AN.A_geo        = 1.0;               %      Account for ecentricity of spherical particle
        AN.gamma_brug   = 1.0;               % [-], Bruggeman pre-exponential multiplier
        AN.alpha_brug   = 1.5;               % [-], Bruggeman exponential factor
        
        % Lithium Foil Properties
        if FLAG.AN_LI_FOIL
            AN.L            = 23.4e-6;       % [m], Direction normal to current flow, thickness of electrode
            AN.r_p          = 50e-6;         % [m], Outer radius of electrode particle
            AN.eps_ed       = 0.99;          % [-], Volume fraction of active material in the composite electrode
            AN.eps_b        = 0.0;           % [-], Volume fraction of binder
            AN.A_geo        = 1.0;           %      Account for ecentricity of spherical particle
            AN.gamma_brug   = 1.0;           % [-], Bruggeman pre-exponential multiplier
            AN.alpha_brug   = 1.5;           % [-], Bruggeman exponential factor
        end

    % ---- Separator ----
        if SIM.cell_geo == 'p'
            SEP.del_z        = sqrt(0.1);    % [m], Direction without tab, length of electrode
            SEP.del_y        = sqrt(0.1);    % [m], Direction with tab, height of electrode
        elseif SIM.cell_geo == 'c'
            SEP.diam         = 1.5e-2;       % [m], coin cell diameter
        end
        SEP.L           = 30e-6;             % [m]  Thickness
        SEP.eps         = 1;                 % [-], Porosity (Volume fraction of electrolyte in sep)
        % SEP.gamma_brug   = 1.0;              % [-], Bruggeman pre-exponential multiplier
        % SEP.alpha_brug   = 1.5;              % [-], Bruggeman exponential factor

    % ---- Cathode ----
        if SIM.cell_geo == 'p'
            CA.del_z        = sqrt(0.1);     % [m], Direction without tab, length of electrode
            CA.del_y        = sqrt(0.1);     % [m], Direction with    tab, height of electrode
        elseif SIM.cell_geo == 'c'
            CA.diam         = 1.5e-2;        % [m], Coin cell diameter
        end
        CA.L            = 100e-6;            % [m], Direction normal to current flow, thickness of electrode
        CA.r_p          = 5E-06;             % [m], Outer radius of electrode particle
        CA.eps_ed       = 0.442994582876546; % [-], Volume fraction of active material in the composite electrode
        CA.eps_b        = 0.259;             % [-], Volume fraction of binder
        CA.A_geo        = 1.0;               %      Account for ecentricity of spherical particle
        CA.gamma_brug   = 1.0;               % [-], Bruggeman pre-exponential multiplier
        CA.alpha_brug   = 1.5;               % [-], Bruggeman exponential factor
        
    % ---- Electrolyte ----
        EL.gamma_brug     = 1.0;               % [-], Bruggeman pre-exponential multiplier
        %EL.alpha_brug   = 1.5;               % [-], Bruggeman exponential factor
        EL.alpha_brug_an  = 2.0; % [-], Bruggeman exponential factor
        EL.alpha_brug_sep = 2.0; % [-], Bruggeman exponential factor
        EL.alpha_brug_ca  = 2.0; % [-], Bruggeman exponential factor


%% Overall Simulation Parameters
% Thermal 
    % Ambient Temperature
    switch FLAG.TempBC
        case 0 % Manual
            SIM.T_inf = 20 + 273.15; % [K], Ambient temperature
        case 1 % Isothermal20
            SIM.T_inf = 20 + 273.15; % [K], Ambient temperature
        case 2 % Anode Cold
            SIM.T_inf = 20 + 273.15; % [K], Ambient temperature
        case 3 % Cathode Cold
            SIM.T_inf = 20 + 273.15; % [K], Ambient temperature
        case 4 % Isothermal18
            SIM.T_inf = 18 + 273.15; % [K], Ambient temperature
        case 5 % Isothermal22
            SIM.T_inf = 22 + 273.15; % [K], Ambient temperature
    end

    % Cell Initial Temperature Profile
    if FLAG.InitialThermalGradient
        switch FLAG.TempBC
            case 0 % Manual
                SIM.AN_Temp_init = 22 + 273.15; % [K], Anode   temperature
                SIM.CA_Temp_init = 18 + 273.15; % [K], Cathode temperature
            case 1 % Isothermal20
                SIM.AN_Temp_init = 20 + 273.15; % [K], Anode   temperature
                SIM.CA_Temp_init = 20 + 273.15; % [K], Cathode temperature
            case 2 % Anode Cold
                SIM.AN_Temp_init = 18 + 273.15; % [K], Anode   temperature
                SIM.CA_Temp_init = 22 + 273.15; % [K], Cathode temperature
            case 3 % Cathode Cold
                SIM.AN_Temp_init = 22 + 273.15; % [K], Anode   temperature
                SIM.CA_Temp_init = 18 + 273.15; % [K], Cathode temperature
            case 4 % Isothermal18
                SIM.AN_Temp_init = 18 + 273.15; % [K], Anode   temperature
                SIM.CA_Temp_init = 18 + 273.15; % [K], Cathode temperature
            case 5 % Isothermal22
                SIM.AN_Temp_init = 22 + 273.15; % [K], Anode   temperature
                SIM.CA_Temp_init = 22 + 273.15; % [K], Cathode temperature
        end
    else
        switch FLAG.TempBC
            case 0 % Manual
                SIM.Temp_init = 20 + 273.15; % [K], Initial temperature of the cell
            case 1 % Isothermal20
                SIM.Temp_init = 20 + 273.15; % [K], Initial temperature of the cell
            case 4 % Isothermal18
                SIM.Temp_init = 18 + 273.15; % [K], Initial temperature of the cell
            case 5 % Isothermal22
                SIM.Temp_init = 22 + 273.15; % [K], Initial temperature of the cell
        end
    end
    
    % Outside (faces normal to y,z direction) Thermal Boundary Condition
    switch FLAG.T_BC_Outside
        case 3 % Insulated
            SIM.h = 0; % [W m^-2 K^-1],  Convection coefficient on the ends of the battery (Forced vs natural)
        case 4 % Convection
            SIM.h = 4; % [W m^-2 K^-1],  Convection coefficient on the ends of the battery (Forced vs natural)
    end
    
    % Anode Thermal Boundary Condition
    switch FLAG.T_BC_AN 
        case 1 % Known Temperature [K]
            switch FLAG.TempBC
                case 0 % Manual
                    SIM.Temp_AN_BC = 25 + 273.15;
                case 1 % Isothermal20
                    SIM.Temp_AN_BC = 20 + 273.15;
                case 2 % Anode Cold
                    SIM.Temp_AN_BC = 18 + 273.15;
                case 3 % Cathode Cold
                    SIM.Temp_AN_BC = 22 + 273.15;
                case 4 % Isothermal18
                    SIM.Temp_AN_BC = 18 + 273.15;
                case 5 % Isothermal22
                    SIM.Temp_AN_BC = 22 + 273.15;
            end
        case 2 % Known heat flux
            SIM.q_AN_BC = -2; % [W m^-2]
        case 3 % Insulated
            SIM.q_AN_BC = 0; % [W m^-2]
        case 4 % Convection
            SIM.h_AN_BC = 4; % [W m^-2 K^-1]
    end

    % Cathode Thermal Boundary Condition
    switch FLAG.T_BC_CA 
        case 1 % Known Temperature [K]
            switch FLAG.TempBC
                case 0 % Manual
                    SIM.Temp_CA_BC = 25 + 273.15;
                case 1 % Isothermal20
                    SIM.Temp_CA_BC = 20 + 273.15;
                case 2 % Anode Cold
                    SIM.Temp_CA_BC = 22 + 273.15;
                case 3 % Cathode Cold
                    SIM.Temp_CA_BC = 18 + 273.15;
                case 4 % Isothermal18
                    SIM.Temp_CA_BC = 18 + 273.15;
                case 5 % Isothermal22
                    SIM.Temp_CA_BC = 22 + 273.15;
            end
        case 2 % Known heat flux
            SIM.q_CA_BC = -1; % [W m^-2]
        case 3 % Insulated
            SIM.q_CA_BC = 0; % [W m^-2]
        case 4 % Convection
            SIM.h_CA_BC = 4; % [W m^-2 K^-1]
    end

% Simulation run time
    if SIM.SimMode == 10 % EIS Ho-Kalman
        SIM.initial_offset = SIM.Tsample; % [s], How long there is an initial zero current
    elseif SIM.SimMode ~= 8 % ~PRBS
        if FLAG.RampThermalGradient
            SIM.initial_offset = SIM.RampThermalGradientTime*2; % [s], How long there is an initial zero current
        else
            SIM.initial_offset = 0;           % [s], How long there is an initial zero current
        end     
    end

% Properties for SOC calcualtion
    if SIM.SimMode ~= 8 % ~PRBS
        % Properties for SOC calcualtion
        %%% Wiley
        SIM.VoltageMax         = 4.2 ;  % [V] 
        SIM.VoltageMin         = 3.4 ;  % [V] 
        SIM.AnodeFormation_X   = 0.00;  % [-] 
        SIM.CathodeFormation_X = 1.00;  % [-] 

        % SIM.OneC_measured      = 22.7; % [A/m^2], Measured cap (used for demand)    (2.14 mAh cm^-2)
        % SIM.VoltageMax         = [] ;   % [V] 
        % SIM.VoltageMin         = [] ;   % [V] 
        % SIM.AnodeFormation_X   = []; % [-] 
        % SIM.CathodeFormation_X = []; % [-] 
        % 
        % SIM.AnodeStoich_SOC0     = 0.0700; % [-] 
        % SIM.CathodeStoich_SOC0   = 0.8900; % [-] 
        % SIM.AnodeStoich_SOC100   = 0.8434; % [-] 
        % SIM.CathodeStoich_SOC100 = 0.3400; % [-]
        
        %%% Wiley Half Cell NMC
        % SIM.VoltageMax         = 5.4 ;  % [V] 
        % SIM.VoltageMin         = 2.9 ;  % [V] 
        % SIM.AnodeFormation_X   = 1.00;  % [-] 
        % SIM.CathodeFormation_X = 0.00;  % [-] 
    else
        SIM.VoltageMax         = 5 ;  % [V] 
        SIM.VoltageMin         = 2 ;  % [V] 
        SIM.AnodeFormation_X   = 0.00;  % [-] 
        SIM.CathodeFormation_X = 1.00;  % [-] 
    end

end


%% OOOLLLDDD Properties
% AN.i_o_a    = 0.5;      % [],             i_o exponent
    % AN.i_o_b    = 0.5;      % [],             i_o exponent
    % AN.i_o_c    = 0.5;      % [],             i_o exponent
    % AN.R_dl     = 0;        % [Ohm],          Double Layer Resistance
    % AN.MW       = 6.941;      % [kg kmol^-1],   Molecular weight of Li
    % AN.MW_Lith  = 6.941;      % [kg kmol^-1],   Molecular weight of Li
    % AN.specCap  = (350)*1e0;  % [Ahr kg^-1],    Specific capacity
            
    % CA.i_o_a    = 0.5;    % i_o exponent
    % CA.i_o_b    = 0.5;    % i_o exponent
    % CA.i_o_c    = 0.5;    % i_o exponent
    % CA.k_ct     = 1e-7;   
    % CA.R_dl     = 0.0005;   % [Ohm],          Double-layer resisitance
    % CA.MW       = (89.61275)*1e0; % [kg kmol^-1],   Molecular weight of NMC
    % CA.MW_Lith  = (96.55375)*1e0; % [kg kmol^-1],   Molecular weight of LNMC
    % CA.specCap  = (170)*1e0;      % [Ahr kg^-1],    Specific capacity
    
    % AN.coat_density = 0.06;    % [kg m^-2], Active material coating density
    % AN.SA_gram      = 4.2;     % [m^2/g], surface area of active material per gram of active material
    
    % CA.coat_density = 0.06;  % [kg m^-2], Active material coating density
    % CA.SA_gram      = 4.2;     % [m^2/g], surface area of active material per gram of active material
    






    % EL.rho      = (2.266*1000)*1e0; % [kg m^-3],      Testing
    % EL.c_p      = (720)*1e0;        % [J kg^-1 K^-1], Testing
    % EL.lambda   = 1.04*1e-5;        % [W m^-1 K^-1],  Testing
    