function [AN,CA,SEP,EL,SIM,N,FLAG] = batt_inputs(SIM)
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

%% Simulation Modes
% 1) Polarization
% 2) Harmonic Perturbation
% 3) State Space EIS
% 4) Known BC Profile Controller
% 5) MOO Controller (MOO: Multi-Objective Optimization)
% 6) Simulink
% 7) Manual Current Profile

% 0) Files that are used strictly for data and not to run a simulation


%% Flags
FLAG.AN_LI_FOIL = 0; % 1 if the anode   is a Li foil
FLAG.CA_LI_FOIL = 0; % 1 if the cathode is a Li foil %%%%%Not IMplemented

FLAG.R_AN   = 1; % 1 if radial diffusion in anode   active material is considered
FLAG.R_CA   = 1; % 1 if radial diffusion in cathode active material is considered

FLAG.COE    = 0; % Cons of Energy (Temperature). 0 if dTdt = 0
    % T_BC
    % 1) Known Temperature  T(0,t)      = T_s
    % 2) Known Heat Flux    -k dTdx|x=0 = q''_s
    % 3) Insulated             dTdx|x=0 = 0
    % 4) Convection         -k dTdx|x=0 = h[T_inf - T(0,t)]
    FLAG.T_BC_AN = 4; 
    FLAG.T_BC_CA = 4; 

    FLAG.CV_CONV           = 0; % Convection is applied to the surface of each CV, not just the end CV %%%%%%%NOT IMPLEMENTED
    FLAG.HEAT_GEN_TOTAL    = 1; % Heat generation total within the control volume
    FLAG.HEAT_GEN_IONIC    = 1; % Heat generation due to ohmic losses from ionic current (el phase)
    FLAG.HEAT_GEN_ELECTRIC = 1; % Heat generation due to ohmic losses from electronic current (ed phase)
    FLAG.HEAT_GEN_RXN      = 1; % Heat generation from current across R_SEI %%% May not be accurate
    FLAG.HEAT_GEN_CHEM_RXN = 0; % Heat generation from a reaction at the SEI %%%%%%%NOT IMPLEMENTED
    
FLAG.V_SEI            = 1; % 1 if the overpotential is calculated using V_SEI
FLAG.SEI_growth       = 0; % 1 if the SEI grows over time (increased R_SEI) 
FLAG.HeatGen          = 0; % 1 if heat is generated by the i_BV across R_SEI

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
%     FLAG.VARIABLE_sigma    = 0;
    FLAG.VARIABLE_kappa    = 1;
    FLAG.VARIABLE_activity = 0;
    FLAG.VARIABLE_D_Liion  = 0;
    FLAG.VARIABLE_tf_num   = 0;
    FLAG.VARIABLE_D_o_AN   = 0; % 1 if the active material diffusion coefficient is dependent on concentration
    FLAG.VARIABLE_D_o_CA   = 0; % 1 if the active material diffusion coefficient is dependent on concentration

%%%% If CONSTANT_PROPS_FROM_HANDLES or VARIABLE_PROPS_FROM_HANDLES are 1
%%%% but none of the individual properties are a 1, then all properties are
%%%% equal to the defaults set below. 
if ( FLAG.CONSTANT_PROPS_FROM_HANDLES && FLAG.VARIABLE_PROPS_FROM_HANDLES)
    warning('Both FLAG.CONSTANT_PROPS_FROM_HANDLES and FLAG.VARIABLE_PROPS_FROM_HANDLES are 1.')
end

FLAG.SaveSolnDiscreteTime = 0; % 1 if evaluate the ode soln at a given sampling rate
	SIM.SaveTimeStep = 1;      % [s], Sampling rate of the ode output

FLAG.SaveSystemForEst = 0; % 1, if save the system to be used in the estimator ONLY FOR SS EIS (SimMode 3)

FLAG.doPostProcessing = 1;   % 1 if the postprocessing function is performed after a simulation completes
    FLAG.ReduceSolnTime = 0; % 1 if the results that are saved don't use all the points produced by t_soln ######NOT IMPLEMENTED YET
FLAG.Plot             = 0;   % 1 if the results plot immediately

FLAG.PRBS_predefinded = 0; % 1 if using a PRBS current that is predefined


%% Numerical Parameters
N.N_CV_AN   = 10;  % Number of x-direction anode     control volumes (CV) %%%Nodes now. Centered for now
N.N_CV_SEP  = 5;   % Number of x-direction seperator control volumes (CV)
N.N_CV_CA   = 10;  % Number of x-direction cathode   control volumes (CV)
N.N_R_AN    = 10;  % Number of r-direction anode     control volumes (CV) %%%No less than 10 should be used
N.N_R_CA    = 10;  % Number of r-direction cathode   control volumes (CV) %%%No less than 10 should be used

% N.N_CV_AN   = 3;  % Number of x-direction anode     control volumes (CV) %%%Nodes now. Centered for now
% N.N_CV_SEP  = 3;   % Number of x-direction seperator control volumes (CV)
% N.N_CV_CA   = 3;  % Number of x-direction cathode   control volumes (CV)
% N.N_R_AN    = 3;  % Number of r-direction anode     control volumes (CV) %%%No less than 10 should be used
% N.N_R_CA    = 3;  % Number of r-direction cathode   control volumes (CV) %%%No less than 10 should be used


%% fsolve options
options = optimoptions('fsolve');
options.FunctionTolerance = 1e-15;
% options.Display = 'iter';
options.Display = 'none';
SIM.fsolve_options = options;


%% Setting Simulation Parameters based on mode of operation
%% ---- Polarization ----
if SIM.SimMode == 1
    SIM.charge_frac         = 1.0;  % How deep do we want to charge/discharge? 
    SIM.t_ramp              = 0;    % [s], Ramp time for load current to go from 0 to i_user
%%%% Input from batch mode 
%     SIM.C_rate              = 2;	% How many charges per hour, ABSOLUTE VALUE
%     SIM.ChargeOrDischarge   = 1;   % -1 if Charge, 1 if Discharge 
%     SIM.SOC_start           = 50;    % [%], Initial state of charge of the cell 
end


%% ---- Harmonic Perturbation ----
if SIM.SimMode == 2
    SIM.C_rate    = 1/5;	  % How many charges per hour, ABSOLUTE VALUE
    SIM.N_cycles  = 10;   % [cycles], number of cycles or full oscillation of the sin wave
    SIM.t_ramp    = 0;
%%%% Input from batch mode 
%     SIM.freq      = 1e-2; % [rad/s], frequency of the sin wave 
%     SIM.SOC_start = 50;   % [%], Initial state of charge 
end


%% ---- State Space EIS ----
if SIM.SimMode == 3
    SIM.C_rate    = 1/20;  % How many charges per hour, ABSOLUTE VALUE
    SIM.t_ramp    = 0;
    
%%%% Input from batch mode 
%     SIM.freq      = 1e-2; % [rad/s], frequency of the sin wave 
%     SIM.SOC_start = 50;   % [%], Initial state of charge
end


%% ---- Known BC Profile Controller ----
if SIM.SimMode == 4
    SIM.DiscreteTimeStep    = 0.001;
    SIM.ControllerHandle    = @Controller_CV; % Function handle of the controller to call
    SIM.ZeroTime            = 0.01; % [s], How long a CV has to be at min i_user before simulation calls it quits 

%%%% Input from CreateProject
%     SIM.SOC_start           = 81.93;    % [%], Initial state of charge of the cell 81.93 ~ 4.0V
%%%% Input from controller
%     SIM.C_rate              = 2;	% How many charges per hour, ABSOLUTE VALUE
%     SIM.ChargeOrDischarge   = 1;   % -1 if Charge, 1 if Discharge 
%     SIM.charge_frac         = 1.0;  % How deep do we want to charge/discharge? 
%     SIM.t_ramp              = 0;    % [s], Ramp time for load current to go from 0 to i_user

end


%% ---- MOO Controller ----
if SIM.SimMode == 5
    SIM.SOC_start           = 81.93;    % [%], Initial state of charge of the cell 81.93 ~ 4.0V

%%%% Input from controller
%     SIM.C_rate              = 2;	% How many charges per hour, ABSOLUTE VALUE
%     SIM.ChargeOrDischarge   = 1;   % -1 if Charge, 1 if Discharge 
%     SIM.charge_frac         = 1.0;  % How deep do we want to charge/discharge? 
%     SIM.t_ramp              = 0;    % [s], Ramp time for load current to go from 0 to i_user

end


%% ---- Manual Profile ----
if SIM.SimMode == 7
    FLAG.Optimize_Profile       = 1; % 1 if the charge current profile is modified until it reaches the specified tolerance
        % ? Function handle if different types of things we are optimizing for?
        FLAG.Save_Current_Profile   = 1; % 1 if saving the final refined profile to its own .mat file
    
    SIM.increase_percent = 1.3; % If profile is being optimized for plating, this is measure for how much the profile increases
    SIM.decrease_percent = 2;   % If profile is being optimized for plating, this is measure for how much the profile decreases
    SIM.SOC_start        = 10;   % [%], Initial state of charge of the cell 
    SIM.C_rate           = 1;  % How many charges per hour, ABSOLUTE VALUE, This is used for initialization of i_user_amp used as a reference current
    SIM.C_rate_min       = 1/20;
    SIM.ramp_time        = 1; % [s],
    SIM.ChargeOrDischarge = -1; % -1 if Charge, 1 if Discharge 
    
    SIM.t_ramp    = 0;
    
%%%% Inputs loaded from makeCurrentProfile.m (SIM.profile_filepath,SIM)
% SIM.tol_Delta_phi
% SIM.max_iterations
% SIM.N_regions

end


%% Battery Chemistry (Cell Performance)
	% ---- Anode ----
    AN.EqPotentialHandle = @E_eqGraphite;
    AN.i_oHandle         = @i_oC6;
    AN.sigmaHandle       = @sigmaC6;
    AN.D_oHandle         = @D_o_Graphite;
    
    AN.k_o      = 3.80589129594505E-09;     % [A m^-2],       Exchange current density Rate constant
    % AN.i_o_a    = 0.5;      % [],             i_o exponent
    % AN.i_o_b    = 0.5;      % [],             i_o exponent
    % AN.i_o_c    = 0.5;      % [],             i_o exponent
    AN.alpha_a  = 0.5;        % [-],            Symmetry factor, annodic 
    AN.alpha_c  = 0.5;        % [-],            Symmetry factor, cathodic
    AN.C_dl     = 3e-6;      % [F/m^2],        Double-layer capacitance
    % AN.R_dl     = 0;        % [Ohm],          Double Layer Resistance
    AN.R_SEI    = 0.00733702042141719;     % [Ohm m^2],      Solid electrolyte interface resistance
    AN.sigma    = 100;              % [S m^-1],       Electrical conductivity (ed phase)
    AN.D_o      = 3E-13;            % [m^2 s^-1],     Solid-state diffusion coefficient
    AN.MW       = (72.0642)*1e0;    % [kg kmol^-1],   Molecular weight of C6
    AN.MW_Lith  = (79.0052)*1e0;    % [kg kmol^-1],   Molecular weight of LiC6
    AN.specCap  = (350)*1e0;        % [Ahr kg^-1],    Specific capacity
    AN.rho      = (2.266*1000)*1e0; % [kg m^-3],      Density of the electrode material without Li
    AN.c_p      = (720)*1e0;        % [J kg^-1 K^-1], Specific heat capacity  of separator material
    % AN.k        = 470;              % [W m^-1 K^-1],  Thermal conductivity
    AN.k        = (0.6)*1e0;        % [W m^-1 K^-1],  Thermal conductivity
    %sciencedirect.com/science/article/pii/S0735193317300179
    AN.C_Li_max    = 28.4851292093828;        % [kmol m^-3], Max concentration of lithium in the active material
    
        % Lithium Foil Properties
        if FLAG.AN_LI_FOIL
            AN.EqPotentialHandle = @E_eqLiFoil;
            AN.i_oHandle         = @i_oLiFoil;
            AN.sigmaHandle       = @sigmaC6;
            AN.D_oHandle         = @D_o_Graphite;
            AN.k_o      = 1E-07;     % [A m^-2],       Exchange current density Rate constant
            AN.alpha_a  = 0.3;        % [-],            Symmetry factor, annodic 
            AN.alpha_c  = 0.7;        % [-],            Symmetry factor, cathodic
            AN.C_dl     = 3e-6;       % [F/m^2],        Double-layer capacitance
            AN.R_SEI    = 1e-2;       % [Ohm m^2],      Solid electrolyte interface resistance
            AN.sigma    = 100;       % [S m^-1],       Electrical conductivity (ed phase)
            AN.D_o      = 3E-13;      % [m^2 s^-1],     Solid-state diffusion coefficient
            AN.MW       = 6.941;      % [kg kmol^-1],   Molecular weight of Li
            AN.MW_Lith  = 6.941;      % [kg kmol^-1],   Molecular weight of Li
            % AN.specCap  = (350)*1e0;  % [Ahr kg^-1],    Specific capacity
            AN.rho      = 0.534*1000; % [kg m^-3],      Density of the electrode material without Li
            AN.c_p      = 3000;       % [J kg^-1 K^-1], Specific heat capacity  of separator material
            AN.k        = 71.2;       % [W m^-1 K^-1],  Thermal conductivity
            AN.C_Li_max    = AN.rho/AN.MW; % [kmol m^-3], Max concentration of lithium in the active material
        end

    % ---- Separator ----
    SEP.rho = (2266)*1e0; % [kg m^-3],      Density of separator material
    SEP.c_p = (720)*1e0;  % [J kg^-1 K^-1], Specific heat capacity 
    % SEP.k   = 470;  % [W m^-1 K^-1], Thermal conductivity of separator material
    SEP.k   = (0.6)*1e0;  % [W m^-1 K^-1],  Thermal conductivity of separator material
    
    % ---- Cathode ----
    CA.EqPotentialHandle = @E_eqNMC;
    CA.i_oHandle         = @i_oC6;  
    CA.sigmaHandle       = @sigmaNMC;
    CA.D_oHandle         = @D_o_NMC532;
    
    CA.k_o      = 9.72997022033729E-10;   % [A m^-2],       Exchange current density Rate constant
    % CA.i_o_a    = 0.5;    % i_o exponent
    % CA.i_o_b    = 0.5;    % i_o exponent
    % CA.i_o_c    = 0.5;    % i_o exponent
    % CA.k_ct     = 1e-7;   
    CA.alpha_a  = 0.5;      % [-],            Symmetry factor, annodic 
    CA.alpha_c  = 0.5;      % [-],            Symmetry factor, cathodic
    CA.C_dl     = 1e-6;     % [F/m^2],        Double-layer capacitance
    % CA.R_dl     = 0.0005;   % [Ohm],          Double-layer resisitance
    CA.R_SEI    = 0.00611911635271628;   % [ohm m^2],      Solid electrolyte interface resistance
    CA.sigma    = 3.8;     % [S m^-1],       Electrical conductivity (am phase)
    CA.D_o      = 5E-13;    % [m^2 s^-1],     Solid-state diffusion coefficient
    CA.MW       = (89.61275)*1e0; % [kg kmol^-1],   Molecular weight of NMC
    CA.MW_Lith  = (96.55375)*1e0; % [kg kmol^-1],   Molecular weight of LNMC
%     CA.specCap  = 155;      % [Ahr kg^-1],    Specific capacity
    CA.specCap  = (170)*1e0;      % [Ahr kg^-1],    Specific capacity
%     CA.rho      = 2222.68;  % [kg m^-3],      Density of the electrode material without Li
    CA.rho      = (2739.6961)*1e0;  % [kg m^-3],      Density of the electrode material without Li
    CA.c_p      = (720)*1e0;      % [J kg^-1 K^-1], Specific heat capacity of cathode material
    % CA.k        = 470;      % [W m^-1 K^-1], Thermal conductivity of cathode material
    CA.k        = (0.6)*1e0;      % [W m^-1 K^-1],  Thermal conductivity of cathode material
    CA.C_Li_max    = 40.544149086751;        % [kmol m^-3], Max concentration of lithium in the active material
    
    % ---- Electrolyte ----
    EL.tf_numHandle       = @transferenceNumber;
    EL.ActivityHandle     = @activity;
    EL.D_o_Li_ionHandle   = @D_oLiion;
    EL.kappaHandle        = @kappa;
    
    EL.D_o_Li_ion   = 7.5E-11; % [m^2 s^-1], Li^+ liquid diffusion coefficient
    EL.kappa        = 0.28;      % [S m^-1], Ionic conductivity
    EL.C            = 1.0;      % [kmol m^-3], Li concentration
    EL.tf_num       = 0.363;    % [-], Transference Number 
    EL.Activity     = 1;        % 
    EL.rho          = 2266;     % [kg m^-3],      Density of the electrolyte %%%%%%%%%%Guess value
    EL.c_p          = 720;      % [J kg^-1 K^-1], Specific heat capacity of electrolyte
    EL.k            = 470;      % [W m^-1 K^-1],  Thermal conductivity of electrolyte
    

%% Cell Geometry
% c - coin
% p - pouch
% s - spiral

SIM.cell_geo = 'p';

    % ---- Anode ----
    if SIM.cell_geo == 'p'
        AN.del_z        = sqrt(0.1); % [m], Direction without tab, length of electrode
        AN.del_y        = sqrt(0.1); % [m], Direction with tab, height of electrode
    elseif SIM.cell_geo == 'c'
        AN.diam         = 1.5e-2;  % [m], coin cell diameter
    end
    AN.L            = 100e-6;   % [m], Direction normal to current flow, thickness of electrode
    AN.r_p          = 12.5E-06;    % [m], outer radius of electrode particle
    AN.eps_ed       = 0.320700992697794;   % [-], volume fraction of active material in the composite electrode
    AN.eps_b        = 0.172;     % [-], volume fraction of binder
    AN.A_geo        = 1.0;     % Account for ecentricity of spherical particle
    AN.gamma_brug   = 1.0;     % [-], Bruggeman pre-exponential multiplier
    AN.alpha_brug   = 1.5;     % [-], Bruggeman exponential factor
    % AN.coat_density = 0.06;    % [kg m^-2], Active material coating density
    % AN.SA_gram      = 4.2;     % [m^2/g], surface area of active material per gram of active material
    
        % Lithium Foil Properties
        if FLAG.AN_LI_FOIL
            AN.L            = 23.4e-6; % [m], Direction normal to current flow, thickness of electrode
            AN.r_p          = 50e-6;  % [m], outer radius of electrode particle
            AN.eps_ed       = 0.99;    % [-], volume fraction of active material in the composite electrode
            AN.eps_b        = 0.0;    % [-], volume fraction of binder
            AN.A_geo        = 1.0;    % Account for ecentricity of spherical particle
            AN.gamma_brug   = 1.0;    % [-], Bruggeman pre-exponential multiplier
            AN.alpha_brug   = 1.5;    % [-], Bruggeman exponential factor
            % AN.coat_density = 0.06;    % [kg m^-2], Active material coating density
            % AN.SA_gram      = 4.2;     % [m^2/g], surface area of active material per gram of active material
        end

    % ---- Separator ----
    if SIM.cell_geo == 'p'
        SEP.del_z        = sqrt(0.1); % [m], Direction without tab, length of electrode
        SEP.del_y        = sqrt(0.1); % [m], Direction with tab, height of electrode
    elseif SIM.cell_geo == 'c'
        SEP.diam         = 1.5e-2;  % [m], coin cell diameter
    end
    SEP.L           = 30e-6;   % [m] thickness
    SEP.eps         = 1;    % [-], Porosity (Volume fraction of electrolyte in sep)
%     SEP.gamma_brug   = 1.0;     % [-], Bruggeman pre-exponential multiplier
%     SEP.alpha_brug   = 1.5;     % [-], Bruggeman exponential factor

    % ---- Cathode ----
    if SIM.cell_geo == 'p'
        CA.del_z        = sqrt(0.1); % [m], Direction without tab, length of electrode
        CA.del_y        = sqrt(0.1); % [m], Direction with tab, height of electrode
    elseif SIM.cell_geo == 'c'
        CA.diam         = 1.5e-2;  % [m], coin cell diameter
    end
    CA.L            = 100e-6; % [m], Direction normal to current flow, thickness of electrode
    CA.r_p          = 5E-06;  % [m], outer radius of electrode particle
    CA.eps_ed       = 0.442994582876546; % [-], volume fraction of active material in the composite electrode
    CA.eps_b        = 0.259;   % [-], volume fraction of binder
    CA.A_geo        = 1.0;     % Account for ecentricity of spherical particle
    CA.gamma_brug   = 1.0;     % [-], Bruggeman pre-exponential multiplier
    CA.alpha_brug   = 1.5;     % [-], Bruggeman exponential factor
    % CA.coat_density = 0.06;  % [kg m^-2], Active material coating density
    % CA.SA_gram      = 4.2;     % [m^2/g], surface area of active material per gram of active material
    
    % ---- Electrolyte ----
    EL.gamma_brug   = 1.0;     % [-], Bruggeman pre-exponential multiplier
    EL.alpha_brug   = 1.5;     % [-], Bruggeman exponential factor


%% Overall Simulation Parameters
%%% Thermal 
SIM.Temp_start     = 25 + 273.15; % [K], Initial temperature of the cell 298.15
SIM.T_inf          = 25 + 273.15; % [K], Ambient temperature
SIM.h              = 1;           % [],  Convection coefficient on the ends of the battery (Forced vs natural)
if FLAG.T_BC_AN == 2 % Known heat flux
    SIM.q_AN_BC = 1; % [W m^-2]
end
if FLAG.T_BC_CA == 2 % Known heat flux
    SIM.q_CA_BC = 1; % [W m^-2]
end

% Simulation run time
SIM.initial_offset = 0;          % [s], How long there is an initial zero current

% Properties for SOC calcualtion
%%% Wiley
SIM.VoltageMax         = 4.2 ;  % [V] 
SIM.VoltageMin         = 3.4 ;  % [V] 
SIM.AnodeFormation_X   = 0.00;  % [-] 
SIM.CathodeFormation_X = 1.00;  % [-] 

%%% Wiley Half Cell NMC
% SIM.VoltageMax         = 5.4 ;  % [V] 
% SIM.VoltageMin         = 2.9 ;  % [V] 
% SIM.AnodeFormation_X   = 1.00;  % [-] 
% SIM.CathodeFormation_X = 0.00;  % [-] 
end