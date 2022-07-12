function [AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS] = batt_init(AN,CA,SEP,EL,SIM,N,FLAG)
%% Batt Initialization
% This function is used to initialize all of the user inputs and calculates
% the necessary parameters and pointers for the model

%% Constants
    CONS.F = 96485338.3; % [C kmol^-1], Faraday's Constant 
    CONS.R = 8314.472;   % [J kmol^-1 K^-1], Gas Constant
    if FLAG.COE == 0
        CONS.T = SIM.Temp_start;     % [K], Temperature
    end

%% Control Volume Modification for Li Foil
% Adjust Number of CV
    if FLAG.AN_LI_FOIL
        N.N_CV_AN = 1; % Single control volume in the x-dir
        FLAG.R_AN = 0; % 1 if radial gradients are considered
        N.N_R_AN  = 1; % Number of radial control volumes
    end
    if FLAG.CA_LI_FOIL
        N.N_CV_CA = 1; % Single control volume in the x-dir
        FLAG.R_CA = 0; % 1 if radial gradients are considered     
        N.N_R_CA  = 1; % Number of radial control volumes
    end

%% Control Volume Modification for Distributed
% If radial concentration gradients are not considered, then set the number
% of radial control volumes equal to 1
if ~FLAG.R_AN
    N.N_R_AN  = 1;
end
if ~FLAG.R_CA
    N.N_R_CA  = 1;
end 

N.N_R_max = max(N.N_R_AN,N.N_R_CA);

%% Pointers for Tracked Variables
    i = 1;
    P.T         = i; i = i + 1;
    P.phi_el    = i; i = i + 1;
    P.phi_ed    = i; i = i + 1;
    P.V_1       = i; i = i + 1;
    P.V_2       = i; i = i + 1;
    P.i_PS      = i; i = i + 1;
    P.C_Liion   = i; i = i + 1;
    P.C_Li      = i; i = i + 1;
    
	P.C_Li_surf_AN  = P.C_Li + N.N_R_AN - 1;
	P.C_Li_surf_CA  = P.C_Li + N.N_R_CA - 1;
	P.C_Li_surf_max = max( P.C_Li_surf_AN , P.C_Li_surf_CA );
    
% Pointers for SEP region
    i = 1;
    P.SEP.T       = i; i = i + 1;
    P.SEP.phi_el  = i; i = i + 1;
    P.SEP.C_Liion = i; i = i + 1;
    
% Pointers for Electrostatic 
    i = 1;
    P.ES.phi_el = i; i = i + 1;
    P.ES.phi_ed = i; i = i + 1;
    P.ES.V_1    = i; i = i + 1;
    P.ES.V_2    = i; i = i + 1;
    P.ES.i_PS   = i; i = i + 1;
    P.ES.i_dl   = i; i = i + 1;
    
    N.N_ES_var  = 6;

% Pointers for State Space EIS
if SIM.SimMode == 3
    i = 1;
    P.SS.omega    = i; i = i + 1;
    P.SS.Z_mag    = i; i = i + 1;
    P.SS.Z_Re     = i; i = i + 1;
    P.SS.Z_Im     = i; i = i + 1;
    P.SS.Z_dB     = i; i = i + 1;
    P.SS.Z_ps_deg = i; i = i + 1;
end

% Number of nonRadial SV
    N.N_SV_nR  = 7;    
    N.N_SV_SEP = 3;
    
% Pointers for properties
    i = 1;
    P.sigma      = i; i = i + 1;
    P.kappa      = i; i = i + 1;
    P.R_SEI      = i; i = i + 1;
    P.k          = i; i = i + 1;
    P.rho        = i; i = i + 1;
    P.c_p        = i; i = i + 1;
    P.D_o_Li_ion = i; i = i + 1;
    P.activity   = i; i = i + 1;
    P.tf_num     = i; i = i + 1;
    P.D_o        = i; i = i + 1;
    
    N.N_prop = P.D_o + N.N_R_max - 1;

%% Region Indexing
% Number of state variables in each CV
    N.N_SV_AN  = N.N_SV_nR + N.N_R_AN;
    N.N_SV_SEP = N.N_SV_SEP;
    N.N_SV_CA  = N.N_SV_nR + N.N_R_CA;

    N.N_SV_AN_tot  = N.N_CV_AN  * N.N_SV_AN;
    N.N_SV_SEP_tot = N.N_CV_SEP * N.N_SV_SEP;
    N.N_SV_CA_tot  = N.N_CV_CA  * N.N_SV_CA;

    N.N_SV_tot = N.N_SV_AN_tot + N.N_SV_SEP_tot + N.N_SV_CA_tot;
    N.N_CV_tot = N.N_CV_AN + N.N_CV_SEP + N.N_CV_CA;

% Regional Indexing 
    N.CV_Region_AN  =                          1 : N.N_CV_AN; % AN CVs
    N.CV_Region_SEP = N.N_CV_AN              + 1 : N.N_CV_AN + N.N_CV_SEP; % SEP CVs
    N.CV_Region_CA  = N.N_CV_AN + N.N_CV_SEP + 1 : N.N_CV_tot; % CA CVs

%% Numerical Discretization
% Del-x of each region
    % ---- Anode ----
    AN.del_x = AN.L/(N.N_CV_AN);
    AN.x_vec = (AN.del_x/2):(AN.del_x):(AN.L-AN.del_x/2);
    AN.x_half_vec = 0 : AN.del_x : AN.L;
    
    % ---- Separator ----
    SEP.del_x = SEP.L/(N.N_CV_SEP);
    SEP.x_vec = (AN.L + SEP.del_x/2):(SEP.del_x):(AN.L+SEP.L-SEP.del_x/2);
    SEP.x_half_vec = AN.L : SEP.del_x : (AN.L+SEP.L);
    
    % ---- Cathode ----
    CA.del_x = CA.L/(N.N_CV_CA);
    CA.x_vec = ((AN.L+SEP.L) + CA.del_x/2):(CA.del_x):((AN.L+SEP.L+CA.L)-CA.del_x/2);
    CA.x_half_vec = (AN.L+SEP.L) : CA.del_x : (AN.L+SEP.L+CA.L);
    
    % Overall X-position Vector
    SIM.x_vec      = [AN.x_vec     , SEP.x_vec            , CA.x_vec            ];
    SIM.x_half_vec = [AN.x_half_vec, SEP.x_half_vec(2:end), CA.x_half_vec(2:end)];
    
% Del-r of each region
    % ---- Anode ----
    AN.del_r = AN.r_p/(N.N_R_AN);
    AN.r_vec = (AN.del_r/2):AN.del_r:AN.r_p-(AN.del_r/2);
    AN.r_half_vec = 0:AN.del_r:AN.r_p;
    
    % ---- Cathode ----
    CA.del_r = CA.r_p/(N.N_R_CA);
    CA.r_vec = (CA.del_r/2):CA.del_r:CA.r_p-(CA.del_r/2);
    CA.r_half_vec = 0:CA.del_r:CA.r_p;

%% Various Needed Calculations
% Initialize PROPS
    PROPS = zeros( N.N_prop , N.N_CV_tot );
    
    PROPS( P.sigma      , N.CV_Region_AN ) = AN.sigma * ones( 1 , N.N_CV_AN );
    PROPS( P.sigma      , N.CV_Region_CA ) = CA.sigma * ones( 1 , N.N_CV_CA );
    PROPS( P.kappa      , : )              = EL.kappa * ones( 1 , N.N_CV_tot );
    PROPS( P.R_SEI      , N.CV_Region_AN ) = AN.R_SEI * ones( 1 , N.N_CV_AN );
    PROPS( P.R_SEI      , N.CV_Region_CA ) = CA.R_SEI * ones( 1 , N.N_CV_CA );
    PROPS( P.k          , : )              = [AN.k  *ones(1,N.N_CV_AN), SEP.k  *ones(1,N.N_CV_SEP), CA.k  *ones(1,N.N_CV_CA)];
    PROPS( P.rho        , : )              = [AN.rho*ones(1,N.N_CV_AN), SEP.rho*ones(1,N.N_CV_SEP), CA.rho*ones(1,N.N_CV_CA)];
    PROPS( P.c_p        , : )              = [AN.c_p*ones(1,N.N_CV_AN), SEP.c_p*ones(1,N.N_CV_SEP), CA.c_p*ones(1,N.N_CV_CA)];
    PROPS( P.D_o_Li_ion , : )              = EL.D_o_Li_ion * ones( 1 , N.N_CV_tot );
    PROPS( P.activity   , : )              = EL.Activity   * ones( 1 , N.N_CV_tot );
    PROPS( P.tf_num     , : )              = EL.tf_num     * ones( 1 , N.N_CV_tot );
    PROPS( P.D_o:end    , N.CV_Region_AN ) = AN.D_o        * ones( N.N_R_AN , N.N_CV_AN );
    PROPS( P.D_o:end    , N.CV_Region_CA ) = CA.D_o        * ones( N.N_R_CA , N.N_CV_CA );
    
% Electrolyte volume fraction
    AN.eps_el  = 1 - AN.eps_ed - AN.eps_b;
    SEP.eps_el =    SEP.eps;               %     SEP.eps_el = 1 - SEP.eps; %%%%%%%%%
    CA.eps_el  = 1 - CA.eps_ed - CA.eps_b;

% Initialize Bruggeman
    if FLAG.Bruggeman
        % Initialize
        CONS.BRUG = ones(size(PROPS));
        
        % Calculate Brug for each region/phase
        if FLAG.BRUG_ED
            tau_fac_AN_ED   = AN.gamma_brug * AN.eps_ed  ^ (1 - AN.alpha_brug);
            tau_fac_CA_ED   = CA.gamma_brug * CA.eps_ed  ^ (1 - CA.alpha_brug);
            tau_fac_SEP_ED  = 1; %%%%%%% Placeholder
            
            tau_fac_ED = [tau_fac_AN_ED*ones(1,N.N_CV_AN) , tau_fac_SEP_ED*ones(1,N.N_CV_SEP) , tau_fac_CA_ED*ones(1,N.N_CV_CA)];
            eps_ED     = [AN.eps_ed    *ones(1,N.N_CV_AN) , 1             *ones(1,N.N_CV_SEP) , CA.eps_ed    *ones(1,N.N_CV_CA)];
            brug_ED    = eps_ED ./ tau_fac_ED;
        end
        if FLAG.BRUG_EL
            tau_fac_AN_EL  = EL.gamma_brug * AN.eps_el  ^ (1 - EL.alpha_brug);
            tau_fac_CA_EL  = EL.gamma_brug * CA.eps_el  ^ (1 - EL.alpha_brug);
            tau_fac_SEP_EL = EL.gamma_brug * SEP.eps_el ^ (1 - EL.alpha_brug);
            
            tau_fac_EL = [tau_fac_AN_EL*ones(1,N.N_CV_AN) , tau_fac_SEP_EL*ones(1,N.N_CV_SEP) , tau_fac_CA_EL*ones(1,N.N_CV_CA)];
            eps_EL     = [AN.eps_el    *ones(1,N.N_CV_AN) , SEP.eps_el    *ones(1,N.N_CV_SEP) , CA.eps_el    *ones(1,N.N_CV_CA)];
            brug_EL    = eps_EL ./ tau_fac_EL;
        end
        
        % Fill Brug Matrix based on FLAG
        if FLAG.BRUG_sigma
            CONS.BRUG(P.sigma , :)      = brug_ED;
        end
        if FLAG.BRUG_kappa
            CONS.BRUG(P.kappa , :)      = brug_EL;
        end
        if FLAG.BRUG_activity
            CONS.BRUG(P.activity , :)   = brug_EL;
        end
        if FLAG.BRUG_D_Liion
            CONS.BRUG(P.D_o_Li_ion , :) = brug_EL;
        end
        if FLAG.BRUG_tf_num
            CONS.BRUG(P.tf_num , :)     = brug_EL;
        end
        if FLAG.BRUG_D_o_AN
            %!!!!!!! Not implemented yet
        end
        if FLAG.BRUG_D_o_CA
            %!!!!!!! Not implemented yet
        end
        
        % Apply Brug to PROPS
        PROPS = PROPS .* CONS.BRUG;
    end
%% Geometry Calcs
% Control Volume - Volume
    % A_c  ~ cross-sectional area in the x-direction
    % Vol  ~ Total volume of the electrode
    % dVol ~ Volume of each control volume
    
    % ---- Anode ----
    if SIM.cell_geo == 'p'
        AN.A_c  = AN.del_z * AN.del_y;
    elseif SIM.cell_geo == 'c'
        AN.A_c  = pi * (AN.diam / 2)^2;
    end
    AN.Vol  = AN.A_c * AN.L;
    AN.dVol = AN.A_c * AN.del_x;
    
    % ---- Separator ----
    if SIM.cell_geo == 'p'
        SEP.A_c  = SEP.del_z * SEP.del_y;
    elseif SIM.cell_geo == 'c'
        SEP.A_c  = pi * (SEP.diam / 2)^2;
    end
    SEP.Vol  = SEP.A_c * SEP.L;
    SEP.dVol = SEP.A_c * SEP.del_x;
    
    % ---- Cathode ----
    if SIM.cell_geo == 'p'
        CA.A_c  = CA.del_z * CA.del_y;
    elseif SIM.cell_geo == 'c'
        CA.A_c  = pi * (CA.diam / 2)^2;
    end
    CA.Vol  = CA.A_c * CA.L;
    CA.dVol = CA.A_c * CA.del_x;
    
    SIM.CV_vec = [AN.dVol*ones(1,N.N_CV_AN),...
                 SEP.dVol*ones(1,N.N_CV_SEP),...
                  CA.dVol*ones(1,N.N_CV_CA)];

% Various Geometry Parameters
    % ---- Anode ----
    AN.V_sphere   = (4/3)*pi*AN.r_p^3;         % [m^3],        Volume of a single spherical particle
    AN.V_ed       = AN.eps_ed * AN.Vol;        % [m^3],        Volume of the electrode
    AN.V_el       = AN.eps_el * AN.Vol;        % [m^3],        Volume of the electrolyte
    AN.dVol_el    = AN.V_el / N.N_CV_AN;       % [m^3 per CV], Volume of electrolyte in each control volume
    AN.Np_tot     = AN.V_ed / AN.V_sphere;     % Number of total spherical particles
    AN.Np_CV      = AN.Np_tot / N.N_CV_AN;     % Number of particles per control volume
    AN.A_surf     = 4*pi*AN.r_p^2;             % [m^2],        Surface area of a particle
    AN.A_surf     = AN.A_surf * AN.A_geo;      % [m^2],        Account for ecentricity of spherical particle
    AN.A_surf_tot = AN.A_surf * AN.Np_tot;     % [m^2],        Total surface area om the electrode region
    AN.A_surf_CV  = AN.A_surf_tot / N.N_CV_AN; % [m^2 per CV], Surface area per control volume
    AN.A_s        = AN.A_surf_tot / AN.Vol;    % [m^2_surf / m^3_CV]  Particle surface area per control volume
    
    if FLAG.AN_LI_FOIL
        AN.A_surf = AN.A_c;
        AN.A_surf_CV = AN.A_surf;
    end

    % ---- Separator ----
    SEP.V_el      = SEP.eps_el * SEP.Vol;
    SEP.dVol_el   = SEP.V_el / N.N_CV_SEP;
    
    % ---- Cathode ----
    CA.V_sphere   = (4/3)*pi*CA.r_p^3;         % [m^3],        Volume of a single spherical particle
    CA.V_ed       = CA.eps_ed * CA.Vol;        % [m^3],        Volume of the electrode
    CA.V_el       = CA.eps_el * CA.Vol;        % [m^3],        Volume of the electrolyte
    CA.dVol_el    = CA.V_el / N.N_CV_CA;       % [m^3 per CV], Volume of electrolyte in each control volume
    CA.Np_tot     = CA.V_ed / CA.V_sphere;     % Number of total spherical particles
    CA.Np_CV      = CA.Np_tot / N.N_CV_CA;     % Number of particles per control volume
    CA.A_surf     = 4*pi*CA.r_p^2;             % [m^2],        Surface area of a particle
    CA.A_surf     = CA.A_surf * AN.A_geo;      % [m^2],        Account for ecentricity of spherical particle
    CA.A_surf_tot = CA.A_surf * CA.Np_tot;     % [m^2],        Total surface area om the electrode region
    CA.A_surf_CV  = CA.A_surf_tot / N.N_CV_CA; % [m^2 per CV], Surface area per control volume
    CA.A_s        = CA.A_surf_tot / CA.Vol;    % [m^2_surf / m^3_CV]  Particle surface area per control volume
    
    if FLAG.CA_LI_FOIL
        CA.A_surf = CA.A_c;
        CA.A_surf_CV = CA.A_surf;
    end
    
    % Matrix of all surface areas
    SIM.A_surf_CV_vec = [AN.A_surf_CV*ones(1,N.N_CV_AN),...
                                    zeros(1,N.N_CV_SEP),...
                         CA.A_surf_CV*ones(1,N.N_CV_CA)];
    
% Radial Volume
    % ---- Anode ----
    AN.dVol_r = zeros(N.N_R_AN,1);
    for i = 1:N.N_R_AN
        AN.dVol_r(i) = (4*pi/3)*(AN.r_half_vec(i+1)^3 - AN.r_half_vec(i)^3); % dV of each control volume
    end
    
    % ---- Cathode ----
    CA.dVol_r = zeros(N.N_R_CA,1);
    for i = 1:N.N_R_CA
        CA.dVol_r(i) = (4*pi/3)*(CA.r_half_vec(i+1)^3 - CA.r_half_vec(i)^3); % dV of each control volume
    end
    
%% Capactiy Calculations 
% Max Concentration
    if ~FLAG.preDefined_C_max
        if ~FLAG.AN_LI_FOIL
            AN.C_Li_max = AN.rho / AN.MW;
        end
        if ~FLAG.CA_LI_FOIL
            CA.C_Li_max = CA.rho / CA.MW;
        end
    end
    
% Capacity Calc
    % AN.Cap = (CONS.F/(AN.MW*3600))*AN.V_ed*AN.rho; % [Ahr], Method 1: More theoretical approach
    AN.Cap = AN.C_Li_max * CONS.F *(1/3600)*AN.V_ed; % [Ahr], Method 2: Based on a modeled parameter of C_max
    % AN.Cap_test3 = AN.specCap*AN.V_ed*AN.rho;      % [Ahr], Method 3: Based on something that is more measureable

    % CA.Cap = (CONS.F/(CA.MW*3600))*CA.V_ed*CA.rho; % [Ahr], Method 1: More theoretical approach
    CA.Cap = CA.C_Li_max * CONS.F *(1/3600)*CA.V_ed; % [Ahr], Method 2: Based on a modeled parameter of C_max
    % CA.Cap = CA.specCap*CA.V_ed*CA.rho;            % [Ahr], Method 3: Based on something that is more measureable
    % CA.Cap = CA.Cap * 0.6;                         % [Ahr], Modified Method 2, assumes only 60 percent of the capacity is being used


% Capacity Ratio
    z = (CA.C_Li_max*CA.V_ed)/(AN.C_Li_max*AN.V_ed);

%% Mole Fraction Calculation
    % x and y are the anode and cathode mole fractions respectively
    % Points that are being determined
    % * F: Mole fractions at formation
    % * A: Mole fraction when x = 0
    % * B: Mole fraction when x = 1
    % * C: Mole fraction at V_max
    % * D: Mole fraction at V_min
    options = optimoptions('lsqnonlin');
    options.Display = 'off';

    F_x = SIM.AnodeFormation_X;
    F_y = SIM.CathodeFormation_X;
    % y-intercept
    y_intcep = F_y + z * F_x;

    % y mole fraction at x limits
    A_x = 0;
    A_y = -z*A_x + y_intcep;
    B_x = 1;
    B_y = -z*B_x + y_intcep;

    % Mole fractions at V_max
    x0 = F_x;
    lb = 0;
    ub = 1;
    V_des = SIM.VoltageMax;
    C_x = lsqnonlin(@(x)XfromDesPotential(x,V_des,z,y_intcep,AN,CA),x0,lb,ub,options);
    C_y = YfromX(C_x,z,y_intcep);

    % Mole fractions at V_min
    x0 = F_x;
    lb = 0;
    ub = 1;
    V_des = SIM.VoltageMin;
    D_x = lsqnonlin(@(x)XfromDesPotential(x,V_des,z,y_intcep,AN,CA),x0,lb,ub,options);
    D_y = YfromX(D_x,z,y_intcep);

    % Limits on x and y
    SIM.x_min = max(A_x,D_x);
    SIM.x_max = min(B_x,C_x);
    SIM.y_min = YfromX(SIM.x_max,z,y_intcep);
    SIM.y_max = YfromX(SIM.x_min,z,y_intcep);

    % Initial Lithiation Fraction
    SIM.x_ini = (SIM.x_max-SIM.x_min)*SIM.SOC_start/100 + SIM.x_min;
    SIM.y_ini = YfromX(SIM.x_ini,z,y_intcep);

    % Determine electrode initial voltage potential
    voltage_ini_an = AN.EqPotentialHandle(SIM.x_ini);
    voltage_ini_ca = CA.EqPotentialHandle(SIM.y_ini);
    voltage_ini_cell = voltage_ini_ca - voltage_ini_an;

    % Calculate electrode initial active material concentration
    concen_ini_an = AN.C_Li_max * SIM.x_ini;
    concen_ini_ca = CA.C_Li_max * SIM.y_ini;
    
%% Load current
% ---- Polarization ----
if SIM.SimMode == 1
    SIM.A_c        = min( AN.A_c, CA.A_c );
    SIM.Cell_Cap   = min( AN.Cap, CA.Cap );
    SIM.I_user_amp = SIM.C_rate*SIM.Cell_Cap*SIM.ChargeOrDischarge; % [A], Load current
    SIM.i_user_amp = SIM.I_user_amp/SIM.A_c; % [A m^-2], Load current density (flux)
% ---- Harmonic Perturbation ----
elseif SIM.SimMode == 2
    SIM.A_c        = min( AN.A_c, CA.A_c );
    SIM.Cell_Cap   = min( AN.Cap, CA.Cap );
    SIM.I_user_amp = SIM.C_rate*SIM.Cell_Cap;
    SIM.i_user_amp = SIM.I_user_amp/SIM.A_c;
% ---- State Space EIS ----
elseif SIM.SimMode == 3
    SIM.A_c        = min( AN.A_c, CA.A_c );
    SIM.Cell_Cap   = min( AN.Cap, CA.Cap );
    SIM.I_user_amp = SIM.C_rate*SIM.Cell_Cap;
    SIM.i_user_amp = SIM.I_user_amp/SIM.A_c;
% ---- Known BC Profile Controller ----
elseif SIM.SimMode == 4 
    SIM.A_c        = min( AN.A_c, CA.A_c );
    SIM.Cell_Cap   = min( AN.Cap, CA.Cap );
% ---- MOO Controller ----
elseif SIM.SimMode == 5
    SIM.A_c        = min( AN.A_c, CA.A_c );
    SIM.Cell_Cap   = min( AN.Cap, CA.Cap );
    %!!!!!!!! May need to do something here for SV initialization with a
    %proper i_user value
% ---- Manual Current Profile ----
elseif SIM.SimMode == 7 
    SIM.A_c        = min( AN.A_c, CA.A_c );
    SIM.Cell_Cap   = min( AN.Cap, CA.Cap );
    SIM.I_user_amp = SIM.C_rate*SIM.Cell_Cap*SIM.ChargeOrDischarge;
    SIM.i_user_amp = SIM.I_user_amp/SIM.A_c;
    SIM.i_user_OG  = SIM.i_user_amp;
end

%% Determine Simulation Time Vector
% ---- Polarization ----
if SIM.SimMode == 1 
    if SIM.C_rate == 0
        t_final = 30; % [s], Final time
    else
        t_final = SIM.charge_frac*3600/SIM.C_rate + SIM.initial_offset + SIM.t_ramp; % [s], Final time
    end
%     SIM.tspan = [0,SIM.initial_offset, t_final];
    SIM.tspan = [0, t_final];
    
% ---- Harmonic Perturbation ----
elseif SIM.SimMode == 2
    SIM.f = SIM.freq / (2*pi);    % [cycles s^-1],  Time to complete a period or a full sinusoid cycle
    SIM.f_s = (2*SIM.f)*20;       % [samples s^-1], Based on Nyquist sampling theory (minimum is 2*f)
    t_sim = SIM.N_cycles / SIM.f; % [s],            Time required for N_cycles of sinusoidal curves
    N_samples = SIM.f_s*t_sim;    % [samples],      Number of samples total
    del_time = t_sim / N_samples; % [s sample^-1],  Time between each sample %%%%%%% is this also 1/f_s??????????
    
    sin_time_vector = 0:del_time:t_sim;                     % Vector of time spaces
    sin_time_vector = sin_time_vector + SIM.initial_offset; % Add the initial offset
    if SIM.initial_offset == 0
        SIM.tspan = sin_time_vector;                        % Overall simulation vector
    else
        SIM.tspan = [0, sin_time_vector];                   % Overall simulation vector
    end
    
% ---- State Space EIS ----
    % No t_span needed for this mode
    
% ---- Known BC Profile Controller ----
elseif SIM.SimMode == 4
    % Determine inside RunSimulation from MO_List
    
% ---- MOO Controller ----
elseif SIM.SimMode == 5
    % Determine inside RunSimulation from MO_List
    
% ---- Manual Current Profile ----
elseif SIM.SimMode == 7 % Manual Current Profile
    % Determine inside RunSimulation from MO_List

end

%% Determine the Output Matrix ( Only Mode 3)
if SIM.SimMode == 3
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
            idx_ed = index_offset + P.phi_ed;
            idx_el = index_offset + P.phi_el;
            SIM.OutputMatrix(j,idx_ed) =  1;
            SIM.OutputMatrix(j,idx_el) = -1;
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
end


%% Determine Initial State Vector
SV_IC = zeros(N.N_SV_tot,1);
% ---- Anode ----
for i = 1:N.N_CV_AN
    index_offset = (i-1)*N.N_SV_AN; 
    % Temp
    SV_IC(index_offset + P.T)       =  SIM.Temp_start;
    % phi_el
    SV_IC(index_offset + P.phi_el)  = -voltage_ini_an;
    % phi_ed
    SV_IC(index_offset + P.phi_ed)  =  0;
    % V_1
    SV_IC(index_offset + P.V_1)     = -voltage_ini_an;
    % V_2
    SV_IC(index_offset + P.V_2)     = -voltage_ini_an;
    % i_PS
    SV_IC(index_offset + P.i_PS)    =  0;
    % C_Li^+
    SV_IC(index_offset + P.C_Liion) =  EL.C; 
    % C_Li
    for j = 1:N.N_R_AN
        SV_IC(index_offset+P.C_Li+j-1)   = concen_ini_an;
    end
end

% ---- Separator ----
for i = 1:N.N_CV_SEP
    index_offset = (i-1)*N.N_SV_SEP + N.N_SV_AN_tot;
    % Temp
    SV_IC(index_offset + P.SEP.T)           =  SIM.Temp_start;
    % phi_el
    SV_IC(index_offset + P.SEP.phi_el)      = -voltage_ini_an;
    % C_Li^+
    SV_IC(index_offset + P.SEP.C_Liion)     =  EL.C;
end

% ---- Cathode ----
for i = 1:N.N_CV_CA
    index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
    % Temp
    SV_IC(index_offset + P.T)       =  SIM.Temp_start;
    % phi_el
    SV_IC(index_offset + P.phi_el)  = -voltage_ini_an;
    % phi_ed
    SV_IC(index_offset + P.phi_ed)  =  voltage_ini_cell;
    % V_1
    SV_IC(index_offset + P.V_1)     = -voltage_ini_an;
    % V_2
    SV_IC(index_offset + P.V_2)     = -voltage_ini_an;
    % i_PS
    SV_IC(index_offset + P.i_PS)    =  0; 
    % C_Li^+    
    SV_IC(index_offset + P.C_Liion) =  EL.C;
    % C_Li
    for j = 1:N.N_R_CA
        SV_IC(index_offset+N.N_SV_nR+j)   = concen_ini_ca;
    end
end

%% Save the Equilibrium Values of the Output (Only Mode 3)
if SIM.SimMode == 3
    SIM.OutputAtEquil = SIM.OutputMatrix*SV_IC;    
end

%% Solve for better Initial conditions for phi
SV = SV1Dto2D(SV_IC , N , P , FLAG);
if ( FLAG.CONSTANT_PROPS_FROM_HANDLES || FLAG.VARIABLE_PROPS_FROM_HANDLES)
    props = getProps( SV , AN , SEP, CA , EL , P , N , CONS , FLAG , PROPS);
    if FLAG.CONSTANT_PROPS_FROM_HANDLES
        % Update PROPS to use initial conditions
        PROPS = props;
    end
else
   props = PROPS; 
end

% Set i_user
    if SIM.SimMode == 3 % SS EIS
        i_user = SIM.i_user_amp; 
        SIM.i_user = i_user; % Needed for the SS impedance function
    elseif SIM.SimMode == 4 % KBCP
        SIM.current_MO_step = 1;
        MO = SIM.Controller_MO_File(SIM.current_MO_step).MO;
        if MO == 2
            i_user = SIM.Cell_Cap/20; %Start fsolve assuming the voltage should be close
        else
            if SIM.Controller_MO_File(SIM.current_MO_step).CorD == 'C'
                i_user = -SIM.Controller_MO_File(SIM.current_MO_step).C_rate * SIM.Cell_Cap / SIM.A_c;
            else
                i_user =  SIM.Controller_MO_File(SIM.current_MO_step).C_rate * SIM.Cell_Cap / SIM.A_c;
            end
        end
    elseif SIM.SimMode == 5 % MOO Controller
        i_user = 0;
    else
        i_user = i_user_calc(0,SIM); 
    end

phi_temp = SV(P.phi_el:P.i_PS, :);
phi_temp(end+1,: ) = zeros(1,N.N_CV_tot ); % Adding a row for the new constraint equation(i_dl)
phi_temp(end,N.CV_Region_SEP) =   NaN(1,N.N_CV_SEP); % No SV in the SEP region for the new constraint
phi_temp = reshape(phi_temp, [],1);

% Remove NaN
    phi_guess = zeros((N.N_ES_var)*N.N_CV_AN + N.N_CV_SEP + (N.N_ES_var)*N.N_CV_CA , 1);
    
    find_non_NaN = ~isnan(phi_temp);
    count = 0;
    for i = 1:length(find_non_NaN)
        if find_non_NaN(i)
            count = count + 1;
            phi_guess(count,1) = phi_temp(i);
        end
    end
    
    phi_guess(end+1) = i_user;
    
    SIM.phi_guess = phi_guess;
    
% Solve for better phi values
phi_soln = fsolve(@(phi) phiFsolveFun(phi,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,i_user,props) , phi_guess , SIM.fsolve_options);
i_user = phi_soln(end);
phi = phi1Dto2D( phi_soln(1:end-1) , N , P , FLAG);

%% Fix phi for Initial State Vector
% ---- Anode ----
for i = 1:N.N_CV_AN
    index_offset = (i-1)*N.N_SV_AN; 
    % phi_el
    SV_IC(index_offset + P.phi_el) =  phi(P.ES.phi_el,i); 
    % phi_ed
    SV_IC(index_offset + P.phi_ed) =  phi(P.ES.phi_ed,i);  
    % V_1
    SV_IC(index_offset + P.V_1)    =  phi(P.ES.V_1   ,i);  
    % V_2
    SV_IC(index_offset + P.V_2)    =  phi(P.ES.V_2   ,i);  
    % i_PS
    SV_IC(index_offset + P.i_PS)   =  phi(P.ES.i_PS  ,i); 
end

% ---- Separator ----
for i = 1:N.N_CV_SEP
    index_offset = (i-1)*N.N_SV_SEP + N.N_SV_AN_tot;
    % phi_el
    SV_IC(index_offset + P.SEP.phi_el) =  phi(P.ES.phi_el,i+N.N_CV_AN); 
end

% ---- Cathode ----
for i = 1:N.N_CV_CA
    index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
    % phi_el
    SV_IC(index_offset + P.phi_el) =  phi(P.ES.phi_el , i+N.N_CV_AN+N.N_CV_SEP); 
    % phi_ed
    SV_IC(index_offset + P.phi_ed) =  phi(P.ES.phi_ed , i+N.N_CV_AN+N.N_CV_SEP);  
    % V_1
    SV_IC(index_offset + P.V_1   ) =  phi(P.ES.V_1    , i+N.N_CV_AN+N.N_CV_SEP);  
    % V_2
    SV_IC(index_offset + P.V_2   ) =  phi(P.ES.V_2    , i+N.N_CV_AN+N.N_CV_SEP);  
    % i_PS
    SV_IC(index_offset + P.i_PS  ) =  phi(P.ES.i_PS   , i+N.N_CV_AN+N.N_CV_SEP); 
end

SIM.SV_IC = SV_IC;

%% Mass Matrix
M = zeros(N.N_SV_tot,N.N_SV_tot);
sim_cap = 0; % A dummy capacitance to help DAE solver by turning algebraics into ODEs. With solver working, this isn't needed anymore

% ---- Anode ----
for i = 1:N.N_CV_AN
    index_offset = (i-1)*N.N_SV_AN;
    % Temp
    M(index_offset+P.T       , index_offset+P.T )        =  1;
    % phi_el
    M(index_offset+P.phi_el  , index_offset+P.phi_el )   =  -AN.C_dl;
    M(index_offset+P.phi_el  , index_offset+P.phi_ed )   =   AN.C_dl;
    % phi_ed
    M(index_offset+P.phi_ed  , index_offset+P.phi_el )   =  -AN.C_dl; 
    M(index_offset+P.phi_ed  , index_offset+P.phi_ed )   =   AN.C_dl; 
    % V_1
    M(index_offset+P.V_1     , index_offset+P.V_1    )   =  sim_cap; 
    % V_2
    M(index_offset+P.V_2     , index_offset+P.V_2    )   =  sim_cap; 
    % i_PS
    M(index_offset+P.i_PS    , index_offset+P.i_PS   )   =  sim_cap;
    % C_Li^+
    M(index_offset+P.C_Liion , index_offset+P.C_Liion)  =  1; 
    % C_Li
    for j = 1:N.N_R_AN
        M(index_offset+P.C_Li+j-1 , index_offset+P.C_Li+j-1) = 1;
    end
end

% ---- Separator ----
for i = 1:N.N_CV_SEP
    index_offset = (i-1)*N.N_SV_SEP + N.N_SV_AN_tot;
    % Temp
    M(index_offset+P.SEP.T       , index_offset+P.SEP.T )        =  1;
    % phi_el
    M(index_offset+P.SEP.phi_el  , index_offset+P.SEP.phi_el )   =  sim_cap;
    % C_Li^+
    M(index_offset+P.SEP.C_Liion , index_offset+P.SEP.C_Liion )  =  1; 
end

% ---- Cathode ----
for i = 1:N.N_CV_CA
    index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
    % Temp
    M(index_offset+P.T       , index_offset+P.T )        =  1;
    % phi_el
    M(index_offset+P.phi_el  , index_offset+P.phi_el )   =  -CA.C_dl;
    M(index_offset+P.phi_el  , index_offset+P.phi_ed )   =   CA.C_dl;
    % phi_ed
    M(index_offset+P.phi_ed  , index_offset+P.phi_el )   =  -CA.C_dl; 
    M(index_offset+P.phi_ed  , index_offset+P.phi_ed )   =   CA.C_dl; 
    % V_1
    M(index_offset+P.V_1     , index_offset+P.V_1    )   =  sim_cap; 
    % V_2
    M(index_offset+P.V_2     , index_offset+P.V_2    )   =  sim_cap; 
    % i_PS
    M(index_offset+P.i_PS    , index_offset+P.i_PS   )   =  sim_cap;
    % C_Li^+
    M(index_offset+P.C_Liion , index_offset+P.C_Liion )  =  1; 
    % C_Li
    for j = 1:N.N_R_CA
        M(index_offset+P.C_Li+j-1 , index_offset+P.C_Li+j-1) = 1;
    end
end

% Fix phi_ed BC eqn @ AN/CC
i = 1;
offset = N.N_SV_AN*(i-1);
M(P.phi_ed + offset , :) = zeros(1,N.N_SV_tot);

SIM.M = M;

%% Set up Jacobian here too (sparse)
%%%%% Later


end
%% Some functions to help with calculations
function y_out = YfromX(x,z,y_intcep)
    y_out = -z*x + y_intcep;
end

function res = XfromDesPotential(x,V_des,z,y_intcep,AN,CA)
    % Solve for y
%     y = -z*x + y_intcep;
    y = YfromX(x,z,y_intcep);
    % Solve An potential
    Eeq_an = AN.EqPotentialHandle(x);
    % Solve Ca potential
    Eeq_ca = CA.EqPotentialHandle(y);
    % Calc residual
    Eeq_cell = Eeq_ca - Eeq_an;
    res = Eeq_cell - V_des;
end