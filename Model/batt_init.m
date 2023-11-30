%% Batt Initialization
% This function is used to initialize all of the user inputs and calculates
% the necessary parameters and pointers for the model

function [AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS] = batt_init(AN,CA,SEP,EL,SIM,N,FLAG)
%% Constants
    CONS.F = 96485338.3; % [C kmol^-1],      Faraday's Constant 
    CONS.R = 8314.472;   % [J kmol^-1 K^-1], Gas Constant


%% Control Volume Modification for Li Foil
% Adjust Number of CV
    % if FLAG.AN_LI_FOIL
    %     N.N_CV_AN = 1; % Single control volume in the x-dir
    %     FLAG.R_AN = 0; % 1 if radial gradients are considered
    %     N.N_R_AN  = 1; % Number of radial control volumes
    % end
    % if FLAG.CA_LI_FOIL
    %     N.N_CV_CA = 1; % Single control volume in the x-dir
    %     FLAG.R_CA = 0; % 1 if radial gradients are considered     
    %     N.N_R_CA  = 1; % Number of radial control volumes
    % end


%% Control Volume Modification for Distributed
    % If radial concentration gradients are not considered, then set the number of radial control volumes equal to 1
    if ~FLAG.R_AN
        N.N_R_AN  = 1;
    end
    if ~FLAG.R_CA
        N.N_R_CA  = 1;
    end 
    
    N.N_R_max = max([N.N_R_AN,N.N_R_CA]);


%% Pointers for Tracked Variables
    i = 1;
    P.T         = i; i = i + 1;
    P.del_phi   = i; i = i + 1; % This row will be del_phi for the AN and CA, then phi_el for the SEP for initialization. Inside fsolve and GovnEqns, phi_el and del_phi are restructured properly.
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

% Extra Pointer for phi_el
    P.phi_el     = P.C_Li_surf_max + 1;
    
% Pointers for Electrostatic 
    i = 1;
    P.ES.del_phi = i; i = i + 1;
    P.ES.phi_ed  = i; i = i + 1;
    P.ES.V_1     = i; i = i + 1;
    P.ES.V_2     = i; i = i + 1;
    P.ES.i_PS    = i; i = i + 1;
    P.ES.i_dl    = i; i = i + 1;
    
    N.N_ES_var   = 6;

% Pointers for State Space EIS
if SIM.SimMode == 3 || SIM.SimMode == 10 || SIM.SimMode == 9
    i = 1;
    P.SS.omega    = i; i = i + 1;
    P.SS.Z_mag    = i; i = i + 1;
    P.SS.Z_Re     = i; i = i + 1;
    P.SS.Z_Im     = i; i = i + 1;
    P.SS.Z_dB     = i; i = i + 1;
    P.SS.Z_ps_deg = i; i = i + 1;
end

% Number of nonRadial SV
    N.N_SV_nR  = P.C_Liion;    
    N.N_SV_SEP = P.SEP.C_Liion;
    
% Pointers for properties
    i = 1;
    P.sigma      = i; i = i + 1;
    P.kappa      = i; i = i + 1;
    P.R_SEI      = i; i = i + 1;
    P.lambda_eff = i; i = i + 1;
    P.rho_eff    = i; i = i + 1;
    P.c_p_eff    = i; i = i + 1;
    P.D_o_Li_ion = i; i = i + 1;
    P.activity   = i; i = i + 1;
    P.tf_num     = i; i = i + 1;
    P.D_o        = i; i = i + 1;
    
    N.N_prop = P.D_o + N.N_R_max - 1;

% Pointers for Output Matrix
    i = 1;
    P.OM.cell_volt  = i; i = i + 1;
    P.OM.delta_phi  = i; i = i + 1;
    P.OM.i_Far      = i; i = i + 1;
    P.OM.eta        = i; i = i + 1;
    P.OM.C_Liion    = i; i = i + 1;
    P.OM.C_Li       = i; i = i + 1;
    P.OM.delta_C_Li = i; i = i + 1;
    P.OM.T          = i; i = i + 1;
    P.OM.V_AN       = i; i = i + 1;
    P.OM.V_SEP      = i; i = i + 1;
    P.OM.V_CA       = i; i = i + 1;
    
    N.N_Out = length(fieldnames(P.OM));


%% Region Indexing
% Number of state variables in each CV
    N.N_SV_AN  = N.N_SV_nR + N.N_R_AN;
    N.N_SV_SEP = N.N_SV_SEP;
    N.N_SV_CA  = N.N_SV_nR + N.N_R_CA;

    N.N_SV_max = max([N.N_SV_AN , N.N_SV_SEP , N.N_SV_CA]);

    N.N_SV_AN_tot  = N.N_CV_AN  * N.N_SV_AN;
    N.N_SV_SEP_tot = N.N_CV_SEP * N.N_SV_SEP;
    N.N_SV_CA_tot  = N.N_CV_CA  * N.N_SV_CA;

    N.N_SV_tot = N.N_SV_AN_tot + N.N_SV_SEP_tot + N.N_SV_CA_tot;
    N.N_CV_tot = N.N_CV_AN + N.N_CV_SEP + N.N_CV_CA;

% Regional Indexing 
    N.CV_Region_AN  =                          1 : N.N_CV_AN;              % AN CVs
    N.CV_Region_SEP = N.N_CV_AN              + 1 : N.N_CV_AN + N.N_CV_SEP; % SEP CVs
    N.CV_Region_CA  = N.N_CV_AN + N.N_CV_SEP + 1 : N.N_CV_tot;             % CA CVs


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
        SIM.x_vec      = [AN.x_vec     , SEP.x_vec            , CA.x_vec            ]; % [m], Position of the control volume's center
        SIM.x_half_vec = [AN.x_half_vec, SEP.x_half_vec(2:end), CA.x_half_vec(2:end)]; % [m], Position of the control volume's edges
        
        diff = SIM.x_vec(2:end) - SIM.x_vec(1:end-1);
        SIM.diff_CV_x_vec  = [AN.del_x/2   , diff                 , CA.del_x/2          ]; % [m], Distance between the control volume's center on either side of the interface
        SIM.del_x_vec = [AN.del_x * ones(1,N.N_CV_AN) , SEP.del_x * ones(1,N.N_CV_SEP) , CA.del_x * ones(1,N.N_CV_CA)];
        SIM.del_x_vec_halved = 0.5*SIM.del_x_vec;

        SIM.interp_x_interface = nan(1,N.N_CV_tot+1);
        for i = 2:N.N_CV_tot
            SIM.interp_x_interface(i) = (SIM.x_half_vec(i) - SIM.x_vec(i-1)) / (SIM.x_vec(i) - SIM.x_vec(i-1));
        end
    
% Del-r of each region
    % ---- Anode ----
        AN.del_r = AN.r_p/(N.N_R_AN);
        AN.r_vec = (AN.del_r/2):AN.del_r:AN.r_p-(AN.del_r/2);
        AN.r_half_vec = 0:AN.del_r:AN.r_p;
        AN.r_half_vec_sq = AN.r_half_vec.^2;
        diff = AN.r_vec(2:end) - AN.r_vec(1:end-1);
        AN.del_CV_r_vec = [nan diff nan];
    
    % ---- Cathode ----
        CA.del_r = CA.r_p/(N.N_R_CA);
        CA.r_vec = (CA.del_r/2):CA.del_r:CA.r_p-(CA.del_r/2);
        CA.r_half_vec = 0:CA.del_r:CA.r_p;
        CA.r_half_vec_sq = CA.r_half_vec.^2;
        diff = CA.r_vec(2:end) - CA.r_vec(1:end-1);
        CA.del_CV_r_vec = [nan diff nan];

    AN_del_CV_r_vec_repmat  = repmat((AN.del_CV_r_vec.^(-1))' , 1 , N.N_CV_AN);
    CA_del_CV_r_vec_repmat  = repmat((CA.del_CV_r_vec.^(-1))' , 1 , N.N_CV_CA);
    SIM.del_CV_r_inv_mat = NaN(N.N_R_max+1 , N.N_CV_tot);
    SIM.del_CV_r_inv_mat(1:N.N_R_AN+1 , N.CV_Region_AN ) = AN_del_CV_r_vec_repmat;
    SIM.del_CV_r_inv_mat(1:N.N_R_CA+1 , N.CV_Region_CA ) = CA_del_CV_r_vec_repmat;


%% Various Needed Calculations
% Electrolyte volume fraction
    AN.eps_el   = 1 - AN.eps_ed - AN.eps_b; % Electrolyte Volume Fraction in the anode
    SEP.eps_el  =    SEP.eps;               % Electrolyte Volume Fraction in the separator
    CA.eps_el   = 1 - CA.eps_ed - CA.eps_b; % Electrolyte Volume Fraction in the cathode
    SEP.eps_sep = 1 - SEP.eps_el;           % Separator material Volume Fraction
    SIM.eps_el_vec = [AN.eps_el * ones(1,N.N_CV_AN) , SEP.eps_el * ones(1,N.N_CV_SEP) , CA.eps_el * ones(1,N.N_CV_CA)];
    SIM.MassFluxPreCalcResistor = SIM.eps_el_vec ./ SIM.del_x_vec_halved;

% Effective thermal properties
    [AN , SEP , CA] = getEffectiveThermalProps(AN , SEP , CA , EL);

% Initialize PROPS
    PROPS = nan( N.N_prop , N.N_CV_tot );
    
    PROPS( P.sigma      , N.CV_Region_AN ) = AN.sigma * ones( 1 , N.N_CV_AN );
    PROPS( P.sigma      , N.CV_Region_CA ) = CA.sigma * ones( 1 , N.N_CV_CA );
    PROPS( P.kappa      , : )              = EL.kappa * ones( 1 , N.N_CV_tot );
    PROPS( P.R_SEI      , N.CV_Region_AN ) = AN.R_SEI * ones( 1 , N.N_CV_AN );
    PROPS( P.R_SEI      , N.CV_Region_CA ) = CA.R_SEI * ones( 1 , N.N_CV_CA );
    PROPS( P.lambda_eff , : )              = [AN.lambda_eff * ones(1,N.N_CV_AN) , SEP.lambda_eff * ones(1,N.N_CV_SEP) , CA.lambda_eff * ones(1,N.N_CV_CA)];
    PROPS( P.rho_eff    , : )              = [AN.rho_eff    * ones(1,N.N_CV_AN) , SEP.rho_eff    * ones(1,N.N_CV_SEP) , CA.rho_eff    * ones(1,N.N_CV_CA)];
    PROPS( P.c_p_eff    , : )              = [AN.c_p_eff    * ones(1,N.N_CV_AN) , SEP.c_p_eff    * ones(1,N.N_CV_SEP) , CA.c_p_eff    * ones(1,N.N_CV_CA)];
    PROPS( P.D_o_Li_ion , : )              = EL.D_o_Li_ion * ones( 1 , N.N_CV_tot );
    PROPS( P.activity   , : )              = EL.Activity   * ones( 1 , N.N_CV_tot );
    PROPS( P.tf_num     , : )              = EL.tf_num     * ones( 1 , N.N_CV_tot );
    PROPS( P.D_o:(P.D_o+N.N_R_AN-1)    , N.CV_Region_AN ) = AN.D_o        * ones( N.N_R_AN , N.N_CV_AN );
    PROPS( P.D_o:(P.D_o+N.N_R_CA-1)    , N.CV_Region_CA ) = CA.D_o        * ones( N.N_R_CA , N.N_CV_CA );

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
            tau_fac_AN_EL  = EL.gamma_brug * AN.eps_el  ^ (1 - EL.alpha_brug_an);
            tau_fac_CA_EL  = EL.gamma_brug * CA.eps_el  ^ (1 - EL.alpha_brug_ca);
            tau_fac_SEP_EL = EL.gamma_brug * SEP.eps_el ^ (1 - EL.alpha_brug_sep);
            
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

% Initial Thermal Gradient
    if FLAG.InitialThermalGradient
        % L_tot        = AN.L+SEP.L+CA.L;
        L_tot        = SIM.x_vec(end) - SIM.x_vec(1);
        Delta_temp   = SIM.CA_Temp_init - SIM.AN_Temp_init;
        Slope        = Delta_temp / L_tot;
        temp_fnc     = @(x) Slope.*x + SIM.AN_Temp_init;
        pos_vec      = SIM.x_vec- SIM.x_vec(1); % Needed for proper BC
        % pos_vec(1)   = 0;
        % pos_vec(end) = AN.L+SEP.L+CA.L;
        Temp_vec     = temp_fnc(pos_vec);
    else
        Temp_vec = SIM.Temp_init * ones( 1 , N.N_CV_tot );
    end
    % % Test Plot
    % plot(SIM.x_vec,Temp_vec-273.15)

% Temperature BC Ramp
    if FLAG.RampThermalGradient % Assuming linear ramping
        % y-intercept
            AN.b_thermal = SIM.Temp_init;
            CA.b_thermal = SIM.Temp_init;

        % slope
            AN.m_thermal = (SIM.Temp_AN_BC - SIM.Temp_init)/SIM.RampThermalGradientTime;
            CA.m_thermal = (SIM.Temp_CA_BC - SIM.Temp_init)/SIM.RampThermalGradientTime;
    end


%% Geometry Calcs
% A_c       ~ cross-sectional area in the x-direction
% A_outside ~ Control Volume surface area
% Vol       ~ Total volume of the electrode
% dVol      ~ Volume of each control volume

% ---- Anode ----
    if SIM.cell_geo == 'p'
        AN.A_c       = AN.del_z * AN.del_y;
        AN.A_outside = (2*AN.del_z + 2*AN.del_y)*AN.del_x;
    elseif SIM.cell_geo == 'c'
        AN.A_c       = pi * (AN.diam / 2)^2;
        AN.A_outside = pi *  AN.diam * AN.del_x;
    end
    AN.Vol  = AN.A_c * AN.L;
    AN.dVol = AN.A_c * AN.del_x;

% ---- Separator ----
    if SIM.cell_geo == 'p'
        SEP.A_c  = SEP.del_z * SEP.del_y;
        SEP.A_outside = (2*SEP.del_z + 2*SEP.del_y)*SEP.del_x;
    elseif SIM.cell_geo == 'c'
        SEP.A_c  = pi * (SEP.diam / 2)^2;
        SEP.A_outside = pi *  SEP.diam * SEP.del_x;
    end
    SEP.Vol  = SEP.A_c * SEP.L;
    SEP.dVol = SEP.A_c * SEP.del_x;

% ---- Cathode ----
    if SIM.cell_geo == 'p'
        CA.A_c  = CA.del_z * CA.del_y;
        CA.A_outside = (2*CA.del_z + 2*CA.del_y)*CA.del_x;
    elseif SIM.cell_geo == 'c'
        CA.A_c  = pi * (CA.diam / 2)^2;
        CA.A_outside = pi *  CA.diam * CA.del_x;
    end
    CA.Vol  = CA.A_c * CA.L;
    CA.dVol = CA.A_c * CA.del_x;

% Vectorize
    SIM.CV_vec = [AN.dVol*ones(1,N.N_CV_AN),...
                 SEP.dVol*ones(1,N.N_CV_SEP),...
                  CA.dVol*ones(1,N.N_CV_CA)];

    SIM.A_outside_vec = [AN.A_outside*ones(1,N.N_CV_AN),...
                        SEP.A_outside*ones(1,N.N_CV_SEP),...
                         CA.A_outside*ones(1,N.N_CV_CA)];

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
        AN.A_surf_tot = AN.A_surf * AN.Np_tot;     % [m^2],        Total surface area in the electrode region
        AN.A_surf_CV  = AN.A_surf_tot / N.N_CV_AN; % [m^2 per CV], Surface area per control volume
        AN.A_s        = AN.A_surf_tot / AN.Vol;    % [m^2_surf / m^3_CV]  Particle surface area per control volume
        
        % if FLAG.AN_LI_FOIL
        %     AN.A_surf    = AN.A_c;
        %     AN.A_surf_CV = AN.A_surf;
        % end

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
    
        % if FLAG.CA_LI_FOIL
        %     CA.A_surf    = CA.A_c;
        %     CA.A_surf_CV = CA.A_surf;
        % end
    
    % Vectorize
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
    

%% Capacity Calculations 
% Max Concentration
    if ~FLAG.preDefined_C_max
    end
    
% Capacity Calc
    AN.Cap = AN.C_Li_max * CONS.F *(1/3600)*AN.V_ed; % [Ahr], Method 2: Based on a modeled parameter of C_max
    CA.Cap = CA.C_Li_max * CONS.F *(1/3600)*CA.V_ed; % [Ahr], Method 2: Based on a modeled parameter of C_max

if isempty(SIM.VoltageMax) 
    %% Mole Fraction Calculation
    % Determine V_min and V_max
        voltage_min_an = AN.EqPotentialHandle(SIM.AnodeStoich_SOC0);
        voltage_min_ca = CA.EqPotentialHandle(SIM.CathodeStoich_SOC0);
        SIM.VoltageMin = voltage_min_ca - voltage_min_an; % Voltage at SOC0
        V_min = voltage_min_ca - voltage_min_an;

        voltage_max_an = AN.EqPotentialHandle(SIM.AnodeStoich_SOC100);
        voltage_max_ca = CA.EqPotentialHandle(SIM.CathodeStoich_SOC100);
        SIM.VoltageMax = voltage_max_ca - voltage_max_an; % Voltage at SOC100
        V_max = voltage_max_ca - voltage_max_an;

    % Set stoich limits
        SIM.x_min = SIM.AnodeStoich_SOC0;
        SIM.x_max = SIM.AnodeStoich_SOC100;
        SIM.y_min = SIM.CathodeStoich_SOC100;
        SIM.y_max = SIM.CathodeStoich_SOC0;
    
    % Initial Lithiation Fraction
        SIM.x_ini = (SIM.x_max-SIM.x_min)*SIM.SOC_start/100 + SIM.x_min;
        SIM.y_ini = (SIM.y_min-SIM.y_max)*SIM.SOC_start/100 + SIM.y_max;
    
    % Determine electrode initial voltage potential
        voltage_ini_an   = AN.EqPotentialHandle(SIM.x_ini);
        voltage_ini_ca   = CA.EqPotentialHandle(SIM.y_ini);
        voltage_ini_cell = voltage_ini_ca - voltage_ini_an;
    
    % Calculate electrode initial active material concentration
        concen_ini_an = AN.C_Li_max * SIM.x_ini;
        concen_ini_ca = CA.C_Li_max * SIM.y_ini;
        
    % % Test OCP
    %     SOC_vec = 0:0.1:100;
    %     xFromSOC = @(SOC) (SIM.x_max-SIM.x_min)*SOC/100 + SIM.x_min;
    %     yFromSOC = @(SOC) (SIM.y_min-SIM.y_max)*SOC/100 + SIM.y_max;
    %     x_vec = xFromSOC(SOC_vec);
    %     y_vec = yFromSOC(SOC_vec);
    %     v_an_vec = AN.EqPotentialHandle(x_vec);
    %     v_ca_vec = CA.EqPotentialHandle(y_vec);
    %     v_vec = v_ca_vec - v_an_vec;
    %     figure
    %     plot(SOC_vec , v_vec , '-k' , 'LineWidth', 2)

else
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
        options         = optimoptions('lsqnonlin');
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
        x0    = F_x;
        lb    = 0;
        ub    = 1;
        V_des = SIM.VoltageMax;
        C_x   = lsqnonlin(@(x)XfromDesPotential(x,V_des,z,y_intcep,AN,CA),x0,lb,ub,options);
        C_y   = YfromX(C_x,z,y_intcep);
    
    % Mole fractions at V_min
        x0    = F_x;
        lb    = 0;
        ub    = 1;
        V_des = SIM.VoltageMin;
        D_x   = lsqnonlin(@(x)XfromDesPotential(x,V_des,z,y_intcep,AN,CA),x0,lb,ub,options);
        D_y   = YfromX(D_x,z,y_intcep);
    
    % Limits on x and y
        SIM.x_min = max(A_x,D_x);
        SIM.x_max = min(B_x,C_x);
        SIM.y_min = YfromX(SIM.x_max,z,y_intcep);
        SIM.y_max = YfromX(SIM.x_min,z,y_intcep);
    
    % Initial Lithiation Fraction
        SIM.x_ini = (SIM.x_max-SIM.x_min)*SIM.SOC_start/100 + SIM.x_min;
        SIM.y_ini = YfromX(SIM.x_ini,z,y_intcep);
    
    % Determine electrode initial voltage potential
        voltage_ini_an   = AN.EqPotentialHandle(SIM.x_ini);
        voltage_ini_ca   = CA.EqPotentialHandle(SIM.y_ini);
        voltage_ini_cell = voltage_ini_ca - voltage_ini_an;
    
    % Calculate electrode initial active material concentration
        concen_ini_an = AN.C_Li_max * SIM.x_ini;
        concen_ini_ca = CA.C_Li_max * SIM.y_ini;
end

 
%% Determine Minimum Parameters
    SIM.A_c      = min( AN.A_c, CA.A_c );
    SIM.Cell_Cap = min( AN.Cap, CA.Cap );


%% User Time Vector and Current
    [SIM] = UserCurrentProfile(SIM,FLAG);


%% Determine the Output Matrix 
%!!!!!!!!!!!!!! Add isfield to handle logic of turning off pointers

    N.N_In  = 1;
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
        % Delta V_AN  (phi_el,ANSEP - 0)
        % Delta V_SEP (phi_el,CASEP - phi_el,ANSEP)
        % Delta V_CA  (phi_ed,CA    - phi_el,CASEP)
    SIM.OutputMatrix = zeros(N.N_Out , N.N_SV_tot);
        % Cell Voltage
            idx_phi_ed_AN = P.phi_ed;

            i = N.N_CV_CA;
            index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
            idx_phi_ed_CA = index_offset + P.phi_ed;

            SIM.OutputMatrix(P.OM.cell_volt,idx_phi_ed_AN) = -1;
            SIM.OutputMatrix(P.OM.cell_volt,idx_phi_ed_CA) =  1;
            N.IDX_CellVoltage = [idx_phi_ed_AN , idx_phi_ed_CA]; % This is used in batt_events

        % @AN/SEP
            i = N.N_CV_AN;
            index_offset = (i-1)*N.N_SV_AN;

        % Delta Phi      @AN/SEP
            idx = index_offset + P.del_phi;
            SIM.OutputMatrix(P.OM.delta_phi,idx)  =  1;

        % i_Far          @AN/SEP
            idx = index_offset + P.i_PS;
            SIM.OutputMatrix(P.OM.i_Far,idx)      =  1;

        % Eta           @AN/SEP
            idx = index_offset + P.V_2;
            SIM.OutputMatrix(P.OM.eta,idx)        =  1;
            idx = index_offset + P.V_1;
            SIM.OutputMatrix(P.OM.eta,idx)        = -1;

        % C_Liion        @AN/SEP
            idx = index_offset + P.C_Liion;
            SIM.OutputMatrix(P.OM.C_Liion,idx)    =  1;

        % C_Li (Surface) @AN/SEP
            idx = index_offset + P.C_Li_surf_AN;
            SIM.OutputMatrix(P.OM.C_Li,idx)       =  1;%/AN.C_Li_max;

        % Delta C_Li     @AN/SEP (Over entire radius of particle)
            idx = index_offset + P.C_Li_surf_AN;     % Surface Node
            SIM.OutputMatrix(P.OM.delta_C_Li,idx) =  1;
            idx = index_offset + P.C_Li;         % Most Interior Node
            SIM.OutputMatrix(P.OM.delta_C_Li,idx) = -1;

        % Temperature    @AN/SEP
            idx = index_offset + P.T;
            SIM.OutputMatrix(P.OM.T,idx)          =  1;

        % Delta V_AN (phi_el,ANSEP - 0)
            i = N.N_CV_AN;
            index_offset = (i-1)*N.N_SV_AN;
            idx = index_offset + P.phi_ed;
            SIM.OutputMatrix(P.OM.V_AN , idx)     =  1;
            
            i = N.N_CV_AN;
            index_offset = (i-1)*N.N_SV_AN;
            idx = index_offset + P.del_phi;
            SIM.OutputMatrix(P.OM.V_AN , idx)     = -1;
            
        % Delta V_SEP (phi_el,CASEP - phi_el,ANSEP)
            i = N.N_CV_SEP;
            index_offset = (i-1)*N.N_SV_SEP + N.N_SV_AN_tot;
            idx = index_offset + P.SEP.phi_el;
            SIM.OutputMatrix(P.OM.V_SEP , idx)    =  1;

            i = 1;
            index_offset = (i-1)*N.N_SV_SEP + N.N_SV_AN_tot;
            idx = index_offset + P.SEP.phi_el;
            SIM.OutputMatrix(P.OM.V_SEP , idx)    = -1;
            
        % Delta V_CA (phi_ed,CA - phi_el,CASEP)
            i = N.N_CV_CA;
            index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
            idx = index_offset + P.phi_ed;
            SIM.OutputMatrix(P.OM.V_CA , idx)     =  1;

            i = 1;
            index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
            idx = index_offset + P.phi_ed;
            SIM.OutputMatrix(P.OM.V_CA , idx)     =  -(1);  

            i = 1;
            index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
            idx = index_offset + P.del_phi;
            SIM.OutputMatrix(P.OM.V_CA , idx)     = -(-1);        


%% Determine Initial State Vector
    SV_IC = zeros(N.N_SV_tot,1);

% ---- Anode ----
    for i = 1:N.N_CV_AN
        index_offset = (i-1)*N.N_SV_AN; 
        CV_offset    = i;
        % Temp
            SV_IC(index_offset + P.T)       =  Temp_vec(CV_offset);
        % delta phi
            SV_IC(index_offset + P.del_phi) =  voltage_ini_an;
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
                SV_IC(index_offset+P.C_Li+j-1) = concen_ini_an;
            end
    end
    
% ---- Separator ----
    for i = 1:N.N_CV_SEP
        index_offset = (i-1)*N.N_SV_SEP + N.N_SV_AN_tot;
        CV_offset    = i + N.N_CV_AN;
        % Temp
            SV_IC(index_offset + P.SEP.T)           =  Temp_vec(CV_offset);
        % phi_el
            SV_IC(index_offset + P.SEP.phi_el)      = -voltage_ini_an;
        % C_Li^+
            SV_IC(index_offset + P.SEP.C_Liion)     =  EL.C;
    end
    
% ---- Cathode ----
    for i = 1:N.N_CV_CA
        index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
        CV_offset    = i + N.N_CV_AN + N.N_CV_SEP;
        % Temp
            SV_IC(index_offset + P.T)       =  Temp_vec(CV_offset);
        % delta phi
            SV_IC(index_offset + P.del_phi) =  voltage_ini_ca;
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
                SV_IC(index_offset+N.N_SV_nR+j) = concen_ini_ca;
            end
    end


%% Save the Equilibrium Values of the Output
    SIM.OutputAtEquil = SIM.OutputMatrix*SV_IC;
    SIM.SV_IC_Static  = SV_IC;


%% Property Inverse Calcs
    CONS.F_inv = CONS.F^-1;
    CONS.R_inv = CONS.R^-1;

    AN.R_SEI_inv = AN.R_SEI^-1;
    CA.R_SEI_inv = CA.R_SEI^-1;

    AN.C_Li_max_inv = AN.C_Li_max^-1;
    CA.C_Li_max_inv = CA.C_Li_max^-1;

     AN.del_x_inv =  AN.del_x^-1;
    SEP.del_x_inv = SEP.del_x^-1;
     CA.del_x_inv =  CA.del_x^-1;

    AN.del_r_inv = AN.del_r^-1;
    CA.del_r_inv = CA.del_r^-1;

    AN.A_surf_CV_inv = AN.A_surf_CV^-1;
    CA.A_surf_CV_inv = CA.A_surf_CV^-1;

    for j = 1:N.N_R_AN
        AN.r_half_vec_diffcubed(j,1) = AN.r_half_vec(j+1)^3-AN.r_half_vec(j)^3;
    end
    AN.r_half_vec_diffcubed_inv = AN.r_half_vec_diffcubed.^-1;

    for j = 1:N.N_R_CA
        CA.r_half_vec_diffcubed(j,1) = CA.r_half_vec(j+1)^3-CA.r_half_vec(j)^3;
    end
    CA.r_half_vec_diffcubed_inv = CA.r_half_vec_diffcubed.^-1;

    EL.C_inv = EL.C^-1;

    SIM.diff_CV_x_vec_inv    = SIM.diff_CV_x_vec.^(-1);
    SIM.del_x_vec_halved_inv = SIM.del_x_vec_halved.^(-1);


%% Solve for better Initial conditions for phi
    SV = SV1Dto2D(SV_IC , N.N_SV_max, N.N_CV_tot, N.N_SV_AN_tot, N.N_SV_SEP_tot, N.N_SV_AN, N.N_SV_SEP, N.N_SV_CA, N.N_CV_AN, N.N_CV_SEP, N.N_CV_CA, N.CV_Region_AN, N.CV_Region_SEP, N.CV_Region_CA, P.T, P.del_phi, P.C_Liion, P.SEP.T, P.SEP.phi_el, P.SEP.C_Liion);
    if ( FLAG.CONSTANT_PROPS_FROM_HANDLES || FLAG.VARIABLE_PROPS_FROM_HANDLES)
        [T, ~, ~, ~, ~, ~, ~, Ce, ~, ~, ~, X_AN, X_CA, ~, ~, ~, ~] = extractSV(SV,P.T, P.del_phi, P.phi_ed, P.phi_el, P.V_1, P.V_2, P.i_PS, P.C_Liion, P.C_Li, P.C_Li_surf_AN, P.C_Li_surf_CA, N.CV_Region_AN, N.CV_Region_CA, N.N_R_max, AN.C_Li_max_inv, CA.C_Li_max_inv, EL.C_inv, CONS.R);
        props = getProps( Ce, T, X_AN,  X_CA, FLAG.VARIABLE_kappa, ...
                        FLAG.VARIABLE_D_Liion, FLAG.VARIABLE_activity, FLAG.VARIABLE_tf_num, FLAG.VARIABLE_D_o_AN, FLAG.VARIABLE_D_o_CA, FLAG.Bruggeman,...
                        P.kappa, P.D_o_Li_ion, P.activity, P.tf_num, P.D_o, ...
                        N.N_R_AN, N.N_R_CA, N.CV_Region_AN, N.CV_Region_CA, ...
                        CONS.BRUG, ...
                        AN.D_oHandle, CA.D_oHandle, ...
                        PROPS, EL.kappaHandle, EL.D_o_Li_ionHandle, EL.ActivityHandle, EL.tf_numHandle); 
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
        if MO == 2 % CV
%             i_user = SIM.Cell_Cap/20; %Start fsolve assuming the voltage should be close
            i_user = SIM.Cell_Cap; %Start fsolve assuming the voltage should be close
        else % CC or Relax
            if SIM.Controller_MO_File(SIM.current_MO_step).CorD == 'C'
                i_user = -SIM.Controller_MO_File(SIM.current_MO_step).C_rate * SIM.Cell_Cap / SIM.A_c;
            else
                i_user =  SIM.Controller_MO_File(SIM.current_MO_step).C_rate * SIM.Cell_Cap / SIM.A_c;
            end
        end
    elseif SIM.SimMode == 5 % MOO Controller
        i_user = 0;
    elseif SIM.SimMode == 8 % ---- PRBS ----
        if FLAG.AddInputNoise
            i_user = SIM.profile_current_Noisy2D(1,1);
        else
            i_user = 0;
        end
    else
        i_user_in = nan;
        % i_user = i_user_calc(0,SIM,FLAG,i_user_in);
        SIM.Amp = 0;
        i_user = i_user_calc(0, i_user_in, SIM.SimMode, SIM.profile_time, SIM.profile_current, SIM.Amp, SIM, FLAG);
    end

% Create phi vector    
    phi_temp                      = SV(P.del_phi:P.i_PS, :);
    phi_temp(end+1,: )            = zeros(1,N.N_CV_tot );    % Adding a row for the new constraint equation(i_dl)
    phi_temp(end,N.CV_Region_SEP) = NaN(1,N.N_CV_SEP);       % No SV in the SEP region for the new constraint
    phi_temp                      = reshape(phi_temp, [],1);

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
    i_user   = phi_soln(end);
    phi      = phi1Dto2D( phi_soln(1:end-1) , N , P , FLAG);


%% Fix phi for Initial State Vector
% ---- Anode ----
    for i = 1:N.N_CV_AN
        index_offset = (i-1)*N.N_SV_AN; 
        % phi_el
            SV_IC(index_offset + P.del_phi) =  phi(P.ES.del_phi,i); 
        % phi_ed
            SV_IC(index_offset + P.phi_ed)  =  phi(P.ES.phi_ed,i);  
        % V_1
            SV_IC(index_offset + P.V_1)     =  phi(P.ES.V_1   ,i);  
        % V_2
            SV_IC(index_offset + P.V_2)     =  phi(P.ES.V_2   ,i);  
        % i_PS
            SV_IC(index_offset + P.i_PS)    =  phi(P.ES.i_PS  ,i); 
    end

% ---- Separator ----
    for i = 1:N.N_CV_SEP
        index_offset = (i-1)*N.N_SV_SEP + N.N_SV_AN_tot;
        % phi_el
            SV_IC(index_offset + P.SEP.phi_el) =  phi(P.ES.del_phi,i+N.N_CV_AN); 
    end

% ---- Cathode ----
    for i = 1:N.N_CV_CA
        index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
        % phi_el
            SV_IC(index_offset + P.del_phi) =  phi(P.ES.del_phi , i+N.N_CV_AN+N.N_CV_SEP); 
        % phi_ed
            SV_IC(index_offset + P.phi_ed)  =  phi(P.ES.phi_ed , i+N.N_CV_AN+N.N_CV_SEP);  
        % V_1
            SV_IC(index_offset + P.V_1   )  =  phi(P.ES.V_1    , i+N.N_CV_AN+N.N_CV_SEP);  
        % V_2
            SV_IC(index_offset + P.V_2   )  =  phi(P.ES.V_2    , i+N.N_CV_AN+N.N_CV_SEP);  
        % i_PS
            SV_IC(index_offset + P.i_PS  )  =  phi(P.ES.i_PS   , i+N.N_CV_AN+N.N_CV_SEP); 
    end

    SIM.SV_IC = SV_IC;


%% SV Index 1D to 2D
    SV = SV1Dto2D(SV_IC , N.N_SV_max, N.N_CV_tot, N.N_SV_AN_tot, N.N_SV_SEP_tot, N.N_SV_AN, N.N_SV_SEP, N.N_SV_CA, N.N_CV_AN, N.N_CV_SEP, N.N_CV_CA, N.CV_Region_AN, N.CV_Region_SEP, N.CV_Region_CA, P.T, P.del_phi, P.C_Liion, P.SEP.T, P.SEP.phi_el, P.SEP.C_Liion);
    N.IDX_1Dto2D = find(~isnan(SV));
    SIM.SV_nan   = nan(N.N_SV_max, N.N_CV_tot);


%% Mass Matrix
    M = zeros(N.N_SV_tot,N.N_SV_tot);
    sim_cap = 0; % A dummy capacitance to help DAE solver by turning algebraics into ODEs. With solver working, this isn't needed anymore
    
% ---- Anode ----
    for i = 1:N.N_CV_AN
        index_offset = (i-1)*N.N_SV_AN;
        % Temp
            M(index_offset+P.T       , index_offset+P.T )        =  AN.rho_eff * AN.c_p_eff * AN.dVol; 
        % del_phi
            M(index_offset+P.del_phi , index_offset+P.del_phi )  =  AN.C_dl;
        % phi_ed 
            M(index_offset+P.phi_ed  , index_offset+P.phi_ed )   =  sim_cap; 
        % V_1
            M(index_offset+P.V_1     , index_offset+P.V_1    )   =  sim_cap; 
        % V_2
            M(index_offset+P.V_2     , index_offset+P.V_2    )   =  sim_cap; 
        % i_PS
            M(index_offset+P.i_PS    , index_offset+P.i_PS   )   =  sim_cap;
        % C_Li^+
            M(index_offset+P.C_Liion , index_offset+P.C_Liion)   =  AN.eps_el; 
        % C_Li
            for j = 1:N.N_R_AN
                M(index_offset+P.C_Li+j-1 , index_offset+P.C_Li+j-1) = 1;
            end
    end

% ---- Separator ----
    for i = 1:N.N_CV_SEP
        index_offset = (i-1)*N.N_SV_SEP + N.N_SV_AN_tot;
        % Temp
            M(index_offset+P.SEP.T       , index_offset+P.SEP.T )        =  SEP.rho_eff * SEP.c_p_eff * SEP.dVol;
        % phi_el
            M(index_offset+P.SEP.phi_el  , index_offset+P.SEP.phi_el )   =  sim_cap;
        % C_Li^+
            M(index_offset+P.SEP.C_Liion , index_offset+P.SEP.C_Liion )  =  SEP.eps_el; 
    end

% ---- Cathode ----
    for i = 1:N.N_CV_CA
        index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
        % Temp
            M(index_offset+P.T       , index_offset+P.T )        =  CA.rho_eff * CA.c_p_eff * CA.dVol;
        % del_phi
            M(index_offset+P.del_phi , index_offset+P.del_phi )  =  CA.C_dl;
        % phi_ed
            M(index_offset+P.phi_ed  , index_offset+P.phi_ed )   =  sim_cap; 
        % V_1
            M(index_offset+P.V_1     , index_offset+P.V_1    )   =  sim_cap; 
        % V_2
            M(index_offset+P.V_2     , index_offset+P.V_2    )   =  sim_cap; 
        % i_PS
            M(index_offset+P.i_PS    , index_offset+P.i_PS   )   =  sim_cap;
        % C_Li^+
            M(index_offset+P.C_Liion , index_offset+P.C_Liion )  =  CA.eps_el; 
        % C_Li
            for j = 1:N.N_R_CA
                M(index_offset+P.C_Li+j-1 , index_offset+P.C_Li+j-1) = 1;
            end
    end

% Fix phi_ed BC eqn @ AN/CC
    i = 1;
    offset = N.N_SV_AN*(i-1);
    M(P.phi_ed + offset , :) = zeros(1,N.N_SV_tot);

% Fix Known Temperature BC
    if FLAG.COE
        if FLAG.T_BC_AN == 1 % 1) Known Temperature
            i = 1;
            index_offset = (i-1)*N.N_SV_AN;
            M(index_offset+P.T , index_offset+P.T ) =  sim_cap; 
        end
        if FLAG.T_BC_CA == 1 % 1) Known Temperature
            i = N.N_CV_CA;
            index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
            M(index_offset+P.T , index_offset+P.T ) =  sim_cap;
        end
    end

    SIM.M = M;


%% Make Mass for just diffEq
% Make indices vector
    idx_diff = [];
    idx_algb = [];
    
% ---- Anode ----
    for i = 1:N.N_CV_AN
        index_offset = (i-1)*N.N_SV_AN;
        % Temp
            if FLAG.COE && FLAG.T_BC_AN == 1 && i == 1
                idx_algb(end+1) = index_offset + P.T;
            else
                idx_diff(end+1) = index_offset + P.T;
            end
        % del_phi
            idx_diff(end+1) = index_offset + P.del_phi;
        % phi_ed
            idx_algb(end+1) = index_offset + P.phi_ed;
        % V_1
            idx_algb(end+1) = index_offset + P.V_1;
        % V_2
            idx_algb(end+1) = index_offset + P.V_2;
        % i_PS
            idx_algb(end+1) = index_offset + P.i_PS;
        % C_Li^+
            idx_diff(end+1) = index_offset + P.C_Liion;
        % C_Li
            for j = 1:N.N_R_AN
                idx_diff(end+1) = index_offset + P.C_Li + j-1;
            end
    end

% ---- Separator ----
    for i = 1:N.N_CV_SEP
        index_offset = (i-1)*N.N_SV_SEP + N.N_SV_AN_tot;
        % Temp
            idx_diff(end+1) = index_offset + P.SEP.T;
        % phi_el
            idx_algb(end+1) = index_offset + P.SEP.phi_el;
        % C_Li^+
            idx_diff(end+1) = index_offset + P.SEP.C_Liion;
    end

% ---- Cathode ----
    for i = 1:N.N_CV_CA
        index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
        % Temp
            if FLAG.COE && FLAG.T_BC_CA == 1 && i == N.N_CV_CA
                idx_algb(end+1) = index_offset + P.T;
            else
                idx_diff(end+1) = index_offset + P.T;
            end
        % del_phi
            idx_diff(end+1) = index_offset + P.del_phi;
        % phi_ed
            idx_algb(end+1) = index_offset + P.phi_ed;
        % V_1
            idx_algb(end+1) = index_offset + P.V_1;
        % V_2
            idx_algb(end+1) = index_offset + P.V_2;
        % i_PS
            idx_algb(end+1) = index_offset + P.i_PS;
        % C_Li^+
            idx_diff(end+1) = index_offset + P.C_Liion;
        % C_Li
            for j = 1:N.N_R_CA
                idx_diff(end+1) = index_offset + P.C_Li + j-1;
            end
    end

    SIM.diff_idx = idx_diff;
    SIM.algb_idx = idx_algb;

    N.N_diff = length(SIM.diff_idx);
    N.N_algb = length(SIM.algb_idx);


%% Set up Jacobian here too (sparse)
    if N.N_SV_tot == 855
        data = load('JPattern_sparse855MatlabVaryProp.mat','JPattern_sparse');
    elseif N.N_SV_tot == 355
        data = load('JPattern_sparse355MatlabVaryProp.mat','JPattern_sparse');%%!!!!!!!!!!!!!!!
    end
    SIM.JPattern = data.JPattern_sparse;


end


%% Some functions to help with calculations
%%%%%%%%%%
function y_out = YfromX(x,z,y_intcep)
    y_out = -1/z*x + y_intcep;
end

%%%%%%%%%%
function res = XfromDesPotential(x,V_des,z,y_intcep,AN,CA)
    % Solve for y: y = -z*x + y_intcep;
    y = YfromX(x,z,y_intcep);
    % Solve An potential
        Eeq_an = AN.EqPotentialHandle(x);

    % Solve Ca potential
        Eeq_ca = CA.EqPotentialHandle(y);

    % Calc residual
        Eeq_cell = Eeq_ca - Eeq_an;
        res = Eeq_cell - V_des;
end