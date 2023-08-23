%% Batt Initialization
% This function is used to initialize all of the user inputs and calculates
% the necessary parameters and pointers for the model

function [AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS] = batt_init(AN,CA,SEP,EL,SIM,N,FLAG)
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
    P.phi_el    = P.C_Li_surf_max + 1;
    
% Pointers for Electrostatic 
    i = 1;
    P.ES.del_phi = i; i = i + 1;
    P.ES.phi_ed  = i; i = i + 1;
    P.ES.V_1     = i; i = i + 1;
    P.ES.V_2     = i; i = i + 1;
    P.ES.i_PS    = i; i = i + 1;
    P.ES.i_dl    = i; i = i + 1;
    
    N.N_ES_var  = 6;

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
    P.k          = i; i = i + 1;
    P.rho        = i; i = i + 1;
    P.c_p        = i; i = i + 1;
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
    
    N.N_Out = length(fieldnames(P.OM));


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
        SIM.x_vec      = [AN.x_vec     , SEP.x_vec            , CA.x_vec            ]; % [m], Position of the control volume's center
        SIM.x_half_vec = [AN.x_half_vec, SEP.x_half_vec(2:end), CA.x_half_vec(2:end)]; % [m], Position of the control volume's edges
    
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
    PROPS( P.D_o:(P.D_o+N.N_R_AN-1)    , N.CV_Region_AN ) = AN.D_o        * ones( N.N_R_AN , N.N_CV_AN );
    PROPS( P.D_o:(P.D_o+N.N_R_CA-1)    , N.CV_Region_CA ) = CA.D_o        * ones( N.N_R_CA , N.N_CV_CA );
    
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

% Initial Thermal Gradient
    if FLAG.InitialThermalGradient
        L_tot = SIM.x_vec(end) - SIM.x_vec(1);
        Delta_temp = SIM.CA_Temp - SIM.AN_Temp;
        Slope = Delta_temp / L_tot;
        temp_fnc = @(x) Slope.*x + SIM.AN_Temp;
        Temp_vec = temp_fnc(SIM.x_vec);
    else
        Temp_vec = SIM.Temp_start * ones( 1 , N.N_CV_tot );
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
    AN.A_surf_tot = AN.A_surf * AN.Np_tot;     % [m^2],        Total surface area in the electrode region
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
    

%% Capacity Calculations 
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
    AN.Cap = AN.C_Li_max * CONS.F *(1/3600)*AN.V_ed; % [Ahr], Method 2: Based on a modeled parameter of C_max
    CA.Cap = CA.C_Li_max * CONS.F *(1/3600)*CA.V_ed; % [Ahr], Method 2: Based on a modeled parameter of C_max

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
        voltage_ini_an   = AN.EqPotentialHandle(SIM.x_ini);
        voltage_ini_ca   = CA.EqPotentialHandle(SIM.y_ini);
        voltage_ini_cell = voltage_ini_ca - voltage_ini_an;

    % Calculate electrode initial active material concentration
        concen_ini_an = AN.C_Li_max * SIM.x_ini;
        concen_ini_ca = CA.C_Li_max * SIM.y_ini;
    
    
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
    SIM.OutputMatrix = zeros(N.N_Out , N.N_SV_tot);
        % Cell Voltage
            idx_phi_ed_AN = P.phi_ed;

            i = N.N_CV_CA(end);
            index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
            idx_phi_ed_CA = index_offset + P.phi_ed;

            SIM.OutputMatrix(P.OM.cell_volt,idx_phi_ed_AN) = -1;
            SIM.OutputMatrix(P.OM.cell_volt,idx_phi_ed_CA) =  1;
        % @AN/SEP
            i = N.N_CV_AN(end);
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
                SV_IC(index_offset+P.C_Li+j-1)   = concen_ini_an;
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
                SV_IC(index_offset+N.N_SV_nR+j)   = concen_ini_ca;
            end
    end


%% Save the Equilibrium Values of the Output
    SIM.OutputAtEquil = SIM.OutputMatrix*SV_IC;
    SIM.SV_IC_Static  = SV_IC;


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
        i_user = i_user_calc(0,SIM,FLAG,i_user_in); 
    end

% Create phi vector    
    phi_temp = SV(P.del_phi:P.i_PS, :);
    phi_temp(end+1,: ) = zeros(1,N.N_CV_tot ); % Adding a row for the new constraint equation(i_dl)
    phi_temp(end,N.CV_Region_SEP) = NaN(1,N.N_CV_SEP); % No SV in the SEP region for the new constraint
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


%% Mass Matrix
M = zeros(N.N_SV_tot,N.N_SV_tot);
sim_cap = 0; % A dummy capacitance to help DAE solver by turning algebraics into ODEs. With solver working, this isn't needed anymore

% ---- Anode ----
    for i = 1:N.N_CV_AN
        index_offset = (i-1)*N.N_SV_AN;
        % Temp
            M(index_offset+P.T       , index_offset+P.T )        =  1;
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
            M(index_offset+P.SEP.T       , index_offset+P.SEP.T )        =  1;
        % phi_el
            M(index_offset+P.SEP.phi_el  , index_offset+P.SEP.phi_el )   =  sim_cap;
        % C_Li^+
            M(index_offset+P.SEP.C_Liion , index_offset+P.SEP.C_Liion )  =  SEP.eps_el; 
    end

% ---- Cathode ----
    for i = 1:N.N_CV_CA
        index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
        % Temp
            M(index_offset+P.T       , index_offset+P.T )        =  1;
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


    SIM.M = M;


%% Make Mass for just diffEq
% if SIM.SimMode == 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%% TestPurposes
    % Make indices vector
    idx_diff = [];
    idx_algb = [];
    % ---- Anode ----
    for i = 1:N.N_CV_AN
        index_offset = (i-1)*N.N_SV_AN;
        % Temp
            idx_diff(end+1) = index_offset + P.T;
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
            idx_diff(end+1) = index_offset + P.T;
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

%     SIM.M_DiffEq = M(idx_diff,idx_diff);
%     SIM.M_DiffEq_Inv = inv(SIM.M_DiffEq);
    SIM.diff_idx = idx_diff;
    SIM.algb_idx = idx_algb;

    N.N_diff = length(SIM.diff_idx);
    N.N_algb = length(SIM.algb_idx);
% end


%% Set up Jacobian here too (sparse)
%%%%% Later

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


%% OOOOOOOLLLLLDDDD Stuff

% AN.Cap = (CONS.F/(AN.MW*3600))*AN.V_ed*AN.rho; % [Ahr], Method 1: More theoretical approach

    % AN.Cap_test3 = AN.specCap*AN.V_ed*AN.rho;      % [Ahr], Method 3: Based on something that is more measureable

    % CA.Cap = (CONS.F/(CA.MW*3600))*CA.V_ed*CA.rho; % [Ahr], Method 1: More theoretical approach

    % CA.Cap = CA.specCap*CA.V_ed*CA.rho;            % [Ahr], Method 3: Based on something that is more measureable
    % CA.Cap = CA.Cap * 0.6;                         % [Ahr], Modified Method 2, assumes only 60 percent of the capacity is being used

    % %%%%%%%Testing
    % x_min = 0.0592;
    % x_max = 0.9815;
    % y_intcep = 0.8612;
    % 
    % % Initial Lithiation Fraction
    % SIM.x_ini = (x_max-x_min)*SIM.SOC_start/100 + x_min;
    % SIM.y_ini = YfromX(SIM.x_ini,z,y_intcep);
    % 
    % % Determine electrode initial voltage potential
    % voltage_ini_an = AN.EqPotentialHandle(SIM.x_ini);
    % voltage_ini_ca = CA.EqPotentialHandle(SIM.y_ini);
    % voltage_ini_cell = voltage_ini_ca - voltage_ini_an;
    % 
    % % Calculate electrode initial active material concentration
    % concen_ini_an = AN.C_Li_max * SIM.x_ini;
    % concen_ini_ca = CA.C_Li_max * SIM.y_ini;
    
    % %%%%%%%Testing
    % SOC = 0:1:100;
    % x_ini = (SIM.x_max-SIM.x_min)*SOC/100 + SIM.x_min;
    % y_ini = YfromX(x_ini,z,y_intcep);
    % voltage_an = AN.EqPotentialHandle(x_ini);
    % voltage_ca = CA.EqPotentialHandle(y_ini);
    % cell_voltage = voltage_ca - voltage_an;
    % 
    % %%%%%%
    % figure
    % title('Cell Voltage and X_{surf} vs SOC')
    % xlabel('SOC (%)')
    % 
    % yyaxis left
    % plot(SOC , cell_voltage , 'Linewidth' , 2 )
    % ylabel('Cell Voltage (V)')
    % 
    % yyaxis right
    % plot(SOC , x_ini, 'Linewidth' , 2 )
    % ylabel('X_surf')
    % 
    % %%%%%
    % figure
    % title('X_{surf} vs SOC')
    % xlabel('SOC (%)')
    % 
    % yyaxis left
    % plot(SOC , x_ini , 'Linewidth' , 2 )
    % ylabel('Anode')
    % 
    % yyaxis right
    % plot(SOC , y_ini, 'Linewidth' , 2 )
    % ylabel('Cathode')
    
% User Time Vector and Current
% % ---- Polarization ----
%     if SIM.SimMode == 1
%         % SIM.A_c        = min( AN.A_c, CA.A_c );
%         % SIM.Cell_Cap   = min( AN.Cap, CA.Cap );
%         SIM.I_user_amp = SIM.C_rate*SIM.Cell_Cap*SIM.ChargeOrDischarge; % [A], Load current
%         SIM.i_user_amp = SIM.I_user_amp/SIM.A_c; % [A m^-2], Load current density (flux)
% 
% % ---- Harmonic Perturbation ----
%     elseif SIM.SimMode == 2
%         % SIM.A_c        = min( AN.A_c, CA.A_c );
%         % SIM.Cell_Cap   = min( AN.Cap, CA.Cap );
%         SIM.I_user_amp = SIM.C_rate*SIM.Cell_Cap;
%         SIM.i_user_amp = SIM.I_user_amp/SIM.A_c;
% 
% % ---- State Space EIS ----
%     elseif SIM.SimMode == 3
%         % SIM.A_c        = min( AN.A_c, CA.A_c );
%         % SIM.Cell_Cap   = min( AN.Cap, CA.Cap );
%         SIM.I_user_amp = SIM.C_rate*SIM.Cell_Cap;
%         SIM.i_user_amp = SIM.I_user_amp/SIM.A_c;
% 
% % ---- Known BC Profile Controller ----
%     elseif SIM.SimMode == 4 
%         % SIM.A_c        = min( AN.A_c, CA.A_c );
%         % SIM.Cell_Cap   = min( AN.Cap, CA.Cap );
% 
% % ---- MOO Controller ----
%     elseif SIM.SimMode == 5
%         % SIM.A_c        = min( AN.A_c, CA.A_c );
%         % SIM.Cell_Cap   = min( AN.Cap, CA.Cap );
%         %!!!!!!!! May need to do something here for SV initialization with a
%         %proper i_user value
% 
% % ---- Manual Current Profile ----
%     elseif SIM.SimMode == 7 
%         % SIM.A_c        = min( AN.A_c, CA.A_c );
%         % SIM.Cell_Cap   = min( AN.Cap, CA.Cap );
%         SIM.I_user_amp = SIM.C_rate*SIM.Cell_Cap*SIM.ChargeOrDischarge;
%         SIM.i_user_amp = SIM.I_user_amp/SIM.A_c;
%         SIM.i_user_OG  = SIM.i_user_amp;
% 
% % ---- PRBS ----
%     elseif SIM.SimMode == 8
%         % SIM.A_c        = min( AN.A_c, CA.A_c );
%         % SIM.Cell_Cap   = min( AN.Cap, CA.Cap );
%         % N_PRBS  = 255;                      % Length of the input           [-]
%         N_PRBS  = 2^9 - 1;                  % Length of the input           [-] % New PRBS (zero-mean)
%         TYPE    = 'PRBS';                   % Create a PRBS signal          [-]
%         BAND    = [0 1];                    % Freq. band                    [s]
%         LEVELS  = [-1 1];                   % PRBS limits                   [-]
%         PRBS_gen= idinput(N_PRBS,TYPE,BAND,LEVELS);         % PRBS Demand   [A/m^2]
%         % t_Demand = ((0:1:SIM.PRBSLength-1)*(SIM.Tswitch))'; % PRBS Demand t [s]
%         % C_Demand = PRBS_gen(1:SIM.PRBSLength)*SIM.PRBSAmp;  % PRBS Demand C [A/m^2]
%         % t_Demand = ((0:1:SIM.PRBSLength-9)*(SIM.Tswitch))'; % PRBS Demand t [s]     % New PRBS (zero-mean)
%         C_Demand = PRBS_gen(9:SIM.PRBSLength)*SIM.PRBSAmp;  % PRBS Demand C [A/m^2] % New PRBS (zero-mean)
% 
%         if SIM.MakeLongPRBSSignal
%             num_repeats = ceil(SIM.DesiredLength/length(C_Demand));
%             C_Demand_single = C_Demand;
%             for i = 1:num_repeats-1
%                 if mod(i,2)==1
%                     C_Demand = [C_Demand ; -C_Demand_single];
%                 else
%                     C_Demand = [C_Demand ;  C_Demand_single];
%                 end
%             end
%             C_Demand = C_Demand(1:SIM.DesiredLength);
%         end
% 
%         t_Demand = ( ((1:length(C_Demand))-1 )*(SIM.Tswitch))' ;
% 
% 
%         % Append a no-current entry
%         if SIM.initial_offset > 0
%             t_Demand = [0; t_Demand+SIM.initial_offset]; % PRBS & entry [s]
%             C_Demand = [0; C_Demand];                    % PRBS & entry [A/m^2]
%         end
% 
%         % Add Rest between some switches to help with stability
%             if SIM.AddIntermediateRelaxTime
%                 % Find Zero-Crossings
%                     change_idx = find( C_Demand(1:end-1)~=C_Demand(2:end) ); % True when the following idx doesn't match
%                     change_idx = change_idx + 1; % Change occurs at the value now
% 
%                 % Drop first change
%                     change_idx = change_idx(2:end);
% 
%                 % Find Insert index
%                     mod_idx    = (SIM.NumZeroCrossingUntilNextRelax : SIM.NumZeroCrossingUntilNextRelax : length(change_idx))';
%                     insert_idx = change_idx(mod_idx);
% 
%                 % Create New Current Demand Vector
%                     C_new     = C_Demand(1:insert_idx(1)-1);
%                     N_changes = length(insert_idx);
%                     for i = 1:N_changes-1
%                         C_new = [C_new ; zeros(SIM.NumTsRelax,1) ; C_Demand(insert_idx(i):insert_idx(i+1)-1)];
%                     end
%                     C_new = [ C_new ; zeros(SIM.NumTsRelax,1) ; C_Demand(insert_idx(N_changes):end) ];
% 
%                 % Create New Time Vector
%                     t_new = (0:1:length(C_new)-1)*SIM.Tswitch;
%             else
%                 C_new = C_Demand;
%                 t_new = t_Demand;
%             end
% 
%         % Add Ramp Times
%         SIM.profile_time    = t_new(1);
%         SIM.profile_current = C_new(1);
% 
%         % Step at t^+
%         for i = 2:length(t_new)
%             % @ k
%             SIM.profile_time(end+1,1)    = t_new(i);
%             SIM.profile_current(end+1,1) = C_new(i-1);
% 
%             % After k
%             SIM.profile_time(end+1,1)    = t_new(i) + SIM.Tswitch * SIM.t_ramp_ratio;
%             SIM.profile_current(end+1,1) = C_new(i);
%         end
%         SIM.profile_time    = SIM.profile_time(1:end-1);
%         SIM.profile_current = SIM.profile_current(1:end-1);
% 
%         % % Test Plot
%         %     t_test = 0 : (SIM.Tswitch/5) : SIM.Tswitch*5;
%         %     t_test = t_test + 2* SIM.Tswitch * SIM.t_ramp_ratio; % Shift samples slightly
%         %     i_user = i_user_calc( t_test , SIM);
%         %     figure
%         %     hold on
%         %     plot(t_test           , i_user              ,'ro','Linewidth',2,'DisplayName','Sampled') % This is what the simulation will actually use
%         %     plot(SIM.profile_time , SIM.profile_current ,'-k','Linewidth',2,'DisplayName','Data Points')
%         %     lgn = legend;
%         %     xlim([0,t_test(end)+SIM.Tswitch])
% 
% % ---- EIS from Stitching PRBS ----
%     elseif SIM.SimMode == 9
% 
% % ---- EIS Ho-Kalman ----
%     elseif SIM.SimMode == 10
%         % SIM.A_c        = min( AN.A_c, CA.A_c );
%         % SIM.Cell_Cap   = min( AN.Cap, CA.Cap );
%         SIM.i_user_amp = 1;
%         SIM.I_user_amp = SIM.i_user_amp * SIM.A_c;
% 
%         t_Demand    = (0:1:SIM.HK_nSamples+3) * SIM.Tsample ; % +3 is for (2 inital relax, 1 pulse)
%         C_Demand    = zeros( size(t_Demand) );
%         C_Demand(3) = SIM.i_user_amp;
% 
%         % Step at t^+
%         SIM.profile_time    = t_Demand(1);
%         SIM.profile_current = C_Demand(1);
% 
%         for i = 2:length(t_Demand)
%             % @ k
%             SIM.profile_time(end+1,1)    = t_Demand(i);
%             SIM.profile_current(end+1,1) = C_Demand(i-1);
% 
%             % After k
%             SIM.profile_time(end+1,1)    = t_Demand(i) + SIM.Tsample * SIM.t_ramp_ratio;
%             SIM.profile_current(end+1,1) = C_Demand(i);
%         end
%         SIM.profile_time    = SIM.profile_time(1:end-1);
%         SIM.profile_current = SIM.profile_current(1:end-1);
% 
%         % % Test Plot
%         %     t_test = 0 : (SIM.Tsample/5) : SIM.Tsample*5;
%         %     t_test = t_test + 2* SIM.Tsample * SIM.t_ramp_ratio; % Shift samples slightly
%         %     i_user = i_user_calc( t_test , SIM);
%         %     figure
%         %     hold on
%         %     plot(t_test           , i_user              ,'ro','Linewidth',2,'DisplayName','Sampled') % This is what the simulation will actually use
%         %     plot(SIM.profile_time , SIM.profile_current ,'-k','Linewidth',2,'DisplayName','Data Points')
%         %     lgn = legend;
%         %     xlim([0,t_test(end)+SIM.Tsample])
%     end

% Determine Simulation Time Vector
% % ---- Polarization ----
%     if SIM.SimMode == 1 
%         if SIM.C_rate == 0
%             t_final = 30; % [s], Final time
%         else
%             t_final = SIM.charge_frac*3600/SIM.C_rate + SIM.initial_offset + SIM.t_ramp; % [s], Final time
%         end
%     %     SIM.tspan = [0,SIM.initial_offset, t_final];
%         SIM.tspan = [0, t_final];
% 
% % ---- Harmonic Perturbation ----
%     elseif SIM.SimMode == 2
%         SIM.f     = SIM.freq / (2*pi);    % [cycles s^-1],  Time to complete a period or a full sinusoid cycle
%         SIM.f_s   = (2*SIM.f)*20;         % [samples s^-1], Based on Nyquist sampling theory (minimum is 2*f)
%         t_sim     = SIM.N_cycles / SIM.f; % [s],            Time required for N_cycles of sinusoidal curves
%         N_samples = SIM.f_s*t_sim;        % [samples],      Number of samples total
%         del_time  = t_sim / N_samples;    % [s sample^-1],  Time between each sample %%%%%%% is this also 1/f_s??????????
% 
%         sin_time_vector = 0 : del_time : (t_sim-del_time);                     % Vector of time spaces
%         sin_time_vector = sin_time_vector + SIM.initial_offset; % Add the initial offset
%         if SIM.initial_offset == 0
%             SIM.tspan = sin_time_vector;                        % Overall simulation vector
%         else
%             SIM.tspan = [0, sin_time_vector];                   % Overall simulation vector
%         end
% 
% % ---- State Space EIS ----
%         % No t_span needed for this mode
% 
% % ---- Known BC Profile Controller ----
%     elseif SIM.SimMode == 4
%         % Determine inside RunSimulation from MO_List
% 
% % ---- MOO Controller ----
%     elseif SIM.SimMode == 5
%         % Determine inside RunSimulation from MO_List
% 
% % ---- Manual Current Profile ----
%     elseif SIM.SimMode == 7 % Manual Current Profile
%         % Determine inside RunSimulation from MO_List
% 
% % ---- PRBS ----
%     elseif SIM.SimMode == 8 
%         SIM.tspan = [0, SIM.profile_time(end)];
% 
% % ---- EIS from Stitching PRBS ----
%     elseif SIM.SimMode == 9
%         % SIM.tspan = [0, SIM.profile_time(end)];
% 
% % ---- EIS Ho-Kalman ----
%     elseif SIM.SimMode == 10
%         SIM.tspan = [0, SIM.profile_time(end)];
%     end