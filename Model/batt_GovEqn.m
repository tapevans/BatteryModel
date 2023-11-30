%% Batt Residual Function
% This function is used to return the time derivative of the governing equations.

function dSVdt = batt_GovEqn(t,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,i_user)
%% Persistent 
    persistent AN__A_c
    persistent AN__A_s
    persistent AN__A_surf_CV_inv
    persistent AN__alpha_a
    persistent AN__alpha_c
    persistent AN__b_thermal
    persistent AN__C_Li_max_inv
    persistent AN__D_oHandle
    persistent AN__del_x_inv
    persistent AN__EqPotentialHandle
    persistent AN__i_0_ref
    persistent AN__i_oHandle
    persistent AN__m_thermal
    persistent AN__r_half_vec_diffcubed_inv
    persistent AN__r_half_vec_sq
    persistent AN__R_SEI_inv

    persistent CA__A_c 
    persistent CA__A_s 
    persistent CA__A_surf_CV_inv
    persistent CA__alpha_a
    persistent CA__alpha_c
    persistent CA__b_thermal
    persistent CA__C_Li_max_inv
    persistent CA__D_oHandle
    persistent CA__del_x_inv
    persistent CA__EqPotentialHandle
    persistent CA__i_0_ref
    persistent CA__i_oHandle
    persistent CA__m_thermal
    persistent CA__r_half_vec_diffcubed_inv
    persistent CA__r_half_vec_sq
    persistent CA__R_SEI_inv

    persistent CONS__BRUG
    persistent CONS__F
    persistent CONS__F_inv
    persistent CONS__R

    persistent EL__ActivityHandle
    persistent EL__C_inv
    persistent EL__D_o_Li_ionHandle
    persistent EL__kappaHandle
    persistent EL__tf_numHandle

    persistent FLAG__Bruggeman
    persistent FLAG__COE
    persistent FLAG__HEAT_GEN_TOTAL 
    persistent FLAG__Newman_i_o
    persistent FLAG__RampThermalGradient
    persistent FLAG__T_BC_AN
    persistent FLAG__T_BC_CA
    persistent FLAG__VARIABLE_activity
    persistent FLAG__VARIABLE_D_Liion
    persistent FLAG__VARIABLE_D_o_AN
    persistent FLAG__VARIABLE_D_o_CA
    persistent FLAG__VARIABLE_kappa
    persistent FLAG__VARIABLE_PROPS_FROM_HANDLES
    persistent FLAG__VARIABLE_tf_num

    persistent N__CV_Region_AN
    persistent N__CV_Region_CA
    persistent N__CV_Region_SEP
    persistent N__IDX_1Dto2D
    persistent N__N_CV_AN
    persistent N__N_CV_CA
    persistent N__N_CV_SEP
    persistent N__N_CV_tot
    persistent N__N_R_AN
    persistent N__N_R_CA
    persistent N__N_R_max
    persistent N__N_SV_AN
    persistent N__N_SV_CA
    persistent N__N_SV_SEP

    persistent P__activity 
    persistent P__C_Li 
    persistent P__C_Li_surf_AN
    persistent P__C_Li_surf_CA
    persistent P__C_Liion 
    persistent P__D_o 
    persistent P__D_o_Li_ion 
    persistent P__del_phi 
    persistent P__i_PS 
    persistent P__kappa 
    persistent P__lambda_eff 
    persistent P__phi_ed 
    persistent P__phi_el
    persistent P__SEP__C_Liion
    persistent P__SEP__phi_el
    persistent P__SEP__T
    persistent P__sigma 
    persistent P__T 
    persistent P__tf_num 
    persistent P__V_1 
    persistent P__V_2 

    persistent SEP__del_x_inv 

    persistent SIM__A_c 
    persistent SIM__A_outside_vec
    persistent SIM__Amp
    persistent SIM__del_CV_r_inv_mat
    persistent SIM__del_x_vec_halved_inv
    persistent SIM__diff_CV_x_vec_inv
    persistent SIM__h
    persistent SIM__h_AN_BC
    persistent SIM__h_CA_BC
    persistent SIM__interp_x_interface
    persistent SIM__MassFluxPreCalcResistor
    persistent SIM__profile_current
    persistent SIM__profile_time
    persistent SIM__q_AN_BC
    persistent SIM__q_CA_BC
    persistent SIM__RampThermalGradientTime
    persistent SIM__SimMode 
    persistent SIM__SV_nan
    persistent SIM__T_inf
    persistent SIM__Temp_AN_BC 
    persistent SIM__Temp_CA_BC 

    % persistent AN__A_surf_CV
    % persistent AN__C_Li_max
    % persistent AN__c_p_eff
    % persistent AN__del_CV_r_vec
    % persistent AN__del_r_vec
    % persistent AN__del_x
    % persistent AN__eps_el
    % persistent AN__lambda_eff
    % persistent AN__r_half_vec
    % persistent AN__R_SEI
    % persistent AN__r_vec
    % persistent AN__rho_eff
    % persistent AN__sigma

    % persistent CA__A_surf_CV
    % persistent CA__C_Li_max
    % persistent CA__c_p_eff
    % persistent CA__del_CV_r_vec
    % persistent CA__del_r_vec
    % persistent CA__del_x 
    % persistent CA__eps_el
    % persistent CA__lambda_eff
    % persistent CA__r_half_vec
    % persistent CA__R_SEI 
    % persistent CA__r_vec
    % persistent CA__rho_eff
    % persistent CA__sigma

    % persistent CONS__R_inv

    % persistent EL__C

    % persistent N__N_prop
    % persistent N__N_SV_AN_tot
    % persistent N__N_SV_max
    % persistent N__N_SV_nR
    % persistent N__N_SV_SEP_tot

    % persistent P__c_p_eff
    % persistent P__R_SEI
    % persistent P__rho_eff

    % persistent SEP__c_p_eff
    % persistent SEP__del_x
    % persistent SEP__eps_el
    % persistent SEP__lambda_eff
    % persistent SEP__rho_eff

    % persistent SIM__del_x_vec_halved
    % persistent SIM__diff_CV_x_vec
    % persistent SIM__eps_el_vec

    if isempty(AN__C_Li_max_inv)
        FLAG__Bruggeman                     = FLAG.Bruggeman;
        FLAG__COE                           = FLAG.COE;
        FLAG__HEAT_GEN_TOTAL                = FLAG.HEAT_GEN_TOTAL;
        FLAG__Newman_i_o                    = FLAG.Newman_i_o;
        FLAG__RampThermalGradient           = FLAG.RampThermalGradient;
        FLAG__T_BC_AN                       = FLAG.T_BC_AN;
        FLAG__T_BC_CA                       = FLAG.T_BC_CA;
        FLAG__VARIABLE_activity             = FLAG.VARIABLE_activity;
        FLAG__VARIABLE_D_Liion              = FLAG.VARIABLE_D_Liion;
        FLAG__VARIABLE_D_o_AN               = FLAG.VARIABLE_D_o_AN;
        FLAG__VARIABLE_D_o_CA               = FLAG.VARIABLE_D_o_CA;
        FLAG__VARIABLE_kappa                = FLAG.VARIABLE_kappa;
        FLAG__VARIABLE_PROPS_FROM_HANDLES   = FLAG.VARIABLE_PROPS_FROM_HANDLES;
        FLAG__VARIABLE_tf_num               = FLAG.VARIABLE_tf_num;

        AN__A_c                     = AN.A_c;
        AN__A_s                     = AN.A_s;
        AN__A_surf_CV_inv           = AN.A_surf_CV_inv;
        AN__alpha_a                 = AN.alpha_a;
        AN__alpha_c                 = AN.alpha_c;
        AN__C_Li_max_inv            = AN.C_Li_max_inv;
        AN__D_oHandle               = AN.D_oHandle;
        AN__del_x_inv               = AN.del_x_inv;
        AN__EqPotentialHandle       = AN.EqPotentialHandle;
        AN__i_0_ref                 = AN.i_0_ref;
        AN__i_oHandle               = AN.i_oHandle;
        AN__r_half_vec_diffcubed_inv= AN.r_half_vec_diffcubed_inv;
        AN__r_half_vec_sq           = AN.r_half_vec_sq;
        AN__R_SEI_inv               = AN.R_SEI_inv;

        CA__A_c                     = CA.A_c;
        CA__A_s                     = CA.A_s;
        CA__A_surf_CV_inv           = CA.A_surf_CV_inv;
        CA__alpha_a                 = CA.alpha_a;
        CA__alpha_c                 = CA.alpha_c;
        CA__C_Li_max_inv            = CA.C_Li_max_inv;
        CA__D_oHandle               = CA.D_oHandle;
        CA__del_x_inv               = CA.del_x_inv;
        CA__EqPotentialHandle       = CA.EqPotentialHandle;
        CA__i_0_ref                 = CA.i_0_ref;
        CA__i_oHandle               = CA.i_oHandle;
        CA__r_half_vec_diffcubed_inv= CA.r_half_vec_diffcubed_inv;
        CA__r_half_vec_sq           = CA.r_half_vec_sq;
        CA__R_SEI_inv               = CA.R_SEI_inv;

        CONS__BRUG                  = CONS.BRUG;
        CONS__F                     = CONS.F;
        CONS__F_inv                 = CONS.F_inv;
        CONS__R                     = CONS.R;

        EL__ActivityHandle          = EL.ActivityHandle;
        EL__C_inv                   = EL.C_inv;
        EL__D_o_Li_ionHandle        = EL.D_o_Li_ionHandle;
        EL__kappaHandle             = EL.kappaHandle;
        EL__tf_numHandle            = EL.tf_numHandle;

        N__CV_Region_AN             = N.CV_Region_AN;
        N__CV_Region_CA             = N.CV_Region_CA;
        N__CV_Region_SEP            = N.CV_Region_SEP;
        N__IDX_1Dto2D               = N.IDX_1Dto2D;
        N__N_CV_AN                  = N.N_CV_AN;
        N__N_CV_CA                  = N.N_CV_CA;
        N__N_CV_SEP                 = N.N_CV_SEP;
        N__N_CV_tot                 = N.N_CV_tot;
        N__N_R_AN                   = N.N_R_AN;
        N__N_R_CA                   = N.N_R_CA;
        N__N_R_max                  = N.N_R_max;
        N__N_SV_AN                  = N.N_SV_AN;
        N__N_SV_CA                  = N.N_SV_CA;
        N__N_SV_SEP                 = N.N_SV_SEP;

        P__activity                 = P.activity; 
        P__C_Li                     = P.C_Li;
        P__C_Li_surf_AN             = P.C_Li_surf_AN;
        P__C_Li_surf_CA             = P.C_Li_surf_CA;
        P__C_Liion                  = P.C_Liion;
        P__D_o                      = P.D_o; 
        P__D_o_Li_ion               = P.D_o_Li_ion; 
        P__del_phi                  = P.del_phi;
        P__i_PS                     = P.i_PS;
        P__kappa                    = P.kappa; 
        P__lambda_eff               = P.lambda_eff; 
        P__phi_ed                   = P.phi_ed; 
        P__phi_el                   = P.phi_el;
        P__SEP__C_Liion             = P.SEP.C_Liion;
        P__SEP__phi_el              = P.SEP.phi_el;
        P__SEP__T                   = P.SEP.T;
        P__sigma                    = P.sigma; 
        P__T                        = P.T;
        P__tf_num                   = P.tf_num; 
        P__V_1                      = P.V_1;
        P__V_2                      = P.V_2;

        SEP__del_x_inv              = SEP.del_x_inv;

        SIM__A_c                    = SIM.A_c;
        SIM__A_outside_vec          = SIM.A_outside_vec;
        SIM__del_CV_r_inv_mat       = SIM.del_CV_r_inv_mat;
        SIM__del_x_vec_halved_inv   = SIM.del_x_vec_halved_inv;
        SIM__diff_CV_x_vec_inv      = SIM.diff_CV_x_vec_inv;
        SIM__h                      = SIM.h;
        SIM__interp_x_interface     = SIM.interp_x_interface;
        SIM__MassFluxPreCalcResistor= SIM.MassFluxPreCalcResistor;
        SIM__profile_current        = SIM.profile_current;
        SIM__profile_time           = SIM.profile_time;
        SIM__q_AN_BC                = SIM.q_AN_BC;
        SIM__q_CA_BC                = SIM.q_CA_BC;
        SIM__RampThermalGradientTime= SIM.RampThermalGradientTime;
        SIM__SimMode                = SIM.SimMode;
        SIM__SV_nan                 = SIM.SV_nan;
        SIM__T_inf                  = SIM.T_inf;

        if isfield(SIM , 'Amp')
            SIM__Amp = SIM.Amp;
        end
        if isfield(SIM , 'Temp_AN_BC')
            SIM__Temp_AN_BC = SIM.Temp_AN_BC;
            SIM__Temp_CA_BC = SIM.Temp_CA_BC;
        end
        if isfield(SIM , 'h_AN_BC')
            SIM__h_AN_BC = SIM.h_AN_BC;
        end
        if isfield(SIM , 'h_CA_BC')
            SIM__h_CA_BC = SIM.h_CA_BC;
        end
        if FLAG__RampThermalGradient
            AN__m_thermal = AN.m_thermal;
            AN__b_thermal = AN.b_thermal;
            CA__m_thermal = CA.m_thermal;
            CA__b_thermal = CA.b_thermal;
        end

        % AN__A_surf_CV = AN.A_surf_CV;
        % AN__C_Li_max = AN.C_Li_max;
        % AN__c_p_eff = AN.c_p_eff;
        % AN__del_CV_r_vec = AN.del_CV_r_vec;
        % AN__del_r_vec = AN.del_r_vec;
        % AN__del_x = AN.del_x;
        % AN__eps_el = AN.eps_el;
        % AN__lambda_eff = AN.lambda_eff;
        % AN__r_half_vec = AN.r_half_vec;
        % AN__R_SEI = AN.R_SEI;
        % AN__r_vec = AN.r_vec;
        % AN__rho_eff = AN.rho_eff;
        % AN__sigma = AN.sigma;

        % CA__A_surf_CV = CA.A_surf_CV;
        % CA__C_Li_max = CA.C_Li_max;
        % CA__c_p_eff = CA.c_p_eff;
        % CA__del_CV_r_vec = CA.del_CV_r_vec;
        % CA__del_r_vec = CA.del_r_vec;
        % CA__del_x = CA.del_x;
        % CA__eps_el = CA.eps_el;
        % CA__lambda_eff = CA.lambda_eff;
        % CA__r_half_vec = CA.r_half_vec;
        % CA__R_SEI = CA.R_SEI; 
        % CA__r_vec = CA.r_vec;
        % CA__rho_eff = CA.rho_eff;
        % CA__sigma = CA.sigma;

        % CONS__R_inv = CONS.R_inv;

        % EL__C = EL.C;

        % N__N_prop = N.N_prop; 
        % N__N_SV_AN_tot = N.N_SV_AN_tot;
        % N__N_SV_max = N.N_SV_max;
        % N__N_SV_nR = N.N_SV_nR;
        % N__N_SV_SEP_tot = N.N_SV_SEP_tot;

        % P__c_p_eff = P.c_p_eff; 
        % P__R_SEI = P.R_SEI; 
        % P__rho_eff = P.rho_eff; 

        % SEP__c_p_eff = SEP.c_p_eff;
        % SEP__del_x = SEP.del_x;
        % SEP__eps_el = SEP.eps_el;
        % SEP__lambda_eff = SEP.lambda_eff;
        % SEP__rho_eff = SEP.rho_eff;
        
        % SIM__del_x_vec_halved = SIM.del_x_vec_halved;
        % SIM__diff_CV_x_vec = SIM.diff_CV_x_vec;
        % SIM__eps_el_vec = SIM.eps_el_vec;
    end


%% Organize (reshape) the SV
    SV = SV1Dto2D_short(SV, SIM__SV_nan, N__IDX_1Dto2D);
    % SV = SV1Dto2D(SV , N__N_SV_max, N__N_CV_tot, N__N_SV_AN_tot, N__N_SV_SEP_tot, N__N_SV_AN, N__N_SV_SEP, N__N_SV_CA, N__N_CV_AN, N__N_CV_SEP, N__N_CV_CA, N__CV_Region_AN, N__CV_Region_SEP, N__CV_Region_CA, P__T, P__del_phi, P__C_Liion, P__SEP__T, P__SEP__phi_el, P__SEP__C_Liion);
    SV = addPhiEl2SV(SV, P__phi_ed, P__del_phi, N__CV_Region_SEP, N__N_CV_SEP);

        
%% Pull out each SV vector
    [T, del_phi, phi_ed, phi_el, V_1, V_2, i_PS, Ce, C_Li, C_Li_AN, C_Li_CA, X_AN, X_CA, Ce_norm, Ce_log, eta, RT_inv_vec] = extractSV(SV,P__T, P__del_phi, P__phi_ed, P__phi_el, P__V_1, P__V_2, P__i_PS, P__C_Liion, P__C_Li, P__C_Li_surf_AN, P__C_Li_surf_CA, N__CV_Region_AN, N__CV_Region_CA, N__N_R_max, AN__C_Li_max_inv, CA__C_Li_max_inv, EL__C_inv, CONS__R);
    

%% Obtain Property Values
    if FLAG__VARIABLE_PROPS_FROM_HANDLES
        props = getProps(  Ce, T, X_AN,  X_CA, FLAG__VARIABLE_kappa, ...
                        FLAG__VARIABLE_D_Liion, FLAG__VARIABLE_activity, FLAG__VARIABLE_tf_num, FLAG__VARIABLE_D_o_AN, FLAG__VARIABLE_D_o_CA, FLAG__Bruggeman,...
                        P__kappa, P__D_o_Li_ion, P__activity, P__tf_num, P__D_o, ...
                        N__N_R_AN, N__N_R_CA, N__CV_Region_AN, N__CV_Region_CA, ...
                        CONS__BRUG, ...
                        AN__D_oHandle, CA__D_oHandle, ...
                        PROPS, EL__kappaHandle, EL__D_o_Li_ionHandle, EL__ActivityHandle, EL__tf_numHandle); 
    else
        props = PROPS;
    end


%% Property Vector
    [sigma_vec, kappa_vec, tf_vec, activity_vec, D_o_Li_ion_vec, lambda_vec, D_o_vec] = extractProps(P__sigma, P__kappa, P__tf_num, P__activity, P__D_o_Li_ion, P__lambda_eff, P__D_o, props);


%% Calculate Value at Interface between CV
    sigma_vec_interface    = getInterfaceXDir(sigma_vec   , SIM__interp_x_interface);
    kappa_vec_interface    = getInterfaceXDir(kappa_vec   , SIM__interp_x_interface);
    tf_vec_interface       = getInterfaceXDir(tf_vec      , SIM__interp_x_interface);
    activity_vec_interface = getInterfaceXDir(activity_vec, SIM__interp_x_interface);
    T_interface            = getInterfaceXDir(T           , SIM__interp_x_interface);

    D_o_vec_interface        = [NaN(1,N__N_CV_tot) ; 0.5*(D_o_vec(1:end-1,:) + D_o_vec(2:end,:)) ; NaN(1,N__N_CV_tot)]; % This assumes constant del_r in both particles

    
%% Calculate Gradient
    [phi_ed_diff, phi_ed_grad] = diffAndGradXCalc(phi_ed, SIM__diff_CV_x_vec_inv);
    [phi_el_diff, phi_el_grad] = diffAndGradXCalc(phi_el, SIM__diff_CV_x_vec_inv);
    [Ce_log_diff, Ce_log_grad] = diffAndGradXCalc(Ce_log, SIM__diff_CV_x_vec_inv);
    [Ce_diff    , Ce_grad    ] = diffAndGradXCalc(Ce    , SIM__diff_CV_x_vec_inv);
    [T_diff     , T_grad     ] = diffAndGradXCalc(T     , SIM__diff_CV_x_vec_inv);

    diff        = C_Li(2:end,:) - C_Li(1:end-1,:);
    C_Li_diff   = [NaN(1 , N__N_CV_tot) ; diff ; NaN(1 , N__N_CV_tot)];


%% Calculate i_user
    if     SIM__SimMode == 3 % State Space EIS
        % Uses the i_user value from the function handle
    elseif SIM__SimMode == 5 % MOO Controller
        % Uses the i_user value from the function handle
    else
        i_user = i_user_calc(t, i_user, SIM__SimMode, SIM__profile_time, SIM__profile_current, SIM__Amp, SIM, FLAG);
    end


%% Calculate All Fluxes
% currentCalc
    [i_ed , i_el ] = currentCalc(  sigma_vec_interface, kappa_vec_interface, T_interface, activity_vec_interface, tf_vec_interface, ...
                                   phi_ed_grad, phi_el_grad, Ce_log_grad, ...
                                   CONS__F_inv, CONS__R, ...
                                   N__CV_Region_AN, N__CV_Region_SEP, N__CV_Region_CA, N__N_CV_tot, i_user);

% iFarCalc
    i_Far = iFarCalc( T, Ce_norm, X_AN, X_CA, eta, RT_inv_vec, AN__i_0_ref , AN__alpha_a, AN__alpha_c, CA__i_0_ref, CA__alpha_a, CA__alpha_c, ...
                        FLAG__Newman_i_o, CONS__F, ...
                        N__CV_Region_AN, N__CV_Region_CA, N__N_CV_SEP, AN__i_oHandle , CA__i_oHandle);

% Species Production
    s_dot = i_Far * CONS__F_inv ;

% JLiionCalc
    J_Liion = JLiionCalc( Ce_diff, D_o_Li_ion_vec, tf_vec_interface, N__CV_Region_AN, N__N_CV_tot, CONS__F_inv, SIM__MassFluxPreCalcResistor, i_el);

% JLiCalc
    J_Li    = JLiCalc( C_Li_diff, D_o_vec_interface, N__N_CV_tot, N__N_R_AN, N__N_R_CA, N__CV_Region_AN, N__CV_Region_CA, SIM__del_CV_r_inv_mat, s_dot);

% COE
    if FLAG__COE
        [q_cond , q_conv , q_gen] = thermalAllVectors( T, T_diff, lambda_vec, N__CV_Region_AN, N__N_CV_tot, SIM__del_x_vec_halved_inv, SIM__h, SIM__T_inf, SIM__q_AN_BC, SIM__q_CA_BC, SIM__h_AN_BC, SIM__h_CA_BC, FLAG__T_BC_AN, FLAG__T_BC_CA, FLAG__HEAT_GEN_TOTAL, i_el, i_ed, i_Far);
    else
        q_cond = zeros(1 , N__N_CV_tot+1);
        q_conv = zeros(1 , N__N_CV_tot  );
        q_gen  = zeros(1 , N__N_CV_tot  );
    end


%% Calculate Equilibrium Voltage
    E_eq_an  = AN__EqPotentialHandle( X_AN(end,:) );
    E_eq_ca  = CA__EqPotentialHandle( X_CA(end,:) );
    
    E_eq_vec = [ E_eq_an , zeros(1,N__N_CV_SEP) , E_eq_ca];


%% Solving for dSVdt
%% ---- Anode ----
    dSVdt_AN = zeros(N__N_SV_AN,N__N_CV_AN);
    for i = 1:N__N_CV_AN
        % Temp
            dSVdt_AN(P__T, i) = -(q_cond(i+1) - q_cond(i)) * SIM__A_c ...
                                - q_conv(i)                * SIM__A_outside_vec(i) ...
                                + q_gen(i) ;
        
        % phi_el                    
            % dSVdt_AN(P__del_phi , i) =  (AN__A_c / AN__A_surf_CV)*(i_el(i+1) - i_el(i)  ) ...
            %                           - (SV(P__V_1,i) - SV(P__phi_el,i))/AN__R_SEI;
            dSVdt_AN(P__del_phi , i) =  (AN__A_c * AN__A_surf_CV_inv)*(i_el(i+1) - i_el(i)  ) ...
                                      - ( V_1(i) - phi_el(i) )*AN__R_SEI_inv;
                            
        % phi_ed                                        
            dSVdt_AN(P__phi_ed , i) =   i_ed(i  ) + i_el(i  )...
                                     - (i_ed(i+1) + i_el(i+1));
        
        % V_1                    
            % dSVdt_AN(P__V_1    , i) = - (SV(P__V_1,i) - SV(P__phi_el,i))/AN__R_SEI ...
            %                           +  i_Far(i);
            dSVdt_AN(P__V_1    , i) = - ( V_1(i) - phi_el(i) )*AN__R_SEI_inv ...
                                      +  i_Far(i);
                            
        % V_2                    
            dSVdt_AN(P__V_2    , i) = - i_Far(i)...
                                      + i_PS(i);
        
	    % i_PS                    
            % dSVdt_AN(P__i_PS   , i) =    SV(P__phi_ed,i) - SV(P__V_2,i) -  E_eq_vec(i);
            dSVdt_AN(P__i_PS   , i) =    phi_ed(i) - V_2(i) -  E_eq_vec(i);
                            
	    % C_Li^+
            % dSVdt_AN(P__C_Liion, i) = -(J_Liion(i+1) - J_Liion(i))/ AN__del_x...
            %                           + s_dot(i) * AN__A_s;
            dSVdt_AN(P__C_Liion, i) = -(J_Liion(i+1) - J_Liion(i)) * AN__del_x_inv...
                                      + s_dot(i) * AN__A_s;
                    
        % C_Li
            for j = 1:N__N_R_AN
                    % dSVdt_AN(P__C_Li+j-1, i) = -3*(AN__r_half_vec(j+1)^2 * J_Li(j+1,i) - AN__r_half_vec(j)^2 * J_Li(j,i)) ...
                    %                            / (AN__r_half_vec(j+1)^3-AN__r_half_vec(j)^3);
                    dSVdt_AN(P__C_Li+j-1, i) = -3*(AN__r_half_vec_sq(j+1) * J_Li(j+1,i) - AN__r_half_vec_sq(j) * J_Li(j,i)) ...
                                               * AN__r_half_vec_diffcubed_inv(j);
            end
    end

    % Fix Boundary Conditions
        % Set Anode to be reference electrode
            i = 1;
            dSVdt_AN(P__phi_ed, i) = -phi_ed(i);

        % If Known Temperature BC
           if FLAG__T_BC_AN == 1 && FLAG__COE 
               i = 1;
               if FLAG__RampThermalGradient && t<= SIM__RampThermalGradientTime
                   Temp_AN_BC = AN__m_thermal * t + AN__b_thermal;
                   dSVdt_AN(P__T, i) = T(i) -     Temp_AN_BC;
               else
                   dSVdt_AN(P__T, i) = T(i) - SIM__Temp_AN_BC;
               end
           end


%% ---- Separator ----
    dSVdt_SEP = zeros(N__N_SV_SEP,N__N_CV_SEP);
    for i = 1:N__N_CV_SEP 
        index_offset = N__N_CV_AN + i;  
        % Temp
            dSVdt_SEP(P__SEP__T, i) = -(q_cond(index_offset+1) - q_cond(index_offset)) * SIM__A_c ...
                                      - q_conv(index_offset)                           * SIM__A_outside_vec(index_offset) ...
                                      + q_gen(index_offset) ;
        
        % phi_el
            dSVdt_SEP(P__SEP__phi_el, i) = -(i_el(index_offset+1)-i_el(index_offset));
        
        % C_Li^+
            % dSVdt_SEP(P__SEP__C_Liion, i)= -(J_Liion(index_offset+1) - J_Liion(index_offset))/SEP__del_x;
            dSVdt_SEP(P__SEP__C_Liion, i)= -(J_Liion(index_offset+1) - J_Liion(index_offset))*SEP__del_x_inv;
    end


%% ---- Cathode ----
    dSVdt_CA = zeros(N__N_SV_CA,N__N_CV_CA);
    for i = 1:N__N_CV_CA
        index_offset = N__N_CV_AN + N__N_CV_SEP + i;  
        % Temp
    	    dSVdt_CA(P__T, i) = -(q_cond(index_offset+1) - q_cond(index_offset)) * SIM__A_c ...
                                - q_conv(index_offset)                           * SIM__A_outside_vec(index_offset) ...
                                + q_gen(index_offset) ;
        
        % phi_el
            % dSVdt_CA(P__del_phi , i) =  (CA__A_c / CA__A_surf_CV)*(i_el(index_offset+1) - i_el(index_offset)  ) ...
            %                             - (SV(P__V_1,index_offset) - SV(P__phi_el,index_offset))/CA__R_SEI;
            dSVdt_CA(P__del_phi , i) =  (CA__A_c * CA__A_surf_CV_inv)*(i_el(index_offset+1) - i_el(index_offset)  ) ...
                                        - ( V_1(index_offset) - phi_el(index_offset) )*CA__R_SEI_inv;
        
	    % phi_ed                    
            dSVdt_CA(P__phi_ed , i)  =   i_ed(index_offset  ) + i_el(index_offset  )...
                                      - (i_ed(index_offset+1) + i_el(index_offset+1));
        
	    % V_1                    
            % dSVdt_CA(P__V_1    , i)  = - (SV(P__V_1,index_offset) - SV(P__phi_el,index_offset))/CA__R_SEI ...
            %                            +  i_Far(index_offset);
            dSVdt_CA(P__V_1    , i)  = - ( V_1(index_offset) - phi_el(index_offset) )*CA__R_SEI_inv ...
                                       +  i_Far(index_offset);
        
	    % V_2                    
            dSVdt_CA(P__V_2    , i)  = -  i_Far(index_offset)...
                                       +  i_PS(index_offset);
        
	    % i_PS                    
            dSVdt_CA(P__i_PS   , i)  =    phi_ed(index_offset) - V_2(index_offset) ...
                                       -  E_eq_vec(index_offset);
        
        % C_Li^+
            % dSVdt_CA(P__C_Liion, i) = -(J_Liion(index_offset+1) - J_Liion(index_offset))/ CA__del_x...
            %                           + s_dot(index_offset) * CA__A_s;
            dSVdt_CA(P__C_Liion, i) = -(J_Liion(index_offset+1) - J_Liion(index_offset)) * CA__del_x_inv...
                                      + s_dot(index_offset) * CA__A_s;
        % C_Li
            for j = 1:N__N_R_CA
                    % dSVdt_CA(P__C_Li+j-1, i) = -3*(CA__r_half_vec(j+1)^2 * J_Li(j+1,index_offset) - CA__r_half_vec(j)^2 * J_Li(j,index_offset)) ...
                    %                            / (CA__r_half_vec(j+1)^3-CA__r_half_vec(j)^3);
                    dSVdt_CA(P__C_Li+j-1, i) = -3*(CA__r_half_vec_sq(j+1) * J_Li(j+1,index_offset) - CA__r_half_vec_sq(j) * J_Li(j,index_offset)) ...
                                               * CA__r_half_vec_diffcubed_inv(j);
            end
    end

    % Fix Boundary Conditions
        % If Known Temperature BC
           if FLAG__T_BC_CA == 1 && FLAG__COE 
               i = N__N_CV_CA;
               index_offset     = N__N_CV_AN + N__N_CV_SEP + i;
               if FLAG__RampThermalGradient && t<= SIM__RampThermalGradientTime
                   Temp_CA_BC = CA__m_thermal * t + CA__b_thermal;
                   dSVdt_CA(P__T, i) = T(index_offset) -      Temp_CA_BC;
               else
                   dSVdt_CA(P__T, i) = T(index_offset) - SIM__Temp_CA_BC;
               end
           end


%% Reshape
% Reshape matrix to a column vector
    dSVdt_AN  = reshape(dSVdt_AN ,[],1);
    dSVdt_SEP = reshape(dSVdt_SEP,[],1);
    dSVdt_CA  = reshape(dSVdt_CA ,[],1);


%% Combine all dSVdt from each region
    dSVdt = [dSVdt_AN ; dSVdt_SEP ; dSVdt_CA];


%% Used for troubleshooting
    % if t>SIM.initial_offset
    %    t;
    % end
    % if t>SIM.t_ramp
    %    t;
    % end
end