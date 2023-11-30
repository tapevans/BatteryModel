%% PersistentFunction
function [AN__C_Li_max, AN__C_Li_max_inv, AN__A_surf_CV_inv, AN__R_SEI, AN__R_SEI_inv, AN__del_x_inv, AN__A_s, AN__r_half_vec_sq, AN__A_c, AN__alpha_a, AN__alpha_c, AN__r_half_vec_diffcubed_inv, AN__i_0_ref, CA__C_Li_max, CA__C_Li_max_inv, CA__A_surf_CV_inv, CA__R_SEI, CA__R_SEI_inv, CA__del_x_inv, CA__A_s , CA__r_half_vec_sq, CA__A_c , CA__alpha_a, CA__alpha_c, CA__r_half_vec_diffcubed_inv, CA__i_0_ref, SEP__del_x_inv , SIM__SimMode , SIM__A_c , SIM__A_outside_vec, SIM__Temp_AN_BC , SIM__h, SIM__T_inf, SIM__h_AN_BC, SIM__q_AN_BC, SIM__h_CA_BC, SIM__q_CA_BC, SIM__Amp, SIM__profile_time, SIM__profile_current, SIM__del_CV_r_inv_mat, P__T , P__del_phi , P__phi_ed, P__V_1, P__V_2 , P__i_PS      , P__C_Liion   , P__C_Li      , P__C_Li_surf_AN, P__C_Li_surf_CA, P__SEP__T, P__SEP__phi_el, P__SEP__C_Liion, P__phi_el, P__sigma      , P__kappa      , P__R_SEI      , P__lambda_eff , P__rho_eff    , P__c_p_eff    , P__D_o_Li_ion , P__activity   , P__tf_num     , P__D_o        , N__N_CV_tot, N__CV_Region_AN, N__CV_Region_SEP, N__N_SV_AN, N__N_CV_AN, N__N_R_AN, N__CV_Region_CA, N__N_SV_CA, N__N_CV_CA, N__N_R_CA, N__N_CV_SEP, N__N_SV_SEP, N__N_R_max, N__N_prop, N__N_SV_AN_tot, N__N_SV_SEP_tot, FLAG__VARIABLE_PROPS_FROM_HANDLES, FLAG__COE, FLAG__HEAT_GEN_TOTAL , FLAG__T_BC_AN, FLAG__T_BC_CA, FLAG__Newman_i_o, FLAG__Bruggeman, FLAG__VARIABLE_kappa, FLAG__VARIABLE_D_Liion, FLAG__VARIABLE_activity, FLAG__VARIABLE_tf_num, FLAG__VARIABLE_D_o_AN, FLAG__VARIABLE_D_o_CA, AN__sigma, CA__sigma, AN__lambda_eff, AN__rho_eff, AN__c_p_eff, EL__C_inv, SEP__lambda_eff, SEP__rho_eff, SEP__c_p_eff, CA__lambda_eff, CA__rho_eff, CA__c_p_eff, CONS__F, CONS__R, CONS__F_inv, CONS__BRUG, SIM__RampThermalGradientTime, AN__m_thermal, AN__b_thermal, CA__m_thermal, CA__b_thermal, SIM__diff_CV_x_vec_inv, SIM__interp_x_interface, SIM__MassFluxPreCalcResistor, SIM__del_x_vec_halved_inv, FLAG__RampThermalGradient] = getPersistentVariables(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG)

%% Persistent 
    persistent AN__C_Li_max
    persistent AN__C_Li_max_inv
    persistent AN__A_surf_CV_inv
    persistent AN__R_SEI
    persistent AN__R_SEI_inv
    persistent AN__del_x_inv
    persistent AN__A_s
    persistent AN__r_half_vec_sq
    persistent AN__A_c
    persistent AN__alpha_a
    persistent AN__alpha_c
    persistent AN__r_half_vec_diffcubed_inv
    persistent AN__i_0_ref
    persistent CA__C_Li_max
    persistent CA__C_Li_max_inv
    persistent CA__A_surf_CV_inv
    persistent CA__R_SEI
    persistent CA__R_SEI_inv
    persistent CA__del_x_inv
    persistent CA__A_s 
    persistent CA__r_half_vec_sq
    persistent CA__A_c 
    persistent CA__alpha_a
    persistent CA__alpha_c
    persistent CA__r_half_vec_diffcubed_inv
    persistent CA__i_0_ref
    persistent SEP__del_x_inv 
    persistent SIM__SimMode 
    persistent SIM__A_c 
    persistent SIM__A_outside_vec
    persistent SIM__Temp_AN_BC 
    persistent SIM__h
    persistent SIM__T_inf
    persistent SIM__h_AN_BC
    persistent SIM__q_AN_BC
    persistent SIM__h_CA_BC
    persistent SIM__q_CA_BC
    persistent SIM__Amp
    persistent SIM__profile_time
    persistent SIM__profile_current
    persistent SIM__del_CV_r_inv_mat
    persistent P__T 
    persistent P__del_phi   
    persistent P__phi_ed    
    persistent P__V_1       
    persistent P__V_2       
    persistent P__i_PS      
    persistent P__C_Liion   
    persistent P__C_Li      
	persistent P__C_Li_surf_AN
	persistent P__C_Li_surf_CA
    persistent P__SEP__T
    persistent P__SEP__phi_el
    persistent P__SEP__C_Liion
    persistent P__phi_el
    persistent P__sigma      
    persistent P__kappa      
    persistent P__R_SEI      
    persistent P__lambda_eff 
    persistent P__rho_eff    
    persistent P__c_p_eff    
    persistent P__D_o_Li_ion 
    persistent P__activity   
    persistent P__tf_num     
    persistent P__D_o        
    persistent N__N_CV_tot
    persistent N__CV_Region_AN
    persistent N__CV_Region_SEP
    persistent N__N_SV_AN
    persistent N__N_CV_AN
    persistent N__N_R_AN
    persistent N__CV_Region_CA
    persistent N__N_SV_CA
    persistent N__N_CV_CA
    persistent N__N_R_CA
    persistent N__N_CV_SEP
    persistent N__N_SV_SEP
    persistent N__N_R_max
    persistent N__N_prop
    persistent N__N_SV_AN_tot
    persistent N__N_SV_SEP_tot
    persistent FLAG__VARIABLE_PROPS_FROM_HANDLES
    persistent FLAG__COE
    persistent FLAG__HEAT_GEN_TOTAL 
    persistent FLAG__T_BC_AN
    persistent FLAG__T_BC_CA
    persistent FLAG__Newman_i_o
    persistent FLAG__Bruggeman
    persistent FLAG__VARIABLE_kappa
    persistent FLAG__VARIABLE_D_Liion
    persistent FLAG__VARIABLE_activity
    persistent FLAG__VARIABLE_tf_num
    persistent FLAG__VARIABLE_D_o_AN
    persistent FLAG__VARIABLE_D_o_CA
    persistent AN__sigma
    persistent CA__sigma
    persistent AN__lambda_eff
    persistent AN__rho_eff
    persistent AN__c_p_eff
    persistent EL__C_inv
    persistent SEP__lambda_eff
    persistent SEP__rho_eff
    persistent SEP__c_p_eff
    persistent CA__lambda_eff
    persistent CA__rho_eff
    persistent CA__c_p_eff
    persistent CONS__F
    persistent CONS__R
    persistent CONS__F_inv
    persistent CONS__BRUG
    persistent SIM__RampThermalGradientTime
    persistent AN__m_thermal
    persistent AN__b_thermal
    persistent CA__m_thermal
    persistent CA__b_thermal
    persistent SIM__diff_CV_x_vec_inv
    persistent SIM__interp_x_interface
    persistent SIM__MassFluxPreCalcResistor
    persistent SIM__del_x_vec_halved_inv
    persistent FLAG__RampThermalGradient
    % persistent SIM__eps_el_vec
    % persistent AN__del_r_vec
    % persistent CA__del_r_vec
    % persistent AN__r_vec
    % persistent AN__eps_el
    % persistent AN__del_CV_r_vec
    % persistent CA__r_vec
    % persistent CA__eps_el
    % persistent CA__del_CV_r_vec
    % persistent SEP__eps_el
    % persistent N__N_SV_nR

    % persistent AN__A_surf_CV
    % persistent AN__del_x
    % persistent CA__A_surf_CV
    % persistent CA__del_x 
    % persistent SEP__del_x
    % persistent EL__C
    % persistent CONS__R_inv
    % persistent SIM__diff_CV_x_vec
    % persistent SIM__del_x_vec_halved
    % persistent AN__r_half_vec
    % persistent CA__r_half_vec

    if isempty(AN__C_Li_max)
        FLAG__T_BC_AN = FLAG.T_BC_AN;
        FLAG__T_BC_CA = FLAG.T_BC_CA;
        FLAG__RampThermalGradient = FLAG.RampThermalGradient;
        FLAG__VARIABLE_PROPS_FROM_HANDLES = FLAG.VARIABLE_PROPS_FROM_HANDLES;
        FLAG__COE = FLAG.COE;
        FLAG__HEAT_GEN_TOTAL  = FLAG.HEAT_GEN_TOTAL;
        FLAG__Newman_i_o = FLAG.Newman_i_o;
        FLAG__Bruggeman = FLAG.Bruggeman;
        FLAG__VARIABLE_kappa = FLAG.VARIABLE_kappa;
        FLAG__VARIABLE_D_Liion = FLAG.VARIABLE_D_Liion;
        FLAG__VARIABLE_activity = FLAG.VARIABLE_activity;
        FLAG__VARIABLE_tf_num = FLAG.VARIABLE_tf_num;
        FLAG__VARIABLE_D_o_AN = FLAG.VARIABLE_D_o_AN;
        FLAG__VARIABLE_D_o_CA = FLAG.VARIABLE_D_o_CA;

        AN__sigma = AN.sigma;
        CA__sigma = CA.sigma;
        AN__lambda_eff = AN.lambda_eff;
        AN__rho_eff = AN.rho_eff;
        AN__c_p_eff = AN.c_p_eff;
        SEP__lambda_eff = SEP.lambda_eff;
        SEP__rho_eff = SEP.rho_eff;
        SEP__c_p_eff = SEP.c_p_eff;
        CA__lambda_eff = CA.lambda_eff;
        CA__rho_eff = CA.rho_eff;
        CA__c_p_eff = CA.c_p_eff;

        AN__C_Li_max = AN.C_Li_max;
        AN__C_Li_max_inv = AN.C_Li_max_inv;
        AN__A_surf_CV_inv = AN.A_surf_CV_inv;
        AN__R_SEI = AN.R_SEI;
        AN__R_SEI_inv = AN.R_SEI_inv;
        AN__del_x_inv = AN.del_x_inv;
        AN__A_s = AN.A_s;
        AN__r_half_vec_sq = AN.r_half_vec_sq;
        AN__A_c = AN.A_c;
        AN__alpha_a = AN.alpha_a;
        AN__alpha_c = AN.alpha_c;
        AN__i_0_ref = AN.i_0_ref;
        AN__r_half_vec_diffcubed_inv = AN.r_half_vec_diffcubed_inv;

        CA__C_Li_max = CA.C_Li_max;
        CA__C_Li_max_inv = CA.C_Li_max_inv;
        CA__A_surf_CV_inv = CA.A_surf_CV_inv;
        CA__R_SEI = CA.R_SEI;
        CA__R_SEI_inv = CA.R_SEI_inv;
        CA__del_x_inv = CA.del_x_inv;
        CA__A_s = CA.A_s;
        CA__r_half_vec_sq = CA.r_half_vec_sq;
        CA__A_c = CA.A_c;
        CA__alpha_a = CA.alpha_a;
        CA__alpha_c = CA.alpha_c;
        CA__i_0_ref = CA.i_0_ref;
        CA__r_half_vec_diffcubed_inv = CA.r_half_vec_diffcubed_inv;

        EL__C_inv = EL.C_inv;
        SEP__del_x_inv = SEP.del_x_inv;
        SIM__SimMode = SIM.SimMode;
        SIM__A_c = SIM.A_c;
        SIM__A_outside_vec = SIM.A_outside_vec;
        if isfield(SIM , 'Amp')
            SIM__Amp = SIM.Amp;
        end
        if isfield(SIM , 'Temp_AN_BC')
            SIM__Temp_AN_BC = SIM.Temp_AN_BC;
        end
        SIM__h = SIM.h;
        SIM__T_inf = SIM.T_inf;
        if isfield(SIM , 'h_AN_BC')
            SIM__h_AN_BC = SIM.h_AN_BC;
        end
        if isfield(SIM , 'h_CA_BC')
            SIM__h_CA_BC = SIM.h_CA_BC;
        end
        SIM__q_CA_BC = SIM.q_CA_BC;
        SIM__q_AN_BC = SIM.q_AN_BC;
        SIM__profile_time = SIM.profile_time;
        SIM__profile_current = SIM.profile_current;
        P__T = P.T;
        P__del_phi   = P.del_phi;
        P__phi_ed    = P.phi_ed;    
        P__V_1       = P.V_1;
        P__V_2       = P.V_2;
        P__i_PS      = P.i_PS;
        P__C_Liion   = P.C_Liion;
        P__C_Li      = P.C_Li;
	    P__C_Li_surf_AN  = P.C_Li_surf_AN;
	    P__C_Li_surf_CA  = P.C_Li_surf_CA;
        P__SEP__T       = P.SEP.T;
        P__SEP__phi_el  = P.SEP.phi_el;
        P__SEP__C_Liion = P.SEP.C_Liion;
        P__phi_el = P.phi_el;
        P__sigma = P.sigma;      
        P__kappa = P.kappa;     
        P__R_SEI = P.R_SEI;      
        P__lambda_eff = P.lambda_eff; 
        P__rho_eff = P.rho_eff;    
        P__c_p_eff = P.c_p_eff;    
        P__D_o_Li_ion = P.D_o_Li_ion; 
        P__activity = P.activity;   
        P__tf_num = P.tf_num;     
        P__D_o = P.D_o;   
        N__N_CV_tot = N.N_CV_tot;
        N__CV_Region_AN = N.CV_Region_AN;
        N__N_SV_AN = N.N_SV_AN;
        N__N_CV_AN = N.N_CV_AN;
        N__N_R_AN = N.N_R_AN;
        N__CV_Region_CA = N.CV_Region_CA;
        N__CV_Region_SEP = N.CV_Region_SEP;
        N__N_SV_CA = N.N_SV_CA;
        N__N_CV_CA = N.N_CV_CA;
        N__N_R_CA = N.N_R_CA;
        N__N_CV_SEP = N.N_CV_SEP;
        N__N_SV_SEP = N.N_SV_SEP;
        N__N_R_max = N.N_R_max;
        N__N_prop = N.N_prop;
        N__N_SV_AN_tot = N.N_SV_AN_tot;
        N__N_SV_SEP_tot = N.N_SV_SEP_tot;
        CONS__F = CONS.F;
        CONS__R = CONS.R;
        CONS__F_inv = CONS.F_inv;
        CONS__BRUG = CONS.BRUG;

        SIM__RampThermalGradientTime = SIM.RampThermalGradientTime;
        if FLAG__RampThermalGradient
            AN__m_thermal = AN.m_thermal;
            AN__b_thermal = AN.b_thermal;
            CA__m_thermal = CA.m_thermal;
            CA__b_thermal = CA.b_thermal;
        end

        SIM__diff_CV_x_vec_inv = SIM.diff_CV_x_vec_inv;
        SIM__interp_x_interface = SIM.interp_x_interface;
        SIM__MassFluxPreCalcResistor = SIM.MassFluxPreCalcResistor;
        SIM__del_x_vec_halved_inv = SIM.del_x_vec_halved_inv;
        SIM__del_CV_r_inv_mat = SIM.del_CV_r_inv_mat;

        % SIM__eps_el_vec = SIM.eps_el_vec;
        % AN__del_r_vec  = AN.del_r_vec;
        % CA__del_r_vec  = CA.del_r_vec;
        % AN__r_vec = AN.r_vec;
        % AN__eps_el = AN.eps_el;
        % AN__del_CV_r_vec = AN.del_CV_r_vec;
        % CA__r_vec = CA.r_vec;
        % CA__eps_el = CA.eps_el;
        % CA__del_CV_r_vec = CA.del_CV_r_vec;
        % SEP__eps_el = SEP.eps_el;
        % N__N_SV_nR = N.N_SV_nR;

        % AN__A_surf_CV = AN.A_surf_CV;
        % AN__del_x = AN.del_x;
        % CA__A_surf_CV = CA.A_surf_CV;
        % CA__del_x = CA.del_x;
        % EL__C = EL.C;
        % SEP__del_x = SEP.del_x;
        % CONS__R_inv = CONS.R_inv;
        % SIM__diff_CV_x_vec = SIM.diff_CV_x_vec;
        % SIM__del_x_vec_halved = SIM.del_x_vec_halved;
        % AN__r_half_vec = AN.r_half_vec;
        % CA__r_half_vec = CA.r_half_vec;
    end


end