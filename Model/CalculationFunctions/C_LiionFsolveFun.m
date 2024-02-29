%% C_LiionFsolveFun
% The goal of this function is to solve for an initial electrolyte
% concentration profile that satisfies the balance between thermal,
% electrostatic, and molar gradients. This function will only be used when
% an initial thermal gradient is applied to the cell.
%
% What is being solved: The system is in "equilibrium" when the molar flux
% of the species comes to zero between all the control volumes. There are
% N_CV+1 fluxes (one at each face of each control volume), where N_CV is
% the number of control volumes used in the model. The boundary condition
% at the current collector specifies that the species flux is zero, leaving
% N_CV-1 fluxes that can vary. There are N_CV concentrations to be found.
% For all but one of the concentrations, the residual is set to be the
% condition that the flux leaving its upstream face has to be zero. For the
% last concentration, the residual will be set such that the sum of all of
% the Li^+ species in the electrolyte doesn't change.
%
% One thing to note, this fsolve is used to obtain a better initial
% condition for the elyte concentration profile. In reality, changes in the
% electrochemical potential of the species would cause echem reactions to
% take place at the active material surface, changing more than just the
% electrolyte concentration.


function Res = C_LiionFsolveFun(CLiion,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,i_user,props,m_0,Vol_el)
%% Initialize 
    CLiion = reshape(CLiion, 1,[]);
    SV(P.C_Liion, :) = CLiion;
    SV = addPhiEl2SV(SV, P.phi_ed, P.del_phi, N.CV_Region_SEP, N.N_CV_SEP);


%% Pull out each SV vector
    [T, T_inv_vec, del_phi, phi_ed, phi_el, V_1, V_2, i_PS, Ce, C_Li, C_Li_AN, C_Li_CA, X_AN, X_CA, Ce_norm, Ce_log, eta, RT_inv_vec , CTRG , N_Particles] = extractSV(SV,P.T, P.del_phi, P.phi_ed, P.phi_el, P.V_1, P.V_2, P.i_PS, P.CTRGrow, P.NPartic, P.C_Liion, P.C_Li, P.C_Li_surf_AN, P.C_Li_surf_CA, N.CV_Region_AN, N.CV_Region_CA, N.N_R_max, AN.C_Li_max_inv, CA.C_Li_max_inv, EL.C_inv, CONS.R);


%% Property Vector
    [sigma_vec, kappa_vec, tf_vec, activity_vec, D_o_Li_ion_vec, lambda_vec, D_o_vec] = extractProps(P.sigma, P.kappa, P.tf_num, P.activity, P.D_o_Li_ion, P.lambda_eff, P.D_o, props);


%% Calculate Value at Interface between CV
    sigma_vec_interface      = getInterfaceXDir(sigma_vec      , SIM.interp_x_interface);
    kappa_vec_interface      = getInterfaceXDir(kappa_vec      , SIM.interp_x_interface);
    tf_vec_interface         = getInterfaceXDir(tf_vec         , SIM.interp_x_interface);
    activity_vec_interface   = getInterfaceXDir(activity_vec   , SIM.interp_x_interface);
    D_o_Li_ion_vec_interface = getInterfaceXDir(D_o_Li_ion_vec , SIM.interp_x_interface);
    T_interface              = getInterfaceXDir(T              , SIM.interp_x_interface);
    T_interface_inv_vec      = getInterfaceXDir(T_inv_vec      , SIM.interp_x_interface);
    Ce_interface             = getInterfaceXDir(Ce             , SIM.interp_x_interface);

    D_o_vec_interface        = [NaN(1,N.N_CV_tot) ; 0.5*(D_o_vec(1:end-1,:) + D_o_vec(2:end,:)) ; NaN(1,N.N_CV_tot)]; % This assumes constant del_r in both particles
    dmudc_vec_interface      = dmudc_Latz( Ce_interface , T_interface, tf_vec_interface , activity_vec_interface);

    
%% Calculate Gradient
    [phi_ed_diff, phi_ed_grad] = diffAndGradXCalc(phi_ed, SIM.diff_CV_x_vec_inv);
    [phi_el_diff, phi_el_grad] = diffAndGradXCalc(phi_el, SIM.diff_CV_x_vec_inv);
    % [Ce_log_diff, Ce_log_grad] = diffAndGradXCalc(Ce_log, SIM.diff_CV_x_vec_inv);
    [Ce_diff    , Ce_grad    ] = diffAndGradXCalc(Ce    , SIM.diff_CV_x_vec_inv);
    [T_diff     , T_grad     ] = diffAndGradXCalc(T     , SIM.diff_CV_x_vec_inv);

    diff        = C_Li(2:end,:) - C_Li(1:end-1,:);
    C_Li_diff   = [NaN(1 , N.N_CV_tot) ; diff ; NaN(1 , N.N_CV_tot)];


%% Calc Currents and Voltages
    [~ , i_el ] = currentCalc(  sigma_vec_interface, kappa_vec_interface, dmudc_vec_interface, tf_vec_interface, ...
                                phi_ed_grad, phi_el_grad, Ce_grad,T_grad                                       , ...
                                CONS.F_inv, EL.Beta                                                            , ...
                                N.CV_Region_AN, N.CV_Region_SEP, N.CV_Region_CA, N.N_CV_tot, i_user                 );
    
    J_Liion = JLiionCalc(  D_o_Li_ion_vec_interface , tf_vec_interface, Ce_interface , ...
                           Ce_grad , T_grad , i_el , T_interface_inv_vec             , ...
                           N.CV_Region_AN, N.N_CV_tot, CONS.F_inv , EL.S_T                );


%% Calculate Residual
    m_calc = CLiion*Vol_el';
    % Res    = [J_Liion(2:end-1)*1e8 , (m_calc - m_0)]'; % Using the upstream flux
    Res    = [(m_calc - m_0) , J_Liion(2:end-1)*1e8]'; % Using the downstream flux

end

    % [SV_out] = addPhiEl2SV(SV_in, P.phi_ed, P.del_phi, N.CV_Region_SEP, N.N_CV_SEP);
    % SV = addPhiEl2SV(SV,P,N);

    % i_Far = iFarCalc( T, Ce_norm, X_AN, X_CA, eta, RT_inv_vec, AN.i_0_ref , AN.alpha_a, AN.alpha_c, CA.i_0_ref, CA.alpha_a, CA.alpha_c, ...
    %                     FLAG.Newman_i_o, CONS.F, ...
    %                     N.CV_Region_AN, N.CV_Region_CA, N.N_CV_SEP, AN.i_oHandle , CA.i_oHandle);
    % 
    % E_eq_an  = AN.EqPotentialHandle( SV(P.C_Li_surf_AN , N.CV_Region_AN ) / AN.C_Li_max );
    % E_eq_ca  = CA.EqPotentialHandle( SV(P.C_Li_surf_CA , N.CV_Region_CA ) / CA.C_Li_max );
    % E_eq_vec = [ E_eq_an , zeros(1,N.N_CV_SEP) , E_eq_ca];

%% Calc Currents and Voltages
    % [~ , i_el ] = currentCalc( SV , AN , SEP , CA , EL , P , N , CONS , FLAG , i_user , props);
    % J_Liion      = JLiionCalc( SV , AN , SEP , CA , EL , P , N , CONS , FLAG , i_el , props);

    % [i_ed , i_el ] = currentCalc(  sigma_vec_interface, kappa_vec_interface, T_interface, activity_vec_interface, tf_vec_interface, ...
    % phi_ed_grad, phi_el_grad, Ce_log_grad, ...
    % CONS__F_inv, CONS__R, ...
    % N__CV_Region_AN, N__CV_Region_SEP, N__CV_Region_CA, N__N_CV_tot, i_user);

    % J_Liion = JLiionCalc( Ce_diff, D_o_Li_ion_vec, tf_vec_interface, N__CV_Region_AN, N__N_CV_tot, CONS__F_inv, SIM__MassFluxPreCalcResistor, i_el);