%% extractSV
function [T, del_phi, phi_ed, phi_el, V_1, V_2, i_PS, Ce, C_Li, C_Li_AN, C_Li_CA, X_AN, X_CA, Ce_norm, Ce_log, eta, RT_inv_vec] = extractSV(SV,P__T, P__del_phi, P__phi_ed, P__phi_el, P__V_1, P__V_2, P__i_PS, P__C_Liion, P__C_Li, P__C_Li_surf_AN, P__C_Li_surf_CA, N__CV_Region_AN, N__CV_Region_CA, N__N_R_max, AN__C_Li_max_inv, CA__C_Li_max_inv, EL__C_inv, CONS__R)
    [r,~] = size(SV);
    T       = SV(P__T       , :);      
    del_phi = SV(P__del_phi , :);
    phi_ed  = SV(P__phi_ed  , :);
    if r < P__phi_el
        phi_el = [];
    else
        phi_el  = SV(P__phi_el  , :);
    end
    V_1     = SV(P__V_1     , :);
    V_2     = SV(P__V_2     , :);  
    i_PS    = SV(P__i_PS    , :);
    Ce      = SV(P__C_Liion , :);
    C_Li    = SV(P__C_Li:P__C_Li+N__N_R_max-1 , : );
    C_Li_AN = SV(P__C_Li:P__C_Li_surf_AN , N__CV_Region_AN );
    C_Li_CA = SV(P__C_Li:P__C_Li_surf_CA , N__CV_Region_CA );
    
    X_AN    = C_Li_AN * AN__C_Li_max_inv;
    X_CA    = C_Li_CA * CA__C_Li_max_inv;

    Ce_norm = Ce * EL__C_inv;
    Ce_log  = log(Ce);

    eta   = V_2 - V_1;

    RT_inv_vec = (CONS__R * T).^(-1);

    % X_AN    = C_Li_AN / AN__C_Li_max;
    % X_CA    = C_Li_CA / CA__C_Li_max;
    % Ce_norm = Ce / EL__C;
end