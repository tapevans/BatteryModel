function [sigma_vec, kappa_vec, tf_vec, activity_vec, D_o_Li_ion_vec, lambda_vec, D_o_vec] = extractProps(P__sigma, P__kappa, P__tf_num, P__activity, P__D_o_Li_ion, P__lambda_eff, P__D_o, props)
    sigma_vec      = props( P__sigma      , : );
    kappa_vec      = props( P__kappa      , : );
    tf_vec         = props( P__tf_num     , : );
    activity_vec   = props( P__activity   , : );
    D_o_Li_ion_vec = props( P__D_o_Li_ion , : );
    lambda_vec     = props( P__lambda_eff , : );
    D_o_vec        = props( P__D_o:end    , : );
end