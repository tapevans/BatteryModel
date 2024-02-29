%% Reshape SV: 1D to 2D
function SV_out = SV1Dto2D(SV_in , N__N_SV_max, N__N_CV_tot, N__N_SV_AN_tot, N__N_SV_SEP_tot, N__N_SV_AN, N__N_SV_SEP, N__N_SV_CA, N__N_CV_AN, N__N_CV_SEP, N__N_CV_CA, N__CV_Region_AN, N__CV_Region_SEP, N__CV_Region_CA, P__T, P__del_phi, P__C_Liion, P__SEP__T, P__SEP__phi_el, P__SEP__C_Liion)
    SV_temp = NaN(N__N_SV_max , N__N_CV_tot);
    
    % ---- Anode ----
        SV_an = reshape(SV_in(1:N__N_SV_AN_tot) , N__N_SV_AN , N__N_CV_AN);
        SV_temp(1:N__N_SV_AN , N__CV_Region_AN) = SV_an;
    
    % ---- Separator ----
        SV_sep = reshape(SV_in(N__N_SV_AN_tot+1:N__N_SV_AN_tot+N__N_SV_SEP_tot) , N__N_SV_SEP , N__N_CV_SEP);
        SV_temp( P__T       , N__CV_Region_SEP ) = SV_sep(P__SEP__T       , :);
        SV_temp( P__del_phi , N__CV_Region_SEP ) = SV_sep(P__SEP__phi_el  , :);
        SV_temp( P__C_Liion , N__CV_Region_SEP ) = SV_sep(P__SEP__C_Liion , :);
    
    % ---- Cathode ----
        SV_ca = reshape(SV_in(N__N_SV_AN_tot+N__N_SV_SEP_tot+1:end) , N__N_SV_CA , N__N_CV_CA);
        SV_temp(1:N__N_SV_CA , N__CV_Region_CA) = SV_ca;

    % Set output
        SV_out = SV_temp;
end