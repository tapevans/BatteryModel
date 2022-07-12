%% Reshape SV: 1D to 2D
function SV_out = SV1Dto2D(SV_in , N , P , FLAG)

r = max(N.N_SV_AN , N.N_SV_CA);
SV_temp = NaN(r , N.N_CV_tot);

% ---- Anode ----
    SV_an = reshape(SV_in(1:N.N_SV_AN_tot) , N.N_SV_AN , N.N_CV_AN);
    SV_temp(1:N.N_SV_AN , N.CV_Region_AN) = SV_an;

% ---- Separator ----
    SV_sep = reshape(SV_in(N.N_SV_AN_tot+1:N.N_SV_AN_tot+N.N_SV_SEP_tot) , N.N_SV_SEP , N.N_CV_SEP);
    SV_temp( P.T       , N.CV_Region_SEP ) = SV_sep(P.SEP.T       , :);
    SV_temp( P.phi_el  , N.CV_Region_SEP ) = SV_sep(P.SEP.phi_el  , :);
    SV_temp( P.C_Liion , N.CV_Region_SEP ) = SV_sep(P.SEP.C_Liion , :);

% ---- Cathode ----
    SV_ca = reshape(SV_in(N.N_SV_AN_tot+N.N_SV_SEP_tot+1:end) , N.N_SV_CA , N.N_CV_CA);
    SV_temp(1:N.N_SV_CA , N.CV_Region_CA) = SV_ca;

% Set output
SV_out = SV_temp;

end