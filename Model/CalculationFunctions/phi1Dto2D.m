%% Reshape phi: 1D to 2D
% Reorganize the phi values
% Insert NaN back into the phi matrix for i_ed in the seperator region
%%
function phi_out = phi1Dto2D( phi_in , N , P , FLAG)
phi_temp = NaN(N.N_ES_var, N.N_CV_tot);

% ---- Anode ----
    phi_an = reshape(phi_in(1:N.N_CV_AN*N.N_ES_var) , N.N_ES_var , N.N_CV_AN); 
    phi_temp( : , N.CV_Region_AN) = phi_an;
    
% ---- Separator ----
    phi_sep = reshape(phi_in(N.N_CV_AN*N.N_ES_var+1 : N.N_CV_AN*N.N_ES_var+N.N_CV_SEP) , 1 , N.N_CV_SEP);
    phi_temp( P.ES.phi_el , N.CV_Region_SEP ) = phi_sep;

% ---- Cathode ----
    phi_ca = reshape(phi_in(N.N_CV_AN*N.N_ES_var+N.N_CV_SEP+1:end) , N.N_ES_var , N.N_CV_CA);
    phi_temp( : , N.CV_Region_CA) = phi_ca;

% Set Output
phi_out = phi_temp;

end