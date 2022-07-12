function q_gen = calcHeatGenTot( SV , AN , SEP, CA , EL , P , N , CONS , SIM , FLAG , props , i_el, i_ed, i_Far)
%% Initialize
q_gen_io       = zeros( 1 , N.N_CV_tot );
q_gen_elec     = zeros( 1 , N.N_CV_tot );
q_gen_rxn      = zeros( 1 , N.N_CV_tot );
q_gen_chem_rxn = zeros( 1 , N.N_CV_tot );

%% CALCULATE HEAT GEN FROM IONIC FLUX
if FLAG.HEAT_GEN_IONIC 
    q_gen_io = ( i_el(1:end-1).^2 + i_el(2:end).^2 ) ...
                  ./ ( 2 * props(P.kappa,:) );
end

%% CALCULATE HEAT GEN FROM ELECTRONIC FLUX 
if FLAG.HEAT_GEN_ELECTRIC  %%%%%%%%%%Check the NaN on this one
    q_gen_elec = ( i_ed(1:end-1).^2 + i_ed(2:end).^2 ) ...
                  ./ ( 2 * props(P.sigma,:) );
    q_gen_elec(N.CV_Region_SEP) = zeros(1,N.N_CV_SEP);
end

%% CALCULATE HEAT GEN FROM CURRENT AT SURFACE RXN
if FLAG.HEAT_GEN_RXN 
    q_gen_rxn = i_Far.^2 .* SIM.A_surf_CV_vec.^2 .* props(P.R_SEI,:)...
                ./ SIM.CV_vec ;
end

%% CALCULATE HEAT GEN FROM CHEM RXN
if FLAG.HEAT_GEN_CHEM_RXN 
%     q_gen_chem_rxn %%%%%%%%%%%%%%%%%%%%%% NOT IMPLEMENTED
end

%% Calculate total volumetric heat generation from all sources
q_gen = q_gen_io + q_gen_elec + q_gen_rxn + q_gen_chem_rxn;


end