%% Function to calculate i_Far
% The ith index is the faradaic current of the ith CV
function i_Far = iFarCalc( SV , AN , CA , P , N , CONS , FLAG , props)
%% Initialize
    i_Far = zeros(1 , N.N_CV_tot);

%% Concentration Dependent Parameters
% % Delta phi
%     del_phi = SV(P.phi_ed , :) - SV(P.phi_el , :);

% eta
    eta_AN   = SV(P.V_2 , N.CV_Region_AN ) - SV(P.V_1 , N.CV_Region_AN );
    eta_CA   = SV(P.V_2 , N.CV_Region_CA ) - SV(P.V_1 , N.CV_Region_CA );
    eta   = [eta_AN  , NaN(1,N.N_CV_SEP) ,eta_CA  ];

% i_o
    if FLAG.Newman_i_o
        i_o_an = CONS.F * AN.k_o ...
                           * SV(P.C_Liion     ,N.CV_Region_AN)  .^AN.alpha_a ...
          .* ( AN.C_Li_max - SV(P.C_Li_surf_AN,N.CV_Region_AN) ).^AN.alpha_a ...
                          .* SV(P.C_Li_surf_AN,N.CV_Region_AN)  .^AN.alpha_c;
        i_o_ca = CONS.F * CA.k_o ...
                           * SV(P.C_Liion     ,N.CV_Region_CA)  .^CA.alpha_a ...
          .* ( CA.C_Li_max - SV(P.C_Li_surf_CA,N.CV_Region_CA) ).^CA.alpha_a ...
                          .* SV(P.C_Li_surf_CA,N.CV_Region_CA)  .^CA.alpha_c;
    else
        i_o_an  = AN.i_oHandle( SV(:,N.CV_Region_AN) , P, AN );
        i_o_ca  = CA.i_oHandle( SV(:,N.CV_Region_CA) , P, CA );
    end
    if FLAG.AN_LI_FOIL
        i_o_an = 10; %%%%%%%%%% Hard coded; pulled from https://www.sciencedirect.com/science/article/pii/S0378775312005423
    end
    i_o     = [i_o_an, NaN(1,N.N_CV_SEP), i_o_ca];

%% i_Far Calc
% Anode
    for i = N.CV_Region_AN
        i_Far(i) = i_o(i) * ( exp( AN.alpha_a*CONS.F*eta(i)/( CONS.R*SV(P.T,i) ))... 
                            - exp(-AN.alpha_c*CONS.F*eta(i)/( CONS.R*SV(P.T,i) )));
    end
 
% Cathode
    for i = N.CV_Region_CA
        i_Far(i) = i_o(i) * ( exp( CA.alpha_a*CONS.F*eta(i)/( CONS.R*SV(P.T,i) ))... 
                            - exp(-CA.alpha_c*CONS.F*eta(i)/( CONS.R*SV(P.T,i) )));
    end
end