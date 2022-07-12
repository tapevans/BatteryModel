function Res = phiFsolveFun(phi,SV,AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,i_user,props)
%% Initialize 
phi = phi1Dto2D( phi , N , P , FLAG);
SV(P.phi_el:P.i_PS , :) = phi(P.ES.phi_el:P.ES.i_PS , :);

%% Calc Currents and Voltages
[i_ed , i_el ] = currentCalc( SV , AN , SEP , CA , EL , P , N , CONS , FLAG , i_user , props);
i_Far          = iFarCalc( SV , AN , CA , P , N , CONS, FLAG , props);

E_eq_an  = AN.EqPotentialHandle( SV(P.C_Li_surf_AN , N.CV_Region_AN ) / AN.C_Li_max );
E_eq_ca  = CA.EqPotentialHandle( SV(P.C_Li_surf_CA , N.CV_Region_CA ) / CA.C_Li_max );
E_eq_vec = [ E_eq_an , zeros(1,N.N_CV_SEP) , E_eq_ca];

%% Calculate Residual
% ---- Anode ----
Res_AN = zeros(N.N_ES_var, N.N_CV_AN);
for i = N.CV_Region_AN
    Res_AN(P.ES.phi_el, i) = (AN.A_c / AN.A_surf)*(i_el(i+1) - i_el(i)  ) ...
                            - (SV(P.V_1,i) - SV(P.phi_el,i))/AN.R_SEI...
                            -  phi(P.ES.i_dl,i);
    Res_AN(P.ES.phi_ed, i) = (AN.A_c / AN.A_surf)*(i_ed(i)   - i_ed(i+1))...
                            -  SV(P.i_PS,i)...
                            -  phi(P.ES.i_dl,i);
    Res_AN(P.ES.V_1   , i) = (SV(P.V_1,i) - SV(P.phi_el,i))/AN.R_SEI ...
                            -  i_Far(i);
    Res_AN(P.ES.V_2   , i) = i_Far(i)...
                            -  SV(P.i_PS,i);
    Res_AN(P.ES.i_PS  , i) = SV(P.phi_ed,i) - SV(P.V_2,i) ...
                            -  E_eq_vec(i);
    Res_AN(P.ES.i_dl  , i) = SV(P.phi_ed,i) - SV(P.phi_el,i) ...
                            -  E_eq_vec(i);
end 
    % Fix BC 
    Res_AN(P.ES.phi_ed,1) = SV(P.phi_ed,1) - 0;

% ---- Separator ----
Res_SEP = zeros(1, N.N_CV_SEP);
for i = 1:N.N_CV_SEP
    CV_offset = N.N_CV_AN;
    Res_SEP(P.ES.phi_el, i) = i_el(i+1+CV_offset) - i_el(i+CV_offset);
end

% ---- Cathode ----
Res_CA = zeros(N.N_ES_var, N.N_CV_CA);
for i = 1:N.N_CV_CA
    CV_offset = N.N_CV_AN + N.N_CV_SEP + i;
    Res_CA(P.ES.phi_el, i) = (CA.A_c / CA.A_surf)*(i_el(CV_offset+1) - i_el(CV_offset)  ) ...
                            - (SV(P.V_1, CV_offset) - SV(P.phi_el,CV_offset))/CA.R_SEI...
                            -  phi(P.ES.i_dl,CV_offset);
    Res_CA(P.ES.phi_ed, i) = (CA.A_c / CA.A_surf)*(i_ed(CV_offset)   - i_ed(CV_offset+1))...
                            -  SV(P.i_PS,CV_offset)...
                            -  phi(P.ES.i_dl,CV_offset);
    Res_CA(P.ES.V_1   , i) = (SV(P.V_1,CV_offset) - SV(P.phi_el,CV_offset))/CA.R_SEI ...
                            -  i_Far(CV_offset);
    Res_CA(P.ES.V_2   , i) = i_Far(CV_offset)...
                            -  SV(P.i_PS,CV_offset);
    Res_CA(P.ES.i_PS  , i) = SV(P.phi_ed,CV_offset) - SV(P.V_2,CV_offset) ...
                            -  E_eq_vec(CV_offset);
    Res_CA(P.ES.i_dl  , i) = SV(P.phi_ed,CV_offset) - SV(P.phi_el,CV_offset) ...
                            -  E_eq_vec(CV_offset);
end

% Reshape
Res_AN  = reshape(Res_AN ,[],1);
Res_SEP = reshape(Res_SEP,[],1);
Res_CA  = reshape(Res_CA ,[],1);

Res = [Res_AN ; Res_SEP ; Res_CA];

end