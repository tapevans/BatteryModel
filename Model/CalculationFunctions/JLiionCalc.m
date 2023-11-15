%% Function to calculate J_Liion
% The ith index is the flux at the left (minus) face of the ith control volume
function J_Liion = JLiionCalc( SV , AN , SEP, CA , EL , P , N , CONS , FLAG , i_el , props)
%% Initialize
    J_Liion = zeros(1 , N.N_CV_tot+1 );

%% Calculate Flux
% ---- Anode ----
    % CC/AN Boundary Condition
        i = N.CV_Region_AN(1);
        J_Liion(i) = 0; % No flux through current collector

    % AN Region
        % if ~FLAG.AN_LI_FOIL
            for i = N.CV_Region_AN(2:end)
                D_o_Li_ion = (props( P.D_o_Li_ion , i ) + props( P.D_o_Li_ion , i-1 ))/2;
                tf         = (props( P.tf_num     , i ) + props( P.tf_num     , i-1 ))/2;
                
                J_Liion(i) = - AN.eps_el * D_o_Li_ion * ( SV(P.C_Liion,i) - SV(P.C_Liion,i-1) ) / AN.del_x...
                             + i_el(i) * tf / CONS.F;
            end
        % end

% ---- Separator ----
    % AN/SEP interface
        i = N.CV_Region_SEP(1);
        D_o_Li_ion = (props( P.D_o_Li_ion , i ) + props( P.D_o_Li_ion , i-1 ))/2;
        tf         = (props( P.tf_num     , i ) + props( P.tf_num     , i-1 ))/2;
        
        R_an  = ( AN.del_x/2)/( AN.eps_el * D_o_Li_ion);
        R_sep = (SEP.del_x/2)/(SEP.eps_el * D_o_Li_ion);
        R_tot = R_an + R_sep;
        
        J_Liion(i) = -( SV(P.C_Liion,i) - SV(P.C_Liion,i-1) ) / ( R_tot )...
                     + i_el(i) * tf / CONS.F;
    
    % SEP Region
        for i = N.CV_Region_SEP(2:end)
            D_o_Li_ion = (props( P.D_o_Li_ion , i ) + props( P.D_o_Li_ion , i-1 ))/2;
            tf         = (props( P.tf_num     , i ) + props( P.tf_num     , i-1 ))/2;
            
            J_Liion(i) = - SEP.eps_el * D_o_Li_ion * ( SV(P.C_Liion,i) - SV(P.C_Liion,i-1) ) / SEP.del_x...
                         + i_el(i) * tf / CONS.F;
        end
 
% ---- Cathode ----
    % SEP/CA interface
        i = N.CV_Region_CA(1) ;
        
        D_o_Li_ion = (props( P.D_o_Li_ion , i ) + props( P.D_o_Li_ion , i-1 ))/2;
        tf         = (props( P.tf_num     , i ) + props( P.tf_num     , i-1 ))/2;
        
        R_ca  = ( CA.del_x/2)/( CA.eps_el * D_o_Li_ion);
        R_sep = (SEP.del_x/2)/(SEP.eps_el * D_o_Li_ion);
        R_tot = R_ca + R_sep;
        
        J_Liion(i) = -( SV(P.C_Liion,i) - SV(P.C_Liion,i-1) ) / ( R_tot )...
                     + i_el(i) * tf / CONS.F;
                
    % CA region
        % if ~FLAG.CA_LI_FOIL
            for i = N.CV_Region_CA(2:end)
                D_o_Li_ion = (props( P.D_o_Li_ion , i ) + props( P.D_o_Li_ion , i-1 ))/2;
                tf         = (props( P.tf_num     , i ) + props( P.tf_num     , i-1 ))/2;
                
                J_Liion(i) = - CA.eps_el * D_o_Li_ion * ( SV(P.C_Liion,i) - SV(P.C_Liion,i-1) ) / CA.del_x...
                             + i_el(i) * tf / CONS.F;
            end
        % end

% Boundary Condition at the CA/CC
    i = N.N_CV_tot + 1;
    J_Liion(i) = 0;    % No species flux through current collector

end        