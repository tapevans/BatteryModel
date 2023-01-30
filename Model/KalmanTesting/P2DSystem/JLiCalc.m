%% Function to calculate J_Li
% The ith index is for the ith CV
% The jth index is for the flux at the left (minus) face of the jth radial control volume
% j = 1 is r = 0

function J_Li = JLiCalc( SV , AN , CA , P , N , s_dot , props)
%% Initialize
J_Li = NaN(N.N_R_max+1 , N.N_CV_tot); % [kmol m^-2 s^-1], Radial molar flux (in electrode)(m^-2 is the cross-sectional area)

%% Calc
% Symmetry BC
    j = 1;
    J_Li(j,:) = zeros( 1 , N.N_CV_tot );

% Internal Region
    % ---- Anode ----
    for i = N.CV_Region_AN
        % Internal CV
        for j = 2:N.N_R_AN
            D_o = (props( P.D_o+(j-1)-1 , i ) + props( P.D_o+(j-1)-1 , i ))/2;
            
            J_Li(j,i) = -D_o * ( SV(N.N_SV_nR+j , i) - SV(N.N_SV_nR+j-1 , i) ) ...
                / (AN.r_vec(j) - AN.r_vec(j-1));
        end
    end

    % ---- Cathode ----
    for i = N.CV_Region_CA
        % Internal CV
        for j = 2:N.N_R_CA
            D_o = (props( P.D_o+(j-1)-1 , i ) + props( P.D_o+(j-1)-1 , i ))/2;
            
            J_Li(j,i) = -D_o * ( SV(N.N_SV_nR+j , i) - SV(N.N_SV_nR+j-1 , i) ) ...
                / (CA.r_vec(j) - CA.r_vec(j-1));
        end
    end

% Surface Reaction BC
    j = N.N_R_AN + 1;
    J_Li(j,N.CV_Region_AN) = s_dot(1,N.CV_Region_AN);

    j = N.N_R_CA + 1;
    J_Li(j,N.CV_Region_CA) = s_dot(1,N.CV_Region_CA);

end