%% Reaction Rate (i_o) NMC
% The input parameters are assumed to have the following units:
    %   - X_Li:         Intercalation fraction at the surface           [-]
    %   - C_Liion:      Li-ion concentration at the surface      [kmol/m^3]    
    %   - T:            Temperature                                     [K]
function i_o_NMC = i_oNMC_HZ(SV , P  , ED, EL)
% i_star_T
    RR      = 8314.472; % [J kmol^-1 K^-1], Gas Constant
    E_Li    = 30e6; % [J kmol^-1], activation energy
    T       = SV(P.T , :);
    T_inv   = T.^-1;
    T_o     = 303.15;
    T_o_inv = T_o^-1;
    i_star_T = ED.i_0_ref * exp(-E_Li * (T_inv - T_o_inv) / RR );

% i_star_Liion
    Ce = SV(P.C_Liion,:);
    Ce_norm = Ce/EL.C;
    i_star_Liion = Ce_norm.^ED.alpha_a;

% i_star_Li
    a_k = [
        -  3.585290065824760
        +  32.49768821737960
        -  94.16571081287610
        +  124.0524690073040
        -  75.23567141488800
        +  16.50452829641290];

    C_Li = SV(P.C_Li_surf_AN,:);
    C_Li_norm = C_Li/ED.C_Li_max;
    
    i_star_Li = zeros( 1 , length(C_Li_norm) );
    for i = 1:length(a_k)
        i_star_Li = i_star_Li + a_k(i) * (C_Li_norm.^(i-1));
    end
    i_star_Li = (1 - C_Li_norm).^ED.alpha_a .* C_Li_norm.^ED.alpha_c;
    
% i_o
    i_o_NMC = i_star_T .* i_star_Liion .* i_star_Li;



    % RR      = 8314.472; % Universal gas constant [J/kmol/K]
    % 
    % T       = SV(P.T , :);
    % T_inv   = T.^-1;
    % T_o     = 303.15;
    % T_o_inv = T_o^-1;
    % 
    % C_Liion = SV( P.C_Liion , :);
    % X_Li    = SV( P.C_Li_surf_CA , : ) / ED.C_Li_max; %%%!!P.C_Li_surf_CA is hardcoded for cathode
    
    % i_o_NMC = 3*...
    %     (  16.50452829641290*X_Li.^5 ...
    %     -  75.23567141488800*X_Li.^4 ...
    %     +  124.0524690073040*X_Li.^3 ...
    %     -  94.16571081287610*X_Li.^2 ...
    %     +  32.49768821737960*X_Li    ...
    %     -  3.585290065824760        )...
    %     .* (C_Liion./1.2).^0.5 ... 
    %     .*exp(-30E6 * (T_inv - T_o_inv) / RR );
end

% This came from Andrew C. modeling code
% i_o_NMC = (...
%    16.50452829641290*X_Li.^5 ...
% -  75.23567141488800*X_Li.^4 ...
% +  124.0524690073040*X_Li.^3 ...
% -  94.16571081287610*X_Li.^2 ...
% +  32.49768821737960*X_Li    ...
% -  3.585290065824760)...
% .* (C_Liion./SIM.C_Liion_init).^0.5 ... % Dividing by 1.2 because the initial concentration is 1.2M (Changed to 1200 to convert between mol and kmol)
% *exp(-30E6 * (1/T - 1/T_o) / RR )*10.0;