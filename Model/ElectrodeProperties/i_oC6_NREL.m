%% i_oC6
% The output exchange current density is in [A m^-2]

function i_o = i_oC6_NREL(T, Ce_norm, C_Li_norm , i_0_ref , alpha_a, alpha_c)
% i_star_T
    RR      = 8314.472; % [J kmol^-1 K^-1], Gas Constant
    E_Li    = 30e6;     % [J kmol^-1],      Activation Energy
    T_inv   = T.^-1;
    T_o     = 303.15;
    T_o_inv = T_o^-1;
    i_star_T = i_0_ref * exp(-E_Li * (T_inv - T_o_inv) / RR );

% i_star_Liion
    i_star_Liion = Ce_norm.^alpha_a;

% i_star_Li
    i_star_Li = (1 - C_Li_norm).^alpha_a .* C_Li_norm.^alpha_c;
    
    i_o = i_star_T .* i_star_Liion .* i_star_Li;
end




% RR  = 8314.472; % [J kmol^-1 K^-1], Gas Constant
% 
% T       = SV(P.T , :);
% T_inv   = T.^-1;
% T_o     = 303.15;
% T_o_inv = T_o^-1;
% 
% i_o = 1.3 *0.27*exp(-30E6 * (T_inv - T_o_inv) / RR )       ...
%          .* (SV(P.C_Liion,:)                   ).^ED.alpha_a ...
%          .* (SV(P.C_Li_surf_AN,:)              ).^ED.alpha_c ...
%          .* ((ED.C_Li_max-SV(P.C_Li_surf_AN,:))).^ED.alpha_a;
% 
% % Avoid Complex
%     if ~isreal(i_o)
%         for i = 1:length(i_o)
%             if ~isreal(i_o(i))
%                 i_o(i) = 0;
%             end
%         end
%     end

% My Old Function
% i_o = 0.27*exp(-30E6 * (T_inv - T_o_inv) / RR )       ...
%          .* (SV(P.C_Liion,:)                ).^ED.alpha_a ...
%          .* (SV(P.C_Li_surf,:)              ).^ED.alpha_c ...
%          .* ((ED.C_Li_max-SV(P.C_Li_surf,:))).^ED.alpha_a;
     
% Peter's Function
% i0 =  1.3.*0.27.*exp( -30E6 * (T_inv - T_o_inv) / RR ) ...
%                   .*Li_ion.^an.alpha_a ...
%                   .*(an.Li_max-x.*an.Li_max).^an.alpha_a ...
%                   .*(an.Li_max*x).^an.alpha_c;