function i_o = i_oLiFoil(SV , P  , ED)
RR  = 8314.472;

T       = SV(P.T , :);
T_inv   = T.^-1;
T_o     = 303.15;
T_o_inv = T_o^-1;

i_o = 0.27*exp(-30E6 * (T_inv - T_o_inv) / RR )       ...
         .* (SV(P.C_Liion,:)                ).^ED.alpha_a ...
         .* (SV(P.C_Li_surf,:)              ).^ED.alpha_c ...
         .* ((ED.C_Li_max-SV(P.C_Li_surf,:))).^ED.alpha_a;
end