%% kappa

function [kappa_out] = kappa_HZ(Ce , T)
    Acoeff_Kappa =                                                         0.0001909446*(T.^2.0) - 0.08038545*T + 9.003410;
    Bcoeff_Kappa = -0.000000028875870*(T.^4.0) + 0.000034836380*(T.^3.0) - 0.0158367700*(T.^2.0) + 3.19529500*T - 241.4638;
    Ccoeff_Kappa =  0.000000016537860*(T.^4.0) - 0.000019987600*(T.^3.0) + 0.0090711550*(T.^2.0) - 1.82806400*T + 138.0976;
    Dcoeff_Kappa = -0.000000002791965*(T.^4.0) + 0.000003377143*(T.^3.0) - 0.0015327070*(T.^2.0) + 0.30900030*T - 23.35671;
    
    kappa_out = Ce .* (Acoeff_Kappa + Bcoeff_Kappa .* Ce + Ccoeff_Kappa .* (Ce.^2.0) + Dcoeff_Kappa .* (Ce.^3.0));

    kappa_out = kappa_out*1.0;
end