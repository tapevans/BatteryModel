% Add phi_el to the SV matrix and put NaN in del_phi for the sep region
function [SV_out] = addPhiEl2SV(SV_in,P,N)
    phi_el = SV_in(P.phi_ed,:) - SV_in(P.del_phi,:);
    phi_el(1,N.CV_Region_SEP) = SV_in(P.del_phi,N.CV_Region_SEP);
    SV_out = [SV_in;phi_el];
    SV_out(P.del_phi,N.CV_Region_SEP) = NaN(1,N.N_CV_SEP);
end