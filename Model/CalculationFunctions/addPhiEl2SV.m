% Add phi_el to the SV matrix and put NaN in del_phi for the sep region
function [SV_out] = addPhiEl2SV(SV_in, P__phi_ed, P__del_phi, N__CV_Region_SEP, N__N_CV_SEP)
    phi_el = SV_in(P__phi_ed,:) - SV_in(P__del_phi,:);
    phi_el(1,N__CV_Region_SEP) = SV_in(P__del_phi,N__CV_Region_SEP);
    SV_out = [SV_in;phi_el];
    SV_out(P__del_phi,N__CV_Region_SEP) = NaN(1,N__N_CV_SEP);
end