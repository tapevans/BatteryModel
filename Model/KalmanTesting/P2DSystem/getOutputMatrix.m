%% Create Output Matrix
function [OutputMatrix,P] = getOutputMatrix(N,P,AN)
% Pointers for Output Matrix
    i = 1;
    P.OM.cell_volt = i; i = i + 1;
    P.OM.del_phi   = i; i = i + 1;
%     P.OM.temp      = i; i = i + 1;
    P.OM.C_Liion   = i; i = i + 1;
    P.OM.X_surf    = i; i = i + 1;
    P.OM.i_Far     = i; i = i + 1;
    P.OM.eta       = i; i = i + 1;
    
    N.N_Out = length(fieldnames(P.OM));

    N.N_In  = 1;
        % I_user
    % Outputs
        % Cell Voltage
        % Delta Phi   @AN/SEP
        % Temperature @AN/SEP
        % C_Liion     @AN/SEP
        % X_surf      @AN/SEP
        % i_Far       @AN/SEP
        % eta         @AN/SEP
    OutputMatrix = zeros(N.N_Out , N.N_SV_tot);
        % Cell Voltage
            idx_phi_ed_AN = P.phi_ed;

            i = N.N_CV_CA(end);
            index_offset = (i-1)*N.N_SV_CA + N.N_SV_AN_tot + N.N_SV_SEP_tot;
            idx_phi_ed_CA = index_offset + P.phi_ed;

            OutputMatrix(P.OM.cell_volt,idx_phi_ed_AN) = -1;
            OutputMatrix(P.OM.cell_volt,idx_phi_ed_CA) =  1;
        % @AN/SEP
            i = N.N_CV_AN(end);
            index_offset = (i-1)*N.N_SV_AN;
        % Delta Phi   @AN/SEP
            idx = index_offset + P.del_phi;
            OutputMatrix(P.OM.del_phi,idx) =  1;
%         % Temperature @AN/SEP
%             idx = index_offset + P.T;
%             OutputMatrix(P.OM.temp,idx) = 1;
        % C_Liion     @AN/SEP
            idx = index_offset + P.C_Liion;
            OutputMatrix(P.OM.C_Liion,idx) = 1;
        % X_surf      @AN/SEP
            idx = index_offset + P.C_Li_surf_AN;
            OutputMatrix(P.OM.X_surf,idx) = 1/AN.C_Li_max;
        % i_Far      @AN/SEP
            idx = index_offset + P.i_PS;
            OutputMatrix(P.OM.i_Far,idx) = 1;
        % Eta      @AN/SEP
            idx = index_offset + P.V_2;
            OutputMatrix(P.OM.eta,idx) = 1;
            idx = index_offset + P.V_1;
            OutputMatrix(P.OM.eta,idx) = -1;

end