%% thermalGenerationTotal
%
%

function [q_gen] = thermalGenerationTotal( SV , AN , SEP, CA , EL , SIM , P , N , CONS , FLAG , props, i_el, i_ed, i_Far)
%% Initialize
    q_gen_io_minusface   = zeros( 1 , N.N_CV_tot );
    q_gen_io_plusface    = zeros( 1 , N.N_CV_tot );
    q_gen_elec_minusface = zeros( 1 , N.N_CV_tot );
    q_gen_elec_plusface  = zeros( 1 , N.N_CV_tot );
    q_gen_rxn            = zeros( 1 , N.N_CV_tot );
    q_gen_chem_rxn       = zeros( 1 , N.N_CV_tot );

    V_SEI = SV( P.V_1 , :) - SV( P.phi_el , :);

    % Phi at the cell face
        del_x_vec = [AN.del_x * ones(1,N.N_CV_AN) , SEP.del_x * ones(1,N.N_CV_SEP), CA.del_x * ones(1,N.N_CV_CA)];
        phi_ed = SV( P.phi_ed , : );
        phi_el = SV( P.phi_el , : );

        phi_ed_half = nan( 1 , N.N_CV_tot+1 );
        phi_el_half = nan( 1 , N.N_CV_tot+1 );

        % Internal Faces
        for i = 2:N.N_CV_tot
            phi_ed_half(i) = ((phi_ed(i)-phi_ed(i-1)) / (SIM.x_vec(i) - SIM.x_vec(i-1))) * (SIM.x_half_vec(i) - SIM.x_vec(i-1)) + phi_ed(i-1);
            phi_el_half(i) = ((phi_el(i)-phi_el(i-1)) / (SIM.x_vec(i) - SIM.x_vec(i-1))) * (SIM.x_half_vec(i) - SIM.x_vec(i-1)) + phi_el(i-1);
        end

        % AN/CC
            i = 1;
            phi_ed_half(i) = ((phi_ed(i)-phi_ed_half(i+1)) / (SIM.x_vec(i) - SIM.x_half_vec(i+1))) * (SIM.x_half_vec(i) - SIM.x_vec(i)) + phi_ed(i);
            phi_el_half(i) = ((phi_el(i)-phi_el_half(i+1)) / (SIM.x_vec(i) - SIM.x_half_vec(i+1))) * (SIM.x_half_vec(i) - SIM.x_vec(i)) + phi_el(i);

        % CA/CC
            i = N.N_CV_tot + 1;
            phi_ed_half(i) = ((phi_ed(i-1)-phi_ed_half(i-1)) / (SIM.x_vec(i-1) - SIM.x_half_vec(i-1))) * (SIM.x_half_vec(i) - SIM.x_vec(i-1)) + phi_ed(i-1);
            phi_el_half(i) = ((phi_el(i-1)-phi_el_half(i-1)) / (SIM.x_vec(i-1) - SIM.x_half_vec(i-1))) * (SIM.x_half_vec(i) - SIM.x_vec(i-1)) + phi_el(i-1);

        % % Test plot
        %     figure(1)
        %     hold on
        %     plot(SIM.x_vec      , phi_ed      , '-k')
        %     plot(SIM.x_half_vec , phi_ed_half , 'ob')
        % 
        %     figure(2)
        %     hold on
        %     plot(SIM.x_vec      , phi_el      , '-k','Linewidth',2)
        %     plot(SIM.x_half_vec , phi_el_half , 'ob')


%% Calculate Heat Generation from Ionic Flux
    if FLAG.HEAT_GEN_IONIC 
        for i = 1:N.N_CV_tot
            q_gen_io_minusface(i) = abs( i_el(i)   * SIM.A_c * ( phi_el_half(i)   - phi_el(i) ) );
            q_gen_io_plusface(i)  = abs( i_el(i+1) * SIM.A_c * ( phi_el_half(i+1) - phi_el(i) ) );
        end
    end
    q_gen_io = q_gen_io_minusface + q_gen_io_plusface;


%% Calculate Heat Generation from Electronic Flux 
    if FLAG.HEAT_GEN_ELECTRIC  
    % Anode Region
        for i = N.CV_Region_AN(1:end-1)
            q_gen_elec_minusface(i) = abs( i_ed(i)   * SIM.A_c * ( phi_ed_half(i)   - phi_ed(i) ) );
            q_gen_elec_plusface(i)  = abs( i_ed(i+1) * SIM.A_c * ( phi_ed_half(i+1) - phi_ed(i) ) );
        end

        % Anode/Seperator Interface
            i = N.CV_Region_AN(end);
            q_gen_elec_minusface(i) = abs( i_ed(i)   * SIM.A_c * ( phi_ed_half(i)   - phi_ed(i) ) );

    % Cathode Region
        % Cathode/Seperator Interface
        i = N.CV_Region_CA(1);
        q_gen_elec_plusface(i)  = abs( i_ed(i+1) * SIM.A_c * ( phi_ed_half(i+1) - phi_ed(i) ) );

        for i = N.CV_Region_CA(2:end)
            q_gen_elec_minusface(i) = abs( i_ed(i)   * SIM.A_c * ( phi_ed_half(i)   - phi_ed(i) ) );
            q_gen_elec_plusface(i)  = abs( i_ed(i+1) * SIM.A_c * ( phi_ed_half(i+1) - phi_ed(i) ) );
        end
    end
    q_gen_elec = q_gen_elec_minusface + q_gen_elec_plusface;
    

%% Calculate Heat Generation from Charge Transfer Reaction
    if FLAG.HEAT_GEN_RXN 
        q_gen_rxn = abs( i_Far .* SIM.A_surf_CV_vec .* V_SEI );
        q_gen_rxn(N.CV_Region_SEP) = zeros(1,N.N_CV_SEP);
    end


%% Calculate Heat Generation from Chemical Reaction
% Heat release of reaction (endo/exo process)
    if FLAG.HEAT_GEN_CHEM_RXN 
        % q_gen_chem_rxn %%%%%%%%%%%%%%%%%%%%%% NOT IMPLEMENTED
    end


%% Calculate total volumetric heat generation from all sources
    q_gen = q_gen_io + q_gen_elec + q_gen_rxn + q_gen_chem_rxn;


end