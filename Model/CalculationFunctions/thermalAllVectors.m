%% thermalAllVectors
%
%

function [q_cond , q_conv , q_gen] = thermalAllVectors( SV , AN , SEP, CA , EL , P , N , CONS , FLAG , props, i_el, i_ed, i_Far)
%% Conduction 
    [q_cond] = thermalFluxConduction( SV , AN , SEP, CA , EL , P , N , CONS , FLAG , props); % [J s^-1 m^-2]


%% Convection
    [q_conv] = thermalFluxConvection( SV , SIM , P ); % [J s^-1 m_outside^-2]


%% Heat Gen
    if FLAG.HEAT_GEN_TOTAL 
        [q_gen] = thermalGenerationTotal( SV , AN , SEP, CA , EL , P , N , CONS , FLAG , props, i_el, i_ed, i_Far); % [J s^-1 ]
    else
        q_gen = zeros( 1 , N.N_CV_tot );
    end

end