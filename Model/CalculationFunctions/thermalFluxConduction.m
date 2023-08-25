%% thermalFluxConduction
%
%

%%%!!!!!!!! Not capable of running a Li foil
%%%!!!!!!!! Not capable of running temperature dependent thermal properties

function [q_cond] = thermalFluxConduction( SV , AN , SEP, CA , EL , P , N , CONS , FLAG , props)
%% Initialize    
    q_cond = nan(1,N.N_CV_tot+1); % ith term is the left face of the ith CV
    T = SV( P.T , :);
    lambda_vec = props(P.lambda_eff,:);
    

%% Anode Boundary Condition
    i = N.CV_Region_AN(1);
    switch FLAG.T_BC_AN
        case 1 % Known Temperature
            q_cond(i) = 0;
        case 2 % Known heat flux
            q_cond(i) = SIM.q_AN_BC; % [W m^-2]
        case 3 % Insulated
            q_cond(i) = SIM.q_AN_BC; % [W m^-2]
        case 4 % Convection
            q_cond(i) = SIM.h_AN_BC * ( T(i) - SIM.T_inf ); % [W m^-2]
    end


%% Anode Internal
    for i = N.CV_Region_AN(2:end)
        q_cond(i) = ( T(i) - T(i-1) ) * ( -lambda_vec(i) / AN.del_x );
    end


%% Anode/Separator BC
    i = N.CV_Region_SEP(1);
    R_AN  = ( -lambda_vec(i-1) / (AN.del_x /2) )^-1 ;
    R_SEP = ( -lambda_vec(i)   / (SEP.del_x/2) )^-1 ;
    R_tot = R_AN + R_SEP;
    q_cond(i) = ( T(i) - T(i-1) ) / R_tot;

    
%% Separator Internal
    for i = N.CV_Region_SEP(2:end)
        q_cond(i) = ( T(i) - T(i-1) ) * ( -lambda_vec(i) / SEP.del_x );
    end


%% Separator/Cathode BC
    i = N.CV_Region_CA(1);
    R_SEP = ( -lambda_vec(i-1) / (SEP.del_x/2) )^-1 ;
    R_CA  = ( -lambda_vec(i)   / (CA.del_x /2) )^-1 ;
    R_tot = R_SEP + R_CA;
    q_cond(i) = ( T(i) - T(i-1) ) / R_tot;


%% Cathode Internal
    for i = N.CV_Region_CA(2:end)
        q_cond(i) = ( T(i) - T(i-1) ) * ( -lambda_vec(i) / CA.del_x );
    end


%% Cathode BC
    i = N.N_CV_tot+1;
    switch FLAG.T_BC_CA
        case 1 % Known Temperature
            q_cond(i) = 0;
        case 2 % Known heat flux
            q_cond(i) = SIM.q_CA_BC; % [W m^-2]
        case 3 % Insulated
            q_cond(i) = SIM.q_CA_BC; % [W m^-2]
        case 4 % Convection
            q_cond(i) = SIM.h_CA_BC * ( T(i-1) - SIM.T_inf ); % [W m^-2]
    end

end