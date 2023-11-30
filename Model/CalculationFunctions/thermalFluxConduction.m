%% thermalFluxConduction
%
%

%%%!!!!!!!! Not capable of running a Li foil

function [q_cond] = thermalFluxConduction( T, T_diff, lambda_vec , SIM__del_x_vec_halved_inv, N__CV_Region_AN, N__N_CV_tot, SIM__T_inf, SIM__q_AN_BC, SIM__q_CA_BC, SIM__h_AN_BC, SIM__h_CA_BC, FLAG__T_BC_AN, FLAG__T_BC_CA)
    % Resistor_vector = (lambda_vec ./ SIM__del_x_vec_halved).^-1;
    Resistor_vector = (lambda_vec .* SIM__del_x_vec_halved_inv).^-1;
    Resistor_vector_tot = [nan Resistor_vector(2:end) + Resistor_vector(1:end-1)  nan];
    q_cond = - Resistor_vector_tot.^-1  .* T_diff; %%%% !!!!!!!!!! NOT TESTED                        


%% Anode Boundary Condition
    i = N__CV_Region_AN(1);
    switch FLAG__T_BC_AN
        case 1 % Known Temperature
            q_cond(i) = 0;
        case 2 % Known heat flux
            q_cond(i) = SIM__q_AN_BC; % [W m^-2]
        case 3 % Insulated
            q_cond(i) = SIM__q_AN_BC; % [W m^-2]
        case 4 % Convection
            q_cond(i) = SIM__h_AN_BC * ( SIM__T_inf - T(i) ); % [W m^-2], 
            % Flux is written as positive in the positive x-direction
            % Therefore, positive convection at the anode occurs when T_inf
            % is greater than T_AN
    end                    


%% Cathode BC
    i = N__N_CV_tot+1;
    switch FLAG__T_BC_CA
        case 1 % Known Temperature
            q_cond(i) = 0;
        case 2 % Known heat flux
            q_cond(i) = SIM__q_CA_BC; % [W m^-2]
        case 3 % Insulated
            q_cond(i) = SIM__q_CA_BC; % [W m^-2]
        case 4 % Convection
            q_cond(i) = SIM__h_CA_BC * ( T(i-1) - SIM__T_inf ); % [W m^-2], 
            % Flux is written as positive in the positive x-direction
            % Therefore, positive convection at the cathode occurs when T_CA
            % is greater than T_inf
    end                    


%% Fix Known Temp BC
    if FLAG__T_BC_AN == 1
        q_cond(1) = q_cond(2);
    end
    if FLAG__T_BC_CA == 1
        q_cond(end) = q_cond(end-1);
    end
end