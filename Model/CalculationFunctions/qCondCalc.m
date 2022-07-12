%% Heat Flux via Conduction
% This function is used to calculate the heat flux due to conduction.
% Index i refers to the flux at the left face of ith control volume 
function q_cond = qCondCalc( SV , AN , SEP, CA , EL , P , N , SIM , CONS , FLAG , props)
%% Initialize
% q_cond = zeros(1 , N.N_CV_tot + 1);

% Half resistance of each cell (It is the same for each cell. See Write-up)
R_half = [AN.del_x * ones( 1 , N.N_CV_AN ) , SEP.del_x * ones( 1 , N.N_CV_SEP ) , CA.del_x * ones( 1 , N.N_CV_CA )]...
         ./ (2 * props( P.k , : ) );

R_tot = [ 1 , R_half(2:end)+R_half(1:end-1) , 1 ];

deltaT = [ 0 , ( SV(P.T,2:end)-SV(P.T,1:end-1) ) , 0 ];

%% Calculate Flux
% Interior Nodes
    q_cond = -deltaT./R_tot;

% CC/AN Boundary Condition
    i = N.CV_Region_AN(1);
    switch(FLAG.T_BC_AN)
        case 1 % Known Temp
            %?
        case 2 % Known Heat Flux
            q_cond(i) = SIM.q_AN_BC;
        case 3 % Insulated
            q_cond(i) = 0; % No flux through current collector
        case 4 % Convection
            q_cond(i) = SIM.h * ( SIM.T_inf - SV(P.T , i) );
        otherwise
            disp('Not a valid AN temp BC')                
    end

% Boundary Condition at the CA/CC
    i = N.N_CV_tot + 1;
    switch(FLAG.T_BC_CA)
        case 1 % Known Temp
            %?
        case 2 % Known Heat Flux
            q_cond(i) = SIM.q_CA_BC;
        case 3 % Insulated
            q_cond(i) = 0; % No flux through current collector
        case 4 % Convection
            q_cond(i) = SIM.h * ( SV(P.T , i-1) - SIM.T_inf );
        otherwise
            disp('Not a valid CA temp BC')
    end
end