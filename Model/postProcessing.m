%% Post-Processing Function
% Take all the data and convert into easy to read vectors
function postProcessing(filename)
%% Load file to workspace
    load(filename,'postProcessComplete')


%% do post-processing if it hasn't already
if ~postProcessComplete 
    load(filename)
%% Make Desired Times Vector
    if FLAG.ReduceSolnTime
        % I want the index of all times before 1 second, every 10 seconds after 1 second up until 500 sec before t_soln(end), All indicies of the last 500 seconds
        FirstTime   = 1;
        TimeFromEnd = 500;
        TimeDelta   = 10;

        t_mid_final = t_soln(end) - TimeFromEnd; % last time of the middle region
        t_vec_mid_des = FirstTime+TimeDelta : TimeDelta : t_mid_final;

        % Find Index for 1 sec
        [~,k_1] = min(abs(t_soln-FirstTime));

        % Find Mid indicies
        k_mid = zeros(1,length(t_vec_mid_des)); % indicies of the desired middle region
        for i = 1:length(t_vec_mid_des)
            [~,k_mid(i)] = min(abs(t_soln-t_vec_mid_des(i)));
        end

        % Assemble Index Vector
        idx = [ (1:1:k_1) , k_mid , (k_mid(end)+1:1:length(t_soln)) ];
        N_t_steps = length(idx);

        t_soln_OG = t_soln;
        t_soln = t_soln(idx);
    else
        N_t_steps = length(t_soln);
        idx = 1:1:N_t_steps;
    end


%% Initialize Variables
    % Ones, Zeros, and NaN vectors for each region
        ones_vec_AN  = ones(1,N.N_CV_AN);
        ones_vec_SEP = ones(1,N.N_CV_SEP);
        ones_vec_CA  = ones(1,N.N_CV_CA);

        zeros_vec_AN  = zeros(1,N.N_CV_AN);
        zeros_vec_SEP = zeros(1,N.N_CV_SEP);
        zeros_vec_CA  = zeros(1,N.N_CV_CA);

        nan_vec_AN  = nan(1,N.N_CV_AN);
        nan_vec_SEP = nan(1,N.N_CV_SEP);
        nan_vec_CA  = nan(1,N.N_CV_CA);
        
    % Individual Variables
        max_SV = max( N.N_SV_AN , N.N_SV_CA );
    
        SV          = zeros(max_SV+1 , N.N_CV_tot, N_t_steps );
        
        TemperatureK     = zeros( N_t_steps , N.N_CV_tot );

        q_cond           = zeros( N_t_steps , N.N_CV_tot + 1 );
        q_conv           = zeros( N_t_steps , N.N_CV_tot );
        q_gen            = zeros( N_t_steps , N.N_CV_tot );

        intr_egy_CV      = zeros( N_t_steps , N.N_CV_tot ); % internal energy for each control volume
        intr_egy_tot     = zeros( N_t_steps , 1 );
        intr_egy_tot_exp = zeros( N_t_steps , 1 ); % (Expected) Captures the changes in internal energy from heat flux at the boundaries 

        %%%% (Later think about chemical potential energy)
        
        del_phi      = zeros( N_t_steps , N.N_CV_tot );
        phi_ed       = zeros( N_t_steps , N.N_CV_tot );
        phi_el       = zeros( N_t_steps , N.N_CV_tot );
        V_1          = zeros( N_t_steps , N.N_CV_tot );
        V_2          = zeros( N_t_steps , N.N_CV_tot );
        % V_SEI        = zeros( N_t_steps , N.N_CV_tot );
        eta          = zeros( N_t_steps , N.N_CV_tot );
        i_Far        = zeros( N_t_steps , N.N_CV_tot );
        i_o          = zeros( N_t_steps , N.N_CV_tot );
        Eq           = zeros( N_t_steps , N.N_CV_tot );
        CTR_growth   = zeros( N_t_steps , N.N_CV_tot );
        N_Particles  = zeros( N_t_steps , N.N_CV_tot );
        
        i_el         = zeros( N_t_steps , N.N_CV_tot + 1);
        i_ed         = zeros( N_t_steps , N.N_CV_tot + 1);
        % Cap          = zeros( N_t_steps , 1 );
        
        C_Liion      = zeros( N_t_steps , N.N_CV_tot );
        C_Li_surf    = NaN(   N_t_steps , N.N_CV_tot );
        X_Li_surf    = NaN(   N_t_steps , N.N_CV_tot );
        C_Li         = NaN(   N.N_R_max , N.N_CV_tot, N_t_steps );
        X_Li         = NaN(   N.N_R_max , N.N_CV_tot, N_t_steps );
        % J_Liion      = zeros( N_t_steps , N.N_CV_tot + 1);
        % J_Li
        % props
        
        total_mass   = zeros( N_t_steps , 1 );
        CoC          = zeros( N_t_steps , N.N_CV_tot );


%% Calculate i_user
        if SIM.SimMode == 4
            i_user = i_user_soln;
            i_user = i_user(idx);    
        elseif SIM.SimMode == 5
            % Pull i_user from RunSimulation loop and make into a vector the
            % same size as t_soln
        
        elseif SIM.SimMode == 6 % Simulink
            % Pull i_user from sim output
        
        elseif SIM.SimMode == 8 % ---- PRBS ----
            i_user = i_user_soln;
            i_user = i_user(idx);
        
        else
            i_user_in = nan; %%%%%% Placeholder for now
            % i_user = i_user_calc(t_soln,SIM,FLAG,i_user_in);
            i_user = i_user_calc(t_soln, i_user_in, SIM.SimMode, SIM.profile_time, SIM.profile_current, SIM.Amp, SIM, FLAG);
        end
        i_user = reshape(i_user,[],1); % Ensure it is a column vector    


%% Pre-Thermal Calcs
% There are 4 phases being considered: ED, EL, SEP, B (Binder)
    % c_p vector for each control volume
        c_p_vec_ED  = [AN.c_p*ones_vec_AN          zeros_vec_SEP   CA.c_p*ones_vec_CA];
        c_p_vec_EL  = [EL.c_p*ones_vec_AN    EL.c_p*ones_vec_SEP   EL.c_p*ones_vec_CA];
        c_p_vec_SEP = [      zeros_vec_AN   SEP.c_p*ones_vec_SEP         zeros_vec_CA]; 
        c_p_vec_B   = [AN.c_p*ones_vec_AN          zeros_vec_SEP   CA.c_p*ones_vec_CA];

    % rho vector for each control volume
        rho_vec_ED  = [AN.rho*ones_vec_AN          zeros_vec_SEP   CA.rho*ones_vec_CA];
        rho_vec_EL  = [EL.rho*ones_vec_AN    EL.rho*ones_vec_SEP   EL.rho*ones_vec_CA];
        rho_vec_SEP = [      zeros_vec_AN   SEP.rho*ones_vec_SEP         zeros_vec_CA]; 
        rho_vec_B   = [AN.rho*ones_vec_AN          zeros_vec_SEP   CA.rho*ones_vec_CA];

    % dVol vector for each control volume's phase
        dVol_vec_ED  = [AN.eps_ed*AN.dVol*ones_vec_AN                       zeros_vec_SEP   CA.eps_ed*CA.dVol*ones_vec_CA];
        dVol_vec_EL  = [AN.eps_el*AN.dVol*ones_vec_AN   SEP.eps_el *SEP.dVol*ones_vec_SEP   CA.eps_el*CA.dVol*ones_vec_CA];
        dVol_vec_SEP = [                 zeros_vec_AN   SEP.eps_sep*SEP.dVol*ones_vec_SEP                    zeros_vec_CA]; 
        dVol_vec_B   = [AN.eps_b *AN.dVol*ones_vec_AN                       zeros_vec_SEP   CA.eps_b *CA.dVol*ones_vec_CA];

        %%%%% Test on volume calculation (Is there sum equal to the total volume)
            % Combined_matrix_dVol = [dVol_vec_ED ; dVol_vec_EL ; dVol_vec_SEP ; dVol_vec_B];
            % % idx = find( isnan(Combined_matrix_dVol) );
            % % zeros_vec = zeros(1, length(idx) );
            % % Combined_matrix_dVol(idx) = zeros_vec;
            % dVol_sum = sum(Combined_matrix_dVol,1);
            % % dVol_geo = [];
            % ANSEPCA_Vol = [sum(dVol_sum(N.CV_Region_AN))  sum(dVol_sum(N.CV_Region_SEP))  sum(dVol_sum(N.CV_Region_CA))]
            % ANSEPCA_Vol_geo = [ AN.L*AN.A_c , SEP.L*SEP.A_c , CA.L*CA.A_c]
            % vol_tot = sum(ANSEPCA_Vol)
            % vol_tot_geo = AN.A_c*(AN.L+SEP.L+CA.L)

    % cp*rho*dV for each phase's CV
        cp_rho_dV_vec_ED  = c_p_vec_ED  .* rho_vec_ED  .* dVol_vec_ED;
        cp_rho_dV_vec_EL  = c_p_vec_EL  .* rho_vec_EL  .* dVol_vec_EL;
        cp_rho_dV_vec_SEP = c_p_vec_SEP .* rho_vec_SEP .* dVol_vec_SEP;
        cp_rho_dV_vec_B   = c_p_vec_B   .* rho_vec_B   .* dVol_vec_B;


%% Perform calcs and save results to their variable name ~ perform calcs that depend on props
    for i = 1:N_t_steps
        % i
        % if i == 556
        %     i
        % end
        % Go through the solution vector and reshape every SV to 2D (3D matrix)
        % SV_temp    = SV1Dto2D( SV_soln( idx(i) , : ) , N , P , FLAG );
        % SV( : , : , i ) = addPhiEl2SV(SV_temp,P,N);
        % SV_temp = SV1Dto2D(SV_soln(idx(i),:) , N.N_SV_max, N.N_CV_tot, N.N_SV_AN_tot, N.N_SV_SEP_tot, N.N_SV_AN, N.N_SV_SEP, N.N_SV_CA, N.N_CV_AN, N.N_CV_SEP, N.N_CV_CA, N.CV_Region_AN, N.CV_Region_SEP, N.CV_Region_CA, P.T, P.del_phi, P.C_Liion, P.SEP.T, P.SEP.phi_el, P.SEP.C_Liion);
        SV_temp = SV1Dto2D_short(SV_soln(idx(i),:), SIM.SV_nan, N.IDX_1Dto2D);
        SV( : , : , i ) = addPhiEl2SV(SV_temp, P.phi_ed, P.del_phi, N.CV_Region_SEP, N.N_CV_SEP);

        % Pull out each SV vector
        [TemperatureK(i,:), ~, del_phi(i,:), phi_ed(i,:), phi_el(i,:), V_1(i,:), V_2(i,:), ~, C_Liion(i,:), C_Li(:,:,i), C_Li_AN, C_Li_CA, X_AN, X_CA, Ce_norm, ~, eta(i,:), RT_inv_vec , CTR_growth(i,:) , N_Particles(i,:)] = extractSV(SV( : , : , i ),P.T, P.del_phi, P.phi_ed, P.phi_el, P.V_1, P.V_2, P.i_PS, P.CTRGrow, P.NPartic, P.C_Liion, P.C_Li, P.C_Li_surf_AN, P.C_Li_surf_CA, N.CV_Region_AN, N.CV_Region_CA, N.N_R_max, AN.C_Li_max_inv, CA.C_Li_max_inv, EL.C_inv, CONS.R);

        % R_SEI_growth( i , : ) = SV( P.R_SEI_grow , : , i );
        % Cap_loss( i , : )     = SV( P.Cap_loss   , : , i );
        % AM_loss( i , : )      = SV( P.AM_loss    , : , i );
        
        C_Li_surf( i , N.CV_Region_AN ) = C_Li_AN(end,:); 
        C_Li_surf( i , N.CV_Region_CA ) = C_Li_CA(end,:); 

        X_Li(:,:,i) =  C_Li(:,:,i);
        X_Li( : , N.CV_Region_AN, i) = X_Li( : , N.CV_Region_AN , i) / AN.C_Li_max;
        X_Li( : , N.CV_Region_CA, i) = X_Li( : , N.CV_Region_CA , i) / CA.C_Li_max;
        
        X_Li_surf( i , N.CV_Region_AN ) = X_Li( N.N_R_AN , N.CV_Region_AN , i);
        X_Li_surf( i , N.CV_Region_CA ) = X_Li( N.N_R_CA , N.CV_Region_CA , i);
        

        % Get Properties
            if FLAG.VARIABLE_PROPS_FROM_HANDLES
                % props = getProps( SV( : , : , i ) , AN , SEP, CA , EL , P , N , CONS , FLAG , PROPS);
                props = getProps( C_Liion(i,:), TemperatureK(i,:), X_AN,  X_CA, FLAG.VARIABLE_kappa                                                             , ...
                                FLAG.VARIABLE_D_Liion, FLAG.VARIABLE_activity, FLAG.VARIABLE_tf_num, FLAG.VARIABLE_D_o_AN, FLAG.VARIABLE_D_o_CA, FLAG.Bruggeman , ...
                                P.kappa, P.D_o_Li_ion, P.activity, P.tf_num, P.D_o                                                                              , ...
                                N.N_R_AN, N.N_R_CA, N.CV_Region_AN, N.CV_Region_CA                                                                              , ...
                                CONS.BRUG                                                                                                                       , ...
                                AN.D_oHandle, CA.D_oHandle                                                                                                      , ...
                                PROPS, EL.kappaHandle, EL.D_o_Li_ionHandle, EL.ActivityHandle, EL.tf_numHandle                                                       ); 
            else
                props = PROPS;
            end
        
        % Property Vector
            [sigma_vec, kappa_vec, tf_vec, activity_vec, D_o_Li_ion_vec, lambda_vec, D_o_vec] = extractProps(P.sigma, P.kappa, P.tf_num, P.activity, P.D_o_Li_ion, P.lambda_eff, P.D_o, props);
        
        % Calculate Value at Interface between CV
            sigma_vec_interface      = getInterfaceXDir(sigma_vec         , SIM.interp_x_interface);
            kappa_vec_interface      = getInterfaceXDir(kappa_vec         , SIM.interp_x_interface);
            tf_vec_interface         = getInterfaceXDir(tf_vec            , SIM.interp_x_interface);
            activity_vec_interface   = getInterfaceXDir(activity_vec      , SIM.interp_x_interface);
            D_o_Li_ion_vec_interface = getInterfaceXDir(D_o_Li_ion_vec    , SIM.interp_x_interface);
            T_interface              = getInterfaceXDir(TemperatureK(i,:) , SIM.interp_x_interface);
            Ce_interface             = getInterfaceXDir(C_Liion(i,:)      , SIM.interp_x_interface);
        
            D_o_vec_interface        = [NaN(1,N.N_CV_tot) ; 0.5*(D_o_vec(1:end-1,:) + D_o_vec(2:end,:)) ; NaN(1,N.N_CV_tot)]; % This assumes constant del_r in both particles
            dmudc_vec_interface      = dmudc_Latz( Ce_interface , T_interface, tf_vec_interface , activity_vec_interface);

        % Calculate Gradient
            [phi_ed_diff, phi_ed_grad] = diffAndGradXCalc(phi_ed(i,:)       , SIM.diff_CV_x_vec_inv);
            [phi_el_diff, phi_el_grad] = diffAndGradXCalc(phi_el(i,:)       , SIM.diff_CV_x_vec_inv);
            % [Ce_log_diff, Ce_log_grad] = diffAndGradXCalc(Ce_log            , SIM.diff_CV_x_vec_inv);
            [Ce_diff    , Ce_grad    ] = diffAndGradXCalc(C_Liion(i,:)      , SIM.diff_CV_x_vec_inv);
            [T_diff     , T_grad     ] = diffAndGradXCalc(TemperatureK(i,:) , SIM.diff_CV_x_vec_inv);
        
            diff        = C_Li(2:end,:,i) - C_Li(1:end-1,:,i);
            C_Li_diff   = [NaN(1 , N.N_CV_tot) ; diff ; NaN(1 , N.N_CV_tot)];

        % Calculate Fluxes
        % i_el, i_ed
            [i_ed(i , :) , i_el(i , :) ] = currentCalc( sigma_vec_interface, kappa_vec_interface, dmudc_vec_interface, tf_vec_interface, ...
                                                          phi_ed_grad, phi_el_grad, Ce_grad,T_grad                                       , ...
                                                          CONS.F_inv, EL.Beta                                                            , ...
                                                          N.CV_Region_AN, N.CV_Region_SEP, N.CV_Region_CA, N.N_CV_tot, i_user(i)              );

        % i_o
            if FLAG.Newman_i_o
                disp('Fix this FLAG.Newman_i_o in post-processing')
                % i_o_an = CONS.F * AN.k_o ...
                %                    * SV(P.C_Liion     ,N.CV_Region_AN , i)  .^AN.alpha_a ...
                %   .* ( AN.C_Li_max - SV(P.C_Li_surf_AN,N.CV_Region_AN , i) ).^AN.alpha_a ...
                %                   .* SV(P.C_Li_surf_AN,N.CV_Region_AN , i)  .^AN.alpha_c;
                % i_o_ca = CONS.F * CA.k_o ...
                %                    * SV(P.C_Liion     ,N.CV_Region_CA , i)  .^CA.alpha_a ...
                %   .* ( CA.C_Li_max - SV(P.C_Li_surf_CA,N.CV_Region_CA , i) ).^CA.alpha_a ...
                %                   .* SV(P.C_Li_surf_CA,N.CV_Region_CA , i)  .^CA.alpha_c;
            else
                i_o_an  = AN.i_oHandle( TemperatureK(i,N.CV_Region_AN), Ce_norm(:,N.CV_Region_AN), X_Li_surf( i , N.CV_Region_AN ) , AN.i_0_ref , AN.alpha_a , AN.alpha_c);
                i_o_ca  = CA.i_oHandle( TemperatureK(i,N.CV_Region_CA), Ce_norm(:,N.CV_Region_CA), X_Li_surf( i , N.CV_Region_CA ) , CA.i_0_ref , CA.alpha_a , CA.alpha_c);
            end
            i_o(i , :) = [i_o_an, NaN(1,N.N_CV_SEP), i_o_ca];
            
        % i_Far
            i_Far(i , :) = iFarCalc( TemperatureK(i,:), Ce_norm, X_AN, X_CA, eta(i,:), RT_inv_vec, AN.i_0_ref , AN.alpha_a, AN.alpha_c, CA.i_0_ref, CA.alpha_a, CA.alpha_c, ...
                              FLAG.Newman_i_o, CONS.F, ...
                              N.CV_Region_AN, N.CV_Region_CA, N.N_CV_SEP, AN.i_oHandle , CA.i_oHandle);
        
        % Equilibrium
            Eq_an = AN.EqPotentialHandle( X_Li_surf( i , N.CV_Region_AN ) , TemperatureK(i,N.CV_Region_AN));
            Eq_ca = CA.EqPotentialHandle( X_Li_surf( i , N.CV_Region_CA ) , TemperatureK(i,N.CV_Region_CA));
            Eq (i , : ) = [ Eq_an , NaN(1,N.N_CV_SEP) , Eq_ca ];

        % Internal Energy
            intr_egy_vec_ED  = cp_rho_dV_vec_ED  .* TemperatureK( i , : );
            intr_egy_vec_EL  = cp_rho_dV_vec_EL  .* TemperatureK( i , : );
            intr_egy_vec_SEP = cp_rho_dV_vec_SEP .* TemperatureK( i , : );
            intr_egy_vec_B   = cp_rho_dV_vec_B   .* TemperatureK( i , : );

            intr_egy_CV( i , : ) = sum([intr_egy_vec_ED ; intr_egy_vec_EL ; intr_egy_vec_SEP ; intr_egy_vec_B] ,1);
            intr_egy_tot( i )    = sum(intr_egy_CV( i , : ));

        % Expected Internal Energy
            if ~isfield(SIM,'h_AN_BC')
                SIM.h_AN_BC = [];
            end
            if ~isfield(SIM,'h_CA_BC')
                SIM.h_CA_BC = [];
            end
            if ~isfield(SIM,'q_AN_BC')
                SIM.q_AN_BC = [];
            end
            if ~isfield(SIM,'q_CA_BC')
                SIM.q_CA_BC = [];
            end

            [q_cond(i,:) , q_conv(i,:) , q_gen(i,:)] = thermalAllVectors( TemperatureK( i , : ), T_diff, lambda_vec, N.CV_Region_AN, N.N_CV_tot, SIM.del_x_vec_halved_inv, SIM.h, SIM.T_inf, SIM.q_AN_BC, SIM.q_CA_BC, SIM.h_AN_BC, SIM.h_CA_BC, FLAG.T_BC_AN, FLAG.T_BC_CA, FLAG.HEAT_GEN_TOTAL, i_el(i , :), i_ed(i , :), i_Far(i , :));
            if i == 1
                intr_egy_tot_exp(i) = intr_egy_tot(i);
            else
                % Added Energy
                    added_egy = 0;

                % AN BC
                    added_egy = added_egy + q_cond(i-1,1) * AN.A_c * (t_soln(i) - t_soln(i-1)); % [J]

                % CA BC
                    added_egy = added_egy - q_cond(i-1,end) * CA.A_c * (t_soln(i) - t_soln(i-1)); % [J]

                % Outer BC
                    added_egy = added_egy - sum(q_conv(i-1,:) .* SIM.A_outside_vec) * (t_soln(i) - t_soln(i-1)); % [J]

                % Internal????????????? Added electrical power?????
                    % added_egy = added_egy + abs(i_user(i-1,:)*SIM.A_c * cell_voltage(i-1,:)) * (t_soln(i) - t_soln(i-1)); % [J]

                intr_egy_tot_exp(i) = intr_egy_tot_exp(i-1) + added_egy;
            end
    end  
    
    % SEI Voltage
        V_SEI = V_1 - phi_el;

    % Normalized Liion Concentration
        X_Liion = C_Liion / EL.C; % Normalized with respect to the initial concentration

    % Cell Voltage
        cell_voltage = phi_ed(:,end) - phi_ed(:,1);

    % Temperature in Celsius
        TemperatureC = TemperatureK - 273.15;
    
    % Cell Capacity
        Cap = -SIM.A_c/3600 * cumtrapz( t_soln , i_user );   
        
    % SOC
        SOC = 100*Cap/SIM.Cell_Cap + SIM.SOC_start;

    % Current Vectors
        I_user = i_user * SIM.A_c;
        I_user_norm_Crate = I_user / SIM.Cell_Cap;

    % s_dot
        s_dot = i_Far / CONS.F ;

    % Incremental Capacity
        DelVoltage = cell_voltage(3:end) - cell_voltage(1:end-2);
        DelCap     =          Cap(3:end) -          Cap(1:end-2);
        dQdV       = DelCap ./ DelVoltage;

    % Conservation of Charge Check
    % Divergence of the sum of current flux should equal 0
    % #### This equation is now used in the solver in place of the old
    % phi_ed equation. The old phi_ed equation could be used here later but
    % for now I'll stick with this
    for i = 1:N_t_steps
        for j = 1:N.N_CV_tot
            CoC(i , j) = (i_el(i,j+1) + i_ed(i,j+1)) - (i_el(i,j) + i_ed(i,j));
        end
    end


%% Calculate mass at each time step
    % Volume Vector
    Vol_el  = [ AN.dVol_el*ones(1 , N.N_CV_AN) , SEP.dVol_el*ones(1 , N.N_CV_SEP) , CA.dVol_el*ones(1 , N.N_CV_CA)];
    % Vol_ed = NaN(N.N_R_max , N.N_CV_tot);
       
    % for i = N.CV_Region_AN
    %     Vol_ed(1:N.N_R_AN,i) = AN.dVol_r * AN.Np_CV;
    % end
    % for i = N.CV_Region_CA
    %     Vol_ed(1:N.N_R_CA,i) = CA.dVol_r * CA.Np_CV;
    % end
    % Vol_vec = [Vol_el ; Vol_ed];
    % 
    % % Mass in each CV
    % mass = zeros(N.N_R_max+1 , N.N_CV_tot, N_t_steps ); %%%% The plus 1 is to include C_Liion
    % for i = 1:N_t_steps
    %     mass(:,:,i) = Vol_vec .* SV(P.C_Liion:P.C_Li_surf_max,:,i);
    %     % Sum of masses
    %     total_mass(i,1) = sum(mass(:,:,i),'all','omitnan');
    % end
    % mass_error = total_mass(:) - total_mass(1);
       
    dVol_ed = NaN(N.N_R_max , N.N_CV_tot);
    for i = N.CV_Region_AN
        dVol_ed(1:N.N_R_AN,i) = AN.dVol_r;
    end
    for i = N.CV_Region_CA
        dVol_ed(1:N.N_R_CA,i) = CA.dVol_r;
    end
    for i = 1:N_t_steps
        for j = 1:N.N_R_max
            Vol_ed_wrt_time(j,:,i) = dVol_ed(j,:) .* N_Particles(i,:);
        end
        Vol_vec(:,:,i) = [Vol_el ; Vol_ed_wrt_time(:,:,i)];
    end
    
    % Mass in each CV
    mass = zeros(N.N_R_max+1 , N.N_CV_tot, N_t_steps ); %%%% The plus 1 is to include C_Liion
    for i = 1:N_t_steps
        mass(:,:,i) = Vol_vec(:,:,i) .* SV(P.C_Liion:P.C_Li_surf_max,:,i);
        % Sum of masses
        total_mass(i,1) = sum(mass(:,:,i),'all','omitnan');
    end
    mass_error = total_mass(:) - total_mass(1);
    

    %% Calculations specific to sinusoidal pertebations
    if SIM.SimMode == 2 % ---- Harmonic Perturbation ----
        %% ID Voltage Section
        Y           = fft(cell_voltage);   % Fourier Transform of response
        L           = length(Y);           % Number of samples
        P2          = abs(Y/L);            % 2-sided spectrum (neg and pos frequencies)
        P1          = P2(1:L/2+1);         % single-sided spectrum
        P1(2:end-1) = 2*P1(2:end-1);       % Multiply everything by 2 execpt the DC (0 Hz)
        f           = SIM.f_s*(0:(L/2))/L; % [Hz], Frequency domain
        phase       = angle(Y);
        
        V_off_fft      = P1(1);
        [Amp_ID , idx] = max(P1(2:end));
        freq_ID_fft    = f(idx+1);
        ps             = phase(idx+1);
        % ps             = ps + pi/2 - pi; % Adds pi/2 because fft identifies cos, Subtract pi because voltage decreases when current increases
        ps             = ps + pi/2; % Adds pi/2 because fft identifies cos, Subtract pi because voltage decreases when current increases

        % % Test Plot
        %     V_Id = Amp_ID * sin( 2*pi*freq_ID_fft*t_soln + ps ) + V_off_fft;
        %     V_Id2 = Amp_ID * sin( SIM.freq*t_soln + ps ) + V_off_fft;
        % 
        %     figure
        %     hold on
        %     plot(t_soln,cell_voltage, '-k', 'LineWidth',2,'DisplayName','Simulation')
        %     plot(t_soln,V_Id        , 'ob', 'LineWidth',2,'DisplayName','ID')
        %     plot(t_soln,V_Id2       , 'or'               ,'DisplayName','ID2')
        %     title('Compare ID')
        %     xlabel('time (s)')
        %     ylabel('V(t)')
        %     legend
        

        %% Impedance Calculation
        Z_mag = Amp_ID / SIM.I_user_amp; % Impedance Magnitude
        Z_Re = -Z_mag * cos(ps);          % Real Impedance Component
        Z_Im = -Z_mag * sin(ps);          % Imaginary Impedance Component
        Z_dB = 20*log10(Z_mag);          % Impedance Magnitude in decibel
        Z_angle_deg = ps * 360 / (2*pi) - 180; % Phase Shift in degrees
        % OLD
            % Z_Re = Z_mag * cos(ps);          % Real Impedance Component
            % Z_Im = Z_mag * sin(ps);          % Imaginary Impedance Component
            % Z_angle_deg = ps * 360 / (2*pi); % Phase Shift in degrees
        
        if Z_angle_deg <= -358
            Z_angle_deg = Z_angle_deg + 360; % Angle wrapping
        end

    elseif SIM.SimMode == 7 % ---- Manual Profile ----
        if FLAG.Optimize_Profile && FLAG.Save_Current_Profile
            profile_save_filepath   = [filename(1:end-4), '_CurrentProfile_Output.mat'];
            
            region_time_vec    = SIM.region_time_vec;
            region_current_vec = SIM.region_current_vec;
            profile_time       = SIM.profile_time;
            profile_current    = SIM.profile_current;
            t_final            = region_time_vec(end);
            
            SIM.SimMode = 0;
            
            save(profile_save_filepath,'region_time_vec','region_current_vec','profile_time','profile_current','t_final','profile_save_filepath','SIM')
        end
        SIM.SimMode = 7;
    end

    
    %% Set the variable for finished post-processing
        postProcessComplete = 1;
    

    %% Resave data to the .mat file
        clearvars i 
        SIZE_SV = whos('SV');
        if SIZE_SV.bytes > 1e9
            save(filename,'-v7.3')
        else
            save(filename)
        end
end
end
