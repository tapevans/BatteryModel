%% Update Properties
% This function will be used to determine the desired properties as a
% function of temperature and concentration
function props = getProps( SV , AN , SEP, CA , EL , P , N , CONS , FLAG , PROPS)
%% Initialize
    props = zeros( N.N_prop , N.N_CV_tot );
    %%%%%% Maybe add these in the function call
    %%%%%% It's less clean but then there's less variables
    Ce   = SV( P.C_Liion , :);
    T    = SV( P.T , :);
    X_AN = SV( P.C_Li:P.C_Li_surf_AN , N.CV_Region_AN ) / AN.C_Li_max;
    X_CA = SV( P.C_Li:P.C_Li_surf_CA , N.CV_Region_CA ) / CA.C_Li_max;

    
%% Evaluate Properties 
    % kappa
        if FLAG.VARIABLE_kappa
            props( P.kappa , : )      = EL.kappaHandle(Ce , T);
        else
            props( P.kappa , : )      = PROPS( P.kappa , :);
        end
    
    % D_o_Liion
        if FLAG.VARIABLE_D_Liion
            props( P.D_o_Li_ion , : ) = EL.D_o_Li_ionHandle(Ce , T);
        else
            props( P.D_o_Li_ion , : ) = PROPS( P.D_o_Li_ion , :);
        end
    
    % Activity
        if FLAG.VARIABLE_activity
            props( P.activity , : )   = EL.ActivityHandle(Ce , T);
        else
            props( P.activity , : )   = PROPS( P.activity , :);
        end
    
    % Transference Number
        if FLAG.VARIABLE_tf_num
            props( P.tf_num , : )     = EL.tf_numHandle(Ce , T);
        else
            props( P.tf_num , : )     = PROPS( P.tf_num , :);
        end
    
    % Active Material Diffusion -- Anode
        if FLAG.VARIABLE_D_o_AN
            props( P.D_o:end , N.CV_Region_AN ) = AN.D_oHandle(X_AN);
        else
            props( P.D_o:end , N.CV_Region_AN ) = PROPS( P.D_o:end , N.CV_Region_AN );
        end
    
    % Active Material Diffusion -- Cathode
        if FLAG.VARIABLE_D_o_CA
            props( P.D_o:end , N.CV_Region_CA ) = CA.D_oHandle(X_CA);
        else
            props( P.D_o:end , N.CV_Region_CA ) = PROPS( P.D_o:end , N.CV_Region_CA);
        end
    
    % Parameters that don't have handles yet
        props( P.sigma      , N.CV_Region_AN ) = AN.sigma * ones( 1 , N.N_CV_AN );
        props( P.sigma      , N.CV_Region_CA ) = CA.sigma * ones( 1 , N.N_CV_CA );
        props( P.R_SEI      , N.CV_Region_AN ) = AN.R_SEI * ones( 1 , N.N_CV_AN );
        props( P.R_SEI      , N.CV_Region_CA ) = CA.R_SEI * ones( 1 , N.N_CV_CA );
        props( P.lambda_eff , : )              = [AN.lambda_eff * ones(1,N.N_CV_AN) , SEP.lambda_eff * ones(1,N.N_CV_SEP), CA.lambda_eff * ones(1,N.N_CV_CA)];
        props( P.rho_eff    , : )              = [AN.rho_eff    * ones(1,N.N_CV_AN) , SEP.rho_eff    * ones(1,N.N_CV_SEP), CA.rho_eff    * ones(1,N.N_CV_CA)];
        props( P.c_p_eff    , : )              = [AN.c_p_eff    * ones(1,N.N_CV_AN) , SEP.c_p_eff    * ones(1,N.N_CV_SEP), CA.c_p_eff    * ones(1,N.N_CV_CA)];
        
    % Account for tortuosity (Bruggeman)
    if FLAG.Bruggeman
        props = props .* CONS.BRUG;
    end
end